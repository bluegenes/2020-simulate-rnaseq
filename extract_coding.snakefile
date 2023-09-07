"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s extract_coding.snakefile --use-conda --configfiles config/QfO.yml
"""
# tiny snake to simulate reads from given genome/gene names + gff

#quest_for_orthologs_data: no gff3/bed - currently just simulate from entire file. Might want to subset somehow instead...

outdir = config("output_directory", "output_translate")
datadir = os.path.join(outdir, "data")
simdir = os.path.join(outdir, "simulate_reads")
logsdir = os.path.join(outdir, "logs")



attribute_info=config.get("attributes")
attribute_targets=[]
if any([not attribute_info, full]):
    full_target = expand(os.path.join(simdir, "{ref}", "sample_01.fasta"), ref=reffiles.keys()), # full sequence target
if attribute_info:
    attributeD = {}
    for attr in attribute_info.keys():
        # build reverse dictionary, so we can get the attribute name from the value
        for val in attribute_info[attr]:
            attributeD[val] = attr
    attribute_targets=expand(os.path.join(simdir, "{ref}", "{gene}", "sample_01.fasta"), gene=attributeD.keys(), ref=reffiles.keys())
    # add flanking sequences (or slop)
    if flank:
        attribute_targets+=expand(os.path.join(simdir, "{ref}", "{gene}", "flank{flank}", "sample_01.fasta"), gene=attributeD.keys(), ref=reffiles.keys(), flank=flank)
    if slop:
        attribute_targets+=expand(os.path.join(simdir, "{ref}", "{gene}", "slop{slop}", "sample_01.fasta"), gene=attributeD.keys(), ref=reffiles.keys(), slop=slop)

rule all:
    input:
        full_target,
        attribute_targets # specified above


rule seqtk_fasta_to_fastq:
    input: lambda w: sampleD[w.sample] # sample name: fasta dictionary
    output: os.path.join(simdir, "{refname}", "{sample}.fq.gz")
    #log: os.path.join(logsdir, "{refname}.simreads.log")
    #benchmark: os.path.join(logsdir, "{refname}.simreads.benchmark")
    ##wildcard_constraints:
    #    refname="\w+", # ~only allow letters,numbers, underscore
    #    gene_name="\w+"
    conda: "polyester-env.yml"
    shell:
    """
    seqtk seq -F '#' {input} | gzip -9 > {output}
    """

moltypeD = {"hp": "hydrophobic-polar", "protein": "protein", "dayhoff": "dayhoff"}
pepRefD = {"sprot": "/home/ntpierce/2020-pep/khtools_testing/uniprot_sprot.fasta.gz", "merc": "/home/ntpierce/2020-pep/khtools_testing/MERC.fasta.gz"}



# translate reads --> proteins
rule extract_coding:
    input:
        reads= rules.seqtk_fasta_to_fastqc
        pep_ref= lambda w: pepRefD[w.pep_ref]
    output:
        coding_prot=os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.codingpep.fa"),
        noncoding_nucl=os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.noncoding.fa"),
        low_complexity_nucl=os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.lowcomplexnucl.fa"),
        csv=os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.csv"),
    log: os.path.join(logs_dir, "khtools", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.extract_coding.log")
    benchmark: os.path.join(logs_dir, "khtools", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.extract_coding.benchmark")
    params:
        molecule=lambda w: moltypeD[w.molecule],
    #wildcard_constraints:
    #    ksize=["\d+"],
    conda: os.path.join(wrappers_dir, "khtools-env.yml")
    shell:
        """
        khtools extract-coding --verbose --molecule {params.molecule} --peptide-ksize {wildcards.ksize} --noncoding-nucleotide-fasta {output.noncoding_nucl} --low-complexity-nucleotide-fasta {output.             low_complexity_nucl} --csv {output.csv} {input.pep_ref} {input.reads} > {output.coding_prot} 2> {log}
        """

rule diamond_makedb_nr:
    input: "/group/ctbrowngrp/pierce/databases/nr.gz"
    output: "/group/ctbrowngrp/pierce/databases/nr.dmnd"
    conda: os.path.join(wrappers_dir, "diamond-env.yml")
    log: os.path.join(logs_dir, "diamond", "nr_makedb.log")
    benchmark: os.path.join(logs_dir, "diamond", "nr_makedb.benchmark")
    shell:
        """
        diamond makedb --in {input} --db {output} 2> {log}
        """

# for orthofinder, output needs to be of form: Blast{species}_{species}.txt.gz
rule diamond_blastp_extract_coding:
    input:
        pep = rules.extract_coding.output.coding_prot,
        db = rules.diamond_makedb_nr.output
    output: os.path.join(out_dir, "preprocess", "khtools", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.codingpep.diamond_nr.out")
    conda: os.path.join(wrappers_dir, "diamond-env.yml")
    log: os.path.join(logs_dir, "diamond", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.codingpep.diamond_nr.log")
    benchmark: os.path.join(logs_dir, "diamond", "{sample}_{molecule}_k{ksize}_ref{pep_ref}.codingpep.diamond_nr.benchmark")
    shell:
        """
        diamond blastp -d {input.db} -q {input.pep} -o {output} 2> {log}
        """
