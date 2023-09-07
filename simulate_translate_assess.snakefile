"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s simreads.snakefile --use-conda --configfiles config/Hsapiens_ensembl.yml
"""
# simulate RNAseq reads from transcripts, translate with sencha, assess translation

import os
import sys
import pandas as pd

# set up dirs
envs_dir = "envs"

outdir = config.get("output_directory","output_simreads")
simdir = os.path.join(outdir, "simulate_reads")
logsdir = os.path.join(outdir, "logs")
sencha_indexdir=config.get("sencha_indexdir", "index")
index_dir = os.path.join(outdir, sencha_indexdir)
translate_dir = os.path.join(outdir, "translate")

reference_info = config["sencha_params"].get("peptide_references")

## polyester creates files called "sample_01.fasta.gz", sample_02.fasta.gz", etc
def generate_replist(number_replicates):
    replist = list(range(number_replicates))
    reps=[]
    for num in replist:
        num+=1 # start at 1
        if num < 10:
            num = "0" + str(num)
        reps.append(str(num))
    return reps

# get polyester parameters from config; use replicates to build rep sample numbers
sim_config = config.get("simulation_params")
replist = generate_replist(sim_config.get("num_replicates", 5))

fastaDF = pd.read_csv(config.get("fasta_csv"), dtype=str, sep=",", header=0, index_col=0)
#sampleDF = pd.read_csv(config.get("samples_csv"), dtype=str, sep=",", header=0, index_col=0)

# build output files based on config
output_targets = []
ref_targets = []
output_extensions = ["codingpep.fa", "noncoding.fa", "lowcomplexnucl.fa", "csv"]

alphabets_to_run = config["sencha_params"]["alphabet"] # dictionary of alphabets to run, corresponding ksizes 
# build final file targets based on config info
for alpha, alpha_info in alphabets_to_run.items():
    output_targets+=expand(os.path.join(translate_dir, "{sample}_{rep}_{alphabet}_k{k}_ref{ref}.{ext}"), sample=list(fastaDF.index), alphabet=alpha, k=alpha_info["ksizes"], ref=reference_info.keys(), ext=output_extensions)
    #ref_targets+=expand(os.path.join(index_dir, "ref{ref}_{alphabet}_k{k}.index"), ref=reference_info.keys(), alphabet=alpha, k=alpha_info["ksizes"])
    ref_targets+=expand(os.path.join(index_dir, "ref{ref}_{alphabet}_k{k}.nodegraph"), ref=reference_info.keys(), alphabet=alpha, k=alpha_info["ksizes"])

# abbreviate "hydrophobic-polar" as hp in filenames. Need full alpha name for khtools code
alpha_abbreviations = {"hp": "hydrophobic-polar", "protein": "protein", "dayhoff": "dayhoff"}

rule all:
    input: 
        expand(os.path.join(simdir, "{run_name}_samples.csv"), run_name = runname),
        output_targets

rule polyester_simreads_full:
    input: lambda w: samplesDF["fasta"][w.refname]
    input: lambda w: fastaDF.loc[w.refname, "read1"],
    output: expand(os.path.join(simdir, "{{refname}}", "sample_{rep}.fasta.gz"), rep = replist)
    params:
        output_dir = lambda w: os.path.abspath(os.path.join(simdir, w.refname)),
        num_reps = sim_config.get("num_replicates", 5),
        read_length = sim_config.get("read_length", 150),
        simulate_paired = sim_config.get("paired", False),
        num_reads_per_transcript=sim_config.get("num_reads_per_transcript", 1000),
    log: os.path.join(logsdir, "{refname}.simreads.log")
    benchmark: os.path.join(logsdir, "{refname}.simreads.benchmark")
    threads: 1
    resources:
      mem_mb=16000, #16GB
      runtime=1000
    wildcard_constraints:
        refname="\w+",
        gene_name="\w+"
    conda: "envs/polyester-env.yml"
    script: "scripts/simulate_reads.R"

rule seqtk_fasta_to_fastq_full:
    input: os.path.join(simdir, "{refname}", "sample_{rep}.fasta.gz")
    output: os.path.join(simdir, "{refname}", "{refname}_{rep}.fq.gz")
    log: os.path.join(logsdir, "seqtk", "{refname}_{rep}.seqtk.log")
    benchmark: os.path.join(logsdir, "seqtk", "{refname}_{rep}.seqtk.benchmark")
    threads: 1
    resources:
      mem_mb=4000, #4GB
      runtime=60 #minutes
    wildcard_constraints:
        refname="\w+",
        gene_name="\w+",
        rep="\d+"
    conda: "envs/seqtk-env.yml"
    shell:
        """
        seqtk seq -F 'I' {input} | gzip -9 > {output}
        """

rule write_samples_csv:
    input: expand(os.path.join(simdir, "{ref}", "{ref}_{rep}.fq.gz"), ref=refLinks.keys(), rep=replist),
    output: os.path.join(simdir, "{run_name}_samples.csv")
    threads: 1
    resources:
      mem_mb=1000,
      runtime=15
    run:
        with open(str(output), "w") as csv:
            csv.write("sample_id" + "," + "read_1"+ "\n")
            for f in input:
                refname = os.path.basename(os.path.dirname(f))
                replicate = f.rsplit(".fq.gz")[0].rsplit("_", 1)[1] # _01.fq.gz
                csv.write(f"{refname}_{replicate}" + "," + os.path.abspath(f) + "\n")


rule sencha_index:
    input: lambda w: reference_info[w.ref] # get file path from reference name, config
    output: os.path.join(index_dir, "ref{ref}_{alphabet}_k{ksize}.nodegraph"),
    log: os.path.join(logs_dir, "sencha_index", "ref{ref}_{alphabet}_k{ksize}.index.log")
    benchmark: os.path.join(logs_dir, "sencha_index", "ref{ref}_{alphabet}_k{ksize}.index.benchmark")
    params:
        alphabet=lambda w: alpha_abbreviations[w.alphabet],
    threads: 1
    resources:
        mem_mb=80000,
        runtime=6000 # ~4 days 
    wildcard_constraints:
        ref="\w+",
        alphabet="\w+",
        ksize="\d+",
    conda: os.path.join(envs_dir, "khtools-env.yml")
    shell:
        """
        khtools index --alphabet {params.alphabet} --peptide-ksize {wildcards.ksize} {input} --save-as {output} 2> {log}
        """

rule sencha_translate:
    input:
        #rules.write_samples_csv.output
        #fastq=lambda w: sampleDF.loc[w.sample, "read1"],
        fastq=os.path.join(simdir, "{refname}", "{refname}_{rep}.fq.gz"),
        index= rules.sencha_index.output
    output:
        coding_prot=os.path.join(translate_dir, "{refname}_{rep}_{alphabet}_k{ksize}_ref{ref}.codingpep.fa"),
        coding_nucl=os.path.join(translate_dir, "{refname}_{rep}_{alphabet}_k{ksize}_ref{ref}.codingnucl.fa"),
        noncoding_nucl=os.path.join(translate_dir, "{refname}_{rep}_{alphabet}_k{ksize}_ref{ref}.noncoding.fa"),
        low_complexity_prot=os.path.join(translate_dir, "{refname}_{rep}_{alphabet}_k{ksize}_ref{ref}.lowcomplexprot.fa"),
        low_complexity_nucl=os.path.join(translate_dir, "{refname}_{rep}_{alphabet}_k{ksize}_ref{ref}.lowcomplexnucl.fa"),
        csv=os.path.join(translate_dir, "{refname}_{rep}_{alphabet}_k{ksize}_ref{ref}.csv"),
    #log: os.path.join(logs_dir, "sencha_translate", "{sample}_{alphabet}_k{ksize}_ref{ref}.translate.log") #2>{log} err ("missing PEPTIDES file")
    benchmark: os.path.join(logs_dir, "sencha_translate", "{refname}_{rep}_{alphabet}_k{ksize}_ref{ref}.translate.benchmark")
    params:
        alphabet=lambda w: alpha_abbreviations[w.alphabet],
    threads: 1
    resources:
        mem_mb=80000, #lambda wildcards, attempt: attempt *40000, #16GB*attempt
        runtime=6000
    wildcard_constraints:
        ref="\w+",
        alphabet="\w+",
        ksize="\d+",
    conda: os.path.join(envs_dir, "khtools-env.yml")
    shell:
        """
        khtools translate --verbose --peptides-are-bloom-filter --alphabet {params.alphabet} --peptide-ksize {wildcards.ksize} --noncoding-nucleotide-fasta {output.noncoding_nucl} --low-complexity-nucleotide-fasta {output.low_complexity_nucl} --coding-nucleotide-fasta {output.coding_nucl} --low-complexity-peptide-fasta {output.low_complexity_prot} --csv {output.csv} {input.index} {input.fastq} > {output.coding_prot}
        """
