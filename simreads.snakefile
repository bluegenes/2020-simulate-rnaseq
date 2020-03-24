# tiny snake to simulate reads from given genome/gene names + gff

#quest_for_orthologs_data: no gff3/bed, so can't use this?
#human genome data:

outdir = "output_simreads"
datadir = os.path.join(outdir, "data")
simdir = os.path.join(outdir, "simulate_reads")
logsdir = os.path.join(outdir, "logs")

reffiles = {}
reffiles["Hsapiens_ensembl"] = {"fasta": "/Users/tessa/dib-lab/simulate_reads/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa",
                                "gff": "/Users/tessa/dib-lab/simulate_reads/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf"}

reference_name="Hsapiens_ensembl"

GENES = ["GNAS", "AGER", "CDKN2A", "CSF2RA", "GNRHR", "XBP1"]

rule all:
    input: expand(os.path.join(simdir, "{ref}", "{gene}", "sample_01.fasta"), gene=GENES, ref="Hsapiens_ensembl")

rule filter_gff:
    input: lambda w: reffiles[w.refname]["gff"]
    output: os.path.join(datadir, "{refname}", "{gene_name}.gff3")
    params:
        attribute_to_filter_on="gene_name",
        attribute_value=lambda w: w.gene_name,
    log: os.path.join(logsdir, "{refname}", "{gene_name}.filter_gff.log")
    conda: "polyester-env.yml"
    script: "filter_gff.R"

# right now we're getting the exact matching range. to do: add flank, slop to get regions around these genes
rule bedtools_getfasta:
    input: 
        fasta= lambda w: reffiles[w.refname]["fasta"],
        gff=rules.filter_gff.output
    output: os.path.join(datadir, "{refname}", "{gene_name}.fasta")
    log: os.path.join(logsdir, "bedtools", "{refname}", "{gene_name}.get_fasta.log")
    benchmark: os.path.join(logsdir, "bedtools", "{refname}", "{gene_name}.get_fasta.benchmark")
    conda: "bedtools-env.yml"
    shell: 
        """
        bedtools getfasta -fi {input.fasta} -bed {input.gff} > {output} 2> {log}
        """

rule polyester_simreads_gene:
    input: rules.bedtools_getfasta.output
    output: os.path.join(simdir, "{refname}", "{gene_name}", "sample_01.fasta")
    params:
        output_dir = lambda w: os.path.abspath(os.path.join(simdir, w.refname, w.gene_name)),
        num_reps = 5,
        read_length = 150,
        simulate_paired = False,
        num_reads_per_transcript=1000,
    log: os.path.join(logsdir, "{refname}_{gene_name}.simreads.log")
    benchmark: os.path.join(logsdir, "{refname}_{gene_name}.simreads.benchmark")
    conda: "polyester-env.yml" 
    script: "simulate_reads.R"

#rule download_QfO_2018:
#    output: os.path.join(datadir, "QfO_release_2018_04.tar.gz")
#    log: os.path.join(logs_dir, "QfO.download")
#    shell:
#        """
#        curl -L ftp://ftp.ebi.ac.uk/pub/databases/reference_proteomes/previous_releases/qfo_release-2018_04/QfO_release_2018_04.tar.gz -O {output}
#        """

