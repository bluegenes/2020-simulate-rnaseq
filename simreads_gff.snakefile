"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s simreads_gff.snakefile --use-conda --configfiles config/Hsapiens_ensembl.yml
"""
# tiny snake to simulate reads from given genome/gene names + gff

#quest_for_orthologs_data: no gff3/bed, so can't use this?

outdir = "output_simreads"
datadir = os.path.join(outdir, "data")
simdir = os.path.join(outdir, "simulate_reads")
logsdir = os.path.join(outdir, "logs")

# use config to specify files, gene names
reffiles = config["reference_files"]

attribute_info=config["attributes"]
attributeD = {}
for attr in attribute_info.keys():
    # build reverse dictionary, so we can get the attribute name from the value
    for val in attribute_info[attr]:
        attributeD[val] = attr

rule all:
    input: expand(os.path.join(simdir, "{ref}", "{gene}", "sample_01.fasta"), gene=attributeD.keys(), ref="Hsapiens_ensembl")

# using R and filtering individually is slow and unnecessary. maybe switch to biopython and grab all attributes with a single passthrough
rule filter_gff:
    input: lambda w: reffiles[w.refname]["gff"]
    output: os.path.join(datadir, "{refname}", "{gene_name}.gff3")
    params:
        attribute_to_filter_on=lambda w: attributeD[w.gene_name],
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

sim_config = config.get("simulation_params")

rule polyester_simreads_gene:
    input: rules.bedtools_getfasta.output
    output: os.path.join(simdir, "{refname}", "{gene_name}", "sample_01.fasta")
    params:
        output_dir = lambda w: os.path.abspath(os.path.join(simdir, w.refname, w.gene_name)),
        num_reps = sim_config.get("num_reps", 5),
        read_length = sim_config.get("read_length", 150),
        simulate_paired = sim_config.get("paired", False),
        num_reads_per_transcript=sim_config.get("num_reads_per_transcript", 1000),
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

