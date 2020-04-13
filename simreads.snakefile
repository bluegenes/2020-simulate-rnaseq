"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s simreads.snakefile --use-conda --configfiles config/Hsapiens_ensembl.yml
"""
# tiny snake to simulate reads from given genome/gene names + gff

#quest_for_orthologs_data: no gff3/bed - currently just simulate from entire file. Might want to subset somehow instead...
outdir = config.get("output_directory","output_simreads")
datadir_default = os.path.join(outdir, "data")
datadir = config.get("data_directory", datadir_default)
simdir = os.path.join(outdir, "simulate_reads")
logsdir = os.path.join(outdir, "logs")

# use config to specify files, gene names
reffiles = config["reference_files"]

flank = config.get("flank_bp")
slop = config.get("slop_bp")

full = config.get("full", False)
full_target=[]

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

attribute_info=config.get("attributes")
attribute_targets=[]
if any([not attribute_info, full]):
    full_target = expand(os.path.join(simdir, "{ref}", "{ref}_{rep}.fq.gz"), ref=reffiles.keys(), rep=replist), # full sequence target
if attribute_info:
    attributeD = {}
    for attr in attribute_info.keys():
        # build reverse dictionary, so we can get the attribute name from the value in rule params
        for val in attribute_info[attr]:
            attributeD[val] = attr
    attribute_targets=expand(os.path.join(simdir, "{ref}", "{gene}", "{gene}_{rep}.fq.gz"), gene=attributeD.keys(), ref=reffiles.keys(), rep=replist)
    # add flanking sequences (or slop)
    if flank:
        attribute_targets+=expand(os.path.join(simdir, "{ref}", "{gene}", "flank{flank}", "{gene}_flank{flank}_{rep}.fq.gz"), gene=attributeD.keys(), ref=reffiles.keys(), flank=flank, rep=replist)
    if slop:
        attribute_targets+=expand(os.path.join(simdir, "{ref}", "{gene}", "slop{slop}", "{gene}_slop{slop}_{rep}.fq.gz"), gene=attributeD.keys(), ref=reffiles.keys(), slop=slop, rep=replist)

rule all:
    input: 
        full_target,
        attribute_targets # specified above

# using R and filtering individually is slow and unnecessary. maybe switch to biopython and grab all attributes with a single passthrough
rule filter_gff:
    input: 
        gff= lambda w: os.path.join(datadir, reffiles[w.refname]["gff"]),
    output: os.path.join(datadir, "{refname}", "filtered", "{gene_name}.gff3")
    params:
        attribute_to_filter_on=lambda w: attributeD[w.gene_name],
        attribute_value=lambda w: w.gene_name,
    log: os.path.join(logsdir, "filter_gff", "{refname}_{gene_name}.filter_gff.log")
    wildcard_constraints:
        refname="\w+", # ~only allow letters,numbers, underscore
        gene_name="\w+"
    conda: "envs/polyester-env.yml"
    script: "scripts/filter_gff.R"

# flank gets the regions around the features, excluding the features themselves. 
rule bedtools_flank:
    input:
        fasta= lambda w: os.path.join(datadir, reffiles[w.refname]["fasta"]),
        genome_info= lambda w: os.path.join(datadir, w.refname, reffiles[w.refname]["genome_info"]),
        gff=rules.filter_gff.output
    output: 
        gff=os.path.join(datadir, "{refname}", "filtered", "{gene_name}.flank{flank}.gff"),
        fasta=os.path.join(datadir, "{refname}", "filtered", "{gene_name}.flank{flank}.fasta")
    log: os.path.join(logsdir, "bedtools", "{refname}", "{gene_name}.flank{flank}.get_fasta.log")
    benchmark: os.path.join(logsdir, "bedtools", "{refname}", "{gene_name}.flank{flank}.get_fasta.benchmark")
    params:
        flank= lambda w: w.flank,
    wildcard_constraints:
        refname="\w+",
        gene_name="\w+"
    conda: "envs/bedtools-env.yml"
    shell:
        """
        bedtools flank -b {params.flank} -i {input.gff} -g {input.genome_info} > {output.gff} 2> {log}
        bedtools getfasta -fi {input.fasta} -bed {output.gff} > {output.fasta} 2>> {log}
        """

# flank/slop options:
#-b	Increase the BED/GFF/VCF entry by the same number base pairs in each direction. Integer.
#-l	The number of base pairs to subtract from the start coordinate. Integer.
#-r	The number of base pairs to add to the end coordinate. Integer.
#-s	Define -l and -r based on strand. For example. if used, -l 500 for a negative-stranded feature, it will add 500 bp to the end coordinate.
#-pct	Define -l and -r as a fraction of the feature’s length. E.g. if used on a 1000bp feature, -l 0.50, will add 500 bp “upstream”. Default = false.

# slop gets the regions around the features INCLUDING the features themselves
rule bedtools_slop:
    input:
        fasta= lambda w: os.path.join(datadir, reffiles[w.refname]["fasta"]),
        genome_info=lambda w: os.path.join(datadir, w.refname, reffiles[w.refname]["genome_info"]),
        gff=rules.filter_gff.output
    output:
        gff=os.path.join(datadir, "{refname}", "filtered", "{gene_name}.slop{slop}.gff"),
        fasta=os.path.join(datadir, "{refname}", "filtered", "{gene_name}.slop{slop}.fasta")
    log: os.path.join(logsdir, "bedtools", "{refname}", "{gene_name}.slop{slop}.get_fasta.log")
    benchmark: os.path.join(logsdir, "bedtools", "{refname}", "{gene_name}.slop{slop}.get_fasta.benchmark")
    params:
        slop= lambda w: w.slop,
        #genome_info= lambda w: os.path.join(datadir, w.refname, reffiles[w.refname]["genome_info"]),
    wildcard_constraints:
        refname="\w+",
        gene_name="\w+"
    conda: "envs/bedtools-env.yml"
    shell:
        """
        bedtools slop -b {params.slop} -i {input.gff} -g {input.genome_info} > {output.gff} 2> {log}
        bedtools getfasta -fi {input.fasta} -bed {output.gff} > {output.fasta} 2>> {log}
        """

rule bedtools_getfasta:
    input: 
        fasta= lambda w: os.path.join(datadir, reffiles[w.refname]["fasta"]),
        #fasta=os.path.join(datadir, "{refname}", "{refname}.fasta"),
        gff=rules.filter_gff.output
    output: os.path.join(datadir, "{refname}", "filtered", "{gene_name}.fasta")
    log: os.path.join(logsdir, "bedtools", "{refname}", "{gene_name}.get_fasta.log")
    benchmark: os.path.join(logsdir, "bedtools", "{refname}", "{gene_name}.get_fasta.benchmark")
    wildcard_constraints:
        refname="\w+",
        gene_name="\w+"
    conda: "envs/bedtools-env.yml"
    shell: 
        """
        bedtools getfasta -fi {input.fasta} -bed {input.gff} > {output} 2> {log}
        """


# simulation rules
rule polyester_simreads_gene:
    input: rules.bedtools_getfasta.output
    output: expand(os.path.join(simdir, "{{refname}}", "{{gene_name}}", "sample_{rep}.fasta.gz"), rep = replist)
    params:
        output_dir = lambda w: os.path.join(simdir, w.refname, w.gene_name),
        num_reps = sim_config.get("num_replicates", 5),
        read_length = sim_config.get("read_length", 150),
        simulate_paired = sim_config.get("paired", False),
        num_reads_per_transcript=sim_config.get("num_reads_per_transcript", 1000),
    wildcard_constraints:
        refname="\w+",
        gene_name="\w+"
    log: os.path.join(logsdir, "{refname}_{gene_name}.simreads.log")
    benchmark: os.path.join(logsdir, "{refname}_{gene_name}.simreads.benchmark")
    conda: "envs/polyester-env.yml" 
    script: "scripts/simulate_reads.R"

rule polyester_simreads_flank:
    input: rules.bedtools_flank.output.fasta 
    output: expand(os.path.join(simdir, "{{refname}}", "{{gene_name}}", "flank{{flank}}", "sample_{rep}.fasta.gz"), rep = replist)
    params:
        output_dir = lambda w: os.path.join(simdir, w.refname, w.gene_name, f"flank{w.flank}"),
        num_reps = sim_config.get("num_replicates", 5),
        read_length = sim_config.get("read_length", 150),
        simulate_paired = sim_config.get("paired", False),
        num_reads_per_transcript=sim_config.get("num_reads_per_transcript", 1000),
    log: os.path.join(logsdir, "{refname}_{gene_name}.flank{flank}.simreads.log")
    benchmark: os.path.join(logsdir, "{refname}_{gene_name}.flank{flank}.simreads.benchmark")
    wildcard_constraints:
        refname="\w+",
        gene_name="\w+"
    conda: "envs/polyester-env.yml"
    script: "scripts/simulate_reads.R"

rule polyester_simreads_slop:
    input: rules.bedtools_slop.output.fasta 
    output: expand(os.path.join(simdir, "{{refname}}", "{{gene_name}}", "slop{{slop}}", "sample_{rep}.fasta.gz"), rep = replist)
    params:
        output_dir = lambda w: os.path.join(simdir, w.refname, w.gene_name, f"slop{w.slop}"),
        num_reps = sim_config.get("num_replicates", 5),
        read_length = sim_config.get("read_length", 150),
        simulate_paired = sim_config.get("paired", False),
        num_reads_per_transcript=sim_config.get("num_reads_per_transcript", 1000),
    log: os.path.join(logsdir, "{refname}_{gene_name}.slop{slop}.simreads.log")
    benchmark: os.path.join(logsdir, "{refname}_{gene_name}.slop{slop}.simreads.benchmark")
    wildcard_constraints:
        refname="\w+",
        gene_name="\w+"
    conda: "envs/polyester-env.yml"
    script: "scripts/simulate_reads.R"

rule polyester_simreads_full:
    input: lambda w: os.path.join(datadir, reffiles[w.refname]["fasta"])
    output: expand(os.path.join(simdir, "{{refname}}", "sample_{rep}.fasta.gz"), rep = replist)
    params:
        output_dir = lambda w: os.path.abspath(os.path.join(simdir, w.refname)),
        num_reps = sim_config.get("num_replicates", 5),
        read_length = sim_config.get("read_length", 150),
        simulate_paired = sim_config.get("paired", False),
        num_reads_per_transcript=sim_config.get("num_reads_per_transcript", 1000),
    log: os.path.join(logsdir, "{refname}.simreads.log")
    benchmark: os.path.join(logsdir, "{refname}.simreads.benchmark")
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
    wildcard_constraints:
        refname="\w+",
        gene_name="\w+",
        rep="\d+"
    conda: "envs/seqtk-env.yml"
    shell:
        """
        seqtk seq -F 'I' {input} | gzip -9 > {output}
        """

rule seqtk_fasta_to_fastq_attribute:
    input: os.path.join(simdir, "{refname}", "{gene_name}", "sample_{rep}.fasta.gz")
    output: os.path.join(simdir, "{refname}", "{gene_name}", "{gene_name}_{rep}.fq.gz")
    log: os.path.join(logsdir, "seqtk", "{refname}_{gene_name}_{rep}.seqtk.log")
    benchmark: os.path.join(logsdir, "seqtk", "{refname}_{gene_name}_{rep}.seqtk.benchmark")
    wildcard_constraints:
        refname="\w+",
        gene_name="\w+",
        rep="\d+"
    conda: "envs/seqtk-env.yml"
    shell:
        """
        seqtk seq -F 'I' {input} | gzip -9 > {output}
        """

rule seqtk_fasta_to_fastq_extrainfo:
    input: os.path.join(simdir, "{refname}", "{gene_name}",  "{extrainfo}", "sample_{rep}.fasta.gz")
    output: os.path.join(simdir, "{refname}","{gene_name}",  "{extrainfo}", "{gene_name}_{extrainfo}_{rep}.fq.gz")
    log: os.path.join(logsdir, "seqtk", "{refname}_{gene_name}_{extrainfo}_{rep}.seqtk.log")
    benchmark: os.path.join(logsdir, "seqtk", "{refname}_{gene_name}_{extrainfo}_{rep}.seqtk.benchmark")
    wildcard_constraints:
        refname="\w+",
        gene_name="\w+",
        extrainfo="\w+",
        rep="\d+"
    conda: "envs/seqtk-env.yml"
    shell:
        """
        seqtk seq -F 'I' {input} | gzip -9 > {output}
        """
