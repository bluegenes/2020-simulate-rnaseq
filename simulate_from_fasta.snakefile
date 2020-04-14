"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s simreads.snakefile --use-conda --configfiles config/Hsapiens_ensembl.yml
"""
# simulate RNAseq reads from transcripts in a fasta file 

import os

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

outdir = config.get("output_directory","output_simreads")
datadir_default = os.path.join(outdir, "data")
datadir = config.get("data_directory", datadir_default)
simdir = os.path.join(outdir, "simulate_reads")
logsdir = os.path.join(outdir, "logs")

# use config to specify files
runname = config["run_name"]
reffiles = config["reference_files"]
# download and simulate these reference types per species
default_reftypes = ["ncRNA", "cds", "cdna", "cdna_abinitio", "qfo_dna"]
reference_types = config.get("reference_types", default_reftypes)

# if using different types of reference files (e.g. noncoding rna vs coding seqs), build a refname: download link dictionary
refLinks={}
for reference in reffiles:
    for reftype, link in reffiles[reference].items():
        refname = f"{reference}_{reftype}"
        refLinks[refname] = link # build reference_name : download link dictionary

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

rule all:
    input: 
        expand(os.path.join(simdir, "{run_name}_samples.csv"), run_name = runname)


rule download_fasta:
    input: lambda w: FTP.remote(refLinks[w.refname], static=True, keep_local=True, immediate_close=True)
    output: os.path.join(datadir,"{refname}.fa.gz")
    log: os.path.join(logsdir,"get_data", "{refname}.log")
    shell: "mv {input} {output} 2> {log}"

rule polyester_simreads_full:
    #input: lambda w: os.path.join(datadir, reffiles[w.refname]["fasta"])
    input: lambda w: os.path.join(datadir,"{refname}.fa.gz") 
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

rule write_samples_csv:
    input: expand(os.path.join(simdir, "{ref}", "{ref}_{rep}.fq.gz"), ref=refLinks.keys(), rep=replist),
    output: os.path.join(simdir, "{run_name}_samples.csv")
    run:
        with open(output, "w") as csv:
            csv.write("sample_id" + "," + "read_1"+ "\n")
            for f in input:
                refname = os.path.dirname(f)
                csv.write(refname + "," + os.path.abspath(f) + "\n")
