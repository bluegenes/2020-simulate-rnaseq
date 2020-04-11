"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s simreads.snakefile --use-conda --configfiles config/Hsapiens_ensembl.yml
"""
# download reference files if necessary

# use config to specify files, gene names
reffiles = config["reference_files"]

# same directories as simreads.snakefile. If running as subworkflow, these will already be present
outdir = config.get("output_directory","output_simreads")
datadir = os.path.join(outdir, "data")
simdir = os.path.join(outdir, "simulate_reads")
logsdir = os.path.join(outdir, "logs")


rule all:
    input: expand(os.path.join(datadir, "{refname}", "{refname}.{filetype}"), refname = reffiles.keys(), filetype= ["fasta", "gff"])


rule download_tar:
    output:
        fasta=os.path.join(datadir, "{refname}", "{refname}.fasta"),
        gff=os.path.join(datadir, "{refname}", "{refname}.gff"),
    params:
        # url and fileformat for full database (optional)
        database_url = lambda w: reffiles[w.refname].get("database_url"),
        database_fileformat = lambda w: reffiles[w.refname].get("database_fileformat", "tar.gz"),
        # relative paths to fasta, gff (gff optional)
        gff_path = lambda w: reffiles[w.refname].get("gff"),
        # independent fasta, gff downloads (optional)
        fasta_url = lambda w: reffiles[w.refname].get("fasta_url"),
        gff_url = lambda w: reffiles[w.refname].get("gff_url"),
        fasta_filetype = lambda w: reffiles[w.refname].get("fasta_fileformat", "gz"),
        gff_filetype = lambda w: reffiles[w.refname].get("gff_fileformat", "gz"),
    log: os.path.join(logsdir, '{refname}.download.log')
    wildcard_constraints:
        refname="\w+", # ~only allow letters,numbers, underscore
        gene_name="\w+"
    script: "scripts/download_wrapper.py"

#rule download_gz:
#    output:

