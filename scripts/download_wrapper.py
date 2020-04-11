__author__ = "N. Tessa Pierce"
__copyright__ = "Copyright 2020 N. Tessa Pierce"
__email__ = "ntpierce@ucdavis.edu"
__license__ = "MIT"


import os

from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)

cmd = ['curl', '-L']

def get_unzip_cmd(fileformat, outfile, zipped=None):
    if fileformat == "gz":
        return ['|', 'gunzip', '-c', '>', outfile]
    elif fileformat == "tar.gz":
        return ["-O", zipped]#, ";" 'tar', '-xzf', zipped, "--directory" , outfile, "--strip-components 1"] ##, '>', outfile]
    else:
        raise ValueError('Valid url_fileformats are "gz" and "tar.gz"')

#potential downloads:
# 1. database containing fasta and gff
# 2. independent fasta (and optional gff)

# figure out what output files are specified
fasta = snakemake.output.get("fasta")
gff = snakemake.output.get("gff") #should be optional

database = os.path.dirname(str(fasta))
zipped_db = os.path.basename(database) + ".tar.gz" # can't give path to curl

# paths to fasta within database
fasta_path = snakemake.params.get("fasta_path")
gff_path = snakemake.params.get("gff_path")

database_url = snakemake.params.get("database_url") #optional
fasta_url = snakemake.params.get("fasta_url") # optional
gff_url = snakemake.params.get("gff_url") #optional

# add check to see if either database_url or fasta_url exists (one is required!)

shell("mkdir -p {database} {log}")

if database_url:
    db_fileformat = snakemake.params.get("database_fileformat", "tar.gz")
    unzip_cmd = get_unzip_cmd(db_fileformat, "{database}", "{zipped_db}")
    db_cmd = cmd + ["{database_url}"] + unzip_cmd + ["{log}"]
    shell(' '.join(db_cmd))
    shell("tar -xzf {zipped_db} --directory {database} --strip-components 1")
    if fasta_path and not fasta_url:
        shell("ln -s {fasta_path} {fasta}")
    if gff_path and not gff_url:
        shell("ln -s {gff_path} {gff}")

if fasta_url: # independent link for fasta
    fasta_fileformat = snakemake.params.get("fasta_fileformat", "gz")
    unzip_cmd = get_unzip_cmd(fasta_fileformat, "{fasta}")
    fasta_cmd = cmd + ["{fasta_url}"] + unzip_cmd + ["{log}"]
    shell(' '.join(fasta_cmd))

if gff_url: # independent link for gff
    gff_fileformat = snakemake.params.get("gff_fileformat", "gz")
    unzip_cmd = get_unzip_cmd(gff_fileformat, "{gff}")
    gff_cmd = cmd + ["{gff_url}"] + unzip_cmd + ["{log}"]
    shell(' '.join(gff_cmd))

if not any([gff_path, gff_url]):
    shell("touch {gff} {log}")

#fileformat = snakemake.params.get("database_fileformat")

#if snakemake.params.get('md5'):
#    cmd.append('&& python -c "assert \'`md5sum {output} | '
#                'awk \'{{print $1}}\'`\' == \'{snakemake.params.md5}\', '
#                '\'MD5sum does not match\'"')


