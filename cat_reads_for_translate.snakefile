import os
import glob

readsdir = config.get("readsdir", "/home/ntpierce/2020-simulate-rnaseq/QfO_vertebrates/simulate_reads")

combine_reftypes = config.get("reftypes", ["ncRNA","qfo_dna"])

species_names = config.get("species_names", ["Homo_sapiens"])
reftype_str = "_".join(combine_reftypes)
#run_name = config.get("run_name", "combined_reads") + "_" + combined_reftypes
outdir = config.get("outdir", os.path.join(readsdir, "combined_reads")) 
NUMREADS = config.get("numreads_in_subsample", 100)

def generate_replist(number_replicates):
    replist = list(range(number_replicates))
    reps=[]
    for num in replist:
        num+=1 # start at 1
        if num < 10:
            num = "0" + str(num)
        reps.append(str(num))
    return reps

replist = generate_replist(config.get("num_replicates", 5))

#Homo_sapiens_qfo_dna_01,/home/ntpierce/2020-simulate-rnaseq/QfO_vertebrates/simulate_reads/Homo_sapiens_qfo_dna/Homo_sapiens_qfo_dna_s100_01.fq.gz


rule all:
    input: expand(os.path.join(outdir, "{species}.{reftypes}.s{numreads}.samples.csv"), species=species_names, reftypes=reftype_str, numreads=NUMREADS)

def aggregate_reads(w):
    reads = []
    for reftype in combine_reftypes:
        reads += expand(os.path.join(readsdir, "{{species}}_{reftype}/{{species}}_{reftype}_s{{numreads}}_{{replicate}}.fq.gz"), reftype=reftype)
    return reads

rule cat_reads:
    input: aggregate_reads
    # reftypes wildcard was not solving properly (hard bc all "_") - so just specify it here
    output: expand(os.path.join(outdir, "{{species}}_{reftypes}_s{{numreads}}_{{replicate}}.fq.gz"), reftypes=reftype_str)
    wildcard_constraints:
        species="\w+",
        reftypes="\w+",
        replicate="\d+",
        numreads="\d+"
    shell:
        """
        cat {input} > {output}
        """
    
rule write_combined_samples_csv:
    input: expand(os.path.join(outdir, "{{species}}_{{reftypes}}_s{{numreads}}_{replicate}.fq.gz"), replicate=replist)
    output: os.path.join(outdir, "{species}.{reftypes}.s{numreads}.samples.csv")
    threads: 1
    resources:
      mem_mb=1000,
      runtime=15
    wildcard_constraints:
        species="\w+",
        reftypes="\w+",
        replicate="\d+",
        numreads="\d+"
    run:
        with open(str(output), "w") as csv:
            csv.write("sample_id" + "," + "read1"+ "\n")
            for f in input:
                #refname = os.path.basename(os.path.dirname(f))
                refname = f"{wildcards.species}_{wildcards.reftypes}"
                replicate = f.rsplit(".fq.gz")[0].rsplit("_", 1)[1] # _01.fq.gz
                csv.write(f"{refname}_{replicate}" + "," + os.path.abspath(f) + "\n")
