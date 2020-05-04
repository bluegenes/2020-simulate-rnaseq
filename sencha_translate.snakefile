import os
import sys
import pandas as pd

sampleDF = pd.read_csv(config.get("samples_csv"), dtype=str, sep=",", header=0, index_col=0)

out_dir = config.get("out_dir", "sencha_out")
index_dir = os.path.join(out_dir, "index")
translate_dir = os.path.join(out_dir, "translate")

data_dir = os.path.join(out_dir, "input_data")
logs_dir = os.path.join(out_dir, "logs")
envs_dir = "envs"

reference_info = config["sencha_params"].get("peptide_references")
refnames = []
default_tablesize = "1e8"

# include tablesize in reference name
for ref, info in reference_info.items():
    tablesize = info.get("tablesize", default_tablesize)
    refnames.append(f"{ref}_t{tablesize}")

# build output files based on config
output_targets = []
ref_targets = []
output_extensions = ["codingpep.fa", "noncoding.fa", "lowcomplexnucl.fa", "csv", "json"]

alphabets_to_run = config["sencha_params"]["alphabet"] # dictionary of alphabets to run, corresponding ksizes 
# build final file targets based on config info
for alpha, alpha_info in alphabets_to_run.items():
    output_targets+=expand(os.path.join(translate_dir, "{sample}_{alphabet}_k{k}_ref{refname}.{ext}"), sample=list(sampleDF.index), alphabet=alpha, k=alpha_info["ksizes"], refname=refnames, ext=output_extensions)
    ref_targets+=expand(os.path.join(index_dir, "ref{refname}_{alphabet}_k{k}.index"), refname=refnames, alphabet=alpha, k=alpha_info["ksizes"])

rule all:
    input: 
        ref_targets,
        output_targets

# abbreviate "hydrophobic-polar" as hp in filenames. Need full alpha name for sencha code
alpha_abbreviations = {"hp": "hydrophobic-polar", "protein": "protein", "dayhoff": "dayhoff"}

rule sencha_index:
    input: lambda w: reference_info[w.ref]["path"] # get file path from reference name, config
    output: os.path.join(index_dir, "ref{ref}_t{tablesize}_{alphabet}_k{ksize}.index"),
    log: os.path.join(logs_dir, "sencha_index", "ref{ref}_t{tablesize}_{alphabet}_k{ksize}.index.log")
    benchmark: os.path.join(logs_dir, "sencha_index", "ref{ref}_t{tablesize}_{alphabet}_k{ksize}.index.benchmark")
    params:
        alphabet=lambda w: alpha_abbreviations[w.alphabet],
    threads: 1
    resources:
        mem_mb=lambda w:  reference_info[w.ref].get("index_memory", 20000),
        runtime=6000 # ~4 days
    wildcard_constraints:
        ref="\w+",
        alphabet="\w+",
        ksize="\d+",
    conda: os.path.join(envs_dir, "sencha-env.yml")
    shell:
        """
        sencha index --alphabet {params.alphabet} --peptide-ksize {wildcards.ksize} {input} --tablesize {wildcards.tablesize} --save-as {output} 2> {log}
        """

rule sencha_translate:
    input:
        fastq=lambda w: sampleDF.loc[w.sample, "read1"],
        index= rules.sencha_index.output
    output:
        coding_prot=os.path.join(translate_dir, "{sample}_{alphabet}_k{ksize}_ref{ref}_t{tablesize}.codingpep.fa"),
        coding_nucl=os.path.join(translate_dir, "{sample}_{alphabet}_k{ksize}_ref{ref}_t{tablesize}.codingnucl.fa"),
        noncoding_nucl=os.path.join(translate_dir, "{sample}_{alphabet}_k{ksize}_ref{ref}_t{tablesize}.noncoding.fa"),
        low_complexity_prot=os.path.join(translate_dir, "{sample}_{alphabet}_k{ksize}_ref{ref}_t{tablesize}.lowcomplexprot.fa"),
        low_complexity_nucl=os.path.join(translate_dir, "{sample}_{alphabet}_k{ksize}_ref{ref}_t{tablesize}.lowcomplexnucl.fa"),
        csv=os.path.join(translate_dir, "{sample}_{alphabet}_k{ksize}_ref{ref}_t{tablesize}.csv"),
        json=os.path.join(translate_dir, "{sample}_{alphabet}_k{ksize}_ref{ref}_t{tablesize}.json"),
    log: os.path.join(logs_dir, "sencha_translate", "{sample}_{alphabet}_k{ksize}_ref{ref}_t{tablesize}.translate.log") #2>{log} err ("missing PEPTIDES file")
    benchmark: os.path.join(logs_dir, "sencha_translate", "{sample}_{alphabet}_k{ksize}_ref{ref}_t{tablesize}.translate.benchmark")
    params:
        alphabet=lambda w: alpha_abbreviations[w.alphabet],
    threads: 1
    resources:
        #mem_mb=80000, #lambda wildcards, attempt: attempt *40000, #16GB*attempt
        #runtime=6000
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=60
    wildcard_constraints:
        ref="\w+",
        alphabet="\w+",
        ksize="\d+",
    conda: os.path.join(envs_dir, "sencha-env.yml")
    shell:
        """
        sencha translate --verbose --peptides-are-bloom-filter --alphabet {params.alphabet} --peptide-ksize {wildcards.ksize} --noncoding-nucleotide-fasta {output.noncoding_nucl} --low-complexity-nucleotide-fasta {output.low_complexity_nucl} --coding-nucleotide-fasta {output.coding_nucl} --low-complexity-peptide-fasta {output.low_complexity_prot} --csv {output.csv} --json-summary {output.json} {input.index} {input.fastq} > {output.coding_prot} 2> {log}
        """
