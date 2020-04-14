import os
import sys
import pandas as pd

sampleDF = pd.read_csv(config.get("samples_csv"), dtype=str, sep=",", header=0, index_col=0)

out_dir = config.get("out_dir", "sencha_out")
index_dir = config.get("out_dir", "sencha_index")
translate_dir = config.get("out_dir", "sencha_translate")
data_dir = os.path.join(out_dir, "input_data")
logs_dir = os.path.join(out_dir, "logs")
envs_dir = "envs"

reference_info = config["sencha_params"].get("peptide_references")

# build output files based on config
output_targets = []
output_extensions = ["codingpep.fa", "noncoding.fa", "lowcomplexnucl.fa", "csv"]

alphabets_to_run = config["sencha_params"]["alphabet"] # dictionary of alphabets to run, corresponding ksizes 
# build final file targets based on config info
for alpha, alpha_info in alphabets_to_run.items():
    output_targets+=expand(os.path.join(translate_dir, "{sample}_{alphabet}_k{k}_ref{ref}.{ext}"), sample=list(sampleDF.index), alphabet=alpha, k=alpha_info["ksizes"], ref=reference_info.keys(), ext=output_extensions)

rule all:
    input: output_targets

# abbreviate "hydrophobic-polar" as hp in filenames. Need full alpha name for khtools code
alpha_abbreviations = {"hp": "hydrophobic-polar", "protein": "protein", "dayhoff": "dayhoff"}

rule sencha_index:
    input: lambda w: reference_info[w.ref] # get file path from reference name, config
    output: os.path.join(index_dir, "ref{ref}_{alphabet}_k{ksize}.index"),
    log: os.path.join(logs_dir, "sencha_index", "ref{ref}_{alphabet}_k{ksize}.index.log")
    benchmark: os.path.join(logs_dir, "sencha_index", "ref{ref}_{alphabet}_k{ksize}.index.benchmark")
    params:
        alphabet=lambda w: alpha_abbreviations[w.alphabet],
    wildcard_constraints:
        ref="\w+",
        alphabet="\w+",
        ksize="\d+",
    conda: os.path.join(envs_dir, "khtools-env.yml")
    shell:
        """
        khtools index --alphabet {wildcards.alphabet} --peptide-ksize {wildcards.ksize} {input} > {output} 2> {log}
        """

rule sencha_translate:
    input:
        fastq=lambda w: sampleDF.loc[w.sample, "read1"],
        index=os.path.join(index_dir, "ref{ref}_{alphabet}_k{ksize}.index"),
    output:
        coding_prot=os.path.join(translate_dir, "{sample}_{alphabet}_k{ksize}_ref{ref}.codingpep.fa"),
        noncoding_nucl=os.path.join(translate_dir, "{sample}_{alphabet}_k{ksize}_ref{ref}.noncoding.fa"),
        low_complexity_nucl=os.path.join(translate_dir, "{sample}_{alphabet}_k{ksize}_ref{ref}.lowcomplexnucl.fa"),
        csv=os.path.join(translate_dir, "{sample}_{alphabet}_k{ksize}_ref{ref}.csv"),
    log: os.path.join(logs_dir, "sencha_translate", "{sample}_{alphabet}_k{ksize}_ref{ref}.translate.log")
    benchmark: os.path.join(logs_dir, "sencha_translate", "{sample}_{alphabet}_k{ksize}_ref{ref}.translate.benchmark")
    params:
        alphabet=lambda w: alpha_abbreviations[w.alphabet],
    wildcard_constraints:
        ref="\w+",
        alphabet="\w+",
        ksize="\d+",
    conda: os.path.join(envs_dir, "khtools-env.yml")
    shell:
        """
        khtools extract-coding --verbose --alphabet {params.alphabet} --peptide-ksize {wildcards.ksize} --noncoding-nucleotide-fasta {output.noncoding_nucl} 
        --low-complexity-nucleotide-fasta {output.low_complexity_nucl} --csv {output.csv} {input.index} {input.fastq} > {output.coding_prot} 2> {log}
        """
