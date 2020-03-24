# tiny sim script: given fasta, simulate some reads using polyester

# run this script from within polyester-env.yml (needs rtracklayer)
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

#library(polyester)
#library(Biostrings)
library(rtracklayer)
#library(BSgenome)

full_gff <- snakemake@input[[1]]
output_gff <- snakemake@output[[1]]

attribute <- snakemake@params[["attribute_to_filter_on"]]
print(attribute)
attribute_value <- snakemake@params[["attribute_value"]]
print(attribute_value)

# 1. Filter full gff and write temp subset gff
# read in gff (uses rtracklayer library)
gff <- import(full_gff)

## right now, only works for gene_name attr
# why does this work, but using attribute _not_ work?
mini_gff <- subset(gff, gene_name == attribute_value)

export(mini_gff, output_gff)

# add other filtering
#filtered_gff <- subset(full_gff, seqnames  %in% seqnames_to_keep)
