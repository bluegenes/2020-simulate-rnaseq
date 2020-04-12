# For human genome: need to keep ONLY chromosomes in the gff

# run this script from within polyester-env.yml (needs rtracklayer)
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

full_gff <- snakemake@input[["gff"]]
output_gff <- snakemake@output[[1]]

library(rtracklayer)
# read in gff (uses rtracklayer library)
h_gtf <- import(full_gff)
seqnames_to_keep = c(1:22, "M", "X", "Y") # keep only seqs assigned to chromosomes
filtered <- subset(h_gtf, seqnames  %in% seqnames_to_keep) # filter all the extra/unassigned/patches/etc entries from the gff file.
export(filtered, output_gff)
