# tiny sim script: given fasta, simulate some FASTQ using rsubread

# run this script from within polyester-env.yml
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(Rsubread)
library(Biostrings)
#library(rtracklayer)
#library(BSgenome)

subset_fasta <- snakemake@input[[1]]
output_dir <- snakemake@params[["output_dir"]]
output_prefix <- snakemake@params[["output_prefix"]]
read_length <- as.integer(snakemake@params[["read_length"]]) # 150
num_reps <- as.integer(snakemake@params[["num_reps"]]) # 5
paired <- snakemake@params[["simulate_paired"]] # true/false
reads_per_tx <- as.integer(snakemake@params[["num_reads_per_transcript"]])

## Modified from Rsubread example

# Scan through the fasta file to get transcript names and lengths
transcripts <- scanFasta(subset_fasta)   #"GENCODE-Human-transcripts.fa.gz")
head(transcripts)
#nsequences <- nrow(transcripts) #- sum(transcripts$Duplicate)

nsequences <- length(readDNAStringSet(subset_fasta))


# Assign reads_per_tx transcript count per sequence. Or counts are relative - set totals via library size?
#counts <- rep(reads_per_tx, nrow(transcripts))
counts <- rep(0, nrow(transcripts))
counts[!transcripts$Duplicate] <- 1
libSize = reads_per_tx * nrow(transcripts)

# Assign a random TPM value to each non-duplicated transcript sequence
#TPMs <- rep(0, nrow(transcripts))
#TPMs[!transcripts$Duplicate] <- rexp(nsequences)

# Generate actual reads.
# The output read file is my-simulated-sample_R1.fastq.gz 
# The true read counts are returned.
#true.counts <- simReads("GENCODE-Human-transcripts.fa.gz", TPMs, "my-simulated-sample")
true.counts <- simReads(subset_fasta, counts, output_prefix, library.size=libSize, read.length=read_length, simulate.sequencing.error=FALSE) #, simulate.sequencing.error=TRUE)
print(true.counts[1:10,])

# enable either coverage or solid number
#coverage <- snakemake@params[["coverage"]]
#coverage = round(20 * width(fastaStrSet) / 100) # 20x coverage

#countmat = matrix(reads_per_tx, nrow=length(readDNAStringSet(subset_fasta)), ncol=num_reps)
#print(countmat)
# simulate reads
#simulate_experiment_countmat(fasta=subset_fasta, readmat=countmat, reportCoverage=TRUE, readlen=read_length, paired=paired, outdir=output_dir, gzip=TRUE)
