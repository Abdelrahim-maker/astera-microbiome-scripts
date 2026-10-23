

##load libraries
# Set the R library location
.libPaths("/fs/ess/PAS0471/jelmer/R/metabar")
dyn.load("/fs/ess/PAS0471/jelmer/software/GLPK/lib/libglpk.so.40", local = FALSE)
install.packages("RCPP")
library(Rcpp)
install.packages("dada2")
library(dada2)
library(phyloseq)
library(tidyverse)
install.packages("ShortRead")
library(ShortRead)
path <- "/fs/ess/PAS1212/Afaf_rich/16srRNA_2"
list.files(path)

# Create the phyloseq object
metadata_file <- "/fs/ess/PAS1212/Afaf_rich/16srRNA_2/metadata_cutadapt.tsv"
metadata_df <- read.table(file = metadata_file, sep = "\t", header = TRUE)

colnames(metadata_df)[1] <- "SampleID"
rownames(metadata_df) <- metadata_df$SampleID

head(metadata_df)
metadata_df$SampleID
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

head(basename(fnFs))

# Extract the sample IDs from the FASTQ file names
# sub() arguments: sub(pattern, replacement, vector)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sampleIDs <- sub("_R.*", "", basename(fnFs))
head(sampleIDs)

# Check whether the sample IDs match
identical(sort(metadata_df$SampleID), sampleIDs)
setdiff(sort(metadata_df$SampleID), sampleIDs)
setdiff(sampleIDs, sort(metadata_df$SampleID))

# Define the output dirs
outdir <- "/fs/ess/PAS1212/Afaf_rich/16srRNA_2/results/dada"                     # Most output
filter_dir <- "/fs/ess/PAS1212/Afaf_rich/16srRNA_2/results/dada/filtered_fastq"  # Filtered FASTQ files

# Create directories if they are not already created:
dir.create(filter_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(filter_dir, showWarnings = FALSE, recursive = TRUE)




# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
#fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
#fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
##sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Extract sample names
sample.names <- sub("_R.*", "", basename(fnFs))
#sample.names <- sapply(strsplit(basename(fnFs), "_"), function(x) paste(x[1:3], collapse = "_"))

# View the extracted sample names
print(sample.names)

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered_2", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered_2", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Specify the output directory
n_cores <- 4

# Now, run filterAndTrim() ##maxEE was relaxed for the reverse reads to allow all reads to pass the filter
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(148, 148),
                     maxN = 0, maxEE = c(2, 5), truncQ = 2, rm.phix = FALSE,
                     compress = FALSE, multithread = n_cores, verbose = TRUE)
head(out)

#this code helped with the error that not all files exist
exists <- file.exists(filtFs)
derepFs <- derepFastq(filtFs[exists], verbose=TRUE)
derepRs <- derepFastq(filtRs[exists], verbose=TRUE)
names(derepFs) <- sample.names[exists]
names(derepRs) <- sample.names[exists]
errF <- learnErrors(derepFs, multithread=TRUE)
errR <- learnErrors(derepRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)


dadaFs <- dada(derepFs, err=errF, pool = FALSE, multithread=n_cores)

dadaRs <- dada(derepRs, err=errR, pool = FALSE, multithread=n_cores)

dadaFs[[1]]


mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

##remove dereplicated objects
# Save the dereplicated objects
saveRDS(derepFs, file = file.path(outdir, "fastqs_derep_F.rds"))
saveRDS(derepRs, file = file.path(outdir, "fastqs_derep_R.rds"))

rm(derepFs, derepRs) # Remove objects from environment

##make the sequence table
seqtab_all <- makeSequenceTable(mergers)
dim(seqtab_all)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_all)))

seqtab <- removeBimeraDenovo(seqtab_all, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab)

sum(seqtab)/sum(seqtab_all)

# Save the sequence table
saveRDS(seqtab, file = file.path(outdir, "seqtab.rds"))

getN <- function(x) sum(getUniques(x))

##try this
# Get the numbers of unique reads at each step
#getN <- function(x) sum(getUniques(x))

denoised_F <- sapply(dadaFs, getN)
denoised_R <- sapply(dadaRs, getN)
merged <- sapply(mergers, getN)

##this worked better
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
#colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
#rownames(track) <- sample.names
#head(track)

# Write the summary table to file
#outdir <- "/fs/ess/PAS1212/Afaf_rich/16srRNA_2/results"
#write.table(track,
            #file = file.path(outdir, "ASVs_summary.txt"),
            #sep = "\t", quote = FALSE, row.names = TRUE)

# Create a summary table
nreads_summary <- data.frame(out,
                             denoised_F,
                             denoised_R,
                             merged,
                             nonchim = rowSums(seqtab),
                             row.names = sampleIDs)
colnames(nreads_summary)[1:2] <- c("input", "filtered")
# Have a look at the first few rows of the summary table
head(nreads_summary)
# Write the summary table to file
write.table(nreads_summary,
            file = file.path(outdir, "nreads_summary.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE)


##taxonomy

# Assign taxonomy to the ASVs
tax_key <- "data/ref/silva_nr99_v138.1_train_set.fa.gz"
tax <- assignTaxonomy(seqtab, tax_key, multi = TRUE, verbose = TRUE)
saveRDS(tax, file.path(outdir, "tax.rds"))

if(!identical(getSequences(tax), getSequences(seqtab))) stop("Taxonomy mismatch.")

# How many ASVs are in each phylum?
table(tax[,"Phylum"], useNA = "ifany")

# Write the output files

# Prepare sequences and headers:
asv_seqs <- colnames(seqtab)
asv_headers <- paste(">ASV", 1:ncol(seqtab), sep = "_")

# Interleave headers and sequences:
asv_fasta <- c(rbind(asv_headers, asv_seqs))

# Write fasta file:
write(asv_fasta, file = file.path(outdir, "ASVs.fa"))

# Fix the sample IDs in the sequence table
rownames(seqtab) <- sub("_R.*", "", rownames(seqtab))



ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE),
               sample_data(metadata_df),
               tax_table(tax))
# Saves the phyloseq object as an .rds file (which can be imported directly by phyloseq):
saveRDS(ps, file = file.path(outdir, "ps.rds"))

