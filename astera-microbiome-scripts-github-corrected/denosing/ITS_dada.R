# Set the R library location
.libPaths("/fs/ess/PAS0471/jelmer/R/metabar")
dyn.load("/fs/ess/PAS0471/jelmer/software/GLPK/lib/libglpk.so.40", local = FALSE)
# Load the packages
library(dada2)
library(phyloseq)
library(tidyverse)
n_cores <- 4

# Define the input files

# Dir with input FASTQ files (post-cutadapt trimmed reads):
indir <- "/fs/ess/PAS1212/Afaf_rich/ITS/results/cutadapt/ITS"



# FASTA file with training data - a database your sequences will be compared to (SILVA, UNITE, etc):
# (Check for an up-to-date version at <https://benjjneb.github.io/dada2/training.html>)
tax_key <- "/fs/ess/PAS1212/Afaf_rich/ITS/sh_general_release_dynamic_18.07.2023_dev.fasta"  

# File with sample metadata (data about your samples, treatments, blocking, etc):
metadata_file <- "/fs/ess/PAS1212/Afaf_rich/ITS/metadata_cutadapt.tsv" 

# Define the output dirs
outdir <- "/fs/ess/PAS1212/Afaf_rich/ITS/results/cutadapt/ITS/results/dada"                     # Most output
filter_dir <- "/fs/ess/PAS1212/Afaf_rich/ITS/results/dada/filtered_fastq"  # Filtered FASTQ files

# Create directories if they are not already created:
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(filter_dir, showWarnings = FALSE, recursive = TRUE)


# Create vectors with FASTQ file names
fastqs_raw_F <- sort(list.files(indir, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fastqs_raw_R <- sort(list.files(indir, pattern = "_R2_001.fastq.gz", full.names = TRUE))

head(fastqs_raw_F)
# Read and prepare the metadata
metadata_df <- read.table(file = metadata_file, sep = "\t", header = TRUE)
# Assuming metadata_df is your metadata dataframe
#metadata_df <- metadata_df[!grepl("R2", metadata_df$SampleID), ]


colnames(metadata_df)[6] <- "SampleID"
rownames(metadata_df) <- metadata_df$SampleID

head(metadata_df)


# Check the sample IDs in the metadata
metadata_df$SampleID

# Check the FASTA file names
# (Note, basename() strips the dir name from the filename)
head(basename(fastqs_raw_F))
# Extract the sample IDs from the FASTQ file names
# sub() arguments: sub(pattern, replacement, vector)
sampleIDs <- sub("_R.*", "", basename(fastqs_raw_F))


#metadata_df$SampleID <- sub("_R.*", "", basename(fastqs_raw_F))
#head(sampleIDs)
#head(metadata_df$SampleID)
metadata_df$SampleID <- sub("_L001_.*$", "", metadata_df$SampleID)
# Remove "_L001" from each element in sampleIDs
sampleIDs <- gsub("_L001", "", sampleIDs)
# Check whether the sample IDs match
identical(sort(metadata_df$SampleID), sampleIDs)
setdiff(sort(metadata_df$SampleID), sampleIDs)

setdiff(sampleIDs, sort(metadata_df$SampleID))

# Define the paths to the filtered (output) FASTQ files
# (filter_dir directory set and created in the "set file paths" step,
# still have seperate names for FWD and RVS reads at this point.)
fastqs_filt_F <- file.path(filter_dir, paste0(sampleIDs, "_F_filt.fastq"))
fastqs_filt_R <- file.path(filter_dir, paste0(sampleIDs, "_R_filt.fastq"))
plotQualityProfile(fastqs_raw_F[1:6])
plotQualityProfile(fastqs_raw_R[1:6])

# Filter the FASTQ files
n_cores <- 7
filter_results <-
  filterAndTrim(fastqs_raw_F, fastqs_filt_F,
                fastqs_raw_R, fastqs_filt_R,
                #truncLen = c(249,245),
                maxN = 0,
                maxEE = c(2,2),
                truncQ = 2,
                rm.phix = TRUE,
                multithread = n_cores,
                compress = FALSE, verbose = TRUE) 


head(filter_results)

# Dereplicate the FASTQ files
#fastqs_derep_F <- derepFastq(fastqs_filt_F, verbose = FALSE)
#fastqs_derep_R <- derepFastq(fastqs_filt_R, verbose = FALSE)

#names(fastqs_derep_F) <- sampleIDs
#names(fastqs_derep_R) <- sampleIDs

# Error training

exists <- file.exists(fastqs_filt_F)
fastqs_derep_F <- derepFastq(fastqs_filt_F[exists], verbose=TRUE)
fastqs_derep_R <- derepFastq(fastqs_filt_R[exists], verbose=TRUE)
names(fastqs_derep_F) <- sample.names[exists]
names(fastqs_derep_R) <- sample.names[exists]
err_F <- learnErrors(fastqs_derep_F, multithread=TRUE)
err_R <- learnErrors(fastqs_derep_R, multithread=TRUE)

#err_F <- learnErrors(fastqs_derep_F, multithread = n_cores, verbose = TRUE)
#err_R <- learnErrors(fastqs_derep_R, multithread = n_cores, verbose = TRUE)

# Plot the error profiles
plotErrors(err_F, nominalQ = TRUE)
plotErrors(err_R, nominalQ = TRUE)

# Infer ASVs
dada_Fs <- dada(fastqs_derep_F, err = err_F, pool = FALSE, multithread = n_cores)
dada_Rs <- dada(fastqs_derep_R, err = err_R, pool = FALSE, multithread = n_cores)

dada_Fs[[1]]

# Merge forward and reverse reads
mergers <- mergePairs(dada_Fs, fastqs_derep_F,
                      dada_Rs, fastqs_derep_R,
                      verbose = TRUE)

# Save the dereplicated objects
saveRDS(fastqs_derep_F, file = file.path(outdir, "fastqs_derep_F.rds"))
saveRDS(fastqs_derep_R, file = file.path(outdir, "fastqs_derep_R.rds"))

rm(fastqs_derep_F, fastqs_derep_R) # Remove objects from environment


# Create the sequence table
seqtab_all <- makeSequenceTable(mergers)

# Check the dimensions of the sequence table
dim(seqtab_all)

# Check the ASV sequence lengths
table(nchar(getSequences(seqtab_all)))

# Identify and remove chimeras
seqtab <- removeBimeraDenovo(seqtab_all,
                             method = "consensus",
                             multithread = n_cores,
                             verbose = TRUE)
# Proportion of retained sequences:
sum(seqtab) / sum(seqtab_all)

# Save the sequence table
saveRDS(seqtab, file = file.path(outdir, "seqtab.rds"))


# Get the numbers of unique reads at each step
getN <- function(x) sum(getUniques(x))
track <- cbind(filter_results[exists,], sapply(dada_Fs, getN), sapply(dada_Rs, getN), sapply(mergers, getN), rowSums(seqtab))
denoised_F <- sapply(dada_Fs, getN)
denoised_R <- sapply(dada_Rs, getN)
merged <- sapply(mergers, getN)

#We’ll join these vectors together with the filter_results dataframe, and the number of non-chimeric reads:

# Create a summary table
nreads_summary <- data.frame(filter_results[exists,],
                             denoised_F,
                             denoised_R,
                             merged,
                             nonchim = rowSums(seqtab),
                             row.names = sampleIDs)
colnames(nreads_summary)[1:2] <- c("input", "filtered")
head(nreads_summary)
percentage <-  round(nreads_summary[, "nonchim"] / nreads_summary[, "input"] * 100)

nreads_summary <- cbind(nreads_summary, percentage = percentage)

# Write the summary table to file
write.table(nreads_summary,
            file = file.path(outdir, "nreads_summary.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE)
# Assign taxonomy to the ASVs
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

# Create the phyloseq object
ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE),
               sample_data(metadata_df),
               tax_table(tax))

# Saves the phyloseq object as an .rds file (which can be imported directly by phyloseq):
saveRDS(ps, file = file.path(outdir, "ps.rds"))

View(sample_data(ps))

