
# Set the R library location
.libPaths("/fs/ess/PAS0471/jelmer/R/metabar")

.libPaths("/fs/ess/PAS1212/abdelrahim")
source(file.path(Sys.getenv("LMOD_PKG"), "init/R"))
module("load",'gnu/11.2.0')
module('load','geos/3.9.1')
module('load','proj/8.1.0')
module('load','gdal/3.3.1')

dyn.load("/fs/ess/PAS0471/jelmer/software/GLPK/lib/libglpk.so.40", local = FALSE)
# Load the packages
install.packages("BiocManager")
library(BiocManager)
BiocManager::install('dada2')
install.packages("dada2")
library(dada2)
install.packages()
library(phyloseq)
library(tidyverse)
library(vegan)
library(microViz)
library(ggplot2)
install.packages("ggpubr")
library(ggpubr)
path <- "/fs/ess/PAS1212/Afaf_rich/AMF_CS_OH/231117_Rich_demultiplxed_fastqs"
list.files(path)


fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_L001_"), `[`, 1)


plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names



out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,255),
                     maxN=0, maxEE=c(2,2), trimLeft = 10, truncQ=2, rm.phix=FALSE,
                     compress=FALSE, multithread=n_cores, verbose = TRUE) # On Windows set multithread=FALSE
head(out)
file_path <- file.path(path, "filter_reads_summary.csv")
write.csv(out, file_path, row.names=FALSE)

exists <- file.exists(filtFs)
derepFs <- derepFastq(filtFs[exists], verbose=TRUE)
derepRs <- derepFastq(filtRs[exists], verbose=TRUE)
names(derepFs) <- sample.names[exists]
names(derepRs) <- sample.names[exists]
errF <- learnErrors(derepFs, multithread=TRUE)
errR <- learnErrors(derepRs, multithread=TRUE)

dadaFs <- dada(derepFs, err=errF, pool = FALSE, multithread=n_cores)

dadaRs <- dada(derepRs, err=errR, pool = FALSE, multithread=n_cores)

dadaFs[[1]]

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


seqtab_all <- makeSequenceTable(mergers)
dim(seqtab_all)

table(nchar(getSequences(seqtab_all)))

seqtab.nochim <- removeBimeraDenovo(seqtab_all, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab_all)


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


file_path <- file.path(path, "densoised_reads_amf.csv")
write.csv(track, file_path, row.names=FALSE)


###Assign taxonomy
tax_key <- "/fs/ess/PAS1212/Afaf_rich/AMF_CS_OH/231117_Rich_demultiplxed_fastqs/maarjam_seqsdb.fasta"
tax <- assignTaxonomy(seqtab_all, tax_key, multi = TRUE, verbose = TRUE)
outdir <- "/fs/ess/PAS1212/Afaf_rich/AMF_CS_OH/231117_Rich_demultiplxed_fastqs"
saveRDS(tax, file.path(outdir, "tax.rds"))

##inspect tax assignment
taxa.print <- tax # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

