

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
path <- "/fs/ess/PAS1212/Afaf_rich/ITS/results/cutadapt/ITS"
list.files(path)

# Create the phyloseq object
metadata_file <- "/fs/ess/PAS1212/Afaf_rich/ITS/metadata_cutadapt.tsv"
metadata_df <- read.table(file = metadata_file, sep = "\t", header = TRUE)


colnames(metadata_df)[1] <- "SampleID"
rownames(metadata_df) <- metadata_df$SampleID

head(metadata_df)
metadata_df$SampleID

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

head(basename(fnFs))
#sample.names <- sub("_R.*", "", basename(fnFs))

# Extract the sample IDs from the FASTQ file names
# sub() arguments: sub(pattern, replacement, vector)
#fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
#fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sampleIDs <- sub("_R.*", "", basename(fnFs))
#metadata_df$SampleID <- sub("_R.*", "", basename(fnFs))


head(sampleIDs)


# Remove "_R1_001.fastq.gz" from SampleID column
#metadata_df$SampleID <- sub("_R1_001.fastq.gz$", "", metadata_df$SampleID)

# Print the modified SampleID column
print(metadata_df$SampleID)
head(metadata_df$SampleID)


#metadata_df$SampleID <- unique(metadata_df$SampleID)

# Check whether the sample IDs match
identical(sort(metadata_df$SampleID), sampleIDs)
setdiff(sort(metadata_df$SampleID), sampleIDs)
setdiff(sampleIDs, sort(metadata_df$SampleID))

# Define the output dirs
outdir <- "/fs/ess/PAS1212/Afaf_rich/ITS/results/dada"                     # Most output
filter_dir <- "/fs/ess/PAS1212/Afaf_rich/ITS/results/dada/filtered_fastq"  # Filtered FASTQ files

# Create directories if they are not already created:
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
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
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = FALSE, minLen = 50,
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


dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

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
getN <- function(x) sum(getUniques(x))

denoised_F <- sapply(dadaFs, getN)
denoised_R <- sapply(dadaRs, getN)
merged <- sapply(mergers, getN)

##this worked better
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Write the summary table to file
outdir <- "/fs/ess/PAS1212/Afaf_rich/ITS/results"
write.table(track,
            file = file.path(outdir, "ASVs_summary.txt"),
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

###Alpha and beta diversity----
# And the assigned taxon names:
taxa_names(ps)[1:3]

sort(sample_sums(ps))
ps <- subset_samples(ps, sample_sums(ps) > 1000)
# The sample data (metadata)
sample_data(ps)[1:10]

# The count table
head(t(otu_table(ps)[1:10]))
# The taxonomy table
head(tax_table(ps))


fungi <- subset_taxa(ps, Kingdom == "Fungi")
fungi

indirr <- "/fs/ess/PAS1212/Afaf_rich/ITS/results"

taxa <- readRDS(file.path(indirr, "tax.rds"))


###network analysis
indir <- "/fs/ess/PAS1212/Afaf_rich/ITS/results"
ps <- readRDS(file.path(indir, "ps.rds"))


vtill_subset <- phyloseq::subset_samples(ps, management != "NT")
NT_subset  <- phyloseq::subset_samples(ps, management == "NT")

# Construct the differential network 
diff_net_const <- netConstruct(vtill_subset, LNT_subset,
                               filtTax = "highestVar",
                               filtTaxPar = list(highestVar = 75),
                               zeroMethod = "pseudo",
                               normMethod = "clr",
                               measure = "sparcc",
                               verbose = 1, 
                               thresh = 0.2)

# Analyze the network 
diff_net_const_aly <- netAnalyze(diff_net_const,
                                 clustMethod = "cluster_fast_greedy",
                                 graphlet = FALSE,
                                 verbose = 0)


# Compare the analysis of the networks 
comp_nem <- netCompare(bac22rot_net_aly,
                       permTest = FALSE, 
                       verbose = FALSE,
                       seed = 12345)


summary(comp_nem, 
        groupNames = c("LNT", "vtill"),
        showCentr = c("degree", "between", "closeness"), 
        numbNodes = 5)

plot(diff_net_const_aly, 
     sameLayout = TRUE, 
     layoutGroup = "union", 
     colorVec = col,
     borderCol = "gray40", 
     nodeSize = "mclr", 
     nodeSizeSpread = 3, 
     labelScale = TRUE,
     edgeTranspLow = 80, 
     edgeTranspHigh = 50,
     groupNames = c("CS", "CSW"), 
     showTitle = TRUE, 
     cexTitle = 1.2,
     mar = c(1,1,3,1), 
     repulsion = 0.9, 
     shortenLabels = "simple", 
     rmSingles = "inboth",
     nodeFilter = "clustMin", 
     nodeFilterPar = 10, 
     nodeTransp = 50, 
     hubTransp = 30)


# Saves the phyloseq object as an .rds file (which can be imported directly by phyloseq):
saveRDS(ps, file = file.path(outdir, "ps.rds"))


###try network














###post ASV analyses

library(vegan)
library(tidyverse)
library(phyloseq)
library(microViz)

# Set the ggplot plot theme
theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())

# Set a random seed for ordination analyses
set.seed(3131)

ps

# And the assigned taxon names:
taxa_names(ps)[1:3]
sort(sample_sums(ps))
# Remove samples with low counts:
ps <- subset_samples(ps, sample_sums(ps) > 1000)
sample_data(ps)[1:10]
head(t(otu_table(ps)[1:10]))
head(tax_table(ps))
fungi <- subset_taxa(ps, Kingdom == "Bacteria")
fungi


fungi_genera <- fungi %>%
  tax_glom(taxrank = "Genus") %>%                      # agglomerate at Genus level
  transform_sample_counts(function(x) {x / sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.03) %>%                         # Filter out low abundance taxa
  arrange(Genus)            

# Relative abundance by sample and rotation treatment
ggplot(bacteria_genera, aes(x = Sample, y = Abundance, fill = Genus)) + 
  facet_wrap(vars(Rotation), nrow = 2, scales = "free_y") +
  coord_flip() +
  geom_col() +
  scale_fill_manual(values = microViz::distinct_palette(pal = "kelly")) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  labs(y = "Relative Abundance (Genera > 2%)",
       title = "Relative Abundance of Bacteria",
       x = NULL)

# plot alpha diversity by sample
plot_richness(ps, measures = "Observed")


# plot multiple indices 
plot_richness(ps,
              x = "management",
              color = "treatment", 
              measures = c("Observed", "Shannon", "Simpson", "InvSimpson")) +
  theme(legend.position = "top")

# Make alpha diversity figure
plot_richness(ps,
              x = "manegement",
              measures = c("Observed", "Shannon"),
              color = "management",
              nrow = 1) +
  geom_boxplot() +
  geom_jitter(width = 0.05) +
  theme(legend.position = "none") +
  labs(x = "management", y = "Diversity")



# Calculate diversity metrics using estimate_richness and store the data using a vector
rich <- estimate_richness(ps, measures = c("Observed", "Shannon"))

# Compare observed diversity between our sample location
wilcox.observed <- pairwise.wilcox.test(rich$Observed, 
                                        sample_data(ps)$management, 
                                        p.adjust.method = "BH")
wilcox.observed


wilcox.shannon <- pairwise.wilcox.test(rich$Shannon, 
                                       sample_data(ps)$treatment, 
                                       p.adjust.method = "BH")
wilcox.shannon


# Non-metric MultiDimensional Scaling (NMDS)
ord.nmds.jsd_slv <- ordinate(ps, method = "NMDS", distance = "jsd")
# Create a "stressplot"
stressplot(ord.nmds.jsd_slv)

plot_ordination(ps, ord.nmds.jsd_slv,
                color = "management", shape = "treatment") +
  geom_point(size = 4) +
  coord_fixed() 
  #stat_ellipse(type = "t") 

# Subset the phyloseq object
ps_subset <- subset_samples(ps, management %in% c("LNT", "SNT", "vtill"), treatment == "MF")

plot_ordination(ps_subset, ord.nmds.jsd_slv,
                color = "management", shape = "treatment") +
  geom_point(size = 4) +
  coord_fixed() +
stat_ellipse(type = "t") 

# Redundancy Analysis
ord.rda <- ordinate(ps,
                    formula = o ~ management + treatment,
                    method = "CAP",
                    distance = "bray")
# plot the RDA
plot_ordination(ps,
                vegan::scores(ord.rda, scaling = 1),
                type = "sites",
                color = "management",
                shape = "treatment") +
  geom_point(size = 5)+
  theme_classic()

# Storing our metadata in a data frame
sampledf <- data.frame(sample_data(ps))

# Selecting our distance method (in this case the same as we used for NMDS,
# Jensen-Shannon divergence, but we could use different options here)
dist.mat <- phyloseq::distance(ps, method = "jsd")

# Selecting the number of permutations to perform
perm <- how(nperm = 999)
# Selecting blocks from the experiment
#setBlocks(perm) <- with(sampledf, Block)

# Run PERMANOVA with the selected distance method and your variables of interest
# (in this case, location and rotation)
adonis2(dist.mat ~ management * treatment,
        data = sampledf, permutations = perm)

library(DESeq2)

# Converting phyloseq object to DESeq2 format 
dds.df <- phyloseq_to_deseq2(ps, ~management + treatment)

# Identifying the groups in our experiment
dds.df$Group <- factor(paste0(dds.df$management, dds.df$treatment))  


# Setting our experimental desing in the DESeq2 object
design(dds.df) <- formula(~Group) 

# Generating our analysis
dds <- DESeq(dds.df, quiet = TRUE) 

results(dds) 
levels(dds$Group)

# Select the groups to contrast 
res <- results(dds, contrast = c("Group", "NTMF", "vtillMF")) 

# First look at the results 
summary(res)

# Number of taxa with significantly higher abundance in NT (how did we know that NT had higher taxa) :
sum(res$padj < 0.1 & res$log2FoldChange > 0, na.rm = TRUE)

# Number of taxa with significantly lower abundance in vtill:
sum(res$padj < 0.1 & res$log2FoldChange < 0, na.rm = TRUE)

res.order <- res[order(res$pvalue), ]
head(res.order)[c(1:2, 6)]

# Create a data frame from your DESeq results 
res.df <- as.data.frame(res)

# Add the taxonomy of the ASVs to the new data frame 
taxmelt <- as.data.frame(tax_table(bac22rot_phy))
res.df$Tax <- paste(rownames(res.df), taxmelt$Family, sep= " ")

# Filter to keep only log2folds > |10| and pvalues <= 0.05 
res05 <- dplyr::filter(res.df, padj<=0.05, log2FoldChange<=-10 | log2FoldChange>=10) 

# Making sure our new table has everything we need 
head(res05)


# Generating a dataframe from our table
res05.df <- as.data.frame(res05)


ggplot(data = res05.df,
       aes(x = Tax, y = log2FoldChange, fill = log2FoldChange)) +
  geom_bar(stat="identity") +
  coord_flip() +
  labs(title = "Taxa with significantly lower and higher abundance in CSW-WARS 
                when contrasted with CS-WARS") +
  theme(plot.title = element_text(hjust = 0.5))



##network
# Load the packages (package startup messages are not printed below)
library(ggnetwork)
library(phyloseq)
library(NetCoMi)


