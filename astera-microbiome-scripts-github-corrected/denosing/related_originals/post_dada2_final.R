# Set the R library location
.libPaths("/fs/ess/PAS0471/jelmer/R/metabar")

dyn.load("/fs/ess/PAS0471/jelmer/software/GLPK/lib/libglpk.so.40", local = FALSE)
# Load the packages
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
install.packages("ggsignif")
library(ggsignif)


install_github("umerijaz/microbiomeSeq") # Install the package
library(microbiomeseq)
install.packages("microbiomeSeq")
# Define a temporary library location
temp_lib <- tempdir()
# Install the package from GitHub to the temporary library
remotes::install_github("umerijaz/microbiomeSeq", lib = temp_lib)
library("microbiomeSeq", lib.loc = temp_lib)
n_cores <- 4

# Define the input files

# Dir with input FASTQ files (post-cutadapt trimmed reads):
indir <- "/fs/ess/PAS1212/Afaf_rich/ITS/results/cutadapt/ITS/results/dada"
# Set the ggplot plot theme
theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())

# Set a random seed for ordination analyses
set.seed(3131)

ps <- readRDS(file.path(indir, "ps.rds"))
sort(sample_sums(ps))
# Remove samples with low counts:
ps <- subset_samples(ps, sample_sums(ps) > 1000) #this removed 2 samples and we have 7789 taxa

##prepare sequence headers (ASV1, ASV2,...)
dna <- Biostrings::DNAStringSet(taxa_names(ps))

names(dna) <- taxa_names(ps)

ps <- merge_phyloseq(ps, dna)

taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

ps


##modify the incorrect vv in the management
ps@sam_data$management <- gsub("vv", "vtill", ps@sam_data$management)
# Modify the values in the "management" column
ps@sam_data$management <- gsub("LNT", "NT", ps@sam_data$management)

ps@sam_data$management <- gsub("SNT", "NT", ps@sam_data$management)
#ps@sam_data$crop <- gsub("soybean_winter_wheat", "winterwheat", ps@sam_data$crop)



###relative abundance of species
fungi <- subset_taxa(ps, Kingdom == "k__Fungi")
fungi

fungal_species <- fungi %>%
  tax_glom(taxrank = "Species") %>%                      # agglomerate at Genus level
  transform_sample_counts(function(x) {x / sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.03) %>%                         # Filter out low abundance taxa
  arrange(Species) 

# Relative abundance by sample and rotation treatment
ggplot(fungal_species, aes(x = Sample, y = Abundance, fill = Species)) + 
  facet_wrap(vars(management), nrow = 2, scales = "free_y") +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +  # Specify width of bars
  coord_flip() +
  scale_fill_discrete() +
  theme(legend.key.width = unit(1, "cm"),  # Adjust legend key width
        legend.key.height = unit(1, "cm"),
        legend.text = element_text(size = 15),# Adjust legend key height
        axis.text.y = element_blank(),      # Remove y-axis labels
        axis.text.x = element_text(size = 25),  # Set x-axis label size
        axis.title = element_text(size = 25),  # Increase axis label size
        strip.text = element_text(size = 25)) +  # Increase facet wrap label size
  labs(y = "Relative Abundance (Species > 2%)",  # Removed y-axis label
       title = "Relative Abundance of Fungal Species")  


###plot alpha diversity
# plot multiple indices 
plot_richness(ps,
              x = "management",
              color = "crop", 
              measures = c("Observed", "Shannon", "Simpson", "InvSimpson")) +
  theme(legend.position = "top")

# Calculate diversity metrics using estimate_richness and store the data using a vector
rich <- estimate_richness(ps, measures = c("Observed", "Shannon"))

# Compare observed diversity between our sample location
wilcox.observed <- pairwise.wilcox.test(rich$Observed, 
                                        sample_data(ps)$treatment, 
                                        p.adjust.method = "BH")
wilcox.observed

wilcox.shannon <- pairwise.wilcox.test(rich$Shannon, 
                                       sample_data(ps)$treatment, 
                                       p.adjust.method = "BH")
wilcox.shannon

# Extract p-values
p_values_shannon <- wilcox.shannon$p.value



# Plot richness
treatment_colors <- c("MF" = "blue",
                       "MS" = "light brown",
                       "SD" = "orange")
plot_richness(ps,
              x = "treatment",
              measures = c("Observed", "Shannon"),
              nrow = 1) +
  scale_fill_manual(values = treatment_colors) +
  geom_boxplot() +
  geom_jitter(width = 0.05) +
  theme(legend.position = "none") +
  labs(x = "Location", y = "Diversity")

# Add horizontal error bars with stars indicating significance
# Add significance stars
p_with_errorbars <- p +
  geom_signif(comparisons = list(c("Group1", "Group2")),
              map_signif_level = TRUE, textsize = 6)

p_with_errorbars


##########################
# 1. Calculate alpha diversity measures
alpha_div <- estimate_richness(ps, measures = c("Shannon", "Observed"))

# 2. Prepare the data for plotting
alpha_div_df <- as.data.frame(alpha_div) %>%
  mutate(Management = sample_data(ps)$management)  # Assuming 'management' is the grouping variable

# 3. Create boxplots with horizontal error bars and significance letters
p <- ggplot(alpha_div_df, aes(x = Management, y = Shannon)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +  # Add jitter for better visualization
  stat_summary(fun.data = mean_cl_normal, geom = "errorbarh", width = 0.1) +  # Add horizontal error bars
  geom_signif(comparisons = list(c("Group1", "Group2")), 
              map_signif_level = TRUE, textsize = 4, vjust = -0.5) +  # Add letters for significance
  labs(x = "Management", y = "Shannon Diversity") +
  theme_minimal()

# Print the plot
print(p)


#############################################################################################
##beta diversity

##subset only MF all depth
ps_mf <- subset_samples(ps, treatment %in% "MF")
ps_mf_nt_vtill <- subset_samples(ps, treatment %in% "MF"& management %in% c("NT", "vtill"))
ord.nmds.jsd_slv <- ordinate(ps_mf_nt_vtill, method = "PCoA", distance = "jsd")
stressplot(ord.nmds.jsd_slv)
# Define colors for each management type
management_colors <- c("forest" = "darkgreen",
                       "meadow" = "brown",
                       "NT" = "#6495ED",
                       "vtill" = "red")

plot_ordination(ps, ord.nmds.jsd_slv,
                color = "management") +
  geom_point(size = 5) +  # Black outline for points
  scale_color_manual(values = management_colors) +  # Adjust colors
  scale_fill_manual(values = management_colors) +  # Transparent fill for points
  coord_fixed() +
  stat_ellipse(level = 0.95, geom = "polygon", alpha = 0.2, fill = "#CCCCCC") +  # Grey fill for ellipses
  theme(legend.position = "bottom",  # Move legend to bottom
        plot.background = element_rect(fill = "transparent"),  # Set background to transparent
        axis.title = element_text(size = 25),  # Adjust axis label size
        axis.text = element_text(size = 25),   # Adjust axis text size
        legend.text = element_text(size = 25),  # Adjust legend text size
        legend.title = element_text(size = 25))


################################################################################3
###permanova

# Storing our metadata in a data frame
sampledf <- data.frame(sample_data(ps_mf_nt_vtill))

# Selecting our distance method (in this case the same as we used for NMDS,
# Jensen-Shannon divergence, but we could use different options here)
dist.mat <- phyloseq::distance(ps_mf_nt_vtill, method = "jsd")

# Selecting the number of permutations to perform
perm <- how(nperm = 999)
# Selecting blocks from the experiment
setBlocks(perm) <- with(sampledf, Block)

# Run PERMANOVA with the selected distance method and your variables of interest
# (in this case, location and rotation)
adonis2(dist.mat ~ management/Sample_UID_code + depth,
        data = sampledf, permutations = perm)

###################################################################################3
library(DESeq2)

# Get group names to parse onto network construction 
ps_m <- subset_samples(ps, management %in% c("NT", "vtill"))
# Collapse taxa to species level
physeq_species <- tax_glom(ps, taxrank = "Species")
# Converting phyloseq object to DESeq2 format 
dds.df <- phyloseq_to_deseq2(physeq_species, ~management)
# Identifying the groups in our experiment
dds.df$Group <- factor(paste0(dds.df$management))  

# Setting our experimental desing in the DESeq2 object
design(dds.df) <- formula(~Group) 

# Generating our analysis
dds <- DESeq(dds.df, quiet = TRUE) 

results(dds) 
levels(dds$Group)
# Select the groups to contrast 
res <- results(dds, contrast = c("Group", "NT", "vtill")) 

# First look at the results 
summary(res) 
# Number of taxa with significantly higher abundance in NT :
sum(res$padj < 0.1 & res$log2FoldChange > 0, na.rm = TRUE)
# Number of taxa with significantly lower abundance in NT:
res.order <- res[order(res$pvalue), ]
head(res.order)[c(1:2, 6)]
# Create a data frame from your DESeq results 
res.df <- as.data.frame(res)

# Add the taxonomy of the ASVs to the new data frame 
taxmelt <- as.data.frame(tax_table(physeq_species))
res.df$Taz <- paste(rownames(res.df), taxmelt$Species, sep= " ")

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
  labs(title = "Taxa with significantly lower and higher abundance in NT 
                when contrasted with vtill") +
  theme(plot.title = element_text(hjust = 0.5))
#########################################################################################

###core microbiome
library(microbiome)
library(ggvenn)

ps.dataset.cm <- microbiome::transform(ps, "compositional")
ps.CS <- subset_samples(ps.dataset.cm, management == "NT")
ps.CSW <- subset_samples(ps.dataset.cm, management == "vtill")
core.taxa.CS <- core_members(ps.CS, detection = 0.001, prevalence = 95/100)
core.taxa.CSW <- core_members(ps.CSW, detection = 0.001, prevalence = 95/100)


cores.list <- list(NT = core.taxa.CS, vtill = core.taxa.CSW)
cores.list #To see the list of core taxa (ASVs) by Rotation


tax.table.ps <- as.data.frame(tax_table(ps))

shared.taxa <-  core.taxa.CS[core.taxa.CS %in% core.taxa.CSW]
tax <- tax.table.ps[which(row.names(tax.table.ps) %in% shared.taxa ),]

unique.taxa.CS <-  core.taxa.CS[!core.taxa.CS %in% core.taxa.CSW]
tax.CS <- tax.table.ps[which(row.names(tax.table.ps) %in% unique.taxa.CS),]

unique.taxa.CSW <-  core.taxa.CSW[!core.taxa.CSW %in% core.taxa.CS]
tax.CSW <- tax.table.ps[which(row.names(tax.table.ps) %in% unique.taxa.CSW),]

plot <- ggvenn(cores.list, fill_color = c("#0073C2FF", "red"),
               stroke_size = 0.5, set_name_size = 5, text_size = 5,
               show_percentage = T)
print(plot)

###approach 2 for microbiome composition
df <- psmelt(ps)


df.RA <- aggregate(df$Abundance, by = list(Genus = df$Genus), FUN=sum) 
df.RA <- dplyr::rename(df.RA, TotalRA = `x`)
quantile.95.abund <- quantile(df.RA$TotalRA, 0.95)

taxa.95 <- subset(df.RA, TotalRA >= quantile.95.abund) 
taxa.95.genus <- taxa.95$Genus

ps.genus <- tax_glom(ps, "Genus")
ASVs.genus <- otu_table(ps.genus)
df.ASVs.genus <- as.data.frame(ASVs.genus)
rownames(df.ASVs.genus) <- as.data.frame(tax_table(ps.genus))$Genus


#############################################################33
##network
library(ggnetwork)
library(NetCoMi)
set.seed(1234)

# Get group names to parse onto network construction 
ps_m <- subset_samples(ps, management %in% c("NT", "vtill") & treatment %in% "MF")
# Collapse taxa to species level
physeq_species <- tax_glom(ps_m, taxrank = "Species")

# Subset to keep only species-level taxa
physeq_species <- subset_taxa(physeq_species, rank = "Species")

group_vec <- get_variable(ps_m, "management")

# Construct the network 
bac22rot_net <- netConstruct(ps_m,
                             #group = group_vec,
                             filtTax = "highestFreq",
                             filtTaxPar = list(highestFreq = 100),
                             filtSamp = "totalReads",
                             filtSampPar = list(totalReads = 500),
                             zeroMethod = "pseudo",
                             normMethod = "clr",
                             measure = "sparcc",
                             sparsMethod = "none",
                             dissFunc = "signed",
                             verbose = 0,
                             seed = 123456)



bac22rot_net_aly <- netAnalyze(bac22rot_net,
                               clustMethod = "cluster_fast_greedy",
                               graphlet = FALSE,
                               verbose = 0)
summary(bac22rot_net_aly)

# Get number of modules to select colors 
bac22rot_nclust <- as.numeric(max(names(table(bac22rot_net_aly$clustering$clust1))))

# Set color palette 
col <- c(topo.colors(bac22rot_nclust), rainbow(6))

# Plot network
plot(bac22rot_net_aly, 
     colorVec = col,
     borderCol = "gray40", 
     nodeSize = "mclr", 
     nodeSizeSpread = 3, 
     labelScale = TRUE,
     edgeTranspLow = 80, 
     edgeTranspHigh = 50,
     groupNames = c("NT", "vtill"), 
     showTitle = TRUE, 
     cexTitle = 2,
     #cexHubs = 2,
     cexHubLabels = 4,
     cexLabels = 3,
     mar = c(1,1,3,1), 
     repulsion = 0.9, 
     shortenLabels = "simple", 
     rmSingles = "inboth",
     nodeFilter = "clustMin", 
     nodeFilterPar = 10, 
     nodeTransp = 50, 
     hubTransp = 30)



# Subsetting our data into their respective Rotations 
NT_subset <- phyloseq::subset_samples(ps_m, management != "vtill")
vtill_subset  <- phyloseq::subset_samples(ps_m, management == "vtill")

# Construct the differential network 
diff_net_const <- netConstruct(NT_subset, vtill_subset,
                               filtTax = "highestVar",
                               filtTaxPar = list(highestVar = 75),
                               zeroMethod = "pseudo",
                               normMethod = "clr",
                               measure = "sparcc",
                               verbose = 1, 
                               thresh = 0.3)
diff_net_const_aly <- netAnalyze(diff_net_const,
                                 clustMethod = "cluster_fast_greedy",
                                 graphlet = FALSE,
                                 verbose = 0)

comp_nem <- netCompare(bac22rot_net_aly,
                       permTest = FALSE, 
                       verbose = FALSE,
                       seed = 12345)
summary(comp_nem, 
        groupNames = c("NT", "vtill"),
        showCentr = c("degree", "between", "closeness"), 
        numbNodes = 10)



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
     groupNames = c("NT", "vtill"), 
     showTitle = TRUE, 
     cexTitle = 1.2,
     #cexhubLabels = 6,
     cexLabels = 6,
     mar = c(1,1,3,1), 
     repulsion = 0.9, 
     shortenLabels = "simple", 
     rmSingles = "inboth",
     nodeFilter = "clustMin", 
     nodeFilterPar = 10, 
     nodeTransp = 50, 
     hubTransp = 30)
#############################################################################3
##
#build network for NT only 

ps_nt <- subset_samples(ps, management %in% "NT" & treatment %in% "MF")

bac22rot_net <- netConstruct(ps_nt,
                             #group = group_vec,
                             filtTax = "highestFreq",
                             filtTaxPar = list(highestFreq = 100),
                             filtSamp = "totalReads",
                             filtSampPar = list(totalReads = 500),
                             zeroMethod = "pseudo",
                             normMethod = "clr",
                             measure = "sparcc",
                             sparsMethod = "none",
                             dissFunc = "signed",
                             verbose = 0,
                             seed = 123456)



bac22rot_net_aly <- netAnalyze(bac22rot_net,
                               clustMethod = "cluster_fast_greedy",
                               graphlet = FALSE,
                               verbose = 0)
summary(bac22rot_net_aly)

# Get number of modules to select colors 
bac22rot_nclust <- as.numeric(max(names(table(bac22rot_net_aly$clustering$clust1))))

# Set color palette 
library(RColorBrewer)

# Number of clusters - replace this with the actual number of clusters if bac22rot_nclust is not numeric
n_clusters <- bac22rot_nclust 

# Create a color-blind friendly palette for nodes based on the number of clusters
node_colors <- brewer.pal(min(n_clusters, 8), "Pastel2") 
if (n_clusters > 8) {
  node_colors <- colorRampPalette(brewer.pal(8, "Pastel2"))(n_clusters)
}

# Now set the positive and negative colors for edges
color_blind_friendly_pos <- "#2c7bb6"  # Dark blue for positive edges
color_blind_friendly_neg <- "#654321""  # Dark orange for negative edges

# Plot the network
plot(bac22rot_net_aly, 
     colorVec = node_colors,  # Use the color-blind friendly node colors
     borderCol = "gray40", 
     nodeSize = "mclr", 
     nodeSizeSpread = 3, 
     labelScale = TRUE,
     edgeTranspLow = 80, 
     edgeTranspHigh = 50,
     groupNames = c("NT", "vtill"), 
     showTitle = TRUE, 
     cexTitle = 2,
     cexHubLabels = 4,
     highlightHubs = TRUE,
     negDiffCol = TRUE,
     posCol = color_blind_friendly_pos,  # Set color for positive edges
     negCol = color_blind_friendly_neg,  # Set color for negative edges
     cexLabels = 3,
     mar = c(1,1,3,1), 
     repulsion = 0.9, 
     shortenLabels = "simple", 
     rmSingles = "inboth",
     nodeFilter = "clustMin", 
     nodeFilterPar = 10, 
     nodeTransp = 50, 
     hubTransp = 30)

############################################################################################3
##
library(igraph)
#build network for vtill

ps_vtill <- subset_samples(ps, management %in% "vtill" & treatment %in% "MF")
bac22rot_net <- netConstruct(ps_vtill,
                             #group = group_vec,
                             filtTax = "highestFreq",
                             filtTaxPar = list(highestFreq = 100),
                             filtSamp = "totalReads",
                             filtSampPar = list(totalReads = 500),
                             zeroMethod = "pseudo",
                             normMethod = "clr",
                             measure = "sparcc",
                             sparsMethod = "none",
                             dissFunc = "signed",
                             verbose = 0,
                             seed = 123456)



bac22rot_net_aly <- netAnalyze(bac22rot_net,
                               clustMethod = "cluster_fast_greedy",
                               graphlet = FALSE,
                               verbose = 0)
summary(bac22rot_net_aly)

# Get number of modules to select colors 
bac22rot_nclust <- as.numeric(max(names(table(bac22rot_net_aly$clustering$clust1))))

# Set color palette 
# Number of clusters - replace this with the actual number of clusters if bac22rot_nclust is not numeric
n_clusters <- bac22rot_nclust 

# Create a color-blind friendly palette for nodes based on the number of clusters
node_colors <- brewer.pal(min(n_clusters, 8), "Set3") 
if (n_clusters > 8) {
  node_colors <- colorRampPalette(brewer.pal(8, "Set3"))(n_clusters)
}

# Now set the positive and negative colors for edges
color_blind_friendly_pos <- "#2c7bb6"  # Dark blue for positive edges
color_blind_friendly_neg <- "#654321"  # Dark orange for negative edges

# Plot the network
plot(bac22rot_net_aly, 
     colorVec = node_colors,  # Use the color-blind friendly node colors
     borderCol = "gray40", 
     nodeSize = "mclr", 
     nodeSizeSpread = 3, 
     labelScale = TRUE,
     edgeTranspLow = 80, 
     edgeTranspHigh = 50,
     groupNames = c("NT", "vtill"), 
     showTitle = TRUE, 
     cexTitle = 2,
     cexHubLabels = 4,
     highlightHubs = TRUE,
     negDiffCol = TRUE,
     posCol = color_blind_friendly_pos,  # Set color for positive edges
     negCol = color_blind_friendly_neg,  # Set color for negative edges
     cexLabels = 3,
     mar = c(1,1,3,1), 
     repulsion = 0.9, 
     shortenLabels = "simple", 
     rmSingles = "inboth",
     nodeFilter = "clustMin", 
     nodeFilterPar = 10, 
     nodeTransp = 50, 
     hubTransp = 30)




#########################core microbiome

# Transform to compositional abundances
library(microbiome)
# Transform the data
pseq.rel <- microbiome::transform(ps_mf_nt_vtill, "compositional")

# Define the core microbiome
pseq.core <- core(pseq.rel, detection = 0.1/100, prevalence = 50/100)

# Prevalence and detection thresholds
prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(5e-3), log10(.2), length = 10), 3)

# Create the heatmap plot
p <- plot_core(pseq.core, plot.type = "heatmap",
               prevalences = prevalences, detections = detections, min.prevalence = 0.5) +
  xlab("Detection Threshold (Relative Abundance)") +
  ylab("Taxa") +
  labs(fill = "Prevalence") +
  theme_minimal() +  # A clean theme often used in publications
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  scale_fill_gradientn(colors = brewer.pal(9, "Blues"))  # Color-blind friendly palette

# Print the plot
print(p)



############core microbiome using a venn diagram 
library(ggvenn)


ps.dataset.cm <- microbiome::transform(ps_mf_nt_vtill, "compositional")
ps.NT <- subset_samples(ps.dataset.cm, management == "NT")
ps.vtill <- subset_samples(ps.dataset.cm, management == "vtill")

core.taxa.NT <- core_members(ps.NT, detection = 0.1/100, prevalence = 50/100)
core.taxa.vtill <- core_members(ps.vtill, detection = 0.1/100, prevalence = 50/100)

cores.list <- list(NT = core.taxa.NT, vtill = core.taxa.vtill)
cores.list

tax.table.ps <- as.data.frame(tax_table(ps_mf_nt_vtill))

shared.taxa <-  core.taxa.NT[core.taxa.NT %in% core.taxa.vtill]
tax <- tax.table.ps[which(row.names(tax.table.ps) %in% shared.taxa ),]


unique.taxa.NT <-  core.taxa.NT[!core.taxa.NT %in% core.taxa.vtill]
tax.NT <- tax.table.ps[which(row.names(tax.table.ps) %in% unique.taxa.NT),]

unique.taxa.vtill <-  core.taxa.vtill[!core.taxa.vtill %in% core.taxa.NT]
tax.vtill <- tax.table.ps[which(row.names(tax.table.ps) %in% unique.taxa.vtill),]


plot <- ggvenn(cores.list, fill_color = c("#0073C2FF", "red"),
               stroke_size = 0.5, set_name_size = 15, text_size = 15,
               show_percentage = T)
print(plot)
# Load necessary libraries
library(phyloseq)
library(ggplot2)
install.packages("ape")
library(ape)
install.packages("WGCNA")
library(WGCNA)
install.packages("igraph")
library(igraph)
install.packages("ggraph")
library(ggraph)

otu_df <- as.data.frame(otu_table(ps))

View(otu_table(ps))
View(otu_df)
transposed <- t(otu_df)
View(transposed)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomformat")
library(biomformat)

otu_biom <- make_biom(transposed, sample_metadata = NULL, observation_metadata = NULL,
                      id = NULL, matrix_element_type = "int")

biom_file <- "/fs/ess/PAS1212/Afaf_rich/ITS/results/cutadapt/ITS/results/dada/otu3_biom.biom"
setwd("/fs/ess/PAS1212/Afaf_rich/ITS/results/cutadapt/ITS/results/dada")
# Write the biom file
write_biom(otu_biom, biom_file)

seqtab <- readRDS(file.path(indir, "seqtab.rds"))
seqs <- colnames(seqtab)
otab <- otu_table(seqtab, taxa_are_rows=FALSE)
colnames(otab) <- paste0("ASV", seq(ncol(otab)))
otab = t(otab)
write.table(seqs, "dada_seqs.txt",quote=FALSE)
write.table(otab, "dada_table.txt",quote=FALSE,sep="\t")

file_content <- readLines("dada_table.txt")

###################################################################################
##correlation heatmap

grouping_column <- "management"
ps_mf_all <- subset_samples(ps, treatment %in% "SD")
b <- taxa.env.correlation(ps_mf_all, grouping_column, method = "spearman",
                          pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 1,
                          num.taxa = 50, select.variables = NULL)

plot_taxa_env(b) +
  theme_pubr() +  # Apply ggpubr theme for publication quality
  labs(title = "Taxa-Environment Correlation",  # Add plot title
       x = "Spearma Correlation Coefficient",  # Label x-axis
       y = "Taxa") +  # Add caption
  theme(plot.title = element_text(size = 20, hjust = 0.5),  # Adjust title size and position
        axis.title.x = element_text(size = 15),  # Adjust x-axis label size
        axis.title.y = element_text(size = 22),  # Adjust y-axis label size
        axis.text.x = element_text(size = 18, angle = 90, hjust = 1),   # Adjust x-axis tick label size
        axis.text.y = element_text(size = 10),   # Adjust y-axis tick label size
        legend.title = element_text(size = 20),  # Adjust legend title size
        legend.text = element_text(size = 12))   # Adjust legend text size
        
        
#############plot UMAP

otu <- otu_table(ps)
samples <- sample_data(ps)
#otu <- t(otu)
umap <- uwot::umap(otu)
umapdata <- as.data.frame(umap) %>%
  rownames_to_column("SampleID") %>%
  rename(UMAP1 = V1, UMAP2 = V2)

umapdata <- left_join(umapdata, samples, by = "SampleID")
ggplot(umapdata, aes(x = UMAP1, y = UMAP2, color = management)) +
  geom_point() +
  theme(aspect.ratio = 1)


###spearman
grouping_column <- "management"
ps_mf_all <- subset_samples(ps, treatment %in% "SD")
b <- taxa.env.correlation(ps_mf_nt_vtill, grouping_column, method = "spearman",
                          pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 1,
                          select.variables = NULL)

plot_taxa_env(b) +
  #theme_pubr() +  # Apply ggpubr theme for publication quality
  labs(title = "Taxa-Environment Correlation",  # Add plot title
       x = "Spearman Correlation Coefficient",  # Label x-axis
       y = "Taxa") +  # Add caption
  theme(plot.title = element_text(size = 20, hjust = 0.5),  # Adjust title size and position
        axis.title.x = element_text(size = 15),  # Adjust x-axis label size
        axis.title.y = element_text(size = 22),  # Adjust y-axis label size
        axis.text.x = element_text(size = 20, angle = 90, hjust = 1),   # Adjust x-axis tick label size
        axis.text.y = element_text(size = 10),   # Adjust y-axis tick label size
        legend.title = element_text(size = 20),  # Adjust legend title size
        legend.text = element_text(size = 12))
        
