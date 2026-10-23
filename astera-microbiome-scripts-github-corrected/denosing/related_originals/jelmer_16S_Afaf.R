###load libraries
# Set the R library to load packages from
.libPaths("/fs/ess/PAS0471/jelmer/R/metabar")
dyn.load("/fs/ess/PAS0471/jelmer/software/GLPK/lib/libglpk.so.40", local = FALSE)
# Load the packages (package startup messages are not printed below)
library(vegan)
library(tidyverse)
library(phyloseq)
library(microViz)
library(Biostrings)
library(ggplot2)

# Set the ggplot plot theme
theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())

# Set a random seed for ordination analyses
set.seed(3131)

# Dir with input files:
indir <- "/fs/ess/PAS1212/Afaf_rich/16srRNA_2"

# Dir for output:
outdir <- "/fs/ess/PAS1212/Afaf_rich/16srRNA_2"
ps <- readRDS(file.path(indir, "ps.rds"))


###renaming the otus: Asv1, ASV2,...

dna <- Biostrings::DNAStringSet(taxa_names(ps))

names(dna) <- taxa_names(ps)

ps <- merge_phyloseq(ps, dna)

taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

ps



# Create ps subsets for chloroplast, mitochondria, and eukaryotes:
chlr <- subset_taxa(ps, Order == "Chloroplast")
mit <- subset_taxa(ps, Family == "Mitochondria")
euk <- subset_taxa(ps, Kingdom == "Eukaryota")

# Get all taxon IDs that we want to remove:
bad_taxa <- c(taxa_names(chlr), taxa_names(mit), taxa_names(euk))

# Get all taxon IDs that we want to keep:
all_taxa <- taxa_names(ps)
good_taxa <- all_taxa[!(all_taxa %in% bad_taxa)]


# Subset phyloseq object:
ps_bac_arc <- prune_taxa(good_taxa, ps)

pre <- sum(sample_sums(ps))
post <- sum(sample_sums(ps))

post / pre

#post <- ps


# Subset phyloseq object:
ps_bac_arc <- prune_taxa(good_taxa, ps)
ps_bac_arc <- ps
# And the assigned taxon names:
taxa_names(ps)[1:3]
sort(sample_sums(ps))

##correct the darr-cod typo
ps@sam_data$site <- gsub("Darr_CoD", "Darr-COD", ps@sam_data$site)

# Remove samples with low counts: (do not use this here bec u lose samples)
ps <- subset_samples(ps, sample_sums(ps) > 500)
#ps <- subset_taxa(ps, taxa_sums(ps) > 0)

# The sample data (metadata)
sample_data(ps)[1:10]

# The taxonomy table
head(tax_table(ps))


bacteria <- subset_taxa(ps, Kingdom == "Bacteria")
bacteria


###extract otu table----

# Extract the OTU table from the phyloseq object
otu_table <- as(otu_table(ps), "matrix")

# View the OTU table
head(otu_table)

###transform otu table for lefse analysis
# Modify the values in the "management" column
ps@sam_data$management <- gsub("LNT", "NT", ps@sam_data$management)
ps@sam_data$management <- gsub("SNT", "NT", ps@sam_data$management)

###subset the phyloseq by MF and SD for stats as well as MS and SD

ps_MF_SD <- subset_samples(ps, treatment %in% c("MF", "SD"))
ps_MS_SD <- subset_samples(ps, treatment %in% c("MS", "SD"))
ps_MF_MS <- subset_samples(ps, treatment %in% c("MF", "MS"))



physeq_all_depth_MF <- subset_samples(ps, treatment %in% "MF" & management %in% c("vtill", "NT"))
otu_table_all_depth_MF <- as(otu_table(physeq_all_depth_MF), "matrix")
otu_table_all_depth_MF <- t(otu_table_all_depth_MF)
write.csv(otu_table_all_depth_MF, file = file.path(indir, "otu_table_all_depth_MF_lefse.csv"))

physeq_10_MF <- subset_samples(ps, treatment %in% "MF" & depth == 10)
otu_table_10_MF <- as(otu_table(physeq_10_MF), "matrix")
otu_table_10_MF <- t(otu_table_10_MF)
write.csv(otu_table_10_MF, file = file.path(indir, "otu_table_10_MF_lefse.csv"))

otu_table_all_data_3_treat <- as(otu_table(ps_filtered_2), "matrix")
otu_table_all_data_3_treat <- t(otu_table_all_data_3_treat)
write.csv(otu_table_all_data_3_treat, file = file.path(indir, "otu_table_all_data_3_treat.csv"))





physeq_subset_treat_lefse <- subset_samples(ps, depth == 10 & management %in% c("vtill", "NT"))
physeq_subset_treat_2_lefse <- subset_samples(ps, management %in% c("vtill", "NT"))
otu_table_treat <- as(otu_table(physeq_subset_treat_lefse), "matrix")
otu_table_treat <- t(otu_table_treat) ##save this otu table to a lefse directory and only keep samples with 3 treatments

otu_table_treat_2 <- as(otu_table(physeq_subset_treat_2_lefse), "matrix")
otu_table_treat_2 <- t(otu_table_treat_2)
write.csv(otu_table_treat_2, file = file.path(indir, "otu_table_treat_all_depth_lefse.csv"))


write.csv(otu_table_treat, file = file.path(indir, "otu_table_treat_lefse.csv"))

physeq_subset_MF10_lefse <- subset_samples(ps, depth == 10 & treatment %in% "MF" & management %in% c("vtill", "NT"))

otu_table_MF10_lefse <- as(otu_table(physeq_subset_MF10_lefse), "matrix")
otu_table_MF10_lefse <- t(otu_table_MF10_lefse)
write.csv(otu_table_MF10_lefse, file = file.path(indir, "otu_table_MF10_lefse.csv"))

###remove the samples with one treatment from physeq_subset_treat_lefse for ordination and stats
samples_to_remove <- c("cs85_SD_10_LNT_S110_L001", "cs85_MS_10_LNT_S179_L001")

# Remove the specified samples
physeq_treatment_ordi_stats <- subset_samples(physeq_subset_treat_lefse, !(sample_names(physeq) %in% samples_to_remove))

##remove samples from phyloseq_subset_treat_2_lefse (all depth)
patterns_to_remove <- c("^cs102_","^cs15_", "^cs20_", "^cs24_", "^cs26_", "^cs43_", 
   "^cs96_", "^cs97_", "^cs57_", "^cs5_", 
  "^cs66_", "^cs60_", "^cs67_", "^cs72_", "^cs93_", "^cs99_",  
   "^cs47_", "^cs48_","^cs51_", "^cs53_", "^cs59_", "^cs6_", "^cs62_","^cs73_",  "^cs77_",  "^cs78_",
  "^cs79_",  "^cs8_",  "^cs84_",  "^cs85_",  "^cs89_",  "^cs90_",  "^cs99_", "^cs40_", "^cs39_", "^cs37_",
  "^cs100_", "^cs16_", "^cs27_", "^cs3_", "^cs38_", "^cs18_", "^cs19_", "^cs22_", "^cs23_", "^cs17_")  
pattern <- paste0("^(", paste(patterns_to_remove, collapse = "|"), ")")

# Subset the phyloseq object to remove samples matching the patterns
ps_filtered_2 <- subset_samples(ps, !grepl(pattern, sample_names(ps)))












bacteria_genera <- bacteria %>%
  tax_glom(taxrank = "Genus") %>%                      # agglomerate at Genus level
  transform_sample_counts(function(x) {x / sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.03) %>%                         # Filter out low abundance taxa
  arrange(Genus)       


# Relative abundance by sample and rotation treatment
library(RColorBrewer)
library(viridis)
palette <- plasma(length(unique(bacteria_genera$Genus)))
ggplot(bacteria_genera, aes(x = Sample, y = Abundance, fill = Genus)) + 
  facet_grid(vars(management), scales = "free_y") +
  coord_flip() +
  geom_col() +
  scale_fill_manual(values = palette) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  labs(y = "Relative Abundance (Genera > 2%)",
       title = "Relative Abundance of Bacteria",
       x = NULL) +
  theme(axis.text.y = element_blank(),  # Remove y-axis labels
        legend.text = element_text(size = 12),  # Adjust legend label size
        legend.title = element_text(size = 15, face = "bold"),  # Adjust legend title size
        axis.text.x = element_text(size = 15), strip.text = element_text(size = 14))

# plot multiple indices 
plot_richness(ps,
              x = "management",
              color = "site", 
              measures = c("Observed", "Shannon", "Simpson", "InvSimpson")) +
  theme(legend.position = "top")

# Make alpha diversity figure
plot_richness(physeq_treatment_ordi_stats,
              x = "treatment",
              measures = c("Observed", "Shannon"),
              color = "management",
              nrow = 1) +
  geom_boxplot() +
  geom_jitter(width = 0.05) +
  theme(legend.position = "none") +
  labs(x = "management", y = "Diversity")

##compare diversity between managements and sites and treatments
# Calculate diversity metrics using estimate_richness and store the data using a vector
rich <- estimate_richness(ps, measures = c("Observed", "Shannon"))

# Compare observed diversity between our sample location
wilcox.observed <- pairwise.wilcox.test(rich$Observed, 
                                        sample_data(physeq_all_depth_MF)$depth, 
                                        p.adjust.method = "BH")
wilcox.observed


# Non-metric MultiDimensional Scaling (NMDS)
ord.nmds.jsd_slv <- ordinate(ps_object_corrected, method = "NMDS", distance = "jsd")
# Create a "stressplot"
stressplot(ord.nmds.jsd_slv)
plot_ordination(ps_object_corrected, ord.nmds.jsd_slv,
                color = "site", shape = "management") +
  geom_point(size = 4) +
  coord_fixed() 
#stat_ellipse(type = "t") 

# Redundancy Analysis
ps_subset_mf <- subset_samples(physeq_all_depth_MF)

# Subset the sample column based on samples with three different treatments
# Assuming 'ps' is your phyloseq object

# Get the sample data
sample_data <- sample_data(physeq_treat_all_depth_stats)



# Modify the values in the "management" column
ps@sam_data$management <- gsub("LNT", "NT", ps@sam_data$management)
ps@sam_data$management <- gsub("SNT", "NT", ps@sam_data$management)

physeq_subset_mf <- subset_samples(ps, treatment %in% "MF")
# Save the modified Phyloseq object to a new name
physeq_modified <- ps

physeq_subset_mf <- subset_samples(physeq_modified, treatment %in% "MF"& depth == 30)




ord.rda <- ordinate(physeq_subset_mf ,
                    formula = o ~ management + site,
                    method = "NMDS",
                    distance = "bray")
# plot the RDA
#ps_subset_1 <- subset_samples(ps, management %in% c('vtill', 'NT'), treatment %in% "MF")
management_colors <- c("vtill" = "brown", "NT" = "blue", "forest" = "forestgreen")
plot_ordination(ps,
                                   vegan::scores(ord.rda, scaling = 1),
                                   type = "sites",
                                   color = "management") +
  geom_point(size = 5) +
  stat_ellipse(level = 0.95, geom = "polygon", alpha = 0.2, fill = "#CCCCCC") +
  scale_color_manual(values = management_colors) +
  theme(axis.text.x = element_text(size = 20),  # Adjust x-axis text size
        axis.text.y = element_text(size = 20),  # Adjust y-axis text size
        legend.text = element_text(size = 30)) +  # Adjust legend text size
  theme(legend.position = "bottom") +  # Position legend at the bottom
  theme(panel.grid = element_blank()) +  # Remove grid lines
  theme(legend.key = element_rect(fill = "transparent")) 


# Print the modified plot
print(ordination_plot)
install.packages("vegan")
library(vegan)
# Storing our metadata in a data frame
ps_mf_10 <- subset_samples(ps, treatment %in% "MF" & depth == 10)
sampledf <- data.frame(sample_data(ps))

# Selecting our distance method (in this case the same as we used for NMDS,
# Jensen-Shannon divergence, but we could use different options here)
dist.mat <- phyloseq::distance(ps, method = "jsd")

# Selecting the number of permutations to perform
perm <- how(nperm = 999)
# Selecting blocks from the experiment
#setBlocks(perm) <- with(sampledf, site)

# Run PERMANOVA with the selected distance method and your variables of interest
# (in this case, location and rotation)
vegan::adonis2(dist.mat ~ management * depth * treatment,
        data = sampledf, permutations = perm)
vegan::adonis2(dist.mat ~ management/site * treatment, data = sampledf, permutations = perm)

####################################################################################################3

####Differential analysis using deseq2

library(DESeq2)
library(microViz)

# Converting phyloseq object to DESeq2 format
ps_mf <- subset_samples(ps, management %in% c("NT", "vtill")& treatment %in% "MF")

kostic = subset_samples(ps_mf)
kostic
diagdds = phyloseq_to_deseq2(kostic, ~ management)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]


sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kostic)[rownames(sigtab), ], "matrix"))
sigtab <- cbind(sigtab, sample_data(ps_mf)[rownames(sigtab), c("management")])
sigtab <- cbind(sigtab, sample_data(ps_mf)[rownames(sigtab), c("management")])
head(sigtab)
dim(sigtab)

# Set theme
theme_set(theme_bw())

# Define function to set color scale using brewer palette
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

# Phylum order
x_phylum <- tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x_phylum <- sort(x_phylum, decreasing = TRUE)
sigtab$Phylum <- factor(as.character(sigtab$Phylum), levels = names(x_phylum))

# Genus order
x_genus <- tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x_genus <- sort(x_genus, decreasing = TRUE)
sigtab$Genus <- factor(as.character(sigtab$Genus), levels = names(x_genus))

# Define color-blind friendly color palette
color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#8DD3C7", "#FFD92F", "#E6194B", "#3CB44B")

# Plot with improved aesthetics
ggplot(sigtab, aes(x = Genus, y = log2FoldChange, color = Phylum)) +
  geom_point(size = 3) + 
  scale_color_manual(values = color_palette) +  # Manually set color palette
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
        legend.title = element_text(size = 18, face = "bold"),  # Adjust legend title appearance
        legend.text = element_text(size = 18),  # Adjust legend text appearance
        legend.key.size = unit(0.5, "cm"))  # Adjust legend key size



#####################################################3333

##Try deseq to get the diff abundant taxa in each treatment---

sample_data(ps_mf)$management <- factor(sample_data(ps_mf)$management)

# DESeq2 conversion and analysis with DIAGNOSIS variable
diagdds <- phyloseq_to_deseq2(ps_mf, ~ management)
diagdds <- DESeq(diagdds, test = "Wald", fitType = "parametric")

# Results Interpretation
res <- results(diagdds, cooksCutoff = FALSE)
alpha <- 0.01
sigtab <- res[which(res$padj < alpha), ]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(ps_mf)[rownames(sigtab), ], "matrix"))

# Set theme
theme_set(theme_bw())

# Define color-blind friendly color palette
color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#8DD3C7", "#FFD92F", "#E6194B", "#3CB44B")

# Plot with facet wrap by DIAGNOSIS
ggplot(sigtab, aes(x = Genus, y = log2FoldChange, color = Phylum)) +
  geom_point(size = 3) +
  scale_color_manual(values = color_palette) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18),
        legend.key.size = unit(0.5, "cm")) +
  facet_wrap(~ management)



################workshop code for deseq (this one is better and it worked perfectly)
library(DESeq2)
set.seed(1234)
ps_mf_all_depth <- subset_samples(ps, treatment %in% "MF" & management %in% c("vtill", "NT"))
otu_table(ps_mf_all_depth) <- otu_table(ps_mf_all_depth) + 1
otu_table(ps) <- otu_table(ps) + 1

dds.df <- phyloseq_to_deseq2(ps, ~management)

dds.df$Group <- factor(paste0(dds.df$management))  

# Setting our experimental desing in the DESeq2 object
design(dds.df) <- formula(~Group) 

# Generating our analysis
dds <- DESeq(dds.df, test="Wald", fitType="parametric") 

results(dds)

levels(dds$Group)

# Select the groups to contrast 
res <- results(dds, contrast = c("Group", "NT", "vtill"), cooksCutoff = FALSE) 

# First look at the results 
summary(res) 
#alpha = 0.01


# Number of taxa with significantly higher abundance in NT :
sum(res$padj < alpha & res$log2FoldChange > 0, na.rm = TRUE)

# Number of taxa with significantly lower abundance in NT:
sum(res$padj < alpha & res$log2FoldChange < 0, na.rm = TRUE)
res.order <- res[order(res$pvalue), ]
head(res.order)[c(1:2, 6)]


# Create a data frame from your DESeq results 
res.df <- as.data.frame(res)

# Add the taxonomy of the ASVs to the new data frame 
taxmelt <- as.data.frame(tax_table(ps)) # , taxmelt$Family, sep= " " include this to add taxonomy
res.df$Tax <- paste(rownames(res.df))

# Filter to keep only log2folds > |10| and pvalues <= 0.05 
res05 <- dplyr::filter(res.df, padj<=0.05, log2FoldChange<=-4 | log2FoldChange>= 4) 
# Generating a dataframe from our table
res05.df <- as.data.frame(res05)
# Define color palette
color_palette <- c("red", "#9AD8F7") ##this blue shade is for MS "#9AD8F7" and this navy blue is for MF "#0072B2" and the brown is for SD "#8B4513"


# Plot
# Plot
ggplot(data = res05.df,
       aes(x = reorder(Tax, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat="identity", width = 0.7) +
  scale_fill_manual(values = color_palette, guide = FALSE) +
  coord_flip() +
  labs(title = "Taxa with significantly lower and higher abundance in MS when contrasted with SD",
       x = "Lineages",
       y = "log2FoldChange") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 5),
        axis.text.x = element_text(size = 10), # Adjust y-axis text size
        axis.text.y.left = element_text(vjust = 0.5),  # Rotate y-axis text
        legend.position = "none") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 25)) 



###########################################################################################

###Core Microbiome
library(microbiome)
library(ggvenn)
# Make sure we use functions from correct package
transform <- microbiome::transform

# Merge rare taxa to speed up examples
pseq <- transform(ps_mf, "compositional")
pseq <- aggregate_rare(pseq, level = "Genus", detection = 1/100, prevalence = 50/100)




# Pick sample subset
library(phyloseq)
pseq2 <- subset_samples(ps_mf)

# Normal western adults
data(atlas1006)
pseq3 <- atlas1006 %>%
  sggregate_taxa(level = "Phylum") %>%  
  microbiome::transform(transform = "compositional")





#############################################################################################
###network
# Load the packages (package startup messages are not printed below)
library(ggnetwork)
library(phyloseq)
library(NetCoMi)

# Get group names to parse onto network construction 
ps_mf_all_depth <- subset_samples(ps_mf, treatment %in% "MF" & depth == 10 )
ps_subset_1 <- subset_samples(ps, treatment  %in% "MF" & management %in% c("vtill", "NT"))
group_vec <- get_variable(ps_subset_1, "management")

# Construct the network 
bac22rot_net <- netConstruct(ps_subset_1,
                             group = group_vec,
                             filtTax = "highestFreq",
                             filtTaxPar = list(highestFreq = NULL),
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
# Load necessary libraries
library(RColorBrewer)
col <- brewer.pal(4, "PuOr")

# Plot network

# Step 1: Extract taxa names from tax_table()

plot(bac22rot_net_aly, 
     colorVec = col,
     borderCol = "gray40", 
     nodeSize = "mclr", 
     nodeSizeSpread = 3, 
     labelScale = TRUE,
     edgeTranspLow = 80, 
     edgeTranspHigh = 50,
     groupNames = c("vtill", "NT"), 
     showTitle = TRUE, 
     cexTitle = 3,
     cexHubLabels = 1.5, 
     cexNodes = 5,
     hubBorderCol = "darkgray",##this can show u where are the hubs in each network
     mar = c(1,1,3,1), 
     repulsion = 0.9, 
     shortenLabels = "simple", 
     rmSingles = "inboth",
     nodeFilter = "clustMin",
     #nodeColor = "feature",
     nodeFilterPar = 10, 
     nodeTransp = 50, 
     hubTransp = 30) 






##differential networks

# Subsetting our data into their respective Rotations 
vtill_subset <- phyloseq::subset_samples(ps_mf, depth == 10 & management != "NT")
NT_subset  <- phyloseq::subset_samples(ps_mf, depth == 10 & management == "NT")


# Construct the differential network 
diff_net_const <- netConstruct(vtill_subset, NT_subset,
                               filtTax = "highestVar",
                               filtTaxPar = list(highestVar = NULL),
                               zeroMethod = "pseudo",
                               normMethod = "clr",
                               measure = "sparcc",
                               verbose = 1, 
                               thresh = 0.3,
                               sparsMethod = "none")

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
        groupNames = c("vtill", "NT"),
        showCentr = c("degree", "between", "closeness"), 
        numbNodes = 5)
plot(diff_net_const_aly, 
     sameLayout = TRUE, 
     layoutGroup = "union", 
     colorVec = col,
     borderCol = "gray40", 
     nodeSize = "mclr", 
     cexNodes = 0.9,
     nodeSizeSpread = 3, 
     labelScale = TRUE,
     edgeTranspLow = 80, 
     edgeTranspHigh = 50,
     groupNames = c("vtill", "NT"), 
     showTitle = TRUE, 
     cexTitle = 2.8,
     cexHubLabels = 6, 
     #cexNodes = 3,
     mar = c(1,1,3,1), 
     repulsion = 0.9, 
     labels = TRUE,
     shortenLabels = "simple", 
     rmSingles = "inboth",
     nodeFilter = "clustMin", 
     nodeFilterPar = 5,  # Adjust the minimum cluster size
     nodeTransp = 50, 
     hubTransp = 30)

##################try to build each network individually
ps_mf_10_vtill <- subset_samples(ps_mf, management %in% "vtill" & depth == 10)
net_diss <- netConstruct(ps_mf_10_vtill,
                         measure = "aitchison",
                         zeroMethod = "multRepl",
                         sparsMethod = "knn",
                         kNeighbor = 3,
                         verbose = 3)
props_diss <- netAnalyze(net_diss,
                         clustMethod = "hierarchical",
                         clustPar = list(method = "average", k = 3),
                         hubPar = "eigenvector")
plot(props_diss, 
     nodeColor = "cluster", 
     nodeSize = "eigenvector",
     hubTransp = 40,
     edgeTranspLow = 60,
     charToRm = "00000",
     shortenLabels = "simple",
     labelLength = 6,
     mar = c(1, 3, 3, 5))
green2 <- colToTransp("#009900", 40)

legend(0.4, 1.1,
       cex = 2.2,
       legend = c("high similarity (low Aitchison distance)",
                  "low similarity (high Aitchison distance)"), 
       lty = 1, 
       lwd = c(3, 1),
       col = c("darkgreen", green2),
       bty = "n")
