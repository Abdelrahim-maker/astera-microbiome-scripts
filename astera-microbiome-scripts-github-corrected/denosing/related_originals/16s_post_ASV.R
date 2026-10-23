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

# Set the ggplot plot theme
theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())

# Set a random seed for ordination analyses
set.seed(3131)

# Dir with input files:
indir <- "/fs/ess/PAS1212/Afaf_rich/16srRNA_2"

# Dir for output:
outdir <- "/fs/ess/PAS1212/Afaf_rich/16srRNA_2"
ps <- readRDS


dna <- Biostrings::DNAStringSet(taxa_names(ps))

names(dna) <- taxa_names(ps)

ps <- merge_phyloseq(ps, dna)

taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

ps

# Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(ps), "matrix")
# transpose if necessary
if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)

write.table(otu_table(ps), file = "otu_table.csv", sep = "\t", row.names = TRUE, col.names = TRUE)



# Create a new phyloseq object with corrected sample data
ps_object_corrected <- ps
sample_data(ps_object_corrected)$site <- ifelse(sample_data(ps)$site == "Darr_CoD", 
                                                "Darr-CoD", 
                                                sample_data(ps)$site)

# Check if the correction has been applied
sample_data(ps_object_corrected)
# Create ps subsets for chloroplast, mitochondria, and eukaryotes:
chlr <- subset_taxa(ps_object_corrected, Order == "Chloroplast")
mit <- subset_taxa(ps_object_corrected, Family == "Mitochondria")
euk <- subset_taxa(ps_object_corrected, Kingdom == "Eukaryota")

# Get all taxon IDs that we want to remove:
bad_taxa <- c(taxa_names(chlr), taxa_names(mit), taxa_names(euk))

# Get all taxon IDs that we want to keep:
all_taxa <- taxa_names(ps_object_corrected)
good_taxa <- all_taxa[!(all_taxa %in% bad_taxa)]
# Subset phyloseq object:
ps_bac_arc <- prune_taxa(good_taxa, ps)

pre <- sum(sample_sums(ps_object_corrected))
post <- sum(sample_sums(ps_bac_arc))
post / pre

#post <- ps


# Subset phyloseq object:
ps_bac_arc <- prune_taxa(good_taxa, ps_object_corrected)
ps_bac_arc <- ps_object_corrected
# And the assigned taxon names:
taxa_names(ps_object_corrected)[1:3]
sort(sample_sums(ps_object_corrected))

# Remove samples with low counts:
ps_object_corrected <- subset_samples(ps_object_corrected, sample_sums(ps) > 1000)

# The sample data (metadata)
sample_data(ps_object_corrected)[1:10]

# The taxonomy table
head(tax_table(ps_object_corrected))


bacteria <- subset_taxa(ps_object_corrected, Kingdom == "Bacteria")
bacteria

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
plot_richness(ps_object_corrected,
              x = "management",
              color = "site", 
              measures = c("Observed", "Shannon", "Simpson", "InvSimpson")) +
  theme(legend.position = "top")

# Make alpha diversity figure
plot_richness(ps_object_corrected,
              x = "management",
              measures = c("Observed", "Shannon"),
              color = "management",
              nrow = 1) +
  geom_boxplot() +
  geom_jitter(width = 0.05) +
  theme(legend.position = "none") +
  labs(x = "management", y = "Diversity")

##compare diversity between managements and sites and treatments
# Calculate diversity metrics using estimate_richness and store the data using a vector
rich <- estimate_richness(ps_object_corrected, measures = c("Observed", "Shannon"))

# Compare observed diversity between our sample location
wilcox.observed <- pairwise.wilcox.test(rich$Observed, 
                                        sample_data(ps)$management, 
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
ps_subset_2 <- subset_samples(ps, treatment %in% "MF")

# Modify the values in the "management" column
ps@sam_data$management <- gsub("LNT", "NT", ps@sam_data$management)
ps@sam_data$management <- gsub("SNT", "NT", ps@sam_data$management)

# Save the modified Phyloseq object to a new name
physeq_modified <- ps



ord.rda <- ordinate(physeq_modified,
                    formula = o ~ treatment + site,
                    method = "CAP",
                    distance = "bray")
# plot the RDA
#ps_subset_1 <- subset_samples(ps, management %in% c('vtill', 'NT'), treatment %in% "MF")
ordination_plot <- plot_ordination(physeq_modified,
                vegan::scores(ord.rda, scaling = 1),
                type = "sites",
                color = "site",
                shape = "treatment") +
  geom_point(size = 5) +
  theme(legend.text = element_text(size = 18)) +
  theme_classic()
  #stat_ellipse(type = "t")

# Modify point size
ordination_plot + geom_point(size = 4) +  # Adjust the size as needed
  
  # Modify legend font size
ordination_plot + theme_classic()

# Print the modified plot
print(ordination_plot)

library(vegan)
# Storing our metadata in a data frame
sampledf <- data.frame(sample_data(physeq_modified))

# Selecting our distance method (in this case the same as we used for NMDS,
# Jensen-Shannon divergence, but we could use different options here)
dist.mat <- phyloseq::distance(physeq_modified, method = "jsd")

# Selecting the number of permutations to perform
perm <- how(nperm = 999)
# Selecting blocks from the experiment
setBlocks(perm) <- with(sampledf, site)

# Run PERMANOVA with the selected distance method and your variables of interest
# (in this case, location and rotation)
adonis2(dist.mat ~ management * site * depth * treatment,
        data = sampledf, permutations = perm)

library(DESeq2)
library(microViz)

# Converting phyloseq object to DESeq2 format
ps_subset <- subset_samples(ps, management %in% c('vtill', 'NT')& treatment %in% "MF" & depth == 10)

dds.df <- phyloseq_to_deseq2(ps_subset, ~management + site)
# Identifying the groups in our experiment
dds.df$Group <- factor(paste0(dds.df$management, dds.df$site))  




###network
# Load the packages (package startup messages are not printed below)
library(ggnetwork)
library(phyloseq)
library(NetCoMi)

# Get group names to parse onto network construction 
ps_subset_1 <- subset_samples(physeq_modified, management %in% c('vtill', 'NT'))
group_vec <- get_variable(ps_subset_1, "management")

# Construct the network 
bac22rot_net <- netConstruct(ps_subset_1,
                             group = group_vec,
                             filtTax = "highestFreq",
                             filtTaxPar = list(highestFreq = 32,000),
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
     cexTitle = 2.8,
     mar = c(1,1,3,1), 
     repulsion = 0.9, 
     shortenLabels = "simple", 
     rmSingles = "inboth",
     nodeFilter = "clustMin", 
     nodeFilterPar = 10, 
     nodeTransp = 50, 
     hubTransp = 30) 


###or this code
# Get phyla names
taxtab <- as(tax_table(physeq_modified), "matrix")
phyla <- as.factor(gsub("p__", "", taxtab[, "Rank2"]))
names(p) <- taxtab[, "Rank6"]
#table(phyla)

# Define phylum colors
phylcol <- c("cyan", "blue3", "red", "lawngreen", "yellow", "deeppink")

plot(bac22rot_net_aly,
     layout = "spring",
     repulsion = 0.84,
     shortenLabels = "none",
     charToRm = "g__",
     labelScale = FALSE,
     rmSingles = TRUE,
     nodeSize = "clr",
     nodeSizeSpread = 4,
     nodeColor = "feature", 
     #featVecCol = phyla, 
     colorVec =  phylcol,
     posCol = "darkturquoise", 
     negCol = "orange",
     edgeTranspLow = 0,
     edgeTranspHigh = 40,
     cexNodes = 2,
     cexLabels = 2,
     cexHubLabels = 2.5,
     title1 = "Network on genus level with Pearson correlations", 
     showTitle = TRUE,
     cexTitle = 2.3)

# Colors used in the legend should be equally transparent as in the plot
phylcol_transp <- colToTransp(phylcol, 60)

legend(-1.2, 1.2, cex = 2, pt.cex = 2.5, title = "Phylum:", 
       legend=levels(phyla), col = phylcol_transp, bty = "n", pch = 16) 

legend(0.7, 1.1, cex = 2.2, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("darkturquoise","orange"), 
       bty = "n", horiz = TRUE)


##differential networks

# Subsetting our data into their respective Rotations 
vtill_subset <- phyloseq::subset_samples(ps_subset_1, management != "NT")
NT_subset  <- phyloseq::subset_samples(ps_subset_1, management == "NT")

# Construct the differential network 
diff_net_const <- netConstruct(vtill_subset, NT_subset,
                               filtTax = "highestVar",
                               filtTaxPar = list(highestVar = 50),
                               zeroMethod = "pseudo",
                               normMethod = "clr",
                               measure = "sparcc",
                               verbose = 1, 
                               thresh = 0.3)

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
        groupNames = c("NT", "vtill"),
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
     groupNames = c("vtill", "LNT"), 
     showTitle = TRUE, 
     cexTitle = 2.8,
     mar = c(1,1,3,1), 
     repulsion = 0.9, 
     labels = TRUE,
     shortenLabels = "simple", 
     rmSingles = "inboth",
     nodeFilter = "clustMin", 
     nodeFilterPar = 5,  # Adjust the minimum cluster size
     nodeTransp = 50, 
     hubTransp = 30)


###differential abundance using DESEq2
library(DESeq2)

# Converting phyloseq object to DESeq2 format 
dds.df <- phyloseq_to_deseq2(physeq_modified, ~management)

# Identifying the groups in our experiment
dds.df$Group <- factor(paste0(dds.df$management))  

# Setting our experimental desing in the DESeq2 object
design(dds.df) <- formula(~Group) 

# Generating our analysis
dds.df <- dds.df[rowSums(counts(dds.df)) > 0, ]
dds <- DESeq(dds.df, quiet = TRUE) 

results(dds) 

##core microbiome
# Load the packages
library(tidyverse)
library(microbiome)
library(phyloseq)
library(ggvenn)

ps.dataset.cm <- microbiome::transform(physeq_modified, "compositional")

ps.vtill <- subset_samples(ps.dataset.cm, management == "vtill")
ps.NT <- subset_samples(ps.dataset.cm, management == "NT")
core.taxa.vtill <- core_members(ps.vtill, detection = 0.001, prevalence = 95/100)
core.taxa.NT <- core_members(ps.NT, detection = 0.001, prevalence = 95/100)

cores.list <- list(vtill = core.taxa.vtill, NT = core.taxa.NT)
cores.list #To see the list of core taxa (ASVs) by Rotation


tax.table.ps <- as.data.frame(tax_table(physeq_modified))

shared.taxa <-  core.taxa.vtill[core.taxa.vtill %in% core.taxa.NT]
tax <- tax.table.ps[which(row.names(tax.table.ps) %in% shared.taxa ),]

unique.taxa.vtill <-  core.taxa.vtill[!core.taxa.vtill %in% core.taxa.NT]
tax.vtill <- tax.table.ps[which(row.names(tax.table.ps) %in% unique.taxa.vtill),]


unique.taxa.NT <-  core.taxa.NT[!core.taxa.NT %in% core.taxa.vtill]
tax.NT <- tax.table.ps[which(row.names(tax.table.ps) %in% unique.taxa.NT),]
plot <- ggvenn(cores.list, fill_color = c("#0073C2FF", "#EFC000FF"),
               stroke_size = 0.5, set_name_size = 5, text_size = 5,
               show_percentage = T)
print(plot)




