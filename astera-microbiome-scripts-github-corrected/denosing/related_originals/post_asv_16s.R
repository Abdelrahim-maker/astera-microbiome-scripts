##post ASV analyses

# Set the R library to load packages from
.libPaths("/fs/ess/PAS0471/jelmer/R/metabar")
dyn.load("/fs/ess/PAS0471/jelmer/software/GLPK/lib/libglpk.so.40", local = FALSE)
# Load the packages (package startup messages are not printed below)
library(vegan)
library(tidyverse)
library(phyloseq)
library(microViz)
# Set the ggplot plot theme
theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())

# Set a random seed for ordination analyses
set.seed(3131)
# Dir with input files:
indir <- "/fs/ess/PAS1212/Afaf_rich/16srRNA_2/results"

# Dir for output:
outdir <- "/fs/ess/PAS1212/Afaf_rich/16srRNA_2/results/post_asv"

ps <- readRDS(file.path(indir, "ps.rds"))

ps

# And the assigned taxon names:
taxa_names(ps)[1:3]

sort(sample_sums(ps))

# Remove samples with low counts:
ps <- subset_samples(ps, sample_sums(ps) > 1000)

# The sample data (metadata)
sample_data(ps)[1:10]
# The count table
head(t(otu_table(ps)[1:10]))

# The taxonomy table
head(tax_table(ps))

bacteria <- subset_taxa(ps, Kingdom == "Bacteria")
bacteria

bacteria_genera <- bacteria %>%
  tax_glom(taxrank = "Genus") %>%                      # agglomerate at Genus level
  transform_sample_counts(function(x) {x / sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.03) %>%                         # Filter out low abundance taxa
  arrange(Genus)                                       # Sort data frame alphabetically by Genus

# Relative abundance by sample and rotation treatment
ggplot(bacteria_genera, aes(x = Sample, y = Abundance, fill = Genus)) + 
  facet_grid(vars(management), scales = "free_y") +
  coord_flip() +
  geom_col() +
  scale_fill_discrete() +
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
  theme(legend.position = "top",  text = element_text(size = 15)) 

plot_richness(ps, x = "management", measures = c("Observed", "Shannon")) +
  geom_boxplot() +
  theme(legend.position = "top", text = element_text(size =15))

plot_richness(ps, x = "treatment", measures = c("Observed", "Shannon")) +
  geom_boxplot()  


# Make alpha diversity figure
# Plot richness with legend
# Plot richness with automatic legend
plot_richness(ps,
              x = "management",
              measures = c("Observed", "Shannon"),
              color = "management",
              nrow = 1) +
  geom_boxplot() +
  geom_jitter(width = 0.05) +
  theme(legend.position = "right") +  # Change legend position as needed
  labs(x = "Treatment", y = "Diversity")  # Axis labels



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
  coord_fixed() +
  stat_ellipse(type = "t") 

# Subsetting based on management, treatment, and depth
ps_subset <- subset_samples(ps, management %in% c("NT", "vtill") & treatment %in% "MF")

plot_ordination(ps_subset, ord.nmds.jsd_slv,
                color = "management") +
  geom_point(size = 4) +
  coord_fixed() +
  stat_ellipse(type = "t") 


wilcox.shannon <- pairwise.wilcox.test(rich$Shannon, 
                                       sample_data(ps_subset)$management, 
                                       p.adjust.method = "BH")
wilcox.shannon


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
  theme(text = element_text(size =20))


# Storing our metadata in a data frame
sampledf <- data.frame(sample_data(ps))

# Selecting our distance method (in this case the same as we used for NMDS,
# Jensen-Shannon divergence, but we could use different options here)
dist.mat <- phyloseq::distance(ps, method = "jsd")

# Selecting the number of permutations to perform
perm <- how(nperm = 999)
# Selecting blocks from the experiment
#setBlocks(perm) <- with(sampledf)

# Run PERMANOVA with the selected distance method and your variables of interest
# (in this case, location and rotation)
adonis_results <- adonis2(dist.mat ~ management* treatment,
        data = sampledf, permutations = perm)





# Load the packages (package startup messages are not printed below)
library(DESeq2)
library(phyloseq)
library(microViz)
library(ggplot2)

# network analysis
# Construct the network 

# Load the packages (package startup messages are not printed below)
library(ggnetwork)
library(phyloseq)
library(NetCoMi)
# Get group names to parse onto network construction 
group_vec <- get_variable(ps, "management")

bac22rot_net <- netConstruct(ps,
                             group = group_vec,
                             filtTax = "highestFreq",
                             filtTaxPar = list(highestFreq = 50),
                             filtSamp = "totalReads",
                             filtSampPar = list(totalReads = 500),
                             zeroMethod = "pseudo",
                             normMethod = "clr",
                             measure = "pearson",
                             sparsMethod = "none",
                             dissFunc = "signed",
                             verbose = 0,
                             seed = 123456)
bac22rot_net_aly <- netAnalyze(bac22rot_net,
                               clustMethod = "cluster_fast_greedy",
                               graphlet = FALSE,
                               verbose = 0)
