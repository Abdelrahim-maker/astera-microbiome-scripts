# Load the packages
install.packages("dada2")
library(dada2)
install.packages("phyloseq")
library(phyloseq)
library(tidyverse)
library(vegan)
library(microViz)
library(ggplot2)
install.packages("ggpubr")
library(ggpubr)
install.packages("ggsignif")
library(ggsignif)
install.pa
library(microbiomeSeq)
library(microbiomeutilities)
library(microeco)
install.packages("dada2")
library(dada2)
library(ggplot2)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq", force = TRUE)

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
ps@sam_data$crop <- gsub("soybean_winter_wheat", "winterwheat", ps@sam_data$crop)

sample_data(ps)

######normalization
#1-normalize OTU abundance
physeq<-normalise_data(ps,norm.method = "randomsubsample")
#physeq <- normalise_data(ps, norm.method = "varstab" ,fitType="local")

# 2- normalization of sample data
physeq <- normalise_data(physeq, norm.method = "relative", norm.meta=T)

# 3- scale sample data
physeq <- normalise_data(physeq, norm.method = "scale", type="sqrt")

####Alpha diversity
# Define the pairwise comparisons
a_my_comparisons <- list(c("MF", "MS"), c("MF", "SD"), c("MS", "SD"))

# Define the color palette for treatments
treatment_colors <- c("MF" = "#1f78b4", "MS" = "#a6cee3", "SD" = "#8c510a")

# Define the cutpoints and symbols for significance levels
cutpoints <- c(0.0001, 0.001, 0.01, 0.05, 1)
symbols <- c("a", "ab", "b", "ns")

# Set the symnum arguments
symnum.args <- list(cutpoints = cutpoints, symbols = symbols)

# Plot richness and add comparisons
plot_richness(ps, x = "treatment", measures = "Shannon", color = "treatment") +
  geom_boxplot(alpha = 0.6, fill = treatment_colors) + 
  scale_color_manual(values = treatment_colors) +
  theme_pubr() + # Apply ggpubr theme for publication quality
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args) +
  # Modify axis labels and text size
  labs(y = "Shannon Index", x = "") +
  theme(axis.title.x = element_text(size = 28), # Increase x-axis label size
        axis.title.y = element_text(size = 28), # Increase y-axis label size
        axis.text.x = element_text(size = 25),  # Increase x-axis text size
        axis.text.y = element_text(size = 25),  # Increase y-axis text size
        legend.text = element_text(size = 25),  # Increase legend text size
        legend.title = element_text(size = 25), # Increase legend title size
        strip.text = element_blank())           # Remove the strip text
###################################################

# Create rarefaction curve
rarecurve_plot <- rarecurve(t(otu_table(ps)), step = 50, cex = 0.5)

# Customize the plot
rarecurve_plot +
  geom_line(color = "steelblue", size = 1, aes(linetype = "Curve")) +
  scale_linetype_manual(values = "solid") +  # Set line type to solid
  labs(x = "Number of Sequences", y = "Number of OTUs", title = "Rarefaction Curve") +
  theme_minimal(base_size = 14) +  # Use a minimal theme with larger base font size
  theme(axis.title = element_text(size = 16),  # Increase axis title size
        axis.text = element_text(size = 14),   # Increase axis text size
        plot.title = element_text(size = 18))  # 
####################################################################################
##
grouping_column <- "management"
b <- taxa.env.correlation(ps, grouping_column, method = "pearson",
                     pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 1,
                     num.taxa = 50, select.variables = NULL)

plot <- plot_taxa_env(b) +
  theme_pubr() +  # Apply ggpubr theme for publication quality
  labs(title = "Taxa-Environment Correlation",  # Add plot title
       x = "Pearson Correlation Coefficient",  # Label x-axis
       y = "Taxa") +  # Add caption
  theme(plot.title = element_text(size = 20, hjust = 0.5),  # Adjust title size and position
        axis.title.x = element_text(size = 15),  # Adjust x-axis label size
        axis.title.y = element_text(size = 22),  # Adjust y-axis label size
        axis.text.x = element_text(size = 18, angle = 90, hjust = 1),   # Adjust x-axis tick label size
        axis.text.y = element_text(size = 10),   # Adjust y-axis tick label size
        legend.title = element_text(size = 20),  # Adjust legend title size
        legend.text = element_text(size = 12))   # Adjust legend text size

print(plot)
# Print the plot
print(plot)

###############################################################################
##plot abundance
   
physeq_till_NT <- subset_samples(ps, management %in% c("NT", "vtill"))


physeq <- taxa_level(physeq_, which_level = "Species")
physeq <- normalise_data(physeq, norm.method = "relative")
plot_taxa(physeq,grouping_column="management",method="hellinge",number.taxa=21,filename=NULL)
# Modify the plot to remove sample names from x-axis
p <- last_plot() + theme(axis.text.x = element_blank(),  # Remove x-axis labels
                         axis.ticks.x = element_blank()) # Remove x-axis ticks
