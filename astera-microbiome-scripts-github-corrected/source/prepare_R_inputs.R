# Run from the extracted package root after installing the R dependencies.
# Reconstructs legacy inputs from shared selected-cohort data; NOT executed/validated here.
library(phyloseq)
allmeta <- read.csv("data/manuscript_metadata.csv", check.names=FALSE)
files <- c(`16S`="combined_16S_ps.rds", ITS="combined_its.rds", AMF="combined_amf.rds")
for (marker in names(files)) {
  counts <- as.matrix(read.csv(gzfile(paste0("data/",marker,"_counts.csv.gz")),row.names=1,check.names=FALSE))
  taxonomy <- as.matrix(read.csv(paste0("data/",marker,"_taxonomy.csv"),row.names=1,check.names=FALSE))
  m <- allmeta[allmeta$marker==marker,,drop=FALSE]
  rownames(m) <- m$sample_id
  m <- m[rownames(counts),,drop=FALSE]
  m$management <- ifelse(m$management_group=="till","vtill","NT")
  m$treatment <- "MF"
  m$depth <- ifelse(m$Depth_cm=="0-10",10,30)
  m$C <- m[[paste0("C_",marker)]]
  ps <- phyloseq(otu_table(counts,taxa_are_rows=FALSE),tax_table(taxonomy),sample_data(m))
  saveRDS(ps,files[[marker]])
}
