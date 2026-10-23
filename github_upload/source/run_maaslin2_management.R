# Optional rerun adapter; not executed in this release (R/Maaslin2 unavailable).
# Run from package root: Rscript source/run_maaslin2_management.R
stopifnot(requireNamespace("Maaslin2", quietly=TRUE))
for (marker in c("16S","ITS","AMF")) {
 counts <- as.matrix(read.csv(gzfile(file.path("data",paste0(marker,"_counts.csv.gz"))),row.names=1,check.names=FALSE))
 meta <- read.csv(file.path("data",paste0(marker,"_samples.csv")),row.names=2,check.names=FALSE)
 meta <- meta[rownames(counts),,drop=FALSE]
 counts <- counts[,colSums(counts>0)>=max(3,ceiling(.1*nrow(counts))),drop=FALSE]
 rel <- counts/rowSums(counts)
 model_meta <- data.frame(trt2=factor(meta$management_group,levels=c("NT","till")),row.names=rownames(counts))
 Maaslin2::Maaslin2(input_data=as.data.frame(rel),input_metadata=model_meta,
  output=file.path("reproduced","maaslin2",marker),fixed_effects=c("trt2"),
  normalization="NONE",transform="LOG",analysis_method="LM",correction="BH",
  max_significance=.25,plot_heatmap=FALSE,plot_scatter=FALSE)
}
writeLines(capture.output(sessionInfo()),"reproduced/maaslin2/sessionInfo.txt")
