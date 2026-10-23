# Optional focused rerun; NOT executed in this release.
# Original archived subsets were not saved as sample-ID manifests. This rerun
# records its selections and does NOT claim bitwise reproduction of historical draws.
# Run from package root: Rscript source/run_balanced_sensitivity.R
for (pkg in c("vegan","Hmisc","igraph","tibble")) stopifnot(requireNamespace(pkg,quietly=TRUE))
library(igraph)
library(tibble)
clr_transform <- function(mat){
  mat <- as.matrix(mat) + 1
  gm <- exp(rowMeans(log(mat)))
  log(mat / gm)
}

build_corr_network <- function(mat, cor_cut = 0.6, q_cut = 0.05, max_taxa = 250){
  prev <- colMeans(mat > 0)
  keep <- which(prev >= 0.1)
  if (length(keep) < 10) keep <- order(colSums(mat), decreasing = TRUE)[seq_len(min(100, ncol(mat)))]
  mat2 <- mat[, keep, drop = FALSE]
  if (ncol(mat2) > max_taxa) {
    top_idx <- order(colSums(mat2), decreasing = TRUE)[seq_len(max_taxa)]
    mat2 <- mat2[, top_idx, drop = FALSE]
  }
  clr <- clr_transform(mat2)
  rc <- Hmisc::rcorr(as.matrix(clr), type = "spearman")
  R <- rc$r
  P <- rc$P
  ut <- which(upper.tri(R), arr.ind = TRUE)
  qvals <- p.adjust(P[ut], method = "BH")
  keep_e <- which(abs(R[ut]) >= cor_cut & qvals < q_cut)
  if (length(keep_e) == 0) return(make_empty_graph())
  ed <- tibble(
    from = colnames(R)[ut[keep_e, 1]],
    to = colnames(R)[ut[keep_e, 2]],
    weight = R[ut[keep_e]]
  )
  graph_from_data_frame(ed, directed = FALSE)
}


set.seed(20260904)
out <- "reproduced/balanced_sensitivity"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
net_rows <- list(); beta_rows <- list(); selected <- list()
for (marker in c("16S","ITS","AMF")) {
 x <- as.matrix(read.csv(gzfile(file.path("data",paste0(marker,"_counts.csv.gz"))),row.names=1,check.names=FALSE))
 m <- read.csv(file.path("data",paste0(marker,"_samples.csv")),row.names=2)
 m <- m[rownames(x),,drop=FALSE]
 nt <- rownames(m)[m$management_group=="NT"]; ti <- rownames(m)[m$management_group=="till"]
 for (iteration in 1:20) {
  ids <- c(sample(nt,length(ti)),ti)
  selected[[length(selected)+1]] <- data.frame(marker,iteration,sample_id=ids,management=m[ids,"management_group"])
  z <- x[ids,,drop=FALSE]; group <- factor(m[ids,"management_group"],levels=c("NT","till"))
  d <- vegan::vegdist(z/rowSums(z),method="bray")
  a <- vegan::adonis2(d~group,permutations=999)
  disp <- vegan::permutest(vegan::betadisper(d,group),permutations=999)
  beta_rows[[length(beta_rows)+1]] <- data.frame(marker,iteration,R2=a$R2[1],p=a$`Pr(>F)`[1],dispersion_p=disp$tab$`Pr(>F)`[1])
  for (g in c("NT","till")) {
   graph <- build_corr_network(z[group==g,,drop=FALSE]);n <- vcount(graph);e <- ecount(graph)
   net_rows[[length(net_rows)+1]] <- data.frame(marker,iteration,management=g,nodes=n,edges=e,density=if(n>1) 2*e/n/(n-1) else NA,mean_degree=if(n>0) 2*e/n else NA)
  }
 }
}
write.csv(do.call(rbind,selected),file.path(out,"selected_samples.csv"),row.names=FALSE)
write.csv(do.call(rbind,beta_rows),file.path(out,"beta.csv"),row.names=FALSE)
write.csv(do.call(rbind,net_rows),file.path(out,"networks.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
