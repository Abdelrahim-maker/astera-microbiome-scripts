#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  options(stringsAsFactors = FALSE)
})

install_if_missing <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

cran_pkgs <- c(
  "phyloseq", "vegan", "dplyr", "tibble", "stringr", "readr", "purrr",
  "tidyr", "igraph", "Hmisc"
)
for (pkg in cran_pkgs) install_if_missing(pkg)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("Maaslin2", quietly = TRUE)) {
  BiocManager::install("Maaslin2", ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(purrr)
  library(tidyr)
  library(igraph)
  library(Hmisc)
})

root_dir <- "."
out_dir <- file.path(root_dir, "analysis_outputs", "imbalance_sensitivity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_iter <- 20L
set.seed(20260904)

rds_files <- c(
  `16S` = "combined_16S_ps.rds",
  ITS = "combined_its.rds",
  AMF = "combined_amf.rds"
)

find_group_levels <- function(meta){
  lower_vals <- lapply(meta, function(x) tolower(as.character(x)))

  mf_score <- vapply(lower_vals, function(v) {
    sum(v %in% c("mf") | str_detect(v, "^mf$|manure|m[-_ ]?fert"), na.rm = TRUE)
  }, numeric(1))
  ntvt_score <- vapply(lower_vals, function(v) {
    sum(str_detect(v, "^nt$|no[-_ ]?till|notill|v[-_ ]?till|vtill|vertical|till"), na.rm = TRUE)
  }, numeric(1))

  soil_col <- names(which.max(mf_score))
  trt_col <- names(which.max(ntvt_score))

  if (mf_score[[soil_col]] == 0 || ntvt_score[[trt_col]] == 0) {
    stop("Could not auto-detect MF and tillage columns")
  }

  soil_vals <- tolower(as.character(meta[[soil_col]]))
  trt_vals <- tolower(as.character(meta[[trt_col]]))
  mf_idx <- which(soil_vals == "mf" | str_detect(soil_vals, "^mf$|manure|m[-_ ]?fert"))
  trt_mf <- trt_vals[mf_idx]

  nt_candidates <- unique(trt_mf[str_detect(trt_mf, "^nt$|no[-_ ]?till|notill")])
  till_candidates <- unique(trt_mf[str_detect(trt_mf, "^till$|v[-_ ]?till|vtill|vertical")])
  if (length(nt_candidates) == 0) nt_candidates <- unique(trt_mf[str_detect(trt_mf, "nt")])
  if (length(till_candidates) == 0) till_candidates <- unique(trt_mf[str_detect(trt_mf, "till") & !str_detect(trt_mf, "nt|no[-_ ]?till|notill")])

  list(mf_idx = mf_idx, nt_levels = nt_candidates, till_levels = till_candidates, soil_col = soil_col, trt_col = trt_col)
}

subset_mf_nt_till <- function(ps){
  meta <- as(sample_data(ps), "data.frame")
  gl <- find_group_levels(meta)

  meta$sample_key <- rownames(meta)
  meta$soil_std <- tolower(as.character(meta[[gl$soil_col]]))
  meta$trt_std <- tolower(as.character(meta[[gl$trt_col]]))

  keep <- meta %>%
    filter((soil_std == "mf" | str_detect(soil_std, "^mf$|manure|m[-_ ]?fert")) &
             (trt_std %in% c(gl$nt_levels, gl$till_levels))) %>%
    mutate(trt2 = ifelse(trt_std %in% gl$nt_levels, "NT", "till"))

  ps2 <- prune_samples(keep$sample_key, ps)
  sample_data(ps2)$trt2 <- factor(keep$trt2[match(sample_names(ps2), keep$sample_key)], levels = c("NT", "till"))
  ps2
}

otu_samples_by_taxa <- function(ps){
  m <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) m <- t(m)
  m
}

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

network_stats <- function(g, label){
  if (vcount(g) == 0 || ecount(g) == 0) {
    return(tibble(label = label, nodes = vcount(g), edges = ecount(g), density = NA_real_, mean_degree = NA_real_, transitivity = NA_real_, modularity = NA_real_))
  }
  cl <- cluster_louvain(g, weights = abs(E(g)$weight))
  tibble(
    label = label,
    nodes = vcount(g),
    edges = ecount(g),
    density = edge_density(g, loops = FALSE),
    mean_degree = mean(degree(g)),
    transitivity = transitivity(g, type = "globalundirected"),
    modularity = modularity(cl)
  )
}

alpha_tests <- function(ps, domain, subset_label){
  adiv <- estimate_richness(ps, measures = c("Observed", "Shannon")) %>%
    rownames_to_column("sample_key")
  meta <- as(sample_data(ps), "data.frame") %>%
    rownames_to_column("sample_key")
  dat <- left_join(adiv, meta, by = "sample_key")
  obs <- suppressWarnings(wilcox.test(Observed ~ trt2, data = dat, exact = FALSE))
  sha <- suppressWarnings(wilcox.test(Shannon ~ trt2, data = dat, exact = FALSE))
  tibble(
    domain = domain,
    subset = subset_label,
    metric = c("Observed", "Shannon"),
    p_value = c(obs$p.value, sha$p.value),
    nt_median = c(median(dat$Observed[dat$trt2 == "NT"], na.rm = TRUE), median(dat$Shannon[dat$trt2 == "NT"], na.rm = TRUE)),
    till_median = c(median(dat$Observed[dat$trt2 == "till"], na.rm = TRUE), median(dat$Shannon[dat$trt2 == "till"], na.rm = TRUE))
  )
}

ordination_tests <- function(ps, domain, subset_label){
  mat <- otu_samples_by_taxa(ps)
  rel <- decostand(mat, method = "total")
  dist_bray <- vegdist(rel, method = "bray")
  meta <- as(sample_data(ps), "data.frame")

  ad <- adonis2(dist_bray ~ trt2, data = meta, permutations = 999)
  ad_df <- as.data.frame(ad)
  ad_df$term <- rownames(ad_df)
  model_row <- ad_df[ad_df$term == "Model", , drop = FALSE]

  bd <- betadisper(dist_bray, meta$trt2)
  pd <- permutest(bd, permutations = 999)

  tibble(
    domain = domain,
    subset = subset_label,
    permanova_r2 = model_row$R2[1],
    permanova_p = model_row$`Pr(>F)`[1],
    permdisp_f = unname(pd$tab[1, "F"]),
    permdisp_p = unname(pd$tab[1, "Pr(>F)"])
  )
}

run_maaslin <- function(ps, domain, subset_label, iter_id = NA_integer_){
  if (!requireNamespace("Maaslin2", quietly = TRUE)) {
    return(list(summary = tibble(domain = domain, subset = subset_label, iter = iter_id, n_sig = NA_integer_), features = tibble()))
  }

  mat <- otu_samples_by_taxa(ps)
  mat <- mat[, colSums(mat > 0) >= max(3, ceiling(0.1 * nrow(mat))), drop = FALSE]
  if (ncol(mat) < 5) {
    return(list(summary = tibble(domain = domain, subset = subset_label, iter = iter_id, n_sig = 0L), features = tibble()))
  }

  rel <- decostand(mat, method = "total")
  feat <- as.data.frame(rel)
  meta <- as(sample_data(ps), "data.frame")
  meta$trt2 <- factor(meta$trt2, levels = c("NT", "till"))
  out_tmp <- tempfile(pattern = paste0("maaslin_", domain, "_"))
  dir.create(out_tmp, showWarnings = FALSE, recursive = TRUE)

  suppressMessages(
    Maaslin2::Maaslin2(
      input_data = feat,
      input_metadata = meta,
      output = out_tmp,
      fixed_effects = c("trt2"),
      normalization = "NONE",
      transform = "LOG",
      analysis_method = "LM",
      correction = "BH",
      max_significance = 0.25,
      plot_heatmap = FALSE,
      plot_scatter = FALSE
    )
  )

  res_file <- file.path(out_tmp, "all_results.tsv")
  if (!file.exists(res_file)) {
    unlink(out_tmp, recursive = TRUE, force = TRUE)
    return(list(summary = tibble(domain = domain, subset = subset_label, iter = iter_id, n_sig = NA_integer_), features = tibble()))
  }

  res <- read_tsv(res_file, show_col_types = FALSE) %>%
    filter(metadata == "trt2", value == "till")
  sig <- res %>% filter(qval <= 0.25) %>% select(feature, coef, qval)
  unlink(out_tmp, recursive = TRUE, force = TRUE)

  list(
    summary = tibble(domain = domain, subset = subset_label, iter = iter_id, n_sig = nrow(sig)),
    features = mutate(sig, domain = domain, subset = subset_label, iter = iter_id)
  )
}

network_compare <- function(ps, domain, subset_label, iter_id = NA_integer_){
  nt_ids <- sample_names(ps)[sample_data(ps)$trt2 == "NT"]
  till_ids <- sample_names(ps)[sample_data(ps)$trt2 == "till"]
  nt_mat <- otu_samples_by_taxa(prune_samples(nt_ids, ps))
  till_mat <- otu_samples_by_taxa(prune_samples(till_ids, ps))

  g_nt <- build_corr_network(nt_mat)
  g_till <- build_corr_network(till_mat)
  st_nt <- network_stats(g_nt, "NT")
  st_till <- network_stats(g_till, "till")

  deg_w <- if (vcount(g_nt) > 2 && vcount(g_till) > 2 && ecount(g_nt) > 0 && ecount(g_till) > 0) suppressWarnings(wilcox.test(degree(g_nt), degree(g_till), exact = FALSE)$p.value) else NA_real_
  deg_ks <- if (vcount(g_nt) > 2 && vcount(g_till) > 2 && ecount(g_nt) > 0 && ecount(g_till) > 0) suppressWarnings(ks.test(degree(g_nt), degree(g_till))$p.value) else NA_real_

  bind_rows(st_nt, st_till) %>%
    mutate(domain = domain, subset = subset_label, iter = iter_id, degree_wilcox_p = deg_w, degree_ks_p = deg_ks)
}

joint_permanova <- function(ps_list, subset_label, iter_id = NA_integer_){
  common_ids <- Reduce(intersect, lapply(ps_list, sample_names))
  if (length(common_ids) < 6) {
    return(tibble(subset = subset_label, iter = iter_id, r2 = NA_real_, p_value = NA_real_))
  }
  joint_list <- list()
  for (domain in names(ps_list)) {
    m <- otu_samples_by_taxa(prune_samples(common_ids, ps_list[[domain]]))
    prev <- colMeans(m > 0)
    m <- m[, prev >= 0.05, drop = FALSE]
    if (ncol(m) > 200) {
      m <- m[, order(colSums(m), decreasing = TRUE)[1:200], drop = FALSE]
    }
    clr <- clr_transform(m)
    colnames(clr) <- paste0(domain, "|", colnames(clr))
    joint_list[[domain]] <- clr
  }
  joint_mat <- do.call(cbind, joint_list)
  meta <- as(sample_data(prune_samples(common_ids, ps_list[[1]])), "data.frame")
  d_joint <- dist(scale(joint_mat))
  ad <- adonis2(d_joint ~ trt2, data = meta, permutations = 999)
  ad_df <- as.data.frame(ad)
  ad_df$term <- rownames(ad_df)
  model_row <- ad_df[ad_df$term == "Model", , drop = FALSE]
  tibble(subset = subset_label, iter = iter_id, r2 = model_row$R2[1], p_value = model_row$`Pr(>F)`[1])
}

all_ps <- lapply(rds_files, readRDS)
mf_ps <- imap(all_ps, ~subset_mf_nt_till(.x))

alpha_obs <- list()
ord_obs <- list()
maaslin_obs <- list()
maaslin_obs_features <- list()
network_obs <- list()

for (domain in names(mf_ps)) {
  ps <- mf_ps[[domain]]
  alpha_obs[[domain]] <- alpha_tests(ps, domain, "observed")
  ord_obs[[domain]] <- ordination_tests(ps, domain, "observed")
  ma <- run_maaslin(ps, domain, "observed", 0L)
  maaslin_obs[[domain]] <- ma$summary
  maaslin_obs_features[[domain]] <- ma$features
  network_obs[[domain]] <- network_compare(ps, domain, "observed", 0L)
}

write_csv(bind_rows(alpha_obs), file.path(out_dir, "alpha_observed.csv"))
write_csv(bind_rows(ord_obs), file.path(out_dir, "ordination_observed_permanova_permdisp.csv"))
write_csv(bind_rows(maaslin_obs), file.path(out_dir, "maaslin_observed_summary.csv"))
write_csv(bind_rows(maaslin_obs_features), file.path(out_dir, "maaslin_observed_significant_features.csv"))
write_csv(bind_rows(network_obs), file.path(out_dir, "network_observed_stats.csv"))

joint_obs <- joint_permanova(mf_ps, "observed", 0L)
write_csv(joint_obs, file.path(out_dir, "joint_observed_permanova.csv"))

joint_common_ids <- Reduce(intersect, lapply(mf_ps, sample_names))
joint_meta <- as(sample_data(prune_samples(joint_common_ids, mf_ps[[1]])), "data.frame")
joint_nt_ids <- rownames(joint_meta)[joint_meta$trt2 == "NT"]
joint_till_ids <- rownames(joint_meta)[joint_meta$trt2 == "till"]
joint_n_each <- min(length(joint_nt_ids), length(joint_till_ids))

alpha_bal <- list()
ord_bal <- list()
maas_bal <- list()
maas_bal_features <- list()
network_bal <- list()
joint_bal <- list()

for (iter_id in seq_len(n_iter)) {
  subs <- list()
  for (domain in names(mf_ps)) {
    ps <- mf_ps[[domain]]
    nt_ids <- sample_names(ps)[sample_data(ps)$trt2 == "NT"]
    till_ids <- sample_names(ps)[sample_data(ps)$trt2 == "till"]
    n_each <- min(length(nt_ids), length(till_ids))
    keep_ids <- c(sample(nt_ids, n_each), sample(till_ids, n_each))
    subs[[domain]] <- prune_samples(keep_ids, ps)

    alpha_bal[[paste(domain, iter_id, sep = "_")]] <- mutate(alpha_tests(subs[[domain]], domain, "balanced"), iter = iter_id)
    ord_bal[[paste(domain, iter_id, sep = "_")]] <- mutate(ordination_tests(subs[[domain]], domain, "balanced"), iter = iter_id)
    ma <- run_maaslin(subs[[domain]], domain, "balanced", iter_id)
    maas_bal[[paste(domain, iter_id, sep = "_")]] <- ma$summary
    maas_bal_features[[paste(domain, iter_id, sep = "_")]] <- ma$features
    network_bal[[paste(domain, iter_id, sep = "_")]] <- network_compare(subs[[domain]], domain, "balanced", iter_id)
  }

  if (joint_n_each >= 3) {
    joint_keep <- c(sample(joint_nt_ids, joint_n_each), sample(joint_till_ids, joint_n_each))
    subs_common <- lapply(mf_ps, function(ps) prune_samples(joint_keep, ps))
    joint_bal[[as.character(iter_id)]] <- joint_permanova(subs_common, "balanced", iter_id)
  }
}

alpha_bal_df <- bind_rows(alpha_bal)
ord_bal_df <- bind_rows(ord_bal)
maas_bal_df <- bind_rows(maas_bal)
maas_bal_features_df <- bind_rows(maas_bal_features)
network_bal_df <- bind_rows(network_bal)
joint_bal_df <- bind_rows(joint_bal)

write_csv(alpha_bal_df, file.path(out_dir, "alpha_balanced_subsamples.csv"))
write_csv(ord_bal_df, file.path(out_dir, "ordination_balanced_subsamples_permanova_permdisp.csv"))
write_csv(maas_bal_df, file.path(out_dir, "maaslin_balanced_subsamples_summary.csv"))
write_csv(maas_bal_features_df, file.path(out_dir, "maaslin_balanced_subsamples_significant_features.csv"))
write_csv(network_bal_df, file.path(out_dir, "network_balanced_subsamples_stats.csv"))
write_csv(joint_bal_df, file.path(out_dir, "joint_balanced_subsamples_permanova.csv"))

alpha_balance_summary <- alpha_bal_df %>%
  group_by(domain, metric) %>%
  summarise(
    median_p = median(p_value, na.rm = TRUE),
    min_p = min(p_value, na.rm = TRUE),
    max_p = max(p_value, na.rm = TRUE),
    frac_p_lt_0_05 = mean(p_value < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

ord_balance_summary <- ord_bal_df %>%
  group_by(domain) %>%
  summarise(
    median_permanova_r2 = median(permanova_r2, na.rm = TRUE),
    median_permanova_p = median(permanova_p, na.rm = TRUE),
    frac_permanova_p_lt_0_05 = mean(permanova_p < 0.05, na.rm = TRUE),
    median_permdisp_p = median(permdisp_p, na.rm = TRUE),
    frac_permdisp_p_lt_0_05 = mean(permdisp_p < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

maaslin_balance_summary <- maas_bal_df %>%
  group_by(domain) %>%
  summarise(
    median_n_sig = median(n_sig, na.rm = TRUE),
    min_n_sig = min(n_sig, na.rm = TRUE),
    max_n_sig = max(n_sig, na.rm = TRUE),
    .groups = "drop"
  )

maaslin_feature_recurrence <- maas_bal_features_df %>%
  group_by(domain, feature) %>%
  summarise(
    times_detected = n(),
    freq_detected = n() / n_iter,
    median_coef = median(coef, na.rm = TRUE),
    best_q = min(qval, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(domain, desc(freq_detected), best_q)

network_balance_summary <- network_bal_df %>%
  select(domain, iter, label, nodes, edges, density, mean_degree, transitivity, modularity, degree_wilcox_p, degree_ks_p) %>%
  group_by(domain, label) %>%
  summarise(
    median_nodes = median(nodes, na.rm = TRUE),
    median_edges = median(edges, na.rm = TRUE),
    median_mean_degree = median(mean_degree, na.rm = TRUE),
    median_density = median(density, na.rm = TRUE),
    median_transitivity = median(transitivity, na.rm = TRUE),
    median_modularity = median(modularity, na.rm = TRUE),
    median_degree_wilcox_p = median(degree_wilcox_p, na.rm = TRUE),
    median_degree_ks_p = median(degree_ks_p, na.rm = TRUE),
    frac_degree_wilcox_p_lt_0_05 = mean(degree_wilcox_p < 0.05, na.rm = TRUE),
    frac_degree_ks_p_lt_0_05 = mean(degree_ks_p < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

joint_balance_summary <- joint_bal_df %>%
  summarise(
    median_r2 = median(r2, na.rm = TRUE),
    median_p = median(p_value, na.rm = TRUE),
    frac_p_lt_0_05 = mean(p_value < 0.05, na.rm = TRUE)
  )

write_csv(alpha_balance_summary, file.path(out_dir, "alpha_balanced_summary.csv"))
write_csv(ord_balance_summary, file.path(out_dir, "ordination_balanced_summary.csv"))
write_csv(maaslin_balance_summary, file.path(out_dir, "maaslin_balanced_summary.csv"))
write_csv(maaslin_feature_recurrence, file.path(out_dir, "maaslin_balanced_feature_recurrence.csv"))
write_csv(network_balance_summary, file.path(out_dir, "network_balanced_summary.csv"))
write_csv(joint_balance_summary, file.path(out_dir, "joint_balanced_summary.csv"))

message("Imbalance sensitivity analysis complete: ", normalizePath(out_dir))