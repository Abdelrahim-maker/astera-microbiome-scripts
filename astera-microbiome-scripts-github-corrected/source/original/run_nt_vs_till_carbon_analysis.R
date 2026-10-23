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
  "phyloseq", "vegan", "ggplot2", "dplyr", "tidyr", "tibble", "stringr",
  "readr", "purrr", "igraph", "Hmisc"
)
for (p in cran_pkgs) install_if_missing(p)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("Maaslin2", quietly = TRUE)) {
  BiocManager::install("Maaslin2", ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(purrr)
  library(igraph)
  library(Hmisc)
})

set.seed(20260904)

root_dir <- "."
out_dir <- file.path(root_dir, "analysis_outputs", "carbon_analysis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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

network_stats <- function(g){
  if (vcount(g) == 0 || ecount(g) == 0) {
    return(tibble(nodes = vcount(g), edges = ecount(g), density = NA_real_, mean_degree = NA_real_, transitivity = NA_real_, modularity = NA_real_))
  }
  cl <- cluster_louvain(g, weights = abs(E(g)$weight))
  tibble(
    nodes = vcount(g),
    edges = ecount(g),
    density = edge_density(g, loops = FALSE),
    mean_degree = mean(degree(g)),
    transitivity = transitivity(g, type = "globalundirected"),
    modularity = modularity(cl)
  )
}

run_maaslin_carbon <- function(ps, domain, out_dir_domain){
  mat <- otu_samples_by_taxa(ps)
  mat <- mat[, colSums(mat > 0) >= max(3, ceiling(0.1 * nrow(mat))), drop = FALSE]
  if (ncol(mat) < 5) return(tibble())

  rel <- decostand(mat, method = "total")
  feat <- as.data.frame(rel)
  meta <- as(sample_data(ps), "data.frame") %>%
    mutate(
      trt2 = factor(trt2, levels = c("NT", "till")),
      C = suppressWarnings(as.numeric(C))
    ) %>%
    filter(!is.na(C))

  feat <- feat[rownames(meta), , drop = FALSE]
  if (nrow(meta) < 10) return(tibble())

  meta$C_by_till <- ifelse(meta$trt2 == "till", meta$C, 0)

  maas_out <- file.path(out_dir_domain, paste0(domain, "_maaslin2_carbon"))
  dir.create(maas_out, showWarnings = FALSE, recursive = TRUE)

  suppressMessages(
    Maaslin2::Maaslin2(
      input_data = feat,
      input_metadata = meta,
      output = maas_out,
      fixed_effects = c("trt2", "C", "C_by_till"),
      normalization = "NONE",
      transform = "LOG",
      analysis_method = "LM",
      correction = "BH",
      max_significance = 0.25,
      plot_heatmap = FALSE,
      plot_scatter = FALSE
    )
  )

  res_file <- file.path(maas_out, "all_results.tsv")
  if (!file.exists(res_file)) return(tibble())

  readr::read_tsv(res_file, show_col_types = FALSE) %>%
    mutate(domain = domain)
}

all_ps <- lapply(rds_files, readRDS)
mf_ps <- imap(all_ps, ~ subset_mf_nt_till(.x))

carbon_meta_summary <- list()
carbon_treatment_tests <- list()
carbon_permanova <- list()
carbon_varpart <- list()
carbon_maaslin_all <- list()
network_boot <- list()
hub_carbon_assoc <- list()

hub_tbl <- read_csv(file.path("analysis_outputs", "hub_nodes_taxonomy_maaslin.csv"), show_col_types = FALSE)

for (domain in names(mf_ps)) {
  domain_label <- domain
  ps <- mf_ps[[domain]]
  domain_dir <- file.path(out_dir, domain)
  dir.create(domain_dir, showWarnings = FALSE, recursive = TRUE)

  meta <- as(sample_data(ps), "data.frame") %>%
    rownames_to_column("sample_key") %>%
    mutate(
      trt2 = factor(trt2, levels = c("NT", "till")),
      C = suppressWarnings(as.numeric(C))
    )

  carbon_meta_summary[[domain]] <- meta %>%
    group_by(trt2) %>%
    summarise(
      domain = domain,
      n = n(),
      n_non_na_c = sum(!is.na(C)),
      mean_c = mean(C, na.rm = TRUE),
      sd_c = sd(C, na.rm = TRUE),
      median_c = median(C, na.rm = TRUE),
      .groups = "drop"
    )

  c_dat <- meta %>% filter(!is.na(C))
  if (n_distinct(c_dat$trt2) == 2) {
    wt <- wilcox.test(C ~ trt2, data = c_dat, exact = FALSE)
    carbon_treatment_tests[[domain]] <- tibble(
      domain = domain,
      p_value = wt$p.value,
      nt_median = median(c_dat$C[c_dat$trt2 == "NT"], na.rm = TRUE),
      till_median = median(c_dat$C[c_dat$trt2 == "till"], na.rm = TRUE)
    )

    p_c <- ggplot(c_dat, aes(x = trt2, y = C, fill = trt2)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.15, alpha = 0.75, size = 1.8) +
      scale_fill_manual(values = c("NT" = "#1b9e77", "till" = "#d95f02")) +
      labs(
        title = paste0(domain, " MF: Soil carbon by treatment"),
        subtitle = paste0("Wilcoxon p = ", format(wt$p.value, digits = 3, scientific = TRUE)),
        x = NULL,
        y = "Carbon (C)"
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "none")
    ggsave(file.path(domain_dir, paste0(domain, "_carbon_by_treatment.png")), p_c, width = 6, height = 4.3, dpi = 320)
  }

  # PERMANOVA with treatment + carbon + interaction
  mat <- otu_samples_by_taxa(ps)
  rel <- decostand(mat, method = "total")
  rel <- rel[match(c_dat$sample_key, rownames(rel)), , drop = FALSE]
  d <- vegdist(rel, method = "bray")

  ad_seq <- adonis2(d ~ trt2 + C + trt2:C, data = c_dat, permutations = 999, by = "terms")
  ad_mar <- adonis2(d ~ trt2 + C + trt2:C, data = c_dat, permutations = 999, by = "margin")

  seq_tbl <- as.data.frame(ad_seq) %>% rownames_to_column("term") %>% mutate(domain = domain, test_type = "sequential")
  mar_tbl <- as.data.frame(ad_mar) %>% rownames_to_column("term") %>% mutate(domain = domain, test_type = "marginal")
  carbon_permanova[[domain]] <- bind_rows(seq_tbl, mar_tbl)

  # Variance partitioning on Hellinger-transformed data
  X <- model.matrix(~ trt2, data = c_dat)[, -1, drop = FALSE]
  Z <- as.matrix(c_dat[, "C", drop = FALSE])
  Y <- decostand(rel, method = "hellinger")
  vp <- varpart(Y, X, Z)
  fract <- as.data.frame(vp$part$indfract)
  carbon_varpart[[domain]] <- tibble(
    domain = domain,
    fraction = rownames(fract),
    adj_r2 = fract$Adj.R.square
  )

  # MaAsLin2 with carbon and treatment-carbon slope difference proxy
  ma <- run_maaslin_carbon(ps, domain, domain_dir)
  if (nrow(ma) > 0) {
    ma2 <- ma %>%
      mutate(
        significance = qval <= 0.25,
        effect_type = case_when(
          metadata == "trt2" & value == "till" ~ "treatment_main",
          metadata == "C" ~ "carbon_main_NT_slope",
          metadata == "C_by_till" ~ "carbon_slope_difference_till_vs_NT",
          TRUE ~ "other"
        )
      )
    write_csv(ma2, file.path(domain_dir, paste0(domain, "_maaslin2_treatment_carbon_results.csv")))
    carbon_maaslin_all[[domain]] <- ma2
  }

  # Hub/enriched ASV carbon-abundance associations by treatment
  mat_counts <- otu_samples_by_taxa(ps)
  keep_taxa <- intersect(colnames(mat_counts), hub_tbl$asv[hub_tbl$domain == domain_label])
  if (length(keep_taxa) > 0) {
    clr <- clr_transform(mat_counts[, keep_taxa, drop = FALSE])
    clr_df <- as.data.frame(clr) %>% rownames_to_column("sample_key") %>%
      pivot_longer(-sample_key, names_to = "asv", values_to = "clr_abund") %>%
      left_join(meta %>% select(sample_key, trt2, C), by = "sample_key") %>%
      left_join(hub_tbl %>% filter(.data$domain == domain_label) %>%
                  mutate(network_treatment = ifelse(grepl("_NT$", network), "NT", ifelse(grepl("_till$", network), "till", NA_character_))) %>%
                  select(asv, network, degree, hub, maaslin_enriched, taxonomy, network_treatment),
                by = "asv") %>%
      filter(!is.na(C), !is.na(network_treatment), trt2 == network_treatment)

    if (nrow(clr_df) > 0) {
      corr_tbl <- clr_df %>%
        group_by(network, trt2, asv, degree, hub, maaslin_enriched, taxonomy) %>%
        summarise(
          n = n(),
          cor_pearson = ifelse(n() >= 5, cor(C, clr_abund, use = "complete.obs", method = "pearson"), NA_real_),
          p_pearson = ifelse(n() >= 5, suppressWarnings(cor.test(C, clr_abund, method = "pearson")$p.value), NA_real_),
          .groups = "drop"
        ) %>%
        mutate(domain = domain_label) %>%
        mutate(q_pearson = p.adjust(p_pearson, method = "BH"))

      hub_carbon_assoc[[domain]] <- corr_tbl
    }
  }

  # Bootstrap network-carbon relationships within each treatment
  n_boot <- 80L
  for (trt in c("NT", "till")) {
    ids <- meta$sample_key[meta$trt2 == trt & !is.na(meta$C)]
    if (length(ids) < 8) next
    boot_n <- min(10L, length(ids))

    for (b in seq_len(n_boot)) {
      keep_ids <- sample(ids, size = boot_n, replace = TRUE)
      mat_sub <- mat_counts[match(keep_ids, rownames(mat_counts)), , drop = FALSE]
      g <- build_corr_network(mat_sub)
      st <- network_stats(g)
      network_boot[[paste(domain, trt, b, sep = "_")]] <- st %>%
        mutate(
          domain = domain,
          trt2 = trt,
          boot = b,
          mean_C = mean(meta$C[match(keep_ids, meta$sample_key)], na.rm = TRUE)
        )
    }
  }
}

write_csv(bind_rows(carbon_meta_summary), file.path(out_dir, "carbon_metadata_summary.csv"))
write_csv(bind_rows(carbon_treatment_tests), file.path(out_dir, "carbon_treatment_wilcox.csv"))
write_csv(bind_rows(carbon_permanova), file.path(out_dir, "carbon_permanova_by_domain.csv"))
write_csv(bind_rows(carbon_varpart), file.path(out_dir, "carbon_variance_partition_by_domain.csv"))

ma_all <- bind_rows(carbon_maaslin_all)
if (nrow(ma_all) > 0) {
  write_csv(ma_all, file.path(out_dir, "carbon_maaslin_all_domains.csv"))
  ma_sum <- ma_all %>%
    filter(effect_type != "other") %>%
    group_by(domain, effect_type) %>%
    summarise(
      n_total = n(),
      n_sig_q025 = sum(qval <= 0.25, na.rm = TRUE),
      .groups = "drop"
    )
  write_csv(ma_sum, file.path(out_dir, "carbon_maaslin_effect_summary.csv"))
}

hub_assoc <- bind_rows(hub_carbon_assoc)
if (nrow(hub_assoc) > 0) {
  write_csv(hub_assoc, file.path(out_dir, "hub_enriched_asv_carbon_correlations.csv"))
  top_assoc <- hub_assoc %>%
    filter(!is.na(q_pearson)) %>%
    arrange(q_pearson) %>%
    slice_head(n = 24)
  write_csv(top_assoc, file.path(out_dir, "hub_enriched_asv_carbon_correlations_top.csv"))

  p_top <- top_assoc %>%
    mutate(label = paste0(domain, " ", asv, " (", trt2, ")")) %>%
    ggplot(aes(x = reorder(label, cor_pearson), y = cor_pearson, fill = trt2)) +
    geom_col(alpha = 0.85) +
    coord_flip() +
    scale_fill_manual(values = c("NT" = "#1b9e77", "till" = "#d95f02")) +
    labs(title = "Top hub/enriched ASV carbon correlations", x = NULL, y = "Pearson correlation (C vs CLR abundance)") +
    theme_bw(base_size = 10)
  ggsave(file.path(out_dir, "hub_enriched_top_carbon_correlations.png"), p_top, width = 9, height = 7, dpi = 320)
}

net_boot <- bind_rows(network_boot)
if (nrow(net_boot) > 0) {
  write_csv(net_boot, file.path(out_dir, "network_bootstrap_stats_with_carbon.csv"))

  net_assoc <- net_boot %>%
    group_by(domain, trt2) %>%
    summarise(
      n_boot = n(),
      cor_mean_degree = cor(mean_C, mean_degree, use = "complete.obs"),
      p_mean_degree = suppressWarnings(cor.test(mean_C, mean_degree)$p.value),
      cor_density = cor(mean_C, density, use = "complete.obs"),
      p_density = suppressWarnings(cor.test(mean_C, density)$p.value),
      cor_modularity = cor(mean_C, modularity, use = "complete.obs"),
      p_modularity = suppressWarnings(cor.test(mean_C, modularity)$p.value),
      .groups = "drop"
    ) %>%
    mutate(
      q_mean_degree = p.adjust(p_mean_degree, method = "BH"),
      q_density = p.adjust(p_density, method = "BH"),
      q_modularity = p.adjust(p_modularity, method = "BH")
    )
  write_csv(net_assoc, file.path(out_dir, "network_metric_carbon_associations.csv"))

  p_net <- ggplot(net_boot, aes(x = mean_C, y = mean_degree, color = trt2)) +
    geom_point(alpha = 0.6, size = 1.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    facet_wrap(~ domain, scales = "free_y") +
    scale_color_manual(values = c("NT" = "#1b9e77", "till" = "#d95f02")) +
    labs(title = "Bootstrap network mean degree vs carbon", x = "Mean carbon of sampled subset", y = "Network mean degree") +
    theme_bw(base_size = 11)
  ggsave(file.path(out_dir, "network_mean_degree_vs_carbon_bootstrap.png"), p_net, width = 9, height = 4.5, dpi = 320)
}

# Joint integrated PERMANOVA with treatment + carbon + interaction
common_ids <- Reduce(intersect, lapply(mf_ps, sample_names))
if (length(common_ids) >= 10) {
  joint_list <- list()
  for (domain in names(mf_ps)) {
    m <- otu_samples_by_taxa(prune_samples(common_ids, mf_ps[[domain]]))
    prev <- colMeans(m > 0)
    m <- m[, prev >= 0.05, drop = FALSE]
    if (ncol(m) > 200) m <- m[, order(colSums(m), decreasing = TRUE)[1:200], drop = FALSE]
    clr <- clr_transform(m)
    colnames(clr) <- paste0(domain, "|", colnames(clr))
    joint_list[[domain]] <- clr
  }
  joint_mat <- do.call(cbind, joint_list)
  joint_meta <- as(sample_data(prune_samples(common_ids, mf_ps[[1]])), "data.frame") %>%
    rownames_to_column("sample_key") %>%
    mutate(C = suppressWarnings(as.numeric(C))) %>%
    filter(!is.na(C))

  joint_mat <- joint_mat[match(joint_meta$sample_key, rownames(joint_mat)), , drop = FALSE]
  d_joint <- dist(scale(joint_mat))

  j_seq <- adonis2(d_joint ~ trt2 + C + trt2:C, data = joint_meta, permutations = 999, by = "terms")
  j_mar <- adonis2(d_joint ~ trt2 + C + trt2:C, data = joint_meta, permutations = 999, by = "margin")

  joint_tbl <- bind_rows(
    as.data.frame(j_seq) %>% rownames_to_column("term") %>% mutate(test_type = "sequential"),
    as.data.frame(j_mar) %>% rownames_to_column("term") %>% mutate(test_type = "marginal")
  )
  write_csv(joint_tbl, file.path(out_dir, "carbon_joint_permanova.csv"))
}

# Compact text summary for manuscript insertion
sum_lines <- c("# Carbon Analysis Summary", "")
if (file.exists(file.path(out_dir, "carbon_treatment_wilcox.csv"))) {
  cw <- read_csv(file.path(out_dir, "carbon_treatment_wilcox.csv"), show_col_types = FALSE)
  sum_lines <- c(sum_lines, "## Carbon Differences Between NT and Till")
  for (i in seq_len(nrow(cw))) {
    r <- cw[i, ]
    sum_lines <- c(sum_lines, sprintf("- %s: NT median C = %.3f, till median C = %.3f (Wilcoxon p = %s)", r$domain, r$nt_median, r$till_median, format(r$p_value, digits = 3, scientific = TRUE)))
  }
  sum_lines <- c(sum_lines, "")
}

if (file.exists(file.path(out_dir, "carbon_permanova_by_domain.csv"))) {
  cp <- read_csv(file.path(out_dir, "carbon_permanova_by_domain.csv"), show_col_types = FALSE) %>%
    filter(test_type == "sequential", term %in% c("trt2", "C", "trt2:C"))
  sum_lines <- c(sum_lines, "## Domain PERMANOVA (sequential terms)")
  for (i in seq_len(nrow(cp))) {
    r <- cp[i, ]
    sum_lines <- c(sum_lines, sprintf("- %s %s: R2 = %.4f, p = %s", r$domain, r$term, r$R2, format(r$`Pr(>F)`, digits = 3, scientific = TRUE)))
  }
  sum_lines <- c(sum_lines, "")
}

if (file.exists(file.path(out_dir, "carbon_maaslin_effect_summary.csv"))) {
  cm <- read_csv(file.path(out_dir, "carbon_maaslin_effect_summary.csv"), show_col_types = FALSE)
  sum_lines <- c(sum_lines, "## MaAsLin2 Carbon-Aware Effects (q <= 0.25)")
  for (i in seq_len(nrow(cm))) {
    r <- cm[i, ]
    sum_lines <- c(sum_lines, sprintf("- %s %s: %d/%d significant", r$domain, r$effect_type, r$n_sig_q025, r$n_total))
  }
  sum_lines <- c(sum_lines, "")
}

if (file.exists(file.path(out_dir, "network_metric_carbon_associations.csv"))) {
  cn <- read_csv(file.path(out_dir, "network_metric_carbon_associations.csv"), show_col_types = FALSE)
  sum_lines <- c(sum_lines, "## Network Metric Associations with Carbon")
  for (i in seq_len(nrow(cn))) {
    r <- cn[i, ]
    sum_lines <- c(sum_lines, sprintf("- %s %s: mean_degree r = %.3f (q = %s), density r = %.3f (q = %s), modularity r = %.3f (q = %s)", r$domain, r$trt2, r$cor_mean_degree, format(r$q_mean_degree, digits = 3, scientific = TRUE), r$cor_density, format(r$q_density, digits = 3, scientific = TRUE), r$cor_modularity, format(r$q_modularity, digits = 3, scientific = TRUE)))
  }
  sum_lines <- c(sum_lines, "")
}

sum_lines <- c(sum_lines, "## Output Files", "- Domain plots and tables are in per-domain subfolders.", "- Cross-domain summaries are in this folder.")
writeLines(sum_lines, file.path(out_dir, "carbon_analysis_summary.md"))

message("Carbon analysis complete. Outputs written to: ", out_dir)
