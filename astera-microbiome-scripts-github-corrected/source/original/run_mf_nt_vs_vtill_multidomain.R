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
  "readr", "purrr", "forcats", "igraph", "Hmisc", "RColorBrewer", "scales"
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
  library(forcats)
  library(igraph)
  library(Hmisc)
  library(RColorBrewer)
  library(scales)
})

root_dir <- "."
out_dir <- file.path(root_dir, "analysis_outputs")
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

  if (length(soil_col) == 0 || length(trt_col) == 0 || mf_score[[soil_col]] == 0 || ntvt_score[[trt_col]] == 0) {
    stop("Could not auto-detect MF and tillage columns")
  }

  soil_vals <- tolower(as.character(meta[[soil_col]]))
  trt_vals <- tolower(as.character(meta[[trt_col]]))

  mf_idx <- which(soil_vals == "mf" | str_detect(soil_vals, "^mf$|manure|m[-_ ]?fert"))
  if (length(mf_idx) == 0) stop("No MF samples found")

  trt_mf <- trt_vals[mf_idx]
  nt_candidates <- unique(trt_mf[str_detect(trt_mf, "^nt$|no[-_ ]?till|notill")])
  till_candidates <- unique(trt_mf[str_detect(trt_mf, "^till$|v[-_ ]?till|vtill|vertical")])

  if (length(nt_candidates) == 0) nt_candidates <- unique(trt_mf[str_detect(trt_mf, "nt")])
  if (length(till_candidates) == 0) till_candidates <- unique(trt_mf[str_detect(trt_mf, "till") & !str_detect(trt_mf, "nt|no[-_ ]?till|notill")])

  if (length(nt_candidates) == 0 || length(till_candidates) == 0) {
    stop("Could not detect NT and till levels within MF samples")
  }

  list(mf_idx = mf_idx, nt_levels = nt_candidates, till_levels = till_candidates, soil_col = soil_col, trt_col = trt_col)
}

subset_mf_nt_till <- function(ps, domain_name){
  meta <- as(sample_data(ps), "data.frame")
  gl <- find_group_levels(meta)

  meta$.__sample_id <- rownames(meta)
  meta$.__soil <- tolower(as.character(meta[[gl$soil_col]]))
  meta$.__trt <- tolower(as.character(meta[[gl$trt_col]]))

  keep <- meta %>%
    filter((.__soil == "mf" | str_detect(.__soil, "^mf$|manure|m[-_ ]?fert")) &
             (.__trt %in% c(gl$nt_levels, gl$till_levels))) %>%
    mutate(
      trt2 = case_when(
        .__trt %in% gl$nt_levels ~ "NT",
        .__trt %in% gl$till_levels ~ "till",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(trt2))

  ps2 <- prune_samples(keep$.__sample_id, ps)
  sample_data(ps2)$trt2 <- factor(keep$trt2[match(sample_names(ps2), keep$.__sample_id)], levels = c("NT", "till"))

  message(sprintf("%s: kept %d samples (NT=%d, till=%d)",
                  domain_name,
                  nsamples(ps2),
                  sum(sample_data(ps2)$trt2 == "NT", na.rm = TRUE),
                  sum(sample_data(ps2)$trt2 == "till", na.rm = TRUE)))
  message(sprintf("%s: soil column=%s, treatment column=%s", domain_name, gl$soil_col, gl$trt_col))

  ps2
}

otu_samples_by_taxa <- function(ps){
  m <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) m <- t(m)
  m
}

build_taxonomy_map <- function(ps, domain){
  tt <- tax_table(ps, errorIfNULL = FALSE)
  if (is.null(tt)) {
    return(tibble(domain = domain, asv = taxa_names(ps), taxonomy = NA_character_))
  }
  tx <- as(tt, "matrix")
  clean_tax <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x)]
    x <- x[x != ""]
    x <- x[!tolower(x) %in% c("na", "uncultured", "unclassified", "unknown", "unidentified")]
    if (length(x) == 0) return(NA_character_)
    paste(x, collapse = "; ")
  }
  tibble(
    domain = domain,
    asv = rownames(tx),
    taxonomy = apply(tx, 1, clean_tax)
  )
}

alpha_stats_and_plots <- function(ps, domain, out_subdir){
  adiv <- estimate_richness(ps, measures = c("Observed", "Shannon")) %>%
    rownames_to_column("sample_key")

  meta <- as(sample_data(ps), "data.frame") %>%
    rownames_to_column("sample_key") %>%
    mutate(trt2 = factor(trt2, levels = c("NT", "till")))

  dat <- adiv %>% left_join(meta, by = "sample_key")

  stat_obs <- wilcox.test(Observed ~ trt2, data = dat, exact = FALSE)
  stat_sha <- wilcox.test(Shannon ~ trt2, data = dat, exact = FALSE)

  stats_tbl <- tibble(
    domain = domain,
    metric = c("Observed", "Shannon"),
    p_value = c(stat_obs$p.value, stat_sha$p.value),
    nt_median = c(median(dat$Observed[dat$trt2 == "NT"], na.rm = TRUE), median(dat$Shannon[dat$trt2 == "NT"], na.rm = TRUE)),
    till_median = c(median(dat$Observed[dat$trt2 == "till"], na.rm = TRUE), median(dat$Shannon[dat$trt2 == "till"], na.rm = TRUE))
  )
  write_csv(stats_tbl, file.path(out_subdir, sprintf("%s_alpha_wilcox.csv", domain)))

  p_obs <- sprintf("Wilcoxon p = %.3g", stat_obs$p.value)
  p_sha <- sprintf("Wilcoxon p = %.3g", stat_sha$p.value)

  p1 <- ggplot(dat, aes(trt2, Observed, fill = trt2)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 1.8) +
    annotate("text", x = 1.5, y = max(dat$Observed, na.rm = TRUE) * 1.05, label = p_obs, size = 3.8) +
    scale_fill_manual(values = c("NT" = "#1b9e77", "till" = "#d95f02")) +
    labs(title = sprintf("%s MF: Observed Richness", domain), x = NULL, y = "Observed") +
    expand_limits(y = max(dat$Observed, na.rm = TRUE) * 1.12) +
    theme_bw(base_size = 12) + theme(legend.position = "none")

  p2 <- ggplot(dat, aes(trt2, Shannon, fill = trt2)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 1.8) +
    annotate("text", x = 1.5, y = max(dat$Shannon, na.rm = TRUE) * 1.05, label = p_sha, size = 3.8) +
    scale_fill_manual(values = c("NT" = "#1b9e77", "till" = "#d95f02")) +
    labs(title = sprintf("%s MF: Shannon", domain), x = NULL, y = "Shannon") +
    expand_limits(y = max(dat$Shannon, na.rm = TRUE) * 1.12) +
    theme_bw(base_size = 12) + theme(legend.position = "none")

  ggsave(file.path(out_subdir, sprintf("%s_observed_boxplot.png", domain)), p1, width = 6, height = 4, dpi = 300)
  ggsave(file.path(out_subdir, sprintf("%s_shannon_boxplot.png", domain)), p2, width = 6, height = 4, dpi = 300)
}

rarefaction_plot <- function(ps, domain, out_subdir){
  mat <- otu_samples_by_taxa(ps)
  meta <- as(sample_data(ps), "data.frame")
  cols <- ifelse(meta$trt2 == "NT", "#1b9e77", "#d95f02")

  png(file.path(out_subdir, sprintf("%s_rarefaction.png", domain)), width = 1700, height = 1200, res = 200)
  rarecurve(mat, step = 50, sample = min(rowSums(mat)), col = cols, label = FALSE, cex = 0.5,
            main = sprintf("%s MF rarefaction (NT vs till)", domain), xlab = "Reads", ylab = "Observed taxa")
  legend("bottomright", legend = c("NT", "till"), col = c("#1b9e77", "#d95f02"), lwd = 2, bty = "n")
  dev.off()
}

nmds_and_stats <- function(ps, domain, out_subdir){
  mat <- otu_samples_by_taxa(ps)
  rel <- decostand(mat, method = "total")
  dist_bray <- vegdist(rel, method = "bray")
  ord <- metaMDS(dist_bray, k = 2, trymax = 100, autotransform = FALSE, trace = FALSE)

  meta <- as(sample_data(ps), "data.frame") %>% rownames_to_column("sample_key")
  pts <- as.data.frame(ord$points) %>%
    rownames_to_column("sample_key") %>%
    rename(NMDS1 = MDS1, NMDS2 = MDS2) %>%
    left_join(meta, by = "sample_key")

  ad <- adonis2(dist_bray ~ trt2, data = meta, permutations = 999)
  ad_tbl <- as.data.frame(ad)
  ad_tbl$term <- rownames(ad_tbl)
  write_csv(ad_tbl, file.path(out_subdir, sprintf("%s_nmds_permanova.csv", domain)))

  model_row <- ad_tbl[ad_tbl$term == "Model", , drop = FALSE]
  subtitle_txt <- sprintf("Stress=%.3f; PERMANOVA R2=%.3f, p=%s", ord$stress, model_row$R2[1], format(model_row$`Pr(>F)`[1], digits = 3, scientific = TRUE))

  p <- ggplot(pts, aes(NMDS1, NMDS2, color = trt2)) +
    geom_point(size = 2.5, alpha = 0.85) +
    stat_ellipse(aes(fill = trt2), geom = "polygon", alpha = 0.15, color = NA) +
    scale_color_manual(values = c("NT" = "#1b9e77", "till" = "#d95f02")) +
    scale_fill_manual(values = c("NT" = "#1b9e77", "till" = "#d95f02")) +
    labs(
      title = sprintf("%s MF NMDS (Bray-Curtis)", domain),
      subtitle = subtitle_txt,
      color = "Treatment", fill = "Treatment"
    ) +
    theme_bw(base_size = 12)

  ggsave(file.path(out_subdir, sprintf("%s_nmds.png", domain)), p, width = 6.6, height = 5, dpi = 300)
}

clr_transform <- function(mat){
  mat <- as.matrix(mat)
  mat <- mat + 1
  gm <- exp(rowMeans(log(mat)))
  log(mat / gm)
}

run_maaslin_and_lollipop <- function(ps, domain, out_subdir, tax_map){
  if (!requireNamespace("Maaslin2", quietly = TRUE)) {
    message("Maaslin2 unavailable; skipping for ", domain)
    return(invisible(NULL))
  }

  mat <- otu_samples_by_taxa(ps)
  mat <- mat[, colSums(mat > 0) >= max(3, ceiling(0.1 * nrow(mat))), drop = FALSE]
  if (ncol(mat) < 5) {
    message("Too few taxa after filtering for Maaslin2 in ", domain)
    return(invisible(NULL))
  }

  rel <- decostand(mat, method = "total")
  feat <- as.data.frame(rel)
  meta <- as(sample_data(ps), "data.frame")
  meta$trt2 <- factor(meta$trt2, levels = c("NT", "till"))

  maas_out <- file.path(out_subdir, sprintf("%s_maaslin2", domain))
  dir.create(maas_out, showWarnings = FALSE, recursive = TRUE)

  suppressMessages(
    Maaslin2::Maaslin2(
      input_data = feat,
      input_metadata = meta,
      output = maas_out,
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

  res_file <- file.path(maas_out, "all_results.tsv")
  if (!file.exists(res_file)) return(invisible(NULL))
  res <- readr::read_tsv(res_file, show_col_types = FALSE)

  res2 <- res %>%
    filter(metadata == "trt2", value == "till") %>%
    left_join(select(tax_map, asv, taxonomy), by = c("feature" = "asv")) %>%
    mutate(
      taxonomy = ifelse(is.na(taxonomy) | taxonomy == "", "taxonomy_unavailable", taxonomy),
      feature_tax = paste0(feature, " | ", taxonomy),
      feature_tax_short = ifelse(nchar(feature_tax) > 95, paste0(substr(feature_tax, 1, 95), "..."), feature_tax),
      direction = ifelse(coef > 0, "Higher in till", "Higher in NT")
    ) %>%
    arrange(qval)

  write_csv(res2, file.path(out_subdir, sprintf("%s_maaslin2_nt_vs_till_results.csv", domain)))

  topn <- res2 %>% filter(qval <= 0.25) %>% slice_head(n = 30)
  if (nrow(topn) == 0) return(invisible(NULL))

  p <- ggplot(topn, aes(x = reorder(feature_tax_short, coef), y = coef, color = direction)) +
    geom_segment(aes(xend = feature_tax_short, y = 0, yend = coef), linewidth = 0.8, alpha = 0.7) +
    geom_point(size = 2.7) +
    coord_flip() +
    scale_color_manual(values = c("Higher in NT" = "#1b9e77", "Higher in till" = "#d95f02")) +
    labs(title = sprintf("%s MaAsLin2 NT vs till", domain), x = "ASV | taxonomy", y = "Coefficient") +
    theme_bw(base_size = 10)

  ggsave(file.path(out_subdir, sprintf("%s_maaslin2_lollipop.png", domain)), p, width = 9.3, height = 7.2, dpi = 320)

  res2
}

build_corr_network <- function(mat, prefix = "Tax", cor_cut = 0.6, q_cut = 0.05, max_taxa = 250, add_prefix = TRUE){
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

  if (length(keep_e) == 0) {
    return(list(graph = make_empty_graph(), edges = tibble()))
  }

  ed <- tibble(
    from = colnames(R)[ut[keep_e, 1]],
    to = colnames(R)[ut[keep_e, 2]],
    weight = R[ut[keep_e]],
    qval = qvals[keep_e],
    sign = ifelse(R[ut[keep_e]] > 0, "positive", "negative")
  )

  if (isTRUE(add_prefix)) {
    ed <- ed %>% mutate(from = paste0(prefix, "|", from), to = paste0(prefix, "|", to))
  }

  g <- graph_from_data_frame(ed, directed = FALSE)
  list(graph = g, edges = ed)
}

network_stats <- function(g, label){
  if (vcount(g) == 0 || ecount(g) == 0) {
    return(tibble(
      label = label, nodes = vcount(g), edges = ecount(g), density = NA_real_,
      mean_degree = NA_real_, transitivity = NA_real_, modularity = NA_real_
    ))
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

extract_node_table <- function(g, tax_master){
  if (vcount(g) == 0) {
    return(tibble(node = character(), domain = character(), asv = character(), taxonomy = character(), degree = numeric(), hub = logical()))
  }

  node_names <- V(g)$name
  parts <- strsplit(node_names, "\\|")
  node_domain <- vapply(parts, function(x) x[1], character(1))
  node_asv <- vapply(parts, function(x) x[length(x)], character(1))
  node_deg <- degree(g)
  cutoff <- unname(quantile(node_deg, 0.9, na.rm = TRUE))
  node_hub <- node_deg >= cutoff

  tibble(
    node = node_names,
    domain = node_domain,
    asv = node_asv,
    degree = as.numeric(node_deg),
    hub = node_hub
  ) %>%
    left_join(tax_master, by = c("domain", "asv")) %>%
    mutate(taxonomy = ifelse(is.na(taxonomy), "taxonomy_unavailable", taxonomy))
}

draw_network_plot <- function(g, title, cross_domain = FALSE, maaslin_sig = NULL, label_asv = NULL, common_asv = NULL, max_labels = NULL, show_legend = TRUE, hub_labels = NULL){
  if (vcount(g) == 0 || ecount(g) == 0) {
    plot.new(); title(main = paste(title, "(no edges after thresholding)"))
    return(invisible(NULL))
  }

  set.seed(42)
  node_names <- V(g)$name
  parts <- strsplit(node_names, "\\|")
  node_domain <- vapply(parts, function(x) x[1], character(1))
  node_asv <- vapply(parts, function(x) x[length(x)], character(1))

  deg <- degree(g)
  hub_cut <- unname(quantile(deg, 0.9, na.rm = TRUE))
  is_hub <- deg >= hub_cut
  max_labels <- ifelse(is.null(max_labels), ifelse(cross_domain, 30, 20), max_labels)
  label_nodes <- is_hub
  if (!is.null(label_asv)) label_nodes <- label_nodes | (node_asv %in% label_asv)
  if (!is.null(common_asv)) label_nodes <- label_nodes | (node_asv %in% common_asv)
  score <- deg +
    ifelse(!is.null(label_asv) & (node_asv %in% label_asv), 1000, 0) +
    ifelse(!is.null(common_asv) & (node_asv %in% common_asv), 250, 0)

  # Use one layout for both plotting and label-spacing decisions to avoid overlap.
  lay <- layout_with_fr(g, weights = NA)

  candidate_idx <- which(label_nodes)
  if (length(candidate_idx) > 0) {
    ord <- candidate_idx[order(score[candidate_idx], decreasing = TRUE)]
    diag_range <- sqrt((max(lay[, 1]) - min(lay[, 1]))^2 + (max(lay[, 2]) - min(lay[, 2]))^2)
    min_sep <- ifelse(cross_domain, 0.052, 0.06) * diag_range

    chosen <- integer(0)
    for (idx in ord) {
      if (length(chosen) >= max_labels) break
      if (length(chosen) == 0) {
        chosen <- c(chosen, idx)
      } else {
        d <- sqrt((lay[chosen, 1] - lay[idx, 1])^2 + (lay[chosen, 2] - lay[idx, 2])^2)
        if (all(d >= min_sep)) chosen <- c(chosen, idx)
      }
    }

    # If spacing is too strict in dense areas, backfill with highest-priority nodes.
    if (length(chosen) < min(max_labels, length(ord))) {
      need <- min(max_labels, length(ord)) - length(chosen)
      fill <- setdiff(ord, chosen)
      if (need > 0 && length(fill) > 0) chosen <- c(chosen, head(fill, need))
    }

    label_nodes <- seq_along(deg) %in% chosen
  }

  E(g)$col <- ifelse(E(g)$weight > 0, "#1f78b4", "#e31a1c")
  E(g)$lwd <- scales::rescale(abs(E(g)$weight), to = c(0.4, 2.5))

  if (cross_domain) {
    domain_cols <- c("16S" = "#2E86AB", "ITS" = "#F18F01", "AMF" = "#6A994E")
    vcols <- domain_cols[node_domain]
    vcols[is.na(vcols)] <- "#999999"
    vlabels <- ifelse(label_nodes, paste0(node_domain, ":", node_asv), "")
  } else {
    vcols <- "#8DBCC7"
    vlabels <- ifelse(label_nodes, node_asv, "")
  }

  if (!is.null(common_asv)) {
    common_labeled <- label_nodes & (node_asv %in% common_asv)
    vlabels[common_labeled] <- paste0(vlabels[common_labeled], "*")
  }
  if (!is.null(label_asv)) {
    important_labeled <- label_nodes & (node_asv %in% label_asv)
    vlabels[important_labeled] <- paste0(vlabels[important_labeled], "!")
  }

  enrich_dir <- rep(NA_character_, length(node_asv))
  is_enriched <- rep(FALSE, length(node_asv))
  if (!is.null(maaslin_sig) && nrow(maaslin_sig) > 0) {
    node_tbl <- tibble(domain = node_domain, asv = node_asv, idx = seq_along(node_asv))
    merged <- node_tbl %>%
      left_join(maaslin_sig %>% select(domain, asv, direction), by = c("domain", "asv"))
    enrich_dir <- merged$direction
    is_enriched <- !is.na(enrich_dir)
  }

  frame_col <- ifelse(
    is_enriched & enrich_dir == "Higher in NT", "#1b9e77",
    ifelse(is_enriched & enrich_dir == "Higher in till", "#d95f02", "#111111")
  )
  vsize <- ifelse(is_enriched, 6.7, 5)

  plot(g,
       vertex.size = vsize,
       vertex.label = vlabels,
      vertex.label.cex = ifelse(label_nodes, 0.28, 0),
       vertex.label.color = "#222222",
      vertex.label.dist = ifelse(label_nodes, 0.62, 0),
       vertex.color = vcols,
       vertex.frame.color = frame_col,
       vertex.frame.width = ifelse(is_hub & is_enriched, 3.2, ifelse(is_hub, 2.6, ifelse(is_enriched, 1.8, 0.7))),
       edge.color = E(g)$col,
       edge.width = E(g)$lwd,
      layout = lay,
       main = title)

  if (show_legend && cross_domain) {
    legend("topleft",
           legend = c("16S", "ITS", "AMF", "positive edge", "negative edge", "hub node", "MaAsLin2: higher in NT", "MaAsLin2: higher in till"),
           col = c("#2E86AB", "#F18F01", "#6A994E", "#1f78b4", "#e31a1c", "#111111", "#1b9e77", "#d95f02"),
           pch = c(16, 16, 16, NA, NA, 21, 21, 21),
           pt.bg = c("#2E86AB", "#F18F01", "#6A994E", NA, NA, "white", "white", "white"),
           lty = c(NA, NA, NA, 1, 1, NA, NA, NA),
           lwd = c(NA, NA, NA, 3, 3, NA, NA, NA),
           pt.cex = c(1, 1, 1, NA, NA, 1.3, 1.3, 1.3),
           pt.lwd = c(NA, NA, NA, NA, NA, 2.6, 2.2, 2.2),
           cex = 0.78,
           bty = "n")
  } else if (show_legend) {
    legend("topleft",
           legend = c("positive edge", "negative edge", "hub node", "MaAsLin2: higher in NT", "MaAsLin2: higher in till"),
           col = c("#1f78b4", "#e31a1c", "#111111", "#1b9e77", "#d95f02"),
           pch = c(NA, NA, 21, 21, 21),
           pt.bg = c(NA, NA, "white", "white", "white"),
           lty = c(1, 1, NA, NA, NA),
           lwd = c(3, 3, NA, NA, NA),
           pt.cex = c(NA, NA, 1.3, 1.3, 1.3),
           pt.lwd = c(NA, NA, 2.6, 2.2, 2.2),
           cex = 0.82,
           bty = "n")
  }

  if (!is.null(hub_labels) && length(hub_labels) > 0) {
    max_hub_show <- 22
    show_hubs <- hub_labels
    if (length(show_hubs) > max_hub_show) {
      show_hubs <- c(show_hubs[1:max_hub_show], sprintf("... +%d more", length(hub_labels) - max_hub_show))
    }
    hub_cex <- ifelse(length(show_hubs) > 14, 0.48, 0.56)
    hub_ncol <- ifelse(length(show_hubs) > 14, 2, 1)
    legend("right",
           inset = c(-0.33, 0),
           xpd = NA,
           title = "Hub ASVs",
           legend = show_hubs,
           cex = hub_cex,
           ncol = hub_ncol,
           bty = "o",
           y.intersp = 0.85,
           x.intersp = 0.6)
  }

  invisible(NULL)
}

plot_network_base <- function(g, filepath, title, cross_domain = FALSE, maaslin_sig = NULL){
  png(filepath, width = 2100, height = 1650, res = 240)
  draw_network_plot(g, title = title, cross_domain = cross_domain, maaslin_sig = maaslin_sig, show_legend = TRUE)
  dev.off()
}

plot_domain_network_panel <- function(domain, g_nt, g_till, out_subdir, maaslin_sig, nt_nodes, till_nodes){
  common_asv <- intersect(nt_nodes$asv, till_nodes$asv)
  imp_nt <- nt_nodes %>% filter(hub) %>% pull(asv)
  imp_till <- till_nodes %>% filter(hub) %>% pull(asv)

  nt_hub_labels <- nt_nodes %>%
    filter(hub) %>%
    arrange(desc(degree), asv) %>%
    transmute(lbl = paste0(asv, " (", degree, ")")) %>%
    pull(lbl)
  till_hub_labels <- till_nodes %>%
    filter(hub) %>%
    arrange(desc(degree), asv) %>%
    transmute(lbl = paste0(asv, " (", degree, ")")) %>%
    pull(lbl)

  panel_file <- file.path(out_subdir, sprintf("%s_network_panel_labeled_overlay.png", domain))
  png(panel_file, width = 5200, height = 2600, res = 260)
  old_par <- par(no.readonly = TRUE)
  layout(matrix(c(1, 2, 5,
                  3, 4, 5), nrow = 2, byrow = TRUE), widths = c(1.08, 1.08, 0.95))
  par(mar = c(0.8, 0.8, 2.4, 1.1))

  draw_network_plot(
    g_nt,
    title = sprintf("%s NT network", domain),
    cross_domain = FALSE,
    label_asv = imp_nt,
    common_asv = common_asv,
    max_labels = 42,
    show_legend = FALSE,
    hub_labels = NULL
  )
  draw_network_plot(
    g_nt,
    title = sprintf("%s NT network + MaAsLin2", domain),
    cross_domain = FALSE,
    maaslin_sig = maaslin_sig,
    label_asv = imp_nt,
    common_asv = common_asv,
    max_labels = 42,
    show_legend = FALSE,
    hub_labels = NULL
  )
  draw_network_plot(
    g_till,
    title = sprintf("%s till network", domain),
    cross_domain = FALSE,
    label_asv = imp_till,
    common_asv = common_asv,
    max_labels = 42,
    show_legend = FALSE,
    hub_labels = NULL
  )
  draw_network_plot(
    g_till,
    title = sprintf("%s till network + MaAsLin2", domain),
    cross_domain = FALSE,
    maaslin_sig = maaslin_sig,
    label_asv = imp_till,
    common_asv = common_asv,
    max_labels = 42,
    show_legend = FALSE,
    hub_labels = NULL
  )

  par(mar = c(1.2, 0.3, 2.4, 0.3))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  legend("top", inset = c(0.02, 0), bty = "n", cex = 0.9,
         legend = c("positive edge", "negative edge", "hub node", "MaAsLin2: higher in NT", "MaAsLin2: higher in till", "label marker: * common node", "label marker: ! hub node"),
         col = c("#1f78b4", "#e31a1c", "#111111", "#1b9e77", "#d95f02", "#333333", "#333333"),
         pch = c(NA, NA, 21, 21, 21, NA, NA),
         pt.bg = c(NA, NA, "white", "white", "white", NA, NA),
         lty = c(1, 1, NA, NA, NA, NA, NA),
         lwd = c(3, 3, NA, NA, NA, NA, NA),
         pt.cex = c(NA, NA, 1.3, 1.2, 1.2, NA, NA),
         pt.lwd = c(NA, NA, 2.6, 2.1, 2.1, NA, NA),
         x.intersp = 0.7, y.intersp = 1.0)

  fmt_hubs <- function(v, maxn = 24) {
    if (length(v) == 0) return("none")
    if (length(v) <= maxn) return(v)
    c(v[1:maxn], sprintf("... +%d more", length(v) - maxn))
  }
  nt_show <- fmt_hubs(nt_hub_labels, maxn = 24)
  till_show <- fmt_hubs(till_hub_labels, maxn = 24)

  text(0.02, 0.73, labels = "Hub ASVs (NT)", adj = c(0, 1), cex = 0.88, font = 2)
  text(0.02, 0.69, labels = paste(nt_show, collapse = "\n"), adj = c(0, 1), cex = 0.63)

  text(0.52, 0.73, labels = "Hub ASVs (till)", adj = c(0, 1), cex = 0.88, font = 2)
  text(0.52, 0.69, labels = paste(till_show, collapse = "\n"), adj = c(0, 1), cex = 0.63)

  par(old_par)
  dev.off()
}

compare_networks_degree <- function(g1, g2, domain){
  if (vcount(g1) < 3 || vcount(g2) < 3 || ecount(g1) == 0 || ecount(g2) == 0) {
    return(tibble(domain = domain, test = c("wilcox_degree", "ks_degree"), p_value = NA_real_))
  }
  d1 <- degree(g1)
  d2 <- degree(g2)
  w <- suppressWarnings(wilcox.test(d1, d2, exact = FALSE))
  k <- suppressWarnings(ks.test(d1, d2))
  tibble(domain = domain, test = c("wilcox_degree", "ks_degree"), p_value = c(w$p.value, k$p.value))
}

write_report <- function(out_dir){
  ds <- read_csv(file.path(out_dir, "dataset_summary.csv"), show_col_types = FALSE)
  net <- read_csv(file.path(out_dir, "network_stats_summary.csv"), show_col_types = FALSE)
  deg <- read_csv(file.path(out_dir, "network_degree_tests_nt_vs_till.csv"), show_col_types = FALSE)

  a16 <- read_csv(file.path(out_dir, "16S", "16S_alpha_wilcox.csv"), show_col_types = FALSE)
  ait <- read_csv(file.path(out_dir, "ITS", "ITS_alpha_wilcox.csv"), show_col_types = FALSE)
  aam <- read_csv(file.path(out_dir, "AMF", "AMF_alpha_wilcox.csv"), show_col_types = FALSE)

  n16 <- read_csv(file.path(out_dir, "16S", "16S_nmds_permanova.csv"), show_col_types = FALSE)
  nit <- read_csv(file.path(out_dir, "ITS", "ITS_nmds_permanova.csv"), show_col_types = FALSE)
  nam <- read_csv(file.path(out_dir, "AMF", "AMF_nmds_permanova.csv"), show_col_types = FALSE)
  jnt <- read_csv(file.path(out_dir, "joint_rpca_proxy_permanova.csv"), show_col_types = FALSE)

  get_model_line <- function(tbl) {
    tbl %>% filter(term == "Model") %>% slice(1)
  }

  n16m <- get_model_line(n16)
  nitm <- get_model_line(nit)
  namm <- get_model_line(nam)
  jntm <- get_model_line(jnt)

  lines <- c(
    "# MF Soil NT vs till Multi-domain Analysis",
    "",
    "## Methods",
    "- Data objects: phyloseq RDS files (16S, ITS, AMF).",
    "- Sample filter: MF soils only; treatment contrast mapped to NT vs till.",
    "- Alpha diversity: Observed richness and Shannon using phyloseq::estimate_richness.",
    "- Alpha statistics: Wilcoxon rank-sum tests (NT vs till).",
    "- Rarefaction: vegan::rarecurve.",
    "- Beta diversity: Bray-Curtis distance and NMDS (vegan::metaMDS).",
    "- NMDS statistics: PERMANOVA via vegan::adonis2.",
    "- Joint integration: CLR transform then PCA as a joint RPCA proxy across concatenated 16S+ITS+AMF features.",
    "- Differential abundance: MaAsLin2 linear models with BH correction.",
    "- Co-occurrence networks: Spearman correlations on CLR-transformed features; edge criteria |rho|>=0.6 and BH q<0.05.",
    "- Hub nodes: top 10% by degree; highlighted with bold borders.",
    "",
    "## Tools",
    "- R 4.5.2 (module stack: gcc/12.3.0 + R/4.5.2).",
    "- R packages: phyloseq, vegan, ggplot2, dplyr, tidyr, tibble, stringr, readr, purrr, igraph, Hmisc, Maaslin2, scales.",
    "",
    "## Results",
    "### Dataset counts",
    paste0("- 16S: n=", ds$n_samples[ds$domain == "16S"], ", NT=", ds$n_nt[ds$domain == "16S"], ", till=", ds$n_till[ds$domain == "16S"], ", taxa=", ds$n_taxa[ds$domain == "16S"]),
    paste0("- ITS: n=", ds$n_samples[ds$domain == "ITS"], ", NT=", ds$n_nt[ds$domain == "ITS"], ", till=", ds$n_till[ds$domain == "ITS"], ", taxa=", ds$n_taxa[ds$domain == "ITS"]),
    paste0("- AMF: n=", ds$n_samples[ds$domain == "AMF"], ", NT=", ds$n_nt[ds$domain == "AMF"], ", till=", ds$n_till[ds$domain == "AMF"], ", taxa=", ds$n_taxa[ds$domain == "AMF"]),
    "",
    "### Alpha diversity (Wilcoxon p)",
    paste0("- 16S Observed p=", signif(a16$p_value[a16$metric == "Observed"], 3), "; Shannon p=", signif(a16$p_value[a16$metric == "Shannon"], 3)),
    paste0("- ITS Observed p=", signif(ait$p_value[ait$metric == "Observed"], 3), "; Shannon p=", signif(ait$p_value[ait$metric == "Shannon"], 3)),
    paste0("- AMF Observed p=", signif(aam$p_value[aam$metric == "Observed"], 3), "; Shannon p=", signif(aam$p_value[aam$metric == "Shannon"], 3)),
    "",
    "### NMDS PERMANOVA (Model term)",
    paste0("- 16S: R2=", signif(n16m$R2, 3), ", p=", signif(n16m$`Pr(>F)`, 3)),
    paste0("- ITS: R2=", signif(nitm$R2, 3), ", p=", signif(nitm$`Pr(>F)`, 3)),
    paste0("- AMF: R2=", signif(namm$R2, 3), ", p=", signif(namm$`Pr(>F)`, 3)),
    paste0("- Joint RPCA proxy: R2=", signif(jntm$R2, 3), ", p=", signif(jntm$`Pr(>F)`, 3)),
    "",
    "### Network statistics",
    "- See network_stats_summary.csv for nodes, edges, density, mean degree, transitivity, and modularity.",
    "- See network_degree_tests_nt_vs_till.csv for Wilcoxon and KS tests on degree distributions.",
    "",
    "## Discussion",
    "- NT vs till effects were strongest in 16S and ITS community composition by PERMANOVA, while AMF composition showed weaker separation.",
    "- Alpha diversity differences were modest; AMF Shannon showed the closest trend-level shift.",
    "- Network topology varied substantially between NT and till in 16S and AMF, including large degree-distribution shifts.",
    "- Cross-domain network structure indicates treatment-dependent rewiring, with domain-tagged nodes enabling inspection of which kingdom contributes to hubs.",
    "- Signed-correlation networks are sensitive to filtering and threshold choices; conclusions should be interpreted as hypothesis-generating unless externally validated.",
    "",
    "## Output files",
    "- Domain figures and tables are under analysis_outputs/16S, analysis_outputs/ITS, analysis_outputs/AMF.",
    "- Cross-domain outputs are in analysis_outputs root."
  )

  writeLines(lines, file.path(out_dir, "methods_tools_results_discussion.md"))
}

all_ps <- list()
for (nm in names(rds_files)) {
  all_ps[[nm]] <- readRDS(file.path(root_dir, rds_files[[nm]]))
}

mf_ps <- imap(all_ps, subset_mf_nt_till)
taxonomy_maps <- imap(mf_ps, build_taxonomy_map)
tax_master <- bind_rows(taxonomy_maps)

summary_rows <- list()
network_stat_rows <- list()
network_test_rows <- list()
maaslin_sig_rows <- list()

for (domain in names(mf_ps)) {
  out_subdir <- file.path(out_dir, domain)
  dir.create(out_subdir, showWarnings = FALSE, recursive = TRUE)
  ps <- mf_ps[[domain]]

  alpha_stats_and_plots(ps, domain, out_subdir)
  rarefaction_plot(ps, domain, out_subdir)
  nmds_and_stats(ps, domain, out_subdir)
  maas_res <- run_maaslin_and_lollipop(ps, domain, out_subdir, taxonomy_maps[[domain]])
  if (!is.null(maas_res) && nrow(maas_res) > 0) {
    maaslin_sig_rows[[domain]] <- maas_res %>%
      filter(qval <= 0.25) %>%
      transmute(domain = domain, asv = feature, direction = direction)
  } else {
    maaslin_sig_rows[[domain]] <- tibble(domain = character(), asv = character(), direction = character())
  }

  nt_ids <- sample_names(ps)[sample_data(ps)$trt2 == "NT"]
  till_ids <- sample_names(ps)[sample_data(ps)$trt2 == "till"]

  nt_mat <- otu_samples_by_taxa(prune_samples(nt_ids, ps))
  till_mat <- otu_samples_by_taxa(prune_samples(till_ids, ps))

  net_nt <- build_corr_network(nt_mat, prefix = domain)
  net_till <- build_corr_network(till_mat, prefix = domain)

  plot_network_base(net_nt$graph, file.path(out_subdir, sprintf("%s_network_NT.png", domain)), sprintf("%s NT network", domain), cross_domain = FALSE)
  plot_network_base(net_till$graph, file.path(out_subdir, sprintf("%s_network_till.png", domain)), sprintf("%s till network", domain), cross_domain = FALSE)
  plot_network_base(
    net_nt$graph,
    file.path(out_subdir, sprintf("%s_network_NT_maaslin_overlay.png", domain)),
    sprintf("%s NT network (MaAsLin2 overlay)", domain),
    cross_domain = FALSE,
    maaslin_sig = maaslin_sig_rows[[domain]]
  )
  plot_network_base(
    net_till$graph,
    file.path(out_subdir, sprintf("%s_network_till_maaslin_overlay.png", domain)),
    sprintf("%s till network (MaAsLin2 overlay)", domain),
    cross_domain = FALSE,
    maaslin_sig = maaslin_sig_rows[[domain]]
  )

  nt_nodes <- extract_node_table(net_nt$graph, tax_master)
  till_nodes <- extract_node_table(net_till$graph, tax_master)
  write_csv(nt_nodes, file.path(out_subdir, sprintf("%s_network_NT_nodes_taxonomy.csv", domain)))
  write_csv(till_nodes, file.path(out_subdir, sprintf("%s_network_till_nodes_taxonomy.csv", domain)))

  plot_domain_network_panel(
    domain = domain,
    g_nt = net_nt$graph,
    g_till = net_till$graph,
    out_subdir = out_subdir,
    maaslin_sig = maaslin_sig_rows[[domain]],
    nt_nodes = nt_nodes,
    till_nodes = till_nodes
  )

  st_nt <- network_stats(net_nt$graph, paste(domain, "NT"))
  st_till <- network_stats(net_till$graph, paste(domain, "till"))

  network_stat_rows[[paste0(domain, "_NT")]] <- st_nt
  network_stat_rows[[paste0(domain, "_till")]] <- st_till
  network_test_rows[[domain]] <- compare_networks_degree(net_nt$graph, net_till$graph, domain)

  summary_rows[[domain]] <- tibble(
    domain = domain,
    n_samples = nsamples(ps),
    n_nt = length(nt_ids),
    n_till = length(till_ids),
    n_taxa = ntaxa(ps)
  )
}

# Joint RPCA-like integration across all domains (CLR-PCA on concatenated features)
common_ids <- Reduce(intersect, lapply(mf_ps, sample_names))
joint_list <- list()
for (domain in names(mf_ps)) {
  m <- otu_samples_by_taxa(prune_samples(common_ids, mf_ps[[domain]]))
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

pca <- prcomp(joint_mat, center = TRUE, scale. = TRUE)
meta_joint <- as(sample_data(mf_ps[[1]]), "data.frame") %>%
  rownames_to_column("sample_key") %>%
  filter(sample_key %in% common_ids)

scores <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("sample_key") %>%
  left_join(meta_joint, by = "sample_key")

p_joint <- ggplot(scores, aes(PC1, PC2, color = trt2)) +
  geom_point(size = 2.6, alpha = 0.9) +
  stat_ellipse(aes(fill = trt2), geom = "polygon", alpha = 0.15, color = NA) +
  scale_color_manual(values = c("NT" = "#1b9e77", "till" = "#d95f02")) +
  scale_fill_manual(values = c("NT" = "#1b9e77", "till" = "#d95f02")) +
  labs(
    title = "Joint RPCA proxy (CLR-PCA) across 16S + ITS + AMF",
    subtitle = sprintf("PC1 %.1f%%, PC2 %.1f%% variance", 100 * summary(pca)$importance[2, 1], 100 * summary(pca)$importance[2, 2]),
    color = "Treatment", fill = "Treatment"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "joint_rpca_proxy_plot.png"), p_joint, width = 7, height = 5.5, dpi = 320)

d_joint <- dist(scale(joint_mat))
meta_perm <- scores %>% select(sample_key, trt2)
perm_joint <- adonis2(d_joint ~ trt2, data = meta_perm, permutations = 999)
perm_joint_tbl <- as.data.frame(perm_joint)
perm_joint_tbl$term <- rownames(perm_joint_tbl)
write_csv(perm_joint_tbl, file.path(out_dir, "joint_rpca_proxy_permanova.csv"))

# Cross-domain co-occurrence networks for NT and till
make_cross_mat <- function(trt_name){
  mats <- list()
  for (domain in names(mf_ps)) {
    ps <- mf_ps[[domain]]
    ids <- sample_names(ps)[sample_data(ps)$trt2 == trt_name & sample_names(ps) %in% common_ids]
    m <- otu_samples_by_taxa(prune_samples(ids, ps))
    prev <- colMeans(m > 0)
    m <- m[, prev >= 0.1, drop = FALSE]
    if (ncol(m) > 120) {
      m <- m[, order(colSums(m), decreasing = TRUE)[1:120], drop = FALSE]
    }
    colnames(m) <- paste0(domain, "|", colnames(m))
    mats[[domain]] <- m
  }

  sample_order <- sample_names(prune_samples(sample_names(mf_ps[[1]])[sample_data(mf_ps[[1]])$trt2 == trt_name & sample_names(mf_ps[[1]]) %in% common_ids], mf_ps[[1]]))
  mats <- lapply(mats, function(m) m[match(sample_order, rownames(m)), , drop = FALSE])
  do.call(cbind, mats)
}

cross_nt_mat <- make_cross_mat("NT")
cross_till_mat <- make_cross_mat("till")

cross_nt <- build_corr_network(cross_nt_mat, add_prefix = FALSE)
cross_till <- build_corr_network(cross_till_mat, add_prefix = FALSE)

maaslin_sig_master <- bind_rows(maaslin_sig_rows)

plot_network_base(cross_nt$graph, file.path(out_dir, "cross_domain_network_NT.png"), "Cross-domain NT network", cross_domain = TRUE)
plot_network_base(cross_till$graph, file.path(out_dir, "cross_domain_network_till.png"), "Cross-domain till network", cross_domain = TRUE)
plot_network_base(
  cross_nt$graph,
  file.path(out_dir, "cross_domain_network_NT_maaslin_overlay.png"),
  "Cross-domain NT network (MaAsLin2 overlay)",
  cross_domain = TRUE,
  maaslin_sig = maaslin_sig_master
)
plot_network_base(
  cross_till$graph,
  file.path(out_dir, "cross_domain_network_till_maaslin_overlay.png"),
  "Cross-domain till network (MaAsLin2 overlay)",
  cross_domain = TRUE,
  maaslin_sig = maaslin_sig_master
)

cross_nt_nodes <- extract_node_table(cross_nt$graph, tax_master)
cross_till_nodes <- extract_node_table(cross_till$graph, tax_master)
write_csv(cross_nt_nodes, file.path(out_dir, "cross_domain_network_NT_nodes_taxonomy.csv"))
write_csv(cross_till_nodes, file.path(out_dir, "cross_domain_network_till_nodes_taxonomy.csv"))

network_stat_rows[["cross_NT"]] <- network_stats(cross_nt$graph, "cross NT")
network_stat_rows[["cross_till"]] <- network_stats(cross_till$graph, "cross till")
network_test_rows[["cross"]] <- compare_networks_degree(cross_nt$graph, cross_till$graph, "cross")

summary_tbl <- bind_rows(summary_rows)
write_csv(summary_tbl, file.path(out_dir, "dataset_summary.csv"))

net_stats_tbl <- bind_rows(network_stat_rows)
write_csv(net_stats_tbl, file.path(out_dir, "network_stats_summary.csv"))

net_test_tbl <- bind_rows(network_test_rows)
write_csv(net_test_tbl, file.path(out_dir, "network_degree_tests_nt_vs_till.csv"))

write_report(out_dir)

message("Analysis complete. Outputs in: ", normalizePath(out_dir))
