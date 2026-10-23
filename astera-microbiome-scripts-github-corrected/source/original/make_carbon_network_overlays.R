#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(igraph)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(Hmisc)
})

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
  if (length(keep) < 10) {
    keep <- order(colSums(mat), decreasing = TRUE)[seq_len(min(100, ncol(mat)))]
  }
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
    return(make_empty_graph())
  }

  ed <- tibble(
    from = colnames(R)[ut[keep_e, 1]],
    to = colnames(R)[ut[keep_e, 2]],
    weight = R[ut[keep_e]]
  )
  graph_from_data_frame(ed, directed = FALSE)
}

all_ps <- lapply(rds_files, readRDS)
mf_ps <- lapply(all_ps, subset_mf_nt_till)

carbon_tbl <- read_csv("analysis_outputs/carbon_analysis/carbon_maaslin_all_domains.csv", show_col_types = FALSE) %>%
  mutate(feature = as.character(feature))

trt_maaslin <- lapply(names(rds_files), function(domain) {
  f <- file.path("analysis_outputs", domain, paste0(domain, "_maaslin2_nt_vs_till_results.csv"))
  if (!file.exists(f)) return(tibble())
  read_csv(f, show_col_types = FALSE) %>%
    transmute(
      domain = domain,
      asv = feature,
      trt_coef = coef,
      trt_q = qval,
      trt_enriched = qval <= 0.25,
      trt_direction = ifelse(coef > 0, "Higher in till", "Higher in NT")
    )
}) %>% bind_rows()

for (domain in names(mf_ps)) {
  ps <- mf_ps[[domain]]
  mat <- otu_samples_by_taxa(ps)

  carbon_domain <- carbon_tbl %>%
    filter(domain == !!domain, qval <= 0.25, effect_type %in% c("carbon_main_NT_slope", "carbon_slope_difference_till_vs_NT")) %>%
    transmute(
      asv = feature,
      carbon_type = case_when(
        effect_type == "carbon_main_NT_slope" ~ "C slope",
        effect_type == "carbon_slope_difference_till_vs_NT" ~ "C slope diff",
        TRUE ~ "none"
      )
    ) %>%
    group_by(asv) %>%
    summarise(
      carbon_type = ifelse(all(carbon_type == "C slope"), "C slope",
                    ifelse(all(carbon_type == "C slope diff"), "C slope diff", "Both")),
      .groups = "drop"
    )

  domain_trt <- trt_maaslin %>% filter(domain == !!domain)

  for (trt in c("NT", "till")) {
    ids <- sample_names(ps)[sample_data(ps)$trt2 == trt]
    mat_sub <- mat[match(ids, rownames(mat)), , drop = FALSE]
    g <- build_corr_network(mat_sub)
    if (vcount(g) == 0 || ecount(g) == 0) next

    deg <- degree(g)
    hub_cut <- as.numeric(quantile(deg, probs = 0.9, na.rm = TRUE, names = FALSE))
    hubs <- names(deg)[deg >= hub_cut]

    nodes <- tibble(asv = V(g)$name) %>%
      left_join(carbon_domain, by = "asv") %>%
      left_join(domain_trt, by = c("asv")) %>%
      mutate(
        carbon_type = ifelse(is.na(carbon_type), "none", carbon_type),
        trt_enriched = ifelse(is.na(trt_enriched), FALSE, trt_enriched),
        trt_direction = ifelse(is.na(trt_direction), "not_enriched", trt_direction),
        degree = deg[asv],
        hub = asv %in% hubs
      )

    carbon_cols <- c(
      "none" = "#bdbdbd",
      "C slope" = "#66c2a4",
      "C slope diff" = "#fdbb84",
      "Both" = "#bcbddc"
    )

    fill_cols <- c(
      "not_enriched" = "#d9d9d9",
      "Higher in NT" = "#2c7fb8",
      "Higher in till" = "#f03b20"
    )

    V(g)$frame.color <- carbon_cols[nodes$carbon_type]
    V(g)$color <- fill_cols[nodes$trt_direction]
    V(g)$size <- ifelse(nodes$hub, 8.5, 5.5)
    V(g)$label <- ifelse(nodes$carbon_type == "none", "", nodes$asv)
    V(g)$label.cex <- 0.38
    V(g)$label.color <- "black"
    V(g)$shape <- ifelse(nodes$trt_enriched, "square", "circle")

    E(g)$color <- ifelse(E(g)$weight >= 0, "#1f78b4", "#e31a1c")
    E(g)$width <- 0.7 + 1.8 * abs(E(g)$weight)

    lay <- layout_with_fr(g, weights = abs(E(g)$weight) + 1e-6)

    out_png <- file.path("analysis_outputs", domain, paste0(domain, "_network_", trt, "_carbon_overlay.png"))
    png(out_png, width = 2000, height = 1600, res = 250)
    par(mar = c(0.3, 0.3, 3.2, 0.3))

    plot(
      g,
      layout = lay,
      vertex.frame.width = ifelse(nodes$carbon_type == "none", 0.8, 2.6),
      vertex.label.family = "sans",
      main = paste0(domain, " ", trt, " network with carbon-associated nodes")
    )

    legend(
      "topleft",
      legend = c("Carbon-associated (C slope)", "Carbon-associated (C slope diff)", "Carbon-associated (both)", "Not carbon-associated", "Treatment enriched (square)", "Hub node (larger)", "Higher in NT fill", "Higher in till fill"),
      pch = c(21, 21, 21, 21, 22, 21, 21, 21),
      pt.bg = c("#d9d9d9", "#d9d9d9", "#d9d9d9", "#d9d9d9", "#d9d9d9", "#d9d9d9", "#2c7fb8", "#f03b20"),
      col = c("#1b9e77", "#d95f02", "#542788", "#6b6b6b", "#6b6b6b", "#6b6b6b", "#6b6b6b", "#6b6b6b"),
      pt.cex = c(1.4, 1.4, 1.4, 1.4, 1.4, 1.8, 1.4, 1.4),
      bty = "n",
      cex = 0.86
    )

    mtext("Border color = carbon association; fill = treatment direction from MaAsLin2", side = 3, line = 0.4, cex = 0.8)
    dev.off()
  }
}

cat("Wrote carbon-overlay networks to analysis_outputs/<domain>/..._carbon_overlay.png\n")
