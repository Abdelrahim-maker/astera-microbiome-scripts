#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(purrr)
})

out_dir <- file.path("analysis_outputs", "carbon_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ma_file <- file.path(out_dir, "carbon_maaslin_all_domains.csv")
perm_file <- file.path(out_dir, "carbon_permanova_by_domain.csv")
effect_file <- file.path(out_dir, "carbon_maaslin_effect_summary.csv")

if (!file.exists(ma_file)) stop("Missing carbon_maaslin_all_domains.csv")
if (!file.exists(perm_file)) stop("Missing carbon_permanova_by_domain.csv")
if (!file.exists(effect_file)) stop("Missing carbon_maaslin_effect_summary.csv")

ma <- read_csv(ma_file, show_col_types = FALSE)
perm <- read_csv(perm_file, show_col_types = FALSE)
effect <- read_csv(effect_file, show_col_types = FALSE)

sig_carbon <- ma %>%
  filter(effect_type %in% c("carbon_main_NT_slope", "carbon_slope_difference_till_vs_NT"), qval <= 0.25) %>%
  transmute(domain, asv = feature, effect_type, coef, qval) %>%
  group_by(domain, asv, effect_type) %>%
  summarise(
    carbon_coef = coef[which.min(qval)],
    carbon_qval = min(qval, na.rm = TRUE),
    .groups = "drop"
  )

node_files <- c(
  file.path("analysis_outputs", "16S", "16S_network_NT_nodes_with_maaslin.csv"),
  file.path("analysis_outputs", "16S", "16S_network_till_nodes_with_maaslin.csv"),
  file.path("analysis_outputs", "ITS", "ITS_network_NT_nodes_with_maaslin.csv"),
  file.path("analysis_outputs", "ITS", "ITS_network_till_nodes_with_maaslin.csv"),
  file.path("analysis_outputs", "AMF", "AMF_network_NT_nodes_with_maaslin.csv"),
  file.path("analysis_outputs", "AMF", "AMF_network_till_nodes_with_maaslin.csv")
)

missing_nodes <- node_files[!file.exists(node_files)]
if (length(missing_nodes) > 0) {
  stop("Missing network node files: ", paste(missing_nodes, collapse = ", "))
}

nodes <- map_dfr(node_files, function(fp) {
  read_csv(fp, show_col_types = FALSE) %>%
    mutate(
      domain = as.character(domain),
      asv = as.character(asv),
      degree = suppressWarnings(as.numeric(degree)),
      hub = as.logical(hub),
      maaslin_enriched = as.logical(maaslin_enriched),
      taxonomy = as.character(taxonomy)
    )
})

node_flags <- nodes %>%
  mutate(maaslin_enriched = ifelse(is.na(maaslin_enriched), FALSE, maaslin_enriched)) %>%
  group_by(domain, asv) %>%
  summarise(
    in_network = n() > 0,
    in_hub = any(hub, na.rm = TRUE),
    enriched_in_network = any(maaslin_enriched, na.rm = TRUE),
    max_degree = max(degree, na.rm = TRUE),
    taxonomy = dplyr::first(taxonomy[!is.na(taxonomy)]),
    .groups = "drop"
  )

overlap <- sig_carbon %>%
  left_join(node_flags, by = c("domain", "asv")) %>%
  mutate(
    in_network = ifelse(is.na(in_network), FALSE, in_network),
    in_hub = ifelse(is.na(in_hub), FALSE, in_hub),
    enriched_in_network = ifelse(is.na(enriched_in_network), FALSE, enriched_in_network),
    max_degree = ifelse(is.na(max_degree), 0, max_degree),
    taxonomy = ifelse(is.na(taxonomy), "not_in_network", taxonomy)
  )

write_csv(overlap, file.path(out_dir, "carbon_taxa_network_enrichment_overlap.csv"))

overlap_summary <- overlap %>%
  group_by(domain, effect_type) %>%
  summarise(
    n_carbon_sig = n_distinct(asv),
    n_in_network = n_distinct(asv[in_network]),
    n_in_hubs = n_distinct(asv[in_hub]),
    n_enriched_in_network = n_distinct(asv[enriched_in_network]),
    frac_in_network = n_in_network / n_carbon_sig,
    frac_in_hubs = n_in_hubs / n_carbon_sig,
    frac_enriched_in_network = n_enriched_in_network / n_carbon_sig,
    .groups = "drop"
  )

write_csv(overlap_summary, file.path(out_dir, "carbon_taxa_network_enrichment_overlap_summary.csv"))

perm_plot_df <- perm %>%
  filter(test_type == "sequential", term %in% c("trt2", "C", "trt2:C")) %>%
  mutate(
    term = factor(term, levels = c("trt2", "C", "trt2:C"), labels = c("Treatment", "Carbon", "Treatment x Carbon")),
    domain = factor(domain, levels = c("16S", "ITS", "AMF")),
    sig_label = case_when(
      `Pr(>F)` <= 0.001 ~ "***",
      `Pr(>F)` <= 0.01 ~ "**",
      `Pr(>F)` <= 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    label = paste0(sprintf("R2=%.3f", R2), "\n", sig_label)
  )

p_perm <- ggplot(perm_plot_df, aes(x = term, y = domain, fill = R2)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 3.2, lineheight = 0.9) +
  scale_fill_gradient(low = "#fef0d9", high = "#b30000") +
  labs(
    title = "Carbon-aware PERMANOVA by domain",
    subtitle = "Sequential terms with effect size (R2) and significance",
    x = NULL,
    y = NULL,
    fill = "R2"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank())

ggsave(file.path(out_dir, "carbon_permanova_r2_heatmap.png"), p_perm, width = 7.2, height = 3.8, dpi = 320)

effect_plot_df <- effect %>%
  filter(effect_type %in% c("carbon_main_NT_slope", "carbon_slope_difference_till_vs_NT", "treatment_main")) %>%
  mutate(
    effect_type = recode(
      effect_type,
      carbon_main_NT_slope = "Carbon slope (NT)",
      carbon_slope_difference_till_vs_NT = "Carbon slope diff (till-NT)",
      treatment_main = "Treatment main"
    ),
    domain = factor(domain, levels = c("16S", "ITS", "AMF"))
  )

p_effect <- ggplot(effect_plot_df, aes(x = effect_type, y = n_sig_q025, fill = domain)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_text(aes(label = paste0(n_sig_q025, "/", n_total)),
            position = position_dodge(width = 0.8), vjust = -0.35, size = 3) +
  scale_fill_manual(values = c("16S" = "#1b9e77", "ITS" = "#7570b3", "AMF" = "#d95f02")) +
  labs(
    title = "Significant features in carbon-aware MaAsLin2",
    subtitle = "Counts at q <= 0.25 (labels show significant/total)",
    x = NULL,
    y = "Number of significant ASVs",
    fill = "Domain"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 10, hjust = 1))

ggsave(file.path(out_dir, "carbon_maaslin_effect_counts.png"), p_effect, width = 8.5, height = 4.5, dpi = 320)

overlap_plot_df <- overlap_summary %>%
  transmute(
    domain,
    effect_type,
    `Carbon-significant taxa` = n_carbon_sig,
    `In network` = n_in_network,
    `In hubs` = n_in_hubs,
    `Treatment-enriched in network` = n_enriched_in_network
  ) %>%
  pivot_longer(
    cols = c(`Carbon-significant taxa`, `In network`, `In hubs`, `Treatment-enriched in network`),
    names_to = "metric",
    values_to = "count"
  ) %>%
  mutate(
    effect_type = recode(
      effect_type,
      carbon_main_NT_slope = "Carbon slope (NT)",
      carbon_slope_difference_till_vs_NT = "Carbon slope diff (till-NT)"
    ),
    metric = factor(metric, levels = c(
      "Carbon-significant taxa",
      "In network",
      "In hubs",
      "Treatment-enriched in network"
    )),
    domain = factor(domain, levels = c("16S", "ITS", "AMF"))
  )

p_overlap <- ggplot(overlap_plot_df, aes(x = metric, y = count, fill = domain)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  facet_wrap(~ effect_type, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("16S" = "#1b9e77", "ITS" = "#7570b3", "AMF" = "#d95f02")) +
  labs(
    title = "Carbon-significant taxa overlap with network and enrichment",
    subtitle = "Taxa are counted once per domain and carbon effect type",
    x = NULL,
    y = "ASV count",
    fill = "Domain"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

ggsave(file.path(out_dir, "carbon_taxa_overlap_network_enrichment.png"), p_overlap, width = 9.5, height = 7.2, dpi = 320)

top_network_carbon <- overlap %>%
  filter(in_network) %>%
  group_by(domain, asv) %>%
  summarise(
    min_q = min(carbon_qval, na.rm = TRUE),
    max_abs_coef = max(abs(carbon_coef), na.rm = TRUE),
    in_hub = any(in_hub),
    enriched_in_network = any(enriched_in_network),
    max_degree = max(max_degree, na.rm = TRUE),
    taxonomy = dplyr::first(taxonomy),
    .groups = "drop"
  ) %>%
  arrange(min_q, desc(max_abs_coef))

write_csv(top_network_carbon, file.path(out_dir, "carbon_top_network_taxa_ranked.csv"))

top_plot <- top_network_carbon %>%
  slice_head(n = 30) %>%
  mutate(label = paste0(domain, ":", asv))

if (nrow(top_plot) > 0) {
  p_top <- ggplot(top_plot, aes(x = reorder(label, -log10(min_q)), y = -log10(min_q), color = enriched_in_network, shape = in_hub, size = pmax(1, max_degree))) +
    geom_point(alpha = 0.9) +
    coord_flip() +
    scale_color_manual(values = c("TRUE" = "#d95f02", "FALSE" = "#1f78b4")) +
    labs(
      title = "Top carbon-associated taxa present in networks",
      subtitle = "Color = treatment-enriched in network, shape = hub status",
      x = NULL,
      y = "-log10(min carbon q)",
      color = "Enriched in network",
      shape = "Hub",
      size = "Max degree"
    ) +
    theme_bw(base_size = 10)

  ggsave(file.path(out_dir, "carbon_top_network_taxa_dotplot.png"), p_top, width = 9.0, height = 7.5, dpi = 320)
}

cat("Carbon figures and overlap tables written to", out_dir, "\n")
