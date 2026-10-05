# ============================================================
# 06_cross_segment_summary.R
# Cross-segment ARI comparison: bar chart of ARI per segment
# (RSCU vs constellation, RSCU vs host, group vs host)
# ============================================================

source("c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/scripts/common_bat.R")
dirs <- setup_dirs()

# ---- Load ARI results ----
ari_path <- file.path(dirs$tabdir, "ari_results.csv")
if (!file.exists(ari_path)) {
  stop("ari_results.csv not found. Run 05_tanglegram_ari.R first.")
}
ari_df <- read.csv(ari_path, stringsAsFactors = FALSE)

# ---- Load tanglegram results if available ----
tang_path <- file.path(dirs$tabdir, "tanglegram_results.csv")
tang_df <- if (file.exists(tang_path)) read.csv(tang_path, stringsAsFactors = FALSE) else NULL

# ---- Reshape ARI for plotting ----
ari_long <- ari_df %>%
  select(segment, ARI_rscu_vs_group, ARI_rscu_vs_host, ARI_group_vs_host) %>%
  pivot_longer(cols = starts_with("ARI_"), names_to = "comparison", values_to = "ARI") %>%
  mutate(
    comparison = recode(comparison,
      "ARI_rscu_vs_group" = "RSCU vs Constellation",
      "ARI_rscu_vs_host" = "RSCU vs Host",
      "ARI_group_vs_host" = "Constellation vs Host"
    ),
    segment = factor(segment, levels = c("VP4", "VP6", "VP7", "NSP4"))
  )

# ---- ARI bar chart ----
p_ari <- ggplot(ari_long, aes(x = segment, y = ARI, fill = comparison)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", ARI)),
            position = position_dodge(width = 0.8), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c(
    "RSCU vs Constellation" = "#0279EE",
    "RSCU vs Host" = "#FF9400",
    "Constellation vs Host" = "#75A025"
  )) +
  labs(title = "Cross-Segment ARI Comparison",
       x = "Segment", y = "Adjusted Rand Index",
       fill = "Comparison") +
  theme_cub() +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(-0.3, 1.0)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 1, linetype = "dotted", color = "gray70")

png(file.path(dirs$figdir, "ari_cross_segment.png"),
    width = 10, height = 6, units = "in", res = 300)
print(p_ari)
dev.off()
cat("Saved: ari_cross_segment.png\n")

# ---- Entanglement + Baker's gamma summary ----
if (!is.null(tang_df)) {
  p_tang <- ggplot(tang_df, aes(x = segment, y = entanglement)) +
    geom_col(fill = "#0279EE", width = 0.6, alpha = 0.8) +
    geom_text(aes(label = sprintf("%.4f", entanglement)),
              vjust = -0.3, size = 3.5) +
    labs(title = "Tanglegram Entanglement (RSCU vs Phylogeny)",
         x = "Segment", y = "Entanglement (0=perfect match, 1=no match)") +
    theme_cub() +
    coord_cartesian(ylim = c(0, 1))

  png(file.path(dirs$figdir, "entanglement_cross_segment.png"),
      width = 8, height = 5, units = "in", res = 300)
  print(p_tang)
  dev.off()
  cat("Saved: entanglement_cross_segment.png\n")
}

# ---- Summary table ----
summary_table <- ari_df %>%
  mutate(
    concordance = ifelse(ARI_rscu_vs_group > 0.5, "Strong",
                  ifelse(ARI_rscu_vs_group > 0.2, "Moderate",
                   ifelse(ARI_rscu_vs_group > 0, "Weak", "None"))),
    host_signal = ifelse(ARI_rscu_vs_host > 0.5, "Strong",
                  ifelse(ARI_rscu_vs_host > 0.2, "Moderate",
                   ifelse(ARI_rscu_vs_host > 0, "Weak", "None")))
  )

write.csv(summary_table, file.path(dirs$tabdir, "cross_segment_summary.csv"),
          row.names = FALSE)
cat("\nSaved: cross_segment_summary.csv\n")
print(summary_table)

# ---- Interpretation ----
cat("\n=== Interpretation ===\n")
for (i in 1:nrow(summary_table)) {
  seg <- summary_table$segment[i]
  ari_g <- summary_table$ARI_rscu_vs_group[i]
  ari_h <- summary_table$ARI_rscu_vs_host[i]
  conc <- summary_table$concordance[i]
  host_sig <- summary_table$host_signal[i]
  cat(sprintf("  %s: RSCU-group ARI=%.3f (%s concordance), RSCU-host ARI=%.3f (%s host signal)\n",
              seg, ari_g, tolower(conc), ari_h, tolower(host_sig)))
}

cat("\n=== Cross-Segment Summary Complete ===\n")
