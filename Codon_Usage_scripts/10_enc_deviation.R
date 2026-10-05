#!/usr/bin/env Rscript
# ============================================================
# Script 10: ENC-GC3s deviation statistics
# ------------------------------------------------------------
# Purpose: Quantify the deviation of observed ENC values from the
#   expected Wright (1990) curve, per segment. The expected ENC uses
#   the original Wright (1990) formula ENC = 2 + s + 29/[s^2 + (1-s)^2]
#   (wright_expected_enc() in common_bat.R), as drawn in the figures.
#
#   For each segment we report:
#     - mean observed ENC and mean expected ENC (at each sequence's GC3s)
#     - mean deviation (observed - expected)
#     - percentage of sequences falling below the expected curve
#     - Wilcoxon signed-rank test (observed vs expected, paired)
#     - Spearman correlation between GC3s and observed ENC
#
# Input:  output/tables/gc_enc_values_{seg}.csv  (columns: accession, GC, GC3s, ENC, ...)
# Output: output/tables/enc_deviation_stats.csv
# ============================================================

source("c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/scripts/common_bat.R")
dirs <- setup_dirs()

SEGMENTS <- c("VP4", "VP6", "VP7", "NSP4")

results <- lapply(SEGMENTS, function(seg) {
  df <- read.csv(file.path(dirs$tabdir, paste0("gc_enc_values_", seg, ".csv")),
                 stringsAsFactors = FALSE)
  # Expected ENC under the original Wright (1990) null curve,
  # identical to the curve drawn in the neutrality figures (script 07):
  #   ENC_exp = 2 + s + 29 / [s^2 + (1-s)^2],  s = GC3s fraction (capped at 61)
  enc_exp <- pmin(wright_expected_enc(df$GC3s / 100), 61)
  enc_obs <- df$ENC

  diff <- enc_obs - enc_exp
  wil <- suppressWarnings(wilcox.test(enc_obs, enc_exp, paired = TRUE, exact = FALSE))
  sp  <- suppressWarnings(cor.test(df$GC3s, enc_obs, method = "spearman", exact = FALSE))

  data.frame(
    segment        = seg,
    n              = nrow(df),
    mean_enc_obs   = round(mean(enc_obs), 2),
    mean_enc_exp   = round(mean(enc_exp), 2),
    mean_diff      = round(mean(diff), 2),
    sd_diff        = round(sd(diff), 2),
    pct_below      = round(100 * mean(diff < 0), 1),
    wilcoxon_p     = signif(wil$p.value, 3),
    spearman_rho   = round(unname(sp$estimate), 3),
    spearman_p     = signif(sp$p.value, 3),
    stringsAsFactors = FALSE
  )
})

stats <- do.call(rbind, results)
write.csv(stats, file.path(dirs$tabdir, "enc_deviation_stats.csv"), row.names = FALSE)
cat("Saved: enc_deviation_stats.csv\n\n")
print(stats, row.names = FALSE)
cat("\n=== ENC deviation analysis complete ===\n")
