# ============================================================
# 08_coa_pr2_composites.R
# 1. Correspondence analysis (COA) of RSCU per segment (ade4::dudi.coa)
# 2. PR2 parity bias plots per segment (4-fold degenerate 3rd positions)
# 3. Multi-panel composite figures assembled from saved v2 PNGs
# ============================================================

source("c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/scripts/common_bat.R")
suppressPackageStartupMessages({ library(ade4); library(cowplot) })
dirs <- setup_dirs()
options(figdir = dirs$figdir)

SEGMENTS  <- c("VP4", "VP6", "VP7", "NSP4")
LEN_THRESH <- c(VP4=2300, VP6=1150, VP7=950, NSP4=500)

# ---- Load segment sequences + metadata (same logic as 04) ----
load_segment_data <- function(segment) {
  fasta_path <- file.path(dirs$fastadir, paste0(segment, "_cds.fasta"))
  meta_path  <- file.path(dirs$metadir, paste0(segment, "_metadata_parsed.csv"))
  seqs <- seqinr::read.fasta(fasta_path, as.string = TRUE, forceDNAtolower = FALSE)
  meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  seq_list <- list(); meta_list <- list()
  for (sid in names(seqs)) {
    seq <- toupper(seqs[[sid]]); seq <- gsub("-", "", seq)
    n <- nchar(seq)
    if (n < LEN_THRESH[segment]) next
    m <- meta[meta$accession == sid, ]
    if (nrow(m) == 0) {
      short <- strsplit(sid, "_NODE_")[[1]][1]
      m <- meta[meta$accession == short, ]
    }
    if (nrow(m) == 0) next
    m <- m[1, ]
    seq_list[[sid]] <- seq
    meta_list[[sid]] <- data.frame(accession=sid, full_name=m$full_name,
                                   host=m$host, GC=m$GC, source=m$source,
                                   segment=segment, length=n, stringsAsFactors=FALSE)
  }
  list(seqs=seq_list, meta=do.call(rbind, meta_list))
}

# ============================================================
# PART 1 — COA of RSCU
# ============================================================
coa_one_segment <- function(segment) {
  rscu_raw <- read.csv(file.path(dirs$tabdir, paste0("rscu_per_sequence_", segment, ".csv")),
                       check.names = FALSE)
  meta <- read.csv(file.path(dirs$metadir, paste0(segment, "_metadata_parsed.csv")),
                   stringsAsFactors = FALSE)
  # RSCU table has an accession column + 59 codon columns; merge on accession
  meta <- merge(meta, rscu_raw[, "accession", drop = FALSE], by = "accession", sort = FALSE)
  rscu <- rscu_raw[match(meta$accession, rscu_raw$accession), ]
  rownames(rscu) <- rscu$accession
  rscu <- rscu[, setdiff(colnames(rscu), "accession")]
  rscu <- as.data.frame(lapply(rscu, as.numeric))
  stopifnot(nrow(rscu) == nrow(meta))

  # Drop codons with zero variance or all-zero (uninformative / undefined in CA)
  keep <- colSums(rscu) > 0 & apply(rscu, 2, sd) > 0
  rscu_f <- rscu[, keep]
  cat(sprintf("  %s: COA on %d sequences x %d codons (dropped %d invariant)\n",
              segment, nrow(rscu_f), ncol(rscu_f), sum(!keep)))

  coa <- dudi.coa(rscu_f, scannf = FALSE, nf = 3)
  inertia <- coa$eig / sum(coa$eig) * 100

  coords <- data.frame(accession = meta$accession,
                       host = meta$host, GC = meta$GC,
                       Ax1 = coa$li[, 1], Ax2 = coa$li[, 2],
                       stringsAsFactors = FALSE)
  write.csv(coords, file.path(dirs$tabdir, paste0("coa_rscu_coords_", segment, ".csv")),
            row.names = FALSE)

  ax_lab <- function(i) sprintf("Axis %d (%.1f%%)", i, inertia[i])

  # Panel A: colored by constellation group
  pa <- ggplot(coords, aes(Ax1, Ax2, color = GC)) +
    geom_point(size = 4.2, alpha = 0.85) +
    scale_color_manual(values = group_colors, name = "Constellation") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.3) +
    labs(title = paste0(segment, " — by constellation group"),
         x = ax_lab(1), y = ax_lab(2)) +
    theme_cub() +
    theme(legend.position = "right")

  # Panel B: colored by host species
  pb <- ggplot(coords, aes(Ax1, Ax2, color = host)) +
    geom_point(size = 4.2, alpha = 0.85) +
    scale_color_manual(values = host_colors, name = "Host") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.3) +
    labs(title = paste0(segment, " — by host species"),
         x = ax_lab(1), y = ax_lab(2)) +
    theme_cub() +
    theme(legend.position = "right",
          legend.text = element_text(size = 9.5, face = "italic"))

  # No inner panel tags: composites add their own A-D labels
  p <- (pa | pb)
  save_fig(p, paste0("coa_rscu_", segment, "_v2"), width = 13, height = 5.5)

  list(coords = coords, inertia = inertia[1:3])
}

# ============================================================
# PART 2 — PR2 parity bias plots (4-fold degenerate 3rd positions)
# ============================================================
# 4-fold degenerate codon boxes (third position fully synonymous)
FF_BOXES <- c("CT", "TC", "CG", "GC", "CC", "AC", "GT", "GG")

# PR2 coordinates derived from standard seqinr::uco codon counts:
# third-base frequencies across the eight four-fold degenerate boxes
# (incl. the four-fold blocks of Leu/Ser/Arg), per Sueoka (1995).
pr2_coords_seq <- function(seq) {
  cnt <- std_counts(seq)
  ff_codons <- CODONS64[substr(CODONS64, 1, 2) %in% FF_BOXES]
  ff <- cnt[ff_codons]
  n_ff <- sum(ff)
  if (n_ff < 10) return(NULL)
  third <- substr(ff_codons, 3, 3)
  A3 <- sum(ff[third == "A"]) / n_ff; T3 <- sum(ff[third == "T"]) / n_ff
  G3 <- sum(ff[third == "G"]) / n_ff; C3 <- sum(ff[third == "C"]) / n_ff
  data.frame(AT_bias = A3 / (A3 + T3), GC_bias = G3 / (G3 + C3),
             n_4fold = n_ff)
}

pr2_one_segment <- function(segment) {
  data <- load_segment_data(segment)
  rows <- lapply(names(data$seqs), function(sid) {
    r <- pr2_coords_seq(data$seqs[[sid]])
    if (is.null(r)) return(NULL)
    cbind(data.frame(accession = sid), r)
  })
  df <- do.call(rbind, rows)
  df <- merge(df, data$meta[, c("accession", "host", "GC")], by = "accession")
  write.csv(df, file.path(dirs$tabdir, paste0("pr2_coords_", segment, ".csv")),
            row.names = FALSE)

  # Distance from parity center (0.5, 0.5) as a bias magnitude summary
  df$dist_from_center <- sqrt((df$AT_bias - 0.5)^2 + (df$GC_bias - 0.5)^2)
  cat(sprintf("  %s: PR2 for %d sequences; mean dist from parity = %.3f\n",
              segment, nrow(df), mean(df$dist_from_center)))

  p <- ggplot(df, aes(GC_bias, AT_bias, color = GC)) +
    geom_point(size = 4.2, alpha = 0.85) +
    scale_color_manual(values = group_colors, name = "Constellation") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.4) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.4) +
    annotate("point", x = 0.5, y = 0.5, shape = 4, size = 5.0, stroke = 1.2) +
    labs(title = paste0(segment, " — PR2 bias plot"),
         x = "G3 / (G3 + C3)", y = "A3 / (A3 + T3)") +
    theme_cub() +
    theme(legend.position = "right")
  save_fig(p, paste0("pr2_bias_", segment, "_v2"), width = 7, height = 5.5)
  df
}

# ============================================================
# PART 3 — Composite figures from saved v2 PNGs
# ============================================================
composite_2x2 <- function(png_paths, out_name, labels = c("A", "B", "C", "D"),
                          width = 14, height = 12) {
  # Top/left margins give the A-D labels a collision-free zone above panel titles
  panels <- lapply(seq_along(png_paths), function(i) {
    ggdraw() + draw_image(png_paths[i]) +
      theme(plot.margin = margin(20, 8, 8, 30))
  })
  comp <- plot_grid(plotlist = panels, ncol = 2, labels = labels,
                    label_size = 16, label_fontface = "bold",
                    label_x = 0.02, label_y = 0.995, hjust = 0, vjust = 1)
  save_fig(comp, out_name, width = width, height = height)
}

# ============================================================
# MAIN
# ============================================================
cat("=== PART 1: COA of RSCU ===\n")
coa_results <- list()
for (seg in SEGMENTS) {
  cat(sprintf("\n--- %s ---\n", seg))
  coa_results[[seg]] <- coa_one_segment(seg)
}

cat("\n=== PART 2: PR2 bias plots ===\n")
pr2_results <- list()
for (seg in SEGMENTS) {
  cat(sprintf("\n--- %s ---\n", seg))
  pr2_results[[seg]] <- pr2_one_segment(seg)
}

cat("\n=== PART 3: Composite figures ===\n")
fig <- function(pattern, seg) file.path(dirs$figdir, sprintf(pattern, seg))

# Comp 1: RSCU heatmaps
composite_2x2(sapply(SEGMENTS, function(s) fig("rscu_heatmap_%s_v2.png", s)),
              "composite_rscu_heatmaps", width = 16, height = 13)
# Comp 2: RSCU dendrograms
composite_2x2(sapply(SEGMENTS, function(s) fig("rscu_dendrogram_combined_%s_v2.png", s)),
              "composite_rscu_dendrograms", width = 16, height = 13)
# Comp 3: ENC-GC3s neutrality
composite_2x2(sapply(SEGMENTS, function(s) fig("enc_gc3_neutrality_%s_v2.png", s)),
              "composite_enc_gc3s", width = 14, height = 11)
# Comp 4: COA
composite_2x2(sapply(SEGMENTS, function(s) fig("coa_rscu_%s_v2.png", s)),
              "composite_coa_rscu", width = 18, height = 9)
# Comp 5: PR2
composite_2x2(sapply(SEGMENTS, function(s) fig("pr2_bias_%s_v2.png", s)),
              "composite_pr2_bias", width = 14, height = 11)
# Comp 6: tanglegrams
composite_2x2(sapply(SEGMENTS, function(s) fig("tanglegram_ggtree_%s_v2.png", s)),
              "composite_tanglegrams", width = 28, height = 24)
# Comp 7: GC3 boxplots
composite_2x2(sapply(SEGMENTS, function(s) fig("gc3_by_group_%s_v2.png", s)),
              "composite_gc3_by_group", width = 14, height = 11)

# ---- Summary table of COA inertia + PR2 bias ----
summary_rows <- lapply(SEGMENTS, function(seg) {
  data.frame(segment = seg,
             coa_axis1_pct = round(coa_results[[seg]]$inertia[1], 2),
             coa_axis2_pct = round(coa_results[[seg]]$inertia[2], 2),
             pr2_mean_dist = round(mean(pr2_results[[seg]]$dist_from_center), 4),
             n_seq = nrow(pr2_results[[seg]]))
})
coa_pr2_summary <- do.call(rbind, summary_rows)
write.csv(coa_pr2_summary, file.path(dirs$tabdir, "coa_pr2_summary.csv"), row.names = FALSE)
print(coa_pr2_summary)

cat("\n=== Script 08 complete ===\n")
