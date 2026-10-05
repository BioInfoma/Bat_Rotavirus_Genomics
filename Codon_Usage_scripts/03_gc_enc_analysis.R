# ============================================================
# 03_gc_enc_analysis.R
# GC composition (overall GC%, GC1, GC2, GC3, GC3s) + ENC (Wright 1990)
# + ENC-GC3s neutrality plot per segment
# ============================================================

source("c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/scripts/common_bat.R")
dirs <- setup_dirs()

SEGMENTS <- c("VP4", "VP6", "VP7", "NSP4")
LEN_THRESH <- c(VP4=2300, VP6=1150, VP7=950, NSP4=500)

# ---- Load sequences per segment (same as 02) ----
load_segment_data <- function(segment) {
  fasta_path <- file.path(dirs$fastadir, paste0(segment, "_cds.fasta"))
  meta_path  <- file.path(dirs$metadir, paste0(segment, "_metadata_parsed.csv"))
  seqs <- seqinr::read.fasta(fasta_path, as.string = TRUE, forceDNAtolower = FALSE)
  meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  
  seq_list <- list()
  meta_list <- list()
  for (sid in names(seqs)) {
    seq <- toupper(seqs[[sid]])
    seq <- gsub("-", "", seq)
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
    meta_list[[sid]] <- data.frame(accession=sid, host=m$host, GC=m$GC,
                                    segment=segment, length=n, stringsAsFactors=FALSE)
  }
  list(seqs=seq_list, meta=do.call(rbind, meta_list))
}

# ---- Compute GC composition and ENC per sequence ----
compute_gc_enc <- function(seq_list) {
  # Composition via seqinr (GC/GC1/GC2/GC3) + documented synonymous-site
  # derivations from standard uco codon counts
  results <- lapply(names(seq_list), function(sid) {
    bc <- std_composition(seq_list[[sid]])
    cbind(accession=sid, bc)
  })
  df <- do.call(rbind, results)
  # ENC via coRdon (Wright 1990) — batch, coRdon needs a multi-row table
  cnt_mat <- t(sapply(names(seq_list), function(sid) std_counts(seq_list[[sid]])))
  df$ENC <- std_enc(cnt_mat)
  df
}

# ---- ENC-GC3s neutrality plot ----
plot_enc_gc3 <- function(gc_enc_df, meta_df, segment, outdir) {
  df <- merge(gc_enc_df, meta_df[, c("accession", "host", "GC")], by="accession")
  df$GC3s <- as.numeric(df$GC3s)
  df$ENC <- as.numeric(df$ENC)
  
  # Expected ENC curve: original Wright (1990) formula
  #   ENC = 2 + s + 29 / [s^2 + (1-s)^2],  s = GC3s as a fraction
  gc3s_range <- seq(0, 1, by=0.01)
  enc_expected <- pmin(wright_expected_enc(gc3s_range), 61)
  
  expected_df <- data.frame(GC3s = gc3s_range * 100, ENC = enc_expected)
  
  p <- ggplot(df, aes(x=GC3s, y=ENC)) +
    # Expected curve
    geom_line(data=expected_df, aes(x=GC3s, y=ENC), color="gray50", linewidth=0.8, linetype="dashed") +
    # Data points colored by constellation group
    geom_point(aes(color=GC), size=2.5, alpha=0.8) +
    scale_color_manual(values=group_colors, name="Constellation") +
    labs(title=paste0(segment, " — ENC vs GC3s Neutrality Plot"),
         x="GC3s (%)", y="Effective Number of Codons (ENC)") +
    theme_cub() +
    theme(legend.position="right") +
    coord_cartesian(ylim=c(20, 61))
  
  png(file.path(outdir, paste0("enc_gc3_neutrality_", segment, ".png")),
      width=8, height=6, units="in", res=300)
  print(p)
  dev.off()
  cat(sprintf("  Saved: enc_gc3_neutrality_%s.png\n", segment))
}

# ---- GC3 by group boxplot ----
plot_gc3_by_group <- function(gc_enc_df, meta_df, segment, outdir) {
  df <- merge(gc_enc_df, meta_df[, c("accession", "host", "GC")], by="accession")
  df$GC3 <- as.numeric(df$GC3)
  df$GC <- factor(df$GC, levels=names(group_colors)[names(group_colors) %in% unique(df$GC)])
  
  p <- ggplot(df, aes(x=GC, y=GC3, fill=GC)) +
    geom_boxplot(alpha=0.7) +
    geom_jitter(width=0.15, size=1.5, alpha=0.6) +
    scale_fill_manual(values=group_colors) +
    labs(title=paste0(segment, " — GC3 by Constellation Group"),
         x="Constellation Group", y="GC3 (%)") +
    theme_cub() +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))
  
  png(file.path(outdir, paste0("gc3_by_group_", segment, ".png")),
      width=8, height=5, units="in", res=300)
  print(p)
  dev.off()
  cat(sprintf("  Saved: gc3_by_group_%s.png\n", segment))
}

# ---- Main analysis loop ----
for (seg in SEGMENTS) {
  cat(sprintf("\n=== %s ===\n", seg))
  data <- load_segment_data(seg)
  gc_enc <- compute_gc_enc(data$seqs)
  
  # Convert numeric columns
  for (cn in names(gc_enc)) {
    if (cn != "accession") gc_enc[[cn]] <- as.numeric(gc_enc[[cn]])
  }
  
  # Save GC + ENC values
  write.csv(gc_enc, file.path(dirs$tabdir, paste0("gc_enc_values_", seg, ".csv")),
            row.names=FALSE)
  
  # Plots
  plot_enc_gc3(gc_enc, data$meta, seg, dirs$figdir)
  plot_gc3_by_group(gc_enc, data$meta, seg, dirs$figdir)
  
  cat(sprintf("  %s: %d sequences analyzed\n", seg, nrow(gc_enc)))
}

cat("\n=== GC/ENC Analysis Complete ===\n")
