# ============================================================
# 04_cai_analysis.R
# CAI plots and Kruskal-Wallis tests per segment.
# CAI values are computed by 04a_cai_biopython.py (Biopython
# CodonAdaptationIndex; Sharp & Li 1987) — this script only reads
# the resulting cai_values_{seg}.csv files and plots/tests them.
# ============================================================

source("c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/scripts/common_bat.R")
dirs <- setup_dirs()

SEGMENTS <- c("VP4", "VP6", "VP7", "NSP4")
LEN_THRESH <- c(VP4=2300, VP6=1150, VP7=950, NSP4=500)
CODON_DIR <- "c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/host_codon_tables"

# ---- Host-to-codon-table mapping ----
host_map <- host_table_map()

# ---- Load sequences per segment ----
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
    meta_list[[sid]] <- data.frame(accession=sid, full_name=m$full_name,
                                    host=m$host, GC=m$GC, source=m$source,
                                    segment=segment, length=n, stringsAsFactors=FALSE)
  }
  list(seqs=seq_list, meta=do.call(rbind, meta_list))
}

# (CAI computation moved to 04a_cai_biopython.py — Biopython)


# ---- CAI by group boxplot ----
plot_cai_by_group <- function(cai_df, meta_df, segment, outdir) {
  df <- merge(cai_df, meta_df[, c("accession", "GC")], by="accession")
  df <- df[!is.na(df$CAI), ]
  df$GC <- factor(df$GC, levels=names(group_colors)[names(group_colors) %in% unique(df$GC)])

  # Kruskal-Wallis across groups (only if >=2 groups with >=2 obs)
  kw_p <- NA
  tab <- table(df$GC)
  if (sum(tab >= 2) >= 2) {
    kw_p <- kruskal.test(CAI ~ GC, data=df)$p.value
  }
  sub <- ifelse(is.na(kw_p), "insufficient groups for test",
                sprintf("Kruskal-Wallis p = %.3g", kw_p))

  # Sample sizes under labels
  nlab <- df %>% count(GC) %>% mutate(lab=paste0(GC, "\nn=", n))
  df$GC <- factor(df$GC, levels=nlab$GC, labels=nlab$lab)

  p <- ggplot(df, aes(x=GC, y=CAI, fill=GC)) +
    geom_boxplot(alpha=0.7, outlier.shape=NA) +
    geom_jitter(width=0.15, size=1.5, alpha=0.6) +
    scale_fill_manual(values=setNames(group_colors[as.character(nlab$GC)], nlab$lab)) +
    labs(title=paste0(segment, " — CAI by constellation group"), subtitle=sub,
         x="Constellation Group", y="CAI") +
    theme_cub() +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))

  save_fig(p, paste0("cai_by_group_", segment, "_v2"), width=8, height=5)
  invisible(kw_p)
}

# ---- CAI by host boxplot ----
plot_cai_by_host <- function(cai_df, meta_df, segment, outdir) {
  # cai_df already carries 'host'; only take GC from metadata to avoid .x/.y suffixes
  df <- merge(cai_df, meta_df[, c("accession", "GC")], by="accession")
  df <- df[!is.na(df$CAI), ]
  # Order hosts by median CAI
  host_order <- names(sort(tapply(df$CAI, df$host, median)))
  df$host <- factor(df$host, levels=host_order)

  # Kruskal-Wallis across hosts (only if >=2 hosts with >=2 obs)
  kw_p <- NA
  tab <- table(df$host)
  if (sum(tab >= 2) >= 2) {
    kw_p <- kruskal.test(CAI ~ host, data=df)$p.value
  }
  sub <- ifelse(is.na(kw_p), "insufficient replication for test",
                sprintf("Kruskal-Wallis p = %.3g", kw_p))

  p <- ggplot(df, aes(x=host, y=CAI, fill=GC)) +
    geom_boxplot(alpha=0.7, outlier.shape=NA) +
    geom_jitter(width=0.15, size=1.5, alpha=0.6) +
    scale_fill_manual(values=group_colors, name="Constellation") +
    labs(title=paste0(segment, " — CAI by host species"), subtitle=sub,
         x="Host Species", y="CAI") +
    theme_cub() +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8, face="italic"))

  save_fig(p, paste0("cai_by_host_", segment, "_v2"), width=10, height=6)
  invisible(kw_p)
}

# ---- Main analysis loop ----
kw_results <- list()
for (seg in SEGMENTS) {
  cat(sprintf("\n=== %s ===\n", seg))
  data <- load_segment_data(seg)

  # Read CAI values computed by 04a_cai_biopython.py (Biopython)
  cai_full <- read.csv(file.path(dirs$tabdir, paste0("cai_values_", seg, ".csv")),
                       stringsAsFactors=FALSE)
  cai_df <- cai_full[, c("accession", "host", "CAI")]

  # Plots (return Kruskal-Wallis p-values)
  kw_group <- plot_cai_by_group(cai_df, data$meta, seg, dirs$figdir)
  kw_host  <- plot_cai_by_host(cai_df, data$meta, seg, dirs$figdir)
  kw_results[[seg]] <- data.frame(segment=seg,
                                  kw_group_p=kw_group, kw_host_p=kw_host)

  # Summary
  valid <- cai_full[!is.na(cai_full$CAI), ]
  cat(sprintf("  %s: %d sequences, %d with CAI values\n", seg, nrow(cai_full), nrow(valid)))
  cat(sprintf("  CAI range: %.4f - %.4f (median: %.4f)\n",
              min(valid$CAI), max(valid$CAI), median(valid$CAI)))
}

# ---- Save KW test results ----
kw_df <- do.call(rbind, kw_results)
write.csv(kw_df, file.path(dirs$tabdir, "cai_kruskal_wallis.csv"), row.names=FALSE)
cat("\nKruskal-Wallis results:\n"); print(kw_df)

# ---- Composite CAI-by-group figure (2x2) ----
suppressPackageStartupMessages(library(cowplot))
panels <- lapply(SEGMENTS, function(s)
  ggdraw() + draw_image(file.path(dirs$figdir, paste0("cai_by_group_", s, "_v2.png"))) +
    theme(plot.margin = margin(20, 8, 8, 30)))
comp <- plot_grid(plotlist=panels, ncol=2, labels=c("A","B","C","D"),
                  label_size=16, label_fontface="bold",
                  label_x=0.02, label_y=0.995, hjust=0, vjust=1)
save_fig(comp, "composite_cai_by_group", width=14, height=10)

cat("\n=== CAI Analysis Complete ===\n")
