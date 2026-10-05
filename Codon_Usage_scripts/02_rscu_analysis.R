# ============================================================
# 02_rscu_analysis.R
# RSCU analysis: per-sequence RSCU, Euclidean distance, UPGMA clustering,
# heatmap + dendrogram per segment (VP4, VP6, VP7, NSP4)
# ============================================================

source("c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/scripts/common_bat.R")
dirs <- setup_dirs()

# ---- Segment-specific CDS length thresholds ----
LEN_THRESH <- c(VP4=2300, VP6=1150, VP7=950, NSP4=500)
SEGMENTS <- c("VP4", "VP6", "VP7", "NSP4")

# ---- Load and filter sequences per segment ----
load_segment_data <- function(segment) {
  fasta_path <- file.path(dirs$fastadir, paste0(segment, "_cds.fasta"))
  meta_path  <- file.path(dirs$metadir, paste0(segment, "_metadata_parsed.csv"))
  
  # Load FASTA
  seqs <- seqinr::read.fasta(fasta_path, as.string = TRUE, forceDNAtolower = FALSE)
  seq_ids <- names(seqs)
  
  # Load metadata
  meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  
  # Match sequences to metadata
  seq_list <- list()
  meta_list <- list()
  
  for (sid in seq_ids) {
    seq <- toupper(seqs[[sid]])
    # Remove gaps
    seq <- gsub("-", "", seq)
    n <- nchar(seq)
    
    # Filter by length threshold
    if (n < LEN_THRESH[segment]) next
    
    # Find metadata
    m <- meta[meta$accession == sid, ]
    if (nrow(m) == 0) {
      # Try matching by short name (for study sequences)
      short <- strsplit(sid, "_NODE_")[[1]][1]
      m <- meta[meta$accession == short, ]
    }
    if (nrow(m) == 0) next
    
    m <- m[1, ]  # Take first match
    
    seq_list[[sid]] <- seq
    meta_list[[sid]] <- data.frame(
      accession = sid,
      full_name = m$full_name,
      host = m$host,
      GC = m$GC,
      source = m$source,
      segment = segment,
      length = n,
      stringsAsFactors = FALSE
    )
  }
  
  meta_df <- do.call(rbind, meta_list)
  cat(sprintf("  %s: %d sequences (after length filter >= %d bp)\n",
              segment, length(seq_list), LEN_THRESH[segment]))
  
  list(seqs = seq_list, meta = meta_df)
}

# ---- Compute RSCU per sequence ----
compute_rscu_matrix <- function(seq_list) {
  # RSCU via seqinr::uco(index = "rscu") — standard implementation
  rscu_mat <- t(sapply(names(seq_list), function(sid) {
    std_rscu(seq_list[[sid]])
  }))
  # Keep only 59 synonymous codons
  rscu_mat <- rscu_mat[, SYN]
  rownames(rscu_mat) <- names(seq_list)
  rscu_mat
}

# ---- RSCU heatmap per segment ----
plot_rscu_heatmap <- function(rscu_mat, meta_df, segment, outdir) {
  # Order rows by constellation group
  ord <- order(meta_df$GC, meta_df$host)
  rscu_ordered <- rscu_mat[ord, ]
  
  # Column annotations (sequences are columns after transpose)
  col_anno <- HeatmapAnnotation(
    Group = meta_df$GC[ord],
    Host = meta_df$host[ord],
    col = list(
      Group = group_colors[unique(meta_df$GC[ord])],
      Host = setNames(rainbow(length(unique(meta_df$host[ord]))),
                      unique(meta_df$host[ord]))
    ),
    annotation_legend_param = list(
      Group = list(title = "Constellation"),
      Host = list(title = "Host species")
    )
  )
  
  # Column order by amino acid family
  aa_order <- unlist(FAM[order(sapply(FAM, length), decreasing = TRUE)])
  col_order <- intersect(aa_order, colnames(rscu_ordered))
  
  png(file.path(outdir, paste0("rscu_heatmap_", segment, ".png")),
      width = 14, height = 10, units = "in", res = 300)
  draw(
    Heatmap(t(rscu_ordered[, col_order]),
            name = "RSCU",
            col = colorRamp2(c(0, 1, 2, 3), c("#2166AC", "#F7F7F7", "#FDDBC7", "#B2182B")),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_column_names = FALSE,
            row_names_gp = gpar(fontsize = 7),
            column_title = paste0(segment, " — RSCU Heatmap"),
            top_annotation = col_anno,
            use_raster = TRUE,
            raster_quality = 3)
  )
  dev.off()
  cat(sprintf("  Saved: rscu_heatmap_%s.png\n", segment))
}

# ---- RSCU dendrogram per segment ----
plot_rscu_dendrogram <- function(rscu_mat, meta_df, segment, outdir) {
  # Euclidean distance on 59 RSCU values
  dist_mat <- dist(rscu_mat, method = "euclidean")
  
  # UPGMA clustering (average linkage)
  hc <- hclust(dist_mat, method = "average")
  dend <- as.dendrogram(hc)
  
  # Color tips by constellation group
  tip_colors <- group_colors[meta_df$GC[order.dendrogram(dend)]]
  labels_colors(dend) <- tip_colors
  
  png(file.path(outdir, paste0("rscu_dendrogram_", segment, ".png")),
      width = 12, height = 8, units = "in", res = 300)
  par(mar = c(7, 4, 4, 2))
  plot(dend,
       main = paste0(segment, " — RSCU UPGMA Dendrogram"),
       ylab = "Euclidean Distance (59 synonymous codons)",
       cex = 0.5,
       las = 2)
  legend("topright",
         legend = names(group_colors)[names(group_colors) %in% unique(meta_df$GC)],
         fill = group_colors[names(group_colors) %in% unique(meta_df$GC)],
         title = "Constellation",
         cex = 0.7,
         ncol = 2)
  dev.off()
  cat(sprintf("  Saved: rscu_dendrogram_%s.png\n", segment))
  
  # Return the hclust object for tanglegram
  list(hc = hc, dist_mat = dist_mat, dend = dend)
}

# ---- Main analysis loop ----
all_rscu_results <- list()
all_rscu_matrices <- list()
all_meta <- list()

for (seg in SEGMENTS) {
  cat(sprintf("\n=== %s ===\n", seg))
  
  # Load data
  data <- load_segment_data(seg)
  all_meta[[seg]] <- data$meta
  
  # Compute RSCU matrix
  rscu_mat <- compute_rscu_matrix(data$seqs)
  all_rscu_matrices[[seg]] <- rscu_mat
  
  # Save RSCU values
  rscu_df <- as.data.frame(rscu_mat)
  rscu_df$accession <- rownames(rscu_df)
  write.csv(rscu_df, file.path(dirs$tabdir, paste0("rscu_per_sequence_", seg, ".csv")),
            row.names = FALSE)
  
  # Save distance matrix
  dist_mat <- dist(rscu_mat, method = "euclidean")
  write.csv(as.matrix(dist_mat),
            file.path(dirs$tabdir, paste0("rscu_distance_matrix_", seg, ".csv")))
  
  # Plot heatmap
  plot_rscu_heatmap(rscu_mat, data$meta, seg, dirs$figdir)
  
  # Plot dendrogram
  dend_result <- plot_rscu_dendrogram(rscu_mat, data$meta, seg, dirs$figdir)
  all_rscu_results[[seg]] <- dend_result
  
  cat(sprintf("  %s: %d sequences, %d codons analyzed\n", seg, nrow(rscu_mat), ncol(rscu_mat)))
}

# ---- Build host × segment summary table ----
cat("\n=== Host × Segment Summary ===\n")
summary_list <- list()
for (seg in SEGMENTS) {
  m <- all_meta[[seg]]
  agg <- m %>%
    group_by(host, GC) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(segment = seg)
  summary_list[[seg]] <- agg
}
host_segment_summary <- do.call(rbind, summary_list) %>%
  pivot_wider(names_from = segment, values_from = n, values_fill = 0)

write.csv(host_segment_summary,
          file.path(dirs$tabdir, "host_segment_summary.csv"),
          row.names = FALSE)
cat("  Saved: host_segment_summary.csv\n")
print(host_segment_summary)

# ---- Save combined metadata for downstream scripts ----
combined_meta <- do.call(rbind, all_meta)
write.csv(combined_meta,
          file.path(dirs$tabdir, "combined_metadata.csv"),
          row.names = FALSE)
cat("\nSaved combined metadata for all segments\n")

cat("\n=== RSCU Analysis Complete ===\n")
