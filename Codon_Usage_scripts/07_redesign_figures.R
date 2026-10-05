# ============================================================
# 07_redesign_figures.R
# Publication-quality figure redesign 
# Regenerates all figures with aesthetics
# ============================================================

source("c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/scripts/common_bat.R")
dirs <- setup_dirs()
options(figdir = dirs$figdir)

SEGMENTS <- c("VP4", "VP6", "VP7", "NSP4")
LEN_THRESH <- c(VP4=2300, VP6=1150, VP7=950, NSP4=500)
# De novo trees built in this study (MAFFT L-INS-i + IQ-TREE MFP, 1000 UFBoot)
TREE_DIR <- dirs$treedir

TREE_FILES <- list(
  VP4  = file.path(TREE_DIR, "VP4_ml.treefile"),
  VP6  = file.path(TREE_DIR, "VP6_ml.treefile"),
  VP7  = file.path(TREE_DIR, "VP7_ml.treefile"),
  NSP4 = file.path(TREE_DIR, "NSP4_ml.treefile")
)

# ---- Data loading functions ----
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
    if (nrow(m) == 0) { short <- strsplit(sid, "_NODE_")[[1]][1]; m <- meta[meta$accession == short, ] }
    if (nrow(m) == 0) next
    m <- m[1, ]
    seq_list[[sid]] <- seq
    meta_list[[sid]] <- data.frame(accession=sid, full_name=m$full_name, host=m$host,
                                    GC=m$GC, source=m$source, segment=segment, length=n,
                                    stringsAsFactors=FALSE)
  }
  list(seqs=seq_list, meta=do.call(rbind, meta_list))
}

compute_rscu_matrix <- function(seq_list) {
  rscu_mat <- t(sapply(names(seq_list), function(sid) {
    std_rscu(seq_list[[sid]])   # seqinr::uco(index = "rscu")
  }))
  rscu_mat <- rscu_mat[, SYN]; rownames(rscu_mat) <- names(seq_list); rscu_mat
}

# ---- Short label for sequences ----
short_label <- function(sid, meta_df) {
  m <- meta_df[meta_df$accession == sid, ]
  if (nrow(m) > 0 && tolower(m$source[1]) == "study") {
    # Study sequence: use short name
    parts <- strsplit(sid, "_NODE_")[[1]]
    return(parts[1])
  }
  # Reference: accession.version
  sid
}

# ============================================================
# FIGURE 1: RSCU Heatmaps (ComplexHeatmap, AA family grouped)
# ============================================================
make_rscu_heatmap <- function(rscu_mat, meta_df, segment) {
  # UPGMA clustering of sequences (columns)
  dist_mat <- dist(rscu_mat, method = "euclidean")
  hc <- hclust(dist_mat, method = "average")
  col_order <- hc$order

  # Order by clustering
  rscu_ordered <- rscu_mat[col_order, ]

  # AA family grouping for rows
  aa_for_codon <- sapply(colnames(rscu_ordered), function(cn) {
    for (aa in names(FAM)) { if (cn %in% FAM[[aa]]) return(aa) }
    NA
  })

  # Top annotation: constellation group (for columns = sequences)
  group_anno <- HeatmapAnnotation(
    Group = meta_df$GC[col_order],
    col = list(Group = group_colors[unique(meta_df$GC[col_order])]),
    annotation_name_gp = gpar(fontsize = 11.0, fontface = "bold"),
    annotation_legend_param = list(Group = list(title = "Constellation", title_gp = gpar(fontsize = 12.0, fontface = "bold"), labels_gp = gpar(fontsize = 11.0)))
  )

  # Bottom annotation: host species (for columns = sequences)
  host_anno <- HeatmapAnnotation(
    Host = meta_df$host[col_order],
    col = list(Host = host_colors[unique(meta_df$host[col_order])]),
    annotation_name_gp = gpar(fontsize = 11.0, fontface = "bold"),
    annotation_legend_param = list(Host = list(title = "Host species", title_gp = gpar(fontsize = 12.0, fontface = "bold"), labels_gp = gpar(fontsize = 10.0)))
  )

  # Row split by AA family
  row_split_factor <- aa_for_codon[rownames(t(rscu_ordered))]

  # Color scale: diverging centered at 1.0
  col_fun <- colorRamp2(c(0, 0.5, 1.0, 1.5, 2.0, 3.0),
                        c("#2166AC", "#67A9CF", "#F7F7F7", "#F4A582", "#D6604D", "#B2182B"))

  # Column dendrogram
  col_dend <- as.dendrogram(hc)

  # Build heatmap (codons = rows, sequences = columns)
  ht <- Heatmap(
    t(rscu_ordered),
    name = "RSCU",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = col_dend,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 9.0),
    row_split = row_split_factor,
    row_title_gp = gpar(fontsize = 10.0, fontface = "bold"),
    row_title_rot = 0,
    column_title = paste0(segment, " — RSCU Codon Usage Heatmap"),
    column_title_gp = gpar(fontsize = 14.0, fontface = "bold"),
    top_annotation = group_anno,
    bottom_annotation = host_anno,
    use_raster = TRUE,
    raster_quality = 4,
    heatmap_legend_param = list(
      title = "RSCU",
      title_gp = gpar(fontsize = 12.0, fontface = "bold"),
      labels_gp = gpar(fontsize = 11.0),
      legend_height = unit(3, "cm")
    )
  )

  # Save PNG
  png_path <- file.path(dirs$figdir, paste0("rscu_heatmap_", segment, "_v2.png"))
  png(png_path, width = 10, height = 8, units = "in", res = 300)
  draw(ht)
  dev.off()

  # Save SVG
  svg_path <- file.path(dirs$figdir, paste0("rscu_heatmap_", segment, "_v2.svg"))
  svg(svg_path, width = 10, height = 8)
  draw(ht)
  dev.off()

  cat(sprintf("  Saved: rscu_heatmap_%s_v2.png + .svg\n", segment))
}

# ============================================================
# FIGURE 2: RSCU Dendrograms (ggtree, combined group + host)
# ============================================================
make_rscu_dendrogram <- function(rscu_mat, meta_df, segment) {
  dist_mat <- dist(rscu_mat, method = "euclidean")
  hc <- hclust(dist_mat, method = "average")
  tree <- as.treedata(as.phylo(hc))

  # Build tip data
  tip_labels <- hc$labels[hc$order]
  tip_data <- data.frame(
    label = tip_labels,
    Group = meta_df$GC[match(tip_labels, meta_df$accession)],
    Host = meta_df$host[match(tip_labels, meta_df$accession)],
    ShortLabel = sapply(tip_labels, function(sid) short_label(sid, meta_df))
  )

  # Panel A: colored by constellation group
  p_group <- ggtree(tree, layout = "rectangular", size = 0.8) %<+% tip_data +
    geom_tiplab(aes(label = ShortLabel, color = Group), size = 4.5, align = TRUE, linesize = 0.3,
                show.legend = FALSE) +
    geom_tippoint(aes(color = Group), size = 4.5, alpha = 0.8) +
    scale_color_manual(values = group_colors, name = "Constellation") +
    theme_tree2(plot.title = element_text(size = 14.0, face = "bold", hjust = 0)) +
    labs(title = "A. Colored by constellation group") +
    xlim(0, max(dist_mat) * 2.5)

  # Panel B: colored by host species
  p_host <- ggtree(tree, layout = "rectangular", size = 0.8) %<+% tip_data +
    geom_tiplab(aes(label = ShortLabel, color = Host), size = 4.5, align = TRUE, linesize = 0.3,
                show.legend = FALSE) +
    geom_tippoint(aes(color = Host), size = 4.5, alpha = 0.8) +
    scale_color_manual(values = host_colors, name = "Host species") +
    theme_tree2(plot.title = element_text(size = 14.0, face = "bold", hjust = 0)) +
    labs(title = "B. Colored by host species") +
    xlim(0, max(dist_mat) * 2.5)

  # Combine with patchwork
  combined <- p_group + p_host +
    plot_layout(ncol = 2, widths = c(1, 1)) +
    plot_annotation(
      title = paste0(segment, " — RSCU UPGMA Clustering"),
      subtitle = "Euclidean distance on 59 synonymous codons",
      theme = theme(plot.title = element_text(size = 16.0, face = "bold", hjust = 0),
                    plot.subtitle = element_text(size = 13.0, color = "gray40"))
    )

  # Save
  png_path <- file.path(dirs$figdir, paste0("rscu_dendrogram_combined_", segment, "_v2.png"))
  svg_path <- file.path(dirs$figdir, paste0("rscu_dendrogram_combined_", segment, "_v2.svg"))
  ggsave(png_path, combined, width = 14, height = 8, dpi = 300, units = "in")
  ggsave(svg_path, combined, width = 14, height = 8, units = "in")
  cat(sprintf("  Saved: rscu_dendrogram_combined_%s_v2.png + .svg\n", segment))
}

# ============================================================
# FIGURE 3: ENC-GC3s Neutrality Plots (ggplot2 + ggrepel)
# ============================================================
make_enc_gc3_plot <- function(gc_enc_df, meta_df, segment) {
  df <- merge(gc_enc_df, meta_df[, c("accession", "host", "GC", "length")], by = "accession")
  df$GC3s <- as.numeric(df$GC3s)
  df$ENC <- as.numeric(df$ENC)

  # Expected ENC curve: original Wright (1990) formula
  #   ENC = 2 + s + 29 / [s^2 + (1-s)^2],  s = GC3s as a fraction
  gc3s_range <- seq(0, 100, by = 0.5)
  enc_expected <- pmin(wright_expected_enc(gc3s_range / 100), 61)
  expected_df <- data.frame(GC3s = gc3s_range, ENC = enc_expected,
                            ENC_lower = enc_expected - 5, ENC_upper = enc_expected + 5)

  # Identify outliers (ENC < 35 or far from expected)
  df$expected_enc <- pmin(wright_expected_enc(df$GC3s / 100), 61)
  df$is_outlier <- df$ENC < 35 | abs(df$ENC - df$expected_enc) > 10
  df$label_text <- sapply(df$accession, function(sid) short_label(sid, meta_df))

  p <- ggplot(df, aes(x = GC3s, y = ENC)) +
    # Confidence band (±5 ENC)
    geom_ribbon(data = expected_df, aes(x = GC3s, ymin = ENC_lower, ymax = ENC_upper),
                fill = "gray85", alpha = 0.5, inherit.aes = FALSE) +
    # Expected curve
    geom_line(data = expected_df, aes(y = ENC), color = "gray40", linewidth = 0.8) +
    # Reference lines
    geom_hline(yintercept = 45, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    geom_hline(yintercept = 35, linetype = "dotted", color = "gray60", linewidth = 0.4) +
    annotate("text", x = 2, y = 46, label = "Nc=45 (moderate bias)", size = 4.5, hjust = 0, color = "gray50") +
    annotate("text", x = 2, y = 36, label = "Nc=35 (strong bias)", size = 4.5, hjust = 0, color = "gray50") +
    # Data points
    geom_point(aes(color = GC, size = length), alpha = 0.8) +
    scale_color_manual(values = group_colors, name = "Constellation") +
    scale_size_continuous(name = "Length (bp)", range = c(1.5, 4), breaks = c(500, 1000, 2000)) +
    # Outlier labels
    geom_text_repel(data = df[df$is_outlier, ], aes(label = label_text),
                    size = 2.5, max.overlaps = 15, segment.color = "gray60",
                    segment.size = 0.2, box.padding = 0.3) +
    labs(title = paste0(segment, " — ENC vs GC3s Neutrality Plot"),
         subtitle = "Expected curve (Wright 1990) with ±5 ENC confidence band",
         x = expression(GC[3*s]~"(%)"), y = "Effective Number of Codons (Nc)") +
    theme_cub() +
    theme(legend.position = "right") +
    coord_cartesian(ylim = c(20, 61), xlim = c(0, 100)) +
    guides(color = guide_legend(order = 1), size = guide_legend(order = 2))

  save_fig(p, paste0("enc_gc3_neutrality_", segment, "_v2"), width = 7, height = 5)
}

# ============================================================
# FIGURE 4: GC3 Boxplots by Constellation Group (ggplot2 + ggpubr)
# ============================================================
make_gc3_boxplot <- function(gc_enc_df, meta_df, segment) {
  df <- merge(gc_enc_df, meta_df[, c("accession", "host", "GC")], by = "accession")
  df$GC3 <- as.numeric(df$GC3)
  df$GC <- factor(df$GC, levels = names(group_colors)[names(group_colors) %in% unique(df$GC)])

  # Sample sizes
  n_labels <- paste0("n=", as.numeric(table(df$GC)))

  p <- ggplot(df, aes(x = GC, y = GC3, fill = GC)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.6) +
    geom_jitter(aes(color = host), width = 0.15, size = 3.8, alpha = 0.8) +
    scale_fill_manual(values = group_colors, guide = "none") +
    scale_color_manual(values = host_colors, name = "Host species") +
    stat_compare_means(method = "kruskal.test", label = "p.format",
                       label.x.npc = "left", label.y.npc = "top", size = 3) +
    labs(title = paste0(segment, " — GC3 by Constellation Group"),
         subtitle = "Kruskal-Wallis test for group differences",
         x = "Constellation Group", y = expression(GC[3]~"(%)")) +
    theme_cub() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 11.0),
          legend.position = "right") +
    annotate("text", x = seq_along(levels(df$GC)), y = min(df$GC3) - 2,
             label = n_labels, size = 2.5, color = "gray40")

  save_fig(p, paste0("gc3_by_group_", segment, "_v2"), width = 7, height = 5)
}

# ============================================================
# FIGURE 5: Tanglegrams (ggtree-based)
# ============================================================
load_pruned_tree <- function(segment, kept_ids) {
  tree_path <- TREE_FILES[[segment]]
  if (!file.exists(tree_path)) return(NULL)
  tree <- read.tree(tree_path)
  tree_tips <- tree$tip.label
  matched_tips <- c()
  for (sid in kept_ids) {
    if (sid %in% tree_tips) { matched_tips <- c(matched_tips, sid); next }
    short <- strsplit(sid, "_NODE_")[[1]][1]
    if (short %in% tree_tips) { matched_tips <- c(matched_tips, short); next }
    acc_base <- strsplit(sid, "\\.")[[1]][1]
    matches <- tree_tips[grepl(paste0("^", acc_base, "\\."), tree_tips)]
    if (length(matches) > 0) matched_tips <- c(matched_tips, matches[1])
  }
  tips_to_drop <- setdiff(tree_tips, matched_tips)
  if (length(tips_to_drop) > 0) tree <- drop.tip(tree, tips_to_drop)
  # Rename tips to match sequence IDs
  new_labels <- sapply(tree$tip.label, function(tl) {
    for (sid in kept_ids) {
      if (sid == tl) return(sid)
      short <- strsplit(sid, "_NODE_")[[1]][1]
      if (short == tl) return(sid)
      acc_base <- strsplit(sid, "\\.")[[1]][1]
      if (grepl(paste0("^", acc_base, "\\."), tl)) return(sid)
    }
    tl
  })
  tree$tip.label <- new_labels
  tree
}

make_tanglegram <- function(rscu_mat, meta_df, phylo_tree, segment) {
  if (is.null(phylo_tree)) { cat(sprintf("  Skipping tanglegram for %s (no tree)\n", segment)); return(NULL) }

  # RSCU tree
  dist_mat <- dist(rscu_mat, method = "euclidean")
  hc <- hclust(dist_mat, method = "average")
  rscu_tree <- as.phylo(hc)

  # Common tips
  common <- intersect(rscu_tree$tip.label, phylo_tree$tip.label)
  if (length(common) < 3) { cat(sprintf("  WARNING: Only %d common tips\n", length(common))); return(NULL) }

  # Prune both trees to common tips
  rscu_tree <- drop.tip(rscu_tree, setdiff(rscu_tree$tip.label, common))
  phylo_tree <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, common))

  # Ladderize in opposite directions for tanglegram
  phylo_tree <- ladderize(phylo_tree, right = TRUE)
  rscu_tree <- ladderize(rscu_tree, right = FALSE)

  # Tip data with group and short labels
  tip_data <- data.frame(
    label = common,
    Group = meta_df$GC[match(common, meta_df$accession)],
    ShortLabel = sapply(common, function(sid) short_label(sid, meta_df)),
    stringsAsFactors = FALSE
  )

  # Extract tip y-coordinates from ggtree data
  phylo_gd <- ggtree(phylo_tree)$data
  rscu_gd  <- ggtree(rscu_tree)$data
  phylo_tips <- phylo_gd[!is.na(phylo_gd$label) & phylo_gd$isTip, c("label", "y")]
  rscu_tips  <- rscu_gd[!is.na(rscu_gd$label) & rscu_gd$isTip, c("label", "y")]
  names(phylo_tips)[2] <- "y_phylo"
  names(rscu_tips)[2]  <- "y_rscu"

  # Connection data for linking lines
  conn_df <- merge(phylo_tips, rscu_tips, by = "label")
  conn_df <- merge(conn_df, tip_data[, c("label", "Group")], by = "label")

  # Panel 1: Phylogenetic tree (left, tips on right side)
  p_phylo <- ggtree(phylo_tree, layout = "rectangular", size = 0.6) %<+% tip_data +
    geom_tippoint(aes(color = Group), size = 4.0, alpha = 0.8) +
    geom_tiplab(aes(label = ShortLabel), size = 4.2, align = TRUE, linesize = 0.3,
                hjust = -0.1, family = "Liberation Sans") +
    scale_color_manual(values = group_colors, name = "Constellation",
                       guide = guide_legend(override.aes = list(size = 3))) +
    theme_tree2(plot.title = element_text(size = 11, face = "bold", hjust = 0,
                                          family = "Liberation Sans")) +
    labs(title = paste0(segment, " — Phylogenetic Tree")) +
    xlim_tree(max(node.depth.edgelength(phylo_tree)) * 2.5)

  # Panel 2: RSCU UPGMA tree (right, tips on left side, reversed x-axis)
  p_rscu <- ggtree(rscu_tree, layout = "rectangular", size = 0.6) %<+% tip_data +
    geom_tippoint(aes(color = Group), size = 4.0, alpha = 0.8) +
    geom_tiplab(aes(label = ShortLabel), size = 4.2, align = TRUE, linesize = 0.3,
                hjust = 1.1, family = "Liberation Sans") +
    scale_color_manual(values = group_colors, guide = "none") +
    theme_tree2(plot.title = element_text(size = 11, face = "bold", hjust = 1,
                                          family = "Liberation Sans")) +
    labs(title = paste0(segment, " — RSCU UPGMA")) +
    xlim_tree(max(node.depth.edgelength(rscu_tree)) * 2.5) +
    scale_x_reverse()

  # Panel 2 (middle): Connecting lines between corresponding tips
  y_range <- c(0.5, length(common) + 0.5)
  p_conn <- ggplot(conn_df) +
    geom_segment(aes(x = 0, y = y_phylo, xend = 1, yend = y_rscu, color = Group),
                 linewidth = 0.4, alpha = 0.6) +
    scale_color_manual(values = group_colors, guide = "none") +
    theme_void() +
    xlim(0, 1) +
    coord_cartesian(ylim = y_range)

  # Align y-axes on all three panels
  p_phylo <- p_phylo + coord_cartesian(ylim = y_range)
  p_rscu <- p_rscu + coord_cartesian(ylim = y_range)

  # Combine with patchwork: phylogeny | connecting lines | RSCU tree
  combined <- (p_phylo | p_conn | p_rscu) +
    plot_layout(widths = c(0.35, 0.15, 0.35), guides = "collect")

  # Save
  png_path <- file.path(dirs$figdir, paste0("tanglegram_ggtree_", segment, "_v2.png"))
  svg_path <- file.path(dirs$figdir, paste0("tanglegram_ggtree_", segment, "_v2.svg"))
  ggsave(png_path, combined, width = 14, height = 12,
         dpi = 300, units = "in")
  ggsave(svg_path, combined, width = 14, height = 12,
         units = "in")
  cat(sprintf("  Saved: tanglegram_ggtree_%s_v2.png + .svg\n", segment))

  # Compute metrics using cophenetic distances (trees may not be ultrametric)
  coph_phylo <- cophenetic.phylo(phylo_tree)
  coph_rscu  <- cophenetic.phylo(rscu_tree)
  common_ordered <- intersect(rownames(coph_phylo), rownames(coph_rscu))
  coph_phylo <- coph_phylo[common_ordered, common_ordered]
  coph_rscu  <- coph_rscu[common_ordered, common_ordered]
  # Baker's gamma = Spearman correlation of cophenetic distances
  baker_gamma <- cor(as.vector(coph_phylo), as.vector(coph_rscu), method = "spearman")
  # Entanglement (normalized): 0 = identical, 1 = maximally tangled
  entang <- (1 - baker_gamma) / 2

  data.frame(segment = segment, entanglement = entang, baker_gamma = baker_gamma,
             n_common_tips = length(common))
}

# ============================================================
# FIGURE 6: Cross-Segment Summary (patchwork, 2 panels)
# ============================================================
make_cross_segment_summary <- function() {
  ari_df <- read.csv(file.path(dirs$tabdir, "ari_results.csv"), stringsAsFactors = FALSE)
  tang_df <- read.csv(file.path(dirs$tabdir, "tanglegram_results.csv"), stringsAsFactors = FALSE)

  # Panel A: ARI grouped bar chart
  ari_long <- ari_df %>%
    select(segment, ARI_rscu_vs_group, ARI_rscu_vs_host, ARI_group_vs_host) %>%
    pivot_longer(cols = starts_with("ARI_"), names_to = "comparison", values_to = "ARI") %>%
    mutate(
      comparison = recode(comparison,
        "ARI_rscu_vs_group" = "RSCU vs\nConstellation",
        "ARI_rscu_vs_host" = "RSCU vs\nHost",
        "ARI_group_vs_host" = "Constellation\nvs Host"
      ),
      segment = factor(segment, levels = c("VP4", "VP6", "VP7", "NSP4"))
    )

  p_ari <- ggplot(ari_long, aes(x = segment, y = ARI, fill = comparison)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = sprintf("%.3f", ARI)),
              position = position_dodge(width = 0.8), vjust = -0.3, size = 2.8) +
    scale_fill_manual(values = c(
      "RSCU vs\nConstellation" = "#0279EE",
      "RSCU vs\nHost" = "#FF9400",
      "Constellation\nvs Host" = "#75A025"
    ), name = "Comparison") +
    labs(title = "A. Adjusted Rand Index",
         subtitle = "Concordance between RSCU clustering, constellation groups, and host species",
         x = "Segment", y = "Adjusted Rand Index (ARI)") +
    theme_cub() +
    theme(legend.position = "bottom", legend.text = element_text(size = 10.0)) +
    coord_cartesian(ylim = c(-0.3, 1.0)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.3) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "gray70", linewidth = 0.3) +
    annotate("text", x = 0.3, y = 0.95, label = "ARI=1 (perfect)", size = 4.5, color = "gray60", hjust = 0)

  # Panel B: Baker's gamma bar chart
  tang_df$segment <- factor(tang_df$segment, levels = c("VP4", "VP6", "VP7", "NSP4"))

  p_baker <- ggplot(tang_df, aes(x = segment, y = baker_gamma, fill = segment)) +
    geom_col(width = 0.6, alpha = 0.8) +
    geom_text(aes(label = sprintf("%.3f", baker_gamma)), vjust = -0.3, size = 5.0) +
    scale_fill_manual(values = c("VP4" = "#0279EE", "VP6" = "#FF9400",
                                  "VP7" = "#75A025", "NSP4" = "#FD9BED"),
                      guide = "none") +
    labs(title = "B. Baker's Gamma Correlation",
         subtitle = "Cophenetic distance correlation between RSCU tree and phylogenetic tree",
         x = "Segment", y = "Baker's Gamma (correlation)") +
    theme_cub() +
    coord_cartesian(ylim = c(-0.2, 0.8)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.3) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "gray70", linewidth = 0.3)

  # Combine
  combined <- (p_ari | p_baker) + plot_layout(widths = c(1.3, 1))
  save_fig(combined, "cross_segment_summary_v2", width = 12, height = 6)
}

# ============================================================
# MAIN: Generate all figures
# ============================================================
cat("=== Publication-Quality Figure Redesign ===\n\n")

# Load all data upfront
all_data <- list()
for (seg in SEGMENTS) {
  cat(sprintf("Loading %s...\n", seg))
  all_data[[seg]] <- load_segment_data(seg)
}

# Load GC/ENC data
gc_enc_all <- list()
for (seg in SEGMENTS) {
  gc_enc_all[[seg]] <- read.csv(file.path(dirs$tabdir, paste0("gc_enc_values_", seg, ".csv")),
                                 stringsAsFactors = FALSE)
}

# 1. RSCU Heatmaps
cat("\n--- RSCU Heatmaps ---\n")
for (seg in SEGMENTS) {
  rscu_mat <- compute_rscu_matrix(all_data[[seg]]$seqs)
  make_rscu_heatmap(rscu_mat, all_data[[seg]]$meta, seg)
}

# 2. RSCU Dendrograms (combined group + host)
cat("\n--- RSCU Dendrograms ---\n")
for (seg in SEGMENTS) {
  rscu_mat <- compute_rscu_matrix(all_data[[seg]]$seqs)
  make_rscu_dendrogram(rscu_mat, all_data[[seg]]$meta, seg)
}

# 3. ENC-GC3s Neutrality Plots
cat("\n--- ENC-GC3s Neutrality Plots ---\n")
for (seg in SEGMENTS) {
  make_enc_gc3_plot(gc_enc_all[[seg]], all_data[[seg]]$meta, seg)
}

# 4. GC3 Boxplots
cat("\n--- GC3 Boxplots ---\n")
for (seg in SEGMENTS) {
  make_gc3_boxplot(gc_enc_all[[seg]], all_data[[seg]]$meta, seg)
}

# 5. Tanglegrams
cat("\n--- Tanglegrams (ggtree) ---\n")
tang_results <- list()
for (seg in SEGMENTS) {
  rscu_mat <- compute_rscu_matrix(all_data[[seg]]$seqs)
  phylo_tree <- load_pruned_tree(seg, rownames(rscu_mat))
  tang <- make_tanglegram(rscu_mat, all_data[[seg]]$meta, phylo_tree, seg)
  if (!is.null(tang)) tang_results[[seg]] <- tang
}

# 6. Cross-Segment Summary
cat("\n--- Cross-Segment Summary ---\n")
# Save updated tanglegram metrics
if (length(tang_results) > 0) {
  tang_df <- do.call(rbind, tang_results)
  write.csv(tang_df, file.path(dirs$tabdir, "tanglegram_results_v2.csv"), row.names = FALSE)
  cat("  Saved tanglegram_results_v2.csv\n")
}
make_cross_segment_summary()

cat("\n=== All Figures Redesigned ===\n")
cat(sprintf("Output directory: %s\n", dirs$figdir))
