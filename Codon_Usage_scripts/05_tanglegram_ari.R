# ============================================================
# 05_tanglegram_ari.R
# Tanglegram comparing RSCU-based clustering vs phylogenetic tree,
# ARI between RSCU clusters and constellation groups / host species.
# Per-segment analysis.
# ============================================================

source("c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/scripts/common_bat.R")
dirs <- setup_dirs()

SEGMENTS <- c("VP4", "VP6", "VP7", "NSP4")
LEN_THRESH <- c(VP4=2300, VP6=1150, VP7=950, NSP4=500)
# De novo trees built in this study (MAFFT L-INS-i + IQ-TREE MFP, 1000 UFBoot)
TREE_DIR <- dirs$treedir

# ---- Tree file paths ----
TREE_FILES <- list(
  VP4  = file.path(TREE_DIR, "VP4_ml.treefile"),
  VP6  = file.path(TREE_DIR, "VP6_ml.treefile"),
  VP7  = file.path(TREE_DIR, "VP7_ml.treefile"),
  NSP4 = file.path(TREE_DIR, "NSP4_ml.treefile")
)

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

# ---- Load and prune phylogenetic tree to bat-only tips ----
load_pruned_tree <- function(segment, kept_ids) {
  tree_path <- TREE_FILES[[segment]]
  if (!file.exists(tree_path)) {
    cat(sprintf("  WARNING: Tree file not found: %s\n", tree_path))
    return (NULL)
  }
  
  tree <- read.tree(tree_path)
  
  # Get tip labels that match our sequences
  # Study sequences use full NODE names, references use accession.version
  tree_tips <- tree$tip.label
  
  # Find matching tips
  matched_tips <- c()
  for (sid in kept_ids) {
    # Direct match
    if (sid %in% tree_tips) {
      matched_tips <- c(matched_tips, sid)
      next
    }
    # Try short name for study sequences
    short <- strsplit(sid, "_NODE_")[[1]][1]
    if (short %in% tree_tips) {
      matched_tips <- c(matched_tips, short)
      next
    }
    # Try matching by accession (without version)
    acc_base <- strsplit(sid, "\\.")[[1]][1]
    matches <- tree_tips[grepl(paste0("^", acc_base, "\\."), tree_tips)]
    if (length(matches) > 0) {
      matched_tips <- c(matched_tips, matches[1])
    }
  }
  
  # Tips to drop (non-bat or non-matching)
  tips_to_drop <- setdiff(tree_tips, matched_tips)
  
  if (length(tips_to_drop) > 0) {
    tree_pruned <- drop.tip(tree, tips_to_drop)
  } else {
    tree_pruned <- tree
  }
  
  cat(sprintf("  Tree: %d tips (pruned from %d)\n", length(tree_pruned$tip.label), length(tree_tips)))
  
  # Rename tips to match our sequence IDs
  new_labels <- sapply(tree_pruned$tip.label, function(tl) {
    for (sid in kept_ids) {
      if (sid == tl) return(sid)
      short <- strsplit(sid, "_NODE_")[[1]][1]
      if (short == tl) return(sid)
      acc_base <- strsplit(sid, "\\.")[[1]][1]
      if (grepl(paste0("^", acc_base, "\\."), tl)) return(sid)
    }
    tl
  })
  tree_pruned$tip.label <- new_labels
  
  tree_pruned
}

# ---- Compute ARI between RSCU clusters and group labels ----
compute_ari <- function(rscu_mat, meta_df, segment) {
  # Euclidean distance + UPGMA
  dist_mat <- dist(rscu_mat, method = "euclidean")
  hc <- hclust(dist_mat, method = "average")
  
  # Number of constellation groups
  n_groups <- length(unique(meta_df$GC))
  
  # Cut tree into k = n_groups clusters
  rscu_clusters <- cutree(hc, k = n_groups)
  
  # Group labels (constellation)
  group_labels <- as.numeric(as.factor(meta_df$GC))
  
  # Host labels
  host_labels <- as.numeric(as.factor(meta_df$host))
  
  # ARI: RSCU clusters vs constellation groups
  ari_group <- adjustedRandIndex(rscu_clusters, group_labels)
  
  # ARI: RSCU clusters vs host species
  ari_host <- adjustedRandIndex(rscu_clusters, host_labels)
  
  # ARI: constellation groups vs host species (baseline)
  ari_group_host <- adjustedRandIndex(group_labels, host_labels)
  
  data.frame(
    segment = segment,
    n_sequences = nrow(rscu_mat),
    n_groups = n_groups,
    n_hosts = length(unique(meta_df$host)),
    ARI_rscu_vs_group = ari_group,
    ARI_rscu_vs_host = ari_host,
    ARI_group_vs_host = ari_group_host
  )
}

# ---- Tanglegram: RSCU dendrogram vs phylogenetic tree ----
plot_tanglegram <- function(rscu_mat, meta_df, phylo_tree, segment, outdir) {
  if (is.null(phylo_tree)) {
    cat(sprintf("  Skipping tanglegram for %s (no tree)\n", segment))
    return (NULL)
  }
  
  # RSCU dendrogram
  dist_mat <- dist(rscu_mat, method = "euclidean")
  hc <- hclust(dist_mat, method = "average")
  dend_rscu <- as.dendrogram(hc)
  
  # Phylogenetic dendrogram — trees are NOT ultrametric, so represent the
  # phylogeny as a dendrogram via hclust of its cophenetic distances
  dend_phylo <- as.dendrogram(hclust(as.dist(cophenetic.phylo(phylo_tree)), method = "average"))
  # Strip stray names attribute on labels (breaks dendextend label matching)
  labels(dend_phylo) <- as.character(labels(dend_phylo))
  
  # Match tip labels between the two dendrograms
  rscu_labels <- labels(dend_rscu)
  phylo_labels <- labels(dend_phylo)
  
  # Find common labels
  common <- intersect(rscu_labels, phylo_labels)
  if (length(common) < 3) {
    cat(sprintf("  WARNING: Only %d common tips for tanglegram\n", length(common)))
    return (NULL)
  }
  
  # Prune both dendrograms to KEEP only common labels
  # (dendextend::prune REMOVES the listed labels)
  dend_rscu_pruned <- prune(dend_rscu, setdiff(labels(dend_rscu), common))
  dend_phylo_pruned <- prune(dend_phylo, setdiff(labels(dend_phylo), common))
  
  # Color tips by constellation group
  meta_common <- meta_df[match(common, meta_df$accession), ]
  group_colors_common <- group_colors[meta_common$GC]
  
  labels_colors(dend_rscu_pruned) <- group_colors_common[order.dendrogram(dend_rscu_pruned)]
  labels_colors(dend_phylo_pruned) <- group_colors_common[order.dendrogram(dend_phylo_pruned)]
  
  # Compute entanglement (dendrogram-based, on pruned dendrograms)
  entang <- entanglement(dend_rscu_pruned, dend_phylo_pruned)

  # Baker's gamma via cophenetic distance correlation (robust to non-ultrametric trees)
  rscu_phy <- drop.tip(as.phylo(hc), setdiff(as.phylo(hc)$tip.label, common))
  phylo_p  <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, common))
  coph_r <- cophenetic.phylo(rscu_phy)[common, common]
  coph_p <- cophenetic.phylo(phylo_p)[common, common]
  baker_gamma <- cor(as.vector(coph_r), as.vector(coph_p), method = "spearman")
  
  # Plot tanglegram
  png(file.path(outdir, paste0("tanglegram_rscu_vs_phylo_", segment, ".png")),
      width = 16, height = 10, units = "in", res = 200)
  tanglegram(dend_phylo_pruned, dend_rscu_pruned,
             main_left = paste0(segment, " — Phylogenetic Tree"),
             main_right = paste0(segment, " — RSCU UPGMA"),
             margin_inner = 8,
             cex_main = 1.2,
             columns_width = c(5, 3, 5),
             lab.cex = 0.4)
  dev.off()
  
  cat(sprintf("  Saved: tanglegram_rscu_vs_phylo_%s.png\n", segment))
  cat(sprintf("  Entanglement: %.4f, Baker's gamma: %.4f\n", entang, baker_gamma))
  
  data.frame(segment = segment, entanglement = entang, baker_gamma = baker_gamma,
             n_common_tips = length(common))
}

# ---- RSCU dendrogram colored by host ----
plot_rscu_dendro_by_host <- function(rscu_mat, meta_df, segment, outdir) {
  dist_mat <- dist(rscu_mat, method = "euclidean")
  hc <- hclust(dist_mat, method = "average")
  dend <- as.dendrogram(hc)
  
  # Color tips by host species
  host_colors <- setNames(rainbow(length(unique(meta_df$host))),
                          unique(meta_df$host))
  tip_host_colors <- host_colors[meta_df$host[order.dendrogram(dend)]]
  labels_colors(dend) <- tip_host_colors
  
  png(file.path(outdir, paste0("rscu_dendrogram_by_host_", segment, ".png")),
      width = 12, height = 8, units = "in", res = 300)
  par(mar = c(7, 4, 4, 2))
  plot(dend,
       main = paste0(segment, " — RSCU UPGMA (colored by host)"),
       ylab = "Euclidean Distance (59 synonymous codons)",
       cex = 0.5, las = 2)
  legend("topright",
         legend = names(host_colors),
         fill = host_colors,
         title = "Host species",
         cex = 0.5, ncol = 2)
  dev.off()
  cat(sprintf("  Saved: rscu_dendrogram_by_host_%s.png\n", segment))
}

# ---- Main analysis loop ----
ari_results <- list()
tanglegram_results <- list()

for (seg in SEGMENTS) {
  cat(sprintf("\n=== %s ===\n", seg))
  data <- load_segment_data(seg)
  
  # Compute RSCU matrix
  rscu_mat <- t(sapply(names(data$seqs), function(sid) {
    std_rscu(data$seqs[[sid]])   # seqinr::uco(index = "rscu")
  }))
  rscu_mat <- rscu_mat[, SYN]
  rownames(rscu_mat) <- names(data$seqs)
  
  # ARI analysis
  ari <- compute_ari(rscu_mat, data$meta, seg)
  ari_results[[seg]] <- ari
  cat(sprintf("  ARI (RSCU vs group): %.4f\n", ari$ARI_rscu_vs_group))
  cat(sprintf("  ARI (RSCU vs host):  %.4f\n", ari$ARI_rscu_vs_host))
  cat(sprintf("  ARI (group vs host): %.4f\n", ari$ARI_group_vs_host))
  
  # Load and prune phylogenetic tree
  phylo_tree <- load_pruned_tree(seg, rownames(rscu_mat))
  
  # Tanglegram
  tang <- plot_tanglegram(rscu_mat, data$meta, phylo_tree, seg, dirs$figdir)
  if (!is.null(tang)) tanglegram_results[[seg]] <- tang
  
  # RSCU dendrogram by host
  plot_rscu_dendro_by_host(rscu_mat, data$meta, seg, dirs$figdir)
}

# ---- Save ARI results ----
ari_df <- do.call(rbind, ari_results)
write.csv(ari_df, file.path(dirs$tabdir, "ari_results.csv"), row.names = FALSE)
cat("\nSaved: ari_results.csv\n")
print(ari_df)

# ---- Save tanglegram results ----
if (length(tanglegram_results) > 0) {
  tang_df <- do.call(rbind, tanglegram_results)
  write.csv(tang_df, file.path(dirs$tabdir, "tanglegram_results.csv"), row.names = FALSE)
  cat("\nSaved: tanglegram_results.csv\n")
  print(tang_df)
}

cat("\n=== Tanglegram + ARI Analysis Complete ===\n")
