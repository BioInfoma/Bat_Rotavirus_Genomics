# =============================================================================
# Simsek et al. (2021) colour scheme | Nigerian OAU sequences (this study)
# Collapsed nodes: 93, 109  |  Output: 22 × 18 in @ 300 DPI
# =============================================================================

library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)

options(dplyr.summarise.inform = FALSE)

# -----------------------------------------------------------------------------
# FILE PATHS
# -----------------------------------------------------------------------------
treefile <- "RVA1_trimed.aligned.fasta.treefile"
metafile <- "VP1_metadata.csv"

# -----------------------------------------------------------------------------
# COLOURS
# -----------------------------------------------------------------------------
gc_colors <- c(
  "Purple"   = "#6B0F9E",
  "Orange"   = "#F39C12",
  "Green"    = "#27AE60",
  "Brown"    = "#8B4513",
  "Yellow"   = "#F1C40F",
  "Blue"     = "#2980B9",
  "DarkGray" = "#7F8C8D",
  "Red"      = "#C0392B",
  "Gray"     = "#BDC3C7"
)

# -----------------------------------------------------------------------------
# STEP 1 — READ
# -----------------------------------------------------------------------------
tree <- read.tree(treefile)
meta <- read.csv(metafile, stringsAsFactors = FALSE, na.strings = c("", "NA"))

# -----------------------------------------------------------------------------
# STEP 2 — ROOT
# -----------------------------------------------------------------------------
tree <- root(tree, outgroup = "AB009629.2", resolve.root = TRUE)
cat("Rooted on AB009629.2\n")

# -----------------------------------------------------------------------------
# STEP 3 — RENAME TIPS
# -----------------------------------------------------------------------------
name_map        <- setNames(meta$full_name, meta$accession)
tree$tip.label  <- ifelse(tree$tip.label %in% names(name_map),
                          name_map[tree$tip.label],
                          tree$tip.label)
meta$tree_label <- meta$full_name
meta            <- meta[meta$tree_label %in% tree$tip.label, ]
cat("Tips:", length(tree$tip.label), "| Metadata rows:", nrow(meta), "\n")

# -----------------------------------------------------------------------------
# STEP 4 — METADATA TABLE
# -----------------------------------------------------------------------------
clean_meta <- meta %>%
  mutate(
    label      = tree_label,
    is_study   = source == "Study",
    GC_clean   = as.character(GC),
    host_clean = case_when(
      host %in% c("Horse","Alpaca","Camelid","Guanaco") ~ "Horse/Camelid",
      host %in% c("Dog","Cat")                          ~ "Dog/Cat",
      host %in% c("Pig","Bovine","Mouse","Raccoon",
                  "Shrew","SugarGlider","Rat")          ~ "Other mammal",
      host == "Simian" ~ "Non-human primate",
      host == "Bird"   ~ "Avian",
      host == "Human"  ~ "Human",
      host == "Bat"    ~ "Bat",
      TRUE              ~ "Unknown"
    ),
    pt_size   = case_when(is_study ~ 5.5, source == "prev" ~ 5.0, TRUE ~ 4.0),
    pt_stroke = case_when(is_study ~ 2.5, TRUE ~ 1.8),
    pt_fill   = gc_colors[GC_clean],
    pt_fill   = ifelse(is.na(pt_fill), "#BDC3C7", pt_fill),
    pt_color  = case_when(is_study ~ "#C0392B",
                          source == "prev" ~ "black",
                          TRUE ~ "gray40"),
    pt_shape  = case_when(
      host_clean == "Bat"               ~ 21,
      host_clean == "Human"             ~ 22,
      host_clean == "Non-human primate" ~ 24,
      host_clean == "Horse/Camelid"     ~ 24,
      host_clean == "Dog/Cat"           ~ 23,
      host_clean == "Avian"             ~ 25,
      TRUE                               ~ 21
    ),
    txt_color = case_when(is_study ~ "#C0392B",
                          source == "prev" ~ "black",
                          TRUE ~ "black"),
    txt_face  = ifelse(is_study, "bold", "plain"),
    txt_size  = case_when(is_study ~ 11, source == "prev" ~ 11, TRUE ~ 11)
  ) %>%
  select(label, is_study, GC_clean, host_clean,
         pt_size, pt_stroke, pt_fill, pt_color, pt_shape,
         txt_color, txt_face, txt_size)

# -----------------------------------------------------------------------------
# BOOTSTRAP EXTRACTOR
# -----------------------------------------------------------------------------
extract_bootstrap <- function(node_label) {
  if (is.na(node_label) || node_label == "" || grepl("^Node", node_label))
    return(NA_real_)
  as.numeric(str_extract(node_label, "(?<=/)[0-9]+"))
}

# =============================================================================
# STEP 5 — BASE TREE
# =============================================================================
p <- ggtree(tree, ladderize = TRUE, size = 1.5, color = "gray10")

# =============================================================================
# STEP 6 — SCALE THEN COLLAPSE  (Human & other mammal RVAs)
# =============================================================================
p <- scaleClade(p, node = 93,  scale = 0.05)
p <- scaleClade(p, node = 109, scale = 0.05)
p <- collapse(p, node = 93,  mode = "mixed", fill = "black", color = "black")
p <- collapse(p, node = 109, mode = "mixed", fill = "black", color = "black")

collapsed_nodes <- c(93, 109)

# =============================================================================
# STEP 7 — ATTACH METADATA
# =============================================================================
p <- p %<+% clean_meta

p$data$bootstrap <- sapply(p$data$label, extract_bootstrap)
p$data$blab <- ifelse(
  !p$data$isTip & !is.na(p$data$bootstrap) & p$data$bootstrap >= 70,
  as.character(round(p$data$bootstrap)), NA
)

x_max    <- max(p$data$x, na.rm = TRUE)
xlim_max <- x_max * 2.8

# =============================================================================
# STEP 8 — COLLAPSED NODE LABELS
# =============================================================================
collapsed_label_df <- data.frame(
  node  = collapsed_nodes,
  label = rep("Human & other mammal RVAs", 2)
)
node_coords <- p$data %>%
  filter(node %in% collapsed_nodes) %>%
  select(node, x, y)
collapsed_label_df <- left_join(collapsed_label_df, node_coords, by = "node")

# =============================================================================
# STEP 9 — BUILD FINAL TREE PLOT

# =============================================================================
p2 <- p +
  
  # Collapsed node labels
  geom_text(data  = collapsed_label_df,
            aes(x = x_max * 1.01, y = y, label = label),
            hjust    = 0,
            size     = 10.0,          
            color    = "gray20",
            fontface = "italic") +
  
  # Bootstrap values
  geom_text(data = subset(p$data, !is.na(blab)),
            aes(x = x, y = y, label = blab),
            hjust = 1.4, vjust = -0.4,
            size  = 10,              
            color = "gray30") +
  
  # Tip points
  geom_point(data  = subset(p$data, isTip),
             aes(x = x, y = y,
                 fill   = pt_fill,
                 shape  = pt_shape,
                 size   = pt_size,
                 stroke = pt_stroke),
             color       = subset(p$data, isTip)$pt_color,
             show.legend = FALSE) +
  
  # Tip labels
  geom_text(data    = subset(p$data, isTip & !is.na(txt_size)),
            aes(x = x, y = y,
                label    = label,
                color    = txt_color,
                fontface = txt_face,
                size     = txt_size),
            hjust   = 0,
            nudge_x = x_max * 0.04) +
  
  geom_treescale(x = 0, y = -2, width = 0.1,
                 fontsize = 7.0,      
                 color = "gray40") +
  
  scale_shape_identity() +
  scale_fill_identity()  +
  scale_color_identity() +
  scale_size_identity()  +
  xlim(0, xlim_max) +
  theme_tree() +
  theme(
    plot.margin   = margin(15, 10, 15, 15),
    plot.title    = element_text(size = 44, face = "bold",
                                 hjust = 0, margin = margin(b = 8)),
    plot.subtitle = element_text(size = 28, color = "gray30",
                                 hjust = 0, margin = margin(b = 14))
  ) +
  labs(
    title    = "VP1 phylogeny of bat rotavirus A strains",
    subtitle = paste0("Nigerian study sequences (bold, red border) | n = ",
                      Ntip(tree),
                      " taxa | Bootstrap \u226570 shown | Rooted on avian PO\u201113 (AB009629.2)")
  )


# =============================================================================
# STEP 10 — LEGEND: SINGLE-COLUMN (FIX 1 — eliminates all overlap)
# =============================================================================

STP <- 1.45   # vertical step between entries
PT  <- 0.0    # point x position
TX  <- 0.9    # text x offset from point

# ── GC Constellation entries (9 items, single column) ──────────────────────
gc_title_y <- 21.0

gc_pts <- data.frame(
  x    = rep(PT, 9),
  y    = gc_title_y - (1:9) * STP,
  fill = unname(gc_colors),
  label = c(
    "Purple \u2013 R8 GC (G3P[2], Gabon bats)",
    "Orange \u2013 R3 GC ( G3 Bulgaria / Zambia bats)",
    "Green  \u2013 R15 GC (G30/G31 Cameroon)",
    "Brown  \u2013 R13 GC ( G20 Costa Rica)",
    "Yellow \u2013 R16 GC (G25, Cameroon / Saudi bats)",
    "Blue   \u2013 R22 GC (G36,Kenya bat )",
    "Darkgray \u2013 R19 GC ( G33, China)",
    "Red    \u2013 (this study)",
    "Gray   \u2013 Non-bat reference sequences"
  ),
  stringsAsFactors = FALSE
)

# ── Host species entries ───────────────────────────
host_title_y <- gc_title_y - 9 * STP - 2.2

host_pts <- data.frame(
  x     = rep(PT, 7),
  y     = host_title_y - (1:7) * STP,
  shape = c(21, 22, 24, 24, 23, 25, 21),
  label = c("Bat", "Human", "Non-human primate",
            "Horse / Camelid", "Dog / Cat", "Avian", "Other mammal"),
  stringsAsFactors = FALSE
)

y_min <- host_title_y - 8 * STP
y_max <- gc_title_y + 2.5

legend_plot <- ggplot() +
  
  # ── GC Title ──────────────────────────────────────────────────────────────
  annotate("text",
           x = 0, y = gc_title_y,
           label    = "Genotype constellations\n(Simsek et al. 2021)",
           hjust    = 0,
           fontface = "bold",
           size     = 11,           # was 10
           color    = "gray10") +
  
  # ── GC Points (filled squares) ────────────────────────────────────────────
  geom_point(data = gc_pts,
             aes(x = x, y = y),
             shape = 22, size = 6.5,
             fill  = gc_pts$fill,
             color = gc_pts$fill) +
  
  # ── GC Labels ─────────────────────────────────────────────────────────────
  geom_text(data = gc_pts,
            aes(x = x + TX, y = y, label = label),
            hjust = 0,
            size  = 10.0,          # was 8.0
            color = "gray10") +
  
  # ── Host Title ────────────────────────────────────────────────────────────
  annotate("text",
           x = 0, y = host_title_y,
           label    = "Host species",
           hjust    = 0,
           fontface = "bold",
           size     = 11,           # was 10
           color    = "gray10") +
  
  # ── Host Points (varied shapes) ───────────────────────────────────────────
  geom_point(data = host_pts,
             aes(x = x, y = y, shape = factor(shape)),
             size  = 6.5,
             fill  = "gray55",
             color = "gray30") +
  scale_shape_manual(
    values = c("21" = 21, "22" = 22, "23" = 23, "24" = 24, "25" = 25),
    guide  = "none"
  ) +
  
  # ── Host Labels ───────────────────────────────────────────────────────────
  geom_text(data = host_pts,
            aes(x = x + TX, y = y, label = label),
            hjust = 0,
            size  = 10.0,          # was 8.0
            color = "gray10") +
  
  xlim(-0.5, 22) +
  ylim(y_min, y_max) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5))


# =============================================================================
# STEP 11 — COMBINE AND SAVE
# =============================================================================
final <- (p2 | legend_plot) +
  plot_layout(widths = c(3.0, 1.8)) &
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave("VP1_tree_thesis.pdf",
       plot  = final,
       width = 22, height = 32,
       units = "in", dpi = 300,
       limitsize = FALSE)

# ── PNG ───────────────────────────────────────────────
ggsave("VP1_tree_thesis.png",
       plot  = final,
       width = 26, height = 32,
       units = "in", dpi = 300,
       limitsize = FALSE,
       bg    = "white")

cat("Done! Files saved: VP1_tree_thesis.pdf / VP1_tree_thesis.png\n")
