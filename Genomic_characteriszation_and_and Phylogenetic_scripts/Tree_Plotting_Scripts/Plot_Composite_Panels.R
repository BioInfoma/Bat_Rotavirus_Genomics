library(magick)
library(cowplot)
library(ggplot2)

print("Starting panel creation...")

# Paths to the generated PDFs (or PNGs)
trees <- c(
  VP1 = "VP1_tree/VP1_tree_thesis.png",
  VP2 = "VP2/VP2_tree_thesis.png",
  VP3 = "VP3_plus_my_seq.trimmed.aligned.fasta/VP3_tree_thesis.png",
  VP4 = "RVA4_trimmed.aligned.fasta/VP4_tree_thesis.png",
  VP6 = "RVA6_trimmed.aligned.fasta/VP6_tree_thesis.png",
  VP7 = "RVA7_trimed.aligned2.fasta/VP7_tree_thesis.png",
  NSP1 = "NSP1_trimmed_aa.aligned.fasta/NSP1_tree_thesis.png",
  NSP2 = "NSP2_trimmed_aa.aligned.fasta/NSP2_tree_thesis.png",
  NSP3 = "NSP3.reseq_plus_mine_aligned.fasta/NSP3_tree_thesis.png",
  NSP4 = "NSP4_trimmed_aa.aligned.fasta/NSP4_tree_thesis.png",
  NSP5 = "NSP5_trimmed_aa.aligned.fasta/NSP5_tree_thesis.png"
)

# Function to read and prepare image
prep_img <- function(path) {
  img <- magick::image_read(path)
  # No manual label added here to avoid overlapping with the tree's built-in title
  p <- ggdraw() + draw_image(img)
  return(p)
}

print("Loading images...")
plots <- list()
for (name in names(trees)) {
  print(paste("Processing", name))
  plots[[name]] <- prep_img(trees[[name]])
}

print("Creating Tree 1 (VP1, VP2, VP3, NSP1)...")
panel1 <- plot_grid(plots[["VP1"]], plots[["VP2"]], plots[["VP3"]], plots[["NSP1"]],
                    ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), label_size = 80, scale = 0.95, label_x = 0.02, label_y = 0.98)
ggsave("Tree1_VP1_VP2_VP3_NSP1_final.png", panel1, width = 30, height = 30, limitsize = FALSE, dpi = 300)

print("Creating Tree 2 (VP4, VP6, VP7)...")
panel2 <- plot_grid(plots[["VP4"]], plots[["VP6"]], plots[["VP7"]], NULL,
                    ncol = 2, nrow = 2, labels = c("A", "B", "C", ""), label_size = 80, scale = 0.95, label_x = 0.02, label_y = 0.98)
ggsave("Tree2_VP4_VP6_VP7_final.png", panel2, width = 30, height = 30, limitsize = FALSE, dpi = 300)

print("Creating Tree 3 (NSP2, NSP3, NSP4, NSP5)...")
panel3 <- plot_grid(plots[["NSP2"]], plots[["NSP3"]], plots[["NSP4"]], plots[["NSP5"]],
                    ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), label_size = 80, scale = 0.95, label_x = 0.02, label_y = 0.98)
ggsave("Tree3_NSP2_NSP3_NSP4_NSP5_final.png", panel3, width = 30, height = 30, limitsize = FALSE, dpi = 300)

print("Done creating all panels!")
