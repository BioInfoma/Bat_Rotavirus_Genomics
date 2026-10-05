# ============================================================
# common_bat.R — Shared functions for bat-only RVA CUB analysis
# 
#
#   All codon-usage indices are computed with standard
#   tools — seqinr (codon counts, RSCU, base composition) and
#   coRdon (ENC, CAI) — instead of in-house implementations.
#   The in-house functions count_cod(), rscu_from_counts(),
#   enc_wright(), cai_calc(), base_comp() and load_ref_rscu()
#   were REMOVED. What remains here are plotting themes, palettes,
#   IO helpers. The ENC-GC3s expected curve uses the original
#   Wright (1990) formula ENC = 2 + s + 29/[s^2 + (1-s)^2].
# ============================================================

# ---- Required packages ----
suppressPackageStartupMessages({
  library(seqinr)
  library(coRdon)
  library(Biostrings)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(vegan)
  library(ape)
  library(dendextend)
  library(mclust)
  library(dplyr)
  library(tidyr)
  library(ggtree)
  library(treeio)
  library(ggrepel)
  library(ggpubr)
  library(patchwork)
  library(scales)
  library(RColorBrewer)
  library(aplot)
})

# ---- Compatibility shim: ggplot2 4.0.x renamed is.waive -> is_waiver ----
# ggtree 3.14.0 has its own empty() that still calls is.waive()
if (exists("empty", envir=getNamespace("ggtree"))) {
  assignInNamespace("empty",
    function(df) is.null(df) || nrow(df) == 0 || ncol(df) == 0 || ggplot2::is_waiver(df),
    ns="ggtree")
}

# ---- Plotting theme  ----
theme_cub <- function() theme_bw(base_size=14, base_family="Liberation Sans") +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_line(size=0.3, color="gray90"),
    plot.title=element_text(size=15, face="bold", hjust=0),
    plot.subtitle=element_text(size=13, color="gray40"),
    axis.text=element_text(color="black", size=12),
    axis.title=element_text(size=13),
    legend.title=element_text(size=12),
    legend.text=element_text(size=11),
    legend.key.size=unit(0.4, "cm"),
    strip.text=element_text(size=13, face="bold"),
    strip.background=element_rect(fill="gray95", color=NA)
  )

# ---- Host species colors (15 species, colorblind-friendly) ----
# Based on Wong (2011) Nature Methods palette + extensions
host_colors <- c(
  "Eidolon helvum"          = "#E69F00",  # orange
  "Rhinolophus hipposideros"= "#56B4E9",  # sky blue
  "Rhinolophus blasii"      = "#009E73",  # bluish green
  "Rhinolophus euryale"     = "#F0E442",  # yellow
  "Rhinolophus simulator"   = "#0072B2",  # blue
  "Hipposideros gigas"      = "#D55E00",  # vermilion
  "Hipposideros pomona"     = "#CC79A7",  # reddish purple
  "Taphozous melanopogon"   = "#882255",  # wine
  "Taphozous mauritianus"   = "#332288",  # indigo
  "Rousettus aegyptiacus"   = "#117733",  # forest green
  "Rousettus leschenaulti"  = "#999933",  # olive
  "Scotophilus kuhlii"      = "#DDCC77",  # sand
  "Molossus molossus"       = "#44AA99",  # teal
  "Glossophaga soricina"    = "#AA4499",  # magenta
  "Carollia perspicillata"  = "#88CCEE",  # pale blue
  "Unknown bat"             = "#888888"   # gray
)

# ---- Helper: save figure as PNG + SVG ----
save_fig <- function(plot, filename, width=7, height=5, res=300) {
  png_path <- file.path(getOption("figdir", "c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/figures"),
                         paste0(filename, ".png"))
  svg_path <- file.path(getOption("figdir", "c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/figures"),
                         paste0(filename, ".svg"))
  if (inherits(plot, "ggplot") || inherits(plot, "patchwork")) {
    ggsave(png_path, plot, width=width, height=height, dpi=res, units="in")
    ggsave(svg_path, plot, width=width, height=height, units="in")
  } else {
    png(png_path, width=width, height=height, units="in", res=res)
    if (is.function(plot)) plot() else print(plot)
    dev.off()
    svg(svg_path, width=width, height=height)
    if (is.function(plot)) plot() else print(plot)
    dev.off()
  }
  cat(sprintf("  Saved: %s + .svg\n", filename))
}

# ---- Helper: add panel label to ggplot ----
add_panel_label <- function(label, x=-Inf, y=Inf) {
  annotate("text", x=x, y=y, label=label, fontface="bold",
           size=7, hjust=-0.3, vjust=1.5)
}

# ---- Constellation group colors (8 groups) ----
group_colors <- c(
  "Purple" = "#7B2D8E",    # SA11-like GC (G3P[2], Gabon bats)
  "Orange" = "#E67E22",    # G3 bat GC (China/Bulgaria/Zambia bats)
  "Green"  = "#27AE60",    # G30/G31 GC (Cameroon/Africa bats)
  "Brown"  = "#8B4513",    # G20 GC (Costa Rica bat)
  "Yellow" = "#F1C40F",    # G25 GC (Cameroon/Saudi bats)
  "Blue"   = "#2980B9",    # G36 GC (Kenya bat)
  "DarkGray"= "#555555",   # G33 GC (China bat)
  "Red"    = "#C0392B"     # B8-BT-OAU-2, this study
)

# ---- Standard genetic code ----
CT <- c(TTT="F",TTC="F",TTA="L",TTG="L",CTT="L",CTC="L",CTA="L",CTG="L",
        ATT="I",ATC="I",ATA="I",ATG="M",GTT="V",GTC="V",GTA="V",GTG="V",
        TCT="S",TCC="S",TCA="S",TCG="S",CCT="P",CCC="P",CCA="P",CCG="P",
        ACT="T",ACC="T",ACA="T",ACG="T",GCT="A",GCC="A",GCA="A",GCG="A",
        TAT="Y",TAC="Y",TAA="*",TAG="*",TGA="*",CAT="H",CAC="H",CAA="Q",CAG="Q",
        AAT="N",AAC="N",AAA="K",AAG="K",GAT="D",GAC="D",GAA="E",GAG="E",
        TGT="C",TGC="C",TGG="W",CGT="R",CGC="R",CGA="R",CGG="R",AGA="R",AGG="R",
        AGT="S",AGC="S",GGT="G",GGC="G",GGA="G",GGG="G")
STOPS <- c("TAA","TAG","TGA")
SYN <- names(CT)[!names(CT) %in% c("ATG","TGG",STOPS)]  # 59 sense codons

# ---- Synonymous codon families ----
FAM <- list(Phe=c("TTT","TTC"),Leu=c("TTA","TTG","CTT","CTC","CTA","CTG"),
            Ile=c("ATT","ATC","ATA"),Val=c("GTT","GTC","GTA","GTG"),
            Ser=c("TCT","TCC","TCA","TCG","AGT","AGC"),Pro=c("CCT","CCC","CCA","CCG"),
            Thr=c("ACT","ACC","ACA","ACG"),Ala=c("GCT","GCC","GCA","GCG"),
            Tyr=c("TAT","TAC"),His=c("CAT","CAC"),Gln=c("CAA","CAG"),
            Asn=c("AAT","AAC"),Lys=c("AAA","AAG"),Asp=c("GAT","GAC"),Glu=c("GAA","GAG"),
            Cys=c("TGT","TGC"),Arg=c("CGT","CGC","CGA","CGG","AGA","AGG"),Gly=c("GGT","GGC","GGA","GGG"))
DEG2 <- c("Phe","Tyr","Cys","His","Gln","Asn","Lys","Asp","Glu")
DEG3 <- "Ile"; DEG4 <- c("Val","Pro","Thr","Ala","Gly"); DEG6 <- c("Leu","Ser","Arg")

# Canonical codon order shared by seqinr::uco and coRdon (alphabetical DNA)
CODONS64 <- sort(names(CT))

# ============================================================
# Standard-tool wrappers (seqinr + coRdon)
# ============================================================

# ---- Documented preprocessing before any index calculation ----
# Uppercase; trim 3' end to the last complete codon (CDS may be
# partial); drop the terminal stop codon. Returns a character string.
std_trim_cds <- function(seq) {
  seq <- toupper(seq); n <- nchar(seq)
  if (n %% 3 != 0) { n <- n - (n %% 3); seq <- substr(seq, 1, n) }
  if (n < 3) return("")
  cd <- substring(seq, seq(1, n - 2, 3), seq(3, n, 3))
  if (length(cd) > 0 && cd[length(cd)] %in% STOPS) cd <- cd[-length(cd)]
  paste(cd, collapse = "")
}

# ---- Codon counts: seqinr::uco(index = "eff") ----
# Returns a named integer vector over the 64 codons (uppercase names,
# canonical alphabetical order matching coRdon).
std_counts <- function(seq) {
  sq <- std_trim_cds(seq)
  if (nchar(sq) < 3) return(setNames(rep(0L, 64), CODONS64))
  cnt <- uco(s2c(tolower(sq)), index = "eff")
  cnt <- setNames(as.integer(cnt), toupper(names(cnt)))
  cnt[CODONS64]
}

# ---- RSCU: seqinr::uco(index = "rscu") ----
# RSCU is undefined (NA) for codon families absent from a sequence;
# these are set to 0 (family not used), matching common practice.
std_rscu <- function(seq) {
  sq <- std_trim_cds(seq)
  r <- setNames(rep(NA_real_, 64), CODONS64)
  if (nchar(sq) >= 3) {
    v <- uco(s2c(tolower(sq)), index = "rscu")
    r[toupper(names(v))] <- as.numeric(v)
  }
  r[is.na(r)] <- 0
  r[CODONS64]
}

# ---- ENC: coRdon::ENC (Wright 1990) ----
# cnt_mat: matrix/data.frame of codon counts, rows = sequences,
# cols = 64 codons (any order; reordered internally to coRdon's).
std_enc <- function(cnt_mat) {
  # coRdon requires an INTEGER matrix; ENC must be computed on a multi-row
  # codonTable (coRdon drops dimensions internally on 1-row tables).
  cnt_mat <- as.matrix(cnt_mat)[, CODONS64, drop = FALSE]
  colnames(cnt_mat) <- CODONS64
  storage.mode(cnt_mat) <- "integer"
  # alt.init = FALSE: use the standard genetic-code degeneracy classes
  # (Wright 1990); coRdon's alt.init=TRUE default reclassifies alternative
  # initiation codons, which is meant for gene-expression prediction, not
  # for codon-bias description.
  as.numeric(coRdon::ENC(coRdon::codonTable(cnt_mat), alt.init = FALSE))
}

# ---- CAI: coRdon::CAI (Sharp & Li 1987) ----
# viral_cnt: named count vector (64 codons) for one viral sequence.
# host_cnt:  named count vector (64 codons) = host reference table.
# coRdon computes reference RSCU from the host table and applies the
# Sharp & Li algorithm (w = RSCU/RSCUmax per synonymous family;
# geometric mean over codons; single-codon families and stops excluded
# via stop.rm/alt.init defaults).
std_cai <- function(viral_cnt, host_cnt) {
  vm <- matrix(as.integer(viral_cnt[CODONS64]), nrow = 1,
               dimnames = list("viral", CODONS64))
  hm <- matrix(as.integer(host_cnt[CODONS64]), nrow = 1,
               dimnames = list("host", CODONS64))
  vct <- coRdon::codonTable(vm)
  hct <- coRdon::codonTable(hm)
  as.numeric(coRdon::CAI(vct, subsets = list(host = hct)))
}

# ---- Load host codon COUNT table from CSV ----
# CSV columns: codon, aa, count, rscu. Returns named count vector.
load_ref_counts <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  cnt <- setNames(as.numeric(df$count), toupper(df$codon))
  cnt[CODONS64]
}

# ---- Base composition ----
# Overall GC and per-position GC1/GC2/GC3 from seqinr on the trimmed
# CDS. Synonymous-site composition (A3s/T3s/G3s/C3s, GC3s) is derived
# by definition from the standard uco codon counts: third-base
# frequencies across the 59 synonymous codons (no published tool
# computes synonymous-site-only composition directly).
std_composition <- function(seq) {
  sq <- std_trim_cds(seq)
  if (nchar(sq) < 3) return(NULL)
  chars <- s2c(tolower(sq))
  cnt <- std_counts(seq)

  # Overall and per-position GC (seqinr)
  GCall <- GC(chars) * 100
  GC1v  <- GC1(chars) * 100
  GC2v  <- GC2(chars) * 100
  GC3v  <- GC3(chars) * 100
  GC12  <- (GC1v + GC2v) / 2

  # Overall base percentages
  ab <- table(factor(chars, levels = c("a", "t", "g", "c")))
  tot <- sum(ab)
  A_pct <- ab[["a"]] / tot * 100; T_pct <- ab[["t"]] / tot * 100
  G_pct <- ab[["g"]] / tot * 100; C_pct <- ab[["c"]] / tot * 100

  # Synonymous third-position bases from codon counts
  syn_cnt <- cnt[SYN]
  n3 <- sum(syn_cnt)
  third <- substr(SYN, 3, 3)
  A3s <- sum(syn_cnt[third == "A"]) / n3 * 100
  T3s <- sum(syn_cnt[third == "T"]) / n3 * 100
  G3s <- sum(syn_cnt[third == "G"]) / n3 * 100
  C3s <- sum(syn_cnt[third == "C"]) / n3 * 100

  data.frame(A_pct = A_pct, T_pct = T_pct, G_pct = G_pct, C_pct = C_pct,
             AT_pct = A_pct + T_pct, GC_pct = G_pct + C_pct,
             A3s = A3s, T3s = T3s, G3s = G3s, C3s = C3s,
             GC1 = GC1v, GC2 = GC2v, GC3 = GC3v, GC12 = GC12,
             GC3s = G3s + C3s, GCall = GCall)
}

# ---- Wright (1990) expected ENC curve (original publication) ----
# ENC = 2 + s + 29 / [s^2 + (1-s)^2], s = GC3s as a fraction.
wright_expected_enc <- function(s) {
  2 + s + 29 / (s^2 + (1 - s)^2)
}

# ---- Setup output directories ----
setup_dirs <- function(base_dir="c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility"){
  figdir <- file.path(base_dir, "figures")
  tabdir <- file.path(base_dir, "output_tables")
  fastadir <- file.path(base_dir, "data", "fastas")
  metadir <- file.path(base_dir, "data", "metadata")
  treedir <- file.path(base_dir, "data", "trees")
  dir.create(figdir, recursive=TRUE, showWarnings=FALSE)
  dir.create(tabdir, recursive=TRUE, showWarnings=FALSE)
  list(figdir=figdir, tabdir=tabdir, fastadir=fastadir, metadir=metadir, treedir=treedir)
}

# ---- Host species name normalization ----
normalize_host <- function(h){
  h <- trimws(h)
  h <- gsub("Eidolom", "Eidolon helvum", h, ignore.case=TRUE)
  h <- gsub("Straw-coloured fruit bat", "Eidolon helvum", h, ignore.case=TRUE)
  h <- gsub("Straw.colored fruit bat", "Eidolon helvum", h, ignore.case=TRUE)
  h <- gsub("Macronycteris gigas", "Hipposideros gigas", h, ignore.case=TRUE)
  h <- gsub("^Bat$", "Unknown bat", h)
  std <- c(
    "R. hipposideros" = "Rhinolophus hipposideros",
    "R. blasii" = "Rhinolophus blasii",
    "R. euryale" = "Rhinolophus euryale",
    "R. simulator" = "Rhinolophus simulator",
    "R. aegyptiacus" = "Rousettus aegyptiacus",
    "R. leschenaulti" = "Rousettus leschenaulti",
    "H. gigas" = "Hipposideros gigas",
    "H. pomona" = "Hipposideros pomona",
    "T. melanopogon" = "Taphozous melanopogon",
    "T. mauritianus" = "Taphozous mauritianus",
    "S. kuhlii" = "Scotophilus kuhlii",
    "M. molossus" = "Molossus molossus",
    "G. soricina" = "Glossophaga soricina",
    "C. perspicillata" = "Carollia perspicillata",
    "E. helvum" = "Eidolon helvum",
    "Eidolon helvum" = "Eidolon helvum"
  )
  for(pat in names(std)){
    if(tolower(h) == tolower(pat)) return(std[pat])
  }
  for(pat in names(std)){
    if(grepl(pat, h, ignore.case=TRUE)) return(std[pat])
  }
  h
}

# ---- Host-to-codon-table mapping ----
host_table_map <- function(){
  data.frame(
    host_species = c(
      "Rhinolophus hipposideros", "Rousettus aegyptiacus", "Molossus molossus",
      "Eidolon helvum",
      "Rhinolophus simulator", "Taphozous melanopogon", "Rousettus leschenaulti",
      "Glossophaga soricina", "Carollia perspicillata",
      "Rhinolophus blasii", "Rhinolophus euryale",
      "Hipposideros gigas", "Hipposideros pomona",
      "Taphozous mauritianus", "Scotophilus kuhlii"
    ),
    codon_table = c(
      "Rhinolophus_hipposideros_codon.csv", "Rousettus_aegyptiacus_codon.csv", "Molossus_molossus_codon.csv",
      "Eidolon_helvum_codon.csv",
      "Rhinolophus_simulator_codon.csv", "Saccopteryx_bilineata_codon.csv", "Rousettus_leschenaulti_codon.csv",
      "Glossophaga_soricina_codon.csv", "Glossophaga_soricina_codon.csv",
      "Rhinolophus_ferrumequinum_codon.csv", "Rhinolophus_ferrumequinum_codon.csv",
      "Hipposideros_armiger_codon.csv", "Hipposideros_armiger_codon.csv",
      "Saccopteryx_bilineata_codon.csv", "Pipistrellus_kuhlii_codon.csv"
    ),
    source = c(
      "Direct genome CDS (GCF_964194185.1)", "Direct genome CDS (GCF_014176215.1)", "Direct genome CDS (GCF_014108415.1)",
      "Figshare annotation (DOI:10.6084/m9.figshare.29117609.v2)",
      "Trinity transcriptome (SRR38065611-13)", "Proxy: S. bilineata (GCF_036850765.1, same family Emballonuridae)", "Trinity transcriptome (SRR835440-42,SRR6796666)",
      "Trinity transcriptome (SRR13259283-84,SRR13300776)", "Proxy: G. soricina (Trinity transcriptome, same family Phyllostomidae)",
      "Proxy: R. ferrumequinum (GCF_004115265.2, same genus)", "Proxy: R. ferrumequinum (GCF_004115265.2, same genus)",
      "Proxy: H. armiger (GCF_001890085.2, same genus)", "Proxy: H. armiger (GCF_001890085.2, same genus)",
      "Proxy: S. bilineata (GCF_036850765.1, same family)", "Proxy: P. kuhlii (GCF_014108245.1, same family)"
    ),
    stringsAsFactors = FALSE
  )
}
