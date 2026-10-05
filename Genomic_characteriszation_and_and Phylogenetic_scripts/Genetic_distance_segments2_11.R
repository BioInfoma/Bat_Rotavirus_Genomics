# Genetic distance calculation for Rotavirus Segments 2-11

library(ape)
library(dplyr)

setwd("C:/Users/USER/Documents/Masters_Project/Corected_msc/Project_Interpretation/Genetic distance")

find_closest_relative <- function(nigerian_seq, distance_matrix, exclude_nigerian = TRUE) {
  
  distances_to_all <- distance_matrix[nigerian_seq, ]
  
  distances_to_all <- distances_to_all[distances_to_all > 0]
  
  if(exclude_nigerian) {
    non_nigerian <- names(distances_to_all)[!grepl("^B[0-9]+-", names(distances_to_all))]
    distances_to_all <- distances_to_all[non_nigerian]
  }
  
  closest_name <- names(distances_to_all)[which.min(distances_to_all)]
  closest_distance <- min(distances_to_all)
  
  return(list(name = closest_name, distance = closest_distance))
}

process_segment <- function(segment_num, protein_name, tree_file, root_accession) {
  
  cat("\n", paste(rep("=", 70), collapse=""), "\n")
  cat(sprintf("SEGMENT %d (%s) - Rooting on: %s\n", segment_num, protein_name, root_accession))
  cat(paste(rep("=", 70), collapse=""), "\n\n")
  
  if(!file.exists(tree_file)) {
    cat("ERROR: Tree file not found:", tree_file, "\n")
    return(NULL)
  }
  tree <- read.tree(tree_file)
  cat("Tree loaded:", tree_file, "\n")
  cat("Number of tips:", length(tree$tip.label), "\n")
  
  root_tip <- tree$tip.label[tree$tip.label == root_accession]
  if(length(root_tip) == 0) {
    root_base <- sub("\\.[0-9]+$", "", root_accession)
    root_tip <- tree$tip.label[grepl(root_base, tree$tip.label)]
    if(length(root_tip) > 0) {
      cat("NOTE: Exact root accession not found. Using partial match:", root_tip[1], "\n")
      root_accession <- root_tip[1]
    } else {
      cat("ERROR: Root accession", root_accession, "not found in tree!\n")
      cat("Available tips:\n")
      print(head(tree$tip.label, 20))
      return(NULL)
    }
  }
  
  tree <- root(tree, outgroup = root_accession, resolve.root = TRUE)
  cat("Tree successfully rooted on:", root_accession, "\n")
  
  cophenetic_distances <- cophenetic.phylo(tree)
  cat("Cophenetic distance matrix calculated:", nrow(cophenetic_distances), "x", ncol(cophenetic_distances), "\n")
  
  nigerian_seqs <- grep("^B[0-9]+-", tree$tip.label, value = TRUE)
  cat("Nigerian sequences found:", length(nigerian_seqs), "\n")
  
  if(length(nigerian_seqs) == 0) {
    cat("WARNING: No Nigerian sequences found in this tree!\n")
    return(NULL)
  }
  
  for(ns in nigerian_seqs) {
    cat("  -", ns, "\n")
  }
  
  cat("\n--- Calculating genetic distances ---\n\n")
  
  results <- data.frame()
  
  for(seq_id in nigerian_seqs) {
    
    closest <- find_closest_relative(seq_id, cophenetic_distances)
    
    genetic_distance <- closest$distance
    calculated_identity <- (1 - genetic_distance) * 100
    
    len_match <- regmatches(seq_id, regexpr("length_[0-9]+", seq_id))
    seq_length <- ifelse(length(len_match) > 0, 
                         as.numeric(sub("length_", "", len_match)), 
                         NA)
    
    estimated_substitutions <- genetic_distance * seq_length
    
    results <- rbind(results, data.frame(
      Segment = segment_num,
      Protein = protein_name,
      Nigerian_Sequence = seq_id,
      Closest_Relative = closest$name,
      Genetic_Distance = round(genetic_distance, 6),
      Calculated_Identity = round(calculated_identity, 2),
      Sequence_Length = seq_length,
      Estimated_Substitutions = round(estimated_substitutions, 1),
      stringsAsFactors = FALSE
    ))
    
    seq_short <- gsub("_NODE.*", "", seq_id)
    cat(sprintf("  %s:\n", seq_short))
    cat(sprintf("    Closest relative: %s\n", closest$name))
    cat(sprintf("    Genetic distance: %.6f substitutions/site\n", genetic_distance))
    cat(sprintf("    Calculated identity: %.2f%%\n", calculated_identity))
    cat(sprintf("    Sequence length: %s nt\n", ifelse(is.na(seq_length), "unknown", seq_length)))
    cat(sprintf("    Estimated substitutions: %.1f\n\n", estimated_substitutions))
  }
  
  cat("\n--- Segment", segment_num, "Summary ---\n")
  cat("Number of Nigerian sequences:", nrow(results), "\n")
  cat("Genetic distance range:", 
      round(min(results$Genetic_Distance), 6), "-", 
      round(max(results$Genetic_Distance), 6), "substitutions/site\n")
  cat("Identity range:", 
      round(min(results$Calculated_Identity), 2), "% -", 
      round(max(results$Calculated_Identity), 2), "%\n")
  cat("Average estimated substitutions:", 
      round(mean(results$Estimated_Substitutions, na.rm = TRUE), 1), "\n")
  
  csv_file <- sprintf("Segment%d_%s_genetic_distances.csv", segment_num, protein_name)
  write.csv(results, csv_file, row.names = FALSE)
  cat("Results exported to:", csv_file, "\n")
  
  return(results)
}


segments <- list(
  list(num = 2,  protein = "VP2",  
       tree = "RVA2_trimmed.aligned.fasta.treefile",
       root = "AB009630.2"),
  
  list(num = 3,  protein = "VP3",  
       tree = "VP3_plus_mine_all_aligned.fasta.treefile",
       root = "AB009631.2"),
  
  list(num = 4,  protein = "VP4",  
       tree = "RVA4_trimmed.aligned.fasta.treefile",
       root = "AB009632.2"),
  
  list(num = 5,  protein = "NSP1", 
       tree = "NSP1_trimmed_aa.aligned.fasta.treefile",
       root = "AB009633.2"),
  
  list(num = 6,  protein = "VP6",  
       tree = "RVA6_trimmed.aligned.fasta.treefile",
       root = "D16329.2"),
  
  list(num = 7,  protein = "NSP3", 
       tree = "NSP3.reseq_plus_mine_aligned.fasta.treefile",
       root = "AB009626.2"),
  
  list(num = 8,  protein = "NSP2", 
       tree = "NSP2_trimmed_aa.aligned.fasta.treefile",
       root = "AB009625.2"),
  
  list(num = 9,  protein = "VP7",  
       tree = "RVA7_trimed.aligned2.fasta.treefile",
       root = "D82979.2"),
  
  list(num = 10, protein = "NSP4", 
       tree = "NSP4.aligned_plus_top_hit_trimmed.fasta.treefile",
       root = "AB009627.1"),
  
  list(num = 11, protein = "NSP5", 
       tree = "NSP5_trimmed_aa.aligned.fasta.treefile",
       root = "AB009628.2")
)


cat("\n")
cat(paste(rep("#", 70), collapse=""), "\n")
cat("#  ROTAVIRUS GENETIC DISTANCE ANALYSIS - SEGMENTS 2-11\n")
cat("#  Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("#", 70), collapse=""), "\n")

all_results <- data.frame()

for(seg in segments) {
  result <- process_segment(seg$num, seg$protein, seg$tree, seg$root)
  if(!is.null(result)) {
    all_results <- rbind(all_results, result)
  }
}

cat("\n\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("COMBINED SUMMARY - ALL SEGMENTS\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

print(all_results)

write.csv(all_results, "All_Segments_Summary.csv", row.names = FALSE)
cat("\n\nCombined results exported to: All_Segments_Summary.csv\n")

cat("\n--- Per-Segment Summary Statistics ---\n\n")
segment_summary <- all_results %>%
  group_by(Segment, Protein) %>%
  summarise(
    N_Sequences = n(),
    Min_Distance = round(min(Genetic_Distance), 6),
    Max_Distance = round(max(Genetic_Distance), 6),
    Mean_Distance = round(mean(Genetic_Distance), 6),
    Min_Identity = round(min(Calculated_Identity), 2),
    Max_Identity = round(max(Calculated_Identity), 2),
    Mean_Identity = round(mean(Calculated_Identity), 2),
    Mean_Substitutions = round(mean(Estimated_Substitutions, na.rm=TRUE), 1),
    .groups = "drop"
  )
print(as.data.frame(segment_summary))

write.csv(segment_summary, "Segment_Summary_Statistics.csv", row.names = FALSE)
cat("\nSummary statistics exported to: Segment_Summary_Statistics.csv\n")

cat("\n--- Overall Statistics (All Segments) ---\n")
cat("Total Nigerian sequences analyzed:", nrow(all_results), "\n")
cat("Overall distance range:", 
    round(min(all_results$Genetic_Distance), 6), "-", 
    round(max(all_results$Genetic_Distance), 6), "\n")
cat("Overall identity range:", 
    round(min(all_results$Calculated_Identity), 2), "% -", 
    round(max(all_results$Calculated_Identity), 2), "%\n")

cat("\n\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("FOR YOUR MANUSCRIPT:\n")
cat(paste(rep("=", 70), collapse=""), "\n")

for(i in 1:nrow(all_results)) {
  seq_short <- gsub("_NODE.*", "", all_results$Nigerian_Sequence[i])
  cat(sprintf("\nSegment %d (%s) - %s:\n", 
              all_results$Segment[i], all_results$Protein[i], seq_short))
  cat(sprintf("  - Closest relative: %s\n", all_results$Closest_Relative[i]))
  cat(sprintf("  - Genetic distance: %.4f substitutions/site\n", all_results$Genetic_Distance[i]))
  cat(sprintf("  - Nucleotide identity: %.2f%%\n", all_results$Calculated_Identity[i]))
  cat(sprintf("  - Estimated substitutions: %.0f (over %s nt)\n", 
              all_results$Estimated_Substitutions[i],
              ifelse(is.na(all_results$Sequence_Length[i]), "?", all_results$Sequence_Length[i])))
}

cat("\n\nAnalysis complete!\n")
cat("Output files:\n")
cat("  - Individual segment CSVs: Segment[X]_[Protein]_genetic_distances.csv\n")
cat("  - Combined summary: All_Segments_Summary.csv\n")
cat("  - Summary statistics: Segment_Summary_Statistics.csv\n")
