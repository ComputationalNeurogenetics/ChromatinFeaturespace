
genome_bins <- tileGenome(seqlengths = seqlengths(neuronal.merged),
                          tilewidth = 100,
                          cut.last.tile.in.chrom = TRUE)

# Create a FeatureMatrix from scATAC fragments
pileup_data_mat <- FeatureMatrix(
  fragments = Fragments(neuronal.merged),
  features = genome_bins,  # e.g., 100 bp bins
  cells = NULL
)

# Convert to a GRanges object with signal column
pileup_data <- GRanges(
  seqnames = seqnames(genome_bins),
  ranges = IRanges(start = start(genome_bins), end = end(genome_bins)),
  signal = rowSums(pileup_data_mat)  # Sum counts per bin
)

# Assuming you have three GRanges objects:
# `clade.peaks_01_bedpe` contains clade-level features
# `feat.5k.filt` contains 5kbp fixed bins
# `pileup_100bp` contains the ATAC-seq fragment pileups using 100bp bins
# `bulk_peaks` contains bulk ATAC-seq peaks as a reference (if available)

library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(GenomicAlignments)
library(Signac)
library(InteractionSet)

# --- 1) Peak Agreement (Jaccard Index) ---
compute_jaccard <- function(features1, features2) {
  intersection <- length(intersect(features1, features2))
  union_size <- length(union(features1, features2))
  jaccard_index <- intersection / union_size
  return(jaccard_index)
}

jaccard_clade_vs_bulk <- compute_jaccard(clade.peaks_01_bedpe, bulk_peaks)
jaccard_100bp_vs_bulk <- compute_jaccard(pileup_100bp, bulk_peaks)
jaccard_5kb_vs_bulk <- compute_jaccard(feat.5k.filt, bulk_peaks)

jaccard_results <- data.frame(
  Method = c("Clade-Based", "100bp Fixed Bins", "5kb Fixed Bins"),
  Jaccard_Index = c(jaccard_clade_vs_bulk, jaccard_100bp_vs_bulk, jaccard_5kb_vs_bulk)
)

p_jaccard <- ggplot(jaccard_results, aes(x = Method, y = Jaccard_Index, fill = Method)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  labs(title = "Jaccard Index Comparison", y = "Jaccard Index", x = "Method") +
  theme_minimal() +
  scale_fill_manual(values = c("blue", "green", "orange"))

print(p_jaccard)

# --- 2) Feature Stability and Consistency ---
compute_feature_retention <- function(features_list) {
  common_features <- Reduce(intersect, features_list)
  retention_score <- length(common_features) / length(unique(unlist(features_list)))
  return(retention_score)
}

feature_retention_clade <- compute_feature_retention(list(clade.peaks_01_bedpe, clade.peaks_02_bedpe)) # Across datasets
feature_retention_100bp <- compute_feature_retention(list(pileup_100bp, pileup_100bp_sample2))
feature_retention_5kb <- compute_feature_retention(list(feat.5k.filt, feat.5k.filt_sample2))

feature_retention_results <- data.frame(
  Method = c("Clade-Based", "100bp Fixed Bins", "5kb Fixed Bins"),
  Retention_Score = c(feature_retention_clade, feature_retention_100bp, feature_retention_5kb)
)

p_retention <- ggplot(feature_retention_results, aes(x = Method, y = Retention_Score, fill = Method)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  labs(title = "Feature Retention Score Comparison", y = "Retention Score", x = "Method") +
  theme_minimal() +
  scale_fill_manual(values = c("blue", "green", "orange"))

print(p_retention)

# --- 3) Sensitivity to Rare Cell-Specific Features ---
rare_clades <- names(clade_counts[clade_counts / sum(clade_counts) < 0.05])
rare_clade_features <- clade.peaks_01_bedpe[grep(paste(rare_clades, collapse="|"), clade.peaks_01_bedpe$peak_called_in)]
rare_100bp_features <- pileup_100bp[which(rowMeans(as.matrix(pileup_100bp$signal)) < quantile(pileup_100bp$signal, 0.05))]
rare_5kb_features <- feat.5k.filt[which(rowMeans(as.matrix(feat.5k.filt$signal)) < quantile(feat.5k.filt$signal, 0.05))]

sensitivity_results <- data.frame(
  Method = c("Clade-Based", "100bp Fixed Bins", "5kb Fixed Bins"),
  Rare_Features = c(length(rare_clade_features), length(rare_100bp_features), length(rare_5kb_features))
)

p_sensitivity <- ggplot(sensitivity_results, aes(x = Method, y = Rare_Features, fill = Method)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  labs(title = "Sensitivity to Rare Cell-Specific Features", y = "Rare Features Count", x = "Method") +
  theme_minimal() +
  scale_fill_manual(values = c("blue", "green", "orange"))

print(p_sensitivity)

# Print results
cat("Jaccard Index Comparison:\n")
print(jaccard_results)

cat("\nFeature Retention Comparison:\n")
print(feature_retention_results)

cat("\nSensitivity to Rare Features Comparison:\n")
print(sensitivity_results)
