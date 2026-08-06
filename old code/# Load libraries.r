# Load libraries
library(DESeq2)
library(readxl)

# Load count data (assume CSV with genes as rows, samples as columns)
counts <- read.csv("GSE217948_combined_expression_matrix(2).csv", row.names=1)

# Convert to matrix
counts_mat <- as.matrix(counts)

# Filter genes: keep genes with count >= 20 in at least 50% samples
min_counts <- 20
min_samples <- floor(ncol(counts_mat) * 0.5)

keep_genes <- rowSums(counts_mat >= min_counts) >= min_samples
filtered_counts <- counts_mat[keep_genes, ]

# Remove samples with extreme total counts (outliers beyond 3 sd)
sample_sums <- colSums(filtered_counts)
mean_sum <- mean(sample_sums)
sd_sum <- sd(sample_sums)
non_outliers <- (sample_sums > (mean_sum - 3*sd_sum)) & (sample_sums < (mean_sum + 3*sd_sum))
filtered_counts <- filtered_counts[, non_outliers]

# Load metadata (assumes you have metadata file with matching sample IDs)
metadata <- read.csv("metadata.csv", row.names=1)

# Filter metadata to samples remaining after filtering counts
metadata_filtered <- metadata[colnames(filtered_counts), ]

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = metadata_filtered,
                              design = ~ condition)  # Adjust 'condition' as needed

# Run DESeq2 normalization and analysis
dds <- DESeq(dds)

# Get normalized counts
norm_counts <- counts(dds, normalized=TRUE)

expr_ids <- unique(gsub("d\\d+", "", colnames(counts)))

metadata_ids <- metadata$Patient_ID

only_in_expr <- setdiff(expr_ids, metadata_ids)
only_in_metadata <- setdiff(metadata_ids, expr_ids)

length(only_in_expr)
length(only_in_metadata)

head(only_in_expr)
head(only_in_metadata)

# Write normalized counts to CSV
# write.csv(norm_counts, "normalized_counts_DESeq2.csv")
