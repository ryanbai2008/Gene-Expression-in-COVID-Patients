library(DESeq2)
library(readxl)
library(dplyr)
library(stringr)

# Load count data (genes x samples)
counts <- read.csv("data/GSE217948_combined_expression_matrix(1).csv", row.names=1)

expr_mat <- as.matrix(counts)

print(class(colnames(expr_mat)))
print("hello")

# Normalize sample names (column names)
normalize_name <- function(name) {
  # Remove 'd' followed by digits at the end of the string
  str_replace(name, "d\\d+$", "")
}

print(colnames(expr_mat)[100:200])
colnames(expr_mat) <- gsub("d\\d+", "", colnames(expr_mat))

print(colnames(expr_mat)[100:200])

# Load metadata
metadata <- read_excel("data/data.xlsx")

metadata <- as.data.frame(metadata)

rownames(metadata) <- metadata$Patient_ID

missing_samples <- setdiff(colnames(expr_mat), rownames(metadata))
print(missing_samples)

filtered_gene_data <- counts[, colnames(counts) %in% rownames(metadata)]
filtered_metadata <- metadata[rownames(metadata) %in% colnames(filtered_gene_data), ]

all(colnames(filtered_gene_data) == rownames(filtered_metadata))  # should be TRUE


# Keep only common samples
common_samples <- intersect(colnames(expr_mat), rownames(metadata))
expr_sub <- expr_mat[, common_samples]
metadata_sub <- metadata[common_samples, ]

# Filter genes: count >= 20 in at least 50% samplesm 4 is 50%, 2 is 80%
min_counts <- 20
min_samples <- floor(ncol(expr_sub) * 0.5)
keep_genes <- rowSums(expr_sub >= min_counts) >= min_samples
filtered_counts <- expr_sub[keep_genes, ]

# Remove samples with total counts outside 3 SD from mean
sample_sums <- colSums(filtered_counts)
mean_sum <- mean(sample_sums)
sd_sum <- sd(sample_sums)
non_outliers <- (sample_sums > (mean_sum - 2.5*sd_sum)) & (sample_sums < (mean_sum + 2.5*sd_sum))
filtered_counts <- filtered_counts[, non_outliers]
metadata_filtered <- metadata_sub[colnames(filtered_counts), ]

# Make sure Status is factor with no NA
metadata_filtered$Status <- factor(metadata_filtered$Status)
if(any(is.na(metadata_filtered$Status))) {
  stop("NA values in Status column in metadata_filtered")
}

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = metadata_filtered,
                              design = ~ Status)

dds$Status <- relevel(dds$Status, ref = "Control")

# Run DESeq2
dds <- DESeq(dds)

# Get normalized counts
norm_counts <- counts(dds, normalized=TRUE)

norm_counts <- log2(norm_counts+1)

# Save normalized counts
write.csv(norm_counts, "normalized_counts_DESeq2(4).csv")
