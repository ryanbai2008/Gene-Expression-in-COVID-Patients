
# Load libraries
library(DESeq2)
library(limma)
library(readxl)
library(stringr)

# === INPUTS ===
# Load raw counts (genes x samples)

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

rownames(metadata) <- str_trim(metadata$Patient_ID)
# Keep only common samples
common_samples <- intersect(colnames(expr_mat), rownames(metadata))
expr_sub <- expr_mat[, common_samples]
colData <- metadata[common_samples, ]

stopifnot(all(colnames(expr_sub) == rownames(colData)))

# Filter genes: count >= 20 in at least 50% samplesm 4 is 50%, 2 is 80%
min_counts <- 20
min_samples <- floor(ncol(expr_sub) * 0.5)
keep_genes <- rowSums(expr_sub >= min_counts) >= min_samples
filtered_counts <- expr_sub[keep_genes, ]


colData$Status <- factor(colData$Status)
if(any(is.na(colData$Status))) {
  stop("NA values in Status column in metadata_filtered")
}

levels(colData$Status)

# Make syntactically valid names for the levels (no spaces, no special chars)
levels(colData$Status) <- make.names(levels(colData$Status))

print(levels(colData$Status))

colData$Status <- relevel(colData$Status, ref = "Control")


# === STEP 1: Create DESeq2 dataset ===
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = colData,
                              design = ~ Status)


# === STEP 2: Variance Stabilizing Transformation (rlog) ===
vst_dds <- vst(dds, blind = TRUE)
vst_matrix <- assay(vst_dds)

print(class(colData$Status))         # Should be "factor"
print(levels(colData$Status))        # Should show your groups (e.g. "Control", "COVID-19", etc)
print(table(colData$Status))         # Count per group, should be >0 for at least 2 groups


# === STEP 3: Limma DEA ===
group <- colData$Status
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)


fit <- lmFit(vst_matrix, design)

# Define contrast, adjust names to match your actual conditions
# For example, if your groups are "Infected" and "Control":
print(colnames(design))

contrast <- makeContrasts(COVID19_vs_Control = COVID.19 - Control, levels = design)

fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# === STEP 4: Extract and filter results ===
res <- topTable(fit2, adjust = "fdr", number = Inf, sort.by = "P")

deg_res <- res[res$adj.P.Val < 0.05 & abs(res$logFC) > 0.58, ]

# === STEP 5: Export results ===
write.csv(res, "limma_all_results.csv")
write.csv(deg_res, "limma_significant_DEGs.csv")

# Extract p-values and statistics
pvalues <- fit2$p.value
logFCs <- fit2$coefficients
tstats <- fit2$t

# Convert to data.frame for exporting or inspection
stats_df <- data.frame(
  Gene = rownames(pvalues),
  pvalue = pvalues[,1],      # Assuming one contrast
  logFC = logFCs[,1],
  t = tstats[,1]
)

write.csv(stats_df, "limma_stats.csv", row.names=FALSE)
