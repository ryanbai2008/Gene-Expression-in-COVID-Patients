stats_df <- read.csv("limma_stats.csv")

head(stats_df)
summary(stats_df$pvalue)

# Plot distribution of p-values
hist(stats_df$pvalue, breaks=50, main="Histogram of raw p-values", xlab="p-value")

# Optionally check adjusted p-values (FDR) if you saved them or recalculate here
stats_df$adj_pvalue <- p.adjust(stats_df$pvalue, method = "fdr")

summary(stats_df$adj_pvalue)
