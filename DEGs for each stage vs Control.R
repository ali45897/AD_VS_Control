# ================================
# Load required libraries
# ================================
if (!requireNamespace("limma", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("limma")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(limma)
library(pheatmap)
library(ggplot2)
library(ggrepel)

# ================================
# Load data
# ================================
metadata <- read.csv("GSE1297_metadata.csv")  # Columns: GSM_ID, Sample_Name, Diagnosis
expr_matrix <- read.csv("GSE1297_RMA_Normalized.csv", row.names = 1, check.names = FALSE)

# Clean column names
colnames(expr_matrix) <- sub("\\.cel$", "", colnames(expr_matrix))
expr_matrix <- expr_matrix[, metadata$GSM_ID]

# Replace spaces in Diagnosis
metadata$Diagnosis <- gsub(" ", "_", metadata$Diagnosis)
group <- factor(metadata$Diagnosis)
levels(group)  # Should show: Control, Incipient_AD, Moderate_AD, Severe_AD

# ================================
# Design matrix and contrasts
# ================================
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

contrast.matrix <- makeContrasts(
  Severe_vs_Control = Severe_AD - Control,
  Moderate_vs_Control = Moderate_AD - Control,
  Incipient_vs_Control = Incipient_AD - Control,
  levels = design
)

# ================================
# Fit model
# ================================
fit <- lmFit(expr_matrix, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
fit2_treat <- treat(fit2, lfc = 1)

# ================================
# Stage-specific DEG analysis
# ================================
stages <- c("Severe_vs_Control", "Moderate_vs_Control", "Incipient_vs_Control")

for (stage in stages) {
  
  # Extract DEGs
  deg <- topTreat(fit2_treat, coef = stage, number = Inf)
  deg$Gene <- rownames(deg)
  
  # Save all DEGs
  write.csv(deg, paste0("DEG_all_", stage, ".csv"))
  
  # Exploratory DEGs: p < 0.05 & |logFC| > 1
  deg_expl <- deg[deg$P.Value < 0.05 & abs(deg$logFC) > 1, ]
  write.csv(deg_expl, paste0("DEG_exploratory_", stage, ".csv"))
  
  # ================================
  # Volcano plot with top 10 genes labeled
  # ================================
  deg$negLogP <- -log10(deg$P.Value)
  deg$Significant <- ifelse(deg$P.Value < 0.05 & abs(deg$logFC) > 1, "Yes", "No")
  
  # Top 10 genes by |logFC|
  top10 <- head(deg[order(-abs(deg$logFC)), ], 10)
  
  volcano <- ggplot(deg, aes(x = logFC, y = negLogP, color = Significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("No"="grey", "Yes"="red")) +
    geom_text_repel(data = top10, aes(label = Gene),
                    size = 3, box.padding = 0.5, max.overlaps = 10) +
    theme_minimal() +
    labs(title = paste0("Volcano plot: ", stage),
         x = "log2 Fold Change", y = "-log10(P-value)")
  
  ggsave(paste0("Volcano_", stage, ".png"), plot = volcano, width = 8, height = 6, dpi = 300)
  
  # ================================
  # Heatmap of top 20 genes by |logFC|
  # ================================
  top20_genes <- head(deg[order(-abs(deg$logFC)), ], 20)
  expr_top20 <- expr_matrix[rownames(top20_genes), ]
  expr_scaled <- t(scale(t(expr_top20)))  # z-score scaling
  
  annotation_col <- data.frame(Group = metadata$Diagnosis)
  rownames(annotation_col) <- metadata$GSM_ID
  
  png(paste0("Heatmap_top20_", stage, ".png"), width = 1000, height = 800)
  pheatmap(expr_scaled, annotation_col = annotation_col,
           main = paste0("Top 20 genes by |logFC|: ", stage),
           cluster_rows = TRUE, cluster_cols = TRUE)
  dev.off()
}

