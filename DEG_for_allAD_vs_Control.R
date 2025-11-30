# ================================
# Load required libraries
# ================================
library(limma)
library(pheatmap)
library(ggplot2)

# ================================
# Load metadata and expression matrix
# ================================
metadata <- read.csv("GSE1297_metadata.csv")  # 3 columns: GSM_ID, Sample_Name, Diagnosis
expr_matrix <- read.csv("GSE1297_RMA_Normalized.csv", row.names = 1, check.names = FALSE)

# Remove .cel suffix from expression matrix column names
colnames(expr_matrix) <- sub("\\.cel$", "", colnames(expr_matrix))

# Align columns with metadata
expr_matrix <- expr_matrix[, metadata$GSM_ID]

# ================================
# Create AD vs Control group
# ================================
metadata$Group <- ifelse(metadata$Diagnosis == "Control", "Control", "AD")
group <- factor(metadata$Group)

# ================================
# Design matrix for limma
# ================================
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# ================================
# Fit linear model and TREAT
# ================================
fit <- lmFit(expr_matrix, design)
contrast.matrix <- makeContrasts(AD_vs_Control = AD - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
fit2_treat <- treat(fit2, lfc = 1)

# ================================
# Extract DEGs
# ================================
deg_treat <- topTreat(fit2_treat, number = Inf)
deg_treat$Gene <- rownames(deg_treat)

# Exploratory DEGs (p < 0.05 & |logFC| > 1)
deg_exploratory <- deg_treat[deg_treat$P.Value < 0.05 & abs(deg_treat$logFC) > 1, ]

# Save results
write.csv(deg_treat, "DEG_all_AD_vs_Control.csv")
write.csv(deg_exploratory, "DEG_exploratory_AD_vs_Control.csv")

# ================================
# Volcano plot
# ================================
deg_treat$negLogP <- -log10(deg_treat$P.Value)
deg_treat$Significant <- ifelse(deg_treat$P.Value < 0.05 & abs(deg_treat$logFC) > 1, "Yes", "No")

volcano <- ggplot(deg_treat, aes(x = logFC, y = negLogP, color = Significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("No"="grey", "Yes"="red")) +
  theme_minimal() +
  labs(title = "Volcano plot: AD vs Control", x = "log2 Fold Change", y = "-log10(P-value)")

# Save volcano plot
ggsave("Volcano_AD_vs_Control.png", plot = volcano, width = 8, height = 6, dpi = 300)

# ================================
# Heatmap of top 20 genes
# ================================
top20_genes <- head(deg_treat[order(-abs(deg_treat$logFC)), ], 20)
expr_top20 <- expr_matrix[rownames(top20_genes), ]

# Scale expression
expr_scaled <- t(scale(t(expr_top20)))

# Sample annotation
annotation_col <- data.frame(Group = metadata$Group)
rownames(annotation_col) <- metadata$GSM_ID

# Save heatmap as PNG
png("Heatmap_top20_genes_AD_vs_Control.png", width = 1000, height = 800)
pheatmap(expr_scaled, annotation_col = annotation_col, 
         main = "Top 20 genes by |logFC|", cluster_rows = TRUE, cluster_cols = TRUE)
dev.off()
