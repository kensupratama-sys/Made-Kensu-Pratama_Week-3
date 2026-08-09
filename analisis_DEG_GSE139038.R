## ============================================================
## Analisis Differential Gene Expression (DEG) dan Functional
## Enrichment pada Dataset GSE139038 (Breast Cancer vs Paired Normal)
## ============================================================

## ---- 1. Load Library ----
library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(gplots)
library(RColorBrewer)
library(gprofiler2)

## ---- 2. Unduh Data dari GEO ----
gse <- getGEO(
  "GSE139038",
  GSEMatrix = TRUE,
  getGPL = FALSE
)

expr_full <- exprs(gse[[1]])
pheno <- pData(gse[[1]])

## ---- 3. Definisikan Grup Sampel ----
group <- ifelse(
  grepl("^paired normal", pheno$title),
  "Paired_Normal",
  ifelse(
    grepl("^breast cancer", pheno$title),
    "Breast_Cancer",
    "Apparent_Normal"
  )
)

# Hanya gunakan sampel Breast_Cancer dan Paired_Normal
keep <- group != "Apparent_Normal"

expr <- expr_full[, keep]
pheno_sub <- pheno[keep, ]

group_sub <- factor(
  group[keep],
  levels = c("Paired_Normal", "Breast_Cancer")
)

stopifnot(all(colnames(expr) == rownames(pheno_sub)))

## ---- 4. Differential Expression Analysis (limma) ----
design <- model.matrix(~0 + group_sub)
colnames(design) <- levels(group_sub)

fit <- lmFit(expr, design)

contrast_matrix <- makeContrasts(
  Paired_Normal - Breast_Cancer,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

deg <- topTable(
  fit2,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Kategorikan status regulasi gen
deg$Status <- "Not Significant"
deg$Status[deg$adj.P.Val < 0.05 & deg$logFC >= 1] <- "Upregulated"
deg$Status[deg$adj.P.Val < 0.05 & deg$logFC <= -1] <- "Downregulated"

table(deg$Status)

## ---- 5. Volcano Plot ----
volcano <- ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Status), alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano Plot GSE139038",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = "Status"
  ) +
  theme_minimal()

ggsave("plots/volcano_plot.png", volcano, width = 8, height = 6, dpi = 300)

## ---- 6. Heatmap Top 50 Probe Signifikan ----
top50 <- deg[deg$Status != "Not Significant", ]
top50 <- top50[order(top50$adj.P.Val), ][1:50, ]

top50_expr <- expr[rownames(top50), ]
top50_scaled <- t(scale(t(top50_expr)))

annotation_col <- data.frame(Group = group_sub)
rownames(annotation_col) <- colnames(top50_scaled)

pheatmap(
  top50_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  scale = "none",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Probes",
  filename = "plots/heatmap_top50.png",
  width = 8,
  height = 8
)

## ---- 7. Anotasi Probe ke Gene Symbol ----
gpl <- getGEO("GPL27630")
gpl_table <- Table(gpl)

annotation_df <- gpl_table[, c("ID", "Gene_Symbol")]
colnames(annotation_df) <- c("Probe", "Gene")

deg$Probe <- rownames(deg)
deg_annotated <- merge(deg, annotation_df, by = "Probe", all.x = TRUE)

## ---- 8. Daftar Gen Signifikan (Up & Down) ----
deg_sig <- deg_annotated[deg_annotated$Status != "Not Significant", ]

up_genes <- unique(deg_sig$Gene[deg_sig$Status == "Upregulated"])
down_genes <- unique(deg_sig$Gene[deg_sig$Status == "Downregulated"])

# Bersihkan gen ambigu (overlap up & down akibat multi-probe)
overlap <- intersect(up_genes, down_genes)
up_genes <- setdiff(up_genes, overlap)
down_genes <- setdiff(down_genes, overlap)

up_genes <- up_genes[!is.na(up_genes) & up_genes != ""]
down_genes <- down_genes[!is.na(down_genes) & down_genes != ""]

gene_list <- unique(c(up_genes, down_genes))

## ---- 9. Simpan Hasil DEG ke CSV ----
write.csv(deg_annotated, "data/DEG_GSE139038_annotated.csv", row.names = FALSE)
write.csv(data.frame(Gene = up_genes), "data/upregulated_genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = down_genes), "data/downregulated_genes.csv", row.names = FALSE)

## ---- 10. GO Enrichment Analysis (gprofiler2) ----
go_up <- gost(
  query = up_genes,
  organism = "hsapiens",
  sources = c("GO:BP", "GO:MF", "GO:CC"),
  correction_method = "fdr",
  significant = TRUE
)

go_down <- gost(
  query = down_genes,
  organism = "hsapiens",
  sources = c("GO:BP", "GO:MF", "GO:CC"),
  correction_method = "fdr",
  significant = TRUE
)

# Plot GO:BP - Upregulated
go_up_bp <- subset(go_up$result, source == "GO:BP")
go_up_bp <- go_up_bp[order(go_up_bp$p_value), ][1:15, ]
go_up_bp$minus_log10_p <- -log10(go_up_bp$p_value)

ggplot(
  go_up_bp,
  aes(x = minus_log10_p, y = reorder(term_name, minus_log10_p), size = intersection_size)
) +
  geom_point() +
  labs(
    title = "GO Biological Process Enrichment - Upregulated Genes",
    x = "-log10(P-value)", y = "GO Biological Process", size = "Gene Count"
  ) +
  theme_bw()

ggsave("plots/GO_Enrichment_Upregulated.png", width = 10, height = 7, dpi = 300)

# Plot GO:BP - Downregulated
go_down_bp <- subset(go_down$result, source == "GO:BP")
go_down_bp <- go_down_bp[order(go_down_bp$p_value), ][1:15, ]
go_down_bp$minus_log10_p <- -log10(go_down_bp$p_value)

ggplot(
  go_down_bp,
  aes(x = minus_log10_p, y = reorder(term_name, minus_log10_p), size = intersection_size)
) +
  geom_point() +
  labs(
    title = "GO Biological Process Enrichment - Downregulated Genes",
    x = "-log10(P-value)", y = "GO Biological Process", size = "Gene Count"
  ) +
  theme_bw()

ggsave("plots/GO_Enrichment_Downregulated.png", width = 10, height = 7, dpi = 300)

## ---- 11. KEGG Pathway Enrichment ----
kegg_up <- gost(
  query = up_genes,
  organism = "hsapiens",
  sources = "KEGG",
  correction_method = "g_SCS"
)

kegg_down <- gost(
  query = down_genes,
  organism = "hsapiens",
  sources = "KEGG",
  correction_method = "g_SCS"
)

# Plot KEGG - Upregulated
kegg_up_plot <- kegg_up$result[order(kegg_up$result$p_value), ][1:20, ]
kegg_up_plot <- kegg_up_plot[!is.na(kegg_up_plot$p_value), ]
kegg_up_plot$term_name <- factor(kegg_up_plot$term_name, levels = rev(kegg_up_plot$term_name))

ggplot(kegg_up_plot, aes(x = -log10(p_value), y = term_name)) +
  geom_point(aes(size = intersection_size)) +
  labs(
    title = "KEGG Pathway Enrichment - Upregulated Genes",
    x = "-log10(p-value)", y = "KEGG Pathway", size = "Gene Count"
  ) +
  theme_minimal()

ggsave("plots/KEGG_Enrichment_Upregulated.png", width = 10, height = 7, dpi = 300)

# Plot KEGG - Downregulated
kegg_down_plot <- kegg_down$result[order(kegg_down$result$p_value), ][1:20, ]
kegg_down_plot <- kegg_down_plot[!is.na(kegg_down_plot$p_value), ]
kegg_down_plot$term_name <- factor(kegg_down_plot$term_name, levels = rev(kegg_down_plot$term_name))

ggplot(kegg_down_plot, aes(x = -log10(p_value), y = term_name)) +
  geom_point(aes(size = intersection_size)) +
  labs(
    title = "KEGG Pathway Enrichment - Downregulated Genes",
    x = "-log10(p-value)", y = "KEGG Pathway", size = "Gene Count"
  ) +
  theme_minimal()

ggsave("plots/KEGG_Enrichment_Downregulated.png", width = 10, height = 7, dpi = 300)

## ============================================================
## Selesai
## ============================================================
