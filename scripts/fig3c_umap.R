# 26c_fig3c.R (v2) - scRNA KLRF1 localization with celltype annotation
suppressMessages({library(Seurat); library(ggplot2)})
OUT <- "H:/data/.openclaw/workspace/sepsis_project/output/03_immunosenescence/20_stitch/figures_v1"
obj <- readRDS("H:/data/.openclaw/workspace/sepsis_project/output/07_single_cell/seurat_integrated.rds")

ann <- read.csv("H:/data/.openclaw/workspace/sepsis_project/output/07_single_cell/cluster_annotations.csv")
map <- setNames(ann$cell_type, as.character(ann$cluster))
obj$cell_type <- unname(map[as.character(obj$seurat_clusters)])
cat("cell_type table:\n"); print(table(obj$cell_type))

# UMAP by celltype
p_dim <- DimPlot(obj, group.by = "cell_type", label = TRUE, repel = TRUE, pt.size = 0.25, label.size = 3.2) +
  ggtitle("GSE167363 peripheral blood (24,000 cells)") +
  theme(plot.title = element_text(size = 10), legend.text = element_text(size = 7.5))
ggsave(file.path(OUT, "fig3c_umap_celltype.png"), p_dim, width = 6.4, height = 4.6, dpi = 300)

# KLRF1 feature plot
p_feat <- FeaturePlot(obj, features = "KLRF1", pt.size = 0.25, cols = c("grey92", "#C0392B"), order = TRUE) +
  ggtitle("KLRF1 expression") + theme(plot.title = element_text(size = 10))
ggsave(file.path(OUT, "fig3c_umap_klrf1.png"), p_feat, width = 5.0, height = 4.6, dpi = 300)

cat("fig3c done\n")
