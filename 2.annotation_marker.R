library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
setwd("/home/park/project/graduate/GSE200639_tuberculosis/")

rm(list = ls())
immune.combined <- readRDS("integrated_datasets.rds")
immune.combined$integrated_snn_res.0.8 <- factor(immune.combined$integrated_snn_res.0.8, levels = rep(1:23))
Idents(immune.combined) <- 'integrated_snn_res.0.8'
DefaultAssay(immune.combined) <- "RNA"
############# Finding marker genes ################
cluster.all.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(cluster.all.markers, "cluster.all.markers.csv")
library(openxlsx)
cluster_names <- unique(cluster.all.markers$cluster)
wb <- createWorkbook()
for (cluster in cluster_names) {
  subset_markers <- cluster.all.markers[cluster.all.markers$cluster == cluster, ]
  sheet_name <- paste("Cluster", cluster)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, subset_markers)
}
saveWorkbook(wb, file = "cluster_markers.xlsx")
################ Marker based manual annotation #####################
marker.cluster.ids <- c("NK", "B", "B", "Cd4+ T", "Cd8+ T", "B", 'Cd4+ T', 'Cd8+ T', 'Cd8+ T', 'Cd4+ T', 'Cd8+ T', 'Cd4+ T', 'Cd8+ T', 'Cd8+ T', 'Cd4+ T', 'B', 'Unknown', 'Myeloid', 'Non-Immune', 'Non-Immune', 'Non-Immune', 'Non-Immune', 'Plasma')
names(marker.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, marker.cluster.ids)
DimPlot(immune.combined, label = T)
immune.combined$celltype.marker <- Idents(immune.combined)
################ Marker based CIPR annotation #####################
Idents(immune.combined) <- 'integrated_snn_res.0.8'
CIPR.cluster.ids <- c("NK", "B", "B", "Cd4+ T", "Cd8+ T", "B", 'Cd4+ T', 'Cd8+ T', 'Cd8+ T', 'Cd4+ T', 'Cd8+ T', 'Cd8+ T', 'Cd8+ T', 'Cd8+ T', 'Cd4+ T', 'B', 'Unknown', 'Myeloid', 'Non-Immune', 'Myeloid', 'Non-Immune', 'Myeloid', 'Myeloid')
names(CIPR.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, CIPR.cluster.ids)
DimPlot(immune.combined, label = T)
immune.combined$celltype.CIPR <- Idents(immune.combined)


saveRDS(immune.combined, "celltype_marker.rds")
