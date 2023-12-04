library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
setwd("/home/park/project/graduate/GSE200639_tuberculosis/")

rm(list = ls())
immune.combined <- readRDS("celltype_marker.rds")
DefaultAssay(immune.combined) <- "RNA"
Idents(immune.combined) <- 'integrated_snn_res.0.8'
# Cell type annotation
# Immune cell marker
p1 <- FeaturePlot(immune.combined, "Ptprc", label = T, repel = T) # C19,21 : non-immune cell
p2 <- VlnPlot(immune.combined, "Ptprc", pt.size = 0)+NoLegend()
p3 <- DotPlot(immune.combined, features = "Ptprc")
p1/p2
p3
Immune.cluster.ids <- c(rep("Immune", 18), "Non-Immune", "Immune", "Non-Immune", "Immune", "Immune")
names(Immune.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, Immune.cluster.ids)
immune.combined$immune.cell <- Idents(immune.combined)
DimPlot(immune.combined, group.by = "immune.cell") + ggtitle("")
# Myeloid cell marker
Idents(immune.combined) <- "integrated_snn_res.0.8"
p1 <- FeaturePlot(immune.combined, "Itgam", label = T, repel = T) # C1,18 : Myeloid cell
p2 <- VlnPlot(immune.combined, "Itgam", pt.size = 0)+NoLegend()
p1/p2
# T cell marker
p1 <- FeaturePlot(immune.combined, c("Cd3d", "Cd3e"), label = T, repel = T) # C4,5,7,8,9,10,11,12,13,14,15,17,20
p2 <- VlnPlot(immune.combined, c("Cd3d", "Cd3e"), pt.size = 0)+NoLegend()
DotPlot(immune.combined, features =c("Cd3d", "Cd3e"))
p1/p2
# CD4+ T cell
p1 <- FeaturePlot(immune.combined, "Cd4", label = T, repel = T) # C4,7,10,12,15,17,20
p2 <- VlnPlot(immune.combined, c("Cd4"), pt.size = 0)+NoLegend()
DotPlot(immune.combined, features =c("Cd4"))
p1/p2
# CD8+ T cell
p1 <- FeaturePlot(immune.combined, c("Cd8b1","Cd8a"), label = T, repel = T) # C5,8,9,11,13,14,17,20
p2 <- VlnPlot(immune.combined, c("Cd8b1","Cd8a"), pt.size = 0)+NoLegend()
DotPlot(immune.combined, features =c("Cd8a", "Cd8b1"))
p1/p2

immune.combined <- FindSubCluster(immune.combined, 17, algorithm = 4, graph.name = "integrated_snn", resolution = 1)
UMAPPlot(immune.combined, group.by = "sub.cluster", label = T, repel = T) + labs(title ="")
immune.combined$seurat_clusters <- immune.combined$integrated_snn_res.0.8
immune.combined$seurat_clusters[immune.combined$sub.cluster=="17_5"] <- 1
immune.combined$seurat_clusters[immune.combined$sub.cluster=="17_6"] <- 1
immune.combined$seurat_clusters[immune.combined$sub.cluster=="17_4"] <- 6
UMAPPlot(immune.combined, group.by = "seurat_clusters", label = T, repel = T) + labs(title ="")
# B cell marker
p1 <- FeaturePlot(immune.combined, c("Cd79a"), label = T, repel = T) # C2,3,6,16,22,23
p2 <- VlnPlot(immune.combined, c("Cd79a"), pt.size = 0)+NoLegend()
p1/p2
# B/Plasma cell
p1 <- FeaturePlot(immune.combined, c("Ms4a1", "Mzb1"), label = T, repel = T) # C23이 plasma
p2 <- VlnPlot(immune.combined, c("Ms4a1", "Mzb1"), pt.size = 0)+NoLegend()
p1/p2
# NK cell marker
p1 <- FeaturePlot(immune.combined, c("Gzmb", "Ncr1"), label = T, repel = T) # C1, 11일부, 20일부
p2 <- VlnPlot(immune.combined, c("Gzmb", "Ncr1"), pt.size = 0)+NoLegend()
p1/p2


immune.combined <- FindSubCluster(immune.combined, 20, algorithm = 4, graph.name = "integrated_snn", resolution = 1)
UMAPPlot(immune.combined, group.by = "sub.cluster", label = T, repel = T) + labs(title ="")
immune.combined$seurat_clusters[immune.combined$sub.cluster=="20_2"] <- 1
UMAPPlot(immune.combined, group.by = "seurat_clusters", label = T, repel = T) + labs(title ="")


immune.combined <- FindSubCluster(immune.combined, 11, algorithm = 4, graph.name = "integrated_snn", resolution = 1)
DimPlot(immune.combined, group.by = "sub.cluster", label = T)
cluster11 <- subset(immune.combined, subset = seurat_clusters == 11)
DimPlot(cluster11, group.by = "sub.cluster", label = T)
VlnPlot(cluster11, c("Cd8b1","Cd8a"), group.by = "sub.cluster", pt.size = 0)+NoLegend()
VlnPlot(cluster11, c("Gzmb", "Ncr1"), group.by = "sub.cluster", pt.size = 0)+NoLegend()

levels(immune.combined$seurat_clusters) <- rep(1:24)
immune.combined$seurat_clusters[immune.combined$sub.cluster=="11_6"] <- as.factor("24")
Idents(immune.combined) <- "seurat_clusters"
DimPlot(immune.combined, label = T)
p1 <- FeaturePlot(immune.combined, c("Cd8b1","Cd8a"), label = T, repel = T) 
p2 <- VlnPlot(immune.combined, c("Cd8b1","Cd8a"), pt.size = 0)+NoLegend()
p1/p2

p1 <- FeaturePlot(immune.combined, c("Gzmb", "Ncr1"), label = T, repel = T)
p2 <- VlnPlot(immune.combined, c("Gzmb", "Ncr1"), pt.size = 0)+NoLegend()
p1/p2


###################### Assignment ###########################
DimPlot(immune.combined, group.by = "seurat_clusters", label = T)+ggtitle("")
# Major type
Idents(immune.combined) <- "seurat_clusters"
major.cluster.ids <- c("NK", "B", "B", "T", "T", "B", "T", "T", "T", "T", "T", "T", "T", "T", "T", "B", "T", "Myeloid", "Non-Immune", "T", "Non-Immune", "B", "B", "NK")
names(major.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, major.cluster.ids)
DimPlot(immune.combined, label = T)+ggtitle("")
immune.combined$celltype.major <- Idents(immune.combined)
# Minor type
Idents(immune.combined) <- "seurat_clusters"
minor.cluster.ids <- c("NK", "B", "B", "Cd4+ T", "Cd8+ T", "B", "Cd4+ T", "Cd8+ T", "Cd8+ T", "Cd4+ T", "Cd8+ T", "Cd4+ T", "Cd8+ T", "Cd8+ T", "Cd4+ T", "B", "Cd4+/Cd8+ doublet", "Myeloid", "Non-Immune", "NK/T doublet", "Non-Immune", "B", "Plasma", "NKT")
names(minor.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, minor.cluster.ids)
DimPlot(immune.combined, label = T)+ggtitle("")
immune.combined$celltype.minor <- Idents(immune.combined)
DimPlot(immune.combined, label = T, group.by = "celltype.minor")+ggtitle("")

saveRDS(immune.combined, "celltype_manual.rds")
