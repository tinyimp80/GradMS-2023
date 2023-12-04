library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
setwd("/home/park/project/graduate/GSE200639_tuberculosis/")

rm(list = ls())
immune.combined <- readRDS("final_annotation.rds")

######## Cd4+ #######
cd4 <- subset(immune.combined, subset = celltype.minor=="Cd4+ T")
# Naive T cell <- C4,10,15,25
p1 <- FeaturePlot(cd4, c("Ccr7", "Sell"), label = T, repel = T)
p2 <- VlnPlot(cd4, c("Ccr7", "Sell"), pt.size = 0)+NoLegend()
p1/p2
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- FindSubCluster(immune.combined, 7, algorithm = 4, resolution = 0.5,graph.name = "integrated_snn")
DimPlot(immune.combined, group.by = "sub.cluster", label = T, repel = T)
levels(immune.combined$seurat_clusters) <- rep(1:25)
immune.combined$seurat_clusters[immune.combined$sub.cluster=="7_6"] <- as.factor("25")
Idents(immune.combined) <- "seurat_clusters"
DimPlot(immune.combined, label = T)
DefaultAssay(immune.combined) <- "RNA"
cd4 <- subset(immune.combined, subset = celltype.minor=="Cd4+ T")
p1 <- FeaturePlot(cd4, c("Ccr7", "Sell"), label = T, repel = T)
p2 <- VlnPlot(cd4, c("Ccr7", "Sell"), pt.size = 0)+NoLegend()
DotPlot(cd4, features =c("Ccr7", "Sell"))
p1/p2

# Active T cell <- C7,12
p1 <- FeaturePlot(cd4, c("Ifng", "Rora"), label = T, repel = T, min.cutoff = "q10")
p2 <- VlnPlot(cd4, c("Ifng", "Rora"), pt.size = 0)+NoLegend()
DotPlot(cd4, features =c("Ifng", "Rora"))
p1/p2

# IFN Cd4+ <- C10
p1 <- FeaturePlot(cd4, c("Igtp", "Stat1"), label = T, repel = T, min.cutoff = "q10")
p2 <- VlnPlot(cd4, c("Igtp", "Stat1"), pt.size = 0)+NoLegend()
DotPlot(cd4, features =c("Igtp", "Stat1"))
p1/p2

######## Cd8+ #######
cd8 <- subset(immune.combined, subset = celltype.minor=="Cd8+ T")
# Effector Cd8+ <- ###############################################
p1 <- FeaturePlot(immune.combined, c("Igtp", "Stat1"), label = T, repel = T)
p2 <- VlnPlot(immune.combined, c("Igtp", "Stat1"), pt.size = 0)+NoLegend()
p1/p2
FeaturePlot(immune.combined, features = "Ighv")
Ighd

cluster7 <- subset(immune.combined, subset = seurat_clusters == 7)
DefaultAssay(cluster7) <- "integrated"
cluster7 <- RunPCA(cluster7, npcs = 30)
cluster7 <- RunUMAP(cluster7, dims = 1:30)
cluster7 <- FindNeighbors(cluster7)
cluster7 <- FindClusters(cluster7, algorithm = 4, resolution = 1)
DimPlot(cluster7, label = T, repel = T, group.by = "seurat_clusters", reduction = "umap")
DimPlot(cluster7, group.by = "sub.cluster", label = T)
DefaultAssay(cluster7) <- "RNA"
FeaturePlot(cluster7, c("Ccr7", "Sell"), label = T, repel = T)









#################### Assignment ##############
Idents(immune.combined) <- "seurat_clusters"
sub.cluster.ids <- c("NK", "B", "B", "Cd4+ T", "Cd8+ T", "B", "Cd4+ T", "Cd8+ T", "Cd8+ T", "Cd4+ T", "Cd8+ T", "Cd4+ T", "Cd8+ T", "Cd8+ T", "Cd4+ T", "B", "Cd4+/Cd8+ doublet", "Myeloid", "Non-Immune", "NK/T doublet", "Non-Immune", "B", "Plasma", "NKT")
names(sub.cluster.ids) <- levels(sub.cluster.ids)
immune.combined <- RenameIdents(immune.combined, sub.cluster.ids)
DimPlot(immune.combined, label = T)+ggtitle("")
immune.combined$celltype.sub <- Idents(immune.combined)