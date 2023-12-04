library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
setwd("/home/park/project/graduate/GSE200639_tuberculosis/")

rm(list = ls())
# Read 10X files and create Seurat objects
for (file in c("Normal_1","Normal_2")){
  data <- Read10X(paste0("GSE200639_RAW/",file))
  data_obj <- CreateSeuratObject(data, min.features = 200, min.cells = 3, project = "Uninfected")
  data_obj[["percent.mt"]] <- PercentageFeatureSet(data_obj, pattern = "^mt-")
  assign(file, data_obj)
}
for (file in c("D50_1", "D50_2", "D50_3")){
  data <- Read10X(paste0("GSE200639_RAW/",file))
  data_obj <- CreateSeuratObject(data, min.features = 200, min.cells = 3, project = "D50")
  data_obj[["percent.mt"]] <- PercentageFeatureSet(data_obj, pattern = "^mt-")
  assign(file, data_obj)
}
for (file in c("D100_1", "D100_2", "D100_3")){
  data <- Read10X(paste0("GSE200639_RAW/",file))
  data_obj <- CreateSeuratObject(data, min.features = 200, min.cells = 3, project = "D100")
  data_obj[["percent.mt"]] <- PercentageFeatureSet(data_obj, pattern = "^mt-")
  assign(file, data_obj)
}
Normal <- merge(x = Normal_1, y = c(Normal_2), add.cell.id = c("Normal_1","Normal_2"))
Normal <- AddMetaData(Normal, metadata = "Uninfected", col.name = "Condition")
D50 <- merge(x = D50_1, y = c(D50_2, D50_3), add.cell.id = c("D50_1","D50_2", "D50_3"))
D50 <- AddMetaData(D50, "Infected", "Condition")
D100 <- merge(x = D100_1, y = c(D100_2, D100_3), add.cell.ids = c("D100_1", "D100_2", "D100_3"))
D100 <- AddMetaData(D100, "Infected", "Condition")


# Data filtering and doublet removal
Normal <- subset(Normal, subset = nFeature_RNA > 1100 & nFeature_RNA < 4500 & percent.mt < 5)
D50 <- subset(D50, subset = nFeature_RNA > 1100 & nFeature_RNA < 4500 & percent.mt < 5)
D100 <- subset(D100, subset = nFeature_RNA > 1100 & nFeature_RNA < 4500 & percent.mt < 5)


datasets <- list(Normal, D50, D100)
datasets <- lapply(X = datasets, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, vars.to.regress = "percent.mt")
  x <- RunPCA(x, npcs = 17, verbose = F)
  x <- RunUMAP(x, dims = 1:17)
  x <- FindNeighbors(x, dims = 1:17)
  x <- FindClusters(x, resolution = 0.8)
})

library(scDblFinder)
Nor <- scDblFinder(GetAssayData(datasets[[1]], slot = "counts"), clusters = Idents(datasets[[1]])) # Threshold - 0.372
datasets[[1]]$scDblFinder.score <- Nor$scDblFinder.score
FeaturePlot(datasets[[1]], features = "scDblFinder.score")
Normal <- subset(datasets[[1]], subset = scDblFinder.score < 0.372)

D50 <- scDblFinder(GetAssayData(datasets[[2]], slot = "counts"), clusters = Idents(datasets[[2]])) # Threshold - 0.387
datasets[[2]]$scDblFinder.score <- D50$scDblFinder.score
FeaturePlot(datasets[[2]], features = "scDblFinder.score")
D50 <- subset(datasets[[2]], subset = scDblFinder.score < 0.387)


D100 <- scDblFinder(GetAssayData(datasets[[3]], slot = "counts"), clusters = Idents(datasets[[3]])) # Threshold - 0.334
datasets[[3]]$scDblFinder.score <- D100$scDblFinder.score
FeaturePlot(datasets[[3]], features = "scDblFinder.score")
D100 <- subset(datasets[[3]], subset = scDblFinder.score < 0.334)



# Data intergration
immune.anchors <- FindIntegrationAnchors(object.list = list(Normal, D50, D100), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
immune.combined <- ScaleData(immune.combined, vars.to.regress = 'percent.mt')

# Dimensional reduction
immune.combined <- RunPCA(immune.combined, npcs = 17)
ElbowPlot(immune.combined)
immune.combined <- RunUMAP(immune.combined, dims = 1:17)

# Clustering
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:17)
immune.combined <- FindClusters(immune.combined, resolution = c(0.8, 1.4), algorithm = 4)

immune.combined$day <- factor(immune.combined$orig.ident, levels = c("Uninfected", "D50", "D100"))

UMAPPlot(immune.combined, group.by = "seurat_clusters", split.by = "day") + labs(title ="")
UMAPPlot(immune.combined, group.by = "day") + labs(title ="")
UMAPPlot(immune.combined, group.by = "seurat_clusters", label = T, repel = T) + labs(title ="")

saveRDS(immune.combined, "integrated_datasets.rds")