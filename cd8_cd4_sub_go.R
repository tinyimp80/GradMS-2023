library(clusterProfiler)
library(org.Mm.eg.db)
cd8 <- subset(subset, subset = celltype=="Cd8+ T")
Idents(cd8) <- "day"
cd8 <- FindVariableFeatures(cd8, nfeatures = 1000)
cd8 <- ScaleData(cd8, vars.to.regress = "percent.mt")
cd8 <- RunPCA(cd8, npcs = 8)
cd8 <- RunUMAP(cd8, dims = 1:8)
cd8 <- FindNeighbors(cd8, dims = 1:8)
cd8 <- FindClusters(cd8, algorithm = 4, resolution = 0.2)
Idents(cd8) <- "day"
setwd("../cd8_subcluster/")
for(i in levels(factor(cd8$seurat_clusters))){
  cd8_sub <- subset(cd8, subset = seurat_clusters==i)
  DEG_un_50 <- FindMarkers(cd8_sub,ident.1 = "D50", ident.2 = "Uninfected" ,logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
  DEG_un_100 <- FindMarkers(cd8_sub,ident.1 = "D100", ident.2 = "Uninfected" ,logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
  Idents(cd8_sub) <- "Condition"
  DEG_un_TB <- FindMarkers(cd8_sub,ident.1 = "Infected", ident.2 = "Uninfected" ,logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
  # DEG_50_100 <- FindMarkers(cd8_sub,ident.1 = "D100", ident.2 = "D50" ,logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
  
  DEG_un_50 <- DEG_un_50[DEG_un_50$p_val_adj<0.05,]
  DEG_un_50$rank <- -log(DEG_un_50$p_val+1e-312) * DEG_un_50$avg_log2FC
  gsea_un_50 <- DEG_un_50$rank
  names(gsea_un_50) <- rownames(DEG_un_50)
    enrichGO_un_50 <- enrichGO(rownames(DEG_un_50), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
  enrichPlot <- dotplot(enrichGO_un_50, showCategory=15)
  ggsave(paste0("cd8_cluster",i,".uninf_vs_D50.GO.png"), enrichPlot, height = 12, width = 6.8)
  
  DEG_un_100 <- DEG_un_100[DEG_un_100$p_val_adj<0.05,]
  DEG_un_100$rank <- -log(DEG_un_100$p_val+1e-312) * DEG_un_100$avg_log2FC
  gsea_un_100 <- DEG_un_100$rank
  names(gsea_un_100) <- rownames(DEG_un_100)
  enrichGO_un_100 <- enrichGO(rownames(DEG_un_100), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
  enrichPlot <- dotplot(enrichGO_un_100, showCategory=15)
  ggsave(paste0("cd8_cluster",i,".uninf_vs_D100.GO.png"), enrichPlot, height = 12, width = 6.8)
  
  DEG_un_TB <- DEG_un_TB[DEG_un_TB$p_val_adj<0.05,]
  DEG_un_TB$rank <- -log(DEG_un_TB$p_val+1e-312) * DEG_un_TB$avg_log2FC
  gsea_un_TB <- DEG_un_TB$rank
  names(gsea_un_TB) <- rownames(DEG_un_TB)
  enrichGO_un_TB <- enrichGO(rownames(DEG_un_TB), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
  enrichPlot <- dotplot(enrichGO_un_TB, showCategory=15)
  ggsave(paste0("cd8_cluster",i,".uninf_vs_TB.GO.png"), enrichPlot, height = 12, width = 6.8)
}


cd4 <- subset(subset, subset = celltype=="Cd4+ T")
Idents(cd4) <- "day"
cd4 <- FindVariableFeatures(cd4, nfeatures = 1000)
cd4 <- ScaleData(cd4, vars.to.regress = "percent.mt")
cd4 <- RunPCA(cd4, npcs = 8)
cd4 <- RunUMAP(cd4, dims = 1:8)
cd4 <- FindNeighbors(cd4, dims = 1:8)
cd4 <- FindClusters(cd4, algorithm = 4, resolution = 0.2)
Idents(cd4) <- "day"
setwd("../cd4_subcluster/")
for(i in levels(factor(cd4$seurat_clusters))){
  cd4_sub <- subset(cd4, subset = seurat_clusters==i)
  DEG_un_50 <- FindMarkers(cd4_sub,ident.1 = "D50", ident.2 = "Uninfected" ,logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
  DEG_un_100 <- FindMarkers(cd4_sub,ident.1 = "D100", ident.2 = "Uninfected" ,logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
  Idents(cd4_sub) <- "Condition"
  DEG_un_TB <- FindMarkers(cd4_sub,ident.1 = "Infected", ident.2 = "Uninfected" ,logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
  # DEG_50_100 <- FindMarkers(cd4_sub,ident.1 = "D100", ident.2 = "D50" ,logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
  
  DEG_un_50 <- DEG_un_50[DEG_un_50$p_val_adj<0.05,]
  DEG_un_50$rank <- -log(DEG_un_50$p_val+1e-312) * DEG_un_50$avg_log2FC
  gsea_un_50 <- DEG_un_50$rank
  names(gsea_un_50) <- rownames(DEG_un_50)
  enrichGO_un_50 <- enrichGO(rownames(DEG_un_50), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
  enrichPlot <- dotplot(enrichGO_un_50, showCategory=15)
  ggsave(paste0("cd4_cluster",i,".uninf_vs_D50.GO.png"), enrichPlot, height = 12, width = 6.8)
  
  DEG_un_100 <- DEG_un_100[DEG_un_100$p_val_adj<0.05,]
  DEG_un_100$rank <- -log(DEG_un_100$p_val+1e-312) * DEG_un_100$avg_log2FC
  gsea_un_100 <- DEG_un_100$rank
  names(gsea_un_100) <- rownames(DEG_un_100)
  enrichGO_un_100 <- enrichGO(rownames(DEG_un_100), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
  enrichPlot <- dotplot(enrichGO_un_100, showCategory=15)
  ggsave(paste0("cd4_cluster",i,".uninf_vs_D100.GO.png"), enrichPlot, height = 12, width = 6.8)
  
  DEG_un_TB <- DEG_un_TB[DEG_un_TB$p_val_adj<0.05,]
  DEG_un_TB$rank <- -log(DEG_un_TB$p_val+1e-312) * DEG_un_TB$avg_log2FC
  gsea_un_TB <- DEG_un_TB$rank
  names(gsea_un_TB) <- rownames(DEG_un_TB)
  enrichGO_un_TB <- enrichGO(rownames(DEG_un_TB), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
  enrichPlot <- dotplot(enrichGO_un_TB, showCategory=15)
  ggsave(paste0("cd4_cluster",i,".uninf_vs_TB.GO.png"), enrichPlot, height = 12, width = 6.8)
}