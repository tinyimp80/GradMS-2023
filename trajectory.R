setwd("/home/park/project/graduate/GSE200639_tuberculosis/")
rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(stringr)
immune.combined <- readRDS("final_annotation.rds")
immune.combined$sample <- sub("^(\\w+_\\d+).*", "\\1", rownames(immune.combined@meta.data))

##################### cd4
cd4 <- subset(subset, subset = celltype=="Cd4+ T")
Idents(cd4) <- "day"
cd4 <- FindVariableFeatures(cd4, nfeatures = 1000)
cd4 <- ScaleData(cd4, vars.to.regress = "percent.mt")
cd4 <- RunPCA(cd4, npcs = 8)
cd4 <- RunUMAP(cd4, dims = 1:8)
cd4 <- FindNeighbors(cd4, dims = 1:8)
cd4 <- FindClusters(cd4, algorithm = 4, resolution = 0.2)
DimPlot(cd4, label = T)+NoLegend()
FeaturePlot(cd4, features = c("Ccr7", "Sell","Stat1", "Igtp"), min.cutoff = "q10")
# trajectory
cds <- SeuratWrappers::as.cell_data_set(cd4)
fData(cds)$gene_short_name <- rownames(fData(cds))
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

cds@clusters@listData[["UMAP"]][["clusters"]] <- cd4@active.ident
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- cd4@reductions$umap@cell.embeddings
cds <- preprocess_cds(cds)
cds <- learn_graph(cds, use_partition = T)
cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "cluster",
           label_cell_groups=FALSE,
           label_leaves=F,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           label_roots = F)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=T,
           label_leaves=F,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           label_roots = F)
cd4_c3 <- subset(cd4, subset = seurat_clusters==2)
Idents(cd4_c3) <- "Condition"
cd4_c3_up <- FindMarkers(cd4_c3, ident.1 = "Infected", ident.2 = "Uninfected", only.pos = T, min.pct = 0.25, test.use = "MAST")
library(clusterProfiler)
library(org.Mm.eg.db)
cd4_c3_up <- cd4_c3_up[cd4_c3_up$p_val_adj<0.05,]
cd4_c3_up_go <- enrichGO(rownames(cd4_c3_up), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
DimPlot(cd4_c3, split.by = "day")
cd4
cd4_marker <- FindAllMarkers(cd4, min.pct = 0.25, only.pos = T)
cluster_names <- unique(cd4_marker$cluster)
cd4_mk <- list()
for (i in cluster_names) {
  subset_markers <- cd4_marker[cd4_marker$cluster == i, ]
  subset_markers <- subset_markers[subset_markers$p_val_adj<0.05,]
  cd4_mk[[i]] <- subset_markers$gene
}
cd4_compare <- compareCluster(geneClusters = cd4_mk, fun = enrichGO, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "SYMBOL")
dotplot(cd4_compare)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")

FeaturePlot(cd4, "Ifng")
#######################
setwd("output/cd4_subcluster/")
Idents(cd4) <- "day"
for(i in levels(factor(cd4$seurat_clusters))){
  cell_type <- str_replace(i, "\\/", "")
  cell <- subset(cd4, subset = seurat_clusters==i)
  Idents(cell) <- "day"
  DEG_un_50 <- FindMarkers(cell,ident.1 = "D50", ident.2 = "Uninfected" ,logfc.threshold = 0.25, test.use = "MAST")
  DEG_un_100 <- FindMarkers(cell,ident.1 = "D100", ident.2 = "Uninfected" ,logfc.threshold = 0.25, test.use = "MAST")
  DEG_50_100 <- FindMarkers(cell,ident.1 = "D50", ident.2 = "D100" ,logfc.threshold = 0.25, test.use = "MAST")
  volcano_un_50 <- EnhancedVolcano(DEG_un_50,
                                   lab = rownames(DEG_un_50),
                                   x = 'avg_log2FC',
                                   y = 'p_val_adj',drawConnectors = T, ylab = bquote(~-Log[10] ~ adj ~P),
                                   title = '', subtitle = '', caption = 'cutoff : |Log2FC|>1, adj.P<0.05', pCutoff = 0.05, legendPosition = "none", labSize = 3)
  ggsave(paste0(cell_type,".uninf_vs_D50.volcano.png"), volcano_un_50, height = 8, width = 6.8)
  volcano_un_100 <- EnhancedVolcano(DEG_un_100,
                                    lab = rownames(DEG_un_100),
                                    x = 'avg_log2FC',
                                    y = 'p_val_adj',drawConnectors = T, ylab = bquote(~-Log[10] ~ adj ~P),
                                    title = '', subtitle = '', caption = 'cutoff : |Log2FC|>1, adj.P<0.05', pCutoff = 0.05, legendPosition = "none", labSize = 3)
  ggsave(paste0(cell_type,".uninf_vs_D100.volcano.png"), volcano_un_100, height = 8, width = 6.8)
  volcano_50_100 <- EnhancedVolcano(DEG_50_100,
                                    lab = rownames(DEG_50_100),
                                    x = 'avg_log2FC',
                                    y = 'p_val_adj',drawConnectors = T, ylab = bquote(~-Log[10] ~ adj ~P),
                                    title = '', subtitle = '', caption = 'cutoff : |Log2FC|>1, adj.P<0.05', pCutoff = 0.05, legendPosition = "none", labSize = 3)
  ggsave(paste0(cell_type,".D50_vs_D100.volcano.png"), volcano_50_100, height = 8, width = 6.8) 
  
}

cell <- subset(cd4, subset = seurat_clusters==2)
DEG_un_50 <- FindMarkers(cell,ident.1 = "D50", ident.2 = "Uninfected" ,logfc.threshold = 0.25, test.use = "MAST", only.pos = T)
DEG_un_50 <- DEG_un_50[DEG_un_50$p_val_adj<0.05,]
DEG_un_50$rank <- -log(DEG_un_50$p_val+1e-312) * DEG_un_50$avg_log2FC
gsea_un_50 <- DEG_un_50$rank
names(gsea_un_50) <- rownames(DEG_un_50)

enrichGO_un_50 <- enrichGO(rownames(DEG_un_50), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
dotplot(enrichGO_un_50, showCategory=15)

fgseaRes_un_50 <- fgsea(pathways = fgsea_H,
                        stats    = gsea_un_50,
                        minSize  = 15,
                        maxSize  = 500)
gseaplot <- ggplot(fgseaRes_un_50, aes(reorder(pathway, NES), NES)) + geom_col(aes(fill =padj<0.05) ) + coord_flip() +
  theme_minimal() + xlab("")
ggsave("2_un_50.gsea.png", gseaplot, bg =  "white", height = 8, width = 6.8)
#####################
table_samples_by_cell_type <- subset@meta.data %>%
  dplyr::group_by(sample, celltype) %>%
  dplyr::summarize(count = n()) %>%
  tidyr::spread(celltype, count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c("sample", "total_cell_count", dplyr::everything()))
table_samples_by_cell_type


table <- table_samples_by_cell_type %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'sample') %>%
  mutate(sample = factor(sample, levels = levels(factor(subset@meta.data$sample))))

rep(temp_labels$n, 3)
table$percentage <- table$value/rep(temp_labels$n, 3)*100
table$day <- rep(c(rep("D100",3), rep("D50",3), rep("Uninfected",2)),6)
table <- table %>%
  mutate(day = factor(day, levels = c("Uninfected", "D50", "D100")))
table %>%
  ggplot( aes(x=variable, y=percentage, fill=day)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  xlab("")

table[table$variable=="Plasma",] %>%
  ggplot( aes(x=variable, y=percentage, fill=day)) +
  geom_boxplot() +
  xlab("")






my_comparisons <- list( c("Uninfected", "D50"), c("Uninfected", "D100"), c("D50", "D100") )
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
ggboxplot(table[table$variable=="B",], x = "variable", y = "percentage",
          color = "day", palette = "npg")+stat_compare_means()

#############################cd8
cd8 <- subset(subset, subset = celltype=="Cd8+ T")
Idents(cd8) <- "day"
cd8 <- FindVariableFeatures(cd8, nfeatures = 1000)
cd8 <- ScaleData(cd8, vars.to.regress = "percent.mt")
cd8 <- RunPCA(cd8, npcs = 8)
cd8 <- RunUMAP(cd8, dims = 1:8)
cd8 <- FindNeighbors(cd8, dims = 1:8)
cd8 <- FindClusters(cd8, algorithm = 4, resolution = 0.2)
DimPlot(cd8, label = T)+NoLegend()

cd8_marker <- FindAllMarkers(cd8, min.pct = 0.25, only.pos = T)
cluster_names <- unique(cd8_marker$cluster)
cd8_mk <- list()
for (i in cluster_names) {
  subset_markers <- cd8_marker[cd8_marker$cluster == i, ]
  subset_markers <- subset_markers[subset_markers$p_val_adj<0.05,]
  cd8_mk[[i]] <- subset_markers$gene
}
cd8_compare <- compareCluster(geneClusters = cd8_mk, fun = enrichGO, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "SYMBOL")
dotplot(cd8_compare)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")

FeaturePlot(cd8, features = c("Ccr7", "Sell", "Cxcr3", "Ctla2a","Cd44" ,"Ifng"), min.cutoff = "q10")
# trajectory
cds <- SeuratWrappers::as.cell_data_set(cd8)
fData(cds)$gene_short_name <- rownames(fData(cds))
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

cds@clusters@listData[["UMAP"]][["clusters"]] <- cd8@active.ident
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- cd8@reductions$umap@cell.embeddings
cds <- preprocess_cds(cds)
cds <- learn_graph(cds, use_partition = T)
cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "cluster",
           label_cell_groups=FALSE,
           label_leaves=F,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           label_roots = F)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=T,
           label_leaves=F,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           label_roots = F)

setwd("../cd8_subcluster/")

cell <- subset(cd8, subset = seurat_clusters==5)
DEG_un_50 <- FindMarkers(cell,ident.1 = "D50", ident.2 = "Uninfected" ,logfc.threshold = 0.25, test.use = "MAST", only.pos = F)
DEG_un_50 <- DEG_un_50[DEG_un_50$p_val_adj<0.05,]
DEG_un_50$rank <- -log(DEG_un_50$p_val+1e-312) * DEG_un_50$avg_log2FC
gsea_un_50 <- DEG_un_50$rank
names(gsea_un_50) <- rownames(DEG_un_50)

enrichGO_un_50 <- enrichGO(rownames(DEG_un_50), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
dotplot(enrichGO_un_50, showCategory=15)

fgseaRes_un_50 <- fgsea(pathways = fgsea_H,
                        stats    = gsea_un_50,
                        minSize  = 15,
                        maxSize  = 500)
gseaplot <- ggplot(fgseaRes_un_50, aes(reorder(pathway, NES), NES)) + geom_col(aes(fill =padj<0.05) ) + coord_flip() +
  theme_minimal() + xlab("")
ggsave("2_un_50.gsea.png", gseaplot, bg =  "white", height = 8, width = 6.8)
############################B
bcell <- subset(immune.combined, subset = celltype=="B")
Idents(bcell) <- "day"
bcell <- FindVariableFeatures(bcell, nfeatures = 1000)
bcell <- ScaleData(bcell, vars.to.regress = "percent.mt")
bcell <- RunPCA(bcell, npcs = 8)
bcell <- RunUMAP(bcell, dims = 1:8)
bcell <- FindNeighbors(bcell, dims = 1:8)
bcell <- FindClusters(bcell, algorithm = 4, resolution = 0.3)
DimPlot(bcell, label = T)
FeaturePlot(bcell, c("H2-Ab1", "H2-Eb1", "H2-Dmb2"))
bcell_marker <- FindAllMarkers(bcell, min.pct = 0.25, only.pos = T)


cluster_names <- unique(bcell_marker$cluster)
b_mk <- list()
for (i in cluster_names) {
  subset_markers <- bcell_marker[bcell_marker$cluster == i, ]
  subset_markers <- subset_markers[subset_markers$p_val_adj<0.05,]
  b_mk[[i]] <- subset_markers$gene
}
b_compare <- compareCluster(geneClusters = b_mk, fun = enrichGO, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "SYMBOL")
dotplot(b_compare)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")
# trajectory
cds <- SeuratWrappers::as.cell_data_set(bcell)
fData(cds)$gene_short_name <- rownames(fData(cds))
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

cds@clusters@listData[["UMAP"]][["clusters"]] <- bcell@active.ident
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- bcell@reductions$umap@cell.embeddings
cds <- preprocess_cds(cds)
cds <- learn_graph(cds, use_partition = T)
cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "cluster",
           label_cell_groups=FALSE,
           label_leaves=F,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           label_roots = F)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=T,
           label_leaves=F,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           label_roots = F)
# geneset
setwd("../b_subcluster/")
Idents(bcell) <- "day"
cell <- subset(bcell, subset = seurat_clusters==1)
DEG_un_50 <- FindMarkers(cell,ident.1 = "D50", ident.2 = "Uninfected" ,logfc.threshold = 0.25, test.use = "MAST", only.pos = F)
DEG_un_50 <- DEG_un_50[DEG_un_50$p_val_adj<0.05,]
DEG_un_50$rank <- -log(DEG_un_50$p_val+1e-312) * DEG_un_50$avg_log2FC
gsea_un_50 <- DEG_un_50$rank
names(gsea_un_50) <- rownames(DEG_un_50)

enrichGO_un_50 <- enrichGO(rownames(DEG_un_50), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
dotplot(enrichGO_un_50, showCategory=15)

fgseaRes_un_50 <- fgsea(pathways = fgsea_H,
                        stats    = gsea_un_50,
                        minSize  = 15,
                        maxSize  = 500)
gseaplot <- ggplot(fgseaRes_un_50, aes(reorder(pathway, NES), NES)) + geom_col(aes(fill =padj<0.05) ) + coord_flip() +
  theme_minimal() + xlab("")
ggsave("2_un_50.gsea.png", gseaplot, bg =  "white", height = 8, width = 6.8)