setwd("/home/park/project/graduate/GSE200639_tuberculosis/")
rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(stringr)
library(msigdbr)
fgsea_H = msigdbr(species = "mouse", category = "H")
library(fgsea)
fgsea_H <- fgsea_H %>% split(x = .$gene_symbol, f = .$gs_name)
library(EnhancedVolcano)

immune.combined <- readRDS("final_annotation.rds")
pal <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
         "#44AA99", "#999933", "#888888")

setwd("/home/park/project/graduate/GSE200639_tuberculosis/output/deg_gsea/")

subset <- subset(immune.combined, subset = celltype=="B"|celltype=="Plasma"|celltype=="Cd4+ T"|celltype=="Cd8+ T"|celltype=="NKT"|celltype=="NK")
subset$sample <- sub("^(\\w+_\\d+).*", "\\1",rownames(subset@meta.data))
DefaultAssay(subset) <- "RNA"
Idents(subset) <- "day"
DEG_un_50 <- FindMarkers(subset,ident.1 = "D50", ident.2 = "Uninfected" ,logfc.threshold = 0.25, test.use = "MAST")
DEG_un_100 <- FindMarkers(subset,ident.1 = "D100", ident.2 = "Uninfected" ,logfc.threshold = 0.25, test.use = "MAST")
DEG_50_100 <- FindMarkers(subset,ident.1 = "D100", ident.2 = "D50" ,logfc.threshold = 0.25, test.use = "MAST")

library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "uninf_vs_D50")
writeData(wb, sheet = "uninf_vs_D50", DEG_un_50, rowNames = T)
addWorksheet(wb, "uninf_vs_D100")
writeData(wb, sheet = "uninf_vs_D100", DEG_un_100, rowNames = T)
addWorksheet(wb, "D100_vs_D50")
writeData(wb, sheet = "D100_vs_D50", DEG_50_100, rowNames = T)
saveWorkbook(wb, file = "all.deg.xlsx")

volcano_un_50 <- EnhancedVolcano(DEG_un_50,
                                 lab = rownames(DEG_un_50),
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj',drawConnectors = T, ylab = bquote(~-Log[10] ~ adj ~P),
                                 title = '', subtitle = '', caption = paste("Number of DEGs:", nrow(DEG_un_50)), pCutoff = 0.05, legendPosition = "none", labSize = 3,FCcutoff = 0.5,xlab = bquote(~Log[2] ~ D50/Uninfected ~Exp.))
ggsave("all.uninf_vs_D50.volcano.png", volcano_un_50, height = 8, width = 6.8)
volcano_un_100 <- EnhancedVolcano(DEG_un_100,
                                  lab = rownames(DEG_un_100),
                                  x = 'avg_log2FC',
                                  y = 'p_val_adj',drawConnectors = T, ylab = bquote(~-Log[10] ~ adj ~P),
                                  title = '', subtitle = '', caption = paste("Number of DEGs:", nrow(DEG_un_100)), pCutoff = 0.05, legendPosition = "none", labSize = 3,FCcutoff = 0.5,xlab = bquote(~Log[2] ~ D100/Uninfected ~Exp.))
ggsave("all.uninf_vs_D100.volcano.png", volcano_un_100, height = 8, width = 6.8)
volcano_50_100 <- EnhancedVolcano(DEG_50_100,
                                  lab = rownames(DEG_50_100),
                                  x = 'avg_log2FC',
                                  y = 'p_val_adj',drawConnectors = T, ylab = bquote(~-Log[10] ~ adj ~P),
                                  title = '', subtitle = '', caption = paste("Number of DEGs:", nrow(DEG_50_100)), pCutoff = 0.05, legendPosition = "none", labSize = 3,FCcutoff = 0.5,xlab = bquote(~Log[2] ~ D100/D50 ~Exp.))
ggsave("all.D50_vs_D100.volcano.png", volcano_50_100, height = 8, width = 6.8)

DEG_un_50 <- DEG_un_50[DEG_un_50$p_val_adj<0.05,]
DEG_un_50$rank <- -log(DEG_un_50$p_val+1e-312) * DEG_un_50$avg_log2FC
gsea_un_50 <- DEG_un_50$rank
names(gsea_un_50) <- rownames(DEG_un_50)

DEG_un_100 <- DEG_un_100[DEG_un_100$p_val_adj<0.05,]
DEG_un_100$rank <- -log(DEG_un_100$p_val+1e-312) * DEG_un_100$avg_log2FC
gsea_un_100 <- DEG_un_100$rank
names(gsea_un_100) <- rownames(DEG_un_100)

DEG_50_100 <- DEG_50_100[DEG_50_100$p_val_adj<0.05,]
DEG_50_100$rank <- -log(DEG_50_100$p_val+1e-312) * DEG_50_100$avg_log2FC
gsea_50_100 <- DEG_50_100$rank
names(gsea_50_100) <- rownames(DEG_50_100)

fgseaRes_un_50 <- fgsea(pathways = fgsea_H,
                        stats    = gsea_un_50,
                        minSize  = 15,
                        maxSize  = 500)
fgseaRes_un_100 <- fgsea(pathways = fgsea_H,
                         stats    = gsea_un_100,
                         minSize  = 15,
                         maxSize  = 500)
fgseaRes_50_100 <- fgsea(pathways = fgsea_H,
                         stats    = gsea_50_100,
                         minSize  = 15,
                         maxSize  = 500)

wb <- createWorkbook()
addWorksheet(wb, "uninf_vs_D50")
writeData(wb, sheet = "uninf_vs_D50", fgseaRes_un_50, rowNames = T)
addWorksheet(wb, "uninf_vs_D100")
writeData(wb, sheet = "uninf_vs_D100", fgseaRes_un_100, rowNames = T)
addWorksheet(wb, "D50_vs_D100")
writeData(wb, sheet = "D50_vs_D100", fgseaRes_50_100, rowNames = T)
saveWorkbook(wb, file = "all.gsea.xlsx")


gseaplot <- ggplot(fgseaRes_un_50, aes(reorder(pathway, NES), NES)) + geom_col(aes(fill =padj<0.05) ) + coord_flip() +
  theme_minimal() + xlab("")
ggsave("all_un_50.gsea.png", gseaplot, bg =  "white", height = 8, width = 6.8)
gseaplot <- ggplot(fgseaRes_un_100, aes(reorder(pathway, NES), NES)) + geom_col(aes(fill =padj<0.05) ) + coord_flip() +
  theme_minimal() + xlab("")
ggsave("all_un_100.gsea.png", gseaplot, bg =  "white", height = 8, width = 6.8)
gseaplot <- ggplot(fgseaRes_50_100, aes(reorder(pathway, NES), NES)) + geom_col(aes(fill =padj<0.05) ) + coord_flip() +
  theme_minimal() + xlab("")
ggsave("all_50_100.gsea.png", gseaplot, bg =  "white", height = 8, width = 6.8)


################################################
for(i in levels(factor(subset$celltype))[c(1,3:6)]){
  cell_type <- str_replace(i, "\\/", "")
  cell <- subset(subset, subset = celltype==i)
  
  DEG_un_50 <- FindMarkers(cell,ident.1 = "D50", ident.2 = "Uninfected" ,logfc.threshold = 0.25, test.use = "MAST")
  DEG_un_100 <- FindMarkers(cell,ident.1 = "D100", ident.2 = "Uninfected" ,logfc.threshold = 0.25, test.use = "MAST")
  DEG_50_100 <- FindMarkers(cell,ident.1 = "D100", ident.2 = "D50" ,logfc.threshold = 0.25, test.use = "MAST")
  
  wb <- createWorkbook()
  addWorksheet(wb, "uninf_vs_D50")
  writeData(wb, sheet = "uninf_vs_D50", DEG_un_50, rowNames = T)
  addWorksheet(wb, "uninf_vs_D100")
  writeData(wb, sheet = "uninf_vs_D100", DEG_un_100, rowNames = T)
  addWorksheet(wb, "D50_vs_D100")
  writeData(wb, sheet = "D50_vs_D100", DEG_50_100, rowNames = T)
  saveWorkbook(wb, file = paste0(cell_type,".deg.xlsx"))
  
  volcano_un_50 <- EnhancedVolcano(DEG_un_50,
                             lab = rownames(DEG_un_50),
                             x = 'avg_log2FC',
                             y = 'p_val_adj',drawConnectors = T, ylab = bquote(~-Log[10] ~ adj ~P),
                             title = '', subtitle = '', caption = paste("Number of DEGs:", nrow(DEG_un_50)), pCutoff = 0.05, legendPosition = "none", labSize = 3,FCcutoff = 0.5, xlab = bquote(~Log[2] ~ D50/Uninfected ~Exp.))
  ggsave(paste0(cell_type,".uninf_vs_D50.volcano.png"), volcano_un_50, height = 8, width = 6.8)
  volcano_un_100 <- EnhancedVolcano(DEG_un_100,
                                   lab = rownames(DEG_un_100),
                                   x = 'avg_log2FC',
                                   y = 'p_val_adj',drawConnectors = T, ylab = bquote(~-Log[10] ~ adj ~P),
                                   title = '', subtitle = '', caption = paste("Number of DEGs:", nrow(DEG_un_100)), pCutoff = 0.05, legendPosition = "none", labSize = 3,FCcutoff = 0.5, xlab = bquote(~Log[2] ~ D100/Uninfected ~Exp.))
  ggsave(paste0(cell_type,".uninf_vs_D100.volcano.png"), volcano_un_100, height = 8, width = 6.8)
  volcano_50_100 <- EnhancedVolcano(DEG_50_100,
                                   lab = rownames(DEG_50_100),
                                   x = 'avg_log2FC',
                                   y = 'p_val_adj',drawConnectors = T, ylab = bquote(~-Log[10] ~ adj ~P),
                                   title = '', subtitle = '', caption = paste("Number of DEGs:", nrow(DEG_50_100)), pCutoff = 0.05, legendPosition = "none", labSize = 3,FCcutoff = 0.5, xlab = bquote(~Log[2] ~ D100/D50 ~Exp.))
  ggsave(paste0(cell_type,".D50_vs_D100.volcano.png"), volcano_50_100, height = 8, width = 6.8)
   
  #GSEA
  DEG_un_50 <- DEG_un_50[DEG_un_50$p_val_adj <0.05,]
  DEG_un_50$rank <- -log(DEG_un_50$p_val+1e-312) * DEG_un_50$avg_log2FC
  gsea_un_50 <- DEG_un_50$rank
  names(gsea_un_50) <- rownames(DEG_un_50)
  
  DEG_un_100 <- DEG_un_100[DEG_un_100$p_val_adj<0.05,]
  DEG_un_100$rank <- -log(DEG_un_100$p_val+1e-312) * DEG_un_100$avg_log2FC
  gsea_un_100 <- DEG_un_100$rank
  names(gsea_un_100) <- rownames(DEG_un_100)
  
  DEG_50_100 <- DEG_50_100[DEG_50_100$p_val_adj<0.05,]
  DEG_50_100$rank <- -log(DEG_50_100$p_val+1e-312) * DEG_50_100$avg_log2FC
  gsea_50_100 <- DEG_50_100$rank
  names(gsea_50_100) <- rownames(DEG_50_100)
  
  fgseaRes_un_50 <- fgsea(pathways = fgsea_H,
                    stats    = gsea_un_50,
                    minSize  = 15,
                    maxSize  = 500)
  fgseaRes_un_100 <- fgsea(pathways = fgsea_H,
                    stats    = gsea_un_100,
                    minSize  = 15,
                    maxSize  = 500)
  fgseaRes_50_100 <- fgsea(pathways = fgsea_H,
                    stats    = gsea_50_100,
                    minSize  = 15,
                    maxSize  = 500)
  
  
  wb <- createWorkbook()
  addWorksheet(wb, "uninf_vs_D50")
  writeData(wb, sheet = "uninf_vs_D50", fgseaRes_un_50, rowNames = T)
  addWorksheet(wb, "uninf_vs_D100")
  writeData(wb, sheet = "uninf_vs_D100", fgseaRes_un_100, rowNames = T)
  addWorksheet(wb, "D50_vs_D100")
  writeData(wb, sheet = "D50_vs_D100", fgseaRes_50_100, rowNames = T)
  saveWorkbook(wb, file = paste0(cell_type,".gsea.xlsx"))
  

  gseaplot <- ggplot(fgseaRes_un_50, aes(reorder(pathway, NES), NES)) + geom_col(aes(fill =padj<0.05) ) + coord_flip() +
    theme_minimal() + xlab("")
  ggsave(paste0(cell_type, "_un_50.gsea.png"), gseaplot, bg =  "white", height = 8, width = 6.8)
  gseaplot <- ggplot(fgseaRes_un_100, aes(reorder(pathway, NES), NES)) + geom_col(aes(fill =padj<0.05) ) + coord_flip() +
    theme_minimal() + xlab("")
  ggsave(paste0(cell_type, "_un_100.gsea.png"), gseaplot, bg =  "white", height = 8, width = 6.8)
  gseaplot <- ggplot(fgseaRes_50_100, aes(reorder(pathway, NES), NES)) + geom_col(aes(fill =padj<0.05) ) + coord_flip() +
    theme_minimal() + xlab("")
  ggsave(paste0(cell_type, "_50_100.gsea.png"), gseaplot, bg =  "white", height = 8, width = 6.8)
}

###########################################################
table_samples_by_cell_type <- B@meta.data %>%
  dplyr::group_by(Condition, celltype) %>%
  dplyr::summarize(count = n()) %>%
  tidyr::spread(celltype, count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c("Condition", "total_cell_count", dplyr::everything()))
table_samples_by_cell_type

temp_labels <- B@meta.data %>%
  group_by(Condition) %>%
  tally()

table <- table_samples_by_cell_type %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'Condition') %>%
  mutate(sample = factor(Condition, levels = levels(B@meta.data$Condition)))

rep(temp_labels$n, 3)
table$percentage <- table$value/rep(temp_labels$n, 3)*100
table$sample <- table$Condition
table$Condition <- rep(c(rep("D100",3), rep("D50",3), rep("Uninfected",2)),6)
table <- table %>%
  mutate(Condition = factor(Condition, levels = c("Uninfected", "D50", "D100")))
table %>%
  ggplot( aes(x=variable, y=percentage, fill=Condition)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  xlab("")

table[table$variable=="Plasma",] %>%
  ggplot( aes(x=variable, y=percentage, fill=Condition)) +
  geom_boxplot() +
  xlab("")






my_comparisons <- list( c("Uninfected", "D50"), c("Uninfected", "D100"), c("D50", "D100") )
ggboxplot(table, x = "Condition", y = "percentage",
          color = "Condition", palette = "npg", add = "jitter", facet.by = "variable")+stat_compare_means(comparisons = my_comparisons, symnum.args=symnum.args)

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
###########################################################