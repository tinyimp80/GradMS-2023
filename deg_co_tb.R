Idents(subset) <- "Condition"
DEG_CO_TB <- FindMarkers(subset,ident.1 = "Infected", ident.2 = "Uninfected" ,logfc.threshold = 0.25, test.use = "MAST")

volcano_CO_TB <- EnhancedVolcano(DEG_CO_TB,
                                 lab = rownames(DEG_CO_TB),
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj',drawConnectors = T, ylab = bquote(~-Log[10] ~ adj ~P),
                                 title = '', subtitle = '', caption = paste("Number of DEGs:", nrow(DEG_CO_TB)), pCutoff = 0.05, legendPosition = "none", labSize = 3, FCcutoff = 0.5, xlab = bquote(~Log[2] ~ TB/Uninfected ~Exp.))
ggsave("all.uninf_vs_TB.volcano.png", volcano_CO_TB, height = 8, width = 6.8)


DEG_CO_TB <- DEG_CO_TB[DEG_CO_TB$p_val_adj<0.05,]
DEG_CO_TB$rank <- -log(DEG_CO_TB$p_val+1e-312) * DEG_CO_TB$avg_log2FC
gsea_CO_TB <- DEG_CO_TB$rank
names(gsea_CO_TB) <- rownames(DEG_CO_TB)
fgseaRes_CO_TB <- fgsea(pathways = fgsea_H,
                        stats    = gsea_CO_TB,
                        minSize  = 15,
                        maxSize  = 500)
gseaplot <- ggplot(fgseaRes_CO_TB, aes(reorder(pathway, NES), NES)) + geom_col(aes(fill =padj<0.05) ) + coord_flip() +
  theme_minimal() + xlab("")
