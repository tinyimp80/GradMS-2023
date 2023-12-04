table_samples_by_cell_type <- subset@meta.data %>%
  dplyr::group_by(sample, celltype) %>%
  dplyr::summarize(count = n()) %>%
  tidyr::spread(celltype, count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c("sample", "total_cell_count", dplyr::everything()))
table_samples_by_cell_type
table_samples_by_cell_type$day <- c(rep("D100", 3),rep("D50", 3),rep("Uninfected", 2))
table_samples_by_cell_type$day <- factor(table_samples_by_cell_type$day, levels=c("Uninfected", "D50", "D100"))
bcell_cluster_cell <- cd8@meta.data %>%
  dplyr::group_by(sample, seurat_clusters) %>%
  dplyr::summarize(count = n()) %>%
  tidyr::spread(seurat_clusters, count, fill = 0) %>%
  dplyr::ungroup()
bcell_cluster_cell <- cbind(table_samples_by_cell_type, bcell_cluster_cell[,-1])
bcell_cluster_cell$percent.1 <- bcell_cluster_cell$'1'/bcell_cluster_cell$total_cell_count
bcell_cluster_cell$percent.2 <- bcell_cluster_cell$'2'/bcell_cluster_cell$total_cell_count
bcell_cluster_cell$percent.3 <- bcell_cluster_cell$'3'/bcell_cluster_cell$total_cell_count
bcell_cluster_cell$percent.4 <- bcell_cluster_cell$'4'/bcell_cluster_cell$total_cell_count
bcell_cluster_cell$percent.5 <- bcell_cluster_cell$'5'/bcell_cluster_cell$total_cell_count
bcell_cluster_cell$percent.6 <- bcell_cluster_cell$'6'/bcell_cluster_cell$total_cell_count
bcell_cluster_cell$percent.7 <- bcell_cluster_cell$'7'/bcell_cluster_cell$total_cell_count


ggboxplot(bcell_cluster_cell, x = "day", y = "percent.5", fill = "day",
          ylab = "Proportion of cells in C5", xlab = "")+
  stat_compare_means(comparisons = list(c("Uninfected", "D50"), c("Uninfected", "D100"), c("D50", "D100")), method = "t.test")+  theme(legend.position="none")
table_samples_by_cell_type
table_samples_by_cell_type$percent.B <- table_samples_by_cell_type$B/table_samples_by_cell_type$total_cell_count
table_samples_by_cell_type$percent.Plasma <- table_samples_by_cell_type$Plasma/table_samples_by_cell_type$total_cell_count
table_samples_by_cell_type$percent.cd4 <- table_samples_by_cell_type$'Cd4+ T'/table_samples_by_cell_type$total_cell_count
table_samples_by_cell_type$percent.cd8 <- table_samples_by_cell_type$'Cd8+ T'/table_samples_by_cell_type$total_cell_count
table_samples_by_cell_type$percent.nkt <- table_samples_by_cell_type$NKT/table_samples_by_cell_type$total_cell_count
table_samples_by_cell_type$percent.nk <- table_samples_by_cell_type$NK/table_samples_by_cell_type$total_cell_count
ggboxplot(table_samples_by_cell_type, x = "day", y = "percent.nkt", fill = "day",
          ylab = "Proportion of cells in C5", xlab = "") +
  stat_compare_means(comparisons = list(c("Uninfected", "D50"), c("Uninfected", "D100"), c("D50", "D100")), method = "t.test")+  theme(legend.position="none")










library(clusterProfiler)
library(org.Mm.eg.db)
Idents(subset) <- "Condition"
all_up <- FindMarkers(subset, ident.1 = "Infected", ident.2 = "Uninfected", test.use = "MAST", only.pos = T, min.pct = 0.25)
all_up <- all_up[all_up$p_val_adj<0.05,]
ego_up <- enrichGO(rownames(all_up), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)

all_dn <- FindMarkers(subset, ident.2 = "Infected", ident.1 = "Uninfected", test.use = "MAST", only.pos = T, min.pct = 0.25)
all_dn <- all_dn[all_dn$p_val_adj<0.05,]
ego_dn <- enrichGO(rownames(all_dn), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
dotplot(ego_up, showCategory = 15)+ggtitle("TB Over-represent")

VlnPlot(immune.combined, features = c("Stat1", "Ly6a"), group.by = "day", pt.size = 0)+stat_compare_means()





vp_case1 <- function(gene_signature, file_name, test_sign){
  plot_case1 <- function(signature, y_max = 6){
    VlnPlot(immune.combined, features = signature,
            pt.size = 0.1,
            group.by = "day",
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + stat_compare_means(comparisons = test_sign, label = "p.format")
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1) )
  }
  cowplot::plot_grid(plotlist = plot_list)
  file_name <- paste0(file_name, "_r.png")
  ggsave(file_name, width = 14, height = 8)
}
vp_case1(c("Stat1", "Ly6a"), "output/sig_deg", list(c("Uninfected", "D50"), c("Uninfected", "D100"), c("D50", "D100")))
