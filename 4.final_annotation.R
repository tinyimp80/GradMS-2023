library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggforce)
library(scales)
setwd("/home/park/project/graduate/GSE200639_tuberculosis/")

rm(list = ls())
pal <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
         "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
immune.combined <- readRDS("celltype_manual.rds")
immune.combined$celltype.minor <- factor(immune.combined$celltype.minor, levels=c('B', 'Plasma', 'Cd4+ T', 'Cd8+ T',  'NKT', 'NK', 'Myeloid', 'Non-Immune','Cd4+/Cd8+ doublet', 'NK/T doublet'))
immune.combined$celltype.marker <- factor(immune.combined$celltype.marker, levels = c('B','Plasma', 'Cd4+ T', 'Cd8+ T', 'NK', 'Myeloid', 'Non-Immune', 'Unknown'))
immune.combined$celltype.CIPR <- factor(immune.combined$celltype.CIPR, levels = c('B', 'Cd4+ T', 'Cd8+ T', 'NK', 'Myeloid', 'Non-Immune', 'Unknown'))

manual <- levels(immune.combined$celltype.minor)
marker <- levels(immune.combined$celltype.marker)
CIPR <- levels(immune.combined$celltype.CIPR)
cluster <- levels(immune.combined$integrated_snn_res.0.8)

# 
# ########################### Manual vs Marker-base
# color_assignments <- setNames(
#   c(pal[1:length(manual)], pal[1:length(marker)]),
#   c(manual,marker)
# )
# 
# 
# data <- immune.combined@meta.data %>%
#   group_by(celltype.minor,celltype.marker) %>%
#   tally() %>%
#   ungroup() %>%
#   gather_set_data(1:2) %>%
#   dplyr::mutate(
#     x = factor(x, levels = unique(x)),
#     y = factor(y, levels = unique(y))
#   )
# 
# 
# data_labels <- tibble(
#   group = c(
#     rep('manual', length(manual)),
#     rep('marker', length(marker))
#   )
# ) %>%
#   mutate(
#     hjust = ifelse(group == 'sample', 1, 0),
#     nudge_x = ifelse(group == 'sample', -0.1, 0.1)
#   )
# 
# 
# 
# p1 <- ggplot(data, aes(x, id = id, split = y, value = n)) +
#   geom_parallel_sets(aes(fill = celltype.marker), alpha = 0.75, axis.width = 0.15) +
#   geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
#   geom_text(
#     aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
#     hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
#   ) +
#   scale_x_discrete(labels = c('Manual','Marker-based')) +
#   scale_fill_manual(values = color_assignments) +
#   theme_bw() +
#   theme(
#     legend.position = 'none',
#     axis.title = element_blank(),
#     axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
#     axis.text.y = element_blank(),
#     axis.ticks = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank()
#   )
# p1
# 
# 
# ########################### Manual vs CIPR
# color_assignments <- setNames(
#   c(pal[1:length(manual)], pal[1:length(CIPR)]),
#   c(manual,CIPR)
# )
# 
# 
# data <- immune.combined@meta.data %>%
#   group_by(celltype.minor,celltype.CIPR) %>%
#   tally() %>%
#   ungroup() %>%
#   gather_set_data(1:2) %>%
#   dplyr::mutate(
#     x = factor(x, levels = unique(x)),
#     y = factor(y, levels = unique(y))
#   )
# 
# 
# data_labels <- tibble(
#   group = c(
#     rep('manual', length(manual)),
#     rep('CIPR', length(CIPR))
#   )
# ) %>%
#   mutate(
#     hjust = ifelse(group == 'sample', 1, 0),
#     nudge_x = ifelse(group == 'sample', -0.1, 0.1)
#   )
# 
# 
# 
# p2 <- ggplot(data, aes(x, id = id, split = y, value = n)) +
#   geom_parallel_sets(aes(fill = celltype.CIPR), alpha = 0.75, axis.width = 0.15) +
#   geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
#   geom_text(
#     aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
#     hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
#   ) +
#   scale_x_discrete(labels = c('Manual','CIPR')) +
#   scale_fill_manual(values = color_assignments) +
#   theme_bw() +
#   theme(
#     legend.position = 'none',
#     axis.title = element_blank(),
#     axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
#     axis.text.y = element_blank(),
#     axis.ticks = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank()
#   )
# p2
# 
# 
# p1+p2
# 

########################### cluster - manual
color_assignments <- setNames(
  c(hue_pal()(23), pal),
  c(cluster,manual, "Unknown", "Multiplet")
)


data <- immune.combined@meta.data %>%
  group_by(integrated_snn_res.0.8,celltype.minor) %>%
  tally() %>%
  ungroup() %>%
  gather_set_data(1:2) %>%
  dplyr::mutate(
    x = factor(x, levels = unique(x)),
    y = factor(y, levels = unique(y))
  )


data_labels <- tibble(
  group = c(
    rep('cluster', length(cluster)),
    rep('manual', length(manual))
  )
) %>%
  mutate(
    hjust = ifelse(group == 'sample', 1, 0),
    nudge_x = ifelse(group == 'sample', -0.1, 0.1)
  )



p3 <- ggplot(data, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = celltype.minor), alpha = 0.75, axis.width = 0.15) +
  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
  geom_text(
    aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
    hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
  ) +
  scale_x_discrete(labels = c('Cluster','Expression-based')) +
  scale_fill_manual(values = color_assignments) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_blank(),
    axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )
p3
########################### cluster - marker
data <- immune.combined@meta.data %>%
  group_by(integrated_snn_res.0.8,celltype.marker) %>%
  tally() %>%
  ungroup() %>%
  gather_set_data(1:2) %>%
  dplyr::mutate(
    x = factor(x, levels = unique(x)),
    y = factor(y, levels = unique(y))
  )


data_labels <- tibble(
  group = c(
    rep('cluster', length(cluster)),
    rep('marker', length(marker))
  )
) %>%
  mutate(
    hjust = ifelse(group == 'sample', 1, 0),
    nudge_x = ifelse(group == 'sample', -0.1, 0.1)
  )



p4 <- ggplot(data, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = celltype.marker), alpha = 0.75, axis.width = 0.15) +
  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
  geom_text(
    aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
    hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
  ) +
  scale_x_discrete(labels = c('Cluster','Marker-based')) +
  scale_fill_manual(values = color_assignments) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_blank(),
    axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )
p4
p3+p4
########################### cluster - CIPR
data <- immune.combined@meta.data %>%
  group_by(integrated_snn_res.0.8,celltype.CIPR) %>%
  tally() %>%
  ungroup() %>%
  gather_set_data(1:2) %>%
  dplyr::mutate(
    x = factor(x, levels = unique(x)),
    y = factor(y, levels = unique(y))
  )


data_labels <- tibble(
  group = c(
    rep('cluster', length(cluster)),
    rep('CIPR', length(CIPR))
  )
) %>%
  mutate(
    hjust = ifelse(group == 'sample', 1, 0),
    nudge_x = ifelse(group == 'sample', -0.1, 0.1)
  )



p5 <- ggplot(data, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = celltype.CIPR), alpha = 0.75, axis.width = 0.15) +
  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
  geom_text(
    aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
    hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
  ) +
  scale_x_discrete(labels = c('Cluster','CIPR')) +
  scale_fill_manual(values = color_assignments) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_blank(),
    axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )
p5
p3+p4+p5


################# Final assignment ##################
immune.combined$celltype <- immune.combined$celltype.minor
immune.combined$celltype[immune.combined$integrated_snn_res.0.8==20] <- as.factor("Non-Immune")
immune.combined$celltype[immune.combined$integrated_snn_res.0.8==22] <- as.factor("Non-Immune")
levels(immune.combined$celltype) <- c(levels(immune.combined$celltype), "Multiplet")
immune.combined$celltype[immune.combined$integrated_snn_res.0.8==17 & immune.combined$celltype.minor=="NK"] <- as.factor("Multiplet")
immune.combined$celltype[immune.combined$celltype.minor=="Cd4+/Cd8+ doublet"] <- as.factor("Multiplet")
DimPlot(immune.combined, group.by = "celltype", cols = color_assignments)+ggtitle("")

immune.combined <- ScaleData(immune.combined, vars.to.regress = 'percent.mt')
saveRDS(immune.combined, "final_annotation.rds")

########################## Composition of samples by cell type
table_samples_by_cell_type <- immune.combined@meta.data %>%
  dplyr::group_by(day, celltype) %>%
  dplyr::summarize(count = n()) %>%
  tidyr::spread(celltype, count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c("day", "total_cell_count", dplyr::everything()))

table_samples_by_cell_type %>% knitr::kable()


temp_labels <- immune.combined@meta.data %>%
  group_by(day) %>%
  tally()

table_samples_by_cell_type %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'day') %>%
  mutate(sample = factor(day, levels = levels(immune.combined@meta.data$day))) %>%
  ggplot(aes(day, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(
      x = day,
      y = Inf,
      label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)),
      vjust = -1
    ),
    color = 'black', size = 2.8
  ) +
  scale_fill_manual(name = 'Cell type', values = pal) +
  scale_y_continuous(
    name = 'Percentage [%]',
    labels = scales::percent_format(),
    expand = c(0.01,0)
  ) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )

####################### Celltype heatmap
Idents(immune.combined) <- "celltype"
cellmarker <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(cellmarker, "cellmarker.csv")
library(openxlsx)
cluster_names <- unique(cellmarker$cluster)
wb <- createWorkbook()
for (cluster in cluster_names) {
  subset_markers <- cellmarker[cellmarker$cluster == cluster, ]
  sheet_name <- paste("Cluster", cluster)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, subset_markers)
}
saveWorkbook(wb, file = "cellmarker.xlsx")

cellmarker %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

DotPlot(immune.combined, features = unique(top5$gene), group.by = 'celltype')+RotatedAxis()