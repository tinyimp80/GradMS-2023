library(CellChat)
library(igraph)
library(ComplexHeatmap)
DB <- CellChatDB.mouse

expr <- subset@assays$RNA@data
meta <- subset@meta.data
table(subset$day)
setwd("/home/park/project/graduate/GSE200639_tuberculosis/output/cellchat/")
# for (i in levels(factor(subset$day))[3]) {
  day <- str_replace(i, "\\/", "")
  inf <- rownames(meta)[meta$day == str_replace(i, "\\/", "")]
  inf_input = expr[, inf]
  inf_meta = meta[inf, ]
  unique(inf_meta$celltype)
  inf_meta$celltype <- factor(inf_meta$celltype)
  cellchat_inf <- createCellChat(object = inf_input, meta = inf_meta, group.by = 'celltype')
  
  cellchat_inf <- addMeta(cellchat_inf, meta = inf_meta)
  cellchat_inf <- setIdent(cellchat_inf, ident.use = 'celltype')
  
  levels(cellchat_inf@meta$celltype)
  
  groupSize_inf <- as.numeric(table(cellchat_inf@idents))
  cellchat_inf@DB <- DB
  
  
  
  
  
  # Preprocessing the expression data for cell-cell communication analysis
  cellchat_inf <- subsetData(cellchat_inf)
  cellchat_inf <- identifyOverExpressedGenes(cellchat_inf)
  cellchat_inf <- identifyOverExpressedInteractions(cellchat_inf)
  cellchat_inf <- projectData(cellchat_inf, PPI.mouse)
  
  # Compute the communication probability and infer cellular communication network
  cellchat_inf <- computeCommunProb(cellchat_inf)
  cellchat_inf <- filterCommunication(cellchat_inf, min.cells = 0)
  
  # Extract the inferred cellular communication network as a data frame
  # df.net_covid <- subsetCommunication(cellchat_inf)
  
  # Infer the cell-cell communication at a signaling pathway level
  cellchat_inf <- computeCommunProbPathway(cellchat_inf)
  
  # Calculate the aggregated cell-cell communication network
  cellchat_inf <- aggregateNet(cellchat_inf)
  
  # Compute centrality
  cellchat_inf <- netAnalysis_computeCentrality(cellchat_inf, slot.name = "netP")
  png(paste0(i, "_cellchat.png"))
  netVisual_circle(cellchat_inf@net$weight, vertex.weight = groupSize_inf, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use = pal)
  dev.off()
  cellchat_inf@netP$pathways
  # Hierarchy plot
  # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
  vertex.receiver = seq(1,4) # a numeric vector. 
  png(paste0(i, "_cellchat_mhc1.png"))
  netVisual_aggregate(cellchat_inf, signaling = "MHC-I",  vertex.receiver = vertex.receiver)
  dev.off()
  png(paste0(i, "_cellchat_mhc2.png"))
  netVisual_aggregate(cellchat_inf, signaling = "MHC-II",  vertex.receiver = vertex.receiver)
  dev.off()
  png(paste0(i, "_cellchat_ccl.png"))
  netVisual_aggregate(cellchat_inf, signaling = "CCL",  vertex.receiver = vertex.receiver)
  dev.off()
  png(paste0(i, "_cellchat_il16.png"))
  netVisual_aggregate(cellchat_inf, signaling = "IL16",  vertex.receiver = vertex.receiver)
  dev.off()
  # png(paste0(i, "_cellchat_cxcl.png"))
  # netVisual_aggregate(cellchat_inf, signaling = "CXCL",  vertex.receiver = vertex.receiver)
  # dev.off()
  png(paste0(i, "_cellchat_inf2.png"))
  netVisual_aggregate(cellchat_inf, signaling = "IFN-II",  vertex.receiver = vertex.receiver)
  dev.off()
  png(paste0(i, "_cellchat_pecam1.png"))
  netVisual_aggregate(cellchat_inf, signaling = "PECAM1",  vertex.receiver = vertex.receiver)
  dev.off()
  png(paste0(i, "_cellchat_TNF.png"))
  netVisual_aggregate(cellchat_inf, signaling = "TNF",  vertex.receiver = vertex.receiver)
  dev.off()
  
  # Heatmap
  # par(mfrow=c(1,1))
  # netVisual_heatmap(cellchat_inf, signaling = "TNF", color.heatmap = "Reds")
  # Do heatmap based on a single object
  
  
  # netVisual_bubble(cellchat_inf, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
  # netVisual_bubble(cellchat_uninf, sources.use = 4, targets.use = c(1:6), signaling = "MHC-I", remove.isolate = FALSE)
  
  
  plotGeneExpression(cellchat_inf, signaling = "MHC-I", color.use = pal)
  png(paste0(i, "_exp_mhc1.png"))
  plotGeneExpression(cellchat_inf, signaling = "MHC-I", color.use = pal)
  dev.off()
  png(paste0(i, "_exp_mhc2.png"))
  plotGeneExpression(cellchat_inf, signaling = "MHC-II", color.use = pal)
  dev.off()
  png(paste0(i, "_exp_ccl.png"))
  plotGeneExpression(cellchat_inf, signaling = "CCL",  color.use = pal)
  dev.off()
  png(paste0(i, "_exp_il16.png"))
  plotGeneExpression(cellchat_inf, signaling = "IL16", color.use = pal)
  dev.off()
  # png(paste0(i, "_exp_cxcl.png"))
  # plotGeneExpression(cellchat_inf, signaling = "CXCL", color.use = pal)
  # dev.off()
  png(paste0(i, "_exp_inf2.png"))
  plotGeneExpression(cellchat_inf, signaling = "IFN-II", color.use = pal)
  dev.off()
  png(paste0(i, "_exp_pecam1.png"))
  plotGeneExpression(cellchat_inf, signaling = "PECAM1", color.use = pal)
  dev.off()
  png(paste0(i, "_exp_TNF.png"))
  plotGeneExpression(cellchat_inf, signaling = "TNF", color.use = pal)
  dev.off()
}
