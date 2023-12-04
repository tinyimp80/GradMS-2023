test <- readRDS("integrated_datasets.rds")
pal <- c("#E5F5F9", "#1D91C0", "#67001F", "#CB181D", "#78C679","#A6CEE3", "#FD8D3C", "#A6D854",
"#D4B9DA", "#6A51A3", "#D9D9D9", "#FFF7BC", "#000000" ,"#F0F0F0", "#C7EAE5", "#003C30",
"#8C6BB1", "#C7E9B4", "#762A83",  "#AE017E", "#DF65B0", "#EF3B2C", "#74C476")


test <- RunTSNE(test, dims = 1:17)
DefaultAssay(test) <- "RNA"
DimPlot(test, reduction = "tsne", label = T, cols = pal)
test.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(test.markers, "test.markers.csv")
library(openxlsx)
# 결과 데이터프레임을 클러스터별로 나눕니다
cluster_names <- unique(test.markers$cluster)

# 엑셀 파일 생성
wb <- createWorkbook()

# 각 클러스터를 시트로 추가
for (cluster in cluster_names) {
  subset_markers <- test.markers[test.markers$cluster == cluster, ]
  sheet_name <- paste("Cluster", cluster)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, subset_markers)
}

# 엑셀 파일 저장
saveWorkbook(wb, file = "cluster_markers.1.4.xlsx")
