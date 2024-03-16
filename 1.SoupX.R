library(SoupX)
library(Seurat)
library(DropletUtils)
# 指定数据文件夹的路径
folder_path <- "F:/Xujianming_cellreport/"

# 获取所有待处理的文件夹名称 (SRR1565XXXX 格式的文件夹名)
folder_names <- list.files(folder_path, full.names = TRUE)[1]

# 循环处理每个文件夹中的数据
for (i in seq_along(folder_names)) {
  cat("\n\nProcessing folder", i, ": ", folder_names[i])
  
  # 加载原始基因表达数据
  tod_file <- file.path(folder_names[i], "raw_feature_bc_matrix")
  tod_mat <- Read10X(data.dir = tod_file)
  
  # 加载过滤后的基因表达数据
  toc_file <- file.path(folder_names[i], "filtered_feature_bc_matrix")
  toc_mat <- Read10X(data.dir = toc_file)
  
  # 运行SoupX算法和Seurat分析流程
  sc <- SoupChannel(tod_mat, toc_mat, calcSoupProfile=TRUE)
  srat <- CreateSeuratObject(sc$toc)
  srat <- ScaleData(srat)
  srat <- FindVariableFeatures(srat)
  srat <- RunPCA(srat,pcs.compute=30)
  srat <- RunUMAP(srat, dims=seq(30))
  srat <- FindNeighbors(srat, reductions.type="umap", dims = seq(30), resolution = 1)
  srat <- FindClusters(srat,reductions.type="umap", dims = seq(30), resolution = 1)
  
  # 提取聚类和互动缩放数据，然后写出处理结果
  metadata <- as.data.frame(srat[['umap']]@cell.embeddings)
  colnames(metadata) <- c('RD1','RD2')
  metadata$Cluster <- factor(srat@meta.data[rownames(metadata),'RNA_snn_res.1'])
  sc <- setClusters(sc, setNames(metadata$Cluster, rownames(metadata)))
  sc <- setDR(sc, metadata[colnames(sc$toc), c("RD1", "RD2")])
  sc <- autoEstCont(sc)
  out <- adjustCounts(sc, roundToInt = TRUE)
  
  output_dir <- file.path(paste0(folder_names[i],"/filtered_feature_bc_matrix_SoupX_out"))
  DropletUtils::write10xCounts(x=out,output_dir)
  
  # 清除中间变量以释放内存
  rm(tod_mat, toc_mat, sc, srat, metadata, out)
  gc()
}