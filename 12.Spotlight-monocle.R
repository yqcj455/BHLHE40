#BiocManager::install("SPOTlight")
#devtools::install_github("https://github.com/MarcElosua/SPOTlight")
## 下面的包必须安装成功。Seurat 4.0.2会报错，建议用4.0.1，4.1.1.4.1.3
## 有些时候需要重装SingleCellExperiment才可成功
library(SPOTlight)
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
#BiocManager::install('SpatialExperiment')
#library(SpatialExperiment)
library(scater)
library(scran)
#BiocManager::install("scran")

## 读取单细胞----------------------
load("Mal8_fib_sub_after_newcelltype.Rdata")
#scRNA=subset(scRNA,tissue_type =='Tumor')  ## 一般要运行这行，此处细胞少不运行
#BiocManager::install('scuttle')

# 注意，这步报错，提示需要更新一下seurat包
sce <- as.SingleCellExperiment(Fib_epi)
sce <- logNormCounts(sce)
dec <- modelGeneVar(sce)
# 计算高变基因
hvg <- getTopHVGs(dec, n = 3000)
# 加上细胞注释信息，请注意修改！！！！！！
colLabels(sce) <- colData(sce)$Fib_Epi_Cluster
# 去掉核糖体和线粒体基因
genes <- !grepl(pattern = "^RP[L|S]|MT", x = rownames(sce))
# 计算marker基因
mgs <- scoreMarkers(sce, subset.row = genes)
# 保留最相关的marker基因
mgs_fil <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  # Filter and keep relevant marker genes, those with AUC > 0.8
  x <- x[x$mean.AUC > 0.7, ] # 这里要根据自己的数据来修改AUC值
  # Sort the genes from highest to lowest weight
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  # Add gene and cluster id to the dataframe
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)

## 读取空转数据--------------------
stRNA=readRDS('./stRNA_colon1.RDS')

library(Seurat)
SpatialDimPlot(stRNA)

# 核心过程 反卷积的第二种方法----------------------
#运行时间较长，记得保存
res <- SPOTlight(
  x = sce,
  y = stRNA,
  groups = sce$Fib_Epi_Cluster,  ## 注意修改！！！！！！
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",  ## 不变
  gene_id = "gene")

save(res,file = "Fib_Epi_res_spotlight_colon1.Rdata")
load("Fib_Epi_res_spotlight_colon1.Rdata")
## 提取结果-----
head(mat <- res$mat)[, seq_len(3)]
mod <- res$NMF


## 细胞相关性图（相关性高的可能的互作价值高），但没有太大分析价值----
plotTopicProfiles(
  x = mod,
  y = sce$Fib_Epi_Cluster,  ##注意修改！！！！
  facet = FALSE,
  min_prop = 0.01,
  ncol = 1) +
  theme(aspect.ratio = 1)
plotCorrelationMatrix(mat)
plotInteractions(mat, "heatmap")
plotInteractions(mat, "network")


## 可视化方案------
stRNA[["SPOTlight"]] <- CreateAssayObject(t(res$mat))
DefaultAssay(stRNA) <- "SPOTlight"
celltypes = rownames(stRNA)
p <- SpatialFeaturePlot(stRNA, features = celltypes[1:3], pt.size.factor = 1.6, ncol = 4, crop = TRUE)
pdf("fib_epi_spotlight_SpatialFeaturePlot.pdf",height = 5,width = 12)
print(p)
dev.off()

DefaultAssay(stRNA) <- 'Spatial' ## 记得回来


### 判别细胞类型,依据最大比例
light=as.data.frame(mat)
light_celltype=c()
for (i in 1:nrow(light)) {
  a=colnames(light)[which(light[i,]==max(light[i,]))]
  light_celltype=c(light_celltype,a)
}

table(light_celltype)

anno=data.frame(mat,check.names = F)
anno$celltype=light_celltype
# 为了stlearn
write.csv(anno,file ='spotlight_decon_mat_for_fib_epi.csv')

stRNA[["celltype"]] <- light_celltype


### 空间Monocle拟时序分析------------------------
pdf("fib_epi_spotlight_Spatial.pdf",height = 5,width = 5)
SpatialDimPlot(stRNA, group.by = "seurat_clusters", label = TRUE)
dev.off()
# 取子集，选择的区域是目的基因高表达区域
stRNAsub=subset(stRNA,seurat_clusters %in% c(0,6,7))
pdf("fib_epi_spotlight_Spatial_sub.pdf",height = 5,width = 5)
SpatialDimPlot(stRNAsub, group.by = "celltype", label = F)
dev.off()

anno2=anno[colnames(data),]
# 为了stlearn
write.csv(anno2,file ='spotlight_decon_mat_for_fib_epi_sub.csv')
## 拟时序分析
library(Seurat)

#没有monocle要先安装 BiocManager::install()
#BiocManager::install('monocle',update = F,ask = F)

#没有的包先安装
library(BiocGenerics)
library(monocle)
library(tidyverse)
library(patchwork)

data=as.matrix(stRNAsub@assays$SCT@counts)

data <- as(data, 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = stRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
## 以下代码一律不得修改
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())

mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=8, relative_expr = TRUE)

##使用monocle选择的高变基因，不修改
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
plot_ordering_genes(mycds)

#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序，报错请用4.1以上的R,请勿用4.2以上，并重装monocle
mycds <- orderCells(mycds)


#State轨迹分布图
plot4 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
plot4
dev.off()
plot1 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
plot1
plot2 <- plot_cell_trajectory(mycds, color_by = "celltype")
plot2

##合并出图
pdf("fib_epi_spotlight_Spatial_sub_monocle.pdf",height = 5,width = 10)
plotc <- plot2|plot4
plotc
dev.off()

library(ggsci)
pdf("fib_epi_sub_monocle_BHLHE40.pdf",height = 5,width = 5)
plot_genes_in_pseudotime(mycds['BHLHE40',],color_by = "celltype")+ scale_color_nejm()
dev.off()

pdf("fib_epi_sub_monocle_RBMS1.pdf",height = 5,width = 5)
plot_genes_in_pseudotime(mycds['RBMS1',],color_by = "celltype")+ scale_color_nejm()
dev.off()