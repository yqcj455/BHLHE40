library(CellChat)
library(tidyverse)
library(ggalluvial)
library(Seurat)
library(data.table)
library(ggsci)
library(msigdbr)
library(tibble)

setwd("G:\\YS_work2/9.CellphoneDB_epi_stro_TM-TnM/")
df1 <- fread("allepi_fib_epi_all_normalized_data_TM.txt")
cell=df1$Gene
df1=as.data.frame(df1)
df1=df1[,-1] 
rownames(df1)=cell
library(Matrix)
# 然后将 dataframe 转换为稀疏矩阵
df1_sparse <- as(as.matrix(df1), "sparseMatrix")

pd1 <- data.table::fread("allepi_fib_epi_all_meta_TM.txt")
pd1=as.data.frame(pd1)
pd1=column_to_rownames(pd1,"Cell")
##提取表达矩阵和细胞分类信息
cellchat_TM <- createCellChat(object = df1_sparse, meta = pd1,  group.by = "cell_type")
#####################################################################
df2 <- fread("allepi_fib_epi_all_normalized_data_TnM.txt")
cell=df2$Gene
df2=as.data.frame(df2)
df2=df2[,-1] 
rownames(df2)=cell
library(Matrix)
# 然后将 dataframe 转换为稀疏矩阵
df2_sparse <- as(as.matrix(df2), "sparseMatrix")

pd2 <- data.table::fread("allepi_fib_epi_all_meta_TnM.txt")
pd2=as.data.frame(pd2)
pd2=column_to_rownames(pd2,"Cell")
cellchat_TnM <- createCellChat(object = df2_sparse, meta = pd2,  group.by = "cell_type")

dir.create("compare")
setwd("compare/")
####可选CellChatDB.human, CellChatDB.mouse
CellChatDB <- CellChatDB.human
##下一步不出图的时候运行 dev.new()
showDatabaseCategory(CellChatDB)

##
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor) #共刺激因子
head(CellChatDB$complex)
head(CellChatDB$geneInfo)

cellchat=cellchat_TnM 
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat <- subsetData(cellchat)
plan("multisession", workers = 16)
options(future.globals.maxSize = 64*1024^3)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.TnM = cellchat
#################################
cellchat=cellchat_TM
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat <- subsetData(cellchat)
plan("multisession", workers = 16)
options(future.globals.maxSize = 64*1024^3)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.TM = cellchat
##############################################
cc.list=list(TnM=cc.TnM,TM=cc.TM)
cellchat=mergeCellChat(cc.list,cell.prefix = T,add.names = names(cc.list))
save(cellchat,file = "YS_fib_epi_cellchat_TMvsTnM.Rdata")
##可视化
##所有细胞群总体观：通讯数量与强度对比
# 比较通讯数量和通讯强度的差异图
pdf("Interactions_diffInteraction_count.pdf")
compareInteractions(cellchat, show.legend = FALSE, group = c(1, 3), measure = "count")
dev.off()
pdf("Interactions_diffInteraction_weight.pdf")
compareInteractions(cellchat, show.legend = FALSE, group = c(1, 3), measure = "weight")
dev.off()
##第一个图展示通讯数量之间的差异，第二个图展示通讯强度之间的差异。 

# 数量与强度差异的网络图
pdf("netVisual_diffInteraction_count.pdf")
netVisual_diffInteraction(cellchat, weight.scale = TRUE)
dev.off()

pdf("netVisual_diffInteraction_weight.pdf")
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
dev.off()
##红色是case相对于control上调的，蓝色是下调的

# 数量与强度差异的热图
pdf("netVisual_heatmap_count.pdf")
netVisual_heatmap(cellchat)
dev.off()

pdf("netVisual_heatmap_weight.pdf")
netVisual_heatmap(cellchat, measure = "weight")
dev.off()
#case(TM)和control(TnM)对比，红色是上调，蓝色是下调
library(ggplot2)
# 保守和特异性信号通路的识别与可视化
pdf("rankNet_comparison_stacked_polished.pdf",height = 5,width = 10)
rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE,do.flip = F)
dev.off()

pdf("rankNet_comparison.pdf")
rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)
dev.off()
##左图最下面多个信号通路是case组独有的

# 细胞互作数量对比网络图
weight.max = getMaxWeight(cc.list, attribute = c("idents", "count"))

pdf("netVisual_circle_TnM.pdf")
netVisual_circle(cc.list[[1]]@net$count, weight.scale = TRUE, label.edge = FALSE,
                 edge.weight.max = weight.max[2], edge.width.max = 12, title.name = "TnM")
dev.off()

pdf("netVisual_circle_TM.pdf")
netVisual_circle(cc.list[[2]]@net$count, weight.scale = TRUE, label.edge = FALSE,
                 edge.weight.max = weight.max[2], edge.width.max = 12, title.name = "TM")
dev.off()

table(pd1)
s.cell=c( "Malignant EPCs8", "SFRP2+ CAFs","CXCL1+ CAFs")
count1=cc.list[[1]]@net$count[s.cell,s.cell]
count2=cc.list[[2]]@net$count[s.cell,s.cell]
pdf("netVisual_circle_TnM_select.pdf")
netVisual_circle(count1,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "TnM" )
dev.off()
pdf("netVisual_circle_TM_select.pdf")
netVisual_circle(count2,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "TM" )
dev.off()

#### 通路 ####

df=cbind(df1,df2)
df_sparse=as(as.matrix(df), "sparseMatrix")
pd=rbind(pd1,pd2)
cellchat <- createCellChat(object = df_sparse, meta = pd,  group.by = "cell_type")
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
table(CellChatDB.use[["interaction"]][["pathway_name"]])
####可选CellChatDB.human, CellChatDB.mouse
# CellChatDB <- CellChatDB.human
##下一步不出图的时候运行 dev.new()
# showDatabaseCategory(CellChatDB)
# CellChatDB.use2 <- subsetDB(CellChatDB, search = "ECM-Receptor")
# table(CellChatDB.use2[["interaction"]][["pathway_name"]])
cellchat@DB <- CellChatDB.use # set the used database in the object
#对表达数据进行预处理
##将信号基因的表达数据进行子集化，以节省计算成本
cellchat <- subsetData(cellchat)
library(future)
plan("multisession", workers = 16)
options(future.globals.maxSize = 64*1024^3)
# 识别过表达基因，计算时间较长（10线程大概5分钟）
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)
##相互作用推断
## 1、计算通信概率推断细胞互作的通信网络
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
###如果特定细胞群中只有少数细胞，则过滤掉细胞间的通信
cellchat <- filterCommunication(cellchat, min.cells = 3)

#提取推断出的细胞互作的通信网络数据框，我们提供了一个subsetCommunication 函数，
#可以方便地访问感兴趣的推断的细胞间通信。

##返回一个数据框，包含所有推断的配体/受体级别的细胞-细胞通信。设置slot.name = "netP"以访问信令路径级别的推断通信
df.net.1 <- subsetCommunication(cellchat,slot.name = "netP")
df.net.2 <- subsetCommunication(cellchat )
##
df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
df.net2 <- subsetCommunication(cellchat, signaling = c("CXCL", "TGFb","PERIOSTIN"))

##2、在信号通路水平上推断细胞间的通讯
cellchat <- computeCommunProbPathway(cellchat)
##汇总通信概率来计算细胞间的聚合通信网络。
cellchat <- aggregateNet(cellchat)
##3、计算聚合细胞互作通信网络
groupSize <- as.numeric(table(cellchat@idents))
# 创建文件夹
if (!dir.exists("cellchat_pdf")) {
  dir.create("cellchat_pdf")
}
pdf('cellchat_pdf/cellchat_fib_epi_interactions_and_weights.pdf', onefile = TRUE)
layout(matrix(c(1,2), 1, 2, byrow = TRUE)) # 创建一个两个图形区域的布局

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
#左图：外周各种颜色圆圈的大小表示细胞的数量，圈越大，细胞数越多。发出箭头的细胞表达配体，
#箭头指向的细胞表达受体。配体-受体对越多，线越粗。
#右图：互作的概率或者强度值（强度就是概率值相加）

dev.off() 


##每个细胞如何跟别的细胞互作（互作的强度或概率图）
pdf('cellchat_pdf/cell_interaction_weights.pdf',height = 18,width =18)
mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
##每个细胞如何跟别的细胞互作（number of interaction图）
pdf('cellchat_pdf/cell_interaction_count.pdf',height = 18,width =18)
mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
##可视化每个信号通路
##查看通路

levels(cellchat@idents)            #查看细胞顺序
vertex.receiver = c(1, 2)          #指定靶细胞的索引
cellchat@netP$pathways             #查看富集到的信号通路
pathways.show <- c("CXCL", "TGFb","PERIOSTIN")            #指定需要展示的通路

##层次图
pdf('cellchat_pdf/CXCL_interaction_heirarchy.pdf')
vertex.receiver = seq(1,2) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show[1],  vertex.receiver = vertex.receiver,layout="hierarchy")
dev.off()
# 在层次图中，实体圆和空心圆分别表示源和目标。圆的大小与每个细胞组的细胞数成比例。线越粗，互作信号越强。
# 左图中间的target是我们选定的靶细胞。右图是选中的靶细胞之外的另外一组放在中间看互作。
##圈图
pdf('cellchat_pdf/CXCL_interaction_circle.pdf')
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling =pathways.show[1], layout = "circle")
dev.off()
##和弦图
pdf('cellchat_pdf/CXCL_interaction_chord.pdf')
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling =pathways.show[1], layout = "chord", vertex.size = groupSize)
dev.off()
##热图
pdf('cellchat_pdf/CXCL_interaction_heatmap.pdf')
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show[1], color.heatmap = "Reds")
dev.off()
##纵轴是发出信号的细胞，横轴是接收信号的细胞，热图颜色深浅代表信号强度。
##上侧和右侧的柱子是纵轴和横轴强度的累积

##层次图
pdf('cellchat_pdf/TGFb_interaction_heirarchy.pdf')
vertex.receiver = seq(1,2) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show[2],  vertex.receiver = vertex.receiver,layout="hierarchy")
dev.off()
# 在层次图中，实体圆和空心圆分别表示源和目标。圆的大小与每个细胞组的细胞数成比例。线越粗，互作信号越强。
# 左图中间的target是我们选定的靶细胞。右图是选中的靶细胞之外的另外一组放在中间看互作。
##圈图
pdf('cellchat_pdf/TGFb_interaction_circle.pdf')
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling =pathways.show[2], layout = "circle")
dev.off()
##和弦图
pdf('cellchat_pdf/TGFb_interaction_chord.pdf')
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling =pathways.show[2], layout = "chord", vertex.size = groupSize)
dev.off()
##热图
pdf('cellchat_pdf/TGFb_interaction_heatmap.pdf')
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show[2], color.heatmap = "Reds")
dev.off()

##层次图
pdf('cellchat_pdf/PST_interaction_heirarchy.pdf')
vertex.receiver = seq(1,2) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show[3],  vertex.receiver = vertex.receiver,layout="hierarchy")
dev.off()
# 在层次图中，实体圆和空心圆分别表示源和目标。圆的大小与每个细胞组的细胞数成比例。线越粗，互作信号越强。
# 左图中间的target是我们选定的靶细胞。右图是选中的靶细胞之外的另外一组放在中间看互作。
##圈图
pdf('cellchat_pdf/PST_interaction_circle.pdf')
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling =pathways.show[3], layout = "circle")
dev.off()
##和弦图
pdf('cellchat_pdf/PST_interaction_chord.pdf')
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling =pathways.show[3], layout = "chord", vertex.size = groupSize)
dev.off()
##热图
pdf('cellchat_pdf/PST_interaction_heatmap.pdf')
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show[3], color.heatmap = "Reds")
dev.off()

