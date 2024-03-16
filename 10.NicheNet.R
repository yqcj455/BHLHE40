#### Mal8 ####
# install.packages("devtools")
#devtools::install_github("saeyslab/nichenetr")
library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)
options(timeout=600)
organism = "human"

lr_network = readRDS("lr_network_human_21122021.rds")
ligand_target_matrix = readRDS("ligand_target_matrix_nsga2r_final.rds")
weighted_networks = readRDS("weighted_networks_nsga2r_final.rds")


library(SeuratDisk)
Convert("YS_adata_fib_epi_after_newcelltype.h5ad", dest = "h5seurat", overwrite = F)
seuratObj=LoadH5Seurat("YS_adata_fib_epi_after_newcelltype.h5seurat",meta.data = F)
meta=read.table("adata_fib_epi_subcelltype_obs_infor.csv",header = T,sep = ",")

table(meta$Index == colnames(seuratObj))
# 将每一列添加到sc_new对象的元数据中
for (col in colnames(meta)) {
  seuratObj@meta.data[[col]] <- meta[[col]]
}
save(seuratObj,file = "Malepi_fib_after_newcelltype.Rdata")
library(dplyr)
seuratObj@meta.data$Fib_Epi_Cluster %>% table() # note that the number of cells of some cell types is very low and should preferably be higher for a real application
## .
# CCL11+ CAFs      CXCL1+ CAFs     CXCL14+ CAFs  Malignant EPCs1 Malignant EPCs10 
# 171              629              484             4267              700 
# Malignant EPCs11  Malignant EPCs2  Malignant EPCs3  Malignant EPCs4  Malignant EPCs5 
# 693             4206             4030             3988             2622 
# Malignant EPCs6  Malignant EPCs7  Malignant EPCs8  Malignant EPCs9      MYH11+ CAFs 
# 1586             1384             1252              915              373 
# RGS5+ CAFs      SFRP2+ CAFs 
# 681             1564 
DimPlot(seuratObj, reduction = "umap",group.by ="Fib_Epi_Cluster")


Idents(seuratObj)=seuratObj@meta.data[["Fib_Epi_Cluster"]]
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  receiver = "Malignant EPCs8", 
  condition_colname = "Metastases", condition_oi = "With Metastases", condition_reference = "Without Metastases", 
  sender = c("SFRP2+ CAFs","CXCL1+ CAFs"), 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)
# save(nichenet_output,file = "nichenet_output_Mal8.Rdata")
# 主要根据pearson列进行排名，该值越高表明ligand的target基因富集到receiver细胞的差异基因集中。
# 如手册中提到，若该值高于0.1，则可认为是显著的富集
nichenet_output$ligand_activities
ligand_activities=as.data.frame(nichenet_output$ligand_activities)
head(ligand_activities)
library(ggplot2)

# 选择前30个观测值并按照 'aupr' 排序
ligand_activities_subset <- ligand_activities %>%
  mutate(test_ligand = reorder(test_ligand, aupr)) %>%
  head(30)
pdf("nichenet_ligand_AUPR.pdf",height = 8,width = 3)
# 创建热图
min_val <- min(ligand_activities_subset$aupr)
max_val <- max(ligand_activities_subset$aupr)
mid_val <- (min_val + max_val) / 2

# 创建热图
ggplot(ligand_activities_subset, aes(x = factor(1), y = test_ligand, fill = aupr)) +
  geom_tile(color = "white") + 
  scale_fill_gradient(low = "white", high = "orange", 
                      breaks = c(min_val, mid_val, max_val),  # 设置图例的断点
                      labels = format(c(min_val, mid_val, max_val), digits = 2)) +  # 设置图例的标签
  labs(x = "", y = "Prioritized ligands", fill = "AUPR") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom")  # 将图例放在底部
dev.off()
nichenet_output$top_ligands
# [1] "ENG"      "PATJ"     "TIMP1"    "ICAM1"    "TGFB1"    "CFD"      "COL5A3"   "JAM2"    
# [9] "COL6A3"   "CDCP1"    "EDIL3"    "EFNA1"    "OCLN"     "TRAF2"    "ANXA1"    "COL18A1" 
# [17] "WNT5A"    "COL5A2"   "A2M"      "PVR"      "SERPING1" "EREG"     "COL3A1"   "MMP9"    
# [25] "JAM3"     "MMP1"     "LGALS3"   "EDN1"     "RARRES1"  "LAMB2"   

## 相关可视化
# ligand在不同细胞类型的表达
pdf("nichenet_dotplot.pdf")
nichenet_output$ligand_expression_dotplot
dev.off()
# ligand在TMvsTnM前后的差异表达
nichenet_output$ligand_differential_expression_heatmap

# ligand-target score
nichenet_output$ligand_target_matrix %>% .[1:10,1:6]
nichenet_output$ligand_target_df %>% head
nichenet_output$top_targets
TF=as.data.frame(nichenet_output[["ligand_target_heatmap"]][["data"]])
TGFB1=TF[TF$y=="TGFB1",]
Ligand=arrange(TGFB1,desc(score))
Top_Ld=Ligand$x[1:20]
Top_Ld=as.character(Top_Ld)
write.table(Top_Ld,"Top_ligand_target.txt",quote = F,sep = "\t",row.names = F)
## 可视化
# ligand-target score热图

nichenet_output$ligand_target_heatma

pdf("nichenet_heatmap_mal8_polished.pdf",height = 8,width = 12)
nichenet_output$ligand_target_heatmap + 
  scale_fill_gradient2(low = "whitesmoke",  high = "royalblue") + 
  xlab("Predicted target genes in Malignant EPCs8") + 
  ylab("Prioritized SFRP2+ CAFs and CXCL1+ CAFs ligands") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8))
dev.off()
# 靶基因在reciver细胞类型中的差异表达
pdf("nichenet_target_dotplot_mal8.pdf",height = 5,width = 10)
DotPlot(seuratObj %>% subset(idents = "Malignant EPCs8"), 
        features = Top_Ld, cols = "RdYlBu",
        split.by = "Metastases") + 
  RotatedAxis()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
#找差异表达基因
Mal8_Seurat=seuratObj %>% subset(idents = "Malignant EPCs8")

library(Seurat)
library(cowplot)

# 使用 VlnPlot 生成一个图形列表
plot_list <- lapply(1:20, function(i) {
  VlnPlot(seuratObj %>% subset(idents = "Malignant EPCs8"),
          features = Top_Ld[i],
          split.by = "Metastases",
          pt.size = 0)
})

# 将20个图形组合成一个5*4的布局
combined_plot <- plot_grid(plotlist = plot_list, ncol = 5, nrow = 4)

# 开始 PDF 设备输出,一定要保证尺寸够大
pdf("nichenet_target_vlnplot.pdf",width = 20,height = 16)
# 打印组合后的图形
print(combined_plot)
# 关闭 PDF 设备
dev.off()
pdf("nichenet_target_vlnplot_BHLHE40_Mal8.pdf",width = 6,height = 4)
VlnPlot(seuratObj %>% subset(idents = "Malignant EPCs8"),
        features = 'BHLHE40',
        split.by = "Metastases",
        pt.size = 0)+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
dev.off()

pdf("nichenet_target_vlnplot_BHLHE40_Mal8.pdf",width = 6,height = 4)
VlnPlot(seuratObj %>% subset(idents = "Malignant EPCs8"),
        features = 'FOS',
        split.by = "Metastases",
        pt.size = 0)+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
dev.off()

FeaturePlot(object = seuratObj, features = "FOS")
FeaturePlot(object = seuratObj, features = "BHLHE40")


pdf("nichenet_target_vlnplot_TGFBR1_Mal8.pdf",width = 6,height = 4)
VlnPlot(seuratObj %>% subset(idents = "Malignant EPCs8"),
        features = 'TGFBR1',
        split.by = "Metastases",
        pt.size = 0)+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
dev.off()

pdf("nichenet_target_vlnplot_TGFBR2_Mal8.pdf",width = 6,height = 4)
VlnPlot(seuratObj %>% subset(idents = "Malignant EPCs8"),
        features = 'TGFBR2',
        split.by = "Metastases",
        pt.size = 0)+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
dev.off()

pdf("nichenet_target_vlnplot_TGFBR3_Mal8.pdf",width = 6,height = 4)
VlnPlot(seuratObj %>% subset(idents = "Malignant EPCs8"),
        features = 'TGFBR3',
        split.by = "Metastases",
        pt.size = 0)+ theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
dev.off()



## 上述的综合可视化
pdf("nichenet_activity_target_heatmap.pdf",width = 40,height =8)
nichenet_output$ligand_activity_target_heatmap
dev.off()
pdf("nichenet_receptor_heatmap_Mal8.pdf",height = 8,width = 9)
nichenet_output$ligand_receptor_heatmap
dev.off()

#### CAF ####
Idents(seuratObj)=seuratObj@meta.data[["Fib_Epi_Cluster"]]
nichenet_output2 = nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  receiver = "CXCL1+ CAFs", 
  condition_colname = "Metastases", condition_oi = "With Metastases", condition_reference = "Without Metastases", 
  sender = c("SFRP2+ CAFs","Malignant EPCs8"), 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)
# save(nichenet_output2,file = "nichenet_output_CXCL1.Rdata")
# 主要根据pearson列进行排名，该值越高表明ligand的target基因富集到receiver细胞的差异基因集中。
# 如手册中提到，若该值高于0.1，则可认为是显著的富集
nichenet_output2$ligand_activities
ligand_activities=as.data.frame(nichenet_output2$ligand_activities)
head(ligand_activities)

library(ggplot2)

# 选择前30个观测值并按照 'aupr' 排序
ligand_activities_subset <- ligand_activities %>%
  mutate(test_ligand = reorder(test_ligand, aupr)) %>%
  head(30)
pdf("nichenet_ligand_AUPR_CXCL1.pdf",height = 8,width = 3)
# 创建热图
min_val <- min(ligand_activities_subset$aupr)
max_val <- max(ligand_activities_subset$aupr)
mid_val <- (min_val + max_val) / 2

# 创建热图
ggplot(ligand_activities_subset, aes(x = factor(1), y = test_ligand, fill = aupr)) +
  geom_tile(color = "white") + 
  scale_fill_gradient(low = "white", high = "orange", 
                      breaks = c(min_val, mid_val, max_val),  # 设置图例的断点
                      labels = format(c(min_val, mid_val, max_val), digits = 2)) +  # 设置图例的标签
  labs(x = "", y = "Prioritized ligands", fill = "AUPR") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom")  # 将图例放在底部
dev.off()
nichenet_output2$top_ligands
# [1] "ENG"      "PATJ"     "TIMP1"    "ICAM1"    "TGFB1"    "CFD"      "COL5A3"   "JAM2"    
# [9] "COL6A3"   "CDCP1"    "EDIL3"    "EFNA1"    "OCLN"     "TRAF2"    "ANXA1"    "COL18A1" 
# [17] "WNT5A"    "COL5A2"   "A2M"      "PVR"      "SERPING1" "EREG"     "COL3A1"   "MMP9"    
# [25] "JAM3"     "MMP1"     "LGALS3"   "EDN1"     "RARRES1"  "LAMB2"   

## 相关可视化
# ligand在不同细胞类型的表达
pdf("nichenet_dotplot2.pdf")
nichenet_output2$ligand_expression_dotplot
dev.off()
# ligand在TMvsTnM前后的差异表达
nichenet_output2$ligand_differential_expression_heatmap

# ligand-target score
nichenet_output2$ligand_target_matrix %>% .[1:10,1:6]
nichenet_output2$ligand_target_df %>% head
nichenet_output2$top_targets
TF=as.data.frame(nichenet_output2[["ligand_target_heatmap"]][["data"]])
TGFB1=TF[TF$y=="TGFB1",]
Ligand=arrange(TGFB1,desc(score))
Top_Ld=Ligand$x[1:20]

## 可视化
# ligand-target score热图

nichenet_output2$ligand_target_heatma

pdf("nichenet_heatmap_CXCL1_polished.pdf",height = 8,width = 12)
nichenet_output2$ligand_target_heatmap + 
  scale_fill_gradient2(low = "whitesmoke",  high = "royalblue") + 
  xlab("Predicted target genes in CXCL1+ CAFs ligands") + 
  ylab("Prioritized SFRP2+ CAFs and Malignant EPCs8")
dev.off()
# 靶基因在reciver细胞类型中的差异表达
pdf("nichenet_target_dotplot_CXCL1.pdf",height = 5,width = 10)
DotPlot(seuratObj %>% subset(idents = "CXCL1+ CAFs"), 
        features = Top_Ld, cols = "RdYlBu",
        split.by = "Metastases") + 
  RotatedAxis()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("nichenet_target_vlnplot_CXCL1.pdf",height = 8,width = 16)
VlnPlot(seuratObj %>% subset(idents = "CXCL1+ CAFs"), 
        features = Top_Ld, 
        split.by = "Metastases",    
        pt.size = 0)
dev.off()

pdf("nichenet_target_vlnplot_ARID5B_CXCl1.pdf")
VlnPlot(seuratObj %>% subset(idents = "CXCL1+ CAFs"),
        features = 'ARID5B',
        split.by = "Metastases",
        pt.size = 0)
dev.off()

pdf("nichenet_target_vlnplot_TGFBR1.pdf",height = 8,width = 16)
VlnPlot(seuratObj %>% subset(idents = "CXCL1+ CAFs"), 
        features = "TGFBR1", 
        split.by = "Metastases",    
        pt.size = 0)
dev.off()

## 上述的综合可视化
nichenet_output2$ligand_activity_target_heatmap
pdf("nichenet_receptor_heatmap_CXCL1.pdf",height = 8,width = 9)
nichenet_output2$ligand_receptor_heatmap
dev.off()
