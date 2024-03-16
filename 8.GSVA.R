#### GSVA stro ####
##BiocManager::install("GSVA")
library('GSEABase')
library(GSVA)
library(msigdbr)
library(data.table)

df <- fread("adata_fib_epi_subcelltype_normalized_data.csv")
df=as.data.frame(df)
rownames(df)=df$index
df=df[,-1]
df=as.data.frame(t(df))
pd <- data.table::fread("adata_fib_epi_subcelltype_obs_infor.csv")

m_df<- msigdbr(species = "human",  category = "H" ) #H代表Hallmark
geneSets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

library(clusterProfiler)
geneSet_local=read.gmt("h.all.v2023.1.Hs.symbols.gmt")
geneSet_local=split(geneSet_local$gene,geneSet_local$term)

table(pd$Fib_Epi_Cluster)
df2=df[,pd$Fib_Epi_Cluster=="Malignant EPCs8"]

GSVA_hall_mal8 <- gsva(expr=as.matrix(df2), #必须是一个表达矩阵，不能是seurat对象
                  gset.idx.list=geneSets, #必须是一个列表
                  mx.diff=T, # 数据为正态分布则T，双峰则F
                  kcdf="Gaussian", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson",必须是正整数
                  parallel.sz=4) # 并行线程数目

#BiocManager::install('limma')
library(limma)
pd2=pd[pd$Fib_Epi_Cluster=="Malignant EPCs8",]
# 设置或导入分组
pd2$meta=ifelse(pd2$Metastases=='Without Metastases','TnM','TM')
group <- factor(pd2$meta, levels = c( 'TnM','TM'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(GSVA_hall_mal8)

compare <- makeContrasts(TM - TnM, levels=design)
fit <- lmFit(GSVA_hall_mal8, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
head(Diff)

## 发散条形图绘制(还可以作热图)
## barplot
library(dplyr)
dat_plot <- data.frame(id = row.names(Diff),
                       t = Diff$t)
# 去掉"HALLMARK_"
library(stringr)
dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","") #删除Hallmark
# 新增一列 根据t阈值分类，为什么GSVA输入的是logtpm，那用logfc就没意义，差异小，所以用t
dat_plot$threshold = factor(ifelse(dat_plot$t  >-1, ifelse(dat_plot$t >= 1 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
# 排序，t值和p值含义差不多，但t值不仅可以反映显著性(越大越显著)，还可以根据大小确定上下调
dat_plot <- dat_plot %>% arrange(t)
# 变成因子类型
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
# 绘制
library(ggplot2)

##install.packages("ggthemes")
library(ggthemes)
# install.packages("ggprism")
library(ggprism)
p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#36638a','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
  geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('t value of GSVA score, TM VS TnM') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p
# 添加标签

# 小于-1的数量
low1 <- dat_plot %>% filter(t < -1) %>% nrow()
# 小于0总数量
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
# 小于1总数量
high0 <- dat_plot %>% filter(t < 1) %>% nrow()
# 总的柱子数量
high1 <- nrow(dat_plot)

# 依次从下到上添加标签
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                   hjust = 0,color = 'black') + # 小于-1的为黑色标签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 大于1的为黑色标签

p

ggsave("Mal8_gsva_bar_noscaled.pdf",p,width = 14,height  = 12)

#### CXCL1 CAFs ####
table(pd$Fib_Epi_Cluster)
df2=df[,pd$Fib_Epi_Cluster=="CXCL1+ CAFs"]

GSVA_hall_cxcl1 <- gsva(expr=as.matrix(df2), #必须是一个表达矩阵，不能是seurat对象
                       gset.idx.list=geneSets, #必须是一个列表
                       mx.diff=T, # 数据为正态分布则T，双峰则F
                       kcdf="Gaussian", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson",必须是正整数
                       parallel.sz=4) # 并行线程数目

#BiocManager::install('limma')
library(limma)
pd2=pd[pd$Fib_Epi_Cluster=="CXCL1+ CAFs",]
# 设置或导入分组
pd2$meta=ifelse(pd2$Metastases=='Without Metastases','TnM','TM')
group <- factor(pd2$meta, levels = c( 'TnM','TM'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(GSVA_hall_cxcl1)

compare <- makeContrasts(TM - TnM, levels=design)
fit <- lmFit(GSVA_hall_cxcl1, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
head(Diff)

## 发散条形图绘制(还可以作热图)
## barplot
library(dplyr)
dat_plot <- data.frame(id = row.names(Diff),
                       t = Diff$t)
# 去掉"HALLMARK_"
library(stringr)
dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","") #删除Hallmark
# 新增一列 根据t阈值分类，为什么GSVA输入的是logtpm，那用logfc就没意义，差异小，所以用t
dat_plot$threshold = factor(ifelse(dat_plot$t  >-1, ifelse(dat_plot$t >= 1 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
# 排序，t值和p值含义差不多，但t值不仅可以反映显著性(越大越显著)，还可以根据大小确定上下调
dat_plot <- dat_plot %>% arrange(t)
# 变成因子类型
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
# 绘制
library(ggplot2)

##install.packages("ggthemes")
library(ggthemes)
# install.packages("ggprism")
library(ggprism)
p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#36638a','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
  geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('t value of GSVA score, TM VS TnM') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p
# 添加标签

# 小于-1的数量
low1 <- dat_plot %>% filter(t < -1) %>% nrow()
# 小于0总数量
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
# 小于1总数量
high0 <- dat_plot %>% filter(t < 1) %>% nrow()
# 总的柱子数量
high1 <- nrow(dat_plot)

# 依次从下到上添加标签
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                   hjust = 0,color = 'black') + # 小于-1的为黑色标签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 大于1的为黑色标签

p

ggsave("CXCL1_CAF_gsva_bar_noscaled.pdf",p,width = 14,height  = 12)

#### SFRP2 CAFs ####
table(pd$Fib_Epi_Cluster)
df2=df[,pd$Fib_Epi_Cluster=="SFRP2+ CAFs"]

GSVA_hall_SFRP2 <- gsva(expr=as.matrix(df2), #必须是一个表达矩阵，不能是seurat对象
                        gset.idx.list=geneSets, #必须是一个列表
                        mx.diff=T, # 数据为正态分布则T，双峰则F
                        kcdf="Gaussian", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson",必须是正整数
                        parallel.sz=4) # 并行线程数目

#BiocManager::install('limma')
library(limma)
pd2=pd[pd$Fib_Epi_Cluster=="SFRP2+ CAFs",]
# 设置或导入分组
pd2$meta=ifelse(pd2$Metastases=='Without Metastases','TnM','TM')
group <- factor(pd2$meta, levels = c( 'TnM','TM'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(GSVA_hall_SFRP2)

compare <- makeContrasts(TM - TnM, levels=design)
fit <- lmFit(GSVA_hall_SFRP2, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
head(Diff)

## 发散条形图绘制(还可以作热图)
## barplot
library(dplyr)
dat_plot <- data.frame(id = row.names(Diff),
                       t = Diff$t)
# 去掉"HALLMARK_"
library(stringr)
dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","") #删除Hallmark
# 新增一列 根据t阈值分类，为什么GSVA输入的是logtpm，那用logfc就没意义，差异小，所以用t
dat_plot$threshold = factor(ifelse(dat_plot$t  >-1, ifelse(dat_plot$t >= 1 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
# 排序，t值和p值含义差不多，但t值不仅可以反映显著性(越大越显著)，还可以根据大小确定上下调
dat_plot <- dat_plot %>% arrange(t)
# 变成因子类型
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
# 绘制
library(ggplot2)

##install.packages("ggthemes")
library(ggthemes)
# install.packages("ggprism")
library(ggprism)
p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#36638a','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
  geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('t value of GSVA score, TM VS TnM') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p
# 添加标签

# 小于-1的数量
low1 <- dat_plot %>% filter(t < -1) %>% nrow()
# 小于0总数量
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
# 小于1总数量
high0 <- dat_plot %>% filter(t < 1) %>% nrow()
# 总的柱子数量
high1 <- nrow(dat_plot)

# 依次从下到上添加标签
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                   hjust = 0,color = 'black') + # 小于-1的为黑色标签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 大于1的为黑色标签

p

ggsave("SFRP2_CAF_gsva_bar_noscaled.pdf",p,width = 14,height  = 12)

