#### GSEA TCGA Nb#### 
#占用内存太大，所以最先跑
library(GSEABase)
library(limma) 
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
# #devtools::install_github("junjunlab/GseaVis")
library(GseaVis)
library(data.table)
gene="BHLHE40"                                                           #确定基因名

load("COAD-mRNA-FPKM-onlyT-480.Rdata")
group <- ifelse(exp[gene,]<=median(as.numeric(exp[gene,])),"l","h")
#差异分析：
group_list <- factor(group,levels = c('l','h'))
design <- model.matrix(~factor( group_list ))
fit=lmFit(exp,design)
fit=eBayes(fit)
options(digits = 4)
topTable(fit,coef=2,adjust='BH')
deg=topTable(fit,coef=2,adjust='BH',number = Inf)
head(deg)

geneset1 <- read.gmt("h.all.v2023.1.Hs.symbols.gmt")
geneList <- deg$logFC
names(geneList) <- toupper(rownames(deg))
geneList <- sort(geneList,decreasing = T)
gsea_results <- GSEA(geneList = geneList,TERM2GENE = geneset1,verbose = F,eps=1e-10)
hall=as.data.frame(gsea_results)
write.table(hall,paste0(gene,"-GSEA-hallmark-TCGA.xlsx"),quote = F,sep="\t",row.names = F)

geneset2 <- read.gmt("c2.all.v2023.1.Hs.symbols.gmt")
gsea_results2 <- GSEA(geneList = geneList,TERM2GENE = geneset2,verbose = F,eps=1e-10)
c2all=as.data.frame(gsea_results2)
write.table(c2all,paste0(gene,"-GSEA-c2-all-TCGA.xlsx"),quote = F,sep="\t",row.names = F)

geneset3 <- read.gmt("c5.all.v7.5.1.symbols.gmt")
gsea_results3 <- GSEA(geneList = geneList,TERM2GENE = geneset3,verbose = F,eps=1e-10)
c5all=as.data.frame(gsea_results3)
write.table(c5all,paste0(gene,"-GSEA-c5-all-TCGA.xlsx"),quote = F,sep="\t",row.names = F)

####GSEA 39582 #####
load("GSE39582exp-pd.Rdata")  # Replace with your actual data path
#exp=log2(data+1) #GSE39582无需标准化
exp=deg
group <- ifelse(exp[gene,]<=median(as.numeric(exp[gene,])),"l","h")
#差异分析：
group_list <- factor(group,levels = c('l','h'))
design <- model.matrix(~factor( group_list ))
fit=lmFit(exp,design)
fit=eBayes(fit)
options(digits = 4)
topTable(fit,coef=2,adjust='BH')
deg=topTable(fit,coef=2,adjust='BH',number = Inf)
deg2=cbind(Gene=rownames(deg),deg)

geneset1 <- read.gmt("h.all.v2023.1.Hs.symbols.gmt")
geneList <- deg$logFC
names(geneList) <- toupper(rownames(deg))
geneList <- sort(geneList,decreasing = T)
gsea_results4 <- GSEA(geneList = geneList,TERM2GENE = geneset1,verbose = F,eps=1e-10)
hall2=as.data.frame(gsea_results4)
write.table(hall,paste0(gene,"-GSEA-hallmark-39582.xlsx"),quote = F,sep="\t",row.names = F)


geneset2 <- read.gmt("c2.all.v2023.1.Hs.symbols.gmt")
gsea_results5 <- GSEA(geneList = geneList,TERM2GENE = geneset2,verbose = F,eps=1e-10)
c2all2=as.data.frame(gsea_results5)
write.table(c2all,paste0(gene,"-GSEA-c2-all-39582.xlsx"),quote = F,sep="\t",row.names = F)

geneset3 <- read.gmt("c5.all.v7.5.1.symbols.gmt")
gsea_results6 <- GSEA(geneList = geneList,TERM2GENE = geneset3,verbose = F,eps=1e-10)
c5all2=as.data.frame(gsea_results6)
write.table(c5all,paste0(gene,"-GSEA-c5-all-39582.xlsx"),quote = F,sep="\t",row.names = F)

#### 生存分析 ####
#COAD
load("N_T_fpkm_COAD_sym.Rdata")
meta=read.table("COADtime.txt",header = T)
colnames(meta)=c("ID","time","event")
library(stringr)
group_list = ifelse(as.numeric(str_sub(colnames(mRNA_exp),14,15)) < 10,'tumor','normal')
group_list = factor(group_list,levels = c("normal","tumor"))
exp=mRNA_exp[,group_list=="tumor"]
meta$time=meta$time/365
exp=t(exp)
rownames(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(exp))
library(limma)
exp=avereps(exp)
rownames(meta)=meta$ID
sameid=intersect(rownames(exp),rownames(meta))
exp=exp[sameid,]
meta=meta[sameid,]
library(dplyr)
clinicalExp=cbind(meta,exp[,gene])
colnames(clinicalExp)[4]=gene

library(survival)
library(survminer)
rt=clinicalExp
outTab=data.frame()
a=rt[,gene]<=median(rt[,gene])
if (TRUE & FALSE %in% a){
  diff=survdiff(Surv(time,event) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outTab=rbind(outTab,cbind(gene=gene,pvalue=pValue))
  fit <- survfit(Surv(time, event) ~ a, data = rt)
  summary(fit)
}
pValue=round(pValue,3)
pValue=paste0("=",pValue)

pdf(file=paste(gene,"-COAD.survival.pdf",sep=""),
    width=6,
    height=6)
p=ggsurvplot(fit,
             data=rt,
             risk.table = T,
             pval = T,
             conf.int = F,
             ggtheme = theme_minimal(),
             palette = c("#FF0000","#2200FF"),
             legend.title=gene,
             legend.lab=c("High","Low"),
             xlab="Time(year)",
             ylab="Overall Survival",
             xlim=c(0,10),
             break.time.by=2)
print(p,newpage = F)
dev.off()

library(survival)
library(survminer)
rt=clinicalExp
res.cut <- surv_cutpoint(rt, #数据集
                         time = "time", #生存状态
                         event = "event", #生存时间
                         variables = colnames(rt)[4:ncol(rt)] #需要计算的数据列名
)

res.cut=summary(res.cut)

outTab=data.frame()
a=rt[,gene]<=as.numeric(res.cut[gene,][1])
if (TRUE & FALSE %in% a){
  diff=survdiff(Surv(time,event) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outTab=rbind(outTab,cbind(gene=gene,pvalue=pValue))
  fit <- survfit(Surv(time, event) ~ a, data = rt)
  summary(fit)
}
pValue=round(pValue,3)
pValue=paste0("=",pValue)

pdf(file=paste(gene,"-COAD-best.survival.pdf",sep=""),
    width=6,
    height=6)
p=ggsurvplot(fit,
             data=rt,
             risk.table = T,
             pval = T,
             conf.int = F,
             ggtheme = theme_minimal(),
             palette = c("#FF0000","#2200FF"),
             legend.title=gene,
             legend.lab=c("High","Low"),
             xlab="Time(year)",
             ylab="Overall Survival",
             xlim=c(0,10),
             break.time.by=2)
print(p,newpage = F)
dev.off()

#GSE17536
load("GSE17536-deg-time-OS.Rdata")
library(dplyr)
clinicalExp=cbind(time,t(deg[gene,]))
library(survival)
library(survminer)
rt=clinicalExp
res.cut <- surv_cutpoint(rt, #数据集
                         time = "time", #生存状态
                         event = "event", #生存时间
                         variables = colnames(rt)[4:ncol(rt)] #需要计算的数据列名
)

res.cut=summary(res.cut)
outTab=data.frame()
a=rt[,gene]<=as.numeric(res.cut[gene,][1])
if (TRUE & FALSE %in% a){
  diff=survdiff(Surv(time,event) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outTab=rbind(outTab,cbind(gene=gene,pvalue=pValue))
  fit <- survfit(Surv(time, event) ~ a, data = rt)
  summary(fit)
}
pValue=round(pValue,3)
pValue=paste0("=",pValue)

pdf(file=paste(gene,"-GSE17536-best.survival.pdf",sep=""),
    width=6,
    height=6)
p=ggsurvplot(fit,
             data=rt,
             risk.table = T,
             pval = T,
             conf.int = F,
             ggtheme = theme_minimal(),
             palette = c("#FF0000","#2200FF"),
             legend.title=gene,
             legend.lab=c("High","Low"),
             xlab="Time(year)",
             ylab="Overall Survival",
             xlim=c(0,10),
             break.time.by=2)
print(p,newpage = F)
dev.off()

rt=clinicalExp
outTab=data.frame()
a=rt[,gene]<=median(rt[,gene])
if (TRUE & FALSE %in% a){
  diff=survdiff(Surv(time,event) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outTab=rbind(outTab,cbind(gene=gene,pvalue=pValue))
  fit <- survfit(Surv(time, event) ~ a, data = rt)
  summary(fit)
}
pValue=round(pValue,3)
pValue=paste0("=",pValue)

pdf(file=paste(gene,"-GSE17536.survival.pdf",sep=""),
    width=6,
    height=6)
p=ggsurvplot(fit,
             data=rt,
             risk.table = T,
             pval = T,
             conf.int = F,
             ggtheme = theme_minimal(),
             palette = c("#FF0000","#2200FF"),
             legend.title=gene,
             legend.lab=c("High","Low"),
             xlab="Time(year)",
             ylab="Overall Survival",
             xlim=c(0,10),
             break.time.by=2)
print(p,newpage = F)
dev.off()

load("GSE39582exp-pd.Rdata")
meta=read.table("GSE39582time.txt",header = T)
colnames(meta)=c("ID","time","event")
meta$time=meta$time/365
exp=t(deg)
rownames(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(exp))
library(limma)
exp=avereps(exp)
rownames(meta)=meta$ID
sameid=intersect(rownames(exp),rownames(meta))
exp=exp[sameid,]
meta=meta[sameid,]
clinicalExp=cbind(meta,exp[,gene])
colnames(clinicalExp)[4]=gene
library(survival)
library(survminer)
rt=clinicalExp
res.cut <- surv_cutpoint(rt, #数据集
                         time = "time", #生存状态
                         event = "event", #生存时间
                         variables = colnames(rt)[4:ncol(rt)] #需要计算的数据列名
)

res.cut=summary(res.cut)
outTab=data.frame()
a=rt[,gene]<=as.numeric(res.cut[gene,][1])
if (TRUE & FALSE %in% a){
  diff=survdiff(Surv(time,event) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outTab=rbind(outTab,cbind(gene=gene,pvalue=pValue))
  fit <- survfit(Surv(time, event) ~ a, data = rt)
  summary(fit)
}
pValue=round(pValue,3)
pValue=paste0("=",pValue)

pdf(file=paste(gene,"-GSE39582-best.survival.pdf",sep=""),
    width=6,
    height=6)
p=ggsurvplot(fit,
             data=rt,
             risk.table = T,
             pval = T,
             conf.int = F,
             ggtheme = theme_minimal(),
             palette = c("#FF0000","#2200FF"),
             legend.title=gene,
             legend.lab=c("High","Low"),
             xlab="Time(year)",
             ylab="Overall Survival",
             xlim=c(0,10),
             break.time.by=2)
print(p,newpage = F)
dev.off()

rt=clinicalExp
outTab=data.frame()
a=rt[,gene]<=median(rt[,gene])
if (TRUE & FALSE %in% a){
  diff=survdiff(Surv(time,event) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outTab=rbind(outTab,cbind(gene=gene,pvalue=pValue))
  fit <- survfit(Surv(time, event) ~ a, data = rt)
  summary(fit)
}
pValue=round(pValue,3)
pValue=paste0("=",pValue)

pdf(file=paste(gene,"-GSE39582.survival.pdf",sep=""),
    width=6,
    height=6)
p=ggsurvplot(fit,
             data=rt,
             risk.table = T,
             pval = T,
             conf.int = F,
             ggtheme = theme_minimal(),
             palette = c("#FF0000","#2200FF"),
             legend.title=gene,
             legend.lab=c("High","Low"),
             xlab="Time(year)",
             ylab="Overall Survival",
             xlim=c(0,10),
             break.time.by=2)
print(p,newpage = F)
dev.off()
