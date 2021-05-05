#----------------------------------------------------------------------------
# Expression patterns in whole embryos, hypoblast and epiblast
# Tim Coorens
# tc16@sanger.ac.uk
# August 2020

options(stringsAsFactors = F)
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(viridis)
library(data.table)

source("logisticRegression.R")
source("similarity.R")

#Read in data
embryo.integrated=readRDS("embryo_integrated_allembryos_filtered.Rdata")

embryo_hypo=subset(embryo.integrated,idents="Hypoblast")
matrix<-GetAssayData(object = embryo_hypo)

day=9
matrix_mod<-as.matrix(matrix[,embryo_hypo@meta.data$Age==day])
gene<-as.numeric(matrix_mod["CER1",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
correlations_pval<-apply(matrix_mod,1,function(x){cor.test(gene,x)$p.value})

nodal=c("LEFTY1","LEFTY2","NCLN","TMEFF1","HHEX")
bmp=c("SOSTDC1","SOST","GREM2","DAND5","GREM1","NOG","CHRD","TWSG1")
wnt=c("SFRP1","SFRP2","SFRP4","SFRP5","FRZB","DKK1","DKK2","DKK3","DKK4","DKKL1")
nodal_dat=round(correlations[nodal],digits=2)
nodal_dat=nodal_dat[order(nodal_dat,decreasing=TRUE)]
bmp_dat=round(correlations[bmp],digits=2)
bmp_dat=bmp_dat[order(bmp_dat,decreasing=TRUE)]
wnt_dat=round(correlations[wnt],digits=2)
wnt_dat=wnt_dat[order(wnt_dat,decreasing=TRUE)]

correlations_qval=p.adjust(correlations_pval,method = "BH")
pvals_vec_d11=correlations_qval[c(names(nodal_dat),names(bmp_dat),names(wnt_dat))]

pvals_vec_d9=correlations_qval[c(names(nodal_dat),names(bmp_dat),names(wnt_dat))]


pval=data.frame(day9=pvals_vec_d9,
                day11=pvals_vec_d11)
write.table(pval,"pval_CER1_correlation.txt",quote=F,sep="\t")

dat=c(nodal_dat,bmp_dat,wnt_dat)
inh=c(nodal,bmp,wnt)
dat=round(correlations[inh],digits=2)
cols=colorRampPalette(c("Blue","White","Red"))(201)
names(cols)=round(seq(-1,1,0.01),digits=2)
col_dat=cols[as.character(dat)]
col_dat[is.na(col_dat)]="grey"

val_labels=dat
val_labels[is.na(val_labels)]="N/A"
pdf(paste0("CER1_coexpression_2020_07_12_day",day,".pdf"),useDingbats = F,width=14,height=3.5) #Fig 3i
plot(0,type='n',axes=FALSE,ann=FALSE,ylim=c(0,2.1),xlim=c(-0.25,length(dat)+0.2))
rect(xleft=seq(0,length(dat)-1,1),xright=seq(1,length(dat),1),ybottom=0,ytop=1,col=col_dat)
text(x=seq(0.35,length(dat)-0.65,1),y=1.05,labels=names(dat),srt=45,pos=4)
text(x=seq(0.5,length(dat)-0.5,1),y=0.5,labels=val_labels)
text(x=-0.5,y=0.5,labels="CER1")
dev.off()

#-----
#Subset into epi- and hypoblast objects
embryo.integrated_epi=subset(embryo.integrated,idents="Epiblast")
DefaultAssay(embryo.integrated_epi) <- "RNA"

embryo.integrated_hypo=subset(embryo.integrated,idents="Hypoblast")
DefaultAssay(embryo.integrated_hypo) <- "RNA"

#-----
#FGF signalling - Fig 2, Extended data fig 6
fgf_genes=rownames(embryo.integrated_epi)[grepl("FGF",rownames(embryo.integrated_epi))]
for(gene in fgf_genes){
  pdf(paste0("FGF/",gene,"_all_featureplot.pdf"),width=6,height=5,useDingbats = F)
  print(FeaturePlot(object = embryo.integrated, features = gene, pt.size = 1.75,reduction = "umap",cols = c("lightgrey", "blue"),order=T))
  dev.off()
  
  pdf(paste0("FGF/",gene,"_epiblast_featureplot.pdf"),width=6,height=5,useDingbats = F)
  print(FeaturePlot(object = embryo.integrated_epi, features = gene, pt.size = 1.75,reduction = "umap",cols = c("grey", "blue"),order=T))
  dev.off()
  
  pdf(paste0("FGF/",gene,"_hypoblast_featureplot.pdf"),width=6,height=5,useDingbats = F)
  print(FeaturePlot(object = embryo.integrated_hypo, features = gene, pt.size = 1.75, reduction = "umap",cols = c("grey", "blue"),order=T))
  dev.off()
  print(gene)
}

#----

embryo.integrated_epi <- FindNeighbors(object = embryo.integrated_epi,k.param=10)
embryo.integrated_epi <- FindClusters(embryo.integrated_epi,res=0.25)
pdf("epiblast_umap.pdf",width=4,height=3,useDingbats = F)
DimPlot(embryo.integrated_epi, reduction = "umap", group.by = "seurat_clusters")
dev.off()

p1 <- DimPlot(embryo.integrated_epi, reduction = "tsne", group.by = "seurat_clusters")
markers_epi=FindAllMarkers(embryo.integrated_epi,test.use="roc")
top30 <- markers_epi %>% group_by(cluster) %>% top_n(30, avg_diff)
#----

#Hypoblast clustering and expression patterns
embryo.integrated_hypo <- FindNeighbors(object = embryo.integrated_hypo)
embryo.integrated_hypo <- FindClusters(embryo.integrated_hypo,res=0.25)

markers_hypo=FindAllMarkers(embryo.integrated_hypo,test.use="roc")
top50 <- markers_hypo %>% group_by(cluster) %>% top_n(50, avg_diff)
write.csv(top50,"hypoblast_cluster_markers_top50.csv")
pch_vec=rep(16,ncol(embryo.integrated_hypo))
pch_vec[embryo.integrated_hypo@meta.data$Age==11]=17


pdf("hypoblast_umap_all_cells.pdf",width=5,height=4,useDingbats = F)
DimPlot(embryo.integrated_hypo, reduction = "umap", pt.size=2, group.by = "seurat_clusters")
dev.off()
pdf("hypoblast_umap_day9.pdf",width=5,height=4,useDingbats = F)
DimPlot(subset(embryo.integrated_hypo,Age==9), reduction = "umap", pt.size=2, group.by = "seurat_clusters")
dev.off()
pdf("hypoblast_umap_day11.pdf",width=5,height=4,useDingbats = F)
DimPlot(subset(embryo.integrated_hypo,Age==11), reduction = "umap", pt.size=2, group.by = "seurat_clusters")
dev.off()

#Expression in hypoblast (Fig. 3, Extended data fig 7,9)
for(gene2 in c("HHEX","LEFTY1","LEFTY2","NOG","DKK4")){
  pdf(paste0("hypoblast_CER1_",gene2,"_all_cells.pdf"),width=11,height=4,useDingbats = F)
  print(FeaturePlot(embryo.integrated_hypo, reduction = "umap", features = c("CER1",gene2),pt.size=2, 
                    cols = c("firebrick", "blue"),blend=T,max.cutoff="q90",order=T))
  dev.off()
  pdf(paste0("hypoblast_CER1_",gene2,"_day9.pdf"),width=11,height=4,useDingbats = F)
  print(FeaturePlot(subset(embryo.integrated_hypo,Age==9), reduction = "umap", features = c("CER1",gene2),pt.size=2, 
                    cols = c("firebrick", "blue"),blend=T,max.cutoff="q90",order=T))
  dev.off()
  pdf(paste0("hypoblast_CER1_",gene2,"_day11.pdf"),width=11,height=4,useDingbats = F)
  print(FeaturePlot(subset(embryo.integrated_hypo,Age==11), reduction = "umap", features = c("CER1",gene2),pt.size=2, 
                    cols = c("firebrick", "blue"),blend=T,max.cutoff="q90",order=T))
  dev.off()
}

ave_genes=c("NOG","LEFTY1","LEFTY2","DKK1","SFRP1","HHEX","FRZB","DKK3")
for(gene in ave_genes){
  pdf(paste0(gene,"_hypoblast_featureplot.pdf"),width=6,height=5,useDingbats = F)
  print(FeaturePlot(object = embryo.integrated_hypo, features = gene, pt.size = 1.75, reduction = "umap",cols = c("lightgrey", "blue"),max.cutoff="q90"))
  dev.off()
  
  pdf(paste0(gene,"_hypoblast_featureplot_day9.pdf"),width=6,height=5,useDingbats = F)
  print(FeaturePlot(object = subset(embryo.integrated_hypo,Age==9), features = gene, pt.size = 1.75, reduction = "umap",cols = c("lightgrey", "blue"),max.cutoff="q90"))
  dev.off()
  
  pdf(paste0(gene,"_hypoblast_featureplot_day11.pdf"),width=6,height=5,useDingbats = F)
  print(FeaturePlot(object = subset(embryo.integrated_hypo,Age==11), features = gene, pt.size = 1.75, reduction = "umap",cols = c("lightgrey", "blue"),max.cutoff="q90"))
  dev.off()
}

for(gene in c("PDGFRA","SPARC","RSPO3","FRZB","NID2")){
  pdf(paste0(gene,"_hypoblast_featureplot.pdf"),width=4,height=3,useDingbats = F)
  print(FeaturePlot(object = embryo.integrated_hypo, features = gene, pt.size = 1.75, reduction = "umap",cols = c("grey", "blue")))
  dev.off()
}
for(gene in c("GATA2","CGA","GATA3","PGF")){
  pdf(paste0(gene,"_hypoblast_featureplot.pdf"),width=4,height=3,useDingbats = F)
  print(FeaturePlot(object = embryo.integrated_hypo, features = gene, pt.size = 1.75, reduction = "umap",cols = c("grey", "blue")))
  dev.off()
}

