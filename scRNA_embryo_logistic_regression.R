#----------------------------------------------------------------------------
# Compare with other data sets and perform logistic regressions
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

#------------------------------
# Compare with Zhou et al. 

library(data.table)
fuchou_tpm=fread("GSE109555_All_Embryo_TPM.txt.gz",data.table = F) #Data available at GEO GSE109555
rownames(fuchou_tpm)=fuchou_tpm$V1
sampleInfo=read.csv("Supplementary Table 2 Sample Information.csv") #Available as part of the Zhou et al. paper

#Select the good embryos
include_zhou=c("hm_D8_E2", "hm_D8_E3", "hm_D8_E5", "ha_D8_E1", "hv_D8_E1", "hv_D8_E2", "hv_D8_E3", "hv_D10_E6", "ha_D10_E1", "ha_D10_E2", "hm_D10_E4", "hm_D10_E9", "hv_D10_E7", "hv_D10_E8", "ha_D12_E1", "hv_D12_E1", "hv_D12_E2")
fuchou_tpm=fuchou_tpm[,-1]

#Set up seurat objet
embryo_fuchou=CreateSeuratObject(counts = fuchou_tpm[,sampleInfo$Sample[sampleInfo$Ori_Day_Emb%in%include_zhou]])
embryo_fuchou <- NormalizeData(object = embryo_fuchou)
embryo_fuchou <- ScaleData(object = embryo_fuchou)

#Merge with own data - note: no batch correction as logistic regression doesn't rely on it
embryo_integrated_zhou=merge(x = embryo.integrated, y = embryo_fuchou, 
                             add.cell.ids=c("Own","Zhou"),merge.data = TRUE)

#Annotate clusters
embryo_integrated_zhou@meta.data$Cell_Type[embryo_integrated_zhou@meta.data$seurat_clusters==0]="TE"
embryo_integrated_zhou@meta.data$Cell_Type[embryo_integrated_zhou@meta.data$seurat_clusters==1]="TE"
embryo_integrated_zhou@meta.data$Cell_Type[embryo_integrated_zhou@meta.data$seurat_clusters==2]="Epi"
embryo_integrated_zhou@meta.data$Cell_Type[embryo_integrated_zhou@meta.data$seurat_clusters==3]="Hypo"
embryo_integrated_zhou@meta.data$Cell_Type[rownames(embryo_integrated_zhou@meta.data)%in%paste0("Zhou_",sampleInfo$Sample[sampleInfo$Lineage=="PE"])]="Zhou_Hypo"
embryo_integrated_zhou@meta.data$Cell_Type[rownames(embryo_integrated_zhou@meta.data)%in%paste0("Zhou_",sampleInfo$Sample[sampleInfo$Lineage=="TE"])]="Zhou_TE"
embryo_integrated_zhou@meta.data$Cell_Type[rownames(embryo_integrated_zhou@meta.data)%in%paste0("Zhou_",sampleInfo$Sample[sampleInfo$Lineage=="EPI"])]="Zhou_Epi"

dat=embryo_integrated_zhou@assays$RNA@counts

select=rownames(embryo_integrated_zhou@meta.data)[grepl("Own",rownames(embryo_integrated_zhou@meta.data))]
zhou_select=rownames(embryo_integrated_zhou@meta.data)[grepl("Zhou",rownames(embryo_integrated_zhou@meta.data))]

#Own data is training data for the logistic regression, data from Zhou et al is regressed into it as test data
training_dat=dat[,select]
test_dat=dat[,zhou_select]
classes=embryo_integrated_zhou@meta.data[select,"Cell_Type"]

fit=trainModel(refMat=training_dat,classes=classes)
saveRDS(fit,"logistic_regression_own_data_fit_zhou.Rdata")

preds=predictSimilarity(fit=fit,tgtData=test_dat,classes=embryo_integrated_zhou@meta.data[zhou_select,"Cell_Type"],logits=F)


wilcox.test(preds[classes=="Zhou_Epi","Epi"],preds[classes=="Zhou_Epi","Hypo"])
wilcox.test(preds[classes=="Zhou_Epi","Epi"],preds[classes=="Zhou_Epi","TE"])

wilcox.test(preds[classes=="Zhou_Hypo","Hypo"],preds[classes=="Zhou_Hypo","Epi"])
wilcox.test(preds[classes=="Zhou_Hypo","Hypo"],preds[classes=="Zhou_Hypo","TE"])

wilcox.test(preds[classes=="Zhou_TE","TE"],preds[classes=="Zhou_TE","Epi"])
wilcox.test(preds[classes=="Zhou_TE","TE"],preds[classes=="Zhou_TE","Hypo"])

preds=preds[,c("Epi","Hypo","TE")]
note_color=rep('black',9)
note_color[t(preds)>0.5]="white"
pdf("zhou_own_logistic_regression.pdf")
heatmap.2(t(as.matrix(preds)),
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          dendrogram="none",
          notecol = "black",
          col=colorRampPalette(brewer.pal(9,"Blues")[2:9])(100),
          Rowv = F,
          Colv = F,
          scale='none',
          key=F,
          cellnote=round(t(as.matrix(preds)),digits=2),
          cexRow=1.5,
          cexCol=1.5,
          notecex=1.0,
          mar=c(10,6.5))
dev.off()

#---------------------------------
# Compare with Xiang et al. 

xiang_fpkm=fread("GSE136447_555-samples-fpkm.txt.gz",data.table = F) #Get from GEO
xiang_genes=xiang_fpkm$"Gene Name"
xiang_genes[duplicated(xiang_genes)]=xiang_fpkm$"Gene ID"[duplicated(xiang_genes)]
xiang_genes[duplicated(xiang_genes)]=xiang_fpkm$"Gene_ID"[duplicated(xiang_genes)]
rownames(xiang_fpkm)=xiang_genes
xiang_fpkm=xiang_fpkm[,-c(1:7)]

sampleInfo=read.table("xiang_sampleinfo.txt",header=T)
include_xiang=c("D8_8", "D9A4", "D10A2", "D10A6", "D12A3", "D12A5", "D14A1")

#Set up seurat objet
embryo_xiang=CreateSeuratObject(counts = xiang_fpkm[,sampleInfo$Sample.ID[sampleInfo$Embryo.ID%in%include_xiang]])
embryo_xiang <- NormalizeData(object = embryo_xiang)
embryo_xiang <- ScaleData(object = embryo_xiang)

#Merge with own data - note: no batch correction as logistic regression doesn't rely on it
embryo_integrated_xiang=merge(x = embryo.integrated, y = embryo_xiang, 
                              add.cell.ids=c("Own","Xiang"),merge.data = TRUE)

#Annotate own data into the four clusters
embryo_integrated_xiang@meta.data$Cell_Type[embryo_integrated_xiang@meta.data$seurat_clusters==0]="STB"
embryo_integrated_xiang@meta.data$Cell_Type[embryo_integrated_xiang@meta.data$seurat_clusters==1]="CTB"
embryo_integrated_xiang@meta.data$Cell_Type[embryo_integrated_xiang@meta.data$seurat_clusters==2]="Epi"
embryo_integrated_xiang@meta.data$Cell_Type[embryo_integrated_xiang@meta.data$seurat_clusters==3]="Hypo"

#
embryo_integrated_xiang@meta.data$Cell_Type[rownames(embryo_integrated_xiang@meta.data)%in%paste0("Xiang_",sampleInfo$Sample.ID[sampleInfo$Group=="PE"])]="Xiang_PrE"
#embryo_integrated_xiang@meta.data$Cell_Type[rownames(embryo_integrated_xiang@meta.data)%in%paste0("Xiang_",sampleInfo$Sample.ID[sampleInfo$Group=="ICM"])]="Xiang_ICM"
embryo_integrated_xiang@meta.data$Cell_Type[rownames(embryo_integrated_xiang@meta.data)%in%paste0("Xiang_",sampleInfo$Sample.ID[sampleInfo$Group=="EPI"])]="Xiang_Epi"
embryo_integrated_xiang@meta.data$Cell_Type[rownames(embryo_integrated_xiang@meta.data)%in%paste0("Xiang_",sampleInfo$Sample.ID[sampleInfo$Group=="STB"])]="Xiang_STB"
embryo_integrated_xiang@meta.data$Cell_Type[rownames(embryo_integrated_xiang@meta.data)%in%paste0("Xiang_",sampleInfo$Sample.ID[sampleInfo$Group=="CTB"])]="Xiang_CTB"
#embryo_integrated_xiang@meta.data$Cell_Type[rownames(embryo_integrated_xiang@meta.data)%in%paste0("Xiang_",sampleInfo$Sample.ID[sampleInfo$Group=="EVT"])]="Xiang_EVT"
embryo_integrated_xiang@meta.data$Cell_Type[rownames(embryo_integrated_xiang@meta.data)%in%paste0("Xiang_",sampleInfo$Sample.ID[sampleInfo$Group=="PSA-EPI"])]="Xiang_PSA-EPI"

dat=embryo_integrated_xiang@assays$RNA@counts

select=rownames(embryo_integrated_xiang@meta.data)[grepl("Own",rownames(embryo_integrated_xiang@meta.data))]
xiang_select=rownames(embryo_integrated_xiang@meta.data)[grepl("Xiang",rownames(embryo_integrated_xiang@meta.data))&!is.na(embryo_integrated_xiang@meta.data$Cell_Type)]

training_dat=dat[,select]
test_dat=dat[,xiang_select]

classes=embryo_integrated_xiang@meta.data[select,"Cell_Type"]

fit=trainModel(refMat=training_dat,classes=classes)
saveRDS(fit,"logistic_regression_own_data_fit_xiang.Rdata")

preds=predictSimilarity(fit=fit,tgtData=test_dat,classes=embryo_integrated_xiang@meta.data[xiang_select,"Cell_Type"],logits=F)

wilcox.test(preds[xiang_classes=="Xiang_Epi","Epi"],preds[xiang_classes=="Xiang_Epi","Hypo"])
wilcox.test(preds[xiang_classes=="Xiang_Epi","Epi"],preds[xiang_classes=="Xiang_Epi","CTB"])
wilcox.test(preds[xiang_classes=="Xiang_Epi","Epi"],preds[xiang_classes=="Xiang_Epi","STB"])

wilcox.test(preds[xiang_classes=="Xiang_PrE","Hypo"],preds[xiang_classes=="Xiang_PrE","Epi"])
wilcox.test(preds[xiang_classes=="Xiang_PrE","Hypo"],preds[xiang_classes=="Xiang_PrE","CTB"])
wilcox.test(preds[xiang_classes=="Xiang_PrE","Hypo"],preds[xiang_classes=="Xiang_PrE","STB"])

wilcox.test(preds[xiang_classes=="Xiang_CTB","CTB"],preds[xiang_classes=="Xiang_CTB","Epi"])
wilcox.test(preds[xiang_classes=="Xiang_CTB","CTB"],preds[xiang_classes=="Xiang_CTB","Hypo"])
wilcox.test(preds[xiang_classes=="Xiang_CTB","CTB"],preds[xiang_classes=="Xiang_CTB","STB"])

wilcox.test(preds[xiang_classes=="Xiang_STB","STB"],preds[xiang_classes=="Xiang_STB","Epi"])
wilcox.test(preds[xiang_classes=="Xiang_STB","STB"],preds[xiang_classes=="Xiang_STB","Hypo"])
wilcox.test(preds[xiang_classes=="Xiang_STB","STB"],preds[xiang_classes=="Xiang_STB","CTB"])

wilcox.test(preds[xiang_classes=="Xiang_PSA-EPI","Epi"],preds[xiang_classes=="Xiang_PSA-EPI","STB"])
wilcox.test(preds[xiang_classes=="Xiang_PSA-EPI","Epi"],preds[xiang_classes=="Xiang_PSA-EPI","Hypo"])
wilcox.test(preds[xiang_classes=="Xiang_PSA-EPI","Epi"],preds[xiang_classes=="Xiang_PSA-EPI","CTB"])

write.table(preds,"preds_xiang.txt")

rownames(preds)=c("Xiang_STB","Xiang_Epi","Xiang_PSA-EPI","Xiang_CTB")
preds=preds[c("Xiang_PSA-EPI","Xiang_Epi","Xiang_CTB","Xiang_STB"),]

pdf("xiang_own_logistic_regression.pdf")
heatmap.2(t(as.matrix(preds)),
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          dendrogram="none",
          notecol = "black",
          col=colorRampPalette(brewer.pal(9,"Blues")[2:9])(100),
          Rowv = F,
          Colv = F,
          scale='none',
          key=F,
          cellnote=round(t(as.matrix(preds)),digits=2),
          cexRow=1.5,
          cexCol=1.5,
          notecex=1.0,
          mar=c(10,6.5))
dev.off()

#----------------------
# Compare integrated epiblast dataset to cell lines

embryo_integrated_stirparo=readRDS("embryo_integrated_stirparo.Rdata")
epi_cells=rownames(embryo_integrated_stirparo@meta.data)[embryo_integrated_stirparo@meta.data$Epi!="0"]
epi_state=embryo_integrated_stirparo@meta.data$Epi[embryo_integrated_stirparo@meta.data$Epi!="0"]
embryo_epi=subset(embryo_integrated_stirparo,cells = epi_cells)

training_dat=embryo_epi@assays$RNA@counts
classes=epi_state

fit=trainModel(refMat=training_dat,classes=classes)
saveRDS(fit,"logistic_regression_epiblast_fit.Rdata")

#Data from Takashima et al, 2014
cell_lines_raw = read.table("Abundance_CellLine.txt")
cell_lines_ensembl=rownames(cell_lines_raw)

#Read in ensembl annotation for gene name conversion
ensembl_geneid=read.table("ensembl_geneid.txt",header=T,sep="\t")
ensembl_geneid=ensembl_geneid[ensembl_geneid$Ensembl.gene.ID!="",]
genes=ensembl_geneid$Approved.symbol
names(genes)=ensembl_geneid$Ensembl.gene.ID

genes_cell_line=genes[cell_lines_ensembl]
genes_cell_line[duplicated(genes_cell_line)] = cell_lines_ensembl[duplicated(genes_cell_line)]
cell_lines_raw=cell_lines_raw[!is.na(genes_cell_line),]
genes_cell_line=genes_cell_line[!is.na(genes_cell_line)]
rownames(cell_lines_raw)=genes_cell_line

test_dat=cell_lines_raw
preds=predictSimilarity(fit=fit,tgtData=test_dat,logits=F,minGeneMatch=0.1)
write.table(preds,"preds_cell_lines_takashima.txt")

sample_names=c("H9_reset_R2","H9_reset_R3","H9_reset_R1","H9_R3","H9_R1","H9_R2")

rownames(preds)=sample_names
preds=preds[order(rownames(preds)),c("ICM","Pre-Epi","Post-Epi")]

note_color=rep('black',24)
note_color[preds>0.5]="white"
#Fig. 1g
pdf("Logistic_regression_cell_line_takashima.pdf",useDingbats = F,height=4.2,width=6)
heatmap.2(t(as.matrix(preds)),
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          dendrogram="none",
          notecol = rev(note_color),
          col=colorRampPalette(brewer.pal(9,"Blues")[2:9])(100),
          Rowv = F,
          Colv = F,
          scale='none',
          key=F,
          cellnote=round(t(as.matrix(preds)),digits=2),
          cexRow=1.5,
          notecex=1.0,
          ColSideColors=c(rep('firebrick',3),rep('steelblue',3)), 
          mar=c(10,6.5))
dev.off()


#Theunissen et al.
#Get data from GEO
theunissen=matrix(0,nrow=63682,ncol=5)
n=1
for(x in paste0("GSM22186",c(60,68:71))){
  data=fread(paste0(x,".pileup.gz"),sep="\t",data.table=F) 
  theunissen[,n]=data[,2]
  n=n+1
}
rownames(theunissen)=data[,1]
colnames(theunissen)=c("Male Primed ES - WIBR1 (XY)","Male Naive ES - WIN1 (XY)","Male Naive ES - WIN1 (XY) - 2","Naive - WIBR2","Naive - WIBR3")
theunissen=theunissen[rownames(theunissen)%in%ensembl_convert$Gene.stable.ID,]

gene_meta=read.table("../../Hypoblast/May_2020/GTEX_BulkGenes.tsv",sep="\t",header=T)
gene_names=gene_meta$geneID
gene_names[duplicated(gene_names)]=gene_meta$geneID[duplicated(gene_names)]
gene_length=gene_meta$bp_length
names(gene_length)=gene_names

gene_length=gene_length[!is.na(gene_length)]

tpm3 <- function(counts,len=gene_length){
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}
theunissen=theunissen[rownames(theunissen)%in%names(gene_length),]

theunissen_tpm=tpm3(theunissen,len=gene_length[rownames(theunissen)])
rownames(theunissen)=ensembl[rownames(theunissen)]
rownames(theunissen_tpm)=ensembl[rownames(theunissen_tpm)]
preds=predictSimilarity(fit=fit,tgtData=theunissen,logits=F,minGeneMatch=0.1)

preds=preds[,c("ICM","Pre-Epi","Peri-Epi","Post-Epi")]

note_color=rep('black',nrow(preds)*ncol(preds))
note_color[preds>0.5]="white"

#Extended data Fig 5a
pdf("Logistic_regression_theunissen_4stg.pdf",useDingbats = F,height=6,width=5)
heatmap.2(t(as.matrix(preds)),
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          dendrogram="none",
          notecol = rev(note_color),
          col=colorRampPalette(brewer.pal(9,"Blues")[2:9])(100),
          Rowv = F,
          Colv = F,
          scale='none',
          key=F,
          breaks=col_breaks,
          cellnote=round(t(as.matrix(preds)),digits=3),
          cexRow=1,
          notecex=0.7,
          #ColSideColors=c(rep('firebrick',3),rep('steelblue',3)), 
          mar=c(20,6.5))
dev.off()
write.csv(preds,"Logistic_regression_theunissen_4stg.csv")

#Rostovskaya et al. 
rostovskaya=fread("GSE123055_counts.txt.gz",sep="\t",data.table = F)
rownames(rostovskaya)=rostovskaya$ID
rostovskaya=rostovskaya[,-1]
rostovskaya=rostovskaya[,colnames(rostovskaya)!="V53"]
rostovskaya=rostovskaya[rownames(rostovskaya)%in%ensembl_convert$Gene.stable.ID,]

gene_length=gene_length[1:(length(gene_length)-2)]
gene_length=gene_length[!is.na(gene_length)]

tpm3 <- function(counts,len=gene_length){
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}
rostovskaya=rostovskaya[rownames(rostovskaya)%in%names(gene_length),]
rostovskaya_tpm=tpm3(rostovskaya,len=gene_length[rownames(rostovskaya)])
rownames(rostovskaya)=ensembl[rownames(rostovskaya)]
rownames(rostovskaya_tpm)=ensembl[rownames(rostovskaya_tpm)]
preds=predictSimilarity(fit=fit,tgtData=rostovskaya_tpm,logits=F,minGeneMatch=0.1)

preds=preds[order(rownames(preds)),c("ICM","Pre-Epi","Peri-Epi","Post-Epi")]

note_color=rep('black',nrow(preds)*ncol(preds))
note_color[preds>0.5]="white"

#Extended data Fig 5b
pdf("Logistic_regression_rostovskaya_4stg.pdf",useDingbats = F,height=6.2,width=25)
heatmap.2(t(as.matrix(preds)),
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          dendrogram="none",
          notecol = rev(note_color),
          col=colorRampPalette(brewer.pal(9,"Blues")[2:9])(100),
          Rowv = F,
          Colv = F,
          scale='none',
          key=F,
          cellnote=round(t(as.matrix(preds)),digits=3),
          cexRow=1.5,
          notecex=0.7,
          breaks=col_breaks,
          #ColSideColors=c(rep('firebrick',3),rep('steelblue',3)), 
          mar=c(20,6.5))
dev.off()
write.csv(preds,"Logistic_regression_rostovskaya_4stg.csv")

#


