#----------------------------------------------------------------------------
# Read in cellranger data, perform batch correction and plot basic figures
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

#These are embryos with all lineages (hypoblast, epiblast, and trophoblasts):
include=c(1:4,7,9,10,12:15)
#To check selection, set:
# include=1:15
embryo_data_list=list()
for (n in 1:11){ #Singleplexed embryos
  embryo_data_list[[n]]=Read10X(data.dir = paste0("Embryo",n,"/filtered_feature_bc_matrix/"))
  print(n)
}

for (n in 12:15){ #Multiplexed embryos
  embryo_data_list[[n]]=Read10X(data.dir = paste0("MPlex",n-11,"/filtered_feature_bc_matrix/"))
  print(n)
}

#Turn cellranger counts into seurat objects
embryo_data_seurat=list()
for (n in 1:length(embryo_data_list)){
  embryo_data_seurat[[n]]=CreateSeuratObject(counts = embryo_data_list[[n]], project = paste0("embryo",n))
  print(n)
}

#Add deconvolution information from souporcell to multiplexed embryos
mplex_include=cbind(c(1,1,1,1,2,2,3,4,4),c(1:4,3,5,2,2,5))
for (n in 12:15){
  clusters=read.table(paste0("MPlex",n-11,"/clusters.tsv"),header=T)
  include_sub=mplex_include[mplex_include[,1]==n-11,2]-1
  embryo_data_seurat[[n]]@meta.data$indiv=clusters$assignment
  embryo_data_seurat[[n]]=subset(embryo_data_seurat[[n]], cells=clusters$barcode[clusters$status=="singlet"&clusters$assignment%in%include_sub])
  print(n)
  print(table(clusters$status))
}


batch_vec=c(1,1,2,2,rep(3,4),rep(4,3),rep(5,2),rep(6,2))

#Merge seurat objects per batch (as defined above)
batched_seurat_list=list()
for(k in 1:length(unique(batch_vec))){
  select=which(batch_vec==k&1:15%in%include)
  if(length(select)>0){
    if(length(select)>1){
      batched_seurat_list[[k]]=merge(x = embryo_data_seurat[[select[1]]],y = c(embryo_data_seurat[select[2:length(select)]]),
                                     add.cell.ids = paste0("embryo",select[1:length(select)]),merge.data = TRUE,project = paste0("embryo_b",k))
    }else{
      batched_seurat_list[[k]]=embryo_data_seurat[[select[1]]]
      batched_seurat_list[[k]]@meta.data$orig.ident=paste0("embryo",select[1])
      batched_seurat_list[[k]]=RenameCells(batched_seurat_list[[k]],add.cell.id=paste0("embryo",select[1]))
    }
    batched_seurat_list[[k]]@meta.data[, "batch"] <- k
  }
  print(k)
}

#Per batch, perform basic filtering
for(k in unique(batch_vec)){
  mito.genes <- grep(pattern = "^MT-", x = rownames(batched_seurat_list[[k]]), value = TRUE)
  percent.mito <- colSums(as.matrix(GetAssayData(object = batched_seurat_list[[k]], slot = "counts")[mito.genes, ]))/colSums(as.matrix(GetAssayData(object = batched_seurat_list[[k]], slot = "counts")))
  batched_seurat_list[[k]] <- AddMetaData(object = batched_seurat_list[[k]], metadata = percent.mito, col.name = "percent.mito")
  batched_seurat_list[[k]] <- subset(batched_seurat_list[[k]], subset=nFeature_RNA>1000&percent.mito<0.2)
  batched_seurat_list[[k]] <- NormalizeData(object = batched_seurat_list[[k]])
  batched_seurat_list[[k]] <- ScaleData(object = batched_seurat_list[[k]])
  batched_seurat_list[[k]] <- FindVariableFeatures(object = batched_seurat_list[[k]])
  print(k)
}
#Select all genes for integration
all_genes=rownames(batched_seurat_list[[1]])

#Alternatively, select variable features for integration
# hvg_all=c()
# for(k in unique(batch_vec)){
#   hvg_all=union(hvg_all,head(VariableFeatures(batched_seurat_list[[k]]),n=5000))
# }

#Perform batch correction
anchors <- FindIntegrationAnchors(object.list = batched_seurat_list, dims = 1:30)
embryo.integrated <- IntegrateData(anchorset = anchors, dims = 1:30,features.to.integrate=all_genes)

#Note: 'integrated' is best for dimensionality reduction, not for plotting of expression values
DefaultAssay(embryo.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
embryo.integrated <- ScaleData(embryo.integrated, features = all_genes,verbose = FALSE)
embryo.integrated <- RunPCA(embryo.integrated, npcs = 30, verbose = FALSE)
embryo.integrated <- FindNeighbors(object = embryo.integrated)
embryo.integrated <- FindClusters(embryo.integrated,res=0.05)
embryo.integrated <- RunTSNE(embryo.integrated, reduction = "pca", dims = 1:30)
embryo.integrated=RenameIdents(object=embryo.integrated,'0'="Syncytiotrophoblasts",'1'="Cytotrophoblasts",'2'="Epiblast",'3'="Hypoblast")

embryo.integrated <- RunUMAP(embryo.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(embryo.integrated, reduction = "umap", group.by = "batch")
p2 <- DimPlot(embryo.integrated, reduction = "umap", label = F, 
              repel = TRUE,cols=c("dodgerblue","olivedrab3","firebrick","darkviolet"))
pdf("UMAP_combined.pdf",width=10,height=6)
p1 + p2
dev.off()

#Switch from Embryo# to MPlex#.# for the multiplexed embryos
embryo.integrated@meta.data$orig.ident[embryo.integrated@meta.data$orig.ident=="embryo12"]=paste0("mplex1.",as.numeric(embryo.integrated@meta.data$indiv[embryo.integrated@meta.data$orig.ident=="embryo12"])+1)
embryo.integrated@meta.data$orig.ident[embryo.integrated@meta.data$orig.ident=="embryo13"]=paste0("mplex2.",as.numeric(embryo.integrated@meta.data$indiv[embryo.integrated@meta.data$orig.ident=="embryo13"])+1)
embryo.integrated@meta.data$orig.ident[embryo.integrated@meta.data$orig.ident=="embryo14"]=paste0("mplex3.",as.numeric(embryo.integrated@meta.data$indiv[embryo.integrated@meta.data$orig.ident=="embryo14"])+1)
embryo.integrated@meta.data$orig.ident[embryo.integrated@meta.data$orig.ident=="embryo15"]=paste0("mplex4.",as.numeric(embryo.integrated@meta.data$indiv[embryo.integrated@meta.data$orig.ident=="embryo15"])+1)

#Annotate with ages
embryo.integrated@meta.data$Age=11
embryo.integrated@meta.data$Age[embryo.integrated@meta.data$orig.ident%in%paste0("embryo",c(1,2,7,12,14))]=9

#Save the output
saveRDS(embryo.integrated,"embryo_integrated_allembryos_filtered.Rdata")

#Read in previous output
embryo.integrated=readRDS("embryo_integrated_allembryos_filtered.Rdata")

#Print reads/cell and genes/cell
idents=unique(embryo.integrated@meta.data$orig.ident)
for(i in idents){
  print(paste0(i,": ",median(embryo.integrated@meta.data$nCount_RNA[embryo.integrated@meta.data$orig.ident==i]),
               " reads/cell and ",median(embryo.integrated@meta.data$nFeature_RNA[embryo.integrated@meta.data$orig.ident==i])," genes/cell"))
}

#Basic plots 
pdf("embryo_integrated_umap_2020_07_12.pdf",width=6,height=4,useDingbats = F)
DimPlot(embryo.integrated, reduction = "umap", label = F, 
        repel = TRUE,cols=c("dodgerblue","olivedrab3","firebrick","darkviolet"))
dev.off()
pdf("embryo_integrated_umap_orig_ident_2020_07_20.pdf",width=5,height=4,useDingbats = F)
DimPlot(embryo.integrated, reduction = "umap", pt.size=.5, group.by = "orig.ident")
dev.off()
pdf("embryo_integrated_umap_age_2020_07_20.pdf",width=5,height=4,useDingbats = F)
DimPlot(embryo.integrated, reduction = "umap",pt.size=.5,group.by = "Age",cols=c("darkblue","lightblue"))
dev.off()
pdf("embryo_integrated_pca_orig_ident_2020_07_20.pdf",width=5,height=4,useDingbats = F)
DimPlot(embryo.integrated, reduction = "pca", pt.size=.5, group.by = "orig.ident")
dev.off()
pdf("embryo_integrated_tsne_orig_ident_2020_07_20.pdf",width=5,height=4,useDingbats = F)
DimPlot(embryo.integrated, reduction = "tsne", pt.size=.5, group.by = "orig.ident")
dev.off()
pdf("embryo_integrated_percent_mito_2020_07_20.pdf",useDingbats = F)
VlnPlot(object = embryo.integrated, features = "percent.mito", group.by = "orig.ident")
dev.off()
pdf("embryo_integrated_ncount_2020_07_20.pdf",useDingbats = F)
VlnPlot(object = embryo.integrated, features = "nCount_RNA", group.by = "orig.ident")
dev.off()
pdf("embryo_integrated_nFeature_2020_07_20.pdf",useDingbats = F)
VlnPlot(object = embryo.integrated, features = "nFeature_RNA", group.by = "orig.ident")
dev.off()

#Identify markers genes
markers.all=FindAllMarkers(embryo.integrated,test.use="roc")
saveRDS(markers.all,"markers_all_2020_07_12.Rdata")
markers.all=markers.all[!grepl("MT-",rownames(markers.all)),]
top10 <- markers.all %>% group_by(cluster) %>% top_n(10, avg_diff)

#Add manually-defined, canonical markers
epi_markers=c("POU5F1","SOX2", "NANOG","TDGF1","BMP4")
TE_markers = c("GATA2","GATA3","TEAD3","TEAD4","TP63","BMP4","ITGA6","HLA-G","LRP5","ELF5","KRT18",
               "CGA","CGB1","SDC1","PSG1","CSH1","INHA","CYP19A1","HSD3B1","KRT7")
primitive_endoderm_markers = c("EOMES", "FRZB", "HNF4A",
                               "PDGFRA","GATA6","GATA4","SOX17","OTX2","BMP6","IGF1")
canonical_markers = unique(c(TE_markers,primitive_endoderm_markers,epi_markers))


genes_to_use = c(top10$gene[1:20],TE_markers,top10$gene[21:30],epi_markers,top10$gene[31:40],primitive_endoderm_markers)

pdf("heatmap_clusters_canonical_with_key_2020_07_12.pdf",width=10,height=8,useDingbats = F)
DoHeatmap(object = embryo.integrated, features = as.character(genes_to_use),group.colors=c("dodgerblue","olivedrab3","firebrick","darkviolet"),size=4)+
  scale_fill_gradientn(colors=c("blue","white","red"))
dev.off()
#---------------------------------------------
#Figure 1c-f
#Read in ensembl annotation for gene name conversion
ensembl_geneid=read.table("ensembl_geneid.txt",header=T,sep="\t")
ensembl_geneid=ensembl_geneid[ensembl_geneid$Ensembl.gene.ID!="",]
genes=ensembl_geneid$Approved.symbol
names(genes)=ensembl_geneid$Ensembl.gene.ID

#Read in data from Stirparo et al and convert to proper gene names
cell_info = read.table("Stirparo_cellinfo.txt",header=T,sep="\t")
embryo_stirparo_raw = read.table("Stirparo_FPKM.txt",header=T,sep="\t")
embryo_stirparo_raw_mat=embryo_stirparo_raw[,3:ncol(embryo_stirparo_raw)]
unique_rownames = embryo_stirparo_raw$Gene
unique_rownames[duplicated(unique_rownames)] = embryo_stirparo_raw$Ensembl.ID[duplicated(unique_rownames)]
embryo_stirparo_raw_mat=apply(embryo_stirparo_raw_mat,2,as.numeric)
rownames(embryo_stirparo_raw_mat)=unique_rownames

#Create Seurat object
embryo_stirparo=CreateSeuratObject(counts = embryo_stirparo_raw_mat)
embryo_stirparo <- NormalizeData(object = embryo_stirparo)
embryo_stirparo <- ScaleData(object = embryo_stirparo)

#Integrate using the anchor method, in effect a batch correction - Note: this needs the own data prior to batch correction
anchors <- FindIntegrationAnchors(object.list = c(batched_seurat_list,embryo_stirparo), dims = 1:30,k.filter=50)
all_genes=union(rownames(embryo.integrated),rownames(embryo_stirparo))
embryo_integrated_stirparo <- IntegrateData(anchorset = anchors, dims = 1:30,features.to.integrate=all_genes)

DefaultAssay(embryo_integrated_stirparo) <- "integrated"

# Run the standard workflow for visualization and clustering
embryo_integrated_stirparo <- ScaleData(embryo_integrated_stirparo, features = all_genes,verbose = FALSE)
embryo_integrated_stirparo <- RunPCA(embryo_integrated_stirparo, npcs = 30, verbose = FALSE)
embryo_integrated_stirparo <- FindNeighbors(object = embryo_integrated_stirparo)
embryo_integrated_stirparo <- FindClusters(embryo_integrated_stirparo,res=0.05)
embryo_integrated_stirparo <- RunTSNE(embryo_integrated_stirparo, reduction = "pca", dims = 1:30)

embryo_integrated_stirparo@meta.data$Age[grepl('embryo',embryo_integrated_stirparo@meta.data$orig.ident)]=11
embryo_integrated_stirparo@meta.data$Age[embryo_integrated_stirparo@meta.data$orig.ident%in%paste0("embryo",c(1,2,7,12,14))]=9

#Annotate the epiblast cells 
embryo_integrated_stirparo@meta.data$Epi="0"
embryo_integrated_stirparo@meta.data$Epi[embryo_integrated_stirparo@meta.data$seurat_clusters==2&embryo_integrated_stirparo@meta.data$Age==11]="Post-Epi"
embryo_integrated_stirparo@meta.data$Epi[embryo_integrated_stirparo@meta.data$seurat_clusters==2&embryo_integrated_stirparo@meta.data$Age==9]="Peri-Epi"
embryo_integrated_stirparo@meta.data$Epi[rownames(embryo_integrated_stirparo@meta.data)%in%cell_info$Sample[cell_info$Stage=="epi"]]="Pre-Epi"
embryo_integrated_stirparo@meta.data$Epi[rownames(embryo_integrated_stirparo@meta.data)%in%cell_info$Sample[cell_info$Stage=="earlyICM"]]="ICM"
epi_cells=rownames(embryo_integrated_stirparo@meta.data)[embryo_integrated_stirparo@meta.data$Epi!="0"]

saveRDS(embryo_integrated_stirparo,"embryo_integrated_stirparo.Rdata")

embryo_epi=subset(embryo_integrated_stirparo,cells = epi_cells)
DimPlot(embryo_epi, reduction = "pca", group.by = "Epi") #Fi

#Compare and plot expression levels

mat = GetAssayData(object = embryo_integrated_stirparo,assay="RNA",slot = "data")[,embryo_integrated_stirparo@meta.data$Epi!="0"]
epi_state=embryo_integrated_stirparo@meta.data$Epi[embryo_integrated_stirparo@meta.data$Epi!="0"]

epi_state=factor(epi_state,levels=c("ICM","Pre-Epi","Peri-Epi","Post-Epi"),ordered = T)
naive_pluripotency = c("PRDM14","TFCP2L1","KLF2","KLF4","KLF5","KLF17","ZFP42","TBX3","ESRRB","DNMT3L","NR0B1","UTF1","SOX15")
primed = c("FGF5","FGF2","OTX2","DNMT3B","NODAL","POU3F1","SOX11","SFRP2","SALL2")
core = c("POU5F1","NANOG","SOX2")

#Fig1d
pdf("naive_pluripotency_epiblast_with_yaxis_2020_07_27_4stages.pdf",width=1.5*length(naive_pluripotency),height=2.2,useDingbats = F)
par(mar=c(2,2,4.1,1),mfrow=c(1,length(naive_pluripotency)))
for (Gene in naive_pluripotency){
  stripchart(mat[Gene,] ~ epi_state, vertical = TRUE,xaxt='n',
             method = "jitter", pch = 20, col = "white",
             xlab="",ylab="",main=Gene,cex.main=2)
  
  grid(col = "lightgray", lty = "dashed")
  stripchart(mat[Gene,] ~ epi_state, vertical = TRUE,xlab="",cex=1.2,
             method = "jitter", add=T,pch = 20, col = c("dodgerblue3","olivedrab3","darkviolet","firebrick"))
  median(mat[Gene,epi_state=="ICM"])
  segments(y0=median(mat[Gene,epi_state=="ICM"]),x0=0.85,x1=1.15,lwd=3)
  segments(y0=median(mat[Gene,epi_state=="Pre-Epi"]),x0=1.85,x1=2.15,lwd=3)
  segments(y0=median(mat[Gene,epi_state=="Peri-Epi"]),x0=2.85,x1=3.15,lwd=3)
  segments(y0=median(mat[Gene,epi_state=="Post-Epi"]),x0=3.85,x1=4.15,lwd=3)
}
dev.off()

#Fig1e
pdf("primed_pluripotency_epiblast_with_yaxis_2020_07_27_4stages.pdf",width=1.5*length(primed),height=2.2,useDingbats = F)
par(mar=c(2,2,4.1,1),mfrow=c(1,length(primed)))
for (Gene in primed){
  if(all(mat[Gene,]==0)){
    stripchart(mat[Gene,] ~ epi_state, vertical = TRUE,xaxt='n',
               method = "jitter", pch = 20, col = "white",ylim=c(0,1),
               xlab="",ylab="",main=Gene,cex.main=2)
  }else{
    stripchart(mat[Gene,] ~ epi_state, vertical = TRUE,xaxt='n',
               method = "jitter", pch = 20, col = "white",
               xlab="",ylab="",main=Gene,cex.main=2)
  }
  
  grid(col = "lightgray", lty = "dashed")
  stripchart(mat[Gene,] ~ epi_state, vertical = TRUE,xlab="",cex=1.2,
             method = "jitter", add=T,pch = 20, col = c("dodgerblue3","olivedrab3","darkviolet","firebrick"))
  median(mat[Gene,epi_state=="ICM"])
  segments(y0=median(mat[Gene,epi_state=="ICM"]),x0=0.85,x1=1.15,lwd=3)
  segments(y0=median(mat[Gene,epi_state=="Pre-Epi"]),x0=1.85,x1=2.15,lwd=3)
  segments(y0=median(mat[Gene,epi_state=="Peri-Epi"]),x0=2.85,x1=3.15,lwd=3)
  segments(y0=median(mat[Gene,epi_state=="Post-Epi"]),x0=3.85,x1=4.15,lwd=3)
  
}
dev.off()

#Fig1f
pdf("core_pluripotency_epiblast_with_yaxis_2020_07_27_4stages.pdf",width=1.5*length(core),height=2.2,useDingbats = F)
par(mar=c(2,2,4.1,1),mfrow=c(1,length(core)))
for (Gene in core){
  stripchart(mat[Gene,] ~ epi_state, vertical = TRUE,xaxt='n',
             method = "jitter", pch = 20, col = "white",
             xlab="",ylab="",main=Gene,cex.main=2)
  
  grid(col = "lightgray", lty = "dashed")
  stripchart(mat[Gene,] ~ epi_state, vertical = TRUE,xlab="",cex=1.2,
             method = "jitter", add=T,pch = 20, col = c("dodgerblue3","olivedrab3","darkviolet","firebrick"))
  median(mat[Gene,epi_state=="ICM"])
  segments(y0=median(mat[Gene,epi_state=="ICM"]),x0=0.85,x1=1.15,lwd=3)
  segments(y0=median(mat[Gene,epi_state=="Pre-Epi"]),x0=1.85,x1=2.15,lwd=3)
  segments(y0=median(mat[Gene,epi_state=="Peri-Epi"]),x0=2.85,x1=3.15,lwd=3)
  segments(y0=median(mat[Gene,epi_state=="Post-Epi"]),x0=3.85,x1=4.15,lwd=3)

}
dev.off()
