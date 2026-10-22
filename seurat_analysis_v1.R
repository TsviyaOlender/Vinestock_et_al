library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(rio)
library(grid)
library(gridExtra)
library(scCustomize)
##
source("function_enrichR.R")
save.dir = "/home/labs/olenderlab/lvzvia/zelzer/Ron_scRNA_Nov22/scRNAseq_2/analysis"
output_seurat = paste(save.dir,"/seurat_obj_jan23.rds",sep="")
#---------------------------------------------------------------------------------
# Upload experiment 1
#---------------------------------------------------------------------------------
data.folder="/home/labs/olenderlab/lvzvia/zelzer/Ron_scRNA_Nov22/scRNAseq_1/cellranger"
experiments = c("Runx1_Control","Runx1_Mutant")
experimentNames = c("Runx1.Control.1","Runx1.Mutant.1")
samples.num=2
####
# load matrices and create seurat objects
for (i in 1:samples.num){
  experimentLocation = paste(data.folder,"/",experiments[i],"/outs/filtered_feature_bc_matrix",sep="")
  experiment.data <- Read10X(data.dir =experimentLocation);
  seurat.obj<- CreateSeuratObject(counts  = experiment.data, project = experimentNames[i], min.cells = 3, min.features = 0)
  seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj , pattern = "^mt-")
  assign(experimentNames[i],seurat.obj)
  
}
#---------------------------------------------------------------------------------
# Upload experiment 2
#---------------------------------------------------------------------------------
#Upload experiment 2
data.folder="/home/labs/olenderlab/lvzvia/zelzer/Ron_scRNA_Nov22/scRNAseq_2/cellranger"
experiments = c("Runx1_Control","Runx1_Mutant")
experimentNames = c("Runx1.Control.2","Runx1.Mutant.2")
samples.num=2
# load matrices and create seurat objects
for (i in 1:samples.num){
  experimentLocation = paste(data.folder,"/",experiments[i],"/outs/filtered_feature_bc_matrix",sep="")
  experiment.data <- Read10X(data.dir =experimentLocation);
  seurat.obj<- CreateSeuratObject(counts  = experiment.data, project = experimentNames[i], min.cells = 3, min.features = 0)
  seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj , pattern = "^mt-")
  assign(experimentNames[i],seurat.obj)
  
}
#---------------------------------------------------------------------------------
# Cell filtration, merge objects
#---------------------------------------------------------------------------------
#filter
Runx1.Control.f.1 <- subset(Runx1.Control.1, subset = percent.mt <8 & nCount_RNA > 200 & nCount_RNA < 80000)
Runx1.Mutant.f.1 <- subset(Runx1.Mutant.1, subset = percent.mt < 8 & nCount_RNA > 200 & nCount_RNA < 80000)
Runx1.Control.f.2 <- subset(Runx1.Control.2, subset = percent.mt <8 & nCount_RNA > 200 & nCount_RNA < 90000)
Runx1.Mutant.f.2 <- subset(Runx1.Mutant.2, subset = percent.mt < 7 & nCount_RNA > 200 & nCount_RNA < 75000)

experimentNames = c("Runx1.Control.1","Runx1.Mutant.1","Runx1.Control.2","Runx1.Mutant.2")
seurat.obj.f <- merge(x=Runx1.Control.f.1,y=list(Runx1.Mutant.f.1,Runx1.Control.f.2,Runx1.Mutant.f.2),
                      add.cell.ids = experimentNames)
#---------------------------------------------------------------------------------
# Normalization
# SCT, using v2
#---------------------------------------------------------------------------------
seurat.obj.f <- SCTransform(seurat.obj.f, verbose = FALSE,vst.flavor = "v2")
# after deep review of the data we have decided to apply RPCA integration
#---------------------------------------------------------------------------------
#https://satijalab.org/seurat/articles/integration_rpca.html
ifnb.list <- SplitObject(seurat.obj.f, split.by = "orig.ident")
# normalize and identify variable features for each sample independently
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
ifnb.list <- lapply(X = ifnb.list, FUN = RunPCA, features = features)
# perform integration
experiment.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, 
                                             reduction = "rpca")
# this command creates an 'integrated' data assay
seurat.obj.corr <- IntegrateData(anchorset = experiment.anchors, normalization.method = "SCT", dims = 1:30)
#
DefaultAssay(seurat.obj.corr) <- "integrated"
#---------------------------------------------------------------------------------
# dimensionality reduction and clustering
#---------------------------------------------------------------------------------
seurat.obj.corr <- RunPCA(seurat.obj.corr, verbose = FALSE,assay="integrated")
ElbowPlot(seurat.obj.corr,ndims=50)#
pca.num=30 # Num of PCA
res.num = 0.6 # resolution
####
seurat.obj.corr <- RunUMAP(seurat.obj.corr, reduction = "pca", dims = 1:30)
# Visualization
p1<-DimPlot(seurat.obj.corr, reduction = "umap", group.by = "orig.ident")
p2<-DimPlot(seurat.obj.corr, reduction = "umap", split.by = "orig.ident",label=F)
p1+p2
#---------------------------------------------------------------------------------
# clustering
#---------------------------------------------------------------------------------
seurat.obj.corr <- FindNeighbors(seurat.obj.corr, dims = 1:pca.num)
seurat.obj.corr <- FindClusters(seurat.obj.corr, resolution = res.num)
seurat.obj.corr <- RunUMAP(seurat.obj.corr, reduction = "pca", dims = 1:pca.num, verbose = FALSE)
p1<-DimPlot(seurat.obj.corr, reduction = "umap",label=FALSE,label.size = 6,group.by="orig.ident")
p2<-DimPlot(seurat.obj.corr, reduction = "umap",label=TRUE,label.size = 4)
p1+p2

UMAP.cluster.file.bySample<-paste(save.dir,"/","UMAP_clusters_",pca.num,"PCA_res_",res.num,"_bySample.png",sep="")
UMAP.clusters.file = paste(save.dir,"/","UMAP_clusters_",pca.num,"PCA_res_",res.num,".png",sep="")
ggplot2::ggsave(filename = UMAP.cluster.file.bySample, plot = p1) 
ggplot2::ggsave(filename = UMAP.clusters.file, plot = p2)
#---------------------------------------------------------------------------------
# Find Markers
#---------------------------------------------------------------------------------
seurat.obj.corr.markers <- FindAllMarkers(seurat.obj.corr, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2,assay="SCT")
#TOP10
seurat.obj.corr.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
#top2
seurat.obj.corr.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) -> top2
HP2<-DoHeatmap(seurat.obj.corr, features = c(top2$gene)) + NoLegend()
# save. file names
markers.file = paste(save.dir,"/","markers_",pca.num,"PCA_res_",res.num,".csv",sep="")
top2.file = paste(save.dir,"/","top2..markers_",pca.num,"PCA_res_",res.num,".csv",sep="")
top10.file = paste(save.dir,"/","top10.markers_",pca.num,"PCA_res_",res.num,".csv",sep="")
top2.HM.file = paste(save.dir,"/","top2.heatmap_",pca.num,"PCA_res_",res.num,".png",sep="")
# write
write.csv(seurat.obj.corr.markers,file=markers.file)
write.csv(top2,file=top2.file)
write.csv(top10,file=top10.file)
ggplot2::ggsave(filename = top2.HM.file, plot = HP2)
#---------------------------------------------------------------------------------
# automatic annotation using PanglaoDB
#---------------------------------------------------------------------------------
seurat.obj.corr<-main_enrichR(markers.file,save.dir,seurat.obj.corr,0)
#verify that the identities are now the
DimPlot_scCustom(seurat.obj.f)

#---------------------------------------------------------------------------------
# save object
#---------------------------------------------------------------------------------
session_file_name=paste(save.dir,"/session_info.txt",sep="")
sessionInfo() %>% capture.output(file=session_file_name)
saveRDS(seurat.obj.corr,file=output_seurat)
