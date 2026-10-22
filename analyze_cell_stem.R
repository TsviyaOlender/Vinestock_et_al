library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(rio)
library(grid)
library(gridExtra)
source("function_enrichR.R")
save.dir = "/home/labs/olenderlab/lvzvia/zelzer/Ron_scRNA_Nov22/scRNAseq_2/analysis/integrated_2experiments/Reanalysis_Nov23/Literature_comparison/cell_stem_cell_2022/"
output_seurat = paste(save.dir,"/seurat_lit_cell_stem.rds",sep="")
#---------------------------------------------------------------------------------
# Upload experiment 
#---------------------------------------------------------------------------------
data.folder="/home/labs/olenderlab/lvzvia/zelzer/Ron_scRNA_Nov22/scRNAseq_2/analysis/integrated_2experiments/Reanalysis_Nov23/Literature_comparison/cell_stem_cell_2022/data"
experiments = c("p11","p18","p56")
experimentNames = c("p11","p18","p56")
samples.num=3
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
# merge and filter
#---------------------------------------------------------------------------------
seurat.obj <- merge(x=p11,y=list(p18,p56),
                    add.cell.ids = experimentNames)
seurat.obj.f<-subset(seurat.obj,percent.mt < 20 & nFeature_RNA > 200 & nCount_RNA < 45000)
VlnPlot(seurat.obj.f  , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,pt.size=0.5)
table(seurat.obj.f$orig.ident)
#---------------------------------------------------------------------------------
# batch correction
#---------------------------------------------------------------------------------
#https://satijalab.org/seurat/articles/integration_rpca.html
ifnb.list <- SplitObject(seurat.obj.f, split.by = "orig.ident")
# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")
# this command creates an 'integrated' data assay
seurat.obj.f <- IntegrateData(anchorset = immune.anchors)
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(seurat.obj.f) <- "integrated"

# Run the standard workflow for visualization and clustering
seurat.obj.f <- ScaleData(seurat.obj.f, verbose = FALSE)
seurat.obj.f <- RunPCA(seurat.obj.f, verbose = FALSE)
seurat.obj.f <- RunUMAP(seurat.obj.f, reduction = "pca", dims = 1:25)
seurat.obj.f <- FindNeighbors(seurat.obj.f, reduction = "pca", dims = 1:25)
seurat.obj.f <- FindClusters(seurat.obj.f, resolution = 0.6)

# Visualization
DimPlot(seurat.obj.f)
markers<-FindAllMarkers(seurat.obj.f,min.pct = 0.4,assay="SCT")
markers<-subset(markers,markers$p_val_adj < 0.05 & markers$avg_log2FC > 0.2)
markers.file = "home/labs/olenderlab/lvzvia/zelzer/Ron_scRNA_Nov22/scRNAseq_2/analysis/integrated_2experiments/Reanalysis_Nov23/Literature_comparison/cell_stem_cell_2022/solution_pca25_res0.6/2_markers/cell_markers_25PCA_res_0.6.csv"
write.csv(markers,file=markers.file)
#---------------------------------------------------------------------------------
# automatic annotation using PanglaoDB
#---------------------------------------------------------------------------------
seurat.obj.f<-main_enrichR(markers.file,save.dir,seurat.obj.f,0)
#---------------------------------------------------------------------------------
# after enrichR annotation
# we extract fibroblasts and chrondrocyte cells
#---------------------------------------------------------------------------------
Idents(seurat.obj.f)<-seurat.obj.f$panglaoDB
sub<-subset(seurat.obj.f,idents = c("Fibroblasts","Chondrocytes"))
#---------------------------------------------------------------------------------
# Recluster, with RNA normalization
#---------------------------------------------------------------------------------
sub <- SCTransform(sub, verbose = FALSE,return.only.var.genes = TRUE)
DefaultAssay(sub)<-"SCT"
sub <- RunPCA(sub, verbose = FALSE)
ElbowPlot(sub,ndims=50)
sub <- RunUMAP(sub, reduction = "pca", dims = 1:25)
sub <- FindNeighbors(sub, reduction = "pca", dims = 1:25)
sub <- FindClusters(sub, resolution = 0.8)
DimPlot_scCustom(sub)
DimPlot_scCustom(sub,label=T,group.by = "seurat_clusters",label.size = 6)
#---------------------------------------------------------------------------------
# Matching the clusters with the manuscript names
#---------------------------------------------------------------------------------
FeaturePlot_scCustom(sub,features = c("Bglap","Dmp1"))
FeaturePlot_scCustom(sub,features = c("Scx", "Col1a1"))
paper_genes=c("Ly6a","Cd34","Pdgfra","Cd44","Cebpd","Scx",
              "Tnmd","Col1a1","Col1a2","Egr1","Dcn","Fmod","Sox9","Acan",
              "Col3a1","Col2a1","Bgn","Pth1r","Alpl","Nfatc1","Bglap",
              "Spp1","Ibsp","Dmp1","Mef2c","Foxa3")
supp_genes=c("Cebpd","Myc","Creb5","Bcl3","Nfia","Egr1","Tcf4","Creb3l1","Sox9",
             "Runx2","Mef2a","Mef2c","Foxa3","Nfatc1")
DoHeatmap(sub,features = paper_genes)
seurat_clusters_backup<-sub$seurat_clusters
#---------------------------------------------------------------------------------
# clusters renaming
#---------------------------------------------------------------------------------
Idents(sub)<-sub$seurat_clusters
sub <- RenameIdents(object = sub, "0" = "enthesoblasts",
                    "1"="enthesis_progenitors",
                    "2"="enthesoblasts",
                    "3"="enthesis_progenitors",
                    "4"="mineralizing_chondrocytes",
                    "5"="pre-enthosoblasts",
                    "6"="mineralizing_chondrocytes",
                    "7"="tenoblasts",
                    "8"="osteocytes")

DimPlot(sub,label=T)
levels(as.factor(sub@active.ident))
desired_order=c("enthesis_progenitors","pre-enthosoblasts","enthesoblasts","mineralizing_chondrocytes","tenoblasts","osteocytes","9","10")
sub@active.ident<- factor(sub@active.ident , levels = desired_order)
DoHeatmap(sub,features = paper_genes)
table(sub@active.ident,sub$orig.ident)
VlnPlot_scCustom(sub,features = c("Ibsp","Spp1"),pt.size=1)
DimPlot_scCustom(sub,group.by = "seurat_clusters")
#---------------------------------------------------------------------------------
# save
#---------------------------------------------------------------------------------
sub$biol_names_0.4<-sub@active.ident
saveRDS(sub,"fibroblasts_clusters_integrated_seuratV4.rds")
#---------------------------------------------------------------------------------
# extract p18
#---------------------------------------------------------------------------------
Idents(sub)<-sub$orig.ident
p18<-subset(sub,idents="p18")
Idents(p18)<-p18$biol_names_0.4
p18<-subset(p18,idents=c("tenoblasts","osteocytes","10"),invert=T)
# I generate a new object only with the counts
DimPlot_scCustom(p18,pt.size=1)
count_matrix <- GetAssayData(p18, slot = "counts")
metadata <- p18@meta.data
p18<- CreateSeuratObject(counts = count_matrix, meta.data = metadata)
#---------------------------------------------------------------------------------
# I normalize p18 with SCTnormalization to fit with Ron's data
#---------------------------------------------------------------------------------
p18<-SCTransform(p18, verbose = FALSE,return.only.var.genes = TRUE)
DefaultAssay(p18)<-'SCT'
p18<- RunPCA(p18)
ElbowPlot(p18,ndims=50)
p18<- RunUMAP(p18, reduction = "pca", dims = 1:20)
p18 <- FindNeighbors(p18, reduction = "pca", dims = 1:20)
p18 <- FindClusters(p18, resolution = 0.4)
DimPlot_scCustom(p18,pt.size=1)
#---------------------------------------------------------------------------------
# label transfer
#---------------------------------------------------------------------------------
ref.obj=readRDS("../../relevantCells_ctrl_allGenes_22Jan25.rds")
DimPlot_scCustom(ref.obj)
Idents(ref.obj)<-ref.obj$biol.names
query.object<-p18
experiment.anchors <- FindTransferAnchors(reference = ref.obj, query = query.object, dims = 1:15,
                                          reference.reduction = "pca",normalization.method = "SCT")#,recompute.residuals=FALSE)
predictions <- TransferData(anchorset = experiment.anchors, refdata = ref.obj@meta.data[["biol.names"]], dims = 1:15)
query.object <- AddMetaData(query.object, metadata = predictions)
Idents(query.object) <- query.object$predicted.id
query.object@reductions[["umap"]]@cell.embeddings[,1]<-query.object@reductions[["umap"]]@cell.embeddings[,1]*-1
query.object@reductions[["umap"]]@cell.embeddings[,2]<-query.object@reductions[["umap"]]@cell.embeddings[,2]*-1
DimPlot_scCustom(query.object,pt.size=1,label=F)
FeaturePlot_scCustom(query.object,features = "prediction.score.max",pt.size=1)
####################################################################
saveRDS(query.object, file="p18_sct.rds")
##################################################################
query.object=readRDS("../../relevantCells_ctrl_allGenes_22Jan25.rds")
p1<-DimPlot_scCustom(ref.obj)
p2<-DimPlot_scCustom(query.object)
p1|p2
Idents(ref.obj)<-ref.obj$biol_names_0.4
experiment.anchors <- FindTransferAnchors(reference = ref.obj, query = query.object, dims = 1:15,
                                          reference.reduction = "pca",normalization.method = "SCT")#,recompute.residuals=FALSE)
predictions <- TransferData(anchorset = experiment.anchors, refdata = ref.obj@meta.data[["biol_names_0.4"]], dims = 1:15)
query.object <- AddMetaData(query.object, metadata = predictions)
Idents(query.object) <- query.object$predicted.id
DimPlot_scCustom(query.object,pt.size=1,label=F)
FeaturePlot_scCustom(query.object,features = "prediction.score.max",pt.size=1)