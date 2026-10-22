library(viridis)
library(pheatmap)
library(Seurat)
library(RColorBrewer)


seurat.obj<-readRDS("relevantCells_ctrl_allGenes.rds")
#---------------------------------------------------------------------------------
# read list of variable genes for the calculation
#---------------------------------------------------------------------------------

variable.ordered<-read.table(file="variableDE_sorted120324.txt",header=T)
#---------------------------------------------------------------------------------
# generate the heatmap
# 1. extract the expression relevant genes from the matrix
# 2. generate a matrix with the average expression for cluster
# 3. calculate the pearson correlation
# 4. plot
#---------------------------------------------------------------------------------
# step 1
subset_seurat <- subset(seurat.obj, features = variable.ordered$gene)#$V1 in unordered
# step2
expression_matrix <- GetAssayData(subset_seurat, slot = "data",assay="integrated")
expr.mat.t<-as.matrix(expression_matrix)
# step 3
correlation_matrix <- cor(t(expr.mat.t), method = "pearson")
# step4
# make nice colors
color_palette <- brewer.pal(9, "PuBuGn")
reordered_matrix <- correlation_matrix[variable.ordered$gene, variable.ordered$gene]
# to get rid from gene names set show_rownames=F
data.order=pheatmap(correlation_matrix, cluster_rows = F, cluster_cols = F, breaks =NA,
                    legend=T, show_rownames=F, show_colnames=F, border_color=NA,fontsize=8,
                    color=color_palette,
                    scale="none")
