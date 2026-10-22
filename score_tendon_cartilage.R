library(Seurat)
library(dplyr)
save.dir = "/home/labs/olenderlab/lvzvia/zelzer/Ron_scRNA_Nov22/scRNAseq_2/analysis"
output_seurat = paste(save.dir,"/seurat_obj_jan23.rds",sep="")
#---------------------------------------------------------------------------------
# upload object
#---------------------------------------------------------------------------------
seurat_obj <-readRDS(output_seurat)
#---------------------------------------------------------------------------------
# define lists
#---------------------------------------------------------------------------------
tendon.genes=c("Tnc","Scx","Col1a1","Tnmd","Thbs4","Mkx","Col3a1","Thbs2","Lum","Bgn","Col5a1",
               "Col5a2","Col5a3","Fmod","Timp2","Col8a1","Dcn")

cartilage.genes=c("Wwp2","Sox5","Sox6","Col11a1","Col9a1","Acan",
                  "Col2a1","Bmp1b","Nog","Fgfr1","Sox9","Runx2","Comp","Matn3",
                  "Cnmd","Bmp6","Ucma","Creb3l2","Enpp1","Ext1","Npr2","Pth1r")
#---------------------------------------------------------------------------------
# Calculate module scores
#---------------------------------------------------------------------------------
seurat_obj <- AddModuleScore(
  object = seurat_obj,
  features = list(tendon.present, cartilage.present),
  name = "ModuleScore"
)

# Rename columns
seurat_obj$TendonScore <- seurat_obj$ModuleScore1
seurat_obj$CartilageScore <- seurat_obj$ModuleScore2
#---------------------------------------------------------------------------------
# plot
#---------------------------------------------------------------------------------
p1<-FeaturePlot_scCustom(output_seurat,features = c("TendonScore","CartilageScore"))
file.name<-paste(save.dir,"/tendon_cartilage_scores.jpg",sep="")
ggsave(file.name, plot = p1, width = 20, height = 18, dpi = 150)
#---------------------------------------------------------------------------------
# save object with scores
#---------------------------------------------------------------------------------
saveRDS(seurat_obj,file=output_seurat)