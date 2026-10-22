# Load required libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(Seurat)

make_gglot<-function(data_summary,my_title,my_color){
  data_summary$cluster_name<-as.numeric(data_summary$cluster_name)
  my_plot<-ggplot(data_summary, aes(x = as.numeric(cluster_name), y = mean_y)) +
    geom_ribbon(aes(ymin = mean_y - sd_y, ymax = mean_y + sd_y), fill = my_color, alpha = 0.5) +
    #geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", color = "black", size = 1.2))+
    geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2, span = 1.3) +  # Smooth line
    geom_point(size = 3) +  # Add points at the mean values
    labs(title = my_title, 
         x = "Cluster", 
         y = "Mean Score ± SD") +
    scale_x_continuous(
      breaks = 3:5,  # Define positions for the ticks
      labels = c("Tnd.Ent", "UFC", "MFC")) +
    theme_classic()+
    theme(
    plot.title = element_text(size = 16,hjust = 0.5,face="bold"),    # Title font size
    axis.title.x = element_text(size = 0), # X-axis label font size
    axis.title.y = element_text(size = 16), # Y-axis label font size
    axis.text.x = element_text(size = 16),  # X-axis tick label font size
    axis.text.y = element_text(size = 16)   # Y-axis tick label font size
  )
 
  return(my_plot)
  
}

#seurat object was first called, 
seurat.file="/home/labs/zelzerlab/Collaboration/Ron_scRNA_Nov22/scRNAseq_2/analysis/integrated_2experiments/Reanalysis_Nov23/relevantCells_ctrl_allGenes_22Jan25.rds"
seurat.obj.f=readRDS(seurat.file)
tendon.genes=c("Tnc","Scx","Col1a1","Tnmd","Thbs4","Mkx","Col3a1","Thbs2","Lum","Bgn","Col5a1",
               "Col5a2","Col5a3","Fmod","Timp2","Col8a1","Dcn")
cartilage.genes=c("Wwp2","Sox5","Sox6","Col11a1","Col9a1","Acan",
                  "Col2a1","Nog","Fgfr1","Sox9","Runx2","Comp","Matn3",
                  "Cnmd","Bmp6","Ucma","Creb3l2","Enpp1","Ext1","Npr2","Pth1r")#"Bmp1b",

mineralization.genes=c("Phex","Gja1","Tnfrsf11b","Col1a1","Sparc",
"Runx2","Spp1","Bglap","Alpl","Col10a1","Ibsp","Dlx5","Col1a2")

seurat.obj.f<-AddModuleScore(seurat.obj.f,features = list(tendon.genes),name="tendon_score")
seurat.obj.f<-AddModuleScore(seurat.obj.f,features = list(cartilage.genes),name="cartilage_score")
seurat.obj.f<-AddModuleScore(seurat.obj.f,features = list(mineralization.genes),name="mineralization_score")


seurat.obj.f<-subset(seurat.obj.f,idents=c("Tnd.Ent","UFC","MFC"))
my_metadata=data.frame(cluster_name=seurat.obj.f$biol.names,
                       tendon.score=as.numeric(seurat.obj.f$tendon_score1),
                       cartilage.score=as.numeric(seurat.obj.f$cartilage_score1),
                       mineralization.score=as.numeric(seurat.obj.f$mineralization_score1))
my_metadata$cluster_name<-as.character(my_metadata$cluster_name)
my_metadata$cluster_name <- ifelse(my_metadata$cluster_name == "Tnd.Ent", "3", my_metadata$cluster_name)
my_metadata$cluster_name <- ifelse(my_metadata$cluster_name == "UFC", "4", my_metadata$cluster_name)
my_metadata$cluster_name <- ifelse(my_metadata$cluster_name == "MFC", "5", my_metadata$cluster_name)
####################

# Summarize data: calculate mean and standard deviation for each group
data_summary <- my_metadata %>%
  group_by(cluster_name) %>%
  summarize(mean_y = mean(tendon.score), sd_y = sd(tendon.score))
p1<-make_gglot(data_summary,"Tendon Score","#DDC7CF")
print(p1)

###########
data_summary <- my_metadata %>%
  group_by(cluster_name) %>%
  summarize(mean_y = mean(cartilage.score), sd_y = sd(cartilage.score))
p2<-make_gglot(data_summary,"Cartilage Score","#EE82F8")
print(p2)
# Plot smooth line with standard deviation as a shaded area

#############
data_summary <- my_metadata %>%
  group_by(cluster_name) %>%
  summarize(mean_y = mean(mineralization.score), sd_y = sd(mineralization.score))
p3<-make_gglot(data_summary,"Mineralization Score","#224760")
print(p3)
# plot all together
grid.arrange(p1, p2, p3, ncol = 1)
