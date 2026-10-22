library(enrichR)
library(gridExtra)
library(openxlsx)
library(ggplot2)
library(stringr)
#################################################################################
main_enrichR<- function(marker_file,save.folder,the.object,perCond){
  #perCond= if perCond==1, it is a test between 2 conditions. selecting markers differently
  markerGenes<-read.csv(file=marker_file,header=TRUE)
  # The number of the last cluster
  markerGenes$cluster<-as.factor(markerGenes$cluster)
  clusters<-levels(markerGenes$cluster)
  clusters<-as.integer(clusters)
  clusterT = max(clusters)
  # subsetting
  # only positive markers
  # select genes with FC > 0.5 and min cell with expression > 0.2
  min_FC <-0.2 # equals to fold change of 2
  min_p_adj = 0.05
  min.pt1 = 0.2
 
  if(perCond==1){
    markerGenes<-subset(markerGenes,abs(markerGenes$avg_log2FC) > min_FC & markerGenes$p_val_adj < min_p_adj)
  }else{
    markerGenes<-subset(markerGenes,markerGenes$avg_log2FC > abs(min_FC) & markerGenes$p_val_adj < min_p_adj & 
                          markerGenes$pct.1 > min.pt1)
  }
  test=nrow(markerGenes)
  
  result_go_enrichR<-data.frame()
  result_pathways_enrichR<-data.frame()
  result_cellType_enrichR<-data.frame()
  top_PanglaoDB <-data.frame()
  # create folders for outputs
  run_folders =c('go_enrichR','pathways_enrichR','cellType_enrichR') 
  for (i in 1:3){
    current_folder=paste(save.folder,"/",run_folders[i],sep="")
    ifelse(!dir.exists(current_folder), dir.create(current_folder, showWarnings = FALSE), FALSE)
    run_folders[i]=current_folder
    
  }
  #run all clusters
  for (i in 0:clusterT){
    cluster_name<-paste(i, sep="")
    selected<-subset(markerGenes,markerGenes$cluster==cluster_name)
    checkSize <-nrow(selected)
    send_to_annotation=0
    
  
    if(checkSize > 20){
      send_to_annotation=1
      tempRun<-run_go_enrichR(selected,i,checkSize,run_folders[1])
      result_go_enrichR<-rbind(result_go_enrichR,tempRun)
      #
      tempRun<-run_pathway_enrichR(selected,i,checkSize,run_folders[2])
      result_pathways_enrichR<-rbind(result_pathways_enrichR,tempRun)
      #
      tempRun<-run_cellType_enrichR(selected,i,checkSize,run_folders[3])
      result_cellType_enrichR<-rbind(result_cellType_enrichR,tempRun)
      
    }
    if(send_to_annotation==1){
      panglaoDB<-subset(tempRun,tempRun$database =="PanglaoDB_Augmented_2021")
      panglaoDB<-panglaoDB[order(panglaoDB$Adjusted.P.value), ]
      top_PanglaoDB<-rbind(top_PanglaoDB,panglaoDB[1,])
    }else{
      top_PanglaoDB<-rbind(top_PanglaoDB,cluster_name)
    }
    
  }
  plot_results(run_folders[1],result_go_enrichR,c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021",
                                                  "GO_Biological_Process_2021"))
  plot_results(run_folders[1],result_pathways_enrichR,c("WikiPathway_2021_Human", "KEGG_2021_Human"))
  # save results
  write.xlsx(result_go_enrichR,file=paste0(run_folders[1],"/EnrichR_GO_terms_enrichment.xlsx"),overwrite=TRUE ) 
  write.xlsx(result_pathways_enrichR,file=paste0(run_folders[2],"/EnrichR_pathways_enrichment.xlsx"),overwrite=TRUE ) 
  write.xlsx(result_cellType_enrichR,file=paste0(run_folders[3],"/EnrichR_cellType_enrichment.xlsx"),overwrite=TRUE )
  #######################################################################################################
  # add top celltype to object
  #####################################################################
 
  clusters<-levels(the.object@active.ident)
  new_idents <- setNames(top_PanglaoDB$Term, as.character(clusters))
   # Rename the identities
  the.object <- RenameIdents(the.object, new_idents)
  the.object$panglaoDB=the.object@active.ident

  p1<-DimPlot_scCustom(the.object,label=F)+ggtitle("PanglaoDB annotation")
  p2<-DimPlot_scCustom(the.object,label=T,group.by = "seurat_clusters")
  print(p1|p2)
  annotation.UMAP.file=paste(save.folder,"/UMAP_by_enrichR_annotation.jpeg",sep="")
  ggplot2::ggsave(filename = annotation.UMAP.file, plot = p1|p2,units ="cm",width=30,height=20)
  top.PanglaoDB.file=paste(save.folder,"/top_PanglaoDB.csv",sep="")
  write.csv(top_PanglaoDB,file=top.PanglaoDB.file)
  Idents(the.object)<-the.object$seurat_clusters
  return(the.object)
}

#################################################################################
# FUNCTIONS GO terms
#################################################################################
# run enrichR
# GO terms enrichment
run_go_enrichR<- function(df,clusterN,num_of_genes,outFolder) {
  setEnrichrSite("Enrichr")
  dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021",
           "GO_Biological_Process_2021")
  enriched <- enrichr(df$gene, dbs)
  # plot
  p1<-plotEnrich(enriched$GO_Molecular_Function_2021, showTerms = 10,  numChar = 40, y = "Count", orderBy = "pvalue",title="GO_Molecular_Function_2021")
  p2<-plotEnrich(enriched$GO_Cellular_Component_2021, showTerms = 10,  numChar = 40, y = "Count", orderBy = "pvalue",title="GO_Cellular_Component_2021")
  p3<-plotEnrich(enriched$GO_Biological_Process_2021,showTerms = 10,   numChar = 40, y = "Count", orderBy = "pvalue",title="GO_Biological_Process_2021")
  figN<-paste(outFolder,"/","GO_terms_enrichment",clusterN,".jpeg",sep="")
  ctitle<-paste("GO-terms enrichment for cluster",clusterN," with ",num_of_genes," genes",sep="")
  png(filename = figN,width = 1000,
      height = 1000)
  grid.arrange(p1,p2,p3)
  dev.off()
  #
  GO1<-subset(enriched[["GO_Molecular_Function_2021"]],enriched[["GO_Molecular_Function_2021"]]$Adjusted.P.value < 0.05)
  GO2<-subset(enriched[["GO_Cellular_Component_2021"]],enriched[["GO_Cellular_Component_2021"]]$Adjusted.P.value < 0.05)
  GO3<-subset(enriched[["GO_Biological_Process_2021"]],enriched[["GO_Biological_Process_2021"]]$Adjusted.P.value < 0.05)
  # add database name to output file
  GO1$database<-rep("GO_Molecular_Function_2021",length(GO1$Adjusted.P.value))
  GO2$database<-rep("GO_Cellular_Component_2021",length(GO2$Adjusted.P.value))
  GO3$database<-rep("GO_Biological_Process_2021",length(GO3$Adjusted.P.value))
  #add cluster
  all_GO_db<-rbind(GO1,GO2,GO3)
  all_GO_db$cluster<-rep(clusterN,length(all_GO_db$P.value))
  return(all_GO_db)
  
}
#################################################################################
# function pathway
#################################################################################
# pathway enrichment
run_pathway_enrichR<- function(df,clusterN,num_of_genes,outFolder) {
  dbs <- c("WikiPathway_2021_Human", "KEGG_2021_Human")
  enriched <- enrichr(df$gene, dbs)
  # plot
  p1<-plotEnrich(enriched$WikiPathway_2021_Human,  numChar = 40, y = "Count", orderBy = "pvalue",title="WikiPathway_2021_Human")
  p2<-plotEnrich(enriched$KEGG_2021_Human,  numChar = 40, y = "Count", orderBy = "pvalue",title="KEGG_2021_Human")
  
  
   figN<-paste(outFolder,"/","cluster",clusterN,"_pathways_enrichment",".jpeg",sep="")
  ctitle<-paste("pathway enrichment for cluster",clusterN," with ",num_of_genes," genes",sep="")
  png(filename = figN,width = 1000,
      height = 1000)
  grid.arrange(p1,p2,top=ctitle)
  dev.off()
  ####
  
  GO1<-subset(enriched[["WikiPathway_2021_Human"]],enriched[["WikiPathway_2021_Human"]]$Adjusted.P.value < 0.05)
  GO2<-subset(enriched[["KEGG_2021_Human"]],enriched[["KEGG_2021_Human"]]$Adjusted.P.value < 0.05)

  # add database name to output file
  GO1$database<-rep("WikiPathways_2021_Human",length(GO1$Adjusted.P.value))
  GO2$database<-rep("KEGG_2021_Human",length(GO2$Adjusted.P.value))
  #add cluster
  all_GO_db<-rbind(GO1,GO2)
  all_GO_db$cluster<-rep(clusterN,length(all_GO_db$P.value))
  
  return(all_GO_db)
 
}
#################################################################################
# cell type enrichment
#################################################################################
run_cellType_enrichR<- function(df,clusterN,num_of_genes,outFolder) {
  
  dbs <- c("PanglaoDB_Augmented_2021", "CellMarker_Augmented_2021")
  enriched <- enrichr(df$gene, dbs)
  # plot
  if(length(df$gene) > 0){
  p1<-plotEnrich(enriched$PanglaoDB_Augmented_2021, showTerms = 10,  numChar = 40, y = "Count", orderBy = "pvalue",title="PanglaoDB_Augmented_2021")
  p2<-plotEnrich(enriched$CellMarker_Augmented_2021,showTerms = 10,   numChar = 40, y = "Count", orderBy = "pvalue",title="CellMarker_Augmented_2021")
  }
  
  figN<-paste(outFolder,"/","cluster",clusterN,"_cellType",".jpeg",sep="")
  ctitle<-paste("cellType enrichment for cluster",clusterN," with ",num_of_genes," genes",sep="")
  png(filename = figN,width = 1000,
      height = 1000)
  grid.arrange(p1,p2,top=ctitle)
  dev.off()
  #
  # save results to file
  GO1<-subset(enriched[["PanglaoDB_Augmented_2021"]],enriched[["PanglaoDB_Augmented_2021"]]$Adjusted.P.value < 0.05)
  GO2<-subset(enriched[["CellMarker_Augmented_2021"]],enriched[["CellMarker_Augmented_2021"]]$Adjusted.P.value < 0.05)
  
  # add database name to output file
  GO1$database<-rep("PanglaoDB_Augmented_2021",length(GO1$Adjusted.P.value))
  GO2$database<-rep("CellMarker_Augmented_2021",length(GO2$Adjusted.P.value))
  #add cluster
  all_GO_db<-rbind(GO1,GO2)
  all_GO_db$cluster<-rep(clusterN,length(all_GO_db$P.value))
  return(all_GO_db)

}
#########################################################################################
singleList_enrichR<- function(genes_to_run){
  setEnrichrSite("Enrichr")
  dbs <- c("GO_Molecular_Function_2021")
  enriched <- enrichr(genes_to_run, dbs)
  # plot
  p1<-plotEnrich(enriched$GO_Molecular_Function_2021, showTerms = 10,  numChar = 40, y = "Count", orderBy = "pvalue",title="GO_Biological_Process_2021")
  
  return(p1)
}
###########################################################################################
plot_results<-function(folder_to_save,df,databases){
#  (run_folders[1],result_go_enrichR,c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021",
#"GO_Biological_Process_2021"))
  
  df$cluster<-as.factor(df$cluster)
  row.count=length(df$Overlap)
  for (i in 1:row.count){
    df$gene_count[i]<-str_count(df$Genes[i], ";")+1
    # remove GO number
    term_len<-nchar(df$Term[i])
    term_end=term_len-13
    df$short_name[i]<- substr(df$Term[i], 1, term_end)
  }
  
  databases_num=length(databases)
  for(i in 1:databases_num){
    BP=subset(df,df$database==databases[i])
    BP <- BP %>%group_by(cluster) %>%arrange(P.value) %>%  # Arrange by p-value in ascending order
      slice_head(n = 5)  # Select the top 10 terms
    p1<-ggplot(BP, aes(x = cluster, y = short_name)) +
      geom_point(aes(size = gene_count, color = -log10(Adjusted.P.value))) +  # Using -log10(p_value) to represent p-value
      scale_size(range = c(3, 10)) +  # Adjust the size range
      scale_color_gradient(low = "red", high = "blue") +  # Adjust color scale as needed
      theme_minimal() +  # Clean theme
      labs(
        x = "Cluster",
        y = "GO Term",
        size = "Gene Count",
        color = "-log10(p-value)"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme_classic()  # Rotate x-axis labels if needed
      out_file=paste(folder_to_save,"_",databases[i],"_summary_plot.jpg",sep="")
      ggsave(out_file, plot = p1)
  }

}