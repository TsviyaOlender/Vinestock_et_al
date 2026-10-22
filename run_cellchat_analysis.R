library(CellChat)
library(patchwork)
library(Seurat)
library(scCustomize)
#https://github.com/jinworks/CellChat?tab=readme-ov-file
##############################################################
seurat.obj<-readRDS("seurat_object_11Sep25.rds")
DimPlot(seurat.obj)
meta<-seurat.obj@meta.data
cellchat <- createCellChat(object = seurat.obj, group.by = "biol.names", assay = "SCT")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "biol.names") # set "labels" as default cell identity

##############################################################
#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# use Secreted Signaling, includes WNT signalling
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 
cellchat@DB <- CellChatDB.use
##############################################################
#Preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
##############################################################
#Part II: Inference of cell-cell communication network
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat,file="cellchat_object.rds")
##############################################################
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
##############################################################
#Part III: Visualization of cell-cell communication network
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
##############################################################
pathways.show <- c("WNT") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
par(mar = c(1, 1, 5, 1))
netVisual_aggregate(cellchat, signaling = pathways.show,show.legend = TRUE)

# Chord diagram
netVisual_aggregate(cellchat, signaling = "WNT", layout = "hierarchy")
netVisual_aggregate(
  cellchat, signaling = "WNT", layout = "hierarchy",
  vertex.receiver = 3:7
)
# Compute the network centrality scores
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
##############################################################
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
netAnalysis_contribution(cellchat, signaling = pathways.show)
##############################################################
# one ligand receptor
pairLR.WNT <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.WNT[1,] # show one ligand-receptor pair
# Hierarchy plot
# a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show)
#> [[1]]
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
##############################################################
saveRDS(cellchat,file="cellchat_object.rds")
