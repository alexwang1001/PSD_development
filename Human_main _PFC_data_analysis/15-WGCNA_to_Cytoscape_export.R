library("WGCNA")
library("DEP")
library("SummarizedExperiment")
library("limma")
library("tidyverse")
library("vsn")

load(file = "./Results/dep_corrected.RData")
datExpr <- t(assay(dep_corrected))

load(file = "Results/WGCNA_net.RData")

TOM = TOMsimilarityFromExpr(datExpr,
                             networkType = "signed",
                             power = 14)

moduleColors <- net$colors
#select modules
modules = c("brown", "blue", "turquoise", "yellow")

#select genes in modules
module = "turquoise"
genes = colnames(datExpr)
inModule = is.finite(match(moduleColors, module))
modgenes = genes[inModule]
modTOM = TOM[inModule,inModule]
dimnames(modTOM) = list(modgenes, modgenes)

#which threshold to choose as coexpression threshold for the network?
dim(modTOM)
#percentage of edges pass threshold
(sum(modTOM>0.1)-dim(modTOM)[1])/(dim(modTOM)[1]^2-dim(modTOM)[1])
#mean edge number
(sum(modTOM>0.1)-dim(modTOM)[1])/2/dim(modTOM)[1]


for (module in modules) {
  #select genes in modules
  genes = colnames(datExpr)
  inModule = is.finite(match(moduleColors, module))
  modgenes = genes[inModule]
  modTOM = TOM[inModule,inModule]
  dimnames(modTOM) = list(modgenes, modgenes)
  # Export the network into edge and node list files so that Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modgenes,
                                 altNodeNames = modgenes,
                                 nodeAttr = moduleColors[inModule])
  
}

# Export the whole network into edge and node list files so that Cytoscape can read
cyt2 = exportNetworkToCytoscape(TOM,
                               edgeFile = paste("CytoscapeInput-edges.txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes.txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = genes,
                               altNodeNames = genes,
                               nodeAttr = moduleColors)
