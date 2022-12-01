library("tidyverse")
library("org.Hs.eg.db")

#import BioGRID PPI data
human_entrez <- read.table(file = "BioGRID_from_Kaifang/BioGRID_Human.txt", sep = "\t")
mouse_entrez <- read.table(file = "BioGRID_from_Kaifang/BioGRID_Mouse.txt", sep = "\t")
PPI_entrez <- rbind(human_entrez, mouse_entrez)

PPI_entrez$fromNode <- mapIds(org.Hs.eg.db, keys=as.character(PPI_entrez$V1), column="SYMBOL", keytype="ENTREZID", multiVals="first")
PPI_entrez$toNode <- mapIds(org.Hs.eg.db, keys=as.character(PPI_entrez$V2), column="SYMBOL", keytype="ENTREZID", multiVals="first")

PPI_entrez <- PPI_entrez[!is.na(PPI_entrez$fromNode) & !is.na(PPI_entrez$toNode),]
#make switch protein1 and protein2 so that the matching with co-expression file can be in any direction
PPI <- PPI_entrez[,3:4]
PPI2 <- PPI_entrez[,4:3]
colnames(PPI2) <- c("fromNode","toNode")
PPI <- rbind(PPI, PPI2)
PPI_unique <- unique(PPI)
PPI_unique$interaction <- "pp"

#import co-expression edges
edges <- read.table(file = "Results/PPI_network/WGCNA_cytoscape_export/CytoscapeInput-edges.txt", header = TRUE, sep = "\t")
edges2 <- left_join(edges, PPI_unique, by = c("fromNode", "toNode"))
nodes <- read.table(file = "Results/PPI_network/WGCNA_cytoscape_export/CytoscapeInput-nodes.txt", header = TRUE, sep = "\t")
nodes <- nodes[,c(1,3)]
colnames(nodes) <- c("fromNode","fromNode_module")
edges3 <- left_join(edges2, nodes, by = "fromNode")
colnames(nodes) <- c("toNode","toNode_module")
edges4 <- left_join(edges3, nodes, by = "toNode")

#only include top 10% co-expressed edges as co-expression
edges_sorted <- edges4[order(edges4$weight, decreasing = TRUE),]
edges_all <- edges_sorted %>% filter(!is.na(interaction))
write.table(edges_all, file = "Results/PPI_network/WGCNA_cytoscape_export/edges_all.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


edges_0.1 <- head(edges_sorted, n=155673)
edges_0.1_all <- edges_0.1 %>% filter(!is.na(interaction))
write.table(edges_0.1_all, file = "Results/PPI_network/WGCNA_cytoscape_export/edges_0.1_all.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

modules <- c("brown", "blue", "turquoise", "yellow")
for (module in modules) {
  edges_0.1_module <- edges_0.1 %>% filter(fromNode_module == module & toNode_module == module) %>% filter(!is.na(interaction))
  write.table(edges_0.1_module, file = paste0("Results/PPI_network/WGCNA_cytoscape_export/edges_0.1_", module, ".txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}
