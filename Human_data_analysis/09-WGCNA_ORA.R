library("clusterProfiler")
library("org.Hs.eg.db")
library("ReactomePA")
library("WGCNA")
library("DEP")
library("tidyverse")
library("enrichplot")



load(file = "./Results/WGCNA_net.RData")
ncx_genes <- read.csv(file = "BrainSpan_background_genes/BrainSpan_NCX_background_ID.txt", header = FALSE)

geneModules <- net$colors
genes <- names(geneModules)
ids <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
genes <- as.data.frame(genes)
colnames(genes) <- "SYMBOL"
geneIDs <- dplyr::left_join(genes, ids, )
geneIDsModules <-  cbind(geneIDs, geneModules)
colnames(geneIDsModules)[3] <- "geneModules"
write.csv(geneIDsModules, file = "./Results/WGCNA_ORA/geneIDsModules.csv")
geneIDsModules2 <- geneIDsModules[geneIDsModules$geneModules != "grey",]
write.csv(geneIDsModules2, file = "./Results/WGCNA_ORA/geneIDsModules2.csv")

#use the union of brain expressed genes and module genes as universe background
ncx_genes2 <- bitr(ncx_genes$V1, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
sum(!(geneIDsModules$ENTREZID %in% ncx_genes2$ENTREZID))
universe <- union(geneIDsModules$ENTREZID,ncx_genes2$ENTREZID)
universe <- na.omit(universe)
#save background gene list
background_ENTREZ <- as.data.frame(universe)
colnames(background_ENTREZ) <- "ENTREZID"
backgroundIDs <- bitr(universe, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
write.csv(backgroundIDs, file = "./Results/WGCNA_ORA/backgroundIDs.csv")

#load MSigDB data
c2cp_gmt <- read.gmt("MSigDB/c2.cp.v7.1.symbols.gmt")
c2cp_gmt$term <- tolower(c2cp_gmt$term)

#list gene names in each module
modules <- unique(geneModules)
modules <- modules[!modules == "grey"]
module_name_for_ORA <- list()
for (i in modules) {
  tmp <- geneIDsModules[geneIDsModules$geneModules == i,]$SYMBOL
  module_name_for_ORA <- c(module_name_for_ORA, list(tmp))
}

names(module_name_for_ORA) <- modules
#note that for purpose above, you can also use split() to do it

#gene names in background

res_MSigDB <- compareCluster(geneCluster = module_name_for_ORA,
                             fun = "enricher", 
                             pvalueCutoff = 0.05, 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05,
                             minGSSize = 10, 
                             maxGSSize = 500,
                             universe = backgroundIDs$SYMBOL,
                             TERM2GENE=c2cp_gmt)

save(module_name_for_ORA, file ="Results/module_name_for_ORA.RData")
# save(res_MSigDB, res_Reactome, res_KEGG, res_GO, res_GO2, file = "Results/WGCNA_ORA.RData")
write.csv(res_MSigDB@compareClusterResult, file = "Results/WGCNA_ORA/res_MSigDB_c2cp_ORA.csv")

