library(tidyverse)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(circlize)
library(seriation)
library(WGCNA)
library(viridis)


#subclass data
pseudobulk_res <- read.csv("single_cell_adult_human/pseudobulk_res_subclass_human_10x.csv", row.names = 1)
geneIDsModules <- read.csv(file = "Results/WGCNA_ORA/geneIDsModules.csv", row.names = 1)
pseudobulk_res_PSD <- pseudobulk_res[rownames(pseudobulk_res) %in% geneIDsModules$SYMBOL,]
pseudobulk_res_PSD_cpm <- sweep(pseudobulk_res_PSD*1000000, MARGIN = 2, STATS = colSums(pseudobulk_res_PSD), FUN = "/")

coldata <- tibble(sample = colnames(pseudobulk_res_PSD_cpm))
coldata <- separate(coldata, 1 ,into = c("class", "subclass"), sep = "_")


rowdata <- tibble(SYMBOL = rownames(pseudobulk_res_PSD_cpm))
rowdata <- left_join(rowdata, geneIDsModules, by = "SYMBOL")

#the result from AggregateExpression is count that has not been log-transformed
assay <- log(pseudobulk_res_PSD_cpm+1,2)
#generate a summarized experiment

se <- SummarizedExperiment(assays = as.matrix(assay), colData = coldata, 
                           rowData = rowdata)
colData(se)

saveRDS(se, file = "Results/pseudo_bulk_subclass_adult_human.rds")

assay2 <- cbind(rowdata, assay)
write.csv(assay2,"Results/single_cell_adult_human/log cpm counts.csv", row.names = F)


#scale each gene
se2 <- se
assay(se2) <- t(scale(t(assay(se2))))

saveRDS(se2, file = "Results/pseudo_bulk_subclass_adult_human_scaled.rds")

yellow <- assay(se2[rowData(se2)$geneModules == "yellow",])
write.csv(yellow, "Results/single_cell_adult_human/yellow.csv")

blue <- assay(se2[rowData(se2)$geneModules == "blue",])
write.csv(blue, "Results/single_cell_adult_human/blue.csv")
