library(tidyverse)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(circlize)
library(seriation)

#subclass data
pseudobulk_res <- read.csv("single_cell_human_development_Dmitry/pseudobulk_res_type_human_development_Dmitry.csv", row.names = 1)
geneIDsModules <- read.csv(file = "Results/WGCNA_ORA/geneIDsModules.csv", row.names = 1)
pseudobulk_res_PSD <- pseudobulk_res[rownames(pseudobulk_res) %in% geneIDsModules$SYMBOL,]
pseudobulk_res_PSD_cpm <- sweep(pseudobulk_res_PSD*1000000, MARGIN = 2, STATS = colSums(pseudobulk_res_PSD), FUN = "/")

coldata <- tibble(sample = colnames(pseudobulk_res_PSD_cpm))
coldata <- separate(coldata, 1 ,into = c("class", "subclass", "stage"), sep = "_")
coldata <- coldata %>% mutate(class = case_when(
  class == "EN" ~ "Glutamatergic",
  class == "IN" ~ "GABAergic"
))
coldata <- coldata %>% mutate(subclass = case_when(
  subclass == "EN.IT" ~ "IT",
  subclass == "EN.non.IT" ~ "non.IT",
  subclass == "IN.CGE" ~ "CGE",
  subclass == "IN.MGE" ~ "MGE"
))
coldata <- coldata %>% mutate(stage = case_when(
  stage == "0.1.years" ~ "0-1 years",
  stage == "1.2.years" ~ "1-2 years",
  stage == "10.20.years" ~ "10-20 years",
  stage == "2.4.years" ~ "2-4 years",
  stage == "2nd.trimester" ~ "2nd trimester",
  stage == "3rd.trimester" ~ "3rd trimester",
  stage == "4.10.years" ~ "4-10 years",
  TRUE ~ stage
))

rowdata <- tibble(SYMBOL = rownames(pseudobulk_res_PSD_cpm))
rowdata <- left_join(rowdata, geneIDsModules, by = "SYMBOL")

#the result from AggregateExpression is count that has not been log-transformed
assay <- log(pseudobulk_res_PSD_cpm+1,2)
#generate a summarized experiment

se <- SummarizedExperiment(assays = as.matrix(assay), colData = coldata, 
                           rowData = rowdata)

saveRDS(se, file = "Results/pseudo_bulk_type_age_human_development_Dmitry.rds")

assay2 <- cbind(rowdata, assay)
write.csv(assay2,"Results/single_cell_human_development_Dmitry/log cpm counts.csv", row.names = F)

#scale each gene
se2 <- se
assay(se2) <- t(scale(t(assay(se2))))

saveRDS(se, file = "Results/pseudo_bulk_type_age_human_development_Dmitry_scaled.rds")

yellow <- assay(se2[rowData(se2)$geneModules == "yellow",])
write.csv(yellow, "Results/single_cell_human_development_Dmitry/yellow.csv")

blue <- assay(se2[rowData(se2)$geneModules == "blue",])
write.csv(blue, "Results/single_cell_human_development_Dmitry/blue.csv")
