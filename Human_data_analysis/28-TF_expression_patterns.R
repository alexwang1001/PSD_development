library(tidyverse)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(WGCNA)
library(circlize)
library(seriation)


#import pseudobulk data
pseudobulk_res <- read.csv("single_cell_adult_human/pseudobulk_res_subclass_human_10x.csv", row.names = 1)
pseudobulk_res_cpm <- sweep(pseudobulk_res*1000000, MARGIN = 2, STATS = colSums(pseudobulk_res), FUN = "/")
assay <- log(pseudobulk_res_cpm+1,2)
coldata <- tibble(sample = colnames(pseudobulk_res_cpm))
coldata <- separate(coldata, 1 ,into = c("class", "subclass"), sep = "_")

#blue module
blue_ChEA3 <- read.table(file = "ChEA3/blue_Integrated_meanRank.tsv", header = T, sep = "\t")
assay_filt_blue <- assay[rownames(assay) %in% blue_ChEA3$TF[1:100],]
write.csv(assay_filt_blue, "Results/ChEA3/Blue module TF expression log CPM.csv")
#remove genes with too many missing values
index <- goodSamplesGenes(t(assay_filt_blue))
assay_filt_blue2 <- assay_filt_blue[index$goodGenes,]
assay_filt_blue2_scaled <- t(scale(t(assay_filt_blue2)))
write.csv(assay_filt_blue2_scaled, "Results/ChEA3/Blue module TF expression scaled.csv")

#yellow module
yellow_ChEA3 <- read.table(file = "ChEA3/yellow_Integrated_meanRank.tsv", header = T, sep = "\t")
assay_filt_yellow <- assay[rownames(assay) %in% yellow_ChEA3$TF[1:100],]
write.csv(assay_filt_yellow, "Results/ChEA3/Yellow module TF expression log CPM.csv")
#remove genes with too many missing values
index <- goodSamplesGenes(t(assay_filt_yellow))
assay_filt_yellow2 <- assay_filt_yellow[index$goodGenes,]
assay_filt_yellow2_scaled <- t(scale(t(assay_filt_yellow2)))
write.csv(assay_filt_yellow2_scaled, "Results/ChEA3/Yellow module TF expression scaled.csv")