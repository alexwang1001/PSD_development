library(tidyverse)
library(SummarizedExperiment)

#Human and Macaque RNA-seq data
#import RNA-seq data
dat2 <- read.delim(file = "Transcriptome/human_and_macaque_RNA_seq/nhp_development_RPKM_rmTechRep.txt", row.names = 1)
#import colData
coldata2 <- read.csv(file = "Transcriptome/human_and_macaque_RNA_seq/codata.csv")
#generate rowData
rowdata2 <- as.data.frame(rownames(dat2)) %>% separate(col = 1, into = c("ENSEMBL", "SYMBOL"), sep = "\\|") %>% separate(col = 1, into = c("ENSEMBL", "Version"), sep = "\\.")
rowdata2$Version <- NULL
#generate assay data
assay2 <- dat2[,coldata2$sample]
b <- colSums(assay2)
assay2 <- sweep(x = assay2*1000000, MARGIN = 2, STATS = b, FUN = "/")
assay2 <- log(assay2+1,2)
colnames(assay2) <- coldata2$sample
rownames(assay2) <- rowdata2$ENSEMBL

#generate summarized experiment
se2 <- SummarizedExperiment(assays = assay2, colData = coldata2, 
                           rowData = rowdata2)

NCX <- c("A1C", "DFC", "IPC", "ITC", "M1C", "MFC", "OFC", "S1C", "STC", "V1C", "VFC")
se2_NCX <- se2[,se2$region %in% NCX]
se2_subset <- se2_NCX[,(se2_NCX$species == "macaque" & se2_NCX$days > 80) | (se2_NCX$species == "human" & se2_NCX$days < 7937 & se2_NCX$days > 111)]

saveRDS(se2_subset, file = "Results/Human_Macaque_NCX_subset_transcriptomics_logTPM.rds")
