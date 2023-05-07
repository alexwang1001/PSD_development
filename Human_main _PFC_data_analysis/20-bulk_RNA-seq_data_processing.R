library("tidyverse")
library("SummarizedExperiment")
library("biomaRt")
library("WGCNA")

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


#import developmental RNA-seq data
se <- readRDS("human_and_macaque_RNA_seq/Human_Macaque_NCX_subset_transcriptomics_logTPM.rds")
rowdata <- as.data.frame(rowData(se))
#import human PSD protein list
geneIDsModules <- read.csv("Results/WGCNA_ORA/geneIDsModules.csv", row.names = 1)

#convert ENTREZID to ENSEMBL ID
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
human <-  useMart(biomart="ensembl", dataset = "hsapiens_gene_ensembl", verbose = TRUE, host = "dec2021.archive.ensembl.org")
ids <-  getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name'), 
              filters = 'entrezgene_id', 
              values = geneIDsModules$ENTREZID, 
              mart = human)
colnames(ids) <- c("ENSEMBL", "ENTREZID", "GENENAME")
ids_filt <- ids[ids$ENSEMBL %in% rowdata$ENSEMBL,]
ids_filt[duplicated(ids_filt$ENTREZID),]
#remove duplicated ENTREZID after conversion
ids_filt2 <- ids_filt[!(ids_filt$ENSEMBL %in% c("ENSG00000243902","ENSG00000015479")), ]
geneIDsMpdules_filt <- left_join(ids_filt2, geneIDsModules)
geneIDsMpdules_filt <- geneIDsMpdules_filt[, -3]
rownames(geneIDsMpdules_filt) <- geneIDsMpdules_filt$SYMBOL

#only keep PSD genes and human samples in the RNA-seq dataset
se_filt <- se[geneIDsMpdules_filt$ENSEMBL,se$species == "human"]
#filter out genes with too many missing values or zero variance
goodsamplesgenes <- goodSamplesGenes(t(assay(se_filt)))
se_filt2 <- se_filt[goodsamplesgenes$goodGenes,goodsamplesgenes$goodSamples]
geneIDsMpdules_filt2 <- geneIDsMpdules_filt[goodsamplesgenes$goodGenes,]
#recalculate TPM values with PSD genes
assay2 <- assay(se_filt2)
assay2 <- 2^assay2-1
b <- colSums(assay2)
assay2 <- sweep(x = assay2*1000000, MARGIN = 2, STATS = b, FUN = "/")
assay2 <- log(assay2+1,2)
assay2 <- as.matrix(assay2)
rownames(assay2) <- geneIDsMpdules_filt2$SYMBOL

coldata2 <- data.frame(colData(se_filt2))
rowdata2 <- geneIDsMpdules_filt2

se2 <- SummarizedExperiment(assays=assay2, rowData =rowdata2, colData=coldata2)
saveRDS(se2, "Results/Human_NCX_subset_transcriptomics_PSD_logTPM.rds")


