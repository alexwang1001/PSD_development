library("DEP")
library("dplyr")
library("tibble")
library("tidyr")
library("SummarizedExperiment")
library("AnnotationDbi")
library("org.Mm.eg.db")
library("limma")
library("purrr")
library("vsn")
library("clusterProfiler")

data_iBAQ <- read.csv(file = "Clean_data/iBAQ.csv", header = T)
ids  <- bitr(data_iBAQ$Protein.IDs, fromType='UNIPROT', toType='ENTREZID', OrgDb = "org.Mm.eg.db")
colnames(ids) <- c("Protein.IDs", "ENTREZID")

#remove multiple duplicated uniprot IDs and keep the first one
ids2 <- dplyr::distinct(ids, Protein.IDs, .keep_all = T)
data_iBAQ2 <- data_iBAQ %>% dplyr::left_join(ids2)
write.csv(data_iBAQ2, file = "Clean_Data/iBAQ2.csv", row.names = FALSE)
#two proteins are products of the same gene in mouse, O54931/Akap2/677884 is removed
#manually annotated unmapped IDs
data_iBAQ3 <- read.csv(file = "Clean_Data/iBAQ3.csv", header = T)
data_iBAQ4 <- data_iBAQ3 %>% drop_na(ENTREZID)
data_iBAQ4$SYMBOL<- mapIds(org.Mm.eg.db, keys=as.character(data_iBAQ4$ENTREZID), column="SYMBOL", keytype="ENTREZID", multiVals="first")
data_iBAQ4$Gene.names <- data_iBAQ4$SYMBOL
data_iBAQ4 <- subset(data_iBAQ4, select = -SYMBOL)
data_iBAQ4$ENSEMBL<- mapIds(org.Mm.eg.db, keys=as.character(data_iBAQ4$ENTREZID), column="ENSEMBL", keytype="ENTREZID", multiVals="first")

###########################################################################
#remove all mitochondria related proteins
mito <- read.csv(file = "./mitogenes/Mitogenes.csv", header = TRUE)
data_iBAQ5 <- data_iBAQ4[!(data_iBAQ4$ENTREZID %in% mito$MouseOrthologGeneID),]
#re-normalize data
a <- colSums(data_iBAQ5[,3:22])
data_iBAQ5[,3:22] <- sweep(x = data_iBAQ5[,3:22]*1000000, MARGIN = 2, STATS = a, FUN = "/")
colSums(data_iBAQ5[,3:22])

data4 <- data_iBAQ5
data_unique <- make_unique(data4, "Gene.names", "Protein.IDs", delim = ";")
head(data_unique)
write.csv(data_unique, file = "Clean_Data/data_unique.csv", row.names = FALSE)

#make a SummerizedExperiment
expdesign <- read.csv(file = "./Sample_Info/experimental_design.csv", header = TRUE)

proteins_unique <- data_unique
columns <- c(3:22)
any(!c("name", "ID") %in% colnames(proteins_unique))
any(!c("label", "condition", "replicate") %in% colnames(expdesign))
any(!apply(proteins_unique[, columns], 2, is.numeric))
rownames(proteins_unique) <- proteins_unique$name
raw <- proteins_unique[, columns]
raw[raw == 0] <- NA
raw <- log2(raw)
expdesign <- mutate(expdesign, ID = label)
rownames(expdesign) <- expdesign$ID
duplicated(expdesign$ID)
raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]
row_data <- proteins_unique[, -columns]
rownames(row_data) <- row_data$name
se <- SummarizedExperiment(assays = as.matrix(raw), colData = expdesign, 
                           rowData = row_data)
data_se <- se
colData(data_se)
head(assay(data_se))
rowData(data_se)
dim(data_se)
save(data_se, file = "./Results/data.RData")
