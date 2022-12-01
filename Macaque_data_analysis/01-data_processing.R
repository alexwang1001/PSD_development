library("biomaRt")
library("DEP")
library("tidyverse")
library("SummarizedExperiment")
library("AnnotationDbi")
library("org.Mmu.eg.db")
library("limma")
library("purrr")
library("vsn")
library("clusterProfiler")

#import data
data_iBAQ <- read.csv(file = "Clean_data/iBAQ.csv", header = T)

#ID conversion
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
head(datasets)
ensembl = useDataset("mmulatta_gene_ensembl",mart=ensembl)
attributes = listAttributes(ensembl)
filters <- listFilters(ensembl)
ids <- getBM(attributes=c('uniprotsptrembl', 'external_gene_name', 'entrezgene_id', 'ensembl_gene_id'), 
      filters = 'uniprotsptrembl', 
      values = data_iBAQ$Protein.IDs, 
      mart = ensembl)
colnames(ids) <- c("Protein.IDs", "SYMBOL", "ENTREZID", "ENSEMBL")

#remove multiple duplicated uniprot IDs
ids2 <- dplyr::distinct(ids, Protein.IDs, .keep_all = T)
data_iBAQ2 <- data_iBAQ %>% dplyr::left_join(ids2)
#which ids are unmapped?
data_iBAQ_unmapped <- data_iBAQ2[!(data_iBAQ2$Protein.IDs %in% ids2$Protein.IDs),]
#save files for manually annotate the unmapped IDs
write.csv(data_iBAQ2, file = "Clean_Data/iBAQ2.csv", row.names = FALSE)

#import the updated annotations
data_iBAQ3 <- read.csv(file = "Clean_Data/iBAQ3.csv", header = T)
#remove proteins/genes without ENTREZID
data_iBAQ4 <- data_iBAQ3 %>% drop_na(ENTREZID)
#make sure SYMBOLs are consistent with org.Mmu.eg.db
data_iBAQ4$SYMBOL <- mapIds(org.Mmu.eg.db, keys=as.character(data_iBAQ4$ENTREZID), column="SYMBOL", keytype="ENTREZID", multiVals="first")
colnames(data_iBAQ4)[2] <- "Gene.names"

###########################################################################
#which genes are from mitochondria? Use the human mitocarta 3.0 and find its homologs in macaque
mito <- read.csv(file = "./mitogenes/Mitogenes3_human.csv", header = TRUE)
mito_human_id <- mito$HumanGeneID
mito_human_id <- mito_human_id[-1136]

human <-  useMart(biomart="ensembl", dataset = "hsapiens_gene_ensembl", verbose = TRUE, host = "dec2021.archive.ensembl.org")
macaque <-  useMart(biomart="ensembl", dataset = "mmulatta_gene_ensembl", verbose = TRUE, host = "dec2021.archive.ensembl.org")
res <- getLDS(attributes = c("external_gene_name","entrezgene_id"),
              filters = "entrezgene_id", values = mito_human_id, mart = human,
              attributesL = c("external_gene_name","entrezgene_id"), martL = macaque, uniqueRows = TRUE)
colnames(res) <- c("Human_SYMBOL", "Human_ENTREZID", "Macaque_SYMBOL", "Macaque_ENTREZID")
sum(data_iBAQ4$ENTREZID %in% res$Macaque_ENTREZID)

#remove all mitocondria related proteins
data_iBAQ5 <- data_iBAQ4[!(data_iBAQ4$ENTREZID %in% res$Macaque_ENTREZID),]

#####################################################################
a <- colSums(data_iBAQ5[,5:25])
data_iBAQ5[,5:25] <- sweep(x = data_iBAQ5[,5:25]*1000000, MARGIN = 2, STATS = a, FUN = "/")
colSums(data_iBAQ5[,5:25])

data4 <- data_iBAQ5
data_unique <- make_unique(data4, "Gene.names", "Protein.IDs", delim = ";")
head(data_unique)
write.csv(data_unique, file = "Clean_Data/data_unique.csv", row.names = FALSE)

#make a SummerizedExperiment
expdesign <- read.csv(file = "./Sample_Info/experimental_design.csv", header = TRUE)

proteins_unique <- data_unique
columns <- c(5:25)
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
