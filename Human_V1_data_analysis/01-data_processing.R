library("DEP")
library("dplyr")
library("tibble")
library("tidyr")
library("SummarizedExperiment")
library("limma")
library("purrr")
library("vsn")
library("clusterProfiler")
library("org.Hs.eg.db")
library("stringr")
library("ComplexHeatmap")

PFC_data <- read.csv("data/data_unique_non-postmortem_PFC.csv")
PFC_data <- PFC_data[,c(1:2,73:76,3:72)]
BA17_data <- read.csv("data/riBAQ_one_uniprot.csv")
BA17_data$Gene.names <- NULL
data <- left_join(PFC_data, BA17_data, by = "Protein.IDs")
#assign 0 to all NA values
data[,77:98][is.na(data[,77:98])] <- 0

#re-normalize
a <- colSums(data[,77:98])
data[,77:98] <- sweep(x = data[,77:98]*1000000, MARGIN = 2, STATS = a, FUN = "/")
colSums(data[,7:98])
#save data
write.csv(data, file = "data/data_unique_all.csv", row.names = FALSE)

#make a SummerizedExperiment of final data
expdesign <- read.csv(file = "Sample_Info/experimental_design_all.csv", header = TRUE)

proteins_unique <- data
columns <- c(7:98)
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
# data_se <- se
# colData(data_se)
# head(assay(data_se))
# rowData(data_se)
# dim(data_se)
saveRDS(se, file = "./Results/data_all.rds")

