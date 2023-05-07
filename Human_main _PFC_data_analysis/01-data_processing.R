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

data_riBAQ <- read.csv(file = "./Clean_Data/riBAQ_one_uniprot.csv", header = T)
colnames(data_riBAQ)
ids  <- bitr_kegg(data_riBAQ$Protein.IDs, fromType='uniprot', toType='ncbi-geneid', organism='hsa')
colnames(ids) <- c("Protein.IDs", "ENTREZID")

data_riBAQ_unmapped <- data_riBAQ[!(data_riBAQ$Protein.IDs %in% ids$Protein.IDs),]

#remove multiple uniprot IDs that link to the same ENTREZID, only retain one
ids2 <- dplyr::distinct(ids, Protein.IDs, .keep_all = T)
data_riBAQ2 <- data_riBAQ %>% dplyr::left_join(ids2)

write.csv(data_riBAQ2, file = "Clean_Data/riBAQ2.csv")

#manually annotate the unmapped IDs
data_riBAQ3 <- read.csv(file = "Clean_Data/riBAQ3_contamination_removed.csv", header = T)
data_riBAQ4 <- data_riBAQ3 %>% drop_na(ENTREZID)
class(data_riBAQ4$ENTREZID)
data_riBAQ4$ENTREZID <- as.character(data_riBAQ4$ENTREZID)

ids2  <- bitr (data_riBAQ4$ENTREZID, fromType='ENTREZID', toType='SYMBOL', OrgDb="org.Hs.eg.db")

data_riBAQ5 <- data_riBAQ4 %>% dplyr::left_join(ids2, by = "ENTREZID")
data_riBAQ5$Gene.names <- data_riBAQ5$SYMBOL
data_riBAQ5 <- subset(data_riBAQ5, select = -SYMBOL)

#import proteins localized to mitochondria
mito <- read.csv(file = "./Mitogenes/Mitogenes3_human.csv", header = TRUE)

#remove all mitocondria related genes
data_riBAQ6 <- data_riBAQ5[!(data_riBAQ5$ENTREZID %in% mito$HumanGeneID),]
any(duplicated(data_riBAQ6$Gene.names))
data_unique <- make_unique(data_riBAQ6, "Gene.names", "Protein.IDs", delim = ";")
data_unique$name %>% duplicated() %>% any()
write.csv(data_unique, file = "Clean_Data/data_unique.csv",row.names = FALSE)

#calculate the number of proteins identified at each developmental stage
#assign all match between runs values to 0
identification_method <- read.csv(file = "Clean_Data/identification_method_contamination_removed.csv")
identification_method_unique <- identification_method %>% dplyr::semi_join(data_unique, by = "Protein.IDs")
matched_idx <- identification_method_unique[,2:71] == "By matching"
data_unique_MSMSonly <- data_unique
data_unique_MSMSonly[,3:72][matched_idx] <- 0

a <- colSums(data_unique_MSMSonly[,3:72])
data_unique_MSMSonly[,3:72] <- sweep(x = data_unique_MSMSonly[,3:72]*1000000, MARGIN = 2, STATS = a, FUN = "/")
# colSums(data_unique_MSMSonly[,3:72])

#create a SummerizedExperiment
expdesign <- read.csv(file = "./Sample_Info/experimental_design.csv", header = TRUE)
proteins_unique <- data_unique_MSMSonly
columns <- c(3:72)
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
data_se_MSMSonly <- se
colData(data_se_MSMSonly)
head(assay(data_se_MSMSonly))
rowData(data_se_MSMSonly)
save(data_se_MSMSonly, file = "./Results/data_MSMSonly.RData")


#proteins identified in each stage, only keep those that are identified in at least half of all samples in that stage
plot_missval(data_se_MSMSonly)
bin_data <- assay(data_se_MSMSonly)
idx <- is.na(assay(data_se_MSMSonly))
bin_data[!idx] <- 1
bin_data[idx] <- 0

keep <- bin_data %>%
  data.frame() %>%
  rownames_to_column() %>%
  gather(ID, value, -rowname) %>% 
  left_join(data.frame(colData(data_se_MSMSonly)), by = "ID") %>% 
  group_by(rowname, condition) %>%
  summarize(miss_val = n() - sum(value)) %>% 
  filter(miss_val <= 3) %>% 
  ungroup()

keep2 <- keep %>% spread(key = condition, value = miss_val) %>% as.data.frame()
#Overlap with previous studies in adult tissue
#2011 study
Bayes <- read.csv(file = "comparion_with_previous_PSD_list/bayes.csv", header = FALSE)
#2017 study
Roy_ACCNUM <- read.csv(file = "comparion_with_previous_PSD_list/roy2.csv", header = FALSE)
acc <- as.data.frame(str_split(Roy_ACCNUM$V1,"\\|", simplify = TRUE))
acc$V4 <- str_split(acc$V4,"\\.", simplify = TRUE)[,1]
x <- AnnotationDbi::select(org.Hs.eg.db, keys = acc$V4, column = c("SYMBOL", "ENTREZID"), keytype = "ACCNUM")
acc[,5:6] <- x[,2:3]
acc <- as.data.frame(acc)
write.csv(acc, file = "comparion_with_previous_PSD_list/Roy_sp_draft.csv", row.names = FALSE)
#manually annotate unmatched gene IDs
Roy <- read.csv(file = "comparion_with_previous_PSD_list/Roy_sp.csv")

#remove all mitochondria related genes
mito <- read.csv(file = "./Mitogenes/Mitogenes3_human.csv", header = TRUE)
Bayes_no_mito <-Bayes$V1[!(Bayes$V1 %in% mito$HumanGeneID)]
Bayes_no_mito <- unique(Bayes_no_mito)
Roy_no_mito <- as.character(Roy$ENTREZID[!(Roy$ENTREZID %in% mito$HumanGeneID)])
Roy_no_mito <- unique(Roy_no_mito)

#compare proteins identified in group Year18_22 to those by Bayes et al. and Roy et al.
mydata <- keep2$rowname[!(is.na(keep2$Year18_22))]
mydata <- mapIds(org.Hs.eg.db, mydata, "ENTREZID", "SYMBOL")

m <- make_comb_mat(Bayes = Bayes_no_mito,
                   Roy = Roy_no_mito,
                   Wang = mydata)
m
UpSet(m, comb_order = order(comb_size(m), decreasing = TRUE), width = unit(3, "in"), height = unit(1.2, "in"))
dev.off()


#only keep proteins identified from non-postmortem samples
#group1: proteins identified from GW18-23 non-postmortem samples
keep3 <- bin_data %>%
  data.frame() %>%
  rownames_to_column() %>%
  gather(ID, value, -rowname) %>% 
  left_join(data.frame(colData(data_se_MSMSonly)), by = "ID") %>% 
  group_by(rowname, condition) %>%
  summarize(miss_val = n() - sum(value)) %>% 
  filter(miss_val <= 7) %>% 
  ungroup()


GW18_23 <- keep3 %>% filter(condition %in% c("GW18_19", "GW22_23"))
GW18_23_proteins <- unique(GW18_23$rowname)
GW18_23_IDs <- mapIds(org.Hs.eg.db,keys = GW18_23_proteins, column = "ENTREZID", keytype = "SYMBOL")

#group2: proteins identified from Bayes 2011 from neurosurgery samples
Bayes <- read.csv(file = "comparion_with_previous_PSD_list/bayes.csv", header = FALSE)
protein_list <- unique(c(Bayes$V1,GW18_23_IDs))

#keep proteins identified from non-postmortem samples (group1+group2)
data_unique2 <- data_unique %>% filter(ENTREZID %in% protein_list)
#re-normalize
a <- colSums(data_unique2[,3:72])
data_unique2[,3:72] <- sweep(x = data_unique2[,3:72]*1000000, MARGIN = 2, STATS = a, FUN = "/")
colSums(data_unique2[,3:72])
#add ENSEMBL ID
data_unique3 <- data_unique2
data_unique3$ENSEMBL <- mapIds(org.Hs.eg.db, keys=as.character(data_unique3$ENTREZID), column="ENSEMBL", keytype="ENTREZID", multiVals="first")
#manually annotate ATP5F1EP2, 432369, ENSG00000180389
data_unique3$ENSEMBL[data_unique3$ENTREZID == "432369"] <- "ENSG00000180389"
write.csv(data_unique3, file = "Clean_Data/data_unique_non-postmortem.csv", row.names = FALSE)

#make a SummerizedExperiment of final data
expdesign <- read.csv(file = "./Sample_Info/experimental_design.csv", header = TRUE)

proteins_unique <- data_unique3
columns <- c(3:72)
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
save(data_se, file = "./Results/data.RData")

#make an UpSet plot for each age after filter for proteins identified in non-postmortem brains
rownames(keep2) <- keep2$rowname
keep2 <- as.matrix(keep2[,-1])
keep2[!is.na(keep2)] <- 1
keep2[is.na(keep2)] <- 0
keep4 <- keep2[rownames(keep2) %in% rowData(data_se)$Gene.names,]
#draw an UpSet plot
m <- make_comb_mat(keep4, mode = "distinct", min_set_size = 50)
m2 <- m[comb_size(m) > 30]
pdf("Results/upSet_MSMSonly.pdf") 
UpSet(m2, set_order = colnames(keep4), comb_order = order(comb_size(m2), decreasing = TRUE),width = unit(3, "in"), height = unit(1.2, "in"))
dev.off()
