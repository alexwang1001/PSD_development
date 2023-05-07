library("clusterProfiler")
library("tidyverse")
library("readxl")

#differentially expressed genes in ASD, SCZ, BPD, and MDD
ASD_diff <- read_excel(path = "Gandal_differential_expression/ASD_diff.xlsx", na = "NA")
ASD_diff <- ASD_diff[!is.na(ASD_diff$ENTREZID),]
diff_ID <- bitr(ASD_diff$ENTREZID,fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
diff_ID$ENTREZID <- as.numeric(diff_ID$ENTREZID)
ASD_diff <- dplyr::left_join(ASD_diff, diff_ID)
ASD_diff <- ASD_diff[!is.na(ASD_diff$SYMBOL),]
ASD_diff_up <- ASD_diff %>% filter(log2FC > 0, FDR<0.05)
ASD_diff_down <- ASD_diff %>% filter(log2FC < 0, FDR<0.05)
ASD_diff_up <- data.frame("term" = "ASD_diff_up", "SYMBOL" = ASD_diff_up$SYMBOL)
ASD_diff_down <- data.frame("term" = "ASD_diff_down", "SYMBOL" = ASD_diff_down$SYMBOL)

BPD_diff <- read_excel(path = "Gandal_differential_expression/BPD_diff.xlsx", na = "NA")
BPD_diff <- BPD_diff[!is.na(BPD_diff$ENTREZID),]
BPD_diff <- dplyr::left_join(BPD_diff, diff_ID)
BPD_diff <- BPD_diff[!is.na(BPD_diff$SYMBOL),]
BPD_diff_up <- BPD_diff %>% filter(log2FC > 0, FDR<0.05)
BPD_diff_down <- BPD_diff %>% filter(log2FC < 0, FDR<0.05)
BPD_diff_up <- data.frame("term" = "BPD_diff_up", "SYMBOL" = BPD_diff_up$SYMBOL)
BPD_diff_down <- data.frame("term" = "BPD_diff_down", "SYMBOL" = BPD_diff_down$SYMBOL)

MDD_diff <- read_excel(path = "Gandal_differential_expression/MDD_diff.xlsx", na = "NA")
MDD_diff <- MDD_diff[!is.na(MDD_diff$ENTREZID),]
MDD_diff <- dplyr::left_join(MDD_diff, diff_ID)
MDD_diff <- MDD_diff[!is.na(MDD_diff$SYMBOL),]
MDD_diff_up <- MDD_diff %>% filter(log2FC > 0, FDR<0.05)
MDD_diff_down <- MDD_diff %>% filter(log2FC < 0, FDR<0.05)
MDD_diff_up <- data.frame("term" = "MDD_diff_up", "SYMBOL" = MDD_diff_up$SYMBOL)
MDD_diff_down <- data.frame("term" = "MDD_diff_down", "SYMBOL" = MDD_diff_down$SYMBOL)

SCZ_diff <- read_excel(path = "Gandal_differential_expression/SCZ_diff.xlsx", na = "NA")
SCZ_diff <- SCZ_diff[!is.na(SCZ_diff$ENTREZID),]
SCZ_diff <- dplyr::left_join(SCZ_diff, diff_ID)
SCZ_diff <- SCZ_diff[!is.na(SCZ_diff$SYMBOL),]
SCZ_diff_up <- SCZ_diff %>% filter(log2FC > 0, FDR<0.05)
SCZ_diff_down <- SCZ_diff %>% filter(log2FC < 0, FDR<0.05)
SCZ_diff_up <- data.frame("term" = "SCZ_diff_up", "SYMBOL" = SCZ_diff_up$SYMBOL)
SCZ_diff_down <- data.frame("term" = "SCZ_diff_down", "SYMBOL" = SCZ_diff_down$SYMBOL)

TERM2GENE <- rbind(ASD_diff_up, ASD_diff_down, BPD_diff_up, BPD_diff_down, MDD_diff_up, MDD_diff_down, SCZ_diff_up, SCZ_diff_down)

#check if all gene names are official symbols
ids <- bitr(TERM2GENE$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
geneIDs <- dplyr::left_join(TERM2GENE, ids, )
geneIDs[is.na(geneIDs$ENTREZID),]
write.csv(geneIDs, file = "Results/differential_expression_in_psychiatric_disorder/NDD_diff_gene_list.csv")

backgroundIDs <- read.csv(file = "./Results/WGCNA_ORA/backgroundIDs.csv")
universe <- backgroundIDs$SYMBOL

load(file = "Results/module_name_for_ORA.RData")
length(unlist(module_name_for_ORA))
geneIDsModules <- read.csv(file = "Results/WGCNA_ORA/geneIDsModules.csv", header = TRUE, row.names = 1)
module_name_for_ORA$all <- geneIDsModules$SYMBOL
length(unlist(module_name_for_ORA))

source(file = "functions/enricher_disorder.R")
source(file = "functions/compareCluster_disorder.R")
NDD_diff_ORA <- compareCluster_disorder(geneClusters = module_name_for_ORA,
                                         fun="enricher_disorder",
                                         universe = universe,
                                         TERM2GENE = TERM2GENE)

NDD_diff_ORA_result <- NDD_diff_ORA@compareClusterResult
write.csv(NDD_diff_ORA_result, file = "Results/differential_expression_in_psychiatric_disorder/NDD_diff_ORA_result.csv", row.names = FALSE)

save(NDD_diff_ORA, file = "Results/NDD_diff_ORA.RData")


