library("tidyverse")
library("clusterProfiler")

#import module files
geneIDsModules <- read.csv(file = "Results/WGCNA_ORA/geneIDsModules.csv", row.names = 1)


#generate representative GO annotation
#load MSigDB data
c2cp_gmt <- read.gmt("MSigDB/c2.cp.v7.1.symbols.gmt")
c2cp_gmt$term <- tolower(c2cp_gmt$term)

GO <- geneIDsModules

reactome_translation <- c2cp_gmt %>% filter(term == "reactome_translation")
kegg_axon_guidance <- c2cp_gmt %>% filter(term == "kegg_axon_guidance")
reactome_signaling_by_rho_gtpases <- c2cp_gmt %>% filter(term == "reactome_signaling_by_rho_gtpases")
reactome_postsynaptic_signal_transmission <- c2cp_gmt %>% filter(term == "reactome_neurotransmitter_receptors_and_postsynaptic_signal_transmission")
reactome_neurexins_and_neuroligins <- c2cp_gmt %>% filter(term == "reactome_neurexins_and_neuroligins")

GO$reactome_translation <- ifelse(GO$SYMBOL %in% reactome_translation$gene, 1,0)
GO$kegg_axon_guidance <- ifelse(GO$SYMBOL %in% kegg_axon_guidance$gene, 1,0)
GO$reactome_signaling_by_rho_gtpases <- ifelse(GO$SYMBOL %in% reactome_signaling_by_rho_gtpases$gene, 1,0)
GO$reactome_postsynaptic_signal_transmission <- ifelse(GO$SYMBOL %in% reactome_postsynaptic_signal_transmission$gene, 1,0)
GO$reactome_neurexins_and_neuroligins <- ifelse(GO$SYMBOL %in% reactome_neurexins_and_neuroligins$gene, 1,0)

GO$interested <- ifelse(GO$reactome_translation | GO$kegg_axon_guidance | GO$reactome_signaling_by_rho_gtpases | GO$reactome_postsynaptic_signal_transmission | GO$reactome_neurexins_and_neuroligins, 1, 0)


#generate table of NDD genes with rare variants
NDD_variants <- read.csv(file = "Results/NDD_mutations/NDD_de_novo_variants_number_in_each_gene_corrected.csv")
NDD_variants <- unite(NDD_variants, "type", 1:2)
NDD_variants <- NDD_variants %>% spread(type, n)
NDD_variants$ENTREZID <- as.numeric(NDD_variants$ENTREZID)
rare_variants <- left_join(geneIDsModules, NDD_variants, by = c("ENTREZID", "SYMBOL"))
rare_variants[is.na(rare_variants)] <- 0

rare_variants$epilepsy_missense <- ifelse(rare_variants$epilepsy_missense>2, 1, 0)
rare_variants$ID_missense <- ifelse(rare_variants$ID_missense>2, 1, 0)
rare_variants$DD_missense <- ifelse(rare_variants$DD_missense>2, 1, 0)
rare_variants$ASD_missense <- ifelse(rare_variants$ASD_missense>2, 1, 0)
rare_variants$SCZ_missense <- ifelse(rare_variants$SCZ_missense>2, 1, 0)
rare_variants$Not_epilepsy_missense <- ifelse(rare_variants$epilepsy_missense, 0, 1)
rare_variants$Not_ID_missense <- ifelse(rare_variants$ID_missens, 0, 1)
rare_variants$Not_DD_missense <- ifelse(rare_variants$DD_missense, 0, 1)
rare_variants$Not_ASD_missense <- ifelse(rare_variants$ASD_missense, 0, 1)
rare_variants$Not_SCZ_missense <- ifelse(rare_variants$SCZ_missense, 0, 1)

rare_variants$epilepsy_LGD <- ifelse(rare_variants$epilepsy_LGD>1, 1, 0)
rare_variants$ID_LGD <- ifelse(rare_variants$ID_LGD>1, 1, 0)
rare_variants$DD_LGD <- ifelse(rare_variants$DD_LGD>1, 1, 0)
rare_variants$ASD_LGD <- ifelse(rare_variants$ASD_LGD>1, 1, 0)
rare_variants$SCZ_LGD <- ifelse(rare_variants$SCZ_LGD>1, 1, 0)

rare_variants$Not_epilepsy_LGD <- ifelse(rare_variants$epilepsy_LGD, 0, 1)
rare_variants$Not_ID_LGD <- ifelse(rare_variants$ID_LGD, 0, 1)
rare_variants$Not_DD_LGD <- ifelse(rare_variants$DD_LGD, 0, 1)
rare_variants$Not_ASD_LGD <- ifelse(rare_variants$ASD_LGD, 0, 1)
rare_variants$Not_SCZ_LGD <- ifelse(rare_variants$SCZ_LGD, 0, 1)

rare_variants$missense <- ifelse(rare_variants$epilepsy_missense | rare_variants$ID_missense | rare_variants$DD_missense | rare_variants$ASD_missense | rare_variants$SCZ_missense, 1, 0)
rare_variants$LGD <- ifelse(rare_variants$epilepsy_LGD | rare_variants$ID_LGD | rare_variants$DD_LGD |rare_variants$ASD_LGD | rare_variants$SCZ_LGD, 1, 0)
rare_variants$Not_missense <- ifelse(rare_variants$missense, 0, 1)
rare_variants$Not_LGD <- ifelse(rare_variants$LGD, 0, 1)
rare_variants$mutation <- ifelse(rare_variants$missense | rare_variants$LGD, 1, 0)


#generate table of neuropsychiatric genes with GWAS pvalue and FDR
ASD_GWAS <- read.csv(file = "Results/PPI_network/annotations/top_GWAS_genes_for_PPI_network/ASD.csv")
ASD_GWAS$ASD_FDR <- p.adjust(ASD_GWAS$P, method = "BH")
colnames(ASD_GWAS)[1:2] <- c("ENTREZID", "ASD_p_value")
GWAS <- left_join(geneIDsModules, ASD_GWAS, by = "ENTREZID")

SCZ_GWAS <- read.csv(file = "Results/PPI_network/annotations/top_GWAS_genes_for_PPI_network/SCZ.csv")
SCZ_GWAS$SCZ_FDR <- p.adjust(SCZ_GWAS$P, method = "BH")
colnames(SCZ_GWAS)[1:2] <- c("ENTREZID", "SCZ_p_value")
GWAS <- left_join(GWAS, SCZ_GWAS, by = "ENTREZID")

BPD_GWAS <- read.csv(file = "Results/PPI_network/annotations/top_GWAS_genes_for_PPI_network/BPD.csv")
BPD_GWAS$BPD_FDR <- p.adjust(BPD_GWAS$P, method = "BH")
colnames(BPD_GWAS)[1:2] <- c("ENTREZID", "BPD_p_value")
GWAS <- left_join(GWAS, BPD_GWAS, by = "ENTREZID")

MDD_GWAS <- read.csv(file = "Results/PPI_network/annotations/top_GWAS_genes_for_PPI_network/MDD.csv")
MDD_GWAS$MDD_FDR <- p.adjust(MDD_GWAS$P, method = "BH")
colnames(MDD_GWAS)[1:2] <- c("ENTREZID", "MDD_p_value")
GWAS <- left_join(GWAS, MDD_GWAS, by = "ENTREZID")

GWAS <- GWAS[,c(1:3,7,9,11)]
GWAS[is.na(GWAS)] <- 1
#set cutoff at 5% FDR
GWAS$BPD_FDR <- ifelse(GWAS$BPD_FDR<0.05, 1, 0)
GWAS$SCZ_FDR <- ifelse(GWAS$SCZ_FDR<0.05, 1, 0)
GWAS$MDD_FDR <- ifelse(GWAS$MDD_FDR<0.05, 1, 0)
colnames(GWAS)[4:6] <- c("BPD_GWAS","SCZ_GWAS", "MDD_GWAS")
GWAS$Not_BPD_GWAS <- ifelse(GWAS$BPD_GWAS, 0, 1)
GWAS$Not_SCZ_GWAS <- ifelse(GWAS$SCZ_GWAS, 0, 1)
GWAS$Not_MDD_GWAS <- ifelse(GWAS$MDD_GWAS, 0, 1)
GWAS$GWAS <- ifelse(GWAS$BPD_GWAS | GWAS$SCZ_GWAS | GWAS$MDD_GWAS, 1, 0)

#generate table of DEGs in neuropsychiatric genes
diff <- read.csv("Results/differential_expression_in_psychiatric_disorder/NDD_diff_gene_list.csv", row.names = 1)
#due to isoform and difference between ENTREZ and ENSEMBL annotation, there are duplication of the results
#we first make them unique
diff <- unique(diff)
diff <- diff %>% mutate(var = 1) %>% spread(term, var)
diff$ENTREZID <- as.numeric(diff$ENTREZID)
diff <- left_join(geneIDsModules, diff, by = c("ENTREZID", "SYMBOL"))
for (col in colnames(diff)) {
  diff[,col] <- replace_na(diff[,col], 0)
}

diff <- diff[,c(1:3,6,10,4)]
diff$Not_BPD_diff_down <- ifelse(diff$BPD_diff_down, 0, 1)
diff$Not_SCZ_diff_down <- ifelse(diff$SCZ_diff_down, 0, 1)
diff$Not_ASD_diff_down <- ifelse(diff$ASD_diff_down, 0, 1)
diff$diff_down <- ifelse(diff$BPD_diff_down | diff$SCZ_diff_down | diff$ASD_diff_down, 1, 0)

#activity-dependent proteins
activity <- read.csv(file = "Results/activity_dependent_proteome_changes/activity_dependent_protein_list.csv")
activity <- activity %>% mutate(var = 1) %>% spread(condition, var)
activity$ENTREZID <- as.numeric(activity$ENTREZID)
activity <- left_join(geneIDsModules, activity, by = c("ENTREZID", "SYMBOL"))
for (col in colnames(activity)) {
  activity[,col] <- replace_na(activity[,col], 0)
}
activity <- activity[,c(1:4,6)]
activity$Not_BIC_24h <- ifelse(activity$BIC_24h, 0, 1)
activity$Not_TTX_24h <- ifelse(activity$TTX_24h, 0, 1)
activity$activity_dependent <- ifelse(activity$BIC_24h | activity$TTX_24h, 1, 0)
activity <- activity[,c(1:3,6,4,7,5,8)]

#generate csv files for all tables for PPI network annotation
write.csv(GO, file = "Results/PPI_network/annotations/GO.csv", row.names = FALSE)
write.csv(rare_variants, file = "Results/PPI_network/annotations/NDD_rare_variants.csv", row.names = FALSE)
write.csv(GWAS, file = "Results/PPI_network/annotations/GWAS.csv", row.names = FALSE)
write.csv(diff, file = "Results/PPI_network/annotations/Differentially_expressed.csv", row.names = FALSE)
write.csv(activity, file = "Results/PPI_network/annotations/activity_dependent.csv", row.names = FALSE)
