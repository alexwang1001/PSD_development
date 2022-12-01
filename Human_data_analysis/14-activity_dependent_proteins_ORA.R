library("tidyverse")
library("clusterProfiler")
library("readxl")

#activity-dependent protein list were downloaded from PMID: 29447110
#proteins without ENTREZID or are located at mitochondria are removed

#put all gene lists together
geneIDs <- read.csv("Results/activity_dependent_proteome_changes/activity_dependent_protein_list.csv")
TERM2GENE <- geneIDs[,c(3,2)] 

backgroundIDs <- read.csv(file = "./Results/WGCNA_ORA/backgroundIDs.csv")
universe <- backgroundIDs$SYMBOL

load(file = "Results/module_name_for_ORA.RData")
geneIDsModules <- read.csv(file = "Results/WGCNA_ORA/geneIDsModules.csv", header = TRUE, row.names = 1)
module_name_for_ORA$all <- geneIDsModules$SYMBOL

source(file = "functions/enricher_disorder.R")
source(file = "functions/compareCluster_disorder.R")
activity_dependent_proteins_ORA <- compareCluster_disorder(geneClusters = module_name_for_ORA,
                                        fun="enricher_disorder",
                                        universe = universe,
                                        TERM2GENE = TERM2GENE)

activity_dependent_proteins_ORA_result <- as.data.frame(activity_dependent_proteins_ORA@compareClusterResult)
write.csv(activity_dependent_proteins_ORA_result, file = "Results/activity_dependent_proteome_changes/activity_dependent_proteins_ORA_result.csv", row.names = FALSE)
save(activity_dependent_proteins_ORA, file = "Results/activity_dependent_proteins_ORA.RData")
