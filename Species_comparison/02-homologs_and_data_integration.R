library("tidyverse")
library("SummarizedExperiment")
library("DEP")
library(biomaRt)

#load data
load(file = "../Human PSD filtered by non post-mortem tissue/Results/dep_corrected.RData")
human_dep <- dep_corrected
human_rowdata <- as.data.frame(rowData(human_dep))
human_rowdata$ENTREZID <- as.numeric(human_rowdata$ENTREZID)

load(file = "../Macaque PSD/Results/dep.RData")
macaque_dep <- dep
macaque_rowdata <- as.data.frame(rowData(macaque_dep))

load(file = "../Mouse PSD/Results/dep.RData")
mouse_dep <- dep
mouse_rowdata <- as.data.frame(rowData(mouse_dep))

############################################################
#identify homologs
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
human <-  useMart(biomart="ensembl", dataset = "hsapiens_gene_ensembl", verbose = TRUE, host = "dec2021.archive.ensembl.org")
macaque <-  useMart(biomart="ensembl", dataset = "mmulatta_gene_ensembl", verbose = TRUE, host = "dec2021.archive.ensembl.org")
mouse <-  useMart(biomart="ensembl", dataset = "mmusculus_gene_ensembl", verbose = TRUE, host = "dec2021.archive.ensembl.org")
View(listAttributes(human))

human_to_macaque <- getLDS(attributes = c("ensembl_gene_id"),
                           filters = "ensembl_gene_id", values = human_rowdata$ENSEMBL, mart = human,
                           attributesL = c("ensembl_gene_id"), martL = macaque, uniqueRows = TRUE)

macaque_to_human <- getLDS(attributes = c("ensembl_gene_id"),
                           filters = "ensembl_gene_id", values = macaque_rowdata$ENSEMBL, mart = macaque,
                           attributesL = c("ensembl_gene_id"), martL = human, uniqueRows = TRUE)

human_to_mouse <- getLDS(attributes = c("ensembl_gene_id"),
                           filters = "ensembl_gene_id", values = human_rowdata$ENSEMBL, mart = human,
                           attributesL = c("ensembl_gene_id"), martL = mouse, uniqueRows = TRUE)

mouse_to_human <- getLDS(attributes = c("ensembl_gene_id"),
                           filters = "ensembl_gene_id", values = mouse_rowdata$ENSEMBL, mart = mouse,
                           attributesL = c("ensembl_gene_id"), martL = human, uniqueRows = TRUE)

macaque_to_mouse <- getLDS(attributes = c("ensembl_gene_id"),
                           filters = "ensembl_gene_id", values = macaque_rowdata$ENSEMBL, mart = macaque,
                           attributesL = c("ensembl_gene_id"), martL = mouse, uniqueRows = TRUE)

mouse_to_macaque <- getLDS(attributes = c("ensembl_gene_id"),
                           filters = "ensembl_gene_id", values = mouse_rowdata$ENSEMBL, mart = mouse,
                           attributesL = c("ensembl_gene_id"), martL = macaque, uniqueRows = TRUE)

#filter to only retain one to one homologs
human_to_macaque_one_to_one <- human_to_macaque[!(human_to_macaque$Gene.stable.ID %in% human_to_macaque[duplicated(human_to_macaque$Gene.stable.ID),]$Gene.stable.ID) &
                                                  !(human_to_macaque$Gene.stable.ID.1 %in% human_to_macaque[duplicated(human_to_macaque$Gene.stable.ID.1),]$Gene.stable.ID.1),]

macaque_to_human_one_to_one <- macaque_to_human[!(macaque_to_human$Gene.stable.ID %in% macaque_to_human[duplicated(macaque_to_human$Gene.stable.ID),]$Gene.stable.ID) &
                                                  !(macaque_to_human$Gene.stable.ID.1 %in% macaque_to_human[duplicated(macaque_to_human$Gene.stable.ID.1),]$Gene.stable.ID.1),]

human_to_mouse_one_to_one <- human_to_mouse[!(human_to_mouse$Gene.stable.ID %in% human_to_mouse[duplicated(human_to_mouse$Gene.stable.ID),]$Gene.stable.ID) &
                                              !(human_to_mouse$Gene.stable.ID.1 %in% human_to_mouse[duplicated(human_to_mouse$Gene.stable.ID.1),]$Gene.stable.ID.1),]

mouse_to_human_one_to_one <- mouse_to_human[!(mouse_to_human$Gene.stable.ID %in% mouse_to_human[duplicated(mouse_to_human$Gene.stable.ID),]$Gene.stable.ID) &
                                              !(mouse_to_human$Gene.stable.ID.1 %in% mouse_to_human[duplicated(mouse_to_human$Gene.stable.ID.1),]$Gene.stable.ID.1),]

macaque_to_mouse_one_to_one <- macaque_to_mouse[!(macaque_to_mouse$Gene.stable.ID %in% macaque_to_mouse[duplicated(macaque_to_mouse$Gene.stable.ID),]$Gene.stable.ID) &
                                              !(macaque_to_mouse$Gene.stable.ID.1 %in% macaque_to_mouse[duplicated(macaque_to_mouse$Gene.stable.ID.1),]$Gene.stable.ID.1),]

mouse_to_macaque_one_to_one <- mouse_to_macaque[!(mouse_to_macaque$Gene.stable.ID %in% mouse_to_macaque[duplicated(mouse_to_macaque$Gene.stable.ID),]$Gene.stable.ID) &
                                                  !(mouse_to_macaque$Gene.stable.ID.1 %in% mouse_to_macaque[duplicated(mouse_to_macaque$Gene.stable.ID.1),]$Gene.stable.ID.1),]

#filter to only retain one to one homolgs across all three species and present in all three datasets
human_one_to_one <- Reduce(intersect, list(human_to_macaque_one_to_one$Gene.stable.ID,
                                           macaque_to_human_one_to_one$Gene.stable.ID.1,
                                           human_to_mouse_one_to_one$Gene.stable.ID,
                                           mouse_to_human_one_to_one$Gene.stable.ID.1)
)

macaque_one_to_one <- Reduce(intersect, list(human_to_macaque_one_to_one$Gene.stable.ID.1,
                                             macaque_to_human_one_to_one$Gene.stable.ID,
                                             macaque_to_mouse_one_to_one$Gene.stable.ID,
                                             mouse_to_macaque_one_to_one$Gene.stable.ID.1)
)

mouse_one_to_one <- Reduce(intersect, list(human_to_mouse_one_to_one$Gene.stable.ID.1,
                                           mouse_to_human_one_to_one$Gene.stable.ID,
                                           macaque_to_mouse_one_to_one$Gene.stable.ID.1,
                                           mouse_to_macaque_one_to_one$Gene.stable.ID)
)

#get linked table for genes with one to one holomogs in all three species
colnames(human_to_macaque_one_to_one) <- c("Human_ENSEMBL", "Macaque_ENSEMBL")
colnames(human_to_mouse_one_to_one) <- c("Human_ENSEMBL", "Mouse_ENSEMBL")
human_macaque_mouse_one_to_one <- human_to_macaque_one_to_one %>%
  left_join(human_to_mouse_one_to_one) %>%
  filter(Human_ENSEMBL %in% human_one_to_one, Macaque_ENSEMBL %in% macaque_one_to_one, Mouse_ENSEMBL %in% mouse_one_to_one)

human_rowdata_to_join <- human_rowdata[,c("Protein.IDs", "Gene.names", "ENTREZID", "ENSEMBL")]
colnames(human_rowdata_to_join) <- c("Human_UNIPROT", "Human_SYMBOL", "Human_ENTREZID", "Human_ENSEMBL")

macaque_rowdata_to_join <- macaque_rowdata[,c("Protein.IDs", "Gene.names", "ENTREZID", "ENSEMBL")]
colnames(macaque_rowdata_to_join) <- c("Macaque_UNIPROT", "Macaque_SYMBOL", "Macaque_ENTREZID", "Macaque_ENSEMBL")

mouse_rowdata_to_join <- mouse_rowdata[,c("Protein.IDs", "Gene.names", "ENTREZID", "ENSEMBL")]
colnames(mouse_rowdata_to_join) <- c("Mouse_UNIPROT", "Mouse_SYMBOL", "Mouse_ENTREZID", "Mouse_ENSEMBL")


human_macaque_mouse_one_to_one <- human_macaque_mouse_one_to_one %>% left_join(human_rowdata_to_join) %>%
  left_join(macaque_rowdata_to_join) %>%
  left_join(mouse_rowdata_to_join)

#subset dep objects to focus on human_macaque_mouse_one_to_one genes
human_dep_subset <- human_dep[human_macaque_mouse_one_to_one$Human_SYMBOL,]
macaque_dep_subset <- macaque_dep[human_macaque_mouse_one_to_one$Macaque_SYMBOL,]
mouse_dep_subset <- mouse_dep[human_macaque_mouse_one_to_one$Mouse_SYMBOL,]

#generate a single SummarizedExperiment
#assay
human_assay <- assay(human_dep_subset)
macaque_assay <- assay(macaque_dep_subset)
mouse_assay <- assay(mouse_dep_subset)
assay <- cbind(human_assay, macaque_assay, mouse_assay)

#coldata
human_coldata <- as.data.frame(colData(human_dep))[,c("label", "condition","replicate", "AgeDays", "Log2AgeDays")]
colnames(human_coldata) <- c("label", "condition","replicate", "age", "log2_age_days")
human_coldata$species <- "human"
macaque_coldata <- as.data.frame(colData(macaque_dep))[,c("label", "condition","replicate", "age", "log2_age_days")]
macaque_coldata$species <- "macaque"
mouse_coldata <- as.data.frame(colData(mouse_dep))[,c("label", "condition","replicate", "age", "Log2AgeDays")]
colnames(mouse_coldata) <- c("label", "condition","replicate", "age", "log2_age_days")
mouse_coldata$species <- "mouse"
coldata <- rbind(human_coldata, macaque_coldata, mouse_coldata)
#rowdata with module information
rowdata <- human_macaque_mouse_one_to_one

human_geneIDsModules <- read.csv("../Human PSD filtered by non post-mortem tissue/Results/WGCNA_ORA/geneIDsModules.csv", row.names = 1)
colnames(human_geneIDsModules) <- c("Human_SYMBOL", "Human_ENTREZID", "Human_Module")

macaque_geneIDsModules <- read.csv("../Macaque PSD/Results/WGCNA_ORA/geneIDsModules.csv", row.names = 1)
colnames(macaque_geneIDsModules) <- c("Macaque_SYMBOL", "Macaque_ENTREZID", "Macaque_ENSEMBL", "Macaque_Module")

mouse_geneIDsModules <- read.csv("../Mouse PSD/Results/WGCNA_ORA/geneIDsModules.csv", row.names = 1)
colnames(mouse_geneIDsModules) <- c("Mouse_SYMBOL", "Mouse_ENTREZID", "Mouse_ENSEMBL", "Mouse_Module")

rowdata <- rowdata %>% left_join(human_geneIDsModules) %>%
  left_join(macaque_geneIDsModules) %>%
  left_join(mouse_geneIDsModules)
#generate SummarizedExperiment
species_comparison <- SummarizedExperiment(assays=assay, rowData =rowdata, colData=coldata)
saveRDS(species_comparison, file = "Results/species_comparison.rds")
