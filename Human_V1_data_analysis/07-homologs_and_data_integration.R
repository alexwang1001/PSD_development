library("tidyverse")
library("SummarizedExperiment")
library("DEP")


#load data
species_comparison <- readRDS("Results/species_comparison.rds")
human_V1 <- readRDS("Results/dep_corrected_BA17.rds")

#subset human_V1 data before combination
rowdata <- as.data.frame(rowData(species_comparison))
human_V1_subset <- human_V1[rowdata$Human_SYMBOL,]



#generate a single SummarizedExperiment
#assay
assay <- assay(species_comparison)
human_V1_assay <- assay(human_V1_subset)
#change names so that there is no same name with the PFC samples
colnames(human_V1_assay) <- paste0(colnames(human_V1_assay),"_V1")
assay <- cbind(assay, human_V1_assay)

#coldata
coldata <- as.data.frame(colData(species_comparison))[,c("label", "condition","replicate", "age", "log2_age_days", "species")]
coldata$species <- case_when(
  coldata$species == "human"~ "Human PFC",
  coldata$species == "macaque"~ "Macaque",
  coldata$species == "mouse" ~ "Mouse")

human_V1_coldata <- as.data.frame(colData(human_V1_subset))[,c("label", "condition","replicate", "AgeDays", "Log2AgeDays")]
colnames(human_V1_coldata) <- c("label", "condition","replicate", "age", "log2_age_days")
rownames(human_V1_coldata) <- colnames(human_V1_assay)
human_V1_coldata$label <- colnames(human_V1_assay)
human_V1_coldata$species <- "Human V1"

coldata <- rbind(coldata, human_V1_coldata)

#rowdata with module information
rowdata <- rowData(species_comparison)

#generate SummarizedExperiment
species_comparison_V1 <- SummarizedExperiment(assays=assay, rowData =rowdata, colData=coldata)
saveRDS(species_comparison_V1, file = "Results/species_comparison_V1.rds")
