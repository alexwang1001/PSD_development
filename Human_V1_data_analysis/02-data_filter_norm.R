library("DEP")
library("SummarizedExperiment")
library("limma")
library("vsn")
library("tidyverse")

se <- readRDS("Results/data_all.rds")

#load PFC data
load("PFC_data_norm.RData")
#subset BA17 data to focus on the same set of proteins from the PFC data
se2 =se[rownames(data_norm),]
#normalize data using the same vsn model used for PFC data
vsn.fit <- readRDS("vsn_reference.rds")
vsn.fit.BA17 <- vsn::vsnMatrix(2 ^ assay(se2), reference = vsn.fit)
data_norm_BA17 <- se2
assay(data_norm_BA17) <- vsn::predict(vsn.fit.BA17, 2 ^ assay(se2))

meanSdPlot(assay(se2))
meanSdPlot(assay(data_norm_BA17))
plotMD(assay(se2))
plotMD(assay(data_norm_BA17))
plot_normalization(se2, data_norm_BA17)

#Imputation
set.seed(6)
data_imp <- impute(data_norm_BA17, fun = "MinProb", q = 0.01)
plot_imputation(data_norm_BA17, data_imp)

saveRDS(data_norm_BA17, file = "Results/data_norm_all.rds")
saveRDS(data_imp, file = "Results/data_imp_all.rds")
