library("Seurat")
library("tidyverse")
library("future")
library("data.table")

setwd("/kriegsteinlab/data2/LiWang/snRNA_seq_for_PSD/human_snRNA_seq_10x")
plan("multiprocess", workers = 16)
options(future.globals.maxSize = 250000 * 1024^2)


mat <- fread("exprMatrix.tsv")
meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)
colnames(mat) = str_replace_all(colnames(mat),"\\.","-")
experiment.aggregate <- CreateSeuratObject(counts = mat, project = "human_snRNA_seq", meta.data=meta)

#add original tsne coordinates
tSNE <- read.table("tMinusSNE.coords.tsv", sep = "\t", header = F, row.names = 1)
tSNE <- as.matrix(tSNE)
colnames(tSNE) <- paste0("tSNE_", 1:2)
experiment.aggregate[["tSNE"]] <- CreateDimReducObject(embeddings = tSNE, key = "tSNE_", assay = DefaultAssay(experiment.aggregate))

#nomrlize data
experiment.aggregate <- NormalizeData(experiment.aggregate, normalization.method = "LogNormalize", scale.factor = 10000)

#save seurat objects
save(experiment.aggregate, file = "all_snRNAseq_Seurat_human_10x.RData")

#subset all data to focus on cortical neurons
experiment.CTX.neurons <- subset(experiment.aggregate, subset = subclass_label %in% c("L2/3 IT", "L5 ET", "L5 IT", "L5/6 NP", "L6 CT", "L6 IT", "L6 IT Car3", "L6b", "Lamp5", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip"))

#save seurat objects
save(experiment.CTX.neurons, file = "Cortical_neuron_snRNAseq_Seurat_human_10x.RData")

#pseudobulk on subclass level and data saving
res <- AggregateExpression(experiment.CTX.neurons, group.by= c("class_label", "subclass_label"), slot = "count")
write.csv(res[[1]],"pseudobulk_res_subclass_human_10x.csv", row.names = TRUE)

#pseudobulk on cluster level and data saving
res2 <- AggregateExpression(experiment.CTX.neurons, group.by= c("class_label", "subclass_label", "cluster_label"), slot = "count")
write.csv(res2[[1]],"pseudobulk_res_cluster_human_10x.csv", row.names = TRUE)





