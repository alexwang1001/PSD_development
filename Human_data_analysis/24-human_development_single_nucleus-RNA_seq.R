library(Seurat)
library(tidyverse)
library(future)
library(data.table)
library(scCustomize)
library(patchwork)

setwd("/kriegsteinlab/data3/LiWang/analysis/PSD_development/")
mat <- fread("single_cell_human_development_Dmitry/exprMatrix.tsv.gz")
meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)
colnames(mat) = str_replace_all(colnames(mat),"\\.","-")

experiment.aggregate <- CreateSeuratObject(counts = mat, project = "human_development_snRNA_seq", meta.data=meta)

#add original UMAP coordinates
UMAP <- read.table("UMAP.coords.tsv", sep = "\t", header = F, row.names = 1)
UMAP <- as.matrix(UMAP)
colnames(UMAP) <- paste0("UMAP_", 1:2)
experiment.aggregate[["UMAP"]] <- CreateDimReducObject(embeddings = UMAP, key = "UMAP_", assay = DefaultAssay(experiment.aggregate))

#nomrlize data
experiment.aggregate <- NormalizeData(experiment.aggregate, normalization.method = "LogNormalize", scale.factor = 10000)

#save seurat objects
saveRDS(experiment.aggregate, file = "results/single_cell_human_development_Dmitry/single_cell_human_development_Dmitry.rds")

#plot metadata
DimPlot_scCustom(experiment.aggregate, group.by = "age")
ggsave(file="results/single_cell_human_development_Dmitry/umap_age.png",width = 8, height = 6, units = "in")
DimPlot_scCustom(experiment.aggregate, group.by = "cluster", label = T)
ggsave(file="results/single_cell_human_development_Dmitry/umap_cluster.png",width = 8, height = 6, units = "in")

#check markers
human_markers <- c('FOXG1', 'OTX2', 'EN1', 'IRX3', 'HOXB4', 'ADCYAP1', 'NHLH1','EOMES', 'PPP1R17', 'RYR2', 'TNC', 'DCX', 'SYT1', 'RBFOX3', 'RELN', 'GAD1', 'GAD2', 'SST', 'SATB2', 'BCL11B', 'CUX2', 'NRGN', 'SLC17A7', 'SLC17A6', 'TLE4', 'FOXP2', 'PCP4', 'RORB', 'ADRA1A',
                   'CCN2', 'HCRTR2', 'SEMA3E', 'SYT6', 'GFRA1', 'LRATD2', 'TRPC7', 'TLL1', 'CCDC168', 'PCDH8', 'PIK3C2G', 'NPNT', 'FAM160A1', 'OLFML2B', 'POU3F1', 'SYP', 'ADAM33', 'THEMIS', 'TSHZ2', 'PRSS12', 'PDGFRA', 'PLP1', 'OLIG2', 'BCAS1', 'MBP', 'EGFR', 'C1QA', 'AIF1', 'P2RY12', 'GFAP', 'S100B', 'AQP4', 'SPARCL1', 'ID1',
                   'CLDN5', 'PECAM1', 'ASPN', 'CD93', 'CSPG4', 'RGS5', 'PDGFRB', 'MKI67', 'ASCL1', 'DLX1', 'DLX2', 'DLX6', 'SOX2', 'PAX6', 'HES5', 'VIM', "NES", 'HOPX', 'CRYAB', 'PDGFD', 'FEZF2', 'ALDH1L1', 'GSX2', 'PROX1', 'NR2F2', 'ADARB2', 'LHX6', 'MAF', 'NKX2-1', 'ISL1', 'MEIS1', 'MEIS2',
                   'ADARB2', 'LAMP5', 'PAX6', 'VIP', 'WIF1', 'SV2C', 'NDNF', 'CALB2', 'SCGN', 'PVALB', 'TAC1', 'ETS1', 'BTBD11', 'ADAMTS5', 'TBR1', 'HMGB2', 'LUM', 'OGN', 'DCN', 'APOE', 'CD38', 'CD74', 'CCR2', 'MS4A7', 'MRC1',
                   'CUX1', 'SOX5', 'NR2F1', 'GRIN2B', 'GRIN1', 'NR4A2', 'TTR', 'FOLR1','SP8', 'ARX', 'ZEB2', 'NXPH1', 'ACKR3', 'ERBB4', 'EPHB3', 'LHX8', 'ETV1', 'PBX3', 'TSHZ1', 'GBX1', 'ZIC1', 'EBF1', 'BEST3', 'WNT8B', 'ZIC2', 'ZIC3', 'GBX2', 'SIX3', 'EMX1', 'FGF17', 'SHH',
                   'BARHL2', 'LHX2', 'DBX1', 'WNT3A', 'TP73', 'IRF7', 'IRF8', 'RSPO3', 'RSPO2', 'LMX1A', 'NEUROD1', 'EPHB1', 'UNC5D', 'OLIG1', 'GJA1', 'ALDOC', 'SLC1A2', 'NWD2', 'HS3ST4', 'GREM2', 'ADAMTSL1', 'MEGF11', 'ADAMTSL1', 'MEGF11', 'SULF1', 'EGR2', 'NFATC2', 'AGT', 'CA2', 'HEPACAM', 'POU6F2',
                   'TFAP2C', 'MYT1L', 'POSTN', 'DLX6-AS1', 'ITGA2', 'FBXO32', 'DACH1', "TFAP2C", "MOXD1", "ITGA2", "ID3", "GRM3", "UNC5C", "NFIX", "SULF2", "NRP1", "IL1RAPL2")
human_markers_to_plot <- human_markers[human_markers %in% rownames(experiment.aggregate)]
for (i in human_markers_to_plot) {
  FeaturePlot_scCustom(experiment.aggregate, reduction = "UMAP", features = i, pt.size = 0.1, raster.dpi = c(1024, 1024))
  ggsave(file=paste0("results/single_cell_human_development_Dmitry/markers/", i,".png"),width = 6, height = 5, units = "in")
}
#find all marker genes
Idents(experiment.aggregate) <- "cluster"
markers <- FindAllMarkers(experiment.aggregate, assay = "RNA")
write.csv(markers, "results/single_cell_human_development_Dmitry/markers.csv")
#add annotation
metadata <- experiment.aggregate[[]]
annotation <- read.csv(file = "results/single_cell_human_development_Dmitry/annotation.csv")
metadata_new <- left_join(metadata, annotation)
rownames(metadata_new) <- rownames(metadata)
experiment.aggregate <- AddMetaData(experiment.aggregate, metadata_new[,7:9])

#plot metadata
DimPlot_scCustom(experiment.aggregate, group.by = "type")
ggsave(file="results/single_cell_human_development_Dmitry/umap_type.png",width = 8, height = 6, units = "in")

umap_age <- DimPlot_scCustom(experiment.aggregate, group.by = "age") + NoAxes()
umap_type <- DimPlot_scCustom(experiment.aggregate, group.by = "type") + NoAxes()
umap_age + umap_type
ggsave(file="results/single_cell_human_development_Dmitry/umap_age_type_combined.png",width = 11, height = 4, units = "in")

#plot marker genes
SATB2 <- FeaturePlot_scCustom(experiment.aggregate, reduction = "UMAP", features = "SATB2", pt.size = 0.1, raster.dpi = c(1024, 1024)) + NoLegend() + NoAxes()
LHX2 <- FeaturePlot_scCustom(experiment.aggregate, reduction = "UMAP", features = "LHX2", pt.size = 0.1, raster.dpi = c(1024, 1024)) + NoLegend() + NoAxes()
TLE4 <- FeaturePlot_scCustom(experiment.aggregate, reduction = "UMAP", features = "TLE4", pt.size = 0.1, raster.dpi = c(1024, 1024)) + NoLegend() + NoAxes()
GAD2 <- FeaturePlot_scCustom(experiment.aggregate, reduction = "UMAP", features = "GAD2", pt.size = 0.1, raster.dpi = c(1024, 1024)) + NoLegend() + NoAxes()
LHX6 <- FeaturePlot_scCustom(experiment.aggregate, reduction = "UMAP", features = "LHX6", pt.size = 0.1, raster.dpi = c(1024, 1024)) + NoLegend() + NoAxes()
PROX1 <- FeaturePlot_scCustom(experiment.aggregate, reduction = "UMAP", features = "PROX1", pt.size = 0.1, raster.dpi = c(1024, 1024)) + NoAxes()
SATB2 + LHX2 + TLE4 + GAD2 + LHX6 + PROX1 + plot_layout(ncol = 3)
ggsave(file="results/single_cell_human_development_Dmitry/selective_markers.pdf",width = 11, height = 8, units = "in")

#save annotated data
saveRDS(experiment.aggregate, file = "results/single_cell_human_development_Dmitry/single_cell_human_development_Dmitry_annotated.rds")
#subset data to only focus on neurons
experiment.aggregate.neurons <- subset(experiment.aggregate, subset = type %in% c("EN.IT", "EN.non.IT", "IN.MGE", "IN.CGE"))

#pseudobulk on the age and type level and data saving
res <- AggregateExpression(experiment.aggregate.neurons, group.by= c("subclass", "type", "age"), slot = "count")
write.csv(res[[1]],"results/single_cell_human_development_Dmitry/pseudobulk_res_type_human_development_Dmitry.csv", row.names = TRUE)
