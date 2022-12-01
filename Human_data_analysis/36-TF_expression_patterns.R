library("tidyverse")
library("SummarizedExperiment")
library("ComplexHeatmap")
library("WGCNA")
library("circlize")
library("seriation")


#import pseudobulk data
pseudobulk_res <- read.csv("single_cell_adult_human/pseudobulk_res_subclass_human_10x.csv", row.names = 1)
pseudobulk_res_cpm <- sweep(pseudobulk_res*1000000, MARGIN = 2, STATS = colSums(pseudobulk_res), FUN = "/")
assay <- log(pseudobulk_res_cpm+1,2)
coldata <- tibble(sample = colnames(pseudobulk_res_cpm))
coldata <- separate(coldata, 1 ,into = c("class", "subclass"), sep = "_")

#blue module
blue_ChEA3 <- read.table(file = "ChEA3/blue_Integrated_meanRank.tsv", header = T, sep = "\t")
assay_filt_blue <- assay[rownames(assay) %in% blue_ChEA3$TF[1:100],]
write.csv(assay_filt_blue, "Results/single_cell_adult_human/Blue module TF expression log CPM.csv")
#remove genes with too many missing values
index <- goodSamplesGenes(t(assay_filt_blue))
assay_filt_blue2 <- assay_filt_blue[index$goodGenes,]



#yellow module
yellow_ChEA3 <- read.table(file = "ChEA3/yellow_Integrated_meanRank.tsv", header = T, sep = "\t")
assay_filt_yellow <- assay[rownames(assay) %in% yellow_ChEA3$TF[1:100],]
write.csv(assay_filt_yellow, "Results/single_cell_adult_human/Yellow module TF expression log CPM.csv")
#remove genes with too many missing values
index <- goodSamplesGenes(t(assay_filt_yellow))
assay_filt_yellow2 <- assay_filt_yellow[index$goodGenes,]

#median value for top 100 TFs in blue and yellow modules
median_blue <- colMedians(as.matrix(assay_filt_blue2))
median_yellow <- colMedians(as.matrix(assay_filt_yellow2))
module_TF_median_mtx <- data.frame("blue" = scale(median_blue), "yellow" = scale(median_yellow))
rownames(module_TF_median_mtx) <- coldata$subclass

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

row_ha <- rowAnnotation(Class = coldata$class,
                        col = list(Class = c(GABAergic = 3, Glutamatergic = 7)),
                        gp = gpar(col = "white"))
ha_col <- HeatmapAnnotation(Modules = c("blue", "yellow"),
                            col = list(Modules = c("blue" = "blue", "yellow" = "gold")),
                            show_legend = FALSE,
                            gp = gpar(col = "white"))

o1 = seriate(dist(module_TF_median_mtx), method = "GW")
ht <- Heatmap(module_TF_median_mtx,
               name = "Scaled median expression of TFs",
               col = col_fun,
               border_gp = gpar(col = "black"),
               cluster_columns = F,
               cluster_rows = o1[[1]],
               row_dend_reorder = F,
               left_annotation = row_ha,
               top_annotation = ha_col,
               width = ncol(module_TF_median_mtx)*unit(8, "mm"), 
               height = nrow(module_TF_median_mtx)*unit(6, "mm")
)
ht <- draw(ht)
pdf(file = "Results/ChEA3/scaled_median_expression_of_TFs.pdf", width = 5, height = 5)
draw(ht, heatmap_legend_side = "right")
dev.off()


#heatmap for top TFs
assay_scaled_blue <- t(scale(t(assay_filt_blue2), scale = F))
colnames(assay_scaled_blue) <- coldata$subclass

col_fun <- circlize::colorRamp2(
  seq(-2.5, 2.5, (2.5/5)),
  rev(RColorBrewer::brewer.pal(11, "RdBu")))

ha_col <- HeatmapAnnotation(Class = coldata$class,
                            col = list(Class = c(GABAergic = 3, Glutamatergic = 7)),
                            show_legend = T,
                            gp = gpar(col = "white"))

o1 = seriate(dist(assay_scaled_blue), method = "OLO")
ht2 <- Heatmap(assay_scaled_blue,
               name = "Scaled expression",
               col = col_fun,
               border_gp = gpar(col = "black"),
               column_order = c(8,11,12,7,13,9,14,10,6,5,3,1,4,2),
               cluster_rows = o1[[1]],
               row_dend_reorder = F,
               show_row_names = F,
               top_annotation = ha_col,
               width = ncol(assay_scaled_blue)*unit(4, "mm"), 
               height = nrow(assay_scaled_blue)*unit(0.3, "mm")
)

pdf(file = "Results/ChEA3/blue_TF_expression.pdf", width = 7, height = 5)
draw(ht2, heatmap_legend_side = "right")
dev.off()


assay_scaled_yellow <- t(scale(t(assay_filt_yellow2), scale = F))
colnames(assay_scaled_yellow) <- coldata$subclass

col_fun <- circlize::colorRamp2(
  seq(-2.5, 2.5, (2.5/5)),
  rev(RColorBrewer::brewer.pal(11, "RdBu")))


o1 = seriate(dist(assay_scaled_yellow), method = "OLO")
ht3 <- Heatmap(assay_scaled_yellow,
               name = "Scaled expression",
               col = col_fun,
               border_gp = gpar(col = "black"),
               column_order = c(8,11,12,7,13,9,14,10,6,5,3,1,4,2),
               cluster_rows = o1[[1]],
               row_dend_reorder = F,
               show_row_names = F,
               width = ncol(assay_scaled_yellow)*unit(4, "mm"), 
               height = nrow(assay_scaled_yellow)*unit(0.3, "mm")
)

pdf(file = "Results/ChEA3/yellow_TF_expression.pdf", width = 7, height = 5)
draw(ht3, heatmap_legend_side = "right")
dev.off()

#heatmap vertical list
ht_list <- ht2 %v% ht3
pdf(file = "Results/ChEA3/blue_yellow_heatmap_list.pdf", width = 6, height = 8)
draw(ht_list, ht_gap = unit(0.2, "cm"))
dev.off()
