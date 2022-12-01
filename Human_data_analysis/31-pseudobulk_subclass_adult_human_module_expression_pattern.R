library("tidyverse")
library("SummarizedExperiment")
library("ComplexHeatmap")
library("circlize")
library("seriation")
library("WGCNA")
library("viridis")


#subclass data
pseudobulk_res <- read.csv("single_cell_adult_human/pseudobulk_res_subclass_human_10x.csv", row.names = 1)
geneIDsModules <- read.csv(file = "Results/WGCNA_ORA/geneIDsModules.csv", row.names = 1)
pseudobulk_res_PSD <- pseudobulk_res[rownames(pseudobulk_res) %in% geneIDsModules$SYMBOL,]

#re-normalize the data
pseudobulk_res_PSD_cpm <- sweep(pseudobulk_res_PSD*1000000, MARGIN = 2, STATS = colSums(pseudobulk_res_PSD), FUN = "/")

coldata <- tibble(sample = colnames(pseudobulk_res_PSD_cpm))
coldata <- separate(coldata, 1 ,into = c("class", "subclass"), sep = "_")

rowdata <- tibble(SYMBOL = rownames(pseudobulk_res_PSD_cpm))
rowdata <- left_join(rowdata, geneIDsModules, by = "SYMBOL")

#log2 transform
assay <- log(pseudobulk_res_PSD_cpm+1,2)

#generate a summarized experiment

se <- SummarizedExperiment(assays = as.matrix(assay), colData = coldata, 
                           rowData = rowdata)
colData(se)

saveRDS(se, file = "Results/pseudo_bulk_subclass_adult_human.rds")

assay2 <- cbind(rowdata, assay)
write.csv(assay2,"Results/single_cell_adult_human/log cpm counts.csv", row.names = F)

#calculate median values of each module and standardize it
get_scaled_median_for_modules <- function(se) {
  rowdata<- rowData(se)
  coldata <- colData(se)
  modules <- unique(rowdata[["geneModules"]])
  for (module in modules) {
    module_se <- se[rowdata$geneModules == module,]
    module_se_expr <- data.frame(assay(module_se))
    module_median <- apply(module_se_expr, 2, median)
    module_median_scaled <- scale(module_median)
    colData(se)[module] <- module_median_scaled
  }
  data.frame(colData(se))
}

module_median <- get_scaled_median_for_modules(se)

#focus on blue and yellow modules as they are strongly preserved between protein and RNA
module_median$grey <- NULL
module_median$brown <- NULL
module_median$turquoise <- NULL
module_median_long <- gather(module_median, key = "module", value = "median", 3:4)

#heatmap
module_median_mtx <- module_median[,-c(1,2)]
module_median_mtx <- module_median_mtx[,c(2,1)]
rownames(module_median_mtx) <- module_median$subclass

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

row_ha <- rowAnnotation(Class = module_median$class,
                        col = list(Class = c(GABAergic = 3, Glutamatergic = 7)),
                        gp = gpar(col = "white"))
ha_col <- HeatmapAnnotation(Modules = c("blue", "yellow"),
                            col = list(Modules = c("blue" = "blue", "yellow" = "gold")),
                            show_legend = FALSE,
                            gp = gpar(col = "white"))

o1 = seriate(dist(module_median_mtx), method = "OLO")
ht1 <- Heatmap(module_median_mtx,
               name = "Scaled median expression",
               col = col_fun,
               border_gp = gpar(col = "black"),
               cluster_columns = F,
               cluster_rows = o1[[1]],
               row_gap = unit(5, "mm",),
               row_dend_reorder = F,
               left_annotation = row_ha,
               top_annotation = ha_col,
               heatmap_legend_param = list(
                 legend_direction = "horizontal"
               ),
               width = ncol(module_median_mtx)*unit(8, "mm"), 
               height = nrow(module_median_mtx)*unit(6, "mm")
)
ht1 <- draw(ht1)
pdf(file = "Results/single_cell_adult_human/heatmap.pdf", width = 4, height = 5)
draw(ht1, heatmap_legend_side = "top")
dev.off()

#plot blue module genes in individual neuronal subtypes
assay_filt_blue <- assay[rowdata$geneModules == "blue",]
#remove genes with too many missing values
index <- goodSamplesGenes(t(assay_filt_blue))
assay_filt_blue2 <- assay_filt_blue[index$goodGenes,]
#heatmap
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
o2 = seriate(dist(t(assay_scaled_blue)), method = "OLO")
ht2 <- Heatmap(assay_scaled_blue,
               name = "Scaled expression",
               col = col_fun,
               border_gp = gpar(col = "black"),
               column_order = row_order(ht1),
               cluster_columns = F,
               cluster_rows = o1[[1]],
               row_gap = unit(5, "mm",),
               row_dend_reorder = F,
               show_row_names = F,
               top_annotation = ha_col,
               width = ncol(assay_scaled_blue)*unit(4, "mm"), 
               height = nrow(assay_scaled_blue)*unit(0.3, "mm")
)
pdf(file = "Results/single_cell_adult_human/blue_module_expression.pdf", width = 6, height = 6)
draw(ht2, heatmap_legend_side = "right")
dev.off()


#plot yellow module genes in individual neuronal subtypes
assay_filt_yellow <- assay[rowdata$geneModules == "yellow",]
#remove genes with too many missing values
index <- goodSamplesGenes(t(assay_filt_yellow))
assay_filt_yellow2 <- assay_filt_yellow[index$goodGenes,]
#heatmap
assay_scaled_yellow <- t(scale(t(assay_filt_yellow2), scale = F))
colnames(assay_scaled_yellow) <- coldata$subclass


o1 = seriate(dist(assay_scaled_yellow), method = "GW")
o2 = seriate(dist(t(assay_scaled_yellow)), method = "OLO")
ht3 <- Heatmap(assay_scaled_yellow,
               name = "Scaled expression",
               col = col_fun,
               border_gp = gpar(col = "black"),
               cluster_columns = F,
               cluster_rows = o1[[1]],
               row_gap = unit(5, "mm",),
               row_dend_reorder = F,
               show_row_names = F,
               width = ncol(assay_scaled_yellow)*unit(4, "mm"), 
               height = nrow(assay_scaled_yellow)*unit(0.3, "mm")
)
pdf(file = "Results/single_cell_adult_human/yellow_module_expression.pdf", width = 6, height = 6)
draw(ht3, heatmap_legend_side = "right")
dev.off()
#heatmap vertical list
ht_list <- ht2 %v% ht3
pdf(file = "Results/single_cell_adult_human/blue_yellow_heatmap_list.pdf", width = 6, height = 8)
draw(ht_list, ht_gap = unit(0.2, "cm"))
dev.off()

