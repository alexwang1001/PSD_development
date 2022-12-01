library("ComplexHeatmap")
library("circlize")
library("scales")
library("tidyverse")
library("flexclust")
library("SummarizedExperiment")

species_comparison <- readRDS("Results/species_comparison.rds")

###################################################################
#human sample
#human sample correlation
human_comparison <- species_comparison[, species_comparison$species == "human"]
human_comparison_coldata <- as.data.frame(colData(human_comparison))
human_comparison_coldata <- human_comparison_coldata %>% arrange(age)
human_comparison <- human_comparison[,human_comparison_coldata$label]
human_cor_mtx <- cor(assay(human_comparison),assay(human_comparison), method = "pearson")
write.csv(human_cor_mtx, "Results/species_comparison/human_cor_mtx.csv")

col_fun <-  colorRamp2(seq(0.4,1,0.6/9), c("white",RColorBrewer::brewer.pal(9, "YlGnBu")))
anno <- data.frame(Human_group = human_comparison$condition)
var <- anno$Human_group %>% unique()
cols <- RColorBrewer::brewer.pal(6, "Set2")
names(cols) <- var

row_ha_human <- rowAnnotation(df = anno,
                              col = list(Human_group = cols),
                              show_legend = FALSE,
                              show_annotation_name = FALSE)
col_ha_human <- HeatmapAnnotation(df = anno,
                                  col = list(Human_group = cols),
                                  show_legend = TRUE,
                                  show_annotation_name = FALSE)


ht_human <- Heatmap(human_cor_mtx,
                    name = "Pearson r",
                    border = TRUE,
                    col = col_fun,
                    column_title = "Human PSD",
                    cluster_rows = F,
                    row_dend_reorder = FALSE,
                    cluster_columns = F,
                    column_dend_reorder = FALSE,
                    show_row_dend = FALSE,
                    show_column_dend = FALSE,
                    show_row_names = FALSE,
                    show_column_names = FALSE,
                    left_annotation = row_ha_human,
                    top_annotation = col_ha_human,
                    heatmap_legend_param = list(border = "grey"),
                    width = unit(8, "cm"),
                    height = unit(8, "cm")
)
ht_human <- draw(ht_human)
pdf("Results/species_comparison/human_psd_correlation.pdf", width = 6, height = 6)
ht_human
dev.off()

#human euclidean distance
human_euclidean_mtx <- dist(t(assay(human_comparison)), upper = T)
human_euclidean_mtx <- as.matrix(human_euclidean_mtx)
write.csv(human_euclidean_mtx, "Results/species_comparison/human_euclidean_mtx.csv")

col_fun <-  colorRamp2(seq(40,100,60/8), rev(RColorBrewer::brewer.pal(9, "YlOrRd")))

anno <- data.frame(Human_group = human_comparison$condition)
var <- anno$Human_group %>% unique()
cols <- RColorBrewer::brewer.pal(6, "Set2")
names(cols) <- var

row_ha_human <- rowAnnotation(df = anno,
                              col = list(Human_group = cols),
                              show_legend = FALSE,
                              show_annotation_name = FALSE)
col_ha_human <- HeatmapAnnotation(df = anno,
                                  col = list(Human_group = cols),
                                  show_legend = TRUE,
                                  show_annotation_name = FALSE)


human_eucledian <- Heatmap(human_euclidean_mtx,
                           name = "Euclidean distance",
                           border = TRUE,
                           col = col_fun,
                           column_title = "Human PSD",
                           cluster_rows = F,
                           row_dend_reorder = FALSE,
                           cluster_columns = F,
                           column_dend_reorder = FALSE,
                           show_row_dend = FALSE,
                           show_column_dend = FALSE,
                           show_row_names = FALSE,
                           show_column_names = FALSE,
                           left_annotation = row_ha_human,
                           top_annotation = col_ha_human,
                           heatmap_legend_param = list(border = "grey"),
                           width = unit(8, "cm"),
                           height = unit(8, "cm")
)
human_eucledian <- draw(human_eucledian)

pdf("Results/species_comparison/human_psd_euclidean.pdf", width = 6, height = 6)
human_eucledian
dev.off()

###################################################################
#macaque sample
#macaque sample correlation
macaque_comparison <- species_comparison[, species_comparison$species == "macaque"]
macaque_comparison_coldata <- as.data.frame(colData(macaque_comparison))
macaque_comparison_coldata <- macaque_comparison_coldata %>% arrange(age)
macaque_comparison <- macaque_comparison[,macaque_comparison_coldata$label]
macaque_cor_mtx <- cor(assay(macaque_comparison),assay(macaque_comparison), method = "pearson")
write.csv(macaque_cor_mtx, "Results/species_comparison/macaque_cor_mtx.csv")

col_fun <-  colorRamp2(seq(0.4,1,0.6/9), c("white",RColorBrewer::brewer.pal(9, "YlGnBu")))
anno <- data.frame(Macaque_group = macaque_comparison$condition)
var <- anno$Macaque_group %>% unique()
cols <- RColorBrewer::brewer.pal(5, "Set2")
names(cols) <- var

row_ha_macaque <- rowAnnotation(df = anno,
                                col = list(Macaque_group = cols),
                                show_legend = FALSE,
                                show_annotation_name = FALSE)
row_ha_macaque2 <- rowAnnotation(df = anno,
                                col = list(Macaque_group = cols),
                                show_legend = TRUE,
                                show_annotation_name = FALSE)
col_ha_macaque <- HeatmapAnnotation(df = anno,
                                    col = list(Macaque_group = cols),
                                    show_legend = TRUE,
                                    show_annotation_name = FALSE)


ht_macaque <- Heatmap(macaque_cor_mtx,
                    name = "Pearson r",
                    border = TRUE,
                    col = col_fun,
                    column_title = "Macaque PSD",
                    cluster_rows = F,
                    row_dend_reorder = FALSE,
                    cluster_columns = F,
                    column_dend_reorder = FALSE,
                    show_row_dend = FALSE,
                    show_column_dend = FALSE,
                    show_row_names = FALSE,
                    show_column_names = FALSE,
                    left_annotation = row_ha_macaque,
                    top_annotation = col_ha_macaque,
                    heatmap_legend_param = list(border = "grey"),
                    width = unit(8, "cm"),
                    height = unit(8, "cm")
)
ht_macaque <- draw(ht_macaque)
pdf("Results/species_comparison/macaque_psd_correlation.pdf", width = 6, height = 6)
ht_macaque
dev.off()

###################################################################
#mouse sample
#mouse sample correlation
mouse_comparison <- species_comparison[, species_comparison$species == "mouse"]
mouse_comparison_coldata <- as.data.frame(colData(mouse_comparison))
mouse_comparison_coldata <- mouse_comparison_coldata %>% arrange(age)
mouse_comparison <- mouse_comparison[,mouse_comparison_coldata$label]
mouse_cor_mtx <- cor(assay(mouse_comparison),assay(mouse_comparison), method = "pearson")
write.csv(mouse_cor_mtx, "Results/species_comparison/mouse_cor_mtx.csv")

col_fun <-  colorRamp2(seq(0.4,1,0.6/9), c("white",RColorBrewer::brewer.pal(9, "YlGnBu")))
anno <- data.frame(Mouse_group = mouse_comparison$condition)
var <- anno$Mouse_group %>% unique()
cols <- RColorBrewer::brewer.pal(5, "Set2")
names(cols) <- var

row_ha_mouse <- rowAnnotation(df = anno,
                                col = list(Mouse_group = cols),
                                show_legend = FALSE,
                                show_annotation_name = FALSE)
row_ha_mouse2 <- rowAnnotation(df = anno,
                              col = list(Mouse_group = cols),
                              show_legend = TRUE,
                              show_annotation_name = FALSE)
col_ha_mouse <- HeatmapAnnotation(df = anno,
                                    col = list(Mouse_group = cols),
                                    show_legend = TRUE,
                                    show_annotation_name = FALSE)


ht_mouse <- Heatmap(mouse_cor_mtx,
                    name = "Pearson r",
                    border = TRUE,
                    col = col_fun,
                    column_title = "Mouse PSD",
                    cluster_rows = F,
                    row_dend_reorder = FALSE,
                    cluster_columns = F,
                    column_dend_reorder = FALSE,
                    show_row_dend = FALSE,
                    show_column_dend = FALSE,
                    show_row_names = FALSE,
                    show_column_names = FALSE,
                    left_annotation = row_ha_mouse,
                    top_annotation = col_ha_mouse,
                    heatmap_legend_param = list(border = "grey"),
                    width = unit(8, "cm"),
                    height = unit(8, "cm")
)
ht_mouse <- draw(ht_mouse)
pdf("Results/species_comparison/mouse_psd_correlation.pdf", width = 6, height = 6)
ht_mouse
dev.off()

###################################################################
#cross species comparison
#correlation between human and macaque
human_macaque_cor_mtx <- t(cor(assay(human_comparison), assay(macaque_comparison), method = "pearson"))
write.csv(human_macaque_cor_mtx, "Results/species_comparison/human_macaque_cor_mtx.csv")
col_fun <-  colorRamp2(seq(0.3,0.9,0.6/9), c("white",RColorBrewer::brewer.pal(9, "YlGnBu")))

ht_human_macaque_compare <- Heatmap(human_macaque_cor_mtx,
                                    name = "Pearson r",
                                    border = TRUE,
                                    col = col_fun,
                                    column_title = "Human PSD",
                                    row_title = "Macaque PSD",
                                    row_order = row_order(ht_macaque),
                                    column_order = column_order(ht_human),
                                    show_row_names = FALSE,
                                    show_column_names = FALSE,
                                    left_annotation = row_ha_macaque2,
                                    top_annotation = col_ha_human,
                                    heatmap_legend_param = list(border = "grey"),
                                    width = unit(8, "cm"),
                                    height = unit(8, "cm")
)
ht_human_macaque_compare <- draw(ht_human_macaque_compare)
pdf("Results/species_comparison/ht_human_macaque_compare.pdf", width = 6, height = 6)
ht_human_macaque_compare
dev.off()

col_fun <-  colorRamp2(c(0,1), c("white", "white"))
ht_human_macaque_compare_binary <- Heatmap(human_macaque_cor_mtx,
                                    name = "Pearson r",
                                    border = TRUE,
                                    col = col_fun,
                                    cell_fun = function(j, i, x, y, width, height, fill) {
                                      if (human_macaque_cor_mtx[i,j] > 0.6) {
                                        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "orange",lwd = 0, fill = "orange"))
                                      }
                                    },
                                    column_title = "Human PSD",
                                    row_title = "Macaque PSD",
                                    row_order = row_order(ht_macaque),
                                    column_order = column_order(ht_human),
                                    show_row_names = FALSE,
                                    show_column_names = FALSE,
                                    left_annotation = row_ha_macaque2,
                                    top_annotation = col_ha_human,
                                    show_heatmap_legend = FALSE,
                                    width = unit(8, "cm"),
                                    height = unit(8, "cm"),
                                    
)
ht_human_macaque_compare_binary <- draw(ht_human_macaque_compare_binary)
pdf("Results/species_comparison/ht_human_macaque_compare_binary.pdf", width = 6, height = 6)
ht_human_macaque_compare_binary
dev.off()



#correlation between human and mouse
human_mouse_cor_mtx <- t(cor(assay(human_comparison), assay(mouse_comparison), method = "pearson"))
write.csv(human_mouse_cor_mtx, "Results/species_comparison/human_mouse_cor_mtx.csv")
col_fun <-  colorRamp2(seq(0.3,0.9,0.6/9), c("white",RColorBrewer::brewer.pal(9, "YlGnBu")))

ht_human_mouse_compare <- Heatmap(human_mouse_cor_mtx,
                                    name = "Pearson r",
                                    border = TRUE,
                                    col = col_fun,
                                    column_title = "Human PSD",
                                    row_title = "Mouse PSD",
                                    row_order = row_order(ht_mouse),
                                    column_order = column_order(ht_human),
                                    show_row_names = FALSE,
                                    show_column_names = FALSE,
                                    left_annotation = row_ha_mouse2,
                                    top_annotation = col_ha_human,
                                    heatmap_legend_param = list(border = "grey"),
                                    width = unit(8, "cm"),
                                    height = unit(8, "cm")
)
ht_human_mouse_compare <- draw(ht_human_mouse_compare)
pdf("Results/species_comparison/ht_human_mouse_compare.pdf", width = 6, height = 6)
ht_human_mouse_compare
dev.off()

col_fun <-  colorRamp2(c(0,1), c("white", "white"))
ht_human_mouse_compare_binary <- Heatmap(human_mouse_cor_mtx,
                                         name = "Pearson r",
                                         border = TRUE,
                                         col = col_fun,
                                         cell_fun = function(j, i, x, y, width, height, fill) {
                                           if (human_mouse_cor_mtx[i,j] > 0.6) {
                                             grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "orange",lwd = 0, fill = "orange"))
                                           }
                                         },
                                         column_title = "Human PSD",
                                         row_title = "Mouse PSD",
                                         row_order = row_order(ht_mouse),
                                         column_order = column_order(ht_human),
                                         show_row_names = FALSE,
                                         show_column_names = FALSE,
                                         left_annotation = row_ha_mouse2,
                                         top_annotation = col_ha_human,
                                         show_heatmap_legend = FALSE,
                                         width = unit(8, "cm"),
                                         height = unit(8, "cm")
)
ht_human_mouse_compare_binary <- draw(ht_human_mouse_compare_binary)
pdf("Results/species_comparison/ht_human_mouse_compare_binary.pdf", width = 6, height = 6)
ht_human_mouse_compare_binary
dev.off()



#correlation between macaque and mouse
macaque_mouse_cor_mtx <- t(cor(assay(macaque_comparison), assay(mouse_comparison), method = "pearson"))
write.csv(macaque_mouse_cor_mtx, "Results/species_comparison/macaque_mouse_cor_mtx.csv")
col_fun <-  colorRamp2(seq(0.3,0.9,0.6/9), c("white",RColorBrewer::brewer.pal(9, "YlGnBu")))

ht_macaque_mouse_compare <- Heatmap(macaque_mouse_cor_mtx,
                                  name = "Pearson r",
                                  border = TRUE,
                                  col = col_fun,
                                  column_title = "Macaque PSD",
                                  row_title = "Mouse PSD",
                                  row_order = row_order(ht_mouse),
                                  column_order = column_order(ht_macaque),
                                  show_row_names = FALSE,
                                  show_column_names = FALSE,
                                  left_annotation = row_ha_mouse2,
                                  top_annotation = col_ha_macaque,
                                  heatmap_legend_param = list(border = "grey"),
                                  width = unit(8, "cm"),
                                  height = unit(8, "cm")
)
ht_macaque_mouse_compare <- draw(ht_macaque_mouse_compare)
pdf("Results/species_comparison/ht_macaque_mouse_compare.pdf", width = 6, height = 6)
ht_macaque_mouse_compare
dev.off()

col_fun <-  colorRamp2(c(0,1), c("white", "white"))
ht_macaque_mouse_compare_binary <- Heatmap(macaque_mouse_cor_mtx,
                                         name = "Pearson r",
                                         border = TRUE,
                                         col = col_fun,
                                         cell_fun = function(j, i, x, y, width, height, fill) {
                                           if (macaque_mouse_cor_mtx[i,j] > 0.6) {
                                             grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "orange",lwd = 0, fill = "orange"))
                                           }
                                         },
                                         column_title = "Macaque PSD",
                                         row_title = "Mouse PSD",
                                         row_order = row_order(ht_mouse),
                                         column_order = column_order(ht_macaque),
                                         show_row_names = FALSE,
                                         show_column_names = FALSE,
                                         left_annotation = row_ha_mouse2,
                                         top_annotation = col_ha_macaque,
                                         show_heatmap_legend = FALSE,
                                         width = unit(8, "cm"),
                                         height = unit(8, "cm")
)
ht_macaque_mouse_compare_binary <- draw(ht_macaque_mouse_compare_binary)
pdf("Results/species_comparison/ht_macaque_mouse_compare_binary.pdf", width = 6, height = 6)
ht_macaque_mouse_compare_binary
dev.off()
