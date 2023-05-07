library(ComplexHeatmap)
library(circlize)
library(scales)
library(tidyverse)
library(flexclust)
library(SummarizedExperiment)

species_comparison_V1 <- readRDS("Results/species_comparison_V1.rds")

###################################################################
#human PFC sample correlation
human_PFC_comparison <- species_comparison_V1[, species_comparison_V1$species == "Human PFC"]
human_PFC_comparison_coldata <- as.data.frame(colData(human_PFC_comparison))
human_PFC_comparison_coldata <- human_PFC_comparison_coldata %>% arrange(age)
human_PFC_comparison <- human_PFC_comparison[,human_PFC_comparison_coldata$label]
human_PFC_cor_mtx <- cor(assay(human_PFC_comparison),assay(human_PFC_comparison), method = "pearson")
write.csv(human_PFC_cor_mtx, "Results/species_comparison/human_PFC_cor_mtx.csv")

col_fun <-  colorRamp2(seq(0.4,1,0.6/9), c("white",RColorBrewer::brewer.pal(9, "YlGnBu")))
anno <- data.frame(Human_group = human_PFC_comparison$condition)
var <- anno$Human_group %>% unique()
cols <- RColorBrewer::brewer.pal(6, "Set2")
names(cols) <- var

row_ha_human_PFC <- rowAnnotation(df = anno,
                              col = list(Human_group = cols),
                              show_legend = FALSE,
                              show_annotation_name = FALSE)
col_ha_human_PFC <- HeatmapAnnotation(df = anno,
                                  col = list(Human_group = cols),
                                  show_legend = TRUE,
                                  show_annotation_name = FALSE)


ht_human_PFC <- Heatmap(human_PFC_cor_mtx,
                    name = "Pearson r",
                    border = TRUE,
                    col = col_fun,
                    column_title = "Human PFC PSD",
                    cluster_rows = F,
                    row_dend_reorder = FALSE,
                    cluster_columns = F,
                    column_dend_reorder = FALSE,
                    show_row_dend = FALSE,
                    show_column_dend = FALSE,
                    show_row_names = FALSE,
                    show_column_names = FALSE,
                    left_annotation = row_ha_human_PFC,
                    top_annotation = col_ha_human_PFC,
                    heatmap_legend_param = list(border = "grey"),
                    width = unit(8, "cm"),
                    height = unit(8, "cm")
)
ht_human_PFC <- draw(ht_human_PFC)
pdf("Results/species_comparison/human_PFC_correlation.pdf", width = 6, height = 6)
ht_human_PFC
dev.off()

####################################
#human V1 sample correlation
human_V1_comparison <- species_comparison_V1[, species_comparison_V1$species == "Human V1"]
human_V1_comparison_coldata <- as.data.frame(colData(human_V1_comparison))
human_V1_comparison_coldata <- human_V1_comparison_coldata %>% arrange(age)
human_V1_comparison <- human_V1_comparison[,human_V1_comparison_coldata$label]
human_V1_cor_mtx <- cor(assay(human_V1_comparison),assay(human_V1_comparison), method = "pearson")
write.csv(human_V1_cor_mtx, "Results/species_comparison/human_V1_cor_mtx.csv")

col_fun <-  colorRamp2(seq(0.4,1,0.6/9), c("white",RColorBrewer::brewer.pal(9, "YlGnBu")))
anno <- data.frame(Human_group = human_V1_comparison$condition)
var <- anno$Human_group %>% unique()
cols <- RColorBrewer::brewer.pal(6, "Set2")
names(cols) <- var

row_ha_human_V1 <- rowAnnotation(df = anno,
                              col = list(Human_group = cols),
                              show_legend = FALSE,
                              show_annotation_name = FALSE)
col_ha_human_V1 <- HeatmapAnnotation(df = anno,
                                  col = list(Human_group = cols),
                                  show_legend = TRUE,
                                  show_annotation_name = FALSE)


ht_human_V1 <- Heatmap(human_V1_cor_mtx,
                        name = "Pearson r",
                        border = TRUE,
                        col = col_fun,
                        column_title = "Human V1 PSD",
                        cluster_rows = F,
                        row_dend_reorder = FALSE,
                        cluster_columns = F,
                        column_dend_reorder = FALSE,
                        show_row_dend = FALSE,
                        show_column_dend = FALSE,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        left_annotation = row_ha_human_V1,
                        top_annotation = col_ha_human_V1,
                        heatmap_legend_param = list(border = "grey"),
                        width = unit(8, "cm"),
                        height = unit(8, "cm")
)
ht_human_V1 <- draw(ht_human_V1)
pdf("Results/species_comparison/human_V1_correlation.pdf", width = 6, height = 6)
ht_human_V1
dev.off()
###################################################################
#macaque sample
#macaque sample correlation
macaque_comparison <- species_comparison_V1[, species_comparison_V1$species == "Macaque"]
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
mouse_comparison <- species_comparison_V1[, species_comparison_V1$species == "Mouse"]
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
#human area comparison
PFC_V1_cor_mtx <- t(cor(assay(human_PFC_comparison), assay(human_V1_comparison), method = "pearson"))
write.csv(PFC_V1_cor_mtx, "Results/species_comparison/PFC_V1_cor_mtx.csv")
col_fun <-  colorRamp2(seq(0.3,0.9,0.6/9), c("white",RColorBrewer::brewer.pal(9, "YlGnBu")))

ht_PFC_V1_compare <- Heatmap(PFC_V1_cor_mtx,
                             name = "Pearson r",
                             border = TRUE,
                             col = col_fun,
                             column_title = "Human PFC PSD",
                             row_title = "Human V1 PSD",
                             row_order = row_order(ht_human_V1),
                             column_order = column_order(ht_human_PFC),
                             show_row_names = FALSE,
                             show_column_names = FALSE,
                             left_annotation = row_ha_human_V1,
                             top_annotation = col_ha_human_PFC,
                             heatmap_legend_param = list(border = "grey"),
                             width = unit(8, "cm"),
                             height = unit(8, "cm")
)
ht_PFC_V1_compare <- draw(ht_PFC_V1_compare)
pdf("Results/species_comparison/ht_human_PFC_V1_compare.pdf", width = 6, height = 6)
ht_PFC_V1_compare
dev.off()

#cross species comparison
#correlation between human V1 and macaque
human_macaque_cor_mtx <- t(cor(assay(human_V1_comparison), assay(macaque_comparison), method = "pearson"))
write.csv(human_macaque_cor_mtx, "Results/species_comparison/human_macaque_cor_mtx.csv")
col_fun <-  colorRamp2(seq(0.3,0.9,0.6/9), c("white",RColorBrewer::brewer.pal(9, "YlGnBu")))

ht_human_macaque_compare <- Heatmap(human_macaque_cor_mtx,
                                    name = "Pearson r",
                                    border = TRUE,
                                    col = col_fun,
                                    column_title = "Human V1 PSD",
                                    row_title = "Macaque PSD",
                                    row_order = row_order(ht_macaque),
                                    column_order = column_order(ht_human_V1),
                                    show_row_names = FALSE,
                                    show_column_names = FALSE,
                                    left_annotation = row_ha_macaque2,
                                    top_annotation = col_ha_human_V1,
                                    heatmap_legend_param = list(border = "grey"),
                                    width = unit(8, "cm"),
                                    height = unit(8, "cm")
)
ht_human_macaque_compare <- draw(ht_human_macaque_compare)
pdf("Results/species_comparison/ht_human_macaque_compare.pdf", width = 6, height = 6)
ht_human_macaque_compare
dev.off()

#correlation between human V1 and mouse
human_mouse_cor_mtx <- t(cor(assay(human_V1_comparison), assay(mouse_comparison), method = "pearson"))
write.csv(human_mouse_cor_mtx, "Results/species_comparison/human_mouse_cor_mtx.csv")
col_fun <-  colorRamp2(seq(0.3,0.9,0.6/9), c("white",RColorBrewer::brewer.pal(9, "YlGnBu")))


ht_human_mouse_compare <- Heatmap(human_mouse_cor_mtx,
                                    name = "Pearson r",
                                    border = TRUE,
                                    col = col_fun,
                                    column_title = "Human V1 PSD",
                                    row_title = "Mouse PSD",
                                    row_order = row_order(ht_mouse),
                                    column_order = column_order(ht_human_V1),
                                    show_row_names = FALSE,
                                    show_column_names = FALSE,
                                    left_annotation = row_ha_mouse2,
                                    top_annotation = col_ha_human_V1,
                                    heatmap_legend_param = list(border = "grey"),
                                    width = unit(8, "cm"),
                                    height = unit(8, "cm")
)
ht_human_mouse_compare <- draw(ht_human_mouse_compare)
pdf("Results/species_comparison/ht_human_mouse_compare.pdf", width = 6, height = 6)
ht_human_mouse_compare
dev.off()
