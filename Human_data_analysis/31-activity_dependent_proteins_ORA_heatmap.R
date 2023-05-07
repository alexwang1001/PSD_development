library("tidyverse")
library("ComplexHeatmap")
library("circlize")

load("Results/activity_dependent_proteins_ORA.RData")
res <- as.data.frame(activity_dependent_proteins_ORA@compareClusterResult)
odds <- res %>% dplyr::select(Cluster, Description, OddsRatio) %>% spread(key = Cluster, value = OddsRatio)
odds2 <- odds[,c(5,3,2,4,6)]
rownames(odds2) <- odds$Description

pvalue <- res %>% mutate (`-log10(p.adj)` = -log10(p.adjust)) %>% dplyr::select(Cluster, Description, `-log10(p.adj)`)
mat <- spread(pvalue, key = Cluster, value = `-log10(p.adj)`)
mat2 <- mat[,c(5,3,2,4)]
rownames(mat2) <- mat$Description

ha <- HeatmapAnnotation(Modules = letters[1:4],
                        col = list(Modules = c("a" = "brown", "b" = "blue", "c" = "turquoise", "d" = "gold")),
                        show_legend = FALSE,
                        annotation_name_side = "left")


RColorBrewer::brewer.pal(9, "Reds")
col_fun = colorRamp2(seq(0,6,6/9), c("white",RColorBrewer::brewer.pal(9, "Reds")))
ht1 <- Heatmap(mat2,
               name = "-log10(p.adj)",
               col = col_fun, 
               border = FALSE,
               rect_gp = gpar(col = "gray", lwd = 1),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (mat2[i,j] > -log10(0.0001)) {
                   grid.text(round(odds2[i, j], digits = 1), x, y, gp = gpar(fontsize = 12, col= "white"))
                 } else {
                   grid.text(round(odds2[i, j], digits = 1), x, y, gp = gpar(fontsize = 12, col= "black"))
                 }
                 if (mat2[i,j] > -log10(0.05)) {
                   grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "steelblue2",lwd = 3, fill = NA))
                 }
               },
               cluster_rows = TRUE,
               cluster_columns = FALSE,
               clustering_distance_rows = "pearson",
               clustering_method_rows = "complete",
               show_column_names = TRUE,
               top_annotation = ha,
               heatmap_legend_param = list(
                 title = "p.adj",
                 title_position = "lefttop",
                 at = c(0,2,4,6),
                 labels = c(1, expression(10^-2), expression(10^-4), expression(10^-6)),
                 border = "gray",
                 legend_width = unit(5, "cm"),
                 legend_direction = "horizontal"
               )
)

mat3 <- as.matrix(mat$all)
rownames(mat3) <- mat$Description
colnames(mat3) <- "all"
ht2 <- Heatmap(mat3,
               col = col_fun, 
               border = FALSE,
               rect_gp = gpar(col = "gray", lwd = 1),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (mat3[i,j] > -log10(0.0001)) {
                   grid.text(round(data.frame(odds$all)[i, j], digits = 1), x, y, gp = gpar(fontsize = 12, col= "white"))
                 } else {
                   grid.text(round(data.frame(odds$all)[i, j], digits = 1), x, y, gp = gpar(fontsize = 12, col= "black"))
                 }
                 if (mat3[i,j] > -log10(0.05)) {
                   grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "steelblue2",lwd = 3, fill = NA))
                 }
               },
               cluster_columns = FALSE,
               clustering_distance_rows = "pearson",
               clustering_method_rows = "complete",
               show_column_names = TRUE,
               show_heatmap_legend = FALSE
)

ht_list <- ht1+ht2

pdf(file = "Results/activity_dependent_proteome_changes/activity_dependent_proteins_ORA.pdf", width = 4.8, height = 2.8)
draw(ht_list, heatmap_legend_side = "top", ht_gap = unit(0.5, "cm"))
dev.off()

