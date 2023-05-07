library("tidyverse")
library("ComplexHeatmap")
library("circlize")
library("RColorBrewer")

load(file = "Results/NDD_diff_ORA.RData")

res <- NDD_diff_ORA@compareClusterResult
odds <- res %>% dplyr::select(Cluster, Description, OddsRatio) %>% spread(key = Cluster, value = OddsRatio)
odds2 <- odds[,c(5,3,2,4)]
rownames(odds2) <- odds$Description

pvalue <- res %>% mutate (`-log10(p.adj)` = -log10(p.adjust)) %>% dplyr::select(Cluster, Description, `-log10(p.adj)`)
mat <- spread(pvalue, key = Cluster, value = `-log10(p.adj)`)
mat2 <- mat[,c(5,3,2,4)]
rownames(mat2) <- c("Down_in_ASD", "Up_in_ASD", "Down_in_BPD", "Up_in_BPD", "Down_in_MDD", "Up_in_MDD", "Down_in_SCZ", "Up_in_SCZ")

ha <- HeatmapAnnotation(Modules = letters[1:4],
                        col = list(Modules = c("a" = "brown", "b" = "blue", "c" = "turquoise", "d" = "gold")),
                        show_legend = FALSE,
                        annotation_name_side = "left")


RColorBrewer::brewer.pal(9, "Reds")
col_fun = colorRamp2(seq(0,9,9/9), c("white",RColorBrewer::brewer.pal(9, "Reds")))
ht1 <- Heatmap(mat2,
        name = "-log10(p.adj)",
        col = col_fun, 
        border = FALSE,
        rect_gp = gpar(col = "gray", lwd = 1),
        cell_fun = function(j, i, x, y, width, height, fill) {
                if (mat2[i,j] > -log10(0.00001)) {
                        grid.text(round(odds2[i, j], digits = 1), x, y, gp = gpar(fontsize = 12, col= "white"))
                } else {
                        grid.text(round(odds2[i, j], digits = 1), x, y, gp = gpar(fontsize = 12, col= "black"))
                }
                if (mat2[i,j] > -log10(0.05)) {
                        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "steelblue2",lwd = 3, fill = NA))
                }
        },
        cluster_columns = FALSE,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "complete",
        show_column_names = TRUE,
        row_split = rep(c("Down", "Up"),4),
        row_gap = unit(0.5, "cm"),
        top_annotation = ha,
        heatmap_legend_param = list(
                title = "p.adj",
                title_position = "lefttop",
                at = c(0,3,6,9),
                labels = c(1, expression(10^-3), expression(10^-6), expression(""<=10^-9)),
                border = "gray",
                legend_width = unit(5, "cm"),
                legend_direction = "horizontal"
        )
)

mat3 <- as.matrix(mat$all)
rownames(mat3) <- c("Down_in_ASD", "Up_in_ASD", "Down_in_BPD", "Up_in_BPD", "Down_in_MDD", "Up_in_MDD", "Down_in_SCZ", "Up_in_SCZ")
colnames(mat3) <- "all"
ht2 <- Heatmap(mat3,
               col = col_fun, 
               border = FALSE,
               rect_gp = gpar(col = "gray", lwd = 1),
               cell_fun = function(j, i, x, y, width, height, fill) {
                       if (mat3[i,j] > -log10(0.00001)) {
                               grid.text(round(data.frame(odds$all)[i, j], digits = 1), x, y, gp = gpar(fontsize = 12, col= "white"))
                       } else {
                               grid.text(round(data.frame(odds$all)[i, j], digits = 1), x, y, gp = gpar(fontsize = 12, col= "black"))
                       }
                       if (mat3[i,j] > -log10(0.05)) {
                               grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "steelblue2",lwd = 3, fill = NA))
                       }
               },
               cluster_columns = FALSE,
               clustering_distance_rows = "euclidean",
               clustering_method_rows = "complete",
               show_column_names = TRUE,
               show_heatmap_legend = FALSE
)

ht_list <- ht1+ht2
draw(ht_list, heatmap_legend_side = "top", ht_gap = unit(0.5, "cm"))

pdf(file = "Results/differential_expression_in_psychiatric_disorder/NDD_diff_ORA_heatmap.pdf", width = 5, height = 4.5)
draw(ht_list, heatmap_legend_side = "top", ht_gap = unit(0.5, "cm"))
dev.off()
