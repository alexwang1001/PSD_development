library("tidyverse")
library("ComplexHeatmap")
library("circlize")
library("RColorBrewer")

res <- read.csv(file = "Results/NDD_mutations/NDD_de_novo_mutations_burden.csv")

res <- res %>% filter(!(disorder %in% c("control_LGD","control_missense"))) %>% mutate (`-log10(p.adj)` = -log10(p.adj)) %>% dplyr::select(module, disorder, odds_ratio,`-log10(p.adj)`)

odds <- res %>% select(module, disorder, odds_ratio) %>% spread(key = module, value = odds_ratio)
odds2 <- odds[,c(4,3,5,6)]
rownames(odds2) <- odds$disorder

mat <- res %>% select(module, disorder, `-log10(p.adj)`) %>% spread(key = module, value = `-log10(p.adj)`)
mat2 <- mat[,c(4,3,5,6)]
rownames(mat2) <- mat$disorder

ha <- HeatmapAnnotation(Modules = letters[1:4],
                        col = list(Modules = c("a" = "brown", "b" = "blue", "c" = "turquoise", "d" = "gold")),
                        show_legend = FALSE,
                        annotation_name_side = "left")

RColorBrewer::brewer.pal(9, "Reds")
col_fun <-  colorRamp2(seq(0,6,6/9), c("white",RColorBrewer::brewer.pal(9, "Reds")))


ht1 <- Heatmap(mat2,
               name = "-log10(p.adj)",
               col = col_fun, 
               row_title = NULL,
               border = FALSE,
               rect_gp = gpar(col = "gray", lwd = 1),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (mat2[i,j] > -log10(0.001)) {
                   grid.text(round(odds2[i, j], digits = 1), x, y, gp = gpar(fontsize = 12, col= "white"))
                 } else {
                   grid.text(round(odds2[i, j], digits = 1), x, y, gp = gpar(fontsize = 12, col= "black"))
                 }
                 if (mat2[i,j] > -log10(0.05)) {
                   grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "steelblue2",lwd = 3, fill = NA))
                 }
               },
               cluster_columns = FALSE,
               clustering_distance_rows = "pearson",
               clustering_method_rows = "complete",
               show_column_names = TRUE,
               top_annotation = ha,
               heatmap_legend_param = list(
                 title = "p.adj",
                 title_position = "lefttop",
                 at = c(0,2,4,6),
                 labels = c(1, expression(10^-2), expression(10^-4), expression(""<=10^-6)),
                 border = "gray",
                 legend_width = unit(5, "cm"),
                 legend_direction = "horizontal"
               )
)

ht1

mat3 <- as.matrix(mat$all)
rownames(mat3) <- mat$disorder
colnames(mat3) <- "all"

ht2 <- Heatmap(mat3,
               col = col_fun, 
               row_title = NULL,
               border = FALSE,
               rect_gp = gpar(col = "gray", lwd = 1),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (mat3[i,j] > -log10(0.001)) {
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
draw(ht_list, heatmap_legend_side = "top", ht_gap = unit(0.5, "cm"))

pdf(file = "Results/NDD_mutations/NDD_mutations_burden_heatmap.pdf", width = 6.5, height = 5.5)
draw(ht_list, heatmap_legend_side = "top", ht_gap = unit(0.5, "cm"))
dev.off()

