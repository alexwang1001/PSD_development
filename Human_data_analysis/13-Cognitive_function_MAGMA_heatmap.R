library("tidyverse")
library("readxl")
library("ComplexHeatmap")
library("circlize")
library("RColorBrewer")
 
res <- read_excel(path = "GWAS/cognitive_function/Cognitive_GWAS_module_magma.xlsx")

res_list <- res %>% group_split(module)
res2 <- c()
for (i in 1:5) {
  tmp <- cbind(res_list[[i]],p.adjust( res_list[[i]]$p_value, method = "BH"))
  res2 <- rbind(res2,tmp)
}

colnames(res2)[5] <- "p.adj"

res3 <- res2 %>% mutate(`-log10(p.adj)` = -log10(p.adj)) %>% dplyr::select(module, domain, `-log10(p.adj)`)
mat <- spread(res3, key = module, value = `-log10(p.adj)`)
mat2 <- mat[,c(4,3,5,6)]
rownames(mat2) <- mat$domain

res4 <- res2 %>% dplyr::select(module, domain, beta)
beta <- spread(res4, key = module, value = beta)
beta2 <- beta[,c(4,3,5,6)]
rownames(beta2) <- beta2$domain


ha <- HeatmapAnnotation(Modules = letters[1:4],
                        col = list(Modules = c("a" = "brown", "b" = "blue", "c" = "turquoise", "d" = "gold")),
                        show_legend = FALSE,
                        annotation_name_side = "left")


RColorBrewer::brewer.pal(9, "Reds")
col_fun = colorRamp2(seq(0,3,3/9), c("white",RColorBrewer::brewer.pal(9, "Reds")))
ht1 <- Heatmap(mat2,
               name = "-log10(p.adj)",
               col = col_fun, 
               border = FALSE,
               rect_gp = gpar(col = "grey", lwd = 1),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (mat2[i,j] > -log10(0.01)) {
                   grid.text(round(beta2[i, j], digits = 2), x, y, gp = gpar(fontsize = 20, col= "white"))
                 } else {
                   grid.text(round(beta2[i, j], digits = 2), x, y, gp = gpar(fontsize = 20, col= "black"))
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
                 at = c(0,1,2,3),
                 labels = c(1, expression(10^-1), expression(10^-2), expression(10^-3)),
                 border = "gray",
                 legend_width = unit(5, "cm"),
                 legend_direction = "horizontal"
               )
)

mat3 <- as.matrix(mat$all)
rownames(mat3) <- mat$domain
colnames(mat3) <- "all"

beta3 <- as.matrix(beta$all)
rownames(beta3) <- beta$domain
colnames(beta3) <- "all"


ht2 <- Heatmap(mat3,
               col = col_fun, 
               border = FALSE,
               rect_gp = gpar(col = "grey", lwd = 1),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (mat3[i,j] > -log10(0.01)) {
                   grid.text(round(beta3[i, j], digits = 2), x, y, gp = gpar(fontsize = 20, col= "white"))
                 } else {
                   grid.text(round(beta3[i, j], digits = 2), x, y, gp = gpar(fontsize = 20, col= "black"))
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
draw(ht_list,heatmap_legend_side = "top", ht_gap = unit(0.5, "cm"))

pdf(file = "Results/cognitive_function/cognitive_function_GWAS.pdf", width = 7, height = 4.5)
draw(ht_list,heatmap_legend_side = "top", ht_gap = unit(0.5, "cm"))
dev.off()
