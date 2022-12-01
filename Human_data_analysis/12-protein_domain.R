library("tidyverse")
library("DEP")
library("ComplexHeatmap")
library("circlize")
library("SummarizedExperiment")

load(file = "Results/dep_corrected_specificity.RData")
write.table(rowData(dep_corrected_specificity)$Protein.IDs, file = "SMART_domain_analysis/uniprot_IDs.txt",sep = "\t", row.names = FALSE, quote = FALSE)

SMART <- read.csv("SMART_domain_analysis/domain_matrix_smart_only.csv")
colnames(SMART)[1] <- "Protein.IDs"

rowdata <- data.frame(rowData(dep_corrected_specificity))[,c(1,2)]
colnames(rowdata)[1] <- "SYMBOL"

geneIDsModules <- read.csv(file = "./Results/WGCNA_ORA/geneIDsModules.csv", row.names = 1)

rowdata <- left_join(rowdata, geneIDsModules, by = "SYMBOL")

SMART2 <- left_join(rowdata, SMART, by = "Protein.IDs")
SMART2$geneModules <- factor(SMART2$geneModules, levels = c("brown", "blue", "turquoise", "yellow", "grey"))

write.csv(SMART2, "Results/protein_domains/SMART_result_summary.csv", row.names = FALSE)

SMART3 <- SMART2 %>% filter(geneModules != "grey")
dat <- SMART3[,5:421]
dat2 <- dat
dat2[dat2 != 0] <- 1
table(colSums(dat2))
#only analyze domains that are present more than 6 times in all PSD proteins
dat3 <- dat2[,colSums(dat2) > 6]
rownames(dat3) <- SMART3$SYMBOL
sum(rowSums(dat3) != 0)

dat4 <- dat3
dat4_modules <- SMART3$geneModules


anno <- data.frame(as.character(dat4_modules))
colnames(anno) <- "Modules"
var <- anno$Modules %>% unique() %>% sort()
cols <- var
names(cols) <- var


col_ha <- HeatmapAnnotation(df = anno,
                            col = list(Modules = cols),
                            show_legend = FALSE,
                            show_annotation_name = TRUE,
                            annotation_name_side = "left")


jaccard_dist = function(x, y) {
  s = 1 - sum(x & y)/sum(x | y)
  if(is.na(s)) {
    return(1)
  } else {
    return(s)
  }
}

ht1 <- Heatmap(t(dat4),
               col = c("gray90", "blue"),
               clustering_distance_columns = jaccard_dist,
               clustering_distance_rows = jaccard_dist,
               column_split = dat4_modules,
               column_title = NULL,
               cluster_column_slices = FALSE,
               show_row_dend = TRUE,
               row_names_side = "left",
               show_column_dend = FALSE,
               show_column_names = FALSE,
               show_heatmap_legend = TRUE,
               border = TRUE,
               column_gap = unit(0, "mm"),
               top_annotation = col_ha,
               heatmap_legend_param = list(
                 title = c("Protein contains the domain"),
                 color_bar = "discrete",
                 at = c(1,0),
                 labels = c("Yes", "No"),
                 legend_width = unit(3, "cm"),
                 legend_direction = "horizontal"),
               width = unit(18, "cm"))

pdf("Results/protein_domains/domain_individual.pdf", width = 12, height = 11)
draw(ht1, heatmap_legend_side = "top",height = unit(22,"cm"))
dev.off()

#module level summary comparison
result <- data.frame()
for (i in c("brown", "blue", "turquoise", "yellow")) {
  sub_dat3 <- dat3[SMART3$geneModules == i,]
  sum <- data.frame(colSums(sub_dat3))
  normalized_sum <- sum/nrow(sub_dat3)*100
  colnames(normalized_sum) <- as.character(i)
  if (ncol(result) == 0) {
    result <- normalized_sum
  } else {
    result <- cbind(result, normalized_sum)
  }
}


ha <- HeatmapAnnotation(Modules = letters[1:4],
                        col = list(Modules = c("a" = "brown", "b" = "blue", "c" = "turquoise", "d" = "gold")),
                        show_legend = FALSE,
                        annotation_name_side = "left")

RColorBrewer::brewer.pal(9, "Reds")
col_fun = colorRamp2(seq(0,6,6/9), c("white",RColorBrewer::brewer.pal(9, "YlGnBu")))

ht2 <- Heatmap(result,
        col = col_fun,
        clustering_distance_rows = "pearson",
        show_row_names = TRUE,
        cluster_columns = FALSE,
        show_column_names = TRUE,
        show_heatmap_legend = TRUE,
        border = TRUE,
        top_annotation = ha,
        heatmap_legend_param = list(
          title = "Percentage of proteins carrying the domain",
          at = c(0,2,4,6),
          border = "gray",
          legend_width = unit(5, "cm"),
          legend_direction = "horizontal",
          legend_height = unit(5,"cm")),
        width = unit(5, "cm")
        )
pdf("Results/protein_domains/domain_summary.pdf", width = 12, height = 13)
draw(ht2, heatmap_legend_side = "top",height = unit(25,"cm"))
dev.off()


