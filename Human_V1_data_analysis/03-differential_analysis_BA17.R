library("DEP")
library("SummarizedExperiment")
library("ComplexHeatmap")
library("limma")
library("vsn")
library("tidyverse")
library("scales")
library("patchwork")
library("seriation")

data_imp_all <- readRDS("Results/data_imp_all.rds")
data_imp_BA17 <- data_imp_all[,data_imp_all$BrainArea %in% c("unknown", "BA 17") & data_imp_all$Disorder == "Control"]
colData(data_imp_BA17)<- droplevels(colData(data_imp_BA17))

#using limma
condition <- as.factor(data_imp_BA17$condition)
degradation <- as.factor(data_imp_BA17$Degradation)
design <- model.matrix(terms(~ condition + degradation, keep.order = TRUE))
colnames(design)
colnames(design)[2:7] <- c("GW22_23","GW28_41", "Year00_01", "Year04_08", "Year18_22", "degradation")
colnames(design)
fit1 <- lmFit(assay(data_imp_BA17), design=design)
fit1 <- eBayes(fit1)
head(coef(fit1))
colnames(coef(fit1))
#save comparison data
cntrst <- colnames(design)[2:6]
retrieve_fun <- function(comp, fit) {
  res <- topTable(fit, sort.by = "p", coef = comp, number = Inf, 
                  confint = TRUE)
  res <- res[!is.na(res$P.Value), ]
  res$comparison <- rep(comp , dim(res)[1])
  res <- rownames_to_column(res)
  return(res)
}

limma_res <- map_df(cntrst, retrieve_fun, fit = fit1)
table <- limma_res %>% dplyr::select(rowname, logFC, CI.L, CI.R, t, 
                              P.Value, adj.P.Val, comparison) %>% mutate(comparison = paste0(comparison,
                                                                                             "_vs_GW18_19")) %>% gather(variable, value, -c(rowname, 
                                                                                                                                            comparison)) %>% mutate(variable = recode(variable, logFC = "diff", 
                                                                                                                                                                                      P.Value = "p.val", adj.P.Val = "p.adj")) %>% unite(temp, comparison, 
                                                                                                                                                                                                                                         variable) %>% spread(temp, value)


rowData(data_imp_BA17) <- merge(rowData(data_imp_BA17, use.names = FALSE), table, 
                            by.x = "name", by.y = "rowname", all.x = TRUE, sort = FALSE)



#make all pairwise comparison
cntrst2 <- apply(utils::combn(colnames(design)[6:2] , 2), 2, paste, 
                 collapse = " - ")
made_contrasts <- makeContrasts(contrasts = cntrst2, levels = design)
contrast_fit <- contrasts.fit(fit1, made_contrasts)
eB_fit <- eBayes(contrast_fit)
limma_res2 <- map_df(cntrst2, retrieve_fun, fit = eB_fit)
table2 <- limma_res2 %>% dplyr::select(rowname, logFC, CI.L, CI.R, t, 
                                P.Value, adj.P.Val, comparison) %>% mutate(comparison = gsub(" - ", "_vs_",
                                                                                             comparison)) %>% gather(variable, value, -c(rowname, 
                                                                                                                                         comparison)) %>% mutate(variable = recode(variable, logFC = "diff", 
                                                                                                                                                                                   P.Value = "p.val", adj.P.Val = "p.adj")) %>% unite(temp, comparison, 
                                                                                                                                                                                                                                      variable) %>% spread(temp, value)


rowData(data_imp_BA17) <- merge(rowData(data_imp_BA17, use.names = FALSE), table2, 
                            by.x = "name", by.y = "rowname", all.x = TRUE, sort = FALSE)


cols <- colnames(rowData(data_imp_BA17, use.names = FALSE) %>% as.data.frame())
cols
cols_p <- grep("_p.adj", colnames(rowData(data_imp_BA17, use.names = FALSE) %>% as.data.frame()))


dep <- add_rejections(data_imp_BA17, alpha = 0.05, lfc = log2(2))
sum(dep@elementMetadata$significant)


data_results <- get_results(dep)
data_results %>% filter(significant) %>% nrow()
colnames(data_results)

#remove the effect of degradation
dep_corrected <- dep
design_corrected <- model.matrix(~condition)
colnames(design_corrected)[2:ncol(design_corrected)] <- c("GW22_23","GW28_41", "Year00_01","Year04_08", "Year18_22")

batch <- factor(design[,7])

#we can use contr.sum(2) to do effect coding, thus after batch effect removal, it will be the average level of
#non degraded and degraded samples. However, here we want to remove the effect of degraded samples and only keep
#non degraded samples levels so the following 4 lines are omiteted from removeBatchEffect function
#The contrast function, contr.sum(), gives orthogonal contrasts where you compare every level to the overall mean.
#contrasts(batch1) <- contr.sum(levels(batch1))
#contrasts(batch2) <- contr.sum(levels(batch2))
#contrasts(batch3) <- contr.sum(levels(batch3))
#contrasts(batch4) <- contr.sum(levels(batch4))

covariates <- model.matrix(~batch)
covariates <- covariates[,-1]

tmp <- removeBatchEffect(assay(dep), covariates = covariates, design =design_corrected)
assay(dep_corrected) <- tmp


#plot PCA based on condition
var <- apply(assay(dep_corrected), 1, sd)
df <- assay(dep_corrected)[order(var, decreasing = TRUE)[seq_len(1000)],]
# Calculate PCA
pca <- prcomp(t(df), scale = FALSE)
pca_df <- pca$x %>%
  data.frame() %>%
  rownames_to_column() %>%
  left_join(., data.frame(colData(dep_corrected)), by = c("rowname" = "ID"))

# Calculate the percentage of variance explained
percent <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df$`Days (log2)` <- dep_corrected$Log2AgeDays
pca_df$PC2 <- -(pca_df$PC2)

# plot PCA based on condition
ggplot(pca_df, aes(get(paste0("PC", 1)), get(paste0("PC", 2)))) +
  labs(title = paste0("PCA plot - top ", 1000, " variable proteins"),
       x = paste0("PC", 1, ": ", percent[1], "%"),
       y = paste0("PC", 2, ": ", percent[2], "%")) +
  theme_DEP1() +
  geom_point(aes(col = condition),
             size = 2) +
  labs(col = "condition") +
  scale_color_manual(values= RColorBrewer::brewer.pal(6, "Set2"))

ggsave(filename = "Results/pca.pdf", width = 6, height = 3, units = "in")

#plot heatmap
#the heatmap only uses top 500 significant genes
dep2 <- add_rejections(data_imp_BA17, alpha =0.000001, lfc = log2(2))
assay(dep2) <- assay(dep_corrected)
sum(dep2@elementMetadata$significant)
#mark genes with new color
source(file = "functions/get_annotation_brewer_pal.R")
source(file = "functions/plot_heatmap_brewer_pal_seriation.R")
ht <- plot_heatmap_brewer_pal_seriation(dep2, type = "centered", kmeans = F, col_limit = 5, show_row_names = FALSE, show_column_names = FALSE,
                                        indicate = c("condition"), clustering_distance = "spearman")
pdf("Results/heatmap.pdf", width = 8, height = 6)
draw(ht)
dev.off()

saveRDS(dep, file = "Results/dep_BA17.rds")
saveRDS(dep_corrected, file = "Results/dep_corrected_BA17.rds")
