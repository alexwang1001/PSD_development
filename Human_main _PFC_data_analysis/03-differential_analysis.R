library("DEP")
library("SummarizedExperiment")
library("ComplexHeatmap")
library("limma")
library("vsn")
library("tidyverse")
library("scales")
library("patchwork")

load(file = "./Results/data_imp.RData")

data_ctrl <- data_imp
colData(data_ctrl)<- droplevels(colData(data_ctrl))


#using limma
condition <- as.factor(data_ctrl$condition)
degradation <- as.factor(data_ctrl$Degradation)
design <- model.matrix(terms(~ condition + degradation, keep.order = TRUE))
colnames(design)
colnames(design)[2:7] <- c("GW22_23","GW28_41", "Year00_01", "Year04_08", "Year18_22", "degradation")
colnames(design)
fit1 <- lmFit(assay(data_ctrl), design=design)
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


rowData(data_ctrl) <- merge(rowData(data_ctrl, use.names = FALSE), table, 
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


rowData(data_ctrl) <- merge(rowData(data_ctrl, use.names = FALSE), table2, 
                            by.x = "name", by.y = "rowname", all.x = TRUE, sort = FALSE)


cols <- colnames(rowData(data_ctrl, use.names = FALSE) %>% as.data.frame())
cols
cols_p <- grep("_p.adj", colnames(rowData(data_ctrl, use.names = FALSE) %>% as.data.frame()))


dep <- add_rejections(data_ctrl, alpha = 0.05, lfc = log2(2))
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
#non-degraded and partially-degraded samples. However, here we want to remove the effect of degradation and only keep
#non-degraded levels so the following 4 lines were omitted before applying the removeBatchEffect function
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
#plot PCA based on age
p1 <- ggplot(pca_df, aes(get(paste0("PC", 1)), get(paste0("PC", 2)))) +
  labs(title = paste0("PCA plot - top ", 1000, " variable proteins"),
       x = paste0("PC", 1, ": ", percent[1], "%"),
       y = paste0("PC", 2, ": ", percent[2], "%")) +
  coord_fixed() +
  theme_DEP1() +
  geom_point(aes(col = `Days (log2)`),
             size = 2) +
  labs(col = "days (log2)") +
  scale_color_gradientn(colours = rainbow(5)) +
  guides(#reverse color order (low value on top)
    color = guide_colorbar(reverse = T))

#plot PCA based on sex
p2 <- ggplot(pca_df, aes(get(paste0("PC", 1)), get(paste0("PC", 2)))) +
  labs(title = paste0("PCA plot - top ", 1000, " variable proteins"),
       x = paste0("PC", 1, ": ", percent[1], "%"),
       y = paste0("PC", 2, ": ", percent[2], "%")) +
  coord_fixed() +
  theme_DEP1() +
  geom_point(aes(col = Sex),
             size = 2) +
  labs(col = "sex")


#plot PCA based on degradation
pca_df$Quality <- ifelse(pca_df$Degradation == 0, "good", "fair")
p3 <- ggplot(pca_df, aes(get(paste0("PC", 1)), get(paste0("PC", 2)))) +
  labs(title = paste0("PCA plot - top ", 1000, " variable proteins"),
       x = paste0("PC", 1, ": ", percent[1], "%"),
       y = paste0("PC", 2, ": ", percent[2], "%")) +
  coord_fixed() +
  theme_DEP1() +
  geom_point(aes(col = Quality),
             size = 2) +
  labs(col = "quality")

#plot PCA based on process date
p4 <- ggplot(pca_df, aes(get(paste0("PC", 1)), get(paste0("PC", 2)))) +
  labs(title = paste0("PCA plot - top ", 1000, " variable proteins"),
       x = paste0("PC", 1, ": ", percent[1], "%"),
       y = paste0("PC", 2, ": ", percent[2], "%")) +
  coord_fixed() +
  theme_DEP1() +
  geom_point(aes(col = ProcessTime),
             size = 2) +
  labs(col = "processing date") +
  guides(col=guide_legend(ncol=2))

p1+p2+p3+p4
ggsave(filename = "Results/pca_other.pdf", width = 12, height = 6, units = "in")
#PCA finished

plot_cor(dep_corrected, significant = F, lower = 0.3, upper = 1, pal = "OrRd")




#the heatmap only uses top 500 significant genes
dep2 <- add_rejections(data_ctrl, alpha = 0.00000088, lfc = log2(2))
assay(dep2) <- assay(dep_corrected)
sum(dep2@elementMetadata$significant)


#mark genes with new col annotation color
source(file = "functions/plot_heatmap_brewer_pal_rev_mark.R")
source(file = "functions/get_annotation_brewer_pal_rev.R")
ht <- plot_heatmap_brewer_rev_mark(dep2, type = "centered", kmeans = F, col_limit = 5, show_row_names = FALSE, 
                                     show_column_names = FALSE, indicate = c("condition"), use_raster = FALSE, clustering_distance = "spearman")
pdf("Results/heatmap_with_mark.pdf", width = 8, height = 6)
draw(ht)
dev.off()

#example proteins
source(file = "functions/plot_single2.R")
p1 <- plot_single2(dep_corrected, proteins = c("GRIN2A", "GRIN2B")) + guides(color = guide_legend(override.aes = list(size = 4))) + theme(legend.position="top")
# ggsave(filename = "Results/protein_examples/GRIN2A-B.pdf", width = 6, height = 5, units = "in")

p2 <- plot_single2(dep_corrected, proteins = c("DLG3", "DLG4")) + guides(color = guide_legend(override.aes = list(size = 4))) + theme(legend.position="top")
ggsave(filename = "Results/protein_examples/DLG3-4.pdf", width = 6, height = 5, units = "in")


#Ribosomes
p3 <- plot_single2(dep_corrected, proteins = c("RPS19", "RPS2", "RPS3", "RPS4X", "RPS7", "RPS9")) + guides(color = guide_legend(override.aes = list(size = 4))) + theme(legend.position="top")
# ggsave(filename = "Results/protein_examples/RPS.pdf", width = 6, height = 5, units = "in")

#CCT chaperones
p4 <- plot_single2(dep_corrected, proteins = c("CCT8", "CCT7", "CCT3", "CCT4", "TCP1", "CCT5", "CCT6A", "CCT2")) + guides(color = guide_legend(override.aes = list(size = 4))) + theme(legend.position="top")
# ggsave(filename = "Results/protein_examples/CCT.pdf", width = 6, height = 5, units = "in")

p5 <- plot_single2(dep_corrected, proteins = c("SHANK1", "SHANK2", "SHANK3")) + guides(color = guide_legend(override.aes = list(size = 4))) + theme(legend.position="top")
# ggsave(filename = "Results/protein_examples/SHANK.pdf", width = 6, height = 5, units = "in")
p1+p2+p3+p4+p5 +  plot_layout(ncol = 2)
ggsave(filename = "Results/protein_examples/FigS2.pdf", width = 10, height = 12, units = "in")

#join 4 figures together using patchwork for Fig.1F
library(patchwork)
p1 <- plot_single2(dep_corrected, proteins = c("RPL3", "RPL12", "RPL10", "RPL13A", "RPL14", "RPL24")) + guides(color = guide_legend(override.aes = list(size = 4))) + theme(legend.position="top")
p2 <- plot_single2(dep_corrected, proteins = c("CTNNA1", "CTNNA2", "CTNNB1", "CTNND1", "CTNND2")) + guides(color = guide_legend(override.aes = list(size = 4), nrow=2)) + theme(legend.position="top")
p3 <- plot_single2(dep_corrected, proteins = c("TRIO", "ABR","ARHGEF2", "PREX1")) + guides(col = guide_legend(nrow=2)) + guides(color = guide_legend(override.aes = list(size = 4))) + theme(legend.position="top")
p4 <- plot_single2(dep_corrected, proteins = c("GRIA1", "GRIA2", "GRIA3", "GRIA4")) + guides(color = guide_legend(override.aes = list(size = 4))) + theme(legend.position="top")
p1+p2+p3+p4
ggsave(filename = "Results/protein_examples/all_GSEA.pdf", width = 8, height = 8, units = "in")



save(dep, file = "./Results/dep.RData")
save(dep_corrected, file = "./Results/dep_corrected.RData")
