library("DEP")
library("SummarizedExperiment")
library("limma")
library("vsn")
library("tidyverse")
library("scales")

load(file = "./Results/data_imp.RData")
data_imp
data_imp$condition <- factor(data_imp$condition, levels = c("E75_85", "E95_110", "Year0_1", "Year1_3", "Year8_10"))

#using limma
condition <- data_imp$condition
design <- model.matrix(terms(~ condition))
colnames(design)
colnames(design)[2:5] <- c("E95_110", "Year0_1", "Year1_3", "Year8_10")
colnames(design)
fit1 <- lmFit(assay(data_imp), design=design)
fit1 <- eBayes(fit1)
head(coef(fit1))
colnames(coef(fit1))
#save comparison data
cntrst <- colnames(design)[2:5]
retrieve_fun <- function(comp, fit) {
  res <- topTable(fit, sort.by = "p", coef = comp, number = Inf, 
                  confint = TRUE)
  res <- res[!is.na(res$P.Value), ]
  res$comparison <- rep(comp , dim(res)[1])
  res <- rownames_to_column(res)
  return(res)
}

limma_res <- map_df(cntrst, retrieve_fun, fit = fit1)
table <- limma_res %>% select(rowname, logFC, CI.L, CI.R, t, 
                              P.Value, adj.P.Val, comparison) %>% mutate(comparison = paste0(comparison,
                                                                                             "_vs_E75_85")) %>% gather(variable, value, -c(rowname, 
                                                                                                                                       comparison)) %>% mutate(variable = recode(variable, logFC = "diff", 
                                                                                                                                                                                 P.Value = "p.val", adj.P.Val = "p.adj")) %>% unite(temp, comparison, 
                                                                                                                                                                                                                                    variable) %>% spread(temp, value)

rowData(data_imp) <- merge(rowData(data_imp, use.names = FALSE), table, 
                           by.x = "name", by.y = "rowname", all.x = TRUE, sort = FALSE)

#make all pairwise comparison
cntrst2 <- apply(utils::combn(colnames(design)[5:2] , 2), 2, paste, 
                 collapse = " - ")
made_contrasts <- makeContrasts(contrasts = cntrst2, levels = design)
contrast_fit <- contrasts.fit(fit1, made_contrasts)
eB_fit <- eBayes(contrast_fit)
limma_res2 <- map_df(cntrst2, retrieve_fun, fit = eB_fit)
table2 <- limma_res2 %>% select(rowname, logFC, CI.L, CI.R, t, 
                                P.Value, adj.P.Val, comparison) %>% mutate(comparison = gsub(" - ", "_vs_",
                                                                                             comparison)) %>% gather(variable, value, -c(rowname, 
                                                                                                                                         comparison)) %>% mutate(variable = recode(variable, logFC = "diff", 
                                                                                                                                                                                   P.Value = "p.val", adj.P.Val = "p.adj")) %>% unite(temp, comparison, 
                                                                                                                                                                                                                                      variable) %>% spread(temp, value)


rowData(data_imp) <- merge(rowData(data_imp, use.names = FALSE), table2, 
                           by.x = "name", by.y = "rowname", all.x = TRUE, sort = FALSE)


cols <- colnames(rowData(data_imp, use.names = FALSE) %>% as.data.frame())
cols
cols_p <- grep("_p.adj", colnames(rowData(data_imp, use.names = FALSE) %>% as.data.frame()))


dep <- add_rejections(data_imp, alpha = 0.01, lfc = log2(2))

sum(dep@elementMetadata$significant)

plot1 <- plot_pca(dep,indicate = "condition", n = 1000, point_size = 3)

#plot PCA with axis flipped
var <- apply(assay(dep), 1, sd)
df <- assay(dep)[order(var, decreasing = TRUE)[seq_len(1000)],]
# Calculate PCA
pca <- prcomp(t(df), scale = FALSE)
pca_df <- pca$x %>%
  data.frame() %>%
  rownames_to_column() %>%
  left_join(., data.frame(colData(dep)), by = c("rowname" = "ID"))
pca_df$PC1 = -(pca_df$PC1)
pca_df$PC2 = -(pca_df$PC2)
# Calculate the percentage of variance explained
percent <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pca_df$condition <- factor(pca_df$condition, levels = c("E75_85", "E95_110", "Year0_1", "Year1_3", "Year8_10"))
plot2 <- ggplot(pca_df, aes(get(paste0("PC", 1)), get(paste0("PC", 2)))) +
  labs(title = paste0("PCA plot - top ", 1000, " variable proteins"),
       x = paste0("PC", 1, ": ", percent[1], "%"),
       y = paste0("PC", 2, ": ", percent[2], "%")) +
  theme_DEP1() +
  geom_point(aes(col = fct_relevel(condition, "E75_85", "E95_110", "Year0_1", "Year1_3", "Year8_10")),
             size = 2) +
  labs(col = "condition") +
  coord_fixed(ratio = 1.2, expand = TRUE) +
  scale_color_manual(values= RColorBrewer::brewer.pal(5, "Set2"))

#PCA finished
ggsave(filename = "Results/pca.pdf", width = 6, height = 3, units = "in")


#the heatmap using significant genes
source(file = "functions/get_annotation_brewer_pal.R")
source(file = "functions/plot_heatmap_brewer_pal_seriation.R")
set.seed(0)
pdf("Results/heatmap.pdf", width = 8, height = 6)
plot_heatmap_brewer_pal_seriation(dep, type = "centered", kmeans = F, col_limit = 3, show_row_names = FALSE, show_column_names = FALSE,
                               indicate = c("condition"), clustering_distance = "spearman")
dev.off()

#save data
save(dep, file = "./Results/dep.RData")
