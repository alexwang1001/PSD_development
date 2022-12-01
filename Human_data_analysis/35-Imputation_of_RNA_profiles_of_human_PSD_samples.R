library("tidyverse")
library("SummarizedExperiment")
library("mgcv")
library("visreg")
library("WGCNA")
library("FSA")

#import RNA-seq and gene module information
RNA <- readRDS("Results/Human_NCX_subset_transcriptomics_PSD_logTPM.rds")
#filter out genes with 0 varaince or too many missing values


rowdata <- data.frame(rowData(RNA))
#import PSD protein data
load(file = "Results/dep_corrected.RData")
#impute RNA profiles at PSD sample ages
results <- data.frame("log2_age_days" = dep_corrected$Log2AgeDays)
for (gene in rowdata$SYMBOL) {
  RNA_subset <- RNA[gene,]
  data <- NULL
  data$log2_age_days <- RNA_subset$log2_age_days
  data$expression <- as.numeric(assay(RNA_subset))
  data <- data.frame(data)
  # Build the model
  model <- gam(expression ~ s(log2_age_days), data = data,  method = "REML")
  # Make predictions
  predictions <- model %>% predict(results)
  results[gene] <- predictions
}
#build a new summarizedexperiment
assay <- t(results[,-1])
coldata <- data.frame(colData(dep_corrected))
rowdata <- rowdata
colnames(assay) <- coldata$label
RNA_imputed <- SummarizedExperiment(assays=assay, rowData =rowdata, colData=coldata)

saveRDS(RNA_imputed, "Results/Human_NCX_imputed_RNA_PSD_logTPM.rds")

assay$SYMBOL <- rownames(assay)
assay2 <- rowdata %>% left_join(assay)
assay3 <- assay2[,c(colnames(assay2)[1:4],sort(colnames(assay2)[5:58]))]
write.csv(assay3, "Results/RNA-seq/Imputed_RNA_profiles.csv", row.names = F)

#calculate correlation between RNA and protein in each module
protein <- dep_corrected[rownames(RNA_imputed),]
protein_expr <- data.frame(t(assay(protein)))
RNA_expr <- data.frame(t(assay(RNA_imputed)))
res <- data.frame("correlation" = mapply(cor, protein_expr, RNA_expr,  method = "spearman"))
res$SYMBOL <- rowdata$SYMBOL
res <- left_join(res, rowdata)

write.csv(res, "Results/RNA-seq/RNA_protein_correlation.csv", row.names = F)

res$geneModules <- factor(res$geneModules, levels = c("brown", "blue", "turquoise", "yellow", "grey"))
res %>%
  ggplot( aes(x=geneModules, y=correlation, color = geneModules)) +
  geom_boxplot(notch = TRUE, outlier.alpha = 0.25) +
  scale_color_manual(values = c("brown", "blue", "turquoise", "gold", "grey")) +
  theme_classic() +
  labs(x = "Category", y = "Protein-RNA Correlation", color = "Category") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position="none")
ggsave("Results/RNA-seq/Protein-RNA_correlation.pdf", width = 2.5,height = 3)

kruskal.test(correlation ~ geneModules, data = res)
dunnTest(correlation ~ geneModules, data=res, method = "bh")
