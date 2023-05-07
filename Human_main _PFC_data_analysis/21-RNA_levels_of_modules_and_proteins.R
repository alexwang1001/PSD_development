library(tidyverse)
library(SummarizedExperiment)
library(GenomicRanges)
library(patchwork)


se <- readRDS("Results/Human_NCX_subset_transcriptomics_PSD_logTPM.rds")

#calculate median values of each module and standardize it
get_scaled_median_for_modules <- function(se) {
  rowdata<- rowData(se)
  coldata <- colData(se)
  modules <- unique(rowdata[["geneModules"]])
  for (module in modules) {
    module_se <- se[rowdata$geneModules == module,]
    module_se_expr <- data.frame(assay(module_se))
    module_median <- apply(module_se_expr, 2, median)
    module_median_scaled <- scale(module_median)
    colData(se)[module] <- module_median_scaled
  }
  data.frame(colData(se))
}

module_median <- get_scaled_median_for_modules(se)
write.csv(module_median, "Results/RNA-seq/scaled_median.csv", row.names = F)

module_median$grey <- NULL
module_median_long <- gather(module_median, key = "module", value = "median", 9:12)


ggplot(data = module_median_long, mapping = aes(x = log2_age_days, y = median, color = module, fill = module)) +
  geom_smooth(se = F, show.legend = T) +
  geom_point(size = 0.8, alpha = 0.8, show.legend = F) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "gray80", linetype = "dashed"),
        panel.border = element_rect(fill = NA),
        panel.background = element_rect(fill = "gray98", color = "black")
  ) +
  theme(axis.text=element_text(size=12, color = "black"),
        axis.text.x = element_text(face = c("plain", "plain", "bold", "plain", "plain", "plain"),hjust = c(0.8,0.2,0.5,0.5,0.5,0.5)),
        axis.title = element_text(size=14) ) +
  geom_vline(xintercept = 8.0552824355) +
  labs(x = "Post-conceptional age (log-transformed)",
       y = "Standardized expression") +
  scale_color_manual(values = c("blue","brown","turquoise","gold")) +
  scale_fill_manual(values = c("blue","brown","turquoise","gold")) +
  scale_x_continuous(breaks = c(6.8073549221,7.1996723448,8.0552824355,9.30149619498,10.7532167492,11.6375305515,12.885315061), labels = c("GW18", "GW23", "Birth", "Year01", "Year04", "Year08", "Year20"))

ggsave("Results/RNA-seq/standardized_median.pdf", width = 6, height = 4)


#several examples showing the difference between RNA and Protein
source("functions/plot_single2_normalized.R")
source("functions/plot_single2_RNA_normalized.R")
load("Results/dep_corrected.RData")
p1 <- plot_single2_normalized(dep_corrected, proteins = c("NGEF","RASGRF2", "DLG4", "DLG1"), scale = T)+ guides(color = guide_legend(override.aes = list(size = 4))) + theme(legend.position="top")
ggsave(filename = "Results/RNA-seq/turquoise_yellow_protein_examples.pdf", width = 5, height = 4, units = "in")
p2 <- plot_single2_RNA_normalized(se, genes = c("NGEF","RASGRF2", "DLG4", "DLG1"), scale = T)+ guides(color = guide_legend(override.aes = list(size = 4))) + theme(legend.position="top")
ggsave(filename = "Results/RNA-seq/turquoise_yellow_RNA_examples.pdf", width = 5, height = 4, units = "in")

p1+p2
ggsave(filename = "Results/RNA-seq/turquoise_yellow_RNA_examples_together.pdf", width = 9, height = 4, units = "in")

