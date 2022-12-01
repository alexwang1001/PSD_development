library("tidyverse")
library("DOSE")

chosen_pathways <- read.csv(file = "Results/age_specificity/chosen_pathways.csv",header = TRUE)
specificity_GESA_res <- read.csv(file = "Results/age_specificity/specificity_GESA_res.csv",header = TRUE)

specificity_GESA_res_filt <- specificity_GESA_res %>% filter(NES > 0, p.adjust < 0.05) %>% filter(ID %in% chosen_pathways$ID)

specificity_GESA_res_filt <- specificity_GESA_res_filt %>% mutate("-log10(p.adj)" = -log10(p.adjust))
specificity_GESA_res_filt$Description <- factor(specificity_GESA_res_filt$Description, levels = rev(unique(chosen_pathways$ID)))

p <- ggplot(data = specificity_GESA_res_filt, aes(x = Cluster, y = Description, size = NES, color = `-log10(p.adj)`)) +
  geom_point() + 
  scale_colour_viridis_c(guide=guide_colorbar(reverse=TRUE))
p <- p + xlab("") + ylab("") + ggtitle("") +
  theme_dose(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave(filename = "Results/age_specificity/age_specificity_dotplot.pdf", width = 9, height = 6, units = "in")
