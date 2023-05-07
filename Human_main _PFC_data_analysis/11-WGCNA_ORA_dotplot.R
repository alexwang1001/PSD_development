library("tidyverse")
library("DOSE")

result <- read.csv(file = "Results/WGCNA_ORA/res_MSigDB_c2cp_ORA.csv", header = TRUE)
#choose top 8 pathways for each developmental stage to plot
chosen_pathways <- read.csv(file = "Results/WGCNA_ORA/chosen_pathways_8.csv",header = TRUE)

gsize <- as.numeric(sub("/\\d+$", "", as.character(result$GeneRatio)))
gcsize <- as.numeric(sub("^\\d+/", "", as.character(result$GeneRatio)))
result$GeneRatio = gsize/gcsize


res_filt <- result %>%
  filter(Description %in% chosen_pathways$Description)

#reorder Cluster based on peak time manually
Cluster_levels <- unique(res_filt$Cluster)
Cluster_levels <- Cluster_levels[c(4,2,1,3)]
res_filt$Cluster <- factor(res_filt$Cluster,levels= Cluster_levels)
res_filt <- res_filt[order(res_filt$Cluster),]

#reorder description based on p-value and first five in each module
pahtway_levels <- chosen_pathways$Description
pahtway_levels <- unique(pahtway_levels)
res_filt$Description <- factor(res_filt$Description, levels=rev(pahtway_levels))

res_filt <- res_filt %>% mutate("-log10(p.adj)" = -log10(p.adjust))
#draw dotplot
p <- ggplot(data = res_filt, aes(x = Cluster, y = Description, size = GeneRatio, color = `-log10(p.adj)`)) +
  geom_point() + 
  scale_colour_viridis_c(guide=guide_colorbar(reverse=TRUE), limits = c(1.3,8), oob=scales::squish)
p <- p + xlab("") + ylab("") + ggtitle("") +
  theme_dose(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave(filename = "Results/WGCNA_ORA/Module_pathway_enrichment.pdf", width = 9, height = 7, units = "in")

