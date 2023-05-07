library("tidyverse")
library("DOSE")

result <- read.csv(file = "Results/SynGO_ORA/syngo_ora.csv", header = TRUE)
chosen_pathways <- read.csv(file = "Results/SynGO_ORA/chosen_pathways.csv",header = TRUE)

res_filt <- result %>%
  filter(term %in% chosen_pathways$term)

#reorder cluster based on peak time manually
cluster_levels <- unique(chosen_pathways$cluster)
res_filt$cluster <- factor(res_filt$cluster,levels= cluster_levels)
res_filt <- res_filt[order(res_filt$cluster),]

#reorder description based on p-value and first five in each module
pahtway_levels <- res_filt$term
pahtway_levels <- unique(pahtway_levels)
res_filt$term <- factor(res_filt$term, levels=rev(pahtway_levels))

res_filt <- res_filt %>% mutate("-log10(p.adj)" = -log10(p.adjust))

#draw dotplot
p <- ggplot(data = res_filt, aes(x = cluster, y = term, color = `-log10(p.adj)`)) +
  geom_point(size = 4) + 
  scale_colour_viridis_c(guide=guide_colorbar(reverse=TRUE), limits = c(1.3,5), oob=scales::squish)
p <- p + xlab("") + ylab("") + ggtitle("") +
  theme_dose(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave(filename = "Results/SynGO_ORA/Module_syngo_enrichment.pdf", width = 10.5, height = 7, units = "in")
