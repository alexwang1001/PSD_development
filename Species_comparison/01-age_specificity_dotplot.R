library("tidyverse")
library("DOSE")
library("ComplexHeatmap")

chosen_pathways_human <- read.csv(file = "Results/age_specificity_GESA_comparison/chosen_pathways.csv",header = TRUE)
specificity_GESA_res_human <- read.csv(file = "Results/age_specificity_GESA_comparison/specificity_GESA_res_human.csv",header = TRUE, row.names = 1)
specificity_GESA_res_human$Species <- "Human"

specificity_GESA_res_macaque <- read.csv(file = "Results/age_specificity_GESA_comparison/specificity_GESA_res_macaque.csv",header = TRUE, row.names = 1)
specificity_GESA_res_macaque$Species <- "Macaque"

specificity_GESA_res_mouse <- read.csv(file = "Results/age_specificity_GESA_comparison/specificity_GESA_res_mouse.csv",header = TRUE, row.names = 1)
specificity_GESA_res_mouse$Species <- "Mouse"

specificity_GESA_res <- rbind(specificity_GESA_res_human, specificity_GESA_res_macaque, specificity_GESA_res_mouse)

specificity_GESA_res_filt <- specificity_GESA_res %>% filter(NES > 0, p.adjust < 0.05) %>% filter(ID %in% chosen_pathways_human$ID)


specificity_GESA_res_filt <- specificity_GESA_res_filt %>% mutate("-log10(p.adj)" = -log10(p.adjust))

specificity_GESA_res_filt$Description <- factor(specificity_GESA_res_filt$Description, levels = rev(unique(chosen_pathways_human$ID)))
specificity_GESA_res_filt$Cluster <- factor(specificity_GESA_res_filt$Cluster, levels = c("GW18_19", "GW22_23", "GW28_41", "Year00_01", "Year04_08", "Year18_22",
                                                                                          "E75_85", "E95_110", "Year0_1", "Year1_3", "Year8_10",
                                                                                          "P0", "P9", "P13", "P18", "P36"))


p <- ggplot(data = specificity_GESA_res_filt, aes(x = Cluster, y = Description, size = NES, color = `-log10(p.adj)`)) +
  geom_point() + 
  scale_colour_viridis_c(guide=guide_colorbar(reverse=TRUE))
p <- p + xlab("") + ylab("") + ggtitle("") +
  theme_dose(12) +
  theme(strip.background = element_rect(color = "beige", fill="beige"),
        strip.text = element_text(size = 12)) +
  guides(x =  guide_axis(angle = 45)) +
  facet_wrap(~Species, scales="free_x")
p

ggsave(filename = "Results/age_specificity_GESA_comparison/dotplot.pdf", width = 12, height = 6, units = "in")


