library("tidyverse")
library("SummarizedExperiment")

species_comparison <- readRDS("Results/species_comparison_humanized_age.rds")

rowdata <- as.data.frame(rowData(species_comparison))

brown <- species_comparison[rowdata$Human_Module == "brown", ]
brown_expr <- data.frame(assay(brown))
brown_median <- apply(brown_expr, 2, median)
brown$brown <- brown_median
df_brown <- data.frame(colData(brown))

blue <- species_comparison[rowdata$Human_Module == "blue", ]
blue_expr <- data.frame(assay(blue))
blue_median <- apply(blue_expr, 2, median)
blue$blue <- blue_median
df_blue <- data.frame(colData(blue))

yellow <- species_comparison[rowdata$Human_Module == "yellow", ]
yellow_expr <- data.frame(assay(yellow))
yellow_median <- apply(yellow_expr, 2, median)
yellow$yellow <- yellow_median
df_yellow <- data.frame(colData(yellow))

turquoise <- species_comparison[rowdata$Human_Module == "turquoise", ]
turquoise_expr <- data.frame(assay(turquoise))
turquoise_median <- apply(turquoise_expr, 2, median)
turquoise$turquoise <- turquoise_median
df_turquoise <- data.frame(colData(turquoise))

df_all <- df_brown %>% left_join(df_blue) %>% left_join(df_turquoise) %>% left_join(df_yellow)
write.csv(df_all, file = "Results/species_comparison/module_median_abundance_across_species.csv")


#plot module median levels in different species
df_all_long <- df_all %>% pivot_longer(cols = 9:12, names_to = "module", values_to = "median")

df_all_long$module <- factor(df_all_long$module, levels = c("brown", "blue", "turquoise", "yellow"))

ggplot(data = df_all_long, mapping = aes(x = humanized_log2_age_days, y = median)) +
  geom_smooth(mapping = aes(color = species) ,se = F, span = 1, show.legend = T) +
  geom_point(size = 0.8, alpha = 0.8, show.legend = F) +
  facet_wrap(facets = species ~ module)

ggplot(data = df_all_long, mapping = aes(x = humanized_log2_age_days, y = median, color = species)) +
  geom_smooth(size = 1, alpha = 0.8, se = F, span = 1, show.legend = T) +
  geom_point(size = 1, alpha = 0.8, show.legend = F) +
  facet_wrap(facets = ~ module) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "gray80", linetype = "dashed"),
        panel.border = element_rect(fill = NA),
        panel.background = element_rect(fill = "gray98", color = "black")
  ) +
  theme(axis.text=element_text(size=12, color = "black"),
        axis.text.x = element_text(face = c("plain", "bold", "plain", "plain", "plain", "plain"),hjust = c(0.5,0.5,0.5,0.5,0.5,0.5)),
        axis.title = element_text(size=14)) +
  theme(strip.background = element_rect(color = "beige", fill="beige"),
        strip.text = element_text(size = 16)) +
  theme(legend.key.size = unit(0.4, 'inch'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "top") +
  geom_vline(xintercept = 8.0552824355) +
  labs(x = "Humanized post-conceptional age (log-transformed)",
       y = "Median abundance") +
  scale_x_continuous(breaks = c(6.8073549221,8.0552824355,9.30149619498,10.7532167492,12.885315061), labels = c("GW18", "Birth", "Year01", "Year04", "Year20"))

ggsave("Results/species_comparison/module_median_levels_across_species.pdf", width = 6, height = 6)
ggsave("Results/species_comparison/module_median_levels_across_species.png", width = 6, height = 6)

#plot RhoGEF proteins abundance in different species
RhoGEFs <- species_comparison[c("ARHGEF7", "ARHGEF2", "NGEF", "RASGRF2", "ABR"),]

df <- data.frame(assay(RhoGEFs)) %>%
  rownames_to_column() %>%
  gather(label, val, -rowname) %>%
  left_join(., data.frame(colData(RhoGEFs)), by = "label")
colnames(df)[1] <- "Protein"
#add PREX1 which is only identified in humans
load(file = "../Human PSD filtered by non post-mortem tissue/Results/dep_corrected.RData")
PREX1 <- dep_corrected["PREX1",]
df_PREX1 <- data.frame(assay(PREX1)) %>%
  rownames_to_column() %>%
  gather(label, val, -rowname) %>%
  left_join(., data.frame(colData(RhoGEFs)), by = "label")
colnames(df_PREX1)[1] <- "Protein"

df_final <- rbind(df, df_PREX1)

p <- ggplot(data = df_final, mapping = aes(x = humanized_log2_age_days, y = val, color = species, fill = species)) +
  geom_smooth(se = F, show.legend = F) +
  geom_point(size = 0.8, alpha = 1, show.legend = T) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "gray80", linetype = "dashed"),
        panel.border = element_rect(fill = NA),
        panel.background = element_rect(fill = "gray98", color = "black")
  ) +
  theme(axis.text=element_text(size=12, color = "black"),
        axis.text.x = element_text(face = c("plain", "bold", "plain", "plain", "plain", "plain"),hjust = c(0.5,0.5,0.5,0.5,0.5,0.5)),
        axis.title = element_text(size=14) ) +
  theme(strip.background = element_rect(color = "beige", fill="beige"),
        strip.text = element_text(size = 14)) +
  geom_vline(xintercept = 8.0552824355) +
  labs(x = "Humanized post-conceptional age (log-transformed)",
       y = expression("Abundance ("~log[2]~")")) +
  scale_x_continuous(breaks = c(6.8073549221,8.0552824355,9.30149619498,10.7532167492,12.885315061), labels = c("GW18", "Birth", "Year01", "Year04", "Year20")) +
  scale_y_continuous(limits = c(-2,12)) +
  guides(color = guide_legend(override.aes = list(size = 4))) + theme(legend.position="top") +
  facet_wrap(~Protein, ncol = 2)
p

ggsave("Results/species_comparison/RhoGEF_abundance_across_species.pdf", width = 6, height = 7)
