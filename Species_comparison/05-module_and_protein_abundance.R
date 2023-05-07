library(tidyverse)
library(SummarizedExperiment)

species_comparison <- readRDS("Results/species_comparison_humanized_age.rds")

rowdata <- as.data.frame(rowData(species_comparison))

brown <- species_comparison[rowdata$Human_Module == "brown", ]
brown_expr <- data.frame(assay(brown))
brown_expr_scaled <- t(scale(t(brown_expr)))
brown_mean <- colMeans(brown_expr_scaled)
brown_sem <- apply(brown_expr_scaled, 2, function(x) sd(x)/sqrt(length(x)))
brown$brown_mean <- brown_mean
brown$brown_sem <- brown_sem
df_brown <- data.frame(colData(brown))

blue <- species_comparison[rowdata$Human_Module == "blue", ]
blue_expr <- data.frame(assay(blue))
blue_expr_scaled <- t(scale(t(blue_expr)))
blue_mean <- colMeans(blue_expr_scaled)
blue_sem <- apply(blue_expr_scaled, 2, function(x) sd(x)/sqrt(length(x)))
blue$blue_mean <- blue_mean
blue$blue_sem <- blue_sem
df_blue <- data.frame(colData(blue))

yellow <- species_comparison[rowdata$Human_Module == "yellow", ]
yellow_expr <- data.frame(assay(yellow))
yellow_expr_scaled <- t(scale(t(yellow_expr)))
yellow_mean <- colMeans(yellow_expr_scaled)
yellow_sem <- apply(yellow_expr_scaled, 2, function(x) sd(x)/sqrt(length(x)))
yellow$yellow_mean <- yellow_mean
yellow$yellow_sem <- yellow_sem
df_yellow <- data.frame(colData(yellow))




turquoise <- species_comparison[rowdata$Human_Module == "turquoise", ]
turquoise_expr <- data.frame(assay(turquoise))
turquoise_expr_scaled <- t(scale(t(turquoise_expr)))
turquoise_mean <- colMeans(turquoise_expr_scaled)
turquoise_sem <- apply(turquoise_expr_scaled, 2, function(x) sd(x)/sqrt(length(x)))
turquoise$turquoise_mean <- turquoise_mean
turquoise$turquoise_sem <- turquoise_sem
df_turquoise <- data.frame(colData(turquoise))



df_all <- df_brown %>% left_join(df_blue) %>% left_join(df_turquoise) %>% left_join(df_yellow)
write.csv(df_all, file = "Results/species_comparison/module_mean_and_sem_across_species.csv")


#plot module median levels in different species
df_all_long <- df_all %>% pivot_longer(cols = 9:16, names_to = c("module", ".value"), names_pattern = "(.*)(_.*)$")
colnames(df_all_long)[c(10,11)] <- c("mean", "sem")

df_all_long$module <- factor(df_all_long$module, levels = c("brown", "blue", "turquoise", "yellow"))

ggplot(data = df_all_long, mapping = aes(x = humanized_log2_age_days, y = mean)) +
  geom_smooth(mapping = aes(color = species) ,se = F, span = 1, show.legend = T) +
  geom_point(size = 0.8, alpha = 0.8, show.legend = F) +
  facet_wrap(facets = species ~ module)

ggplot(data = df_all_long, mapping = aes(x = humanized_log2_age_days, y = mean, color = species)) +
  geom_smooth(linewidth = 1, alpha = 0.8, se = F, span = 1, show.legend = T) +
  geom_point(size = 1, alpha = 0.8, show.legend = F) +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0.1, alpha = 0.8) +
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
       y = "Scaled abundance") +
  scale_x_continuous(breaks = c(6.8073549221,8.0552824355,9.30149619498,10.7532167492,12.885315061), labels = c("GW18", "Birth", "Year01", "Year04", "Year20"))

ggsave("Results/species_comparison/module_abundance_across_species.pdf", width = 6, height = 6)
ggsave("Results/species_comparison/module_abundance_across_species.png", width = 6, height = 6)

#plot RhoGEF proteins abundance in different species
RhoGEFs <- species_comparison[c("ARHGEF2", "ARHGEF7", "DOCK3", "NGEF", "RASGRF2"),]

df <- data.frame(assay(RhoGEFs)) %>%
  rownames_to_column() %>%
  gather(label, val, -rowname) %>%
  left_join(., data.frame(colData(RhoGEFs)), by = "label")
colnames(df)[1] <- "Protein"
#add PREX1
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

ggsave("Results/species_comparison/RhoGEF_abundance_across_species.pdf", width = 5.5, height = 7)
