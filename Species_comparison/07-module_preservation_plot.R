library(tidyverse)
library(ggrepel)

#Zsummary and medianRank.pres comparison between humans and macaques
stats <- read.csv(file = "Results/species_comparison/modulePreservationStatistics.csv")
stats_macaque <- stats[c(1,2,4:6),]
stats_macaque$medianRank.pres <- c(2,3,5,4,1)

ggplot(data = stats_macaque, mapping = aes(x = medianRank.pres, y = Zsummary.pres, col = moduleLabel,  label = moduleLabel)) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_label_repel(size =4, point.padding = 0.3, show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "gray90", size = 0.25),
        panel.grid.major.x = element_line(colour = "gray90", size = 0.25),
        panel.border = element_rect(fill = NA),
        panel.background = element_rect(fill = "gray98", color = "black"),
        axis.text=element_text(size=12, color = "black"),
        axis.title = element_text(size=14)
  ) +
  labs(x = "Preservation statistic median rank",
       y = expression("Preservation "~Z[summary])) +
  geom_hline(yintercept = 2, linetype = 2, col = "red") +
  scale_x_reverse(limits = c(6,0), breaks = c(0,2,4,6)) +
  scale_y_continuous(limits = c(-1, 20), breaks = c(0,2,5,10,15,20)) +
  scale_color_manual(values = c("blue", "brown","grey60","turquoise","gold"))

ggsave(filename = "Results/species_comparison/module_preservation_rank_macaque.pdf", width = 4, height = 3, units = "in")

#Zsummary and medianRank.pres comparison between humans and mice
stats_mouse <- stats[c(7,8,10:12),]
stats_mouse$medianRank.pres <- c(2,4,5,3,1)

ggplot(data = stats_mouse, mapping = aes(x = medianRank.pres, y = Zsummary.pres, col = moduleLabel,  label = moduleLabel)) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_label_repel(size =4, point.padding = 0.3, show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "gray90", size = 0.25),
        panel.grid.major.x = element_line(colour = "gray90", size = 0.25),
        panel.border = element_rect(fill = NA),
        panel.background = element_rect(fill = "gray98", color = "black"),
        axis.text=element_text(size=12, color = "black"),
        axis.title = element_text(size=14)
  ) +
  labs(x = "Preservation statistic median rank",
       y = expression("Preservation "~Z[summary])) +
  geom_hline(yintercept = 2, linetype = 2, col = "red") +
  scale_x_reverse(limits = c(6,0), breaks = c(0,2,4,6)) +
  scale_y_continuous(limits = c(-2, 20), breaks = c(0,2,5,10,15,20)) +
  scale_color_manual(values = c("blue", "brown","grey60","turquoise","gold"))

ggsave(filename = "Results/species_comparison/module_preservation_rank_mouse.pdf", width = 4, height = 3, units = "in")
