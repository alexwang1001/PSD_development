library("tidyverse")


load(file = "./Results/WGCNA_net.RData")
load(file = "./Results/dep_corrected.RData")

MEs <- net$MEs
write.csv(MEs, "Results/WGCNA_module_abundance_patterns/MEs.csv")
MEs <- MEs[,colnames(MEs) != "MEgrey"]
colnames(MEs) <- c("blue", "brown", "turquoise", "yellow")
cols <- as.data.frame(colData(dep_corrected))
MEs <- cbind(MEs, age = cols$Log2AgeDays, sample = rownames(MEs))
MEs_col <- gather(MEs, key = "module", value = "ME", 1:4)

ggplot(data = MEs_col, mapping = aes(x = age, y = ME, color = module, fill = module)) +
  geom_smooth(se = F, show.legend = F) +
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
    y = "Scaled abundance") +
  scale_color_manual(values = c("blue", "brown","turquoise","gold")) +
  scale_fill_manual(values = c("blue", "brown","turquoise","gold")) +
  scale_x_continuous(breaks = c(6.8073549221,7.1996723448,8.0552824355,9.30149619498,10.7532167492,11.6375305515,12.885315061), labels = c("GW18", "GW23", "Birth", "Year01", "Year04", "Year08", "Year20"))

ggsave(filename = "Results/WGCNA_module_abundance_patterns/protein_modules_abundance_patterns.pdf", width = 6, height = 4, units = "in")
