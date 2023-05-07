library("clusterProfiler")
library("tidyverse")
library("readxl")

#differentially expressed genes in ASD, SCZ, BPD, and MDD
ASD_diff <- read_excel(path = "Gandal_differential_expression/ASD_diff.xlsx", na = "NA")
ASD_diff <- ASD_diff[!is.na(ASD_diff$ENTREZID),]
diff_ID <- bitr(ASD_diff$ENTREZID,fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
diff_ID$ENTREZID <- as.numeric(diff_ID$ENTREZID)
ASD_diff <- dplyr::left_join(ASD_diff, diff_ID)
ASD_diff <- ASD_diff[!is.na(ASD_diff$SYMBOL),]
ASD_diff <- ASD_diff %>% filter(!is.na(FDR)) %>% mutate(`-log10(p.adj)` = -log(FDR,10), disorder = "ASD")


BPD_diff <- read_excel(path = "Gandal_differential_expression/BPD_diff.xlsx", na = "NA")
BPD_diff <- BPD_diff[!is.na(BPD_diff$ENTREZID),]
BPD_diff <- dplyr::left_join(BPD_diff, diff_ID)
BPD_diff <- BPD_diff[!is.na(BPD_diff$SYMBOL),]
BPD_diff <- BPD_diff %>% filter(!is.na(FDR)) %>% mutate(`-log10(p.adj)` = -log(FDR,10), disorder = "BPD")


MDD_diff <- read_excel(path = "Gandal_differential_expression/MDD_diff.xlsx", na = "NA")
MDD_diff <- MDD_diff[!is.na(MDD_diff$ENTREZID),]
MDD_diff <- dplyr::left_join(MDD_diff, diff_ID)
MDD_diff <- MDD_diff[!is.na(MDD_diff$SYMBOL),]
MDD_diff <- MDD_diff %>% filter(!is.na(FDR)) %>% mutate(`-log10(p.adj)` = -log(FDR,10), disorder = "MDD")


SCZ_diff <- read_excel(path = "Gandal_differential_expression/SCZ_diff.xlsx", na = "NA")
SCZ_diff <- SCZ_diff[!is.na(SCZ_diff$ENTREZID),]
SCZ_diff <- dplyr::left_join(SCZ_diff, diff_ID)
SCZ_diff <- SCZ_diff[!is.na(SCZ_diff$SYMBOL),]
SCZ_diff <- SCZ_diff %>% filter(!is.na(FDR)) %>% mutate(`-log10(p.adj)` = -log(FDR,10), disorder = "SCZ")

diff <- rbind(ASD_diff, BPD_diff, MDD_diff, SCZ_diff)
diff <- diff %>% mutate(significant = FDR < 0.05)
diff$disorder <- fct_relevel(diff$disorder, c("SCZ", "ASD", "BPD", "MDD"))

geneIDsModules <- read.csv(file = "Results/WGCNA_ORA/geneIDsModules.csv", header = TRUE, row.names = 1)
brown<- geneIDsModules %>% filter(geneModules == "brown") %>% .$SYMBOL
blue <- geneIDsModules %>% filter(geneModules == "blue") %>% .$SYMBOL
turquoise <- geneIDsModules %>% filter(geneModules == "turquoise") %>% .$SYMBOL
yellow <- geneIDsModules %>% filter(geneModules == "yellow") %>% .$SYMBOL

diff_brown <- diff %>% filter(SYMBOL %in% brown)
diff_blue <- diff %>% filter(SYMBOL %in% blue)
diff_turquoise <- diff %>% filter(SYMBOL %in% turquoise)
diff_yellow <- diff %>% filter(SYMBOL %in% yellow)

ggplot(diff_brown, aes(x = log2FC, y = `-log10(p.adj)`)) +
  xlim(-0.6, 0.6) +
  ylim(0,6) +
  geom_point(mapping = aes(color = significant), alpha = 0.5, size = 0.5) +
  facet_grid(cols = vars(disorder)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -log(0.05,10), linetype = "dashed") +
  xlab("Fold change (log2)") +
  ylab("p.adj (-log10)") +
  scale_color_manual(values = c("gray", "red")) +
  ggtitle("brown module") +
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(color = "beige", fill="beige"))

ggsave(filename = "Results/differential_expression_in_psychiatric_disorder/brown_volcano.jpeg", width = 6, height = 2, dpi = 300, units = "in")

ggplot(diff_blue, aes(x = log2FC, y = `-log10(p.adj)`)) +
  xlim(-0.6, 0.6) +
  ylim(0,6) +
  geom_point(mapping = aes(color = significant), alpha = 0.5, size = 0.5) +
  facet_grid(cols = vars(disorder)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -log(0.05,10), linetype = "dashed") +
  xlab("Fold change (log2)") +
  ylab("p.adj (-log10)") +
  scale_color_manual(values = c("gray", "red")) +
  ggtitle("blue module") +
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(color = "beige", fill="beige"))

ggsave(filename = "Results/differential_expression_in_psychiatric_disorder/blue_volcano.jpeg", width = 6, height = 2, dpi = 300, units = "in")

ggplot(diff_turquoise, aes(x = log2FC, y = `-log10(p.adj)`)) +
  xlim(-0.6, 0.6) +
  ylim(0,6) +
  geom_point(mapping = aes(color = significant), alpha = 0.5, size = 0.5) +
  facet_grid(cols = vars(disorder)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -log(0.05,10), linetype = "dashed") +
  xlab("Fold change (log2)") +
  ylab("p.adj (-log10)") +
  scale_color_manual(values = c("gray", "red")) +
  ggtitle("turquoise module") +
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(color = "beige", fill="beige"))

ggsave(filename = "Results/differential_expression_in_psychiatric_disorder/turquoise_volcano.jpeg", width = 6, height = 2, dpi = 300, units = "in")

ggplot(diff_yellow, aes(x = log2FC, y = `-log10(p.adj)`)) +
  xlim(-0.6, 0.6) +
  ylim(0,6) +
  geom_point(mapping = aes(color = significant), alpha = 0.5, size = 0.5) +
  facet_grid(cols = vars(disorder)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -log(0.05,10), linetype = "dashed") +
  xlab("Fold change (log2)") +
  ylab("p.adj (-log10)") +
  scale_color_manual(values = c("gray", "red")) +
  ggtitle("yellow module") +
  theme_bw() +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(color = "beige", fill="beige"))

ggsave(filename = "Results/differential_expression_in_psychiatric_disorder/yellow_volcano.jpeg", width = 6, height = 2, dpi = 300, units = "in")


ggplot(diff, aes(x = log2FC, y = `-log10(p.adj)`)) +
  xlim(-0.6, 0.6) +
  ylim(0,6) +
  geom_point(mapping = aes(color = significant), alpha = 0.5, size = 0.2) +
  facet_grid(cols = vars(disorder)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -log(0.05,10), linetype = "dashed") +
  xlab("Fold change (log2)") +
  ylab("p.adj (-log10)") +
  scale_color_manual(values = c("gray", "red")) +
  ggtitle("All genes") +
  theme_bw() +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(color = "beige", fill="beige"))

ggsave(filename = "Results/differential_expression_in_psychiatric_disorder/all_volcano.jpeg", width = 6, height = 2,dpi = 300, units = "in")
