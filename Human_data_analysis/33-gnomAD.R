library("tidyverse")
library("biomaRt")
library("patchwork")
library("FSA")

gnomAD <- read.table(file = "gnomAD/gnomad.v2.1.1.lof_metrics.by_gene.txt", sep = "\t", header = T)

listMarts()
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)
attributes = listAttributes(ensembl)

res <- getBM(attributes = c("ensembl_transcript_id", "entrezgene_id"),
      filters = "ensembl_transcript_id", 
      values = gnomAD$transcript, 
      mart = ensembl)
colnames(res) <- c("transcript", "ENTREZID")

gnomAD_2 <- left_join(gnomAD, res, by = "transcript")
gnomAD_2 <- gnomAD_2[,c(78,1:77)]
gnomAD_3 <- gnomAD_2[!is.na(gnomAD_2$ENTREZID),]


geneIDmodules <- read.csv(file = "Results/WGCNA_ORA/geneIDsModules.csv", row.names = 1)
geneIDmodules_gnomAD <- left_join(geneIDmodules, gnomAD_3, by = "ENTREZID")
geneIDmodules_gnomAD$SYMBOL[duplicated(geneIDmodules_gnomAD$ENTREZID)]
#correct some genes which ENTREZID is not matched with emsembl transcript and those duplicated ("NEBL","PALM2AKAP2", and "SLC4A1")
write.csv(geneIDmodules_gnomAD, file = "Results/gnomAD/geneIDmodules_gnomAD.csv", row.names = F)
geneIDmodules_gnomAD2 <- read.csv("Results/gnomAD/geneIDmodules_gnomAD_2.csv")
geneIDmodules_gnomAD_all <- geneIDmodules_gnomAD2
geneIDmodules_gnomAD_all$geneModules <- "PSD_all"
geneIDmodules_gnomAD3 <- rbind(geneIDmodules_gnomAD2, geneIDmodules_gnomAD_all)

#add data for all background genes expressed in the cortex
backgroundIDs <- read.csv("Results/WGCNA_ORA/backgroundIDs.csv", row.names = 1)
backgroundIDs$gene <- backgroundIDs$SYMBOL
backgroundIDs <- backgroundIDs[,c(2,1,3)]
backgroundIDs_gnomAD <- left_join(backgroundIDs, gnomAD, by="gene")
backgroundIDs_gnomAD <- backgroundIDs_gnomAD[order(backgroundIDs_gnomAD$exp_mis, decreasing = T),]
backgroundIDs_gnomAD <- distinct(backgroundIDs_gnomAD, SYMBOL, .keep_all = T)                               
backgroundIDs_gnomAD$geneModules <- "background"
backgroundIDs_gnomAD2 <- backgroundIDs_gnomAD[,c(1,2,80,3:79)]

gnomAD_final <- rbind(geneIDmodules_gnomAD3, backgroundIDs_gnomAD2)
gnomAD_final <- gnomAD_final %>% drop_na(geneModules, oe_lof_upper)
#boxplot
#LOEUF & mis3
gnomAD_final$geneModules <- factor(gnomAD_final$geneModules, levels = c("brown", "blue", "turquoise", "yellow", "grey", "PSD_all", "background"))

p1 <- ggplot(data = gnomAD_final, mapping = aes(x = geneModules, y = oe_lof_upper, color = geneModules)) +
  geom_boxplot(notch = TRUE, outlier.alpha = 0.25) +
  geom_hline(yintercept = median(unlist(gnomAD_final %>% filter(geneModules == "background") %>% dplyr::select(oe_lof_upper))), linetype="dashed") +
  scale_colour_manual(values = c("brown", "blue", "turquoise","gold","grey","red","black")) +
  ylim(0, 2.5) +
  theme_classic() +
  labs(x = "Category", y = "LOEUF score", color = "Category") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

p2 <- ggplot(data = gnomAD_final, mapping = aes(x = geneModules, y = mis_z, color = geneModules)) +
  geom_boxplot(notch = TRUE, outlier.alpha = 0.25) +
  geom_hline(yintercept = median(unlist(gnomAD_final %>% filter(geneModules == "background") %>% dplyr::select(mis_z))), linetype="dashed") +
  scale_colour_manual(values = c("brown", "blue", "turquoise","gold","grey","red","black")) +
  ylim(-5, 15) +
  theme_classic() +
  labs(x = "Category", y = "missense Z-score", color = "Category") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

p1+p2
ggsave("Results/gnomAD/LOEUF_mis_z.pdf", width = 8,height = 4)

#synonymous
ggplot(data = gnomAD_final, mapping = aes(x = geneModules, y = syn_z, color = geneModules)) +
  geom_boxplot(notch = TRUE, outlier.alpha = 0.25) +
  geom_hline(yintercept = median(unlist(gnomAD_final %>% filter(geneModules == "background") %>% dplyr::select(syn_z))), linetype="dashed") +
  scale_colour_manual(values = c("brown", "blue", "turquoise","gold","grey","red","black")) +
  ylim(-10, 5) +
  theme_classic() +
  labs(x = "Category", y = "synonymous Z-score", color = "Category") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("Results/gnomAD/LOEUF_syn_z.pdf", width = 4,height = 4)

kruskal.test(oe_lof_upper ~ geneModules, data = gnomAD_final)
dunnTest(oe_lof_upper ~ geneModules, data=gnomAD_final, method = "bh")  
pairwise.wilcox.test(gnomAD_final$oe_lof_upper, gnomAD_final$geneModules,
                     p.adjust.method = "BH")

kruskal.test(mis_z ~ geneModules, data = gnomAD_final)
dunnTest(mis_z ~ geneModules, data=gnomAD_final, method = "bh") 
pairwise.wilcox.test(gnomAD_final$mis_z, gnomAD_final$geneModules,
                     p.adjust.method = "BH")

kruskal.test(syn_z ~ geneModules, data = gnomAD_final)
dunnTest(syn_z ~ geneModules, data=gnomAD_final, method = "bh") 
pairwise.wilcox.test(gnomAD_final$syn_z, gnomAD_final$geneModules,
                     p.adjust.method = "BH")
