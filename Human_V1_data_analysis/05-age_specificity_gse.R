library("clusterProfiler")
library("org.Hs.eg.db")
library("ReactomePA")
library("DEP")
library("enrichplot")
library("tidyverse")


#MSigDB
library(msigdbr)
msigdbr_species()
m_df = msigdbr(species = "Homo sapiens", category = "C2")
table(m_df$gs_subcat)
m_t2g <- c()
for (i in c("CP:BIOCARTA", "CP:KEGG", "CP:PID", "CP:REACTOME", "CP")) {
  m_df = msigdbr(species = "Homo sapiens", category = "C2", subcategory = i)
  tmp = m_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
  m_t2g = rbind(m_t2g, tmp)
}

m_t2g$gs_name <- tolower(m_t2g$gs_name)

dep_corrected_specificity = readRDS(file = "Results/dep_corrected_specificity.rds")


for (i in c("GW18_19", "GW22_23","GW28_41", "Year00_01", "Year04_08", "Year18_22")) {
  specificity <- dep_corrected_specificity@elementMetadata[,paste0(i , "_specificity_t")]
  names(specificity) <- dep_corrected_specificity@elementMetadata$name
  specificity <- sort(specificity, decreasing = T)
  
  res <- GSEA(specificity,
              exponent = 1,
              minGSSize = 10,
              maxGSSize = 500,
              pvalueCutoff = 1,
              pAdjustMethod = "BH",
              TERM2GENE = m_t2g,
              verbose = TRUE,
              seed = 0)
  
  assign(paste0(i, "_GESA"), res)
}

dotplot(GW18_19_GESA, split = ".sign") + facet_grid(.~.sign)
dotplot(GW22_23_GESA, split = ".sign") + facet_grid(.~.sign)
dotplot(GW28_41_GESA, split = ".sign") + facet_grid(.~.sign)
dotplot(Year00_01_GESA, split = ".sign") + facet_grid(.~.sign)
dotplot(Year04_08_GESA, split = ".sign") + facet_grid(.~.sign)
dotplot(Year18_22_GESA, split = ".sign") + facet_grid(.~.sign)

GESA <- list(GW18_19_GESA, GW22_23_GESA, GW28_41_GESA, Year00_01_GESA, Year04_08_GESA, Year18_22_GESA)
names(GESA) <- c("GW18_19", "GW22_23","GW28_41", "Year00_01", "Year04_08", "Year18_22")
GESA_res<- merge_result(GESA)
write.csv(GESA_res@compareClusterResult, "Results/age_specificity_GESA_comparison/specificity_GESA_res_human_BA17.csv")
save(GW18_19_GESA, GW22_23_GESA, GW28_41_GESA, Year00_01_GESA, Year04_08_GESA, Year18_22_GESA, file = "Results/age_specificity_GSEA_BA17.RData")

