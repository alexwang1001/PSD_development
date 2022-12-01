library("clusterProfiler")
library("ReactomePA")
library("DEP")
library("enrichplot")
library("tidyverse")
library("msigdbr")


m_df = msigdbr(species = "Macaca mulatta", category = "C2")
table(m_df$gs_subcat)
m_t2g <- c()
for (i in c("CP:BIOCARTA", "CP:KEGG", "CP:PID", "CP:REACTOME", "CP")) {
  m_df = msigdbr(species = "Macaca mulatta", category = "C2", subcategory = i)
  tmp = m_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
  m_t2g = rbind(m_t2g, tmp)
}

m_t2g$gs_name <- tolower(m_t2g$gs_name)

load(file = "./Results/dep_specificity.RData")


for (i in c("E75_85", "E95_110", "Year0_1", "Year1_3", "Year8_10")) {
  specificity <- dep_specificity@elementMetadata[,paste0(i , "_specificity_t")]
  names(specificity) <- dep_specificity@elementMetadata$name
  specificity <- sort(specificity, decreasing = T)
  
  res <- GSEA(specificity,
              exponent = 1,
              minGSSize = 10,
              maxGSSize = 500,
              eps = 0,
              pvalueCutoff = 1,
              pAdjustMethod = "BH",
              TERM2GENE = m_t2g,
              verbose = TRUE,
              seed = 0)
  
  assign(paste0(i, "_GESA"), res)
}

GESA <- list(E75_85_GESA, E95_110_GESA, Year0_1_GESA, Year1_3_GESA, Year8_10_GESA)
names(GESA) <- c("E75_85", "E95_110", "Year0_1", "Year1_3", "Year8_10")
GESA_res<- merge_result(GESA)
write.csv(GESA_res@compareClusterResult, "Results/age_specificity/specificity_GESA_res.csv")
save(E75_85_GESA, E95_110_GESA, Year0_1_GESA, Year1_3_GESA, Year8_10_GESA, file = "Results/age_specificity_GSEA.RData")
