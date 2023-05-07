library("tidyverse")

load("NDD_mutations/NDD_de_novo_variants_all.RData")
geneIDsModules2 <- read.csv(file = "Results/WGCNA_ORA/geneIDsModules2.csv", header = TRUE, row.names = 1)
geneIDsModules <- read.csv(file = "Results/WGCNA_ORA/geneIDsModules.csv", header = TRUE, row.names = 1)
all <- geneIDsModules[,1:2] %>% mutate ("geneModules" = "all")
geneIDsModules3 <- rbind(geneIDsModules2, all)
write.csv(geneIDsModules3, file = "Results/NDD_mutations/geneIDsModules3.csv", row.names = FALSE)

NDD_variants$ENTREZID <- as.numeric(NDD_variants$ENTREZID)
NDD_variants_modules <- left_join(NDD_variants, geneIDsModules, by = c("SYMBOL","ENTREZID"))
write.csv(NDD_variants_modules, file = "Results/NDD_mutations/NDD_variants_modules.csv", row.names = FALSE)

res <- c()

for (color in unique(geneIDsModules3$geneModules)) {
  genes <- geneIDsModules3$SYMBOL[geneIDsModules3$geneModules == color]
  in_module <- NDD_variants %>% group_by(PrimaryPhenotype, MutationType) %>%
    filter((SYMBOL %in% genes)) %>% summarise(in_module_count = sum(n)) %>%
    ungroup()
  not_in_module <- NDD_variants %>% group_by(PrimaryPhenotype, MutationType) %>%
    filter(!(SYMBOL %in% genes)) %>% summarise(not_in_module_count = sum(n)) %>%
    ungroup()
  module_count <- full_join(in_module, not_in_module, by= c("PrimaryPhenotype", "MutationType")) %>%
    complete(PrimaryPhenotype, MutationType) %>%
    mutate(module = color) %>% 
    dplyr::select(module, everything())
  module_count[is.na(module_count)] <- 0
  res <- rbind(res, module_count)
}
module_gene_count <- res

#complete() will make NA value explicit which is very important for the following fisher exact test
#alternatively we can use expand.grid from base R in the for loop:
# module_count <- expand.grid(unique(NDD_variants$PrimaryPhenotype), unique(NDD_variants$MutationType))
# colnames(module_count) <- c("PrimaryPhenotype", "MutationType")
# module_count <- module_count %>% left_join(in_module, by = c("PrimaryPhenotype", "MutationType")) %>% 
#   left_join(not_in_module, by = c("PrimaryPhenotype", "MutationType")) %>% 
#   mutate(module = color) %>% 
#   select(module, everything())
# module_count[is.na(module_count)] <- 0

#fisher exact test
test_FET_enrich <- function(module, PrimaryPhenotype, MutationType) {
  df <- matrix(as.numeric(c(module_gene_count[module_gene_count$module == module & module_gene_count$PrimaryPhenotype == PrimaryPhenotype & module_gene_count$MutationType == MutationType,][4:5],
                            module_gene_count[module_gene_count$module == module & module_gene_count$PrimaryPhenotype == "control" & module_gene_count$MutationType == MutationType,][4:5])), nrow = 2, byrow = TRUE)
  df[is.na(df)] <- 0
  res <- fisher.test(df, alternative = "greater")
  res2 <- c("module" = module, 
            "PrimaryPhenotype" = PrimaryPhenotype, 
            "MutationType" = MutationType, 
            "p_value" = res$p.value, 
            "odds_ratio" = res$estimate,
            "count_in_module_case" = df[1,1],
            "count_not_in_module_case" = df[1,2],
            "count_in_module_control" = df[2,1],
            "count_not_in_module_control" = df[2,2],
            "total_in_case" = df[1,1] + df[1,2],
            "total_in_control" = df[2,1] + df[2,2],
            "fraction_in_case" = df[1,1]/(df[1,1] + df[1,2]),
            "fraction_in_control" = df[2,1]/(df[2,1] + df[2,2])
  )
}


param <- list("module" = as.character(module_gene_count$module), "PrimaryPhenotype" = as.character(module_gene_count$PrimaryPhenotype), "MutationType" = as.character(module_gene_count$MutationType))
res <- pmap(param, test_FET_enrich)
res_df <- data.frame(matrix(unlist(res), nrow=length(res), byrow=T))

colnames(res_df) <- c("module", "PrimaryPhenotype", "MutationType", "p_value", "odds_ratio", "count_in_module_case", "count_not_in_module_case", "count_in_module_control", "count_not_in_module_control", "total_in_case", "total_in_control", "fraction_in_case", "fraction_in_control")
res_df$p_value <- as.numeric(as.character(res_df$p_value))
res_df$odds_ratio <- as.numeric(as.character(res_df$odds_ratio))
res_df <- unite(res_df, disorder, 2:3)

res_list <- res_df %>% group_split(module)
res_list <- lapply(res_list, as.data.frame)
for (i in 1:5) {
  res_list[[i]]$p.adj <- as.character(p.adjust(res_list[[i]]$p_value, method = "BH"))
}

res_list <- bind_rows(res_list)

#save data
write.csv(res_list, file = "Results/NDD_mutations/NDD_de_novo_mutations_burden.csv", row.names = FALSE)
