library("DEP")
library("SummarizedExperiment")
library("limma")
library("vsn")
library("tidyverse")

dep_corrected_specificity <- readRDS("Results/dep_corrected_BA17.rds")

for (i in c("GW18_19", "GW22_23", "GW28_41", "Year00_01", "Year04_08", "Year18_22")) {
  #using limma
  condition <- factor(dep_corrected_specificity$condition == i)
  design <- model.matrix(terms(~ condition,keep.order = TRUE))
  colnames(design)[2] <- i
  fit <- lmFit(assay(dep_corrected_specificity), design=design)
  fit <- eBayes(fit)
  #save comparison data
  cntrst <- colnames(design)[2]
  retrieve_fun <- function(comp, fit) {
    res <- topTable(fit, sort.by = "p", coef = comp, number = Inf, 
                    confint = TRUE)
    res <- res[!is.na(res$P.Value), ]
    res$comparison <- rep(comp , dim(res)[1])
    res <- rownames_to_column(res)
    return(res)
  }
  
  limma_res <- map_df(cntrst, retrieve_fun, fit = fit)
  table <- limma_res %>% select(rowname, logFC, CI.L, CI.R, t, 
                                P.Value, adj.P.Val, comparison) %>% mutate(comparison = paste0(comparison,
                                                                                               "_specificity")) %>% gather(variable, value, -c(rowname, 
                                                                                                                                               comparison)) %>% mutate(variable = recode(variable, logFC = "diff", 
                                                                                                                                                                                         P.Value = "p.val", adj.P.Val = "p.adj")) %>% unite(temp, comparison, 
                                                                                                                                                                                                                                            variable) %>% spread(temp, value)
  
  
  rowData(dep_corrected_specificity) <- merge(rowData(dep_corrected_specificity, use.names = FALSE), table, 
                                              by.x = "name", by.y = "rowname", all.x = TRUE, sort = FALSE)
}

rowdata <- as.data.frame(rowData(dep_corrected_specificity))
rowdata <- rowdata[,sort(c(grep("_diff", colnames(rowdata)), grep("_adj", colnames(rowdata)), c(1,2,3,4,114), grep("specificity_t", colnames(rowdata))))]

saveRDS(dep_corrected_specificity, file = "Results/dep_corrected_specificity.rds")
