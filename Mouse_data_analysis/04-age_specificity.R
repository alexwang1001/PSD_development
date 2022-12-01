library("DEP")
library("SummarizedExperiment")
library("limma")
library("vsn")
library("tidyverse")

load(file = "./Results/dep.RData")
dep_specificity <- dep


for (i in c("P0", "P9", "P13", "P18", "P36")) {
  #using limma
  condition <- factor(dep_specificity$condition == i)
  design <- model.matrix(terms(~ condition, keep.order = TRUE))
  colnames(design)[2] <- i
  fit <- lmFit(assay(dep_specificity), design=design)
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
  
  
  rowData(dep_specificity) <- merge(rowData(dep_specificity, use.names = FALSE), table, 
                                              by.x = "name", by.y = "rowname", all.x = TRUE, sort = FALSE)
}

save(dep_specificity, file = "./Results/dep_specificity.RData")
