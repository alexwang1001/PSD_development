library("DEP")
library("SummarizedExperiment")
library("limma")
library("vsn")
library("tidyverse")

load(file = "./Results/data.RData")
#how many proteins are found at each stage? These include by-matching identifications.
plot_missval(data_se)
bin_data <- assay(data_se)
idx <- is.na(assay(data_se))
bin_data[!idx] <- 1
bin_data[idx] <- 0

#only include proteins that are missing in less than 3 samples in any given stage
keep <- bin_data %>% data.frame() %>% rownames_to_column() %>% 
  gather(ID, value, -rowname) %>% left_join(., data.frame(colData(data_se)), 
                                            by = "ID") %>% group_by(rowname, condition) %>% summarize(miss_val = n() - 
                                                                                                        sum(value)) %>% filter(miss_val <= 3) %>% ungroup()

#remove ASD samples
groupNumber <- keep %>% dplyr::count(condition) %>% filter(condition != c("Year04_08_ASD", "Year18_22_ASD"))

ggplot(data = groupNumber, mapping = aes (x = condition, y = n)) +
  geom_col() +
  theme_bw() +
  theme(text=element_text(size=20)) +
  labs(x = "Developmental Stages",
       y = "Protein Identified")

data_half_filt <- data_se[unique(keep$rowname), ]
plot_missval(data_half_filt)
plot_frequency(data_half_filt)
plot_numbers(data_half_filt)
plot_coverage(data_half_filt)

vsn.fit <- vsn::vsnMatrix(2^assay(data_half_filt))
saveRDS(vsn.fit, file = "./Results/vsn_reference.rds")

#normalize data
data_norm <- normalize_vsn(data_half_filt)
meanSdPlot(assay(data_half_filt))
meanSdPlot(assay(data_norm))
plotMD(assay(data_half_filt))
plotMD(assay(data_norm))
plot_normalization(data_half_filt, data_norm)


#Imputation
set.seed(6)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
plot_imputation(data_norm, data_imp)


save(data_norm, file = "./Results/data_norm.RData")
save(data_imp, file = "./Results/data_imp.RData")


