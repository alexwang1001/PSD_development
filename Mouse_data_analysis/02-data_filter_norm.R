library("DEP")
library("SummarizedExperiment")
library("limma")
library("vsn")
library("tidyverse")
library("ComplexHeatmap")

load(file = "./Results/data.RData")
#how many proteins are found at each stage?
#assign all match between runs values to 0
data_unique <- read.csv(file = "Clean_Data/data_unique.csv")
identification_method <- read.csv(file = "Clean_Data/Ident_method.csv")
identification_method_unique <- identification_method %>% dplyr::semi_join(data_unique, by = "Protein.IDs")
matched_idx <- identification_method_unique[,3:22] == "By matching"

bin_data <- assay(data_se)
idx <- is.na(assay(data_se))
bin_data[!idx] <- 1
bin_data[idx] <- 0
sum(bin_data)
bin_data[matched_idx] <- 0

keep <- bin_data %>% data.frame() %>% rownames_to_column() %>% 
  gather(ID, value, -rowname) %>% left_join(., data.frame(colData(data_se)), 
                                            by = "ID") %>% group_by(rowname, condition) %>% summarize(miss_val = n() - 
                                                                                                        sum(value)) %>% filter(miss_val <= 2) %>% ungroup()

groupNumber <- keep %>% dplyr::count(condition)
groupNumber$condition <- factor(groupNumber$condition,levels = c("P0", "P9", "P13", "P18", "P36"))

plot1 <- ggplot(data = groupNumber, mapping = aes (x = condition, y = n)) +
  geom_col() +
  theme_bw() +
  theme(text=element_text(size=15)) +
  labs(x = "Developmental Stages",
       y = "Protein Number Identified")
ggsave(filename = "Results/protein_number_identified.png", width = 5, height = 6, units = "in")

keep2 <- keep %>% spread(key = condition, value = miss_val) %>% as.data.frame()
rownames(keep2) <- keep2$rowname
keep2 <- as.matrix(keep2[,-1])
keep2[!is.na(keep2)] <- 1
keep2[is.na(keep2)] <- 0
condition <- colnames(keep2)[c(1,5,2,3,4)]

#draw an UpSet plot
m <- make_comb_mat(keep2, mode = "distinct", min_set_size = 50)
m2 <- m[comb_size(m) > 30]

pdf("Results/upSet_MSMSonly.pdf") 
UpSet(m2, set_order = condition, comb_order = order(comb_size(m2), decreasing = TRUE), width = unit(2.5, "in"), height = unit(1, "in"))
dev.off()

#Next we performed similar filtering but include by-matching identifications.
#only keep proteins that are detected in at least half of all samples in a given stage
data_filt <- filter_missval(data_se, thr = 2)
dim(data_filt)
plot_frequency(data_filt)
plot_numbers(data_filt)
plot_coverage(data_filt)

#normalize data
data_norm <- normalize_vsn(data_filt)

#Imputation
set.seed(4)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
plot_imputation(data_norm, data_imp)

save(data_norm, file = "./Results/data_norm.RData")
save(data_imp, file = "./Results/data_imp.RData")

