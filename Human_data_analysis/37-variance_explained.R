library("variancePartition")

load(file = "./Results/data_imp.RData")
data_ctrl <- data_imp[,data_imp$Disorder == "Control"]
colData(data_ctrl)<- droplevels(colData(data_ctrl))
protein_matrix <- assay(data_ctrl)
info <- data.frame(colData(data_ctrl))
info$Degradation <- as.character(info$Degradation)

form <- ~ AgeGroup + Degradation + Sex + ProcessTime
varPart <- fitExtractVarPartModel(protein_matrix, form, info)
vp <- sortCols(varPart)
plotPercentBars(vp[1:10,])
colnames(vp) <- c("Age group", "Processing batch", "Quality", "Sex", "Residuals")
median(vp$age_group)
plotVarPart(vp)
ggsave("Results/variancePartition.pdf", width = 3, height = 2.4)
