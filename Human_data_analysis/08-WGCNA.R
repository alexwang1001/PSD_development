library("WGCNA")
library("DEP")
library("SummarizedExperiment")
library("limma")
library("tidyverse")
library("vsn")

load(file = "./Results/dep_corrected.RData")
# The following setting is important, do not omit.
options(stringsAsFactors = F)
datExpr <- t(assay(dep_corrected))

# Choose a set of soft-threshold powers
powers = c(1:40)
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, networkType = "signed", powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
#The negative "signed R^2" is negative when your network has more genes with
#high connectivity than ones with low connectivity
#(i.e., the regression line for the fit log(n(k))~log(k) has a positive slope).
#It means your network shows a topology in some ways opposite
#(more high connectivity than low connectivity genes) to what is normally expected
#(a lot of low-connetivity gens and fewer high connectivity genes).
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")




#single step module detection
cor <- WGCNA::cor
set.seed(0)
net <- blockwiseModules(datExpr, power = 20,
                        corType = "pearson", # or use robust correlation (bicor)
                        networkType = "signed", deepSplit = 2, minModuleSize = 45,
                        reassignThreshold = 1e-6, mergeCutHeight = 0.15,
                        numericLabels = F, pamRespectsDendro = F,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "TOM",
                        verbose = 3)

table(net$colors)
cor<-stats::cor
# Convert labels to colors for plotting
mergedColors = net$colors
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


MEs = net$MEs
MET <- orderMEs(MEs)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

#save kME info
kME <- signedKME(datExpr,MEs)
kME$SYMBOL <- rownames(kME)
kME$geneModules <- net$colors
write.csv(kME, "Results/WGCNA_module_abundance_patterns/kME.csv", row.names = F)
save(net,file = "./Results/WGCNA_net.RData")
