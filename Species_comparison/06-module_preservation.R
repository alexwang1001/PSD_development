library("WGCNA")
library("SummarizedExperiment")

species_comparison <- readRDS("Results/species_comparison_humanized_age.rds")

dim(assay(species_comparison))
datExpr=data.frame(t(assay(species_comparison)))

# Number of data sets that we work with
nSets = 3
# Object that will contain the expression data
multiExpr = list()
multiExpr[[1]] = list(data = datExpr[1:54,])
multiExpr[[2]] = list(data = datExpr[55:75,])
multiExpr[[3]] = list(data = datExpr[76:95,])
# Names for the two sets
setLabels = c("Human", "Macaque", "Mouse")
# Important: components of multiExpr must carry identificating names
names(multiExpr) = setLabels
# Display the dimensions of the expression data (if you are confused by this construct, ignore it):
lapply(multiExpr, lapply, dim)

# Create an object (list) holding the module labels for each set:
rowData <- data.frame(rowData(species_comparison))
colorList = list(rowData$Human_Module, rowData$Human_Module, rowData$Human_Module);
# Components of the list must be named so that the names can be matched to the names of multiExpr
names(colorList) = setLabels

system.time( {
  mp = modulePreservation(multiExpr, colorList,
                          networkType = "signed",
                          referenceNetworks = 1,
                          loadPermutedStatistics = FALSE,
                          nPermutations = 500,
                          verbose = 3)
} )

# This variable will contain the summary table
summaryTable = NULL
# Loop over all combinations of reference and tests sets
ref = 1;
test = c(2,3);
for (ref in c(1)) for (test in c(2:3)) if (ref!=test)
{
  modules = rownames(mp$preservation$Z[[ref]][[test]]);
  nMods = length(modules);
  sizes = mp$preservation$Z[[ref]][[test]][, 1];
  acc = matrix(NA, nMods, 3);
  if (test!=4)
  {
    acc[match(rownames(mp$accuracy$observed[[ref]][[test]]), modules), ] =
      mp$accuracy$observed[[ref]][[test]][, -1, drop = FALSE];
    colnames(acc) = colnames(mp$accuracy$observed[[ref]][[test]])[-1];
    accZ = mp$accuracy$Z[[ref]][[test]][, -1, drop = FALSE];
    acc.log.p = mp$accuracy$log.p[[ref]][[test]][, -1, drop = FALSE];
    acc.log.pBonf = mp$accuracy$log.pBonf[[ref]][[test]][, -1, drop = FALSE];
  } else {
    accZ = matrix(NA, nMods, 3);
    acc.log.p = matrix(NA, nMods, 3);
    acc.log.pBonf = matrix(NA, nMods, 3);
    colnames(acc) = colnames(mp$accuracy$observed[[1]][[2]])[-1];
    colnames(accZ) = colnames(mp$accuracy$Z[[1]][[2]])[-1];
    colnames(acc.log.p) = colnames(mp$accuracy$log.p[[1]][[2]])[-1];
    colnames(acc.log.pBonf) = colnames(mp$accuracy$log.pBonf[[1]][[2]])[-1];
  }
  # Table of results for this reference-test combination
  tab = cbind(referenceSet = rep(setLabels[ref], nMods),
              testSet = rep(setLabels[test], nMods),
              moduleLabel = modules,
              moduleSize = sizes,
              mp$quality$observed[[ref]][[test]][, -1, drop = FALSE],
              mp$preservation$observed[[ref]][[test]][, -1, drop = FALSE],
              acc,
              mp$referenceSeparability$observed[[ref]][[test]][, -1, drop = FALSE],
              mp$testSeparability$observed[[ref]][[test]][, -1, drop = FALSE],
              mp$quality$Z[[ref]][[test]][, -1, drop = FALSE],
              mp$quality$log.p[[ref]][[test]][, -1, drop = FALSE],
              mp$quality$log.pBonf[[ref]][[test]][, -1, drop = FALSE],
              mp$preservation$Z[[ref]][[test]][, -1, drop = FALSE],
              mp$preservation$log.p[[ref]][[test]][, -1, drop = FALSE],
              mp$preservation$log.pBonf[[ref]][[test]][, -1, drop = FALSE],
              accZ,
              acc.log.p,
              acc.log.pBonf,
              mp$referenceSeparability$Z[[ref]][[test]][, -1, drop = FALSE],
              mp$referenceSeparability$log.p[[ref]][[test]][, -1, drop = FALSE],
              mp$referenceSeparability$log.pBonf[[ref]][[test]][, -1, drop = FALSE],
              mp$testSeparability$Z[[ref]][[test]][, -1, drop = FALSE],
              mp$testSeparability$log.p[[ref]][[test]][, -1, drop = FALSE],
              mp$testSeparability$log.pBonf[[ref]][[test]][, -1, drop = FALSE]
  )
  # Add the table to the main table.
  if (is.null(summaryTable)) summaryTable = tab else summaryTable = rbind(summaryTable, tab)
}
# Save the table in csv format.
write.table(summaryTable, file = "Results/species_comparison/modulePreservationStatistics.csv", row.names = FALSE,
            sep = ",", quote = FALSE)
save(setLabels, multiExpr, colorList, mp, file = "Results/module_preservation.RData")
