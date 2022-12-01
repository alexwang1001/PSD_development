library("SummarizedExperiment")
library("tidyverse")
library("ggrepel")
library("WGCNA")

load("Results/dep_corrected.RData")
protein <- dep_corrected
protein_Expr <- data.frame(t(assay(protein)))

RNA <- readRDS("Results/Human_NCX_subset_transcriptomics_PSD_logTPM.rds")
RNA_Expr <- data.frame(t(assay(RNA)))



# Number of data sets that we work with
nSets = 2
# Object that will contain the expression data
multiExpr = list()
multiExpr[[1]] = list(data = protein_Expr)
multiExpr[[2]] = list(data = RNA_Expr)
# Names for the two sets
setLabels = c("Protein", "RNA")
# Important: components of multiExpr must carry identificating names
names(multiExpr) = setLabels
# Display the dimensions of the expression data (if you are confused by this construct, ignore it):
lapply(multiExpr, lapply, dim)

# Create an object (list) holding the module labels for each set:
#protein
rowData_protein <- data.frame(rowData(protein))
rowData_protein$ENTREZID <- as.numeric(rowData_protein$ENTREZID)
geneIDsModules <- read.csv(file = "Results/WGCNA_ORA/geneIDsModules.csv", row.names = 1)
rowData_protein <- left_join(rowData_protein, geneIDsModules)
rowData_RNA <- data.frame(rowData(RNA_filt))
colorList = list(rowData_protein$geneModules, rowData_RNA$geneModules);
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
test = 2;
for (ref in c(1)) for (test in c(2)) if (ref!=test)
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
write.table(summaryTable, file = "Results/RNA-seq/modulePreservationStatistics.csv", row.names = FALSE,
            sep = ",", quote = FALSE)
save(setLabels, multiExpr, colorList, mp, file = "Results/module_preservation_between_protein_and_RNA.RData")

#plot preservation median.rank and Z.summary
stats <- summaryTable[c(1,2,4:6),]
stats$medianRank.pres <- c(2,4,5,3,1)

ggplot(data = stats, mapping = aes(x = medianRank.pres, y = Zsummary.pres, col = moduleLabel,  label = moduleLabel)) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_label_repel(size =4, point.padding = 0.3, show.legend = FALSE) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "gray90", size = 0.25),
        panel.grid.major.x = element_line(colour = "gray90", size = 0.25),
        panel.border = element_rect(fill = NA),
        panel.background = element_rect(fill = "gray98", color = "black"),
        axis.text=element_text(size=12, color = "black"),
        axis.title = element_text(size=14)
  ) +
  labs(x = "Preservation statistic median rank",
       y = expression("Preservation "~Z[summary])) +
  geom_hline(yintercept = 2, linetype = 2, col = "red") +
  scale_x_reverse(limits = c(6,0), breaks = c(0,2,4,6)) +
  scale_y_continuous(limits = c(-1, 20), breaks = c(0,2,5,10,15,20)) +
  scale_color_manual(values = c("blue", "brown","grey60","turquoise","gold"))

ggsave(filename = "Results/RNA-seq/module_preservation_rank.pdf", width = 4, height = 3, units = "in")
