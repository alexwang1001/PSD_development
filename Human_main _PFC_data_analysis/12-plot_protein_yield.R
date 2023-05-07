library("plotrix")

protein_yield <- read.csv("Results/protein_yield/protein_yield.csv", header = T)
protein_yield <- protein_yield[protein_yield$Disorder == "Control",]

lx <- protein_yield$Log2AgeDays
ly <- protein_yield$total.protein.yield

synapse_number <- read.csv("Results/protein_yield/PFC_synapse_number.csv", header = T)
rx <- synapse_number$Log2AgeDays
ry <- synapse_number$Synapses.100.um3

pdf("Results/protein_yield/protein_yield.pdf", width = 7, height = 5)
twoord.plot(lx, ly, rx, ry, type = c("p", "p"), xlim = c(6.5,15), lylim = c(0,2.8), rylim = c(0,70),
            xlab = "Age (log2-transformed)", ylab = "PSD yield (ug/mg tissue)", rylab = "Synapse number / 100 um3",
            xtickpos = c(6.8073549221, 8.0552824355, 9.30149619498, 10.7532167492, 12.885315061), xticklab = c("GW18", "Birth", "Year01", "Year04", "Year20"),
            lwd = 1.5, cex = 1.5)
legend("bottomright", c("PSD Yield", "Synapse Number"), pch = c(1,2), col = c(1,2),cex =1, pt.lwd = 1.5, text.width = 2.5,text.col = c(1,2))
dev.off()

