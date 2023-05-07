library("tidyverse")

###########################################################################

## This code is to calculate the shortest path lengths within the axon pathway proteins as well as those between the axon pathway proteins and non-pathway proteins in the blue module network.

## Extract the axon pathway proteins from the blue module.

library(igraph)

data_pathway <- read.table("blue_pathway.csv", sep = ",", colClasses = "character", skip = 1)

data_pathway <- data_pathway[data_pathway[,3]=="1",]

data_protein <- data_pathway[,1]

data_protein <- unique(data_protein)

data_protein <- sort(data_protein)

write.table(data_protein, file = "axon.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

## Generate the largest connected component of the blue module network.

data_network <- read.table("blue_network.csv", sep = ",", colClasses = "character")

data_network <- unique(data_network)

data_network <- data_network[data_network[,1]!=data_network[,2],]

data_network_1 <- data_network[data_network[,1]<data_network[,2],]

data_network_2 <- data_network[data_network[,1]>data_network[,2],]

data_network_3 <- data_network_2[,2:1]

write.table(data_network_1, file = "blue_network_1.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(data_network_3, file = "blue_network_1.txt", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

data_network_4 <- read.table("blue_network_1.txt", sep = "\t", colClasses = "character")

data_network_4 <- unique(data_network_4)

data_network_4 <- data_network_4[ do.call(order, data_network_4) ,]

write.table(data_network_4, file = "blue_network_2.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

data_network <- read.table("blue_network_2.txt", sep = "\t", colClasses = "character")

data_network_1 <- graph.data.frame(data_network, directed = FALSE)

data_cluster <- clusters(data_network_1)

data_LCC <- which(which.max(data_cluster$csize)==data_cluster$membership)

data_network_2 <- induced.subgraph(data_network_1, data_LCC)

write.graph(data_network_2, "blue_network_3.txt", format = "ncol")

## Calculate the shortest path lengths within the axon pathway proteins as well as those between the axon pathway proteins and non-pathway proteins in the largest connected component.

data_network <- read.table("blue_network_3.txt", sep = " ", colClasses = "character")

data_network_1 <- graph.data.frame(data_network, directed = FALSE)

data_protein_network <- scan("blue_network_3.txt", what = "character", sep = " ")

data_protein_pathway <- scan("axon.txt", what = "character", sep = "\t")

data_protein_pathway <- intersect(data_protein_pathway, data_protein_network)

data_protein_pathway <- unique(data_protein_pathway)

data_protein_pathway <- sort(data_protein_pathway)

data_protein_network <- setdiff(data_protein_network, data_protein_pathway)

data_protein_network <- unique(data_protein_network)

data_protein_network <- sort(data_protein_network)

data_SPL <- shortest.paths(data_network_1, v=data_protein_pathway, to=data_protein_pathway)

write.table(data_SPL, file = "SPL_axon_axon.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

data_SPL_1 <- shortest.paths(data_network_1, v=data_protein_pathway, to=data_protein_network)

write.table(data_SPL_1, file = "SPL_axon_others.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

###############################################################################################
#brown module translation
brown_translation <- read.table(file = "Results/PPI_network/shortest_path_length/brown_SPL_translation_translation.txt",
                                sep = "\t",
                                header = T,
                                row.names = 1)
brown_translation_list <- as.vector(as.matrix(brown_translation[lower.tri(brown_translation)]))
brown_translation_others <- read.table(file = "Results/PPI_network/shortest_path_length/brown_SPL_translation_others.txt",
                           sep = "\t",
                           header = T,
                           row.names = 1)
brown_translation_others_list <- as.vector(as.matrix(brown_translation_others))
brown_translation_res <- wilcox.test(brown_translation_list, brown_translation_others_list, alternative = "less")
brown_translation_res$mean1 <- mean(brown_translation_list)
brown_translation_res$mean2 <- mean(brown_translation_others_list)
brown_translation_res <- unlist(brown_translation_res)
#significant

#blue module translation
blue_translation <- read.table(file = "Results/PPI_network/shortest_path_length/blue_SPL_translation_translation.txt",
                               sep = "\t",
                               header = T,
                               row.names = 1)
blue_translation_list <- as.vector(as.matrix(blue_translation[lower.tri(blue_translation)]))
blue_translation_others <- read.table(file = "Results/PPI_network/shortest_path_length/blue_SPL_translation_others.txt",
                                       sep = "\t",
                                       header = T,
                                       row.names = 1)
blue_translation_others_list <- as.vector(as.matrix(blue_translation_others))
blue_translation_res <- wilcox.test(blue_translation_list, blue_translation_others_list, alternative = "less")
blue_translation_res$mean1 <- mean(blue_translation_list)
blue_translation_res$mean2 <- mean(blue_translation_others_list)
blue_translation_res <- unlist(blue_translation_res)
#significant

#blue module rho
blue_rho <- read.table(file = "Results/PPI_network/shortest_path_length/blue_SPL_rho_rho.txt",
                        sep = "\t",
                        header = T,
                        row.names = 1)
blue_rho_list <- as.vector(as.matrix(blue_rho[lower.tri(blue_rho)]))
blue_rho_others <- read.table(file = "Results/PPI_network/shortest_path_length/blue_SPL_rho_others.txt",
                               sep = "\t",
                               header = T,
                               row.names = 1)
blue_rho_others_list <- as.vector(as.matrix(blue_rho_others))
blue_rho_res <- wilcox.test(blue_rho_list, blue_rho_others_list, alternative = "less")
blue_rho_res$mean1 <- mean(blue_rho_list)
blue_rho_res$mean2 <- mean(blue_rho_others_list)
blue_rho_res <- unlist(blue_rho_res)
#significant

#turquoise module rho
turquoise_rho <- read.table(file = "Results/PPI_network/shortest_path_length/turquoise_SPL_rho_rho.txt",
                               sep = "\t",
                               header = T,
                               row.names = 1)
turquoise_rho_list <- as.vector(as.matrix(turquoise_rho[lower.tri(turquoise_rho)]))
turquoise_rho_others <- read.table(file = "Results/PPI_network/shortest_path_length/turquoise_SPL_rho_others.txt",
                                      sep = "\t",
                                      header = T,
                                      row.names = 1)
turquoise_rho_others_list <- as.vector(as.matrix(turquoise_rho_others))
turquoise_rho_res <- wilcox.test(turquoise_rho_list, turquoise_rho_others_list, alternative = "less")
turquoise_rho_res$mean1 <- mean(turquoise_rho_list)
turquoise_rho_res$mean2 <- mean(turquoise_rho_others_list)
turquoise_rho_res <- unlist(turquoise_rho_res)
#significant

#turquoise module postsynaptic
turquoise_postsynaptic <- read.table(file = "Results/PPI_network/shortest_path_length/turquoise_SPL_postsynaptic_postsynaptic.txt",
                            sep = "\t",
                            header = T,
                            row.names = 1)
turquoise_postsynaptic_list <- as.vector(as.matrix(turquoise_postsynaptic[lower.tri(turquoise_postsynaptic)]))
turquoise_postsynaptic_others <- read.table(file = "Results/PPI_network/shortest_path_length/turquoise_SPL_postsynaptic_others.txt",
                                   sep = "\t",
                                   header = T,
                                   row.names = 1)
turquoise_postsynaptic_others_list <- as.vector(as.matrix(turquoise_postsynaptic_others))
turquoise_postsynaptic_res <- wilcox.test(turquoise_postsynaptic_list, turquoise_postsynaptic_others_list, alternative = "less")
turquoise_postsynaptic_res$mean1 <- mean(turquoise_postsynaptic_list)
turquoise_postsynaptic_res$mean2 <- mean(turquoise_postsynaptic_others_list)
turquoise_postsynaptic_res <- unlist(turquoise_postsynaptic_res)
#significant

#yellow module postsynaptic
yellow_postsynaptic <- read.table(file = "Results/PPI_network/shortest_path_length/yellow_SPL_postsynaptic_postsynaptic.txt",
                                     sep = "\t",
                                     header = T,
                                     row.names = 1)
yellow_postsynaptic_list <- as.vector(as.matrix(yellow_postsynaptic[lower.tri(yellow_postsynaptic)]))
yellow_postsynaptic_others <- read.table(file = "Results/PPI_network/shortest_path_length/yellow_SPL_postsynaptic_others.txt",
                                            sep = "\t",
                                            header = T,
                                            row.names = 1)
yellow_postsynaptic_others_list <- as.vector(as.matrix(yellow_postsynaptic_others))
yellow_postsynaptic_res <- wilcox.test(yellow_postsynaptic_list, yellow_postsynaptic_others_list, alternative = "less")
yellow_postsynaptic_res$mean1 <- mean(yellow_postsynaptic_list)
yellow_postsynaptic_res$mean2 <- mean(yellow_postsynaptic_others_list)
yellow_postsynaptic_res <- unlist(yellow_postsynaptic_res)
#significant

#yellow module neurexins
yellow_neurexins <- read.table(file = "Results/PPI_network/shortest_path_length/yellow_SPL_neurexins_neurexins.txt",
                                  sep = "\t",
                                  header = T,
                                  row.names = 1)
yellow_neurexins_list <- as.vector(as.matrix(yellow_neurexins[lower.tri(yellow_neurexins)]))
yellow_neurexins_others <- read.table(file = "Results/PPI_network/shortest_path_length/yellow_SPL_neurexins_others.txt",
                                         sep = "\t",
                                         header = T,
                                         row.names = 1)
yellow_neurexins_others_list <- as.vector(as.matrix(yellow_neurexins_others))
yellow_neurexins_res <- wilcox.test(yellow_neurexins_list, yellow_neurexins_others_list, alternative = "less")
yellow_neurexins_res$mean1 <- mean(yellow_neurexins_list)
yellow_neurexins_res$mean2 <- mean(yellow_neurexins_others_list)
yellow_neurexins_res <- unlist(yellow_neurexins_res)
#significant

#summarize
res <- rbind(brown_translation_res, blue_translation_res, blue_axon_res, blue_rho_res, turquoise_rho_res, turquoise_postsynaptic_res, yellow_postsynaptic_res, yellow_neurexins_res)
write.csv(res, "Results/PPI_network/shortest_path_length/res.csv")
