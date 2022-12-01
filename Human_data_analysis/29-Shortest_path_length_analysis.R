library("tidyverse")

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
