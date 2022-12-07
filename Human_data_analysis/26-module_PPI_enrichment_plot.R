library("tidyverse")

## This code is to calculate the number of interactions in the random blue module PPI networks.

## Calculate the degree of the background PPI network.

library(igraph)

data_network <- read.table("PPIN_background.txt", sep = "\t", colClasses = "character")

data_network_1 <- graph.data.frame(data_network, directed = FALSE)

data_degree <- degree(data_network_1)

write.table(data_degree, file = "degree_background.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)

## Generate the degree distribution of the background PPI network for randomization.

data_degree <- read.table("degree_background.txt", sep = "\t")

data_degree[data_degree[,2]>=200,2] <- 200

data_degree[(data_degree[,2]>=100)&(data_degree[,2]<200),2] <- 100

data_degree[(data_degree[,2]>=70)&(data_degree[,2]<100),2] <- 70

data_degree[(data_degree[,2]>=50)&(data_degree[,2]<70),2] <- 50

data_degree[(data_degree[,2]>=40)&(data_degree[,2]<50),2] <- 40

data_degree[(data_degree[,2]>=30)&(data_degree[,2]<40),2] <- 30

data_degree[(data_degree[,2]>=20)&(data_degree[,2]<30),2] <- 20

data_degree[(data_degree[,2]>=10)&(data_degree[,2]<20),2] <- 10

data_degree[data_degree[,2]<10,2] <- 1

write.table(data_degree, file = "degree_background_1.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

## Calculate the degree of the blue module proteins in the background PPI network.

data_degree <- read.table("degree_background_1.txt", sep = "\t", colClasses = "character")

data_protein <- read.table("gene_blue.txt", sep = "\t", colClasses = "character")

data_degree_1 <- merge(data_degree, data_protein, by = 1)

write.table(data_degree_1[,2], file = "degree.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

## Generate the degree distribution of the blue module proteins for randomization.

data_degree <- scan("degree.txt", sep = "\t")

data_degree_1 <- unique(data_degree)

data_degree_1 <- sort(data_degree_1)

data_degree_2 <- numeric(length = length(data_degree_1))

for (i in 1:length(data_degree_1))
{
  data_degree_3 <- data_degree[data_degree==data_degree_1[i]]
  
  data_degree_2[i] <- length(data_degree_3)
}

data_degree_4 <- data.frame(data_degree_1, data_degree_2)

write.table(data_degree_4, file = "degree_1.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

## Calculate the random numbers of interactions of the blue module PPI network.

data_network <- read.table("PPIN_background.txt", sep = "\t", colClasses = "character")

data_degree <- read.table("degree_background_1.txt", sep = "\t", colClasses = "character")

data_degree_1 <- read.table("degree_1.txt", sep = "\t", colClasses = "character")

data_edge <- numeric(length = 10000)

set.seed(100000000)

for (i in 1:100000)
{
  for (j in 1:nrow(data_degree_1))
  {
    data_degree_2 <- data_degree[data_degree[,2]==data_degree_1[j,1],]
    
    data_protein <- sample(data_degree_2[,1], data_degree_1[j,2])
    
    write.table(data_protein, file = "protein_random.txt", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  
  data_protein_1 <- read.table("protein_random.txt", sep = "\t", colClasses = "character")
  
  data_network_1 <- merge(data_network, data_protein_1, by = 1)
  
  data_network_2 <- data_network_1[,2:1]
  
  data_network_3 <- merge(data_network_2, data_protein_1, by = 1)
  
  data_edge[i] <- nrow(data_network_3)
  
  file.remove("protein_random.txt")
}

write.table(data_edge, file = "edge_number_random.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

## Calculate the observed number of interactions of the blue module PPI network.

data_network <- read.table("PPIN_background.txt", sep = "\t", colClasses = "character")

data_protein <- read.table("gene_blue.txt", sep = "\t", colClasses = "character")

data_network_1 <- merge(data_network, data_protein, by = 1)

data_network_2 <- data_network_1[,2:1]

data_network_3 <- merge(data_network_2, data_protein, by = 1)

data_edge <- nrow(data_network_3)

write.table(data_edge, file = "edge_number_blue.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)




##########################################################
#brown module
brown_random <- read.table(file = "Results/PPI_network/Module_PPI_edge_number_enrichment_permutation/brown_edge_number_random.txt",
                         header = F)
colnames(brown_random) <- "Edge_number"
module_edge <- 1270
p <- (sum(brown_random$Edge_number > module_edge)+1)/100001
ggplot(brown_random, aes(x=Edge_number)) +
  geom_density(fill="lightblue") +
  annotate("segment", x = module_edge, xend = module_edge, y = 0, yend = 0.006, colour = "red") +
  xlim(c(800,1350)) +
  annotate("text", x = module_edge, y = 0.0065, label = expression("p<1.00 x"~10^-5)) +
  labs(x= "Edge number", y= "Density") +
  theme_bw()
ggsave("Results/PPI_network/results/enrichment/brown_PPI_enrichment.pdf", width = 6, height = 4, units = "in")

#blue module
blue_random <- read.table(file = "Results/PPI_network/Module_PPI_edge_number_enrichment_permutation/blue_edge_number_random.txt",
                          header = F)
colnames(blue_random) <- "Edge_number"
module_edge <- 864
p <- (sum(blue_random$Edge_number > module_edge)+1)/100001
ggplot(blue_random, aes(x=Edge_number)) +
  geom_density(fill="lightblue") +
  annotate("segment", x = module_edge, xend = module_edge, y = 0, yend = 0.008, colour = "red") +
  xlim(c(580,900)) +
  annotate("text", x = module_edge, y = 0.0085, label = expression("p<1.00 x"~10^-5)) +
  labs(x= "Edge number", y= "Density") +
  theme_bw()
ggsave("Results/PPI_network/results/enrichment/blue_PPI_enrichment.pdf", width = 6, height = 4, units = "in")

#turquoise module
turquoise_random <- read.table(file = "Results/PPI_network/Module_PPI_edge_number_enrichment_permutation/turquoise_edge_number_random.txt",
                          header = F)
colnames(turquoise_random) <- "Edge_number"
module_edge <- 1083
p <- (sum(turquoise_random$Edge_number > module_edge)+1)/100001
ggplot(turquoise_random, aes(x=Edge_number)) +
  geom_density(fill="lightblue") +
  annotate("segment", x = module_edge, xend = module_edge, y = 0, yend = 0.007, colour = "red") +
  xlim(c(650,1200)) +
  annotate("text", x = module_edge, y = 0.0075, label = expression("p<1.00 x"~10^-5)) +
  labs(x= "Edge number", y= "Density") +
  theme_bw()
ggsave("Results/PPI_network/results/enrichment/turquoise_PPI_enrichment.pdf", width = 6, height = 4, units = "in")

#yellow module
yellow_random <- read.table(file = "Results/PPI_network/Module_PPI_edge_number_enrichment_permutation/yellow_edge_number_random.txt",
                               header = F)
colnames(yellow_random) <- "Edge_number"
module_edge <- 464
p <- (sum(yellow_random$Edge_number > module_edge)+1)/100001
ggplot(yellow_random, aes(x=Edge_number)) +
  geom_density(fill="lightblue") +
  annotate("segment", x = module_edge, xend = module_edge, y = 0, yend = 0.013, colour = "red") +
  xlim(c(150,500)) +
  annotate("text", x = module_edge, y = 0.0135, label = expression("p<1.00 x"~10^-5)) +
  labs(x= "Edge number", y= "Density") +
  theme_bw()
ggsave("Results/PPI_network/results/enrichment/yellow_PPI_enrichment.pdf", width = 6, height = 4, units = "in")

#grey module
grey_random <- read.table(file = "Results/PPI_network/Module_PPI_edge_number_enrichment_permutation/grey_edge_number_random.txt",
                            header = F)
colnames(grey_random) <- "Edge_number"
module_edge <- 2614
p <- (sum(grey_random$Edge_number > module_edge)+1)/100001
ggplot(grey_random, aes(x=Edge_number)) +
  geom_density(fill="lightblue") +
  annotate("segment", x = module_edge, xend = module_edge, y = 0, yend = 0.007, colour = "red") +
  xlim(c(2250,2850)) +
  annotate("text", x = module_edge, y = 0.0075, label = expression("p=0.15640")) +
  labs(x= "Edge number", y= "Density") +
  theme_bw()
ggsave("Results/PPI_network/results/enrichment/grey_PPI_enrichment.pdf", width = 6, height = 4, units = "in")

