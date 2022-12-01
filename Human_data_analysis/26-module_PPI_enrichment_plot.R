library("tidyverse")

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




