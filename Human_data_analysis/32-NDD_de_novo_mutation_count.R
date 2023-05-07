library("tidyverse")
library("clusterProfiler")
library("limma")

ssc <- read_tsv(file = "NDD_mutations/denovo_db/denovo-db.ssc-samples.variants.tsv", 
                col_names = TRUE, 
                comment = "##",
                col_types = cols(Chr = col_character())
                )
nssc <- read_tsv(file = "NDD_mutations/denovo_db/denovo-db.non-ssc-samples.variants.tsv", 
                 col_names = TRUE, 
                 comment = "##",
                 col_types = cols(Chr = col_character())
)

denovo_db <- rbind(ssc, nssc)
colnames(denovo_db)[1] <- "SampleID"

unique(denovo_db$PrimaryPhenotype)
denovo_db2 <- denovo_db %>%
  mutate(PrimaryPhenotype = case_when(
    PrimaryPhenotype == "autism" ~ "ASD",
    PrimaryPhenotype == "intellectualDisability" ~ "ID",
    PrimaryPhenotype == "developmentalDisorder" ~ "DD",
    TRUE ~ PrimaryPhenotype)
  )
#LGD mutations include "stop-gained", "start-lost", stop-gained-near-splice", "frameshift","frameshift-near-splice", "splice-donor", "splice-acceptor"
denovo_db_filt <- denovo_db2 %>%
  filter(!is.na(Gene)) %>%
  filter(PrimaryPhenotype %in% c("control","ASD","ID", "DD")) %>%
  filter(FunctionClass %in% c("stop-gained", "start-lost", "stop-gained-near-splice", "frameshift","frameshift-near-splice", "splice-donor", "splice-acceptor", "missense", "missense-near-splice"))

by_gene_individual_count <- denovo_db_filt %>%
  mutate(MutationType = case_when(
    FunctionClass %in% c("stop-gained","start-lost","stop-gained-near-splice", "frameshift","frameshift-near-splice", "splice-donor", "splice-acceptor") ~ "LGD",
    FunctionClass %in% c("missense","missense-near-splice") ~ "missense")
  ) %>%
  group_by(PrimaryPhenotype, MutationType, SampleID) %>%
  dplyr::count(Gene) %>%
  ungroup()

by_gene_count <- by_gene_individual_count %>%
  group_by(PrimaryPhenotype, MutationType) %>%
  dplyr::count(Gene) %>%
  ungroup()

#SCZ
SCZ <- read.csv(file = "NDD_mutations/Schizophrenia/SCZ.csv", header = TRUE)
unique(SCZ$FunctionClass)
SCZ_filt <- SCZ %>%
  filter(!is.na(Gene)) %>%
  filter(FunctionClass %in% c("ptv", "missense", "Frameshift", "stop_gained","start_lost","Missense_SNV", "Splice_SNV"))
SCZ_by_gene_count <- SCZ_filt %>% 
  mutate(MutationType = case_when(
    FunctionClass %in% c("ptv", "Frameshift", "stop_gained", "start_lost", "Splice_SNV") ~ "LGD",
    FunctionClass %in% c("missense", "Missense_SNV") ~ "missense",
    TRUE ~ "error")
  ) %>%
  group_by(PrimaryPhenotype, MutationType) %>%
  dplyr::count(Gene) %>%
  ungroup()
unique(SCZ_by_gene_count$MutationType)

#epilepsy
epilepsy <- read.csv(file = "NDD_mutations/Heyne 2018 NDD with epilepsy/epilepsy.csv", header = TRUE)
unique(epilepsy$FunctionClass)
epilepsy_filt <- epilepsy %>%
  filter(!is.na(Gene)) %>%
  filter(FunctionClass %in% c("lof", "mis"))
epilepsy_by_gene_count <- epilepsy_filt %>% 
  mutate(MutationType = case_when(
    FunctionClass == "lof" ~ "LGD",
    FunctionClass == "mis" ~ "missense",
    TRUE ~ "error")
  ) %>%
  group_by(PrimaryPhenotype, MutationType) %>%
  dplyr::count(Gene) %>%
  ungroup()

#merge all data
NDD_variants <- rbind(by_gene_count, SCZ_by_gene_count, epilepsy_by_gene_count)
#change alias to official gene symbol
SYMBOL <- limma::alias2SymbolTable(NDD_variants$Gene, species = "Hs")
NDD_variants$SYMBOL <- SYMBOL
NDD_variants <- drop_na(NDD_variants, SYMBOL)
NDD_variants <- NDD_variants[,-3]
ids <- bitr(NDD_variants$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
colnames(ids) <- c("SYMBOL","ENTREZID")
NDD_variants <- dplyr::left_join(NDD_variants, ids)
NDD_variants <- NDD_variants[,c(1,2,3,5,4)]

#save results
write.csv(NDD_variants, file = "Results/NDD_mutations/NDD_de_novo_variants_number_in_each_gene.csv" , row.names = FALSE)
#note that NDD_variants have some duplication due to alias of gene symbol. They are corrected and saved in
#"Results/NDD_mutations/NDD_de_novo_variants_number_in_each_gene_corrected.csv"
#this does not affect the results of FET
