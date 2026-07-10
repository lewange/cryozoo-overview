library(readxl)
library(readr)


NCBI_Genome_Annotations <- read_excel("data/NCBI_Genome_Annotations.xlsx",
                                      sheet = "Annotations_without_links",.name_repair = "universal")

ENSEMBL_annotations <- read_csv("data/ENSEMBL_annotations.csv",name_repair = "universal")


full_table<-full_join(NCBI_Genome_Annotations,ENSEMBL_annotations,by=c("Species"="Scientific.name"))

saveRDS(full_table,"data/Ensembl_NCBI_annotations.rds")
