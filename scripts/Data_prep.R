# Data preparation
library(biomaRt)
library(org.Hs.eg.db)
library(stringr)
str_split_fixed(before$type, "_and_", 2)

SAAV <- read_excel("Data/Melanoma proteome/SAAVs identification.xlsx", sheet = "summary")
SAAV$UNIPROT = str_split_fixed(SAAV$Master.Protein.Accession, "-", 2)[,1]
annotation <- select(org.Hs.eg.db, keys=unique(SAAV$UNIPROT), columns=c('UNIPROT', 'SYMBOL', 'ENTREZID'), keytype="UNIPROT")
SAAV.out = merge(SAAV, annotation, by="UNIPROT")

write.csv(SAAV.out, 'Data/Melanoma proteome/SAAV.csv', row.names = FALSE)
