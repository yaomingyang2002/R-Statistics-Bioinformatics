#yaoming yang
#2017-11-27
#Genomic Analysis: comparative Genomics
#compara the gene contents between two species.
#Biomart lets you access ENSEMBLE genomes directly and download relevant information 

#libraries
# required for biomaRT
install.packages("stringi",type="win.binary")

source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library("stringi", lib.loc="~/R/win-library/3.4") # this is required for import biomaRT!
library(biomaRt)
library(Biostrings)
library(seqinr)

#see the list of available databases in ENSEMBLE genomes
listMarts()

#Pick up the ensemble database and see available datasets 
ensembl=useMart("ensembl")

listDatasets(ensembl) 

#retrieve and load the data for the Chimpanzee and the Human datasets.
Chimpanzee <- useDataset("ptroglodytes_gene_ensembl", mart = ensembl)
Human <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

#pull the list of available attributes down for each dataset 
ChimpanzeeAttributes <- listAttributes(Chimpanzee)

HumanAttributes <- listAttributes(Human)
HumanAttributes

#use getBM command to retrieve the data from the ensembl servers, using the mart ID
ChimpanzeeGenes <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"), mart=Chimpanzee)
HumanGenes <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"), mart=Human)

#parse the gene names and gene types into two separate vectors for each species
ChimpanzeeGeneNames <- ChimpanzeeGenes[[1]]
ChimpanzeeGeneTypes <- ChimpanzeeGenes[[2]]

HumanGeneNames <- HumanGenes[[1]]
HumanGeneTypes <- HumanGenes[[2]]

#summarize the Gene Types using the TABLE command
ChimpanzeeTypesTable <- table(ChimpanzeeGeneTypes)
HumanTypesTable <- table(HumanGeneTypes)

ChimpanzeeTypesTable
HumanTypesTable

#compare the number of protein coding genes in humans vs chimpanzees
ChimpanzeeProteins <- ChimpanzeeTypesTable["protein_coding"]
HumanProteins <- HumanTypesTable["protein_coding"]

#can compare other attributes in the ensemble Attributes list

#identify human orthologs in Chimpanzees: "ens_hs_gene" 
ChimpHSNum <- getBM(attributes = "ens_hs_gene", mart = Chimpanzee)

#retrieves the list of both Chimpanzee gene ID along with their human orthologs (if available) in a dataframe 
ChimpanzeeHS <- getBM(attributes = c("ensembl_gene_id","ens_hs_gene"), mart = Chimpanzee)

#do the reverse now using the Human mart 
HumanCNum <- getBM(attributes = "ptroglodytes_homolog_ensembl_gene", mart = Human)
HumanC <- getBM(attributes = c("ensembl_gene_id","ptroglodytes_homolog_ensembl_gene"), mart = Human)

#BioMart introduction:https://www.youtube.com/watch?v=QvGT2G0-hYA
#http://www.ensembl.org/index.html