#open Rstudio (right click to chose) as a 'system administrator' to install all packages to avoid some privilege issue

# install some bioinformatics add-ons and packages: 
install.packages("seqinr")

install.packages("ggplot2")

# install bioConductor
source("https://bioconductor.org/biocLite.R")
biocLite()

# Fetch latest bioConductor packages for Seq Alignment
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")

#for multi seq alignment need msa
source("https://bioconductor.org/biocLite.R")
biocLite("msa")


#lib for extract file
install.packages("reshape2")

#for retrieve DNA seq from Genebank
install.packages("ape")

#retrieve DNA seq from GenNank with rentrez
install.packages("rentrez")

# lib for protein analysis
install.packages("Peptides") 

#for phylogeny Maximum Parsimony analysis
install.packages("phangorn")


#for comparative genomics
# required for biomaRT
install.packages("stringi",type="win.binary")
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")

#for mass-spec analysis software
install.packages("gplots")


#for protein Linear Discriminate Analysis, or LDA.
install.packages("timeSeries")
install.packages("MASS")
install.packages("rgl")
install.packages("ggplot2")

#for protein-protein interaction analysis, PPI
#install & import libraries
source("http://bioconductor.org/biocLite.R")
biocLite("graph")
biocLite("Rgraphviz")
biocLite("RBGL")
biocLite("yeastExpData")

# packages and libraries for Linkage Analysis - GenomeWideAssociationStudies GWAS
install.packages("GenABEL")

#for linkage Analysis -Expression Quantitative Traits Loci, or eQTLs
install.packages("MatrixEQTL")

# For Proteomics- Microarray analysis
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("affy")
biocLite("limma")
library(affy)
library(limma)

# for RNASeq Analysis

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("edgeR")
biocLite("limma")
biocLite("Glimma")
biocLite("org.Mm.eg.db")

install.packages("gplots")
install.packages("RColorBrewer")