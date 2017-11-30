#yaoming yang
#2017-11-27
#Protein Struction,Expression, and Interaction
#Mass Spec Analysis

library(gplots)

#1.distribution plot of peptides in the peptidefrags.txt file. 

#reads in the file 
peptides.txt <- read.table("peptidefrags.txt", header=FALSE)

#convert values into a vector peptides 
peptides <-as.vector(peptides.txt$V1)

#generates a histogram with bins of size 400 or 200
hist(peptides,breaks=400) 
hist(peptides,breaks=200) 

#2 compare three different versions of mass-spec analysis software
#2.1 
#install and import install.packages("gplots")
#load and convert protein ID text files into vector
mascot.txt <- read.table("mascot.txt", header=FALSE)
mascot <-as.vector(mascot.txt$V1)

xtandem.txt <- read.table("xtandem.txt", header=FALSE)
xtandem <-as.vector(xtandem.txt$V1)

protpro.txt <- read.table("protpro.txt", header=FALSE)
protpro <-as.vector(protpro.txt$V1)

#create a list
combinedMSdata <- list(Mascot=mascot, XTandem=xtandem, ProtPro=protpro)

#create a Venn diagram using command from the gplots library
venn(combinedMSdata)
