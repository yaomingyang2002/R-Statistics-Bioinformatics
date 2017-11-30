#R test script
#yaoming yang
#2017-11-26
#Sequence Alignment and Analysis

# Fetch latest bioConductor packages
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")

#libraries
library(Biostrings)
library(msa) 
library(seqinr)
library(phangorn) 


#1. DNA seq pairwise alignment
#read DNA file, with 6 seqs
prokaryotes <- read.fasta(file = "prok.fasta", seqtype = "DNA")

#split file and convert into two simple text string seq for pairwise alignment
seq1 <-as.character(prokaryotes[[1]])
seq1 = paste(seq1, collapse = "")

seq2 <-as.character(prokaryotes[[2]])
seq2 = paste(seq2, collapse = "")

#align seq1 and seq2 using the defaults settings of Biostrings
pairalign <- pairwiseAlignment(pattern = seq2, subject = seq1)

# summary statistics
pairalign
summary(pairalign) 

#convert the alignment to strings
pairalignString = BStringSet( c( toString( subject(pairalign) ), toString(pattern(pairalign))))

#export and write as a FASTA file
writeXStringSet(pairalignString, "aligned.txt", format="FASTA")

#2. Pair protein sequences: Dotplots
#read protein seq using seqinr's read.fasta 
coxgenes <- read.fasta(file = "cox1multi.fasta", seqtype="AA")

#parses the first two sequences out of this multi-sequence FASTA 
cox1 <- as.character(coxgenes[[1]])
cox2 <- as.character(coxgenes[[2]])

#generate the simplest dotplot:
dotPlot(cox1, cox2, main = "Human vs Mouse Cox1 Dotplot")

#Add window parameters:wsize=windowsize nmatch=num nt/bt matched<=wsize
dotPlot(cox1, cox2, wsize = 3, wstep = 3, nmatch = 3, main = "Human vs Mouse Cox1 Dotplot\nwsize = 3, wstep = 3, nmatch = 3")

#compare fragement dotplot, first 100 residues
dotPlot(cox1[1:100], cox2[1:100], wsize = 3, wstep = 3, nmatch = 3, main = "Human vs Mouse Cox1 first 100 AADotplot\nwsize = 3, wstep = 3, nmatch = 3")

#3. multiple seq alignment
# need to install and import library(msa)
#3.1 simple MSA. read fasta files as character sets use the Biostrings readAAStringSet and readDNAStringSet
coxAA <- readAAStringSet("cox1multi.fasta")
prokDNA <- readDNAStringSet("prok.fasta")

#aligns the cox1 protein sequences using CLUSTALW and the default alignment parameters.
coxAligned <- msa(coxAA)
coxAligned
#align the DNA sequences similarly
prokAligned <- msa(prokDNA)
prokAligned

#displays complete alignments:
print(prokAligned, show="complete")
print(coxAligned, show="complete")

#3.2 Multiple Sequence Alignments (MSA) Algorithms:
#using default substitution matrix, CLUSTAL 2.1
prokAligned <-msa(prokDNA, "ClustalW")  
prokAligned

#using Gonnet,ClustalOmega 1.2.0
prokAligned <-msa(prokDNA, "ClustalOmega")
prokAligned

#MUSCLE 3.8.31
prokAligned <-msa(prokDNA, "Muscle") 
prokAligned

#3.3 Multiple Sequence Alignments (MSA) parameters:cluster, gapOpening,
#gapExtension,Maxiters,  SubstitutionMatrix, Order(aligned, input), type(=dna, rna, protein), Verbose
#can be added in a comma-separated fashion to the msa command
#https://www.bioconductor.org/packages/devel/bioc/vignettes/msa/inst/doc/msa.pdf

coxAligned <- msa(coxAA, cluster="upgma")
coxAligned

# 3.4.1 msa output: convert to stringset ->save to fasta format
# can use SeaView: http://doua.prabi.fr/software/seaview
prokAlignStr = as(prokAligned, "DNAStringSet")
writeXStringSet(prokAlignStr, file="prokAligned.fasta")

coxAlignStr = as(coxAligned, "AAStringSet")
writeXStringSet(coxAlignStr, file="coxAligned.fasta")

#3.4.2 Exporting MSAs - PHYLIP Format, open with BBEdit, Sublime Text, or Seaview
#did't have to convert the alignment to a stringset
write.phylip(coxAligned, "coxAligned.phylip")
write.phylip(prokAligned, "prokAligned.phylip")

# 4. Phylogenetic Reconstruction: distance, maximum parsimony, and maximum likelihood 
#4.1 Distance Phylogeny
#convert prokaryotic alignment to seqinr format 
prokAligned2 <- msaConvert(prokAligned, type="seqinr::alignment")

#4.1.1 generate a distance matrix using seqinr with the ape package
#identity: relative to identical
prokdist <- dist.alignment(prokAligned2, "identity")
prokdist

#4.1.2 use this matrix to construct a basic distance phylogenetic tree using the neighbor-joining method
library(ape)
prokTree <- nj(prokdist)
plot(prokTree) 

#4.2 maximum parsimony
# need to install and import install.packages("phangorn")
#4.2.1 generate the phylogenetic data using phangorn's PhyDat and can also use msa  
prokAligned3 <- msaConvert(prokAligned, type="phangorn::phyDat")

#4.2.2 generates parsimony trees with pratchet
ParsTree <- pratchet(prokAligned3)

#4.2.3 plot the ParsTree
plot(ParsTree) 

#4.3 Maximum Likelihood
#4.3.1 calculate the likelihood using the pml function
fit <- pml(prokTree, prokAligned3)

#4.3.2 optimize this tree using optim.pml with Jukes-Cantor model K80, or the Kimura two-parameter model
fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")
fitK80 <- optim.pml(fit, model = "K80", rearrangement = "stochastic")

#4.3.3 plot the maximum likelihood tree 
plot(fitJC) 
plot(fitK80) 

#4.4 Bootstrapping (stastistics)
#using the Jukes-Cantor Maximum-Likelihood tree 
bootstrapped <- bootstrap.pml (fitJC, bs=100, optNni=TRUE, multicore=FALSE, control = pml.control(trace=0))

#plot and visualize it
plotBS(midpoint(fitJC$tree), bootstrapped, p = 50, type="p")

