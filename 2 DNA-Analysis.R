#R test script
#yaoming yang
#2017-11-25
#working with DNA analyses-DNA Sequence and protein seq analyses

#import library
library(seqinr)
library(ape)
library(rentrez) 
library(Peptides) 

#1. load a FASTA (sigle or multiple DNA sequence)file: seqtype = "AA" for amino acid, default DNA
cox1<- read.fasta(file = "cox1.fasta", seqtype = "AA")

#how mangy squences in the file?
length(cox1) 

#get 1st sequence
seq1<- cox1[1] 

#2. retrieve DNA seq from GenNank
#install.packages("ape").
# import library(ape)

#read a DNA seq from GenBank
AB003468 <- read.GenBank("AB003468", as.character = "TRUE")

#write DNA seq into fasta format file
write.dna(AB003468, file ="AB003468.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = " ", colw = 10)

#3. retrieve DNA seq from GenNank with rentrez
# https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
#install.packages("rentrez") and import it
entrez_search(db="nucleotide", term="human superoxide dismutase")

#4. FASTA DNA Analysis
#strip the variable AB003468’s first (and only) sequence 
CloningVector <- AB003468[[1]]

#simple nucleotide count, "word count" analysis 
count <- count(CloningVector,1) #count single nucleotide
count

count <- count(CloningVector,2) #count dinucleotide combinations 
count

count <- count(CloningVector,3) #count trinucleotide  combinations 
count

#compute the GC content
GC <- GC(CloningVector)
GC

# GC Analysis
#compute GCwindow of the size 200
GCwindow <- seq(1, length(CloningVector)-200, by = 200)
GCwindow

#find the length of the GCwindow
n <- length(GCwindow)
Chunks <- numeric(n)

#For loop to compute GC per chunk
for (i in 1:n){
  chunk <- CloningVector[GCwindow[i]:(GCwindow[i]+199)]
  chunkGC <- GC(chunk)
  print(chunkGC)
  Chunks[i] <- chunkGC
}

Chunks

#plot GC Windows
plot(GCwindow,Chunks,type="b",xlab="Nucleotide start position",ylab="GC content")

#5. GC Custom-Window Plot Function
slidingwindowGCplot <- function(windowsize,inputseq){
  #compute GCwindow num with the windowsize
  GCwindow <- seq(1, length(inputseq)-windowsize, by = windowsize)
  #find the length of the GCwindow
  n <- length(GCwindow)
  #make a vector the same length but "blank" for us to fill
  Chunks <- numeric(n)
  #For loop to compute GC per chunk
  for (i in 1:n){
    chunk <- inputseq[GCwindow[i]:(GCwindow[i]+windowsize-1)]
    chunkGC <- GC(chunk)
    print(chunkGC)
    Chunks[i] <- chunkGC
  }
  #plot GC Windows, main=title, main=paste("title", variable), paste combine items
  plot(GCwindow,Chunks,type="b",xlab="Nucleotide start position",ylab="GC content",main=paste("GC Plot with windowsize", windowsize))
}

#invoke function, need to run function first
slidingwindowGCplot(200,CloningVector)

slidingwindowGCplot(175,CloningVector)

#6. protein sequence statistics
#install and import library(Peptides) 
# aa component of seq 1
aaComp(cox1[1])

# aa component of all 4 seq in cox1
aaComp(cox1)

# the aliphatic index of the protein sequence 
#(an indicator of thermostability of globular proteins) 
aIndex(cox1)

#predicts the net charge of the protein charge(seq, pH = 7, pKscale = "Lehninger")
charge(cox1) 

#predicts the net charge of a specified sequence 
charge(seq="FLPVLAG", pH=7, pKscale="EMBOSS")

#hydrophobicity: hydrophobicity(seq)
hydrophobicity(cox1[1]) 
