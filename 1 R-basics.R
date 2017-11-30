#R test script
#yaoming yang
#2017-11-25
#This project is try to demonstrate how to use R statistics in bioinformatics

# install some bioinformatics add-ons and packages: 
install.packages("seqinr")

install.packages("ggplot2")

# install bioConductor
source("https://bioconductor.org/biocLite.R")
biocLite()

#1. R basics
count <-0

#vector
primes <- c(1, 3, 5, 7, 11)
Truth <- c(TRUE, FALSE, TRUE)
Names <- c("Bob", "Ted", "Carol", "Alice") 

count <- 1
primes[2]

#vector and dataframe->table/data
organism <- c("Human", "Chimpanzee", "Yeast")
chromosomes <- c(23, 24, 16)
multicellular <- c(TRUE, TRUE, FALSE)
OrganismTable <- data.frame(organism, chromosomes, multicellular)

#select item in table
OrganismTable$organism 
OrganismTable$organism[2] 

#write and read table/data
write.table(OrganismTable, file = "MyData.csv",row.names=FALSE,
            na="", col.names=FALSE, sep=",")

NewDataFrame <- read.csv("MyData.csv", header=FALSE, sep=",")
NewDataFrame1 <- read.csv("MyData.csv", header=TRUE, sep=",")

#for-loop
count <- 0
for (val in OrganismTable$chromosomes){
  if(val >20)
    count = count+1
}
print (count)

count <- 0
for (i in 1:100){
  if(i >20)
    count = count+1
}
print (count)

#2. simple bar-plotting:
barplot(OrganismTable$chromosomes)

#3. SimplePlot:
#import lib
library(ggplot2)
library(reshape2) 

#read csv data
rawdata <- read.csv("plotdata.csv", header=TRUE)

#A Simple Scatter Plot: aes=aesthetics
ggplot(rawdata, aes(x=Subject, y=a)) + geom_point()

#for a better scatter plot: install  in console and import 
#install.packages("reshape2") 

#Melting the Data

melted = melt(rawdata, id.vars ="Subject", measure.vars = c("a","c","d","e","f","g","j","k"))

ggplot(rawdata, aes(x=Subject, y=a)) + geom_point()
myPlot <- ggplot(melted, aes(x=variable, y=value, col=Subject, group=Subject)) + geom_point() + geom_line()+
  xlab("Samlpe")+
  ylab("#Observed")+
  ggtitle("Some Observations I made in the lab")

#save the file
ggsave(filename="MyPlot.pdf", plot=myPlot)
