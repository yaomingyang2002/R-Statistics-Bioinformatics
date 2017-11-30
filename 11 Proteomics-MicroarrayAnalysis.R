#yaoming yang
#2017-11-28
#Proteomics: 1. Microarray analysis, 2 Differentially Expressed Genes 3. Visualization 

#1. Microarray analysis

#1.1 packages and libraries
source("http://bioconductor.org/biocLite.R")
biocLite()

biocLite("affy")
biocLite("limma")
library(affy)
library(limma)

#1.2 download, unzip and put files directly into working directory
#.CEL files are a standard output from an Affymetrix array experiment 
#Su_CELs1.zip, Su_CELs2.zip , Su_CELs3.zip , Su_CELs4.zip, Brain.fetalbrain.2color.dat.txt 

#load all of the *.CEL files 
affy.data <- ReadAffy()

#1.3 normalize and correct dataset
#This dataset was created using the mas5 probe set, need to normalize the data to compensate for differences between runs 
#Affymetrix provides standardization files for all of their probe sets, included in the affy library
eset.mas5 <- mas5(affy.data)

#to retrieve the normalized expression matrix and store it, can see the data in environment: exprSet.nologs 
exprSet.nologs <- exprs(eset.mas5)

#rename columns using the colnames command, then can see the change in data in environment: exprSet.nologs 
colnames(exprSet.nologs) <- c("brain.1", "brain.2", "fetal.brain.1", "fetal.brain.2", "fetal.liver.1",
                              "fetal.liver.2", "liver.1", "liver.2")
#1.4 determine differences in expression
#log transform of the data using the log command to express differences in expression in discrete single digits, e.g. +1/-1, to represent fold change in expression.
#use log base 2 transforms, we end up with a binary representation. then see diference between exprSet and exprSet.nologs
exprSet <- log(exprSet.nologs, 2)

#1.5 output to a tab-separated file, which can easily be imported into Excel (as can comma-separted, or .csv).
write.table(exprSet, file="Su_mas5_matrix.txt", quote=F, sep="\t")

#1.6 "Absent/Present" check for signals using the mas5calls command since we are using mas5 probes.
#generates a vector with the A/P values for the dataset
data.mas5calls <- mas5calls(affy.data)
#converts this into an expression matrix containing A or P for each tissue and gene combination
data.mas5calls.calls <- exprs(data.mas5calls)

#1.7 the processed image data: containing the color/intensity data for each spot on the array
#brain.fetalbrain.2color.data.tx
#now load that using the command from the limma package:
brain.fetalbrain.2color <- read.maimages("brain.fetalbrain.2color.data.txt", columns=list(G="brain.1",
                                                                                          R="fetal.brain.1", 
                                                                                          Gb="bg1", 
                                                                                          Rb="bg2"))
#normalize the color data using another command from the limma package 
brain.fetalbrain.2color.loess <- normalizeWithinArrays(brain.fetalbrain.2color, method="loess")

#plot a two-panel graph to see the difference for the dataset and the end-result  
par(mfrow=c(1,2))
plotMA(brain.fetalbrain.2color)
plotMA(brain.fetalbrain.2color.loess)

#1.8 calculate the log base 2 ratios of expression to compare those "binary" expression levels between conditions
#1.8.1 mean from exprSet: exprSet contains two columns for each tissue type, e.g. brain.1 and brain.2 the results of two separate chip runs
brain.mean <- apply(exprSet[, c("brain.1", "brain.2")], 1, mean)
fetal.brain.mean <- apply(exprSet[, c("fetal.brain.1", "fetal.brain.2")], 1, mean)
liver.mean <- apply(exprSet[, c("liver.1", "liver.2")], 1, mean)
fetal.liver.mean <- apply(exprSet[, c("fetal.liver.1", "fetal.liver.2")], 1, mean)

#1.8.2 calculate the ratios: log(A / B) = log(A) - log(B) for log-transformations of the expression values
brain.fetal.to.adult <- fetal.brain.mean - brain.mean
liver.fetal.to.adult <- fetal.liver.mean - liver.mean

#1.8.3 output
#make a new dataframe with all the data
all.data <- cbind(exprSet, brain.mean, fetal.brain.mean, liver.mean, fetal.liver.mean, 
                  brain.fetal.to.adult,
                  liver.fetal.to.adult)
# export all.data
write.table(all.data, file = "Microarray_ALL.txt", quote=F, sep="\t")

#2 Differentially Expressed Genes
#2.1 T test on our adult vs fetal brain data.
#create two new datasets containing the concantenated values for adult and fetal brain tissues.
dataset.1 <- exprSet[1, c("brain.1", "brain.2")]
dataset.2 <- exprSet[1, c("fetal.brain.1", "fetal.brain.2")]

#perform the T test: a two-sided, or two-tailed T test – meaning we are looking for outliers on both ends of the spectrum (the outer 5% of both ends of the normal distribution). 
t.test.gene.1 <- t.test(dataset.1, dataset.2, "two.sided")

#apply a T test across the relevant columns to compare adult brain expression values to fetal.
brain.p.value.all.genes <- apply(exprSet, 1, function(x){ t.test(x[1:2], x[3:4]) $p.value } )

#do exactly the same for the liver values 
liver.p.value.all.genes <- apply(exprSet, 1, function(x) { t.test(x[5:6], x[7:8]) $p.value } )

#A/P test from before and filter out genes that are uninformative.
#paste together all of the A/P data for all of the chip data into a single vector
AP <- apply(data.mas5calls.calls, 1, paste, collapse="")
AP

#create a variable that lists the genes that are informative 
genes.present = names(AP[AP != "AAAAAAAA"])

#see how many genes are left 
length(genes.present) 

#create a new dataset by taking the subset of exprSet that has data 
exprSet.present <- exprSet[genes.present,]

exprSet 
exprSet.present 

#FDR, or false discovery rate, correction 
#look at our p values from the t test – but we need to adjust for the FDR, or false
#start by isolating the p values from our t test for only genes that are informative, again using our genes.present value
brain.raw.pvals.present <- brain.p.value.all.genes[genes.present]
liver.raw.pvals.present <- liver.p.value.all.genes[genes.present]

#FDR correction is a simple call already built into the stats methods of R, using the p.adjust command.
brain.fdr.pvals.present <- p.adjust(brain.raw.pvals.present, method="fdr")
liver.fdr.pvals.present <- p.adjust(liver.raw.pvals.present, method="fdr")

#sort genes by p-values (a new dataframe) and look at them in increasing order of significance (meaning that there are significant differences between tissues).
brain.fdr.pvals.present.sorted <- brain.fdr.pvals.present[order(brain.fdr.pvals.present)]
liver.fdr.pvals.present.sorted <- liver.fdr.pvals.present[order(liver.fdr.pvals.present)]

#look at the first 10 in terms of lowest p-values which means highest significance.
brain.fdr.pvals.present.sorted[1:10]
liver.fdr.pvals.present.sorted[1:10]

#output files 
#create a new dataset that contains only the genes that meet that 0.01 p value criteria
brain.DE.probesets <- names(brain.raw.pvals.present[brain.raw.pvals.present < 0.01])
liver.DE.probesets <- names(liver.raw.pvals.present[liver.raw.pvals.present < 0.01])

#pull out the log base 2 ratios for the genes in those datasets
brain.DE.log2.ratios <- all.data[brain.DE.probesets, c("brain.fetal.to.adult", "liver.fetal.to.adult")]
liver.DE.log2.ratios <- all.data[liver.DE.probesets, c("brain.fetal.to.adult", "liver.fetal.to.adult")]

#use write.table and output the significant genes 
write.table(brain.DE.log2.ratios, "brain.DE.log2.ratios.txt", sep="\t", quote=F)
write.table(liver.DE.log2.ratios, "liver.DE.log2.ratios.txt", sep="\t", quote=F)

#3. Visualization
#3.1 "reset" plotting panel back to generate single plots.
# par(mfrow=c(1,2)), a plot area with one row and two columns for plots
# par(mfrow=c(3,3)), a 3 by 3 grid of plots 
par(mfrow=c(1,1))

#3.2 scatter plot
#3.2.1 create plot data to make plot call simpler
x.data <- all.data[, "brain.mean"]
y.data <- all.data[, "fetal.brain.mean"]

#3.2.2 generate the scatter plot
plot(x.data, 
     y.data, 
     main = "Log2 expression in fetal brain (n=2) vs adult brain (n=2)", 
     xlab="brain", 
     ylab="fetal brain", 
     col="blue", 
     cex=0.5)
#3.2.3 add lines with abline command based on slope (1) and intercept (0)
abline(0,1)

#3.3 generate a Ratio Intensity Plot, also known as an MA or Bland-Altman plot.
#MA Plots commonly used for genomic data to display the difference between measurements between two samples by transforming the data onto a log ratio (or M) scale and an average (or A) scale; used as X and Y values.
#3.3.1 calculate our average (A) and log ratio (M) values
brain = all.data[, "brain.mean"]
fetal.brain = all.data[, "fetal.brain.mean"]


#3.4 Volcano Plot: the p-values on the Y axis and the fold-difference on the X.
#plots significance (p value) versus fold-change on the y and x axes.
#3.4.1 generate a matrix that contains the raw and FDR corrected p-values for all of the "present" genes.
expression.plus.pvals = cbind(exprSet.present, brain.raw.pvals.present, brain.fdr.pvals.present,
                              liver.raw.pvals.present, liver.fdr.pvals.present)
#3.4.2 calculate the log base 2 ratios by subtracting the already-log-transformed values.
log2.ratios = expression.plus.pvals[, "brain.1"] -  expression.plus.pvals[, "fetal.brain.1"]
#3.4.3 adult brain expression values
p.values = expression.plus.pvals[, "brain.raw.pvals.present"]

#3.4.4 generate the Volcano plot
plot(log2.ratios, -log(p.values, 10) )

#3.5 put together all of the plotting:
par(mfrow=c(1,2))
#the simple log ratios vs p-values plot
plot(log2.ratios, p.values)
#the log transformed  plot
plot(log2.ratios, -log(p.values, 10) )
