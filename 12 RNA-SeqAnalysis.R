#yaoming yang
#2017-11-28
#RNASeq: 1 RNAseq Data Processing 2. RNAseq Quality Control 3. Analysis -MDS Dimensionality

#install and import Bioconductor libraries:

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("edgeR")
biocLite("limma")
biocLite("Glimma")
biocLite("org.Mm.eg.db")

install.packages("gplots")
install.packages("RColorBrewer")

library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)

library(gplots)
library(RColorBrewer)

#1 RNAseq Data Processing:
#1.1 download and place GSE60450_Lactation-GenewiseCounts.txt, SampleInfo.txt  in working directory
#1.2 read the RNAseq count data: this data has already been processed on the instrument used to capture the data
seqdata <- read.delim("GSE60450_Lactation-GenewiseCounts.txt", stringsAsFactors = FALSE)
sampleinfo <- read.delim("SampleInfo.txt")

#1.3 remove the first two columns from seqdata 
countdata <- seqdata[,-(1:2)]

#1.4 add the Entrez Gene Id's (the first column of seqdata) as the row names for countdata.
rownames(countdata) <- seqdata[,1]

#1.5 extract just a portion of the characters from each of the column names in seqdata 
# start the substr at 1 and end it at 7.
colnames(countdata) <- substr(colnames(countdata),start=1,stop=7)

#1.6 filter the data to eliminate genes expressed at very low levels.
# using Counts Per Million(CPM) values
#start by extracting the CPMs and examining them.
myCPM <- cpm(countdata)

#1.7 filter out anything with a CPM of less than 0.5 which is our threshold.
#creating a dataframe called thresh 
thresh <- myCPM > 0.5

#1.8 see how many genes passed through our filter using the the rowSums command
table(rowSums(thresh))

#1.9 filter out the lowest values and show TRUE in zero or 1 samples
keep <- rowSums(thresh) >= 2
counts.keep <- countdata[keep,]

#1.10 see if our data have the correct number of CPM's by plotting them.
plot(myCPM[,1],countdata[,1])

#1.11 plot it again, truncating the values of the Y and X axis to display:
#only see values under 3.0 on the X axis and under 50 on the Y axis.  
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))

#1.12 add a line at the 0.5 on the X axis ->set of 0.5 CPM using abline(v=0.5).
abline(v=0.5)
       
#2 RNASeq Quality Control
#2.1 convert the truncated data from counts.keep to the dgelist format required for the EdgeR library
y <- DGEList(counts.keep)

#2.2 Quality Control (QC) checks on data:
#2.2.1 barplot to see if there are any discrepancies between the samples 
par(mfrow=c(1,1)) # back to single panel
barplot(y$samples$lib.size,names=colnames(y),las=2) #las=2 command rotates axis text
title("Barplot of library sizes") 

#2.2.2 boxplot of the log of the counts
#generate the log of the CPM 
logcounts <- cpm(y, log=TRUE) #get the log base10 of the values
#boxplot of the logcounts variable 
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2, cex = 0.5)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

#2.2.3 Multidimensional Scaling plot(MDSplot):
#side-by-side MDS plots of the data where we color the samples by cell type on one side and by condition on the other
#set our plot area for dual-plot mode
par(mfrow=c(1,2))

#2.2.3.1 MDS plots of the data color the samples by cell type
#find out how many cell types we have in our samples 
levels(sampleinfo$CellType) #two cell types in our data

#define a color palette to help us distinguish these two cells 
col.cell <- c("red","blue")[sampleinfo$CellType]

#generate an MDS plot by invoking plotMDS and using col.cell as the color element. 
plotMDS(y,col=col.cell)
legend("topleft",fill=c("blue","red"),legend=levels(sampleinfo$CellType))
title("Cell type")

#2.2.3.2 MDS plots of the data color the samples by condition
#find out how many by condition (or status) we have in our samples 
levels(sampleinfo$Status)  #three conditions in our data

#define a color palette to help us distinguish these three conditions 
col.status <- c("blue","red","dark green")[sampleinfo$Statuse]

#generate an MDS plot by invoking plotMDS and using col.status as the color element. 
plotMDS(y, col=col.status)
legend("topleft",fill=c("blue","red","dark green"), legend=levels(sampleinfo$Status),
       cex=0.8)
title("Status")

#2.2.4 using corrected data re-do plot MDS
# take a look at the corrected data SampleInfo_Corrected.txt. 
#load the corrected data
sampleinfo <- read.delim("SampleInfo_Corrected.txt")

#re-do plot MDS
par(mfrow=c(1,2))

col.cell <- c("blue","red")[sampleinfo$CellType]
col.status <- c("blue","red","dark green")[sampleinfo$Status]

plotMDS(y,col=col.cell)
legend("topleft",fill=c("blue","red"),legend=levels(sampleinfo$CellType))
title("Cell type")

plotMDS(y,col=col.status)
legend("topleft",fill=c("blue","red","dark green"), legend = levels(sampleinfo$Status), cex=0.8)
title("Status")

#3. RNAseq Analysis -MDS Dimensionality
#MDS plots for alternative (other) dimensionality
#set it back to single plot mode 
par(mfrow=c(1,1))

#3.1 generate an MDS plot using the default dimensions, but differentiating both the status and cell type in one graph
#generate a character set to plot for the cell type 
char.celltype <- c("X", "O") [sampleinfo$CellType]

#repeat plotMDS include the pch command to plot the cell types by the characters
plotMDS(y,col=col.status,pch=char.celltype,cex=2)
legend("topright",legend=levels(sampleinfo$Status),col=col.status,pch=16)
legend("bottomright",legend=levels(sampleinfo$CellType),pch=c(1,4))
title("MDS Plot for first two dimension")

#Analysis - MDS Dimensionality 
#add the dim command and plot alternate dimensions
par(mfrow=c(1,1))
plotMDS(y,dim=c(3,4), col=col.status,pch=char.celltype,cex=2)
legend("topright",legend=levels(sampleinfo$Status),col=col.status,pch=16)
legend("bottomright",legend=levels(sampleinfo$CellType),pch=c(1,4))
title("MDS Plot for dimensions 3 and 4") 

#dim command should come immediately after the Y 
#we've done with the MDS plots is really an example of cluster analysis 

#4. cluster analysis  – Hierarchical Clustering
#For large datasets the best way to visualize the results of clustering are heatmaps
#gplots contains a command, heatmap.2, that will both perform a hierarchical
# clustering and generate a heatmap.

#4.1 use heatmap.2 to generate a hierarchical clustering heatmap
#comput the variance in each row of our logcounts variable
var_genes <- apply(logcounts, 1, var)

#4.2 limit the data to just the 500 most variant genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]

#4.3 subset logcounts by choosing only the genes that show up in select_var 
highly_variable_lcpm <- logcounts[select_var,]

#4.4 see a summary of the contents of highly_variable_lcpm with the dim command 
dim(highly_variable_lcpm) #[1] 500  12

#4.5 choosing the color palette a heatmap with brewer.pal from the RColorBrewer library
#11 colors from the pre-defined set contained in RColorBrewer called RdYlBu
mypalette <- brewer.pal(11,"RdYlBu") #?brewer.pal 

#4.6 use the colorRampPalette command from the GRDevices
#library to generate a series of interpolated colors 
morecols <- colorRampPalette(mypalette)

#4.7 the clustering and generate the heatmap – again, heatmap.2 from gplots\
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most
variable genes across samples", ColSideColors = col.cell,scale="row")
#trace="none" command instructs heatmap.2 not to draw connecting lines through columns or rows 

#4.8 save our heatmap to a file:
#4.8.1 could do this directly through Rstudio, using the Export button found at the top of the PLOT

#4.8.2 a plot can be redirected to an image file (e.g. png or pdf
png(file="High_var_genes.heatmap.png", width = 600, height = 800 )
# follow by regenerate heatmap, the output instead went into your PNG :
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most
variable genes across samples",ColSideColors=col.cell,scale="row")
#turns OFF the redirect.  
dev.off()

#4.8.3. execute the dev.copy2pdf(fileName, width = px, height = px, horizontal = FALSE, etc) command after ANY plot (higher resolution )
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most
variable genes across samples", ColSideColors = col.cell,scale="row")

dev.copy2pdf(file="High_var_genes.heatmap.pdf")

