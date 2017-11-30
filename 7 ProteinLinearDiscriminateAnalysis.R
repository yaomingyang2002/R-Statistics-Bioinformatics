#yaoming yang
#2017-11-27
#Protein Struction,Expression, and Interaction
#Linear Discriminate Analysis, or LDA.
#LDA is similar to ANOVA and regression analysis, and attempts to find a linear combination of features that can be used to separate two or more classes of objects
#Here:3 different groups/conditions/177 separate samples/ mass spec -> protein quantity values.

#install & import libraries
library(timeSeries)
library(MASS)
library(rgl)
library(ggplot2)

#load ms.csv dataset 
Dataset <- read.csv("ms.csv", 
                    header=TRUE, na.strings="NA", dec=".", strip.white=TRUE)

#copy the protein data from the dataset
RawData <- Dataset[,2:14]

#eliminate "empty" columns by colSds where an empty column will have a 0 standard deviation. 
filledcols = colSds(RawData) != 0.0
RawData <- RawData[,filledcols]

#The lda command from the MASS libr do LDA analysis to creates the LDA probabilities 
#over the specified the class/condition/group ($X1 ~ .) and the data to test against
test1.lda <- lda(Dataset$X1 ~ . , data=Dataset)

#generate the actual LDA values by testing these probabilities against the data itself
test1.lda.values <- predict(test1.lda,Dataset)

#chose different columns data to plot
x <- test1.lda.values$x[,1]
y <- test1.lda.values$x[,2]

#generate a plot 
class <- Dataset$X1

#create a dataframe 
plotdata <- data.frame(class, x, y)

#calculate The center of a cluster, Centroid, use the aggregate mean in scatterplot
centroids <- aggregate(cbind(x,y)~class,plotdata,mean)

#calculate the distance for each centroid
CentroidDistances <- dist(centroids, method = "euclidean", diag = TRUE, upper = FALSE, p = 2)
CentroidDistances

#add the label value to CentroidDistances using  attr
attr(CentroidDistances, "Labels") <- centroids$class

#use ggplot2 to plot
plot1 <- ggplot(plotdata,aes(x,y,color=factor(class))) + geom_point(size=3)
plot1

#add in the centroids more geom_points
ggplot(plotdata,aes(x,y,color=factor(class))) + geom_point(size=3)+ geom_point(data=centroids,size=7)

#add labels to the points using geom_text
ggplot(plotdata,aes(x,y,color=factor(class))) + geom_point(size=3)+ geom_point(data=centroids,size=7) +geom_text(data=centroids, size=7, label=centroids$class, colour="black")

#add a title using the ggtitle command
ggplot(plotdata,aes(x,y,color=factor(class))) + geom_point(size=3) + geom_point(data=centroids,size=7) +geom_text(data=centroids, size=7, label=centroids$class, colour="black") + ggtitle("LDA of Conditions 1-3")

#label each data point with the condition 
ggplot(plotdata,aes(x,y,color=factor(class))) + geom_point(size=3) + geom_point(data=centroids,size=7) +geom_text(data=centroids, size=7, label=centroids$class, colour="black") + ggtitle("LDA of Conditions 1-3") +geom_text(aes(label=Dataset$X1),hjust=0, vjust=0, colour="black")

#output:
plot1 <- ggplot(plotdata,aes(x,y,color=factor(class))) + geom_point(size=3) + geom_point(data=centroids,size=7) +geom_text(data=centroids, size=7, label=centroids$class, colour="black") + ggtitle("LDA of Conditions 1-3") +geom_text(aes(label=Dataset$X1),hjust=0, vjust=0, colour="black")
plot1
print(plot1)
ggsave(filename="LDAplot1.pdf", plot = plot1)

#output the CentroidDistances as a csv 
#write.csv(as.matrix(CentroidDistances, file = "centroiddistances.csv"))
CentroidDistMatrix <- as.matrix(CentroidDistances) 
CentroidDistMatrix
write.csv(CentroidDistMatrix, file = "centroiddistances.csv")
