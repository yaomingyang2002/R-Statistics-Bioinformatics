#yaoming yang
#2017-11-27
#Protein Struction,Expression, and Interaction
#protein-protein interaction (PPI)

#install & import libraries
source("http://bioconductor.org/biocLite.R")
biocLite("graph")
biocLite("Rgraphviz")
biocLite("RBGL")
biocLite("yeastExpData")
library(graph)
library(Rgraphviz)
library(RBGL)
library(yeastExpData)

#use yeast expression data pre-loaded in BioConductor
#https://www.nature.com/articles/ng776
#The litG library contains interaction data represented as a graph, where the vertices, or nodes, represent proteins, and edges between nodes indicate interaction
data("litG")

#retrieve the names of all of the proteins(vertices/nodes) and store them in a vector.   
#The nodes command is from the graph library
litGnodes <- nodes(litG)
 
#query the 25th protein
litGnodes[25] 

# query the data for proteins that interact with yeast protein YFL039C
adj(litG, "YFL039C")

#extract the subgraphs, or connected components using the connectedComp routine in the RBGL library
connectedComponents <- connectedComp(litG)

#see all members of a subgraph 
connectedComponents[[1]]
connectedComponents[[3]]

#save all nodes of subgraph3 to a variable component3
component3 <- connectedComponents[[3]]

#use the subGraph command to pull out the nodes and edges from component3 
#subGraph checks this against the full dataset and gives us the full subgraph
subgraph3 <- subGraph(component3, litG)

#see the full stats on subgraph3
subgraph3

#use the Rgraphviz package to visualize the graph
#for built in Help:the command -> HELP tab-> the search window.
subgraph3plot <- layoutGraph(subgraph3, layoutType="neato")
renderGraph(subgraph3plot)

#look for the degree distribution for the PPI dataset using the degree command from the graph library
numdegrees <- graph::degree(litG)

#sort the values to plot the histogram
numdegreessorted <- sort(numdegrees)

#calculate the mean
meandeg <- mean(numdegreessorted)

#histo plot
hist(numdegreessorted, col="green", main = paste("Degree Distribution - Protein Interactions in litG with a mean of " ,meandeg))

hist(numdegrees, col="red", main = paste("Degree Distribution - Protein Interactions in litG with a mean of " ,meandeg))
