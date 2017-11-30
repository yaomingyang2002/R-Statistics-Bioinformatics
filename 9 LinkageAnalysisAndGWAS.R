#yaoming yang
#2017-11-28
#linkage Analysis & Genome-wide Association Studies (GWAS)

# packages and libraries
install.packages("GenABEL")
library(GenABEL)

#load both the genotype and the phenotype files
data <- load.gwaa.data(phenofile="phenotype.dat", genofile="genotype.raw")

#browse data: genotype and phenotype
#gtdata(data)
#phdata(data)

#access the genotype at markers 3 through 5 for the first individual in the data 
gtdata(data)[1,3:5]

# see chromosome that markers 3 through 5 
chromosome(gtdata(data)[,3:5])

# want the name of marker (SNP) 17000
snpnames(gtdata(data)[ ,17000])

# see the genotypes at marker 1177 for the first four individuals
as.character(gtdata(data)[1:4, 1177])

#see the reference ("wild type") allele 
refallele(gtdata(data)[ ,1177])

#look at the summary
summary.snp.data(gtdata(data))

#look at the distribution of the continuous trait (ct)
ct <- phdata(data)$ct

#plot the distribution of the continuous variable
hist(ct, col="slateblue", main="Distribution of ct")

# see if the continuous trait depends on sex
boxplot(ct ~ phdata(data)$sex, col=c("blue", "red"), xlab="sex", names=c("F", "M"), main="ct by sex")

# see the full description and instructions  check.marker 
?check.marker 

#do some quality control (QC) on the dataset with check.marker 
#call:95% interval, perid: an individual called marker set to 0.95, maf: Minor allele frenquency, p.lev: p value for HWE set 10e-08 
qc <- check.marker(data, call = 0.95, perid.call = 0.95, maf = 1e-08, p.lev = 1e-08) 

#create a new dataset containing only individuals and snps that are "ok"
data.qc <- data[qc$idok, qc$snpok]

#see if the continuous variable (ct) is associated with any of the markers in our dataset.
#using the Cochran-Armitage Trend Test with the qtscore command.
an <- qtscore(ct~1, data=data.qc, trait="gaussian")

#plot the result: cex command scaled the font of the labels down by 50%.
plot(an, col=c("olivedrab","slateblue"), cex=.5, main="Manhattan plot")

#make sure that our p-values for each of the many tests that occurred are really "significant"
#by looking at lambda, which is the genomic inflation factor, and can result in lower than expected p values due to cryptic relatedness or population structure.
#estimate is 1.18 – since this is above 1, we do have a small amount of genomic inflation
estlambda(an[,"P1df"], plot=T)

#correct it, QC data also contains genomic control information, stored in Pc1df
#the resulting lambda is much closer to one, so we should plot that value instead of P1df. 
estlambda(an[,"Pc1df"], plot=T) 

# plot the result with df ="Pc1df"
plot(an, col=c("olivedrab","slateblue"), cex=.5, main="Manhattan plot", df="Pc1df")

#Bonferroni correction
pval.threshold <- 0.05
bonferroni <- -log10(pval.threshold / nids(data.qc))

#re-plot the corrected plot
plot(an, col=c("olivedrab","slateblue"), cex=.5, main="Manhattan plot", df="Pc1df")
abline(h=bonferroni, lty=3, color="red")
