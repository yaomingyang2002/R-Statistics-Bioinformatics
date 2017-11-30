#yaoming yang
#2017-11-28
#linkage Analysis -Expression Quantitative Traits Loci, or eQTLs.

# packages and libraries
install.packages("MatrixEQTL")
library(MatrixEQTL)

#using the following 3 files for the  eQTL analysis 
# SNP.txt, GE.txt, Covariates.txt

#1. set up statistical model for eQTL analysis
useModel = modelLINEAR # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

#provide a significance threshold –  use 10-2: 
#only SNPs that have an affect on gene expression with a significance of 10-2 or greater will be retained
pvOutputThreshold = 1e-2

#choose a covariance matrix, setting all values to identity
#numeric() covariance is equivalent to identity
errorCovariance = numeric()

#load SNP, GE, Covariate data
snps = SlicedData$new()        # init a new and blank snp dataframe
snps$fileDelimiter = "\t"      #value are separated by tabs
snps$fileOmitCharacters = "NA" #ignore values that are NA 
snps$fileSkipRows = 1          #first row is label not data
snps$fileSkipColumns = 1       #first columen is a label
snps$fileSliceSize = 2000      #to read chunks of 2000 rows
snps$LoadFile("SNP.txt")       #the fileName

gene = SlicedData$new()
gene$fileDelimiter = "\t"      
gene$fileOmitCharacters = "NA" 
gene$fileSkipRows = 1          
gene$fileSkipColumns = 1       
gene$fileSliceSize = 2000      
gene$LoadFile("GE.txt")

cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      
cvrt$fileOmitCharacters = "NA" 
cvrt$fileSkipRows = 1          
cvrt$fileSkipColumns = 1       
cvrt$LoadFile("Covariates.txt")

#start analysis with outfile, note that we configured variables use the parameter name need to insert into the analysis 
me = Matrix_eQTL_engine(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt, 
  output_file_name = "outfile.txt", 
  pvOutputThreshold = pvOutputThreshold, 
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  pvalue.hist = FALSE, 
  min.pv.by.genesnp = FALSE, 
  noFDRsaveMemory = FALSE)

#output file is a list of SNP-gene combinations that exceeded our thresholds – indicating that those SNPs are somehow affecting the expression levels of the indicated genes. 
#For example – SNP_11 changes the expression levels of Gene_06

#have the report written to our console 
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n')
cat('Detected eQTLs:', '\n')
show(me$all$eqtls)

# other kind of output
write.csv(me$all$eqtls, file = "myeQTLResult.txt") # in a file
print(me$all$eqtls) #in the screen


#2. Determine if the location of a SNP (close, distant) has effect to our gene 
# download and add snpsloc.txt and geneloc.txt 

#add temporary files to store the cis (in same CH) vs trans (in diff CH) information
output_file_name_cis = tempfile()
output_file_name_tra = tempfile()

#set different significance thresholds for our cis vs trans relationships
pvOutputThreshold_cis = 1e-2
pvOutputThreshold_tra = 1e-2

#establish a maximium cis distance 
cisDist = 1e6

#read in the actual position files 
snpspos = read.table("snpsloc.txt", header = TRUE, stringsAsFactors = FALSE)
genepos = read.table("geneloc.txt", header = TRUE, stringsAsFactors = FALSE)

# run the analysis:
me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt, 
  output_file_name = output_file_name_tra, 
  pvOutputThreshold = pvOutputThreshold_tra, 
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  output_file_name_cis = output_file_name_cis,
  pvOutputThreshold_cis =pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE, 
  noFDRsaveMemory = FALSE)  
#noFDRsaveMemory = write the output directly to file
# min.pv.by.genesnp = save all p values even if they do not meet the threshold if set for TRUE