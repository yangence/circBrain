"%&%"=function(a,b)paste0(a,b)
library(MatrixEQTL)
options(stringsAsFactors = F)
chr=commandArgs(trailingOnly = T)[1]
CMCdir="/media/data3/CMC/data/"
basedir="/media/data3/circCMC/data/"
phenotype_file_name=basedir%&%"qtl_C/CIRI.RATIO.mat.zero.Res.Caucasian465.txt"
phenotypeinfo_file=basedir%&%"STA.pos"
genepos =read.table(phenotypeinfo_file,header=T)
pvOutputThreshold_cis = 0.05
pvOutputThreshold_tra = 0
cisDist = 1e5
covariates_file_name=0
#covariates_file_name=basedir%&%"qtl/covariateDx.txt"

useModel=modelLINEAR
## Load RNA phenotype Level data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(phenotype_file_name)

output_file_name_cis =basedir%&%"qtl_C/RatioZero/resultCis/cis.chr"%&%chr%&%".txt"
output_file_name_tra =""
SNP_file_name =CMCdir%&%"genotype/Gaucasian465/SNPmatrix/chr"%&%chr%&%".genotype.MAF5.mat"
snpspos= read.table(CMCdir%&%"genotype/Gaucasian465/SNPmatrix/chr"%&%chr%&%".genotype.MAF5.pos",header = T)

## Load covariates
#cvrt = SlicedData$new();
#cvrt$fileDelimiter = "\t";      # the TAB character
#cvrt$fileOmitCharacters = "NA"; # denote missing values;
#cvrt$fileSkipRows = 1;          # one row of column labels
#cvrt$fileSkipColumns = 1;       # one column of row labels
#if(length(covariates_file_name)>0) {
#  cvrt$LoadFile(covariates_file_name);
#  }
  # Load genotype data
  snps = SlicedData$new();
  snps$fileDelimiter = " ";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 0;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name)
  snps$columnNames=readLines("/media/data3/CMC/data/genotype/Gaucasian465/Gaucasian.id")


  ## Run the analysis

me = Matrix_eQTL_main(
snps = snps, 
gene = gene, 
#cvrt = cvrt,
output_file_name      = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel, 
#errorCovariance = errorCovariance, 
verbose = TRUE, 
output_file_name.cis  = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos, 
genepos = genepos,
cisDist = cisDist,
pvalue.hist = FALSE,
noFDRsaveMemory = FALSE,
min.pv.by.genesnp=TRUE)
write.table(me$cis$min.pv.gene,basedir%&%"qtl_C/RatioZero/minP/minpvgenesnp.chr"%&%chr%&%".txt",sep="\t",col.names = F,quote = F)


