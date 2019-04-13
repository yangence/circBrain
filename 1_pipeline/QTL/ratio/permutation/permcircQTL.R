"%&%"=function(a,b)paste0(a,b)
library(MatrixEQTL)
options(stringsAsFactors = F)
numArgs=as.numeric(commandArgs(trailingOnly = T)[1])
chr=numArgs %% 22 + 1
thread=numArgs %/% 22 +1
CMCdir="/data/users/lzl/CMC/data/"
basedir="/home/lzl/project/circCMC/script/sm_permQTL_ratio/"
phenotype_file_name=CMCdir%&%"CIRI.RATIO.IM.Res.Cau465.txt"
phenotypeinfo_file=CMCdir%&%"STA.pos"
genepos =read.table(phenotypeinfo_file,header=T)
phenotype.mat=read.table(phenotype_file_name,header=T,stringsAsFactors=F,row.names=1)
pvOutputThreshold_cis = 1e-50
pvOutputThreshold_tra = 0
cisDist = 1e5
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

output_file_name_cis =basedir%&%"qtl/permutation/temp/cis.chr"%&%chr%&%".txt"
output_file_name_tra =""
SNP_file_name =CMCdir%&%"genotype/SNPmatrix/chr"%&%chr%&%".genotype.MAF5.mat"
snpspos= read.table(CMCdir%&%"genotype/SNPmatrix/chr"%&%chr%&%".genotype.MAF5.pos",header = T)

## Load covariates
#cvrt = SlicedData$new();
#cvrt$fileDelimiter = "\t";      # the TAB character
#cvrt$fileOmitCharacters = "NA"; # denote missing values;
#cvrt$fileSkipRows = 1;          # one row of column labels
#cvrt$fileSkipColumns = 1;       # one column of row labels
#if(length(covariates_file_name)>0) {
#  cvrt$LoadFile(covariates_file_name);
#}
# Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = " ";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 0;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name)
snps$columnNames=readLines(basedir%&%"qtl/Gaucasian.id")


## Run the analysis
pertimes=500
permchrALL=as.data.frame(matrix(nrow = nrow(phenotype.mat),ncol = pertimes,data = pvOutputThreshold_tra))
rownames(permchrALL)=rownames(phenotype.mat)
indx=rownames(permchrALL)
#chrcvrt=as.matrix(read.table(covariates_file_name,sep = "\t",header = T,row.names = 1))
for(i in 1:pertimes){
print(i%&%" times")
ranorder=order(runif(ncol(phenotype.mat)))
tmp=as.matrix(phenotype.mat[,ranorder])
#tmpcov=chrcvrt
#tmpcov[1,]=tmpcov[1,ranorder]
#cvrt$CreateFromMatrix(tmpcov)
gene$CreateFromMatrix(tmp)
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
  permchrALL[,i]=as.numeric(me$cis$min.pv.gene[indx])
}
  write.table(t(permchrALL),basedir%&%"qtl/permutation/result/permchr"%&%chr%&%"_thread_"%&%thread%&%".txt",col.names = T,row.names =F,quote = F,sep = " " )
