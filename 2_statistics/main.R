#######################
# date: 8/8/2018
# author: Zelin
# decription: statistical analysis part
######################
stop()
options(stringsAsFactors = F)
library(edgeR)
library(psych)
library(reshape2)
library(stats)
library(sva)
library(ggpubr)
library(missForest)
library(scatterplot3d)
library(RColorBrewer)
library(limma)
"%&%"=function(a,b)paste0(a,b)
baseAlp=c("A","T","G","C")
calcResiduals <- function(geneBySampleValues, samplesByCovariates, varsToAddBackIn=NULL, sampleWeights=NULL) {
  #################################################################################
  # "Recursively" use the calcResiduals() code after this section (in a for loop):
  #################################################################################
  if (is.matrix(sampleWeights)) {
    residualizedMat = matrix(NA, nrow=nrow(geneBySampleValues), ncol=ncol(geneBySampleValues), dimnames=dimnames(geneBySampleValues))
    for (gInd in 1:nrow(geneBySampleValues)) {
      gRow = calcResiduals(geneBySampleValues[gInd, , drop=FALSE], samplesByCovariates, varsToAddBackIn, sampleWeights[gInd, ])
      residualizedMat[gInd, ] = gRow
    }
    return(residualizedMat)
  }
  #################################################################################
  
  #result.lm = lsfit(x=samplesByCovariates, y=t(geneBySampleValues), wt=sampleWeights, intercept=FALSE)
  
  # Formula of "y ~ 0 + x" means no intercept:
  result.lm = lm(t(geneBySampleValues) ~ 0 + samplesByCovariates, weights=sampleWeights)
  covarNames = colnames(samplesByCovariates)
  
  coef = result.lm$coefficients
  isMatrixForm = is.matrix(coef)
  if (isMatrixForm) {
    rownames(coef) = covarNames
  }
  else {
    names(coef) = covarNames
  }
  
  allVarsToAddBack = '(Intercept)'
  allVarsToAddBack = c(allVarsToAddBack, varsToAddBackIn)
  allVarsToAddBack = intersect(allVarsToAddBack, covarNames)
  
  residualizedMat = result.lm$residuals
  if (isMatrixForm) {
    multCoef = coef[allVarsToAddBack, , drop=FALSE]
  }
  else {
    multCoef = coef[allVarsToAddBack]
  }
  residualizedMat = residualizedMat + samplesByCovariates[, allVarsToAddBack, drop=FALSE] %*% multCoef
  
  residualizedMat = t(residualizedMat)
  rownames(residualizedMat) = rownames(geneBySampleValues)
  colnames(residualizedMat) = colnames(geneBySampleValues)
  
  return(residualizedMat)
}
#################Hypothesis###############
CMC.geneSet=read.delim("/media/data3/CMC/data/network/CMC_GeneSets_Hypothesis-driven-for-Enrichement.tsv")
CMC.geneSet.list=list()
for (i in 1:nrow(CMC.geneSet)){
  CMC.geneSet.list[[CMC.geneSet[i,1]]]=strsplit(CMC.geneSet[i,"symBeforeOverlap"],"\\|")
}
CMC.geneSet.list.all=unique(c(unlist(CMC.geneSet.list[["GWAS:PGC2 SCZ"]]),
                              unlist(CMC.geneSet.list[["CNV:SCZ"]]),
                              unlist(CMC.geneSet.list[["De novos: SCZ NS"]]),
                              unlist(CMC.geneSet.list[["De novos:SCZ LoF"]])))

CMC.geneSet.list.GWASaCNV=unique(c(unlist(CMC.geneSet.list[["GWAS:PGC2 SCZ"]]),
                                   unlist(CMC.geneSet.list[["CNV:SCZ"]])))
CMC.geneSet.list.GWAS=unique(unlist(CMC.geneSet.list[["GWAS:PGC2 SCZ"]]))
CMC.geneSet.list.GWAS=CMC.geneSet.list.GWAS[-which(CMC.geneSet.list.GWAS=="NA")]
PRIMARY_VARS=c("DxControl","DxSCZ")
#read another result
DE.PBSI=read.delim("DE.paper",sep="\t")

#read mergeCirc Junction file
CIRI.JUNCTION=read.delim("/media/data3/circCMC/data/combineCIRI.junction.txt")
rownames(CIRI.JUNCTION)=CIRI.JUNCTION[,1]
CIRI.JUNCTION.mat=CIRI.JUNCTION[,-1]
rm(CIRI.JUNCTION)
CIRI.JUNCTION.mat[is.na(CIRI.JUNCTION.mat)]=0
CIRI.JUNCTION.read=rowMeans(CIRI.JUNCTION.mat>=MIN_CIRC_READ)
CIRI.JUNCTION.depth=colSums(CIRI.JUNCTION.mat)
CIRI.JUNCTION.depth.M1=colSums(CIRI.JUNCTION.mat*(CIRI.JUNCTION.mat>1))
CIRI.JUNCTION.depth.M2=colSums(CIRI.JUNCTION.mat*(CIRI.JUNCTION.mat>2))
CIRI.JUNCTION.depth.M3=colSums(CIRI.JUNCTION.mat*(CIRI.JUNCTION.mat>3))
CIRI.JUNCTION.depth.M4=colSums(CIRI.JUNCTION.mat*(CIRI.JUNCTION.mat>4))
CIRI.JUNCTION.depth.M5=colSums(CIRI.JUNCTION.mat*(CIRI.JUNCTION.mat>5))

# sum of circRNA RPKM
# circRNA RPKM = junction reads/(circRNA length × total mapped reads)
# circRNA length = (length of reads-20bp) × 2
CIRI.JUNCTION.sumRPKM=CIRI.JUNCTION.depth
hist(log(CIRI.JUNCTION.depth,10))
CIRI.JUNCTION.SampleNum=rowSums(CIRI.JUNCTION.mat>0)

#read mergeCirc Ratio file
CIRI.RATIO=read.delim("/media/data3/circCMC/data/combineCIRI.ratio.txt")
rownames(CIRI.RATIO)=CIRI.RATIO[,1]
CIRI.RATIO.mat=CIRI.RATIO[,-1]
rm(CIRI.JUNCTION)
CIRI.RATIO.mat[is.na(CIRI.RATIO.mat)]=0
CIRI.RATIO.10559.mat=CIRI.RATIO.mat[STA.paper$circRNA.ID,]
#read basic information for these circRNA
STA.paper=read.delim("STA.paper.txt",header = T)

#exon site and length revised
circ.annot=read.table("./network/circRNA.txt",header=F,sep="\t")
rownames(circ.annot)=circ.annot$V1%&%":"%&%circ.annot$V2%&%"|"%&%circ.annot$V3
circ.annot=circ.annot[STA.paper$circRNA.ID,]
circLen=vector(length = nrow(circ.annot))
for (i in 1:nrow(circ.annot)) {
  tmpend=as.numeric(unlist(strsplit(circ.annot[i,'V6'],",")))
  tmpstart=as.numeric(unlist(strsplit(circ.annot[i,'V5'],",")))
  circLen[i]=sum(tmpend-tmpstart)+length(tmpend)
}
circ.annot$len=circLen
exonIdx=which(STA.paper$region=="exon")
STA.paper[exonIdx,'exon.start.site']=circ.annot[exonIdx,'V5']
STA.paper[exonIdx,'exon.end.site']=circ.annot[exonIdx,'V6']
write.table(STA.paper,"STA.paper.txt",col.names = T,row.names = F,sep = "\t",quote=F)

STA.paper.exon=subset(STA.paper,region=="exon")
STA.paper.exon$LenPerExon=c(STA.paper.exon$predicted.length/STA.paper.exon$exon.number)
STA.paper.gene=STA.paper$gene
STA.paper.gene=unique(unlist(strsplit(STA.paper$gene,",")))
STA.paper.gene=STA.paper.gene[-which(STA.paper.gene=="n/a")]
STA.NatureS3=subset(STA.paper,gene%in%Nature.paper.S3$Ensembl.Gene.ID)
STA.noNatureS3=subset(STA.paper,!(gene%in%Nature.paper.S3$Ensembl.Gene.ID))
STA.paper.intergenic=subset(STA.paper,region=="intergenic_region")
STA.paper.exon=subset(STA.paper,region=="exon")
STA.paper.intron=subset(STA.paper,region=="intron")

write.table(STA.paper.gene,"./enrichment_C/STA.paper.gene",row.names = F,col.names = F,quote = F)

STA.pos=do.call("rbind",strsplit(STA.paper$circRNA.ID,":"))
tmp=do.call("rbind",strsplit(STA.pos[,2],"\\|"))
STA.pos=data.frame(geneID=STA.paper$circRNA.ID,chr=STA.pos[,1],left=tmp[,1],right=tmp[,2])
STA.pos$left=as.numeric(STA.pos$left)
STA.pos$right=as.numeric(STA.pos$right)
rownames(STA.pos)=STA.pos$geneID
STA.pos$strand=STA.paper$strand
write.table(STA.pos,"STA.pos",sep="\t",row.names = F,col.names = T,quote = FALSE)

#circRNA from one gene
STA.paper.clear=STA.paper[nchar(STA.paper$gene)==15,]
STA.paper.clearGene=unique(STA.paper.clear$gene)
writeLines(STA.paper.clearGene,"./enrichment_C/STA.paper.clearGene")

####################circRNA detail type################
STA.bed=data.frame(STA.pos[,c(2,3,4,1)],0,STA.paper[,'strand'])
STA.bed[,2]=STA.bed[,2]-1
write.table(STA.bed,"AnnotationCirc/STA.bed",col.names = F,row.names = F,quote=F,sep="\t")
STA.v19=read.delim("./AnnotationCirc/STA.gencodeV19.txt",header=F)
STA.v19$gene=substr(STA.v19$V15,start=9,stop=23)
STA.v19$transcript=substr(STA.v19$V15,start=42,stop=56)
STA.type=data.frame(matrix(nrow=nrow(STA.paper),ncol=9))
colnames(STA.type)=c("circRNA.ID","CDS","UTR5","UTR3","Intronic","Antisense",
                     "Pseudogene","LincRNA","Intergenic")
STA.type$circRNA.ID=STA.paper$circRNA.ID
rownames(STA.type)=STA.paper$circRNA.ID
#exon type

for (i in 1:nrow(STA.type)) {
  tmpcirc=STA.paper[i,'circRNA.ID']
  if(STA.paper[i,'region']=="exon"){
    tmpv19=subset(STA.v19,V4==tmpcirc&V13==V6)
    tmpv19.trans=subset(tmpv19,V9=="transcript")
    tmpv19.gene=subset(tmpv19,V9=="gene")
    
    tmpChoice=unique(tmpv19.trans[which(tmpv19.trans$V16==max(tmpv19.trans$V16)),'transcript'])
    if(length(tmpChoice)>0){
      tmpv19=subset(tmpv19,transcript%in%tmpChoice)
    }else{
      tmpChoice=unique(tmpv19.gene[which(tmpv19.gene$V16==max(tmpv19.gene$V16)),'gene'])
      tmpv19=subset(tmpv19,gene%in%tmpChoice)
    }
    tmp=unique(tmpv19[,'V9'])
    if(length(grep("pseudogene",tmpv19$V15))>0){
      STA.type[i,'Pseudogene']=1
    }
    if(length(grep("lincRNA",tmpv19$V15))>0){
      STA.type[i,'LincRNA']=1
    }
    if(length(grep("antisense",tmpv19$V15))>0){
      STA.type[i,'Antisense']=1
    }
    
    if("CDS"%in%tmp){
      STA.type[i,'CDS']=1
    }
    if("UTR3"%in%tmp){
      STA.type[i,'UTR3']=1
    }
    if("UTR5"%in%tmp){
      STA.type[i,'UTR5']=1
    }
  }else if(STA.paper[i,'region']=="intron"){
    STA.type[i,'Intronic']=1
  }else{
    tmpv19=subset(STA.v19,V4==tmpcirc)
    tmpv19Strand=unique(tmpv19$V13)
    if(length(tmpv19Strand)>1){
      STA.type[i,'Antisense']=1
    }else if(length(tmpv19Strand)==1){
      tmpStrand=STA.paper[i,'strand']
      if(tmpv19Strand==tmpStrand){
        STA.type[i,'Intergenic']=1
      }else if(tmpv19Strand=="."){
        STA.type[i,'Intergenic']=1
      }else{
        STA.type[i,'Antisense']=1
      }
      
    }else{
      STA.type[i,'Intergenic']=1
    }
    
  }
}
STA.type[is.na(STA.type)]=0
which(rowSums(STA.type[,2:9])==0)
STA.type.exon=STA.type[STA.paper.exon$circRNA.ID,c("CDS","UTR5","UTR3","Antisense","LincRNA","Pseudogene")]
for (i in 1:nrow(STA.type.exon)) {
  tmp=STA.type.exon[i,c("CDS","UTR5","UTR3")]
  if(sum(tmp)>0){
    STA.type.exon[i,'type']=paste(c("CDS","UTR5","UTR3")[which(tmp==1)],collapse = "_")
  }else{
    tmp=STA.type.exon[i,c("Antisense","LincRNA","Pseudogene")]
    STA.type.exon[i,'type']=paste(c("Antisense","LincRNA","Pseudogene")[which(tmp==1)],collapse = "_")
  }
}
STA.type.exon.statistic=table(STA.type.exon$type)
names(STA.type.exon.statistic)[which(names(STA.type.exon.statistic)=="")]="others"
STA.type.exon.statistic['CDS']=STA.type.exon.statistic['CDS']+STA.type.exon.statistic['others']
STA.type.exon.statistic['LincRNA']=STA.type.exon.statistic['LincRNA']+STA.type.exon.statistic['Antisense_LincRNA']
STA.type.exon.statistic['CDS_UTR5_UTR3']=STA.type.exon.statistic['CDS_UTR5_UTR3']+STA.type.exon.statistic['UTR5_UTR3']
pie.plot=STA.type.exon.statistic[c('CDS','CDS_UTR5','CDS_UTR3','Antisense','CDS_UTR5_UTR3','LincRNA','Pseudogene','UTR3','UTR5')]
pie(pie.plot)
###############################filter circRNA###################################
CIRI.JUNCTION.PASS.idx=match(STA.paper$circRNA.ID,rownames(CIRI.JUNCTION.mat))
CIRI.JUNCTION.PASS.mat=CIRI.JUNCTION.mat[CIRI.JUNCTION.PASS.idx,]
CIRI.JUNCTION.PASS.depth=colSums(CIRI.JUNCTION.PASS.mat)
rm(CIRI.JUNCTION.mat)

CIRI.SAMPLE.Observe=rep(0,each=nrow(CIRI.JUNCTION.PASS.mat))
for (i in 1:nrow(CIRI.JUNCTION.PASS.mat)) {
  CIRI.SAMPLE.Observe[i]=sum(CIRI.JUNCTION.PASS.mat[i,]!=0)
}
CIRI.SAMPLE.Observe=data.frame(count=CIRI.SAMPLE.Observe,type=STA.paper$region)
CIRI.SAMPLE.Observe.Cum=data.frame(x=sort(unique(CIRI.SAMPLE.Observe$count)),
                                   y=cumsum(table(CIRI.SAMPLE.Observe$count))/nrow(CIRI.SAMPLE.Observe))
CIRI.SAMPLE.Observe.all=data.frame(CIRI.SAMPLE.Observe,
                                   cum=CIRI.SAMPLE.Observe.Cum[as.character(CIRI.SAMPLE.Observe$count),2])
CIRI.SAMPLE.Observe.all$cum=CIRI.SAMPLE.Observe.all$cum*1250
{
  pdf("./pdf_C/number_sample.pdf",width = 5,height=5)
  par(mar=c(5,5,2,2))
  tmp=ggplot(CIRI.SAMPLE.Observe.all)+
    geom_histogram(aes(count,fill=type,color=type),binwidth =10)+theme_pubr()+scale_fill_manual(values = ggplot2::alpha(c("#F8766D","#00BA38","#619CFF"),.7))+
    geom_line(aes(count,cum),color="black",linetype="dashed")+
    theme(legend.title = element_text(size = 0))+labs(x="Samples",y="Number of circRNAs")+ 
    geom_hline(aes(yintercept=625),colour="grey",linetype="dashed")+
    geom_vline(aes(xintercept=495),colour="grey",linetype="dashed")+
    scale_y_continuous(
      sec.axis = sec_axis(~ . * 1 / 1250 , name = "Cumulative fraction"), 
      limits = c(0, 1250))
  print(tmp)
  dev.off()
}
#########################read mergeCirc circ-to-linear file###################
CIRI.C2I=read.delim("/media/data3/circCMC/data/combineCIRI.ratio.txt")
CIRI.C2I.mat=CIRI.C2I[CIRI.JUNCTION.PASS.idx,-1]
rownames(CIRI.C2I.mat)=rownames(CIRI.JUNCTION.PASS.mat)
rm(CIRI.C2I)

#read mergeCirc linear file
CIRI.Linear=read.delim("/media/data3/circCMC/data/combineCIRI.linear.txt")
CIRI.Linear.depth=colSums(CIRI.Linear[,-1],na.rm=T)
CIRI.Linear.mat=CIRI.Linear[CIRI.JUNCTION.PASS.idx,-1]
CIRI.Linear.depth.pass=colSums(CIRI.Linear.mat,na.rm=T)
rownames(CIRI.Linear.mat)=rownames(CIRI.JUNCTION.PASS.mat)
rm(CIRI.Linear)

#######################read covariates##################
covariates=read.delim("/media/data3/CMC/data/covariate.txt")
covariates[which(covariates$Dx=="BP"),'Dx']="AFF"
covariates=covariates[match(colnames(CIRI.JUNCTION.PASS.mat),covariates$id),]
rownames(covariates)=covariates$id
covariates.SCZ=covariates[covariates$Dx!="AFF",]
covariates.AFF.SCZ=covariates[covariates$Dx!="Control",]
covariates.AFF=covariates[covariates$Dx!="SCZ",]
rownames(covariates.SCZ)=covariates.SCZ$id

names(covariates)
##############################design matrix################################
Circ.design.all=model.matrix(as.formula("~ . "),covariates[,-1])
rownames(Circ.design.all)=covariates$id
Circ.design.SCZwithControl=model.matrix(as.formula("~ . "),covariates.SCZ[,-1])
Circ.design.sig.SCZ.Group=model.matrix(~as.factor(covariates.SCZ$Dx))
Circ.design.AFFwithSCZ=model.matrix(as.formula("~ ."),covariates.AFF.SCZ[,-1])
Circ.design.AFFwithControl=model.matrix(as.formula("~ ."),covariates.AFF[,-1])

######################count normalize ALL#######################
CIRI.JUNCTION.PASS.DEG=DGEList(CIRI.JUNCTION.PASS.mat[,covariates$id],
                               lib.size = CIRI.JUNCTION.depth[covariates$id],
                               genes=rownames(CIRI.JUNCTION.PASS.mat))
CIRI.JUNCTION.PASS.CPM=cpm(CIRI.JUNCTION.PASS.DEG)
CPM.mean.SCZ=rowMeans(CIRI.JUNCTION.PASS.CPM[,which(covariates$Dx=='SCZ')])
CPM.mean.Control=rowMeans(CIRI.JUNCTION.PASS.CPM[,which(covariates$Dx=='Control')])
logCPM.mean=log2(CPM.mean.SCZ/CPM.mean.Control)
#sdlogCPM=apply(log2(CIRI.JUNCTION.PASS.CPM[,which(covariates$Dx%in%c('Control','SCZ'))]),1,sd)
CIRI.JUNCTION.PASS.voom=voom(CIRI.JUNCTION.PASS.DEG,design=Circ.design.all)
CIRI.JUNCTION.PASS.voom.log=CIRI.JUNCTION.PASS.voom$E
CIRI.JUNCTION.PASS.voom.weight=CIRI.JUNCTION.PASS.voom$weights
CIRI.JUNCTION.PASS.voom.log.Mean=rowMeans(CIRI.JUNCTION.PASS.voom.log)
write.table(CIRI.JUNCTION.PASS.voom.log.Mean,"CIRI.JUNCTION.PASS.voom.log.Mean",sep="\t")
#################SCZ and Control filter separately####
CIRI.JUNCTION.SCZ.mat=CIRI.JUNCTION.mat[,rownames(covariates.SCZ)[which(covariates.SCZ$Dx=="SCZ")]]
CIRI.JUNCTION.Control.mat=CIRI.JUNCTION.mat[,rownames(covariates.SCZ)[which(covariates.SCZ$Dx=="Control")]]

CIRI.RATIO.SCZ.mat=CIRI.RATIO.mat[,rownames(covariates.SCZ)[which(covariates.SCZ$Dx=="SCZ")]]
CIRI.RATIO.Control.mat=CIRI.RATIO.mat[,rownames(covariates.SCZ)[which(covariates.SCZ$Dx=="Control")]]

CIRI.JUNCTION.SCZ.mean=rowMeans(CIRI.JUNCTION.SCZ.mat>=1)
CIRI.JUNCTION.Control.mean=rowMeans(CIRI.JUNCTION.Control.mat>=1)

CIRI.RATIO.SCZ.mean=rowMeans(CIRI.RATIO.SCZ.mat)
CIRI.RATIO.Control.mean=rowMeans(CIRI.RATIO.Control.mat)

idx.SCZ=which(CIRI.JUNCTION.SCZ.mean>=0.5 & CIRI.RATIO.SCZ.mean>0.05)
idx.Control=which(CIRI.JUNCTION.Control.mean>=0.5 & CIRI.RATIO.Control.mean>0.05)
idx.dff=unique(setdiff(c(idx.SCZ,idx.Control),intersect(idx.SCZ,idx.Control)))
id.dff=rownames(CIRI.RATIO.SCZ.mat)[idx.dff]
tmp=intersect(STA.paper$circRNA.ID,id.dff)
id.dff=setdiff(id.dff,tmp)
CIRI.JUNCTION.SPECIFIC.DEG=DGEList(CIRI.JUNCTION.mat[id.dff,covariates.SCZ$id],
                               lib.size = CIRI.JUNCTION.depth[covariates.SCZ$id],
                               genes=id.dff)
CIRI.JUNCTION.SPECIFIC.voom=voom(CIRI.JUNCTION.SPECIFIC.DEG,design=Circ.design.SCZwithControl)
CIRI.JUNCTION.SPECIFIC.voom.log=CIRI.JUNCTION.SPECIFIC.voom$E
CIRI.JUNCTION.SPECIFIC.voom.lmFit=lmFit(CIRI.JUNCTION.SPECIFIC.voom,design=Circ.design.SCZwithControl)
CIRI.JUNCTION.SPECIFIC.voom.eBayes=eBayes(CIRI.JUNCTION.SPECIFIC.voom.lmFit)
CIRI.JUNCTION.SPECIFIC.voom.DGE=toptable(CIRI.JUNCTION.SPECIFIC.voom.eBayes,coef="DxSCZ",n=Inf,confint=T,sort.by = "p")

#####################Cumulative fraction plot##################

{
  pdf("./pdf_C/cumulative_fraction_circRNA.pdf",width = 5,height = 5)
  par(mar=c(5,5,1,1),mfrow=c(1,1))
  plot(ecdf(CIRI.JUNCTION.PASS.voom.log.Mean[which(STA.paper$region=="exon")]),ylab="Cumulative fraction",xlab="Average circRNA expression (log2 CPM)",col="#F8766D",main="",do.points=F,verticals=T,xlim=c(0,10),las=1)
  par(new=T)
  plot(ecdf(CIRI.JUNCTION.PASS.voom.log.Mean[which(STA.paper$region=="intergenic_region")]),ylab="",xlab="",col="#00BA38",main="",axes=F,do.points=F,verticals=T,xlim=c(0,10))
  par(new=T)
  plot(ecdf(CIRI.JUNCTION.PASS.voom.log.Mean[which(STA.paper$region=="intron")]),ylab="",xlab="",col="#619CFF",main="",axes=F,do.points=F,verticals=T,xlim=c(0,10))
  text(8,0.07,labels = "Exonic",col = "#F8766D", cex=1,pos=4)
  text(8,0.12,labels = "Intergenic",col = "#00BA38", cex=1,pos=4)
  text(8,0.17,labels = "Intronic",col = "#619CFF", cex=1,pos=4)
  dev.off()
  
}

CIRI.expression.distribution=data.frame(exp=CIRI.JUNCTION.PASS.voom.log.Mean,type=STA.paper$region)
CIRI.expression.distribution.cum=data.frame(x=sort(unique(CIRI.JUNCTION.PASS.voom.log.Mean)),
                                            y=cumsum(table(CIRI.JUNCTION.PASS.voom.log.Mean))/10559)
CIRI.expression.distribution.all=data.frame(CIRI.expression.distribution,
                                            cum=CIRI.expression.distribution.cum[as.character(CIRI.expression.distribution$exp),2])
CIRI.expression.distribution.all$cum=CIRI.expression.distribution.all$cum * 1720
{
  pdf("./pdf_C/Distribution_circRNA.pdf")
  par(mar=c(5,5,2,2))
  tmp=ggplot(CIRI.expression.distribution.all)+
    geom_histogram(aes(exp,fill=type,color=type),binwidth =0.4)+theme_pubr()+scale_fill_manual(values = ggplot2::alpha(c("#F8766D","#00BA38","#619CFF"),.7))+
    geom_line(aes(exp,cum),color="black",linetype="dashed")+
    theme(legend.title = element_text(size = 0))+labs(x="Average circRNA expression (log2 CPM)",y="Frequence")+
    scale_y_continuous(
      sec.axis = sec_axis(~ . * 1 / 1720 , name = "Cumulative fraction"), 
      limits = c(0, 1740))
  print(tmp)
  dev.off()
}

#####################SCZ normalization########################
CIRI.JUNCTION.PASS.SCZ.voom=voom(CIRI.JUNCTION.PASS.DEG[,covariates.SCZ$id,keep.lib.sizes=T],design=Circ.design.SCZwithControl,plot=T)
CIRI.JUNCTION.PASS.SCZ.voom.log=CIRI.JUNCTION.PASS.SCZ.voom$E
CIRI.JUNCTION.PASS.SCZ.voom.weight=CIRI.JUNCTION.PASS.SCZ.voom$weights
CIRI.JUNCTION.PASS.SCZ.voom.log.rowMeans=rowMeans(CIRI.JUNCTION.PASS.SCZ.voom.log)

#####################read  Tophat2 raw quantify counts file###############
Tophat2.raw.counts=read.delim("/media/data3/CMC/data/expression/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_geneExpressionRaw.tsv.gz")
Tophat2.raw.counts=Tophat2.raw.counts[,covariates$id]
Tophat2.raw.counts.LIB_SIZE=colSums(Tophat2.raw.counts)

####################Tophat2 counts to log(cpm)##########################
Tophat2.DGEList=DGEList(Tophat2.raw.counts,lib.size = Tophat2.raw.counts.LIB_SIZE,genes=rownames(Tophat2.raw.counts))
Tophat2.DGEList.cpm=cpm(Tophat2.DGEList)
idx=which(rowMeans(Tophat2.DGEList.cpm-1>0)>=0.5)
Tophat2.DGEList.sub=Tophat2.DGEList[idx,,keep.lib.sizes=T]
writeLines(Tophat2.DGEList.sub$genes$genes,"./enrichment_C/Tophat2.HE.gene")
Tophat2.DGEList.voom=voom(Tophat2.DGEList.sub,design = Circ.design.all,plot = T)
Tophat2.DGEList.voom.log=Tophat2.DGEList.voom$E
Tophat2.DGEList.voom.weight=Tophat2.DGEList.voom$weights
Tophat2.DGEList.sub.genes=rownames(Tophat2.DGEList.voom.log)

Tophat2.DGEList.SCZ=Tophat2.DGEList.sub[,covariates.SCZ$id,keep.lib.sizes=T]
Tophat2.DGEList.SCZ.voom=voom(Tophat2.DGEList.SCZ,design = Circ.design.SCZwithControl,plot = T)
Tophat2.DGEList.SCZ.voom.log=Tophat2.DGEList.SCZ.voom$E
###############read isoform expression######################
Isoform.log.cpm=read.delim("/media/data3/CMC/data/expression/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_isoform-adjustedNoSVA-dataNormalization-includeAncestry-adjustedLogCPM.tsv.gz")
rownames(Isoform.log.cpm)=Isoform.log.cpm$IsoformFeature
Isoform.log.cpm=Isoform.log.cpm[,-1]
Isoform.log.cpm=Isoform.log.cpm[,covariates$id]
#######################read GencodeV19#########################
gtfV19=read.delim("/media/data3/annotation/GENCODE/gencode.v19.annotation.gtf.gz",comment.char = "#",header=F)
gtfV19.exon=subset(gtfV19,V3=="exon")
gtfV19.gene=subset(gtfV19,V3=="gene")
gtfV19.gene$gene=substr(gtfV19.gene$V9,9,23)
gtfV19.gene=gtfV19.gene[!duplicated(gtfV19.gene$gene),]
gtfV19.transcript=subset(gtfV19,V3=="transcript")
gtfV19.transcript.detail=do.call("rbind",strsplit(gtfV19.transcript$V9,split = ";"))
gtfV19.transcript.detail=as.data.frame(gtfV19.transcript.detail)
gtfV19.transcript.detail$transID=substr(do.call("rbind",strsplit(gtfV19.transcript.detail[,2]," "))[,3],1,15)
gtfV19.transcript.detail$transType=do.call("rbind",strsplit(gtfV19.transcript.detail[,5]," "))[,3]
gtfV19.transcript.detail$geneID=substr(do.call("rbind",strsplit(gtfV19.transcript.detail[,1]," "))[,2],1,15)

###################Depth statistic#################################
Depth.merge=data.frame(Gene=Tophat2.raw.counts.LIB_SIZE+CIRI.JUNCTION.depth,
                       Circ=CIRI.JUNCTION.depth,
                       Linear=CIRI.Linear.depth)
colMeans(Depth.merge[which(covariates$Dx=="Control"),])
apply(Depth.merge[which(covariates$Dx=="Control"),], 2, sd)
colMeans(Depth.merge[which(covariates$Dx=="SCZ"),])
apply(Depth.merge[which(covariates$Dx=="SCZ"),], 2, sd)
colMeans(Depth.merge[which(covariates$Dx=="AFF"),])
apply(Depth.merge[which(covariates$Dx=="AFF"),], 2, sd)


Depth.merge=Depth.merge[order(Depth.merge$Gene,decreasing = T),]

Depth.merge.long=data.frame(id=rep(1:nrow(Depth.merge),3),depth=c(Depth.merge$Gene,Depth.merge$Linear,Depth.merge$Circ),
                            type=rep(c("a","b","c"),each=nrow(Depth.merge)))
Depth.merge.long$col=c("#F8766D","#00BA38","#619CFF")[as.numeric(as.factor(Depth.merge.long$type))]
Depth.merge.long=subset(Depth.merge.long,type!="b")
{
  pdf("./pdf_C/sample_depth_distribution.pdf")
  plot(log(depth,10)~id,Depth.merge.long,type="h",col=Depth.merge.long$col,
       las=1,bty="l",xlab="samples",ylab="reads",xaxt="n")
  dev.off()
}


##################Linear counts normalize###########################
CIRI.Linear.noNA.mat=CIRI.Linear.mat
CIRI.Linear.noNA.mat[is.na(CIRI.Linear.mat)]=0
CIRI.Linear.DEG=DGEList(CIRI.Linear.noNA.mat,
                        lib.size = CIRI.Linear.depth,
                        genes=rownames(CIRI.Linear.mat))
CIRI.Linear.voom=voom(CIRI.Linear.DEG,design = Circ.design.all,plot = T)
CIRI.Linear.voom.log=CIRI.Linear.voom$E
CIRI.Linear.voom.NA.log=CIRI.Linear.voom.log
CIRI.Linear.voom.NA.log[is.na(CIRI.Linear.mat)]=NA

###################Focus1 DxSCZ#######################
##########################Dx from limma##################
# Tophat2.DGEList.voom.lmFit=lmFit(Tophat2.DGEList.voom,design=Circ.design.all)
# Tophat2.DGEList.voom.eBays=eBayes(Tophat2.DGEList.voom.lmFit)
# for (i in 1:ncol(Circ.design.all)) {
#   Tophat2.DGEList.voom.DEG=toptable(Tophat2.DGEList.voom.eBays,coef=i,n=Inf,confint=T,sort.by = "p")
#   print(paste(colnames(Circ.design.all)[i],sum(Tophat2.DGEList.voom.DEG$adj.P.Val<0.05)))
# }


Tophat2.DGEList.SCZ.voom.lmFit=lmFit(Tophat2.DGEList.SCZ.voom,design=Circ.design.SCZwithControl)
Tophat2.DGEList.SCZ.voom.eBays=eBayes(Tophat2.DGEList.SCZ.voom.lmFit)
Tophat2.DGEList.SCZ.voom.list=list()
for (i in 1:ncol(Circ.design.SCZwithControl)) {
  Tophat2.DGEList.SCZ.voom.DEG=toptable(Tophat2.DGEList.SCZ.voom.eBays,coef=i,n=Inf,confint=T,sort.by = "p")
  Tophat2.DGEList.SCZ.voom.list[[colnames(Circ.design.SCZwithControl)[i]]]=rownames(Tophat2.DGEList.SCZ.voom.DEG)[which(Tophat2.DGEList.SCZ.voom.DEG$adj.P.Val<0.05)]
  print(paste(colnames(Circ.design.SCZwithControl)[i],sum(Tophat2.DGEList.SCZ.voom.DEG$adj.P.Val<0.05)))
}
#Dx
Gene.DEG=toptable(Tophat2.DGEList.SCZ.voom.eBays,coef="DxSCZ",n=Inf,confint=T,sort.by = "p")
Gene.DEG.sig=subset(Gene.DEG,adj.P.Val<0.05)
length(intersect(Nature.paper.S3$Ensembl.Gene.ID,rownames(Gene.DEG.sig)))
length(Tophat2.DGEList.SCZ.voom.list[["DxSCZ"]])
#Institution
length(unique(union(Tophat2.DGEList.SCZ.voom.list[["InstitutionPenn"]],
                    Tophat2.DGEList.SCZ.voom.list[["InstitutionPitt"]])))
#LIB
length(unique(c(Tophat2.DGEList.SCZ.voom.list[["libBatchA"]],
                Tophat2.DGEList.SCZ.voom.list[["libBatchB"]],
                Tophat2.DGEList.SCZ.voom.list[["libBatchC"]],
                Tophat2.DGEList.SCZ.voom.list[["libBatchD"]],
                Tophat2.DGEList.SCZ.voom.list[["libBatchE"]],
                Tophat2.DGEList.SCZ.voom.list[["libBatchF"]],
                Tophat2.DGEList.SCZ.voom.list[["libBatchG"]],
                Tophat2.DGEList.SCZ.voom.list[["libBatchH"]])))
#Gender
length(Tophat2.DGEList.SCZ.voom.list[["SexMale"]])
#AOD
length(Tophat2.DGEList.SCZ.voom.list[["AOD"]])
Tophat2.DGEList.SCZ.voom.AOD=toptable(Tophat2.DGEList.SCZ.voom.eBays,coef="AOD",n=Inf,confint=T,sort.by = "p")

#PMI
length(Tophat2.DGEList.SCZ.voom.list[["PMI"]])
#RIN
length(unique(union(Tophat2.DGEList.SCZ.voom.list[["RIN"]],
                    Tophat2.DGEList.SCZ.voom.list[["RIN2"]])))
#Ethnicity
length(unique(c(Tophat2.DGEList.SCZ.voom.list[["EV1"]],
                Tophat2.DGEList.SCZ.voom.list[["EV2"]],
                Tophat2.DGEList.SCZ.voom.list[["EV3"]],
                Tophat2.DGEList.SCZ.voom.list[["EV4"]],
                Tophat2.DGEList.SCZ.voom.list[["EV5"]])))

CIRI.JUNCTION.PASS.SCZ.voom.lmFit=lmFit(CIRI.JUNCTION.PASS.SCZ.voom,design=Circ.design.SCZwithControl)
CIRI.JUNCTION.PASS.SCZ.voom.eBayes=eBayes(CIRI.JUNCTION.PASS.SCZ.voom.lmFit)
CIRI.JUNCTION.PASS.SCZ.voom.list=list()
for (i in 1:ncol(Circ.design.SCZwithControl)) {
  CIRI.JUNCTION.PASS.SCZ.voom.DGE=toptable(CIRI.JUNCTION.PASS.SCZ.voom.eBayes,
                                           coef=i,n=Inf,confint=T,sort.by = "p")
  CIRI.JUNCTION.PASS.SCZ.voom.list[[colnames(Circ.design.SCZwithControl)[i]]]=rownames(CIRI.JUNCTION.PASS.SCZ.voom.DGE)[which(CIRI.JUNCTION.PASS.SCZ.voom.DGE$adj.P.Val<0.05)]
  print(paste(colnames(Circ.design.SCZwithControl)[i],sum(CIRI.JUNCTION.PASS.SCZ.voom.DGE$adj.P.Val<0.05)))
}
names(CIRI.JUNCTION.PASS.SCZ.voom.list)
#Dx
CIRI.JUNCTION.PASS.SCZ.voom.list[["DxSCZ"]]
#Institution
length(unique(union(CIRI.JUNCTION.PASS.SCZ.voom.list[["InstitutionPenn"]],
      CIRI.JUNCTION.PASS.SCZ.voom.list[["InstitutionPitt"]])))
#LIB
length(unique(c(CIRI.JUNCTION.PASS.SCZ.voom.list[["libBatchA"]],
                    CIRI.JUNCTION.PASS.SCZ.voom.list[["libBatchB"]],
                    CIRI.JUNCTION.PASS.SCZ.voom.list[["libBatchC"]],
                    CIRI.JUNCTION.PASS.SCZ.voom.list[["libBatchD"]],
                    CIRI.JUNCTION.PASS.SCZ.voom.list[["libBatchE"]],
                    CIRI.JUNCTION.PASS.SCZ.voom.list[["libBatchF"]],
                    CIRI.JUNCTION.PASS.SCZ.voom.list[["libBatchG"]],
                    CIRI.JUNCTION.PASS.SCZ.voom.list[["libBatchH"]])))
#Gender
length(CIRI.JUNCTION.PASS.SCZ.voom.list[["SexMale"]])
#AOD
length(CIRI.JUNCTION.PASS.SCZ.voom.list[["AOD"]])
#PMI
length(CIRI.JUNCTION.PASS.SCZ.voom.list[["PMI"]])
#RIN
length(unique(union(CIRI.JUNCTION.PASS.SCZ.voom.list[["RIN"]],
                    CIRI.JUNCTION.PASS.SCZ.voom.list[["RIN2"]])))
#Ethnicity
length(unique(c(CIRI.JUNCTION.PASS.SCZ.voom.list[["EV1"]],
                CIRI.JUNCTION.PASS.SCZ.voom.list[["EV2"]],
                CIRI.JUNCTION.PASS.SCZ.voom.list[["EV3"]],
                CIRI.JUNCTION.PASS.SCZ.voom.list[["EV4"]],
                CIRI.JUNCTION.PASS.SCZ.voom.list[["EV5"]])))


CIRI.JUNCTION.PASS.voom.lmFit=lmFit(CIRI.JUNCTION.PASS.voom,design=Circ.design.all)
CIRI.JUNCTION.PASS.voom.eBayes=eBayes(CIRI.JUNCTION.PASS.voom.lmFit)
for (i in 1:ncol(Circ.design.all)) {
  CIRI.JUNCTION.PASS.voom.DEG=toptable(CIRI.JUNCTION.PASS.voom.eBayes,
                                       coef=i,n=Inf,confint=T,sort.by = "p")
  print(paste(colnames(Circ.design.all)[i],sum(CIRI.JUNCTION.PASS.voom.DEG$adj.P.Val<0.05)))
}
CIRI.JUNCTION.PASS.MF.voom.DEG=toptable(CIRI.JUNCTION.PASS.MF.voom.eBayes,
                                        coef="AOD",n=Inf,confint=T,sort.by = "p")

#####################################Schizophrenia#################
SCZ.DGE=toptable(CIRI.JUNCTION.PASS.SCZ.voom.eBayes,
                                         coef="DxSCZ",n=Inf,confint=T,sort.by = "p")
SCZ.DGE.P005=subset(SCZ.DGE,P.Value<0.05)
SCZ.DGE.P005.STA=STA.paper[rownames(SCZ.DGE.P005),]
SCZ.DGE.P005.STA.gene=na.omit(unique(unlist(strsplit(SCZ.DGE.P005.STA$geneSymbol,","))))
SCZ.DGE.P005.STA.gene=SCZ.DGE.P005.STA.gene[-which(SCZ.DGE.P005.STA.gene=="NA")]
m=length(unique(intersect(CMC.geneSet.list.GWAS,SCZ.DGE.P005.STA.gene)))
m
#wilcox test
Wt.P=vector(length = nrow(STA.paper))
for (i in 1:nrow(STA.paper)) {
  Wt.P[i]=wilcox.test(CIRI.JUNCTION.PASS.SCZ.voom.log[i,]~covariates.SCZ$Dx)$p.value
}
names(Wt.P)=STA.paper$circRNA.ID
Wt.fdr=p.adjust(Wt.P,"fdr")
hist(-log(Wt.P,10))
hist(Wt.P)
sum(Wt.fdr<0.05)

Wt.fdr005.STA=STA.paper[which(Wt.fdr<0.05),]

Wt.fdr005.DEG=SCZ.DGE[rownames(Wt.fdr005.STA),]
Wt.fdr005.DEG.sig=subset(Wt.fdr005.STA.DEG,P.Value<0.05)
Wt.fdr005.DEG.sig.STA=STA.paper[rownames(Wt.fdr005.DEG.sig),]

Gene.DEG.PSig=subset(Gene.DEG,P.Value<0.05)
idxGene=intersect(rownames(Gene.DEG.noSig),Wt.fdr005.DEG.sig.STA$gene)
Wt.fdr005.DEG.sig.STA.noGenesig=Wt.fdr005.DEG.sig.STA[-which(Wt.fdr005.DEG.sig.STA$gene%in%idxGene),]

#by chance
Wt.fdr005.DEG.sig.STA.noGenesig.gene=na.omit(unique(unlist(strsplit(Wt.fdr005.DEG.sig.STA.noGenesig$geneSymbol,","))))
Wt.fdr005.DEG.sig.STA.noGenesig.gene=Wt.fdr005.DEG.sig.STA.noGenesig.gene[-which(Wt.fdr005.DEG.sig.STA.noGenesig.gene=="NA")]
m=length(unique(intersect(CMC.geneSet.list.GWAS,Wt.fdr005.DEG.sig.STA.noGenesig.gene)))
tmpV=rep(0,10000)
n=nrow(Wt.fdr005.DEG.sig.STA.noGenesig)
for (i in 1:10000) {
  randA=na.omit(unique(unlist(strsplit(sample(STA.noNatureS3[,'geneSymbol'],n),","))))
  n1=length(intersect(CMC.geneSet.list.GWAS,randA))
  tmpV[i]=n1  
}
table(tmpV)
sum(tmpV>=m)
#the DE gene
idxGene=intersect(CMC.geneSet.list.GWAS,Wt.fdr005.DEG.sig.STA.noGenesig.gene)
SCZ.DEG.GWAS=
  Wt.fdr005.DEG.sig.STA.noGenesig[match(idxGene,Wt.fdr005.DEG.sig.STA.noGenesig$geneSymbol),]
boxplot(as.numeric(CIRI.JUNCTION.PASS.SCZ.voom.log[SCZ.DEG.GWAS$circRNA.ID[1],])~covariates.SCZ$Dx)
boxplot(as.numeric(CIRI.JUNCTION.PASS.SCZ.voom.log[SCZ.DEG.GWAS$circRNA.ID[2],])~covariates.SCZ$Dx)
boxplot(as.numeric(CIRI.JUNCTION.PASS.SCZ.voom.log[SCZ.DEG.GWAS$circRNA.ID[3],])~covariates.SCZ$Dx)
######################volcano plot#################################
vol.df=data.frame(FoldChange=as.numeric(SCZ.DGE$logFC))
vol.df$limmaP=-log(SCZ.DGE$P.Value,10)
vol.df$wilcoxonFDR=-log(Wt.fdr[match(rownames(SCZ.DGE),STA.paper$circRNA.ID)],10)
vol.df$color=rgb(30,30,30,maxColorValue = 255,alpha = 100)
vol.df[match(Wt.fdr005.DEG.sig.STA.noGenesig$circRNA.ID,rownames(SCZ.DGE)),'color']=rgb(215,48,31,maxColorValue = 255,alpha = 250)
vol.df[match(SCZ.DEG.GWAS$circRNA.ID,rownames(SCZ.DGE)),'color']="blue"
vol.df$size=1
vol.df[match(Wt.fdr005.DEG.sig.STA.noGenesig$circRNA.ID,rownames(SCZ.DGE)),'size']=1.3
write.table(vol.df,"vol.df.txt",row.names = F,col.names = T,quote=F,sep="\t")
{
  pdf("./pdf_C/SCZ_DE.pdf")
  s3d=scatterplot3d(y=vol.df$wilcoxonFDR,z=vol.df$limmaP,x=vol.df$FoldChange,color = vol.df$color,
                    xlab="FoldChange",
                    ylab="Wilcoxon test FDR",
                    zlab="Limma P-value",
                    box=F,grid = F,pch="")
  addgrids3d(y=vol.df$wilcoxonFDR,z=vol.df$limmaP,x=vol.df$FoldChange,grid = c( "yz","xz","xy"))
  for (i in 1:nrow(vol.df)) {
    s3d$points3d(y=vol.df$wilcoxonFDR[i],z=vol.df$limmaP[i],x=vol.df$FoldChange[i],
                 highlight.3d = "p",col = vol.df$color[i],pch=20)
    
  }
  s3d$box3d()
  dev.off()
}
circname='chr7:137148235|137308300'
SCZ.df=data.frame(value=as.numeric(CIRI.JUNCTION.PASS.SCZ.voom.log[circname,]),group=covariates.SCZ$Dx)
{
  pdf("pdf_C/SCZ_example.pdf",width=6)
  draw.ggbox=ggplot(SCZ.df,aes(x=group,y=value))+geom_boxplot(aes(fill=group))+
    theme_bw(base_size = 12)+theme(legend.position="none")+labs(title =paste(STA.paper[circname,"geneSymbol"] ," ",circname))
  print(draw.ggbox)
  dev.off()
}

#combine DE gene
Nature.S3.incircRNA=intersect(Nature.paper.S3$Ensembl.Gene.ID,STA.paper.clear$gene)
Nature.S3.incircRNA.STA=subset(STA.paper.clear,gene%in%Nature.S3.incircRNA)
SCZ.DGE.inNature=SCZ.DGE[Nature.S3.incircRNA.STA$circRNA.ID,]
SCZ.DGE.inNature$fdr=p.adjust(SCZ.DGE.inNature$P.Value,'fdr')
######################################LIB######################################
libBatchA.DEG=subset(toptable(CIRI.JUNCTION.PASS.SCZ.voom.eBayes,
                 coef="libBatchA",n=Inf,confint=T,sort.by = "p"),adj.P.Val<=0.05)
libBatchB.DEG=subset(toptable(CIRI.JUNCTION.PASS.SCZ.voom.eBayes,
                       coef="libBatchB",n=Inf,confint=T,sort.by = "p"),adj.P.Val<=0.05)
libBatchC.DEG=subset(toptable(CIRI.JUNCTION.PASS.SCZ.voom.eBayes,
                       coef="libBatchC",n=Inf,confint=T,sort.by = "p"),adj.P.Val<=0.05)
libBatchD.DEG=subset(toptable(CIRI.JUNCTION.PASS.SCZ.voom.eBayes,
                       coef="libBatchD",n=Inf,confint=T,sort.by = "p"),adj.P.Val<=0.05)
libBatchE.DEG=subset(toptable(CIRI.JUNCTION.PASS.SCZ.voom.eBayes,
                       coef="libBatchE",n=Inf,confint=T,sort.by = "p"),adj.P.Val<=0.05)
libBatchF.DEG=subset(toptable(CIRI.JUNCTION.PASS.SCZ.voom.eBayes,
                       coef="libBatchF",n=Inf,confint=T,sort.by = "p"),adj.P.Val<=0.05)
libBatchG.DEG=subset(toptable(CIRI.JUNCTION.PASS.SCZ.voom.eBayes,
                       coef="libBatchG",n=Inf,confint=T,sort.by = "p"),adj.P.Val<=0.05)
libBatchH.DEG=subset(toptable(CIRI.JUNCTION.PASS.SCZ.voom.eBayes,
                       coef="libBatchH",n=Inf,confint=T,sort.by = "p"),adj.P.Val<=0.05)
libBatchAll.DEG=unique(c(rownames(libBatchA.DEG),rownames(libBatchB.DEG),rownames(libBatchC.DEG),
                  rownames(libBatchD.DEG),rownames(libBatchE.DEG),rownames(libBatchF.DEG),
                  rownames(libBatchG.DEG),rownames(libBatchH.DEG)))
nolibBatchALL.DEG=setdiff(STA.paper$circRNA.ID,libBatchAll.DEG)
t.test(CIRI.JUNCTION.PASS.SCZ.voom.log.rowMeans[libBatchAll.DEG],
  CIRI.JUNCTION.PASS.SCZ.voom.log.rowMeans[nolibBatchALL.DEG])
t.test(STA.paper[libBatchAll.DEG,'predicted.length'],
       STA.paper[nolibBatchALL.DEG,'predicted.length'])
t.test(STA.paper[libBatchAll.DEG,'exon.number'],
       STA.paper[nolibBatchALL.DEG,'exon.number'])

libBatch.mean.exp.random=vector(length = 10000)
for (i in 1:10000) {
  libBatch.mean.exp.random[i]=mean(sample(CIRI.JUNCTION.PASS.SCZ.voom.log.rowMeans,4143))
}
mean(CIRI.JUNCTION.PASS.SCZ.voom.log.rowMeans[libBatchAll.DEG])
range(libBatch.mean.exp.random)

libBatch.median.exp.random=vector(length = 10000)
for (i in 1:10000) {
  libBatch.median.exp.random[i]=median(sample(CIRI.JUNCTION.PASS.SCZ.voom.log.rowMeans,4143))
}
median(CIRI.JUNCTION.PASS.SCZ.voom.log.rowMeans[libBatchAll.DEG])
range(libBatch.median.exp.random)
{
  pdf("./pdf_C/LIB_ExprandSample.pdf",height=6)
  hist(libBatch.median.exp.random,main="",breaks = 10,
       xlim=c(min(libBatch.median.exp.random),median(CIRI.JUNCTION.PASS.SCZ.voom.log.rowMeans[libBatchAll.DEG])))
  par(new=T)
  plot(median(STA.paper[libBatchAll.DEG,'predicted.length']),0,pch=16,axes=F,ylim=c(0,2600),xlab="",ylab="",
       xlim=c(663,max(libBatch.median.length.random)))
  dev.off()
}


libBatch.median.length.random=vector(length = 10000)
for (i in 1:10000) {
  libBatch.median.length.random[i]=median(sample(STA.paper$predicted.length,4143))
}
median(STA.paper[libBatchAll.DEG,'predicted.length'])
range(libBatch.median.length.random)
{
  pdf("./pdf_C/LIB_randSample.pdf",height=6)
  hist(libBatch.median.length.random,xaxt="n",main="",breaks = 10,
       xlim=c(663,max(libBatch.median.length.random)))
  axis(side=1,at=c(663,seq(750,950,100)),
       labels =c(663,seq(750,950,100)))
  par(new=T)
  plot(median(STA.paper[libBatchAll.DEG,'predicted.length']),0,pch=16,axes=F,ylim=c(0,2600),xlab="",ylab="",
       xlim=c(663,max(libBatch.median.length.random)))
  dev.off()
}

#####################################Institution###########################
InstitutionPenn.DEG=subset(toptable(CIRI.JUNCTION.PASS.SCZ.voom.eBayes,
                              coef="InstitutionPenn",n=Inf,confint=T,sort.by = "p"),adj.P.Val<=0.05)
InstitutionPitt.DEG=subset(toptable(CIRI.JUNCTION.PASS.SCZ.voom.eBayes,
                                    coef="InstitutionPitt",n=Inf,confint=T,sort.by = "p"),adj.P.Val<=0.05)
InstitutionAll.DEG=unique(c(rownames(InstitutionPenn.DEG),rownames(InstitutionPitt.DEG)))

Institution.mean.exp.random=vector(length = 10000)
for (i in 1:10000) {
  Institution.mean.exp.random[i]=mean(sample(CIRI.JUNCTION.PASS.SCZ.voom.log.rowMeans,300))
}
mean(CIRI.JUNCTION.PASS.SCZ.voom.log.rowMeans[InstitutionAll.DEG])
range(Institution.mean.exp.random)

Institution.median.exp.random=vector(length = 10000)
for (i in 1:10000) {
  Institution.median.exp.random[i]=median(sample(CIRI.JUNCTION.PASS.SCZ.voom.log.rowMeans,300))
}
median(CIRI.JUNCTION.PASS.SCZ.voom.log.rowMeans[InstitutionAll.DEG])
range(Institution.median.exp.random)

Institution.median.length.random=vector(length = 10000)
for (i in 1:10000) {
  Institution.median.length.random[i]=median(sample(STA.paper$predicted.length,300))
}
median(STA.paper[InstitutionAll.DEG,'predicted.length'])
range(Institution.median.length.random)

###################################Age of Death##############################
AOD.DEG=toptable(CIRI.JUNCTION.PASS.SCZ.voom.eBayes,
                                         coef="AOD",n=Inf,confint=T,sort.by = "p")
AOD.DEG.sig=subset(AOD.DEG,adj.P.Val<=0.05)
AOD.DEG.sig$Exp=CIRI.JUNCTION.PASS.SCZ.voom.log.rowMeans[rownames(AOD.DEG.sig)]
AOD.DEG.sig$type=STA.paper[rownames(AOD.DEG.sig),'region']

write.table(AOD.DEG.sig,"AOD.DEG.sig.txt",col.names=F,row.names = T,sep="\t",quote=F)

AOD.DEG.circID=rownames(AOD.DEG.sig)

CDR1as=CIRI.JUNCTION.PASS.SCZ.voom.log["chrX:139865340|139866824",]
plot(CDR1as~covariates.SCZ$AOD)
summary(lm(CDR1as~covariates.SCZ$AOD))

{
  pdf("./pdf_C/CDR1as_AOD.pdf",height = 6)
  tmpname="CDR1as"
  cdr1as.df=data.frame(Exp=CDR1as,
                       AOD=covariates.SCZ$AOD)
  tmpcor=cor(cdr1as.df$Exp,cdr1as.df$AOD)
  tmp=ggplot(cdr1as.df,aes(x=AOD,y = Exp))+
    geom_point(size = 2,colour=rgb(215,48,31,maxColorValue = 255,alpha = 200)) + 
    stat_smooth(method="lm",colour=rgb(66,146,198,maxColorValue = 255,alpha = 250))+
    labs(title = tmpname,
         x="Age of death",y="CircRNA expression (log2 CPM)")+ theme_bw(base_size = 12)
  print(tmp)
  dev.off()
}

STA.AOD.DEG=STA.paper[AOD.DEG.circID,]
AOD.DEG.clear=intersect(unique(STA.AOD.DEG$gene),STA.paper.clearGene)
STA.paper.clearGene.geneSymbol=STA.paper[match(STA.paper.clearGene,STA.paper$gene),'geneSymbol']
AOD.DEG.clear.geneSymbol=STA.paper[match(AOD.DEG.clear,STA.paper$gene),'geneSymbol']
writeLines(AOD.DEG.clear,"./enrichment_C/AOD.DEG.gene")

AOD.DEG.POS.circID=rownames(subset(AOD.DEG,adj.P.Val<=0.05 & logFC>0))
AOD.DEG.NEG.circID=rownames(subset(AOD.DEG,adj.P.Val<=0.05 & logFC<0))
AOD.DEG.POS.gene=intersect(unique(STA.paper[AOD.DEG.POS.circID,'gene']),STA.paper.clearGene)
AOD.DEG.NEG.gene=intersect(unique(STA.paper[AOD.DEG.NEG.circID,'gene']),STA.paper.clearGene)
writeLines(AOD.DEG.POS.gene,"./enrichment_C/AOD.DEG.POS.gene")
writeLines(AOD.DEG.NEG.gene,"./enrichment_C/AOD.DEG.NEG.gene")

#
AOD.inTophat2=Tophat2.DGEList.SCZ.voom.AOD[STA.AOD.DEG$gene,]
AOD.inTophat2.merge=data.frame(cbind(AOD.inTophat2,AOD.DEG.sig))
AOD.inTophat2.merge$circID=AOD.DEG.circID
plot(x = -log(AOD.inTophat2.merge$P.Value), y=-log(AOD.inTophat2.merge$P.Value.1),pch=20)
AOD.inTophat2.merge.noSigGene=subset(AOD.inTophat2.merge,P.Value>0.05)

AOD.DEG.sig$parentalSig=0
AOD.DEG.sig[AOD.inTophat2.merge.noSigGene$circID,'parentalSig']=1
#enrichment with gwas loci
GWAS.trait.table=table(GWAS.list$MAPPED_TRAIT)
tmpTrait=names(which(GWAS.trait.table>10))
GWAS.list.sub=subset(GWAS.list,MAPPED_TRAIT%in%tmpTrait)
AOD.DEG.GWAS.df=data.frame(matrix(nrow=length(tmpTrait),ncol=4))
rownames(AOD.DEG.GWAS.df)=tmpTrait
for (i in 1:length(tmpTrait)) {
  tmpDis=tmpTrait[i]
  tmpGWAS.gene=subset(GWAS.list.sub,MAPPED_TRAIT%in%tmpDis)$MAPPED_GENE
  tmpGWAS.gene=unique(gsub("([ ])", "", unlist(strsplit(unlist(strsplit(tmpGWAS.gene,"-")),","))))
  m=length(intersect(tmpGWAS.gene,AOD.DEG.clear.geneSymbol))
  n=length(intersect(tmpGWAS.gene,STA.paper.clearGene.geneSymbol))
  ft=fisher.test(matrix(c(m,length(AOD.DEG.clear)-m,n,length(STA.paper.clearGene)-n),nrow=2,ncol=2))
  AOD.DEG.GWAS.df[i,]=c(m,n,ft$estimate,ft$p.value)
}
rownames(AOD.DEG.GWAS.df)=tmpTrait
AOD.DEG.GWAS.df$FDR=p.adjust(AOD.DEG.GWAS.df$X4,"fdr")
GWAS.disease.table.threshold
#plot AOD fold change with mean circRNA expression
AOD.DEG.color=rep(rgb(40,40,40,200,maxColorValue = 255),nrow(AOD.DEG))
AOD.DEG.color[which(AOD.DEG$adj.P.Val<0.05)]="red"
AOD.plot.df=data.frame(FC=AOD.DEG$logFC,
                       Exp=CIRI.JUNCTION.PASS.SCZ.voom.log.rowMeans[rownames(AOD.DEG)],
                       col=AOD.DEG.color)

{
  pdf("./pdf_C/AOD_FC.pdf",height = 5)
  tmp=ggplot(AOD.plot.df,aes(x=Exp,y = FC))+
    geom_point(size = 2,aes(colour=col))+ scale_color_manual(values = c("#282828C8","red"))+
    labs(x="Mean circRNA expression (log2 CPM)",y="FoldChange")+ theme_bw(base_size = 12)+
    theme(legend.position="none")+geom_hline(aes(yintercept=0),colour=rgb(215,60,41,200,maxColorValue = 255),size=1.5)
  print(tmp)
  dev.off()
}

###################################PMI##############################
PMI.DEG=toptable(CIRI.JUNCTION.PASS.SCZ.voom.eBayes,
                 coef="PMI",n=Inf,confint=T,sort.by = "p")
PMI.DEG.circID=rownames(subset(PMI.DEG,adj.P.Val<=0.05))
STA.PMI.DEG=STA.paper[PMI.DEG.circID,]
PMI.DEG.clear=intersect(unique(STA.PMI.DEG$gene),STA.paper.clearGene)

################################RIN##################################
Gene.RIN.DEG=toptable(Tophat2.DGEList.SCZ.voom.eBays,coef="RIN",n=Inf,confint=T,sort.by = "p")
Gene.RIN2.DEG=toptable(Tophat2.DGEList.SCZ.voom.eBays,coef="RIN2",n=Inf,confint=T,sort.by = "p")
plot(Tophat2.DGEList.SCZ.voom.log[rownames(Gene.RIN.DEG)[9],]~covariates.SCZ$RIN)
abline(lm(Tophat2.DGEList.SCZ.voom.log[rownames(Gene.RIN.DEG)[9],]~covariates.SCZ$RIN))
RIN.DEG=toptable(CIRI.JUNCTION.PASS.SCZ.voom.eBayes,
                 coef="RIN",n=Inf,confint=T,sort.by = "p")
RIN2.DEG=toptable(CIRI.JUNCTION.PASS.SCZ.voom.eBayes,
                  coef="RIN2",n=Inf,confint=T,sort.by = "p")
RIN.DEG.circID=unique(c(rownames(subset(RIN.DEG,adj.P.Val<=0.05)),rownames(subset(RIN2.DEG,adj.P.Val<=0.05))))
STA.RIN.DEG=STA.paper[RIN.DEG.circID,]
#evaluate length by randomly sample 28 exonic circRNAs
RIN.length.median.random=vector(length = 10000)
for (i in 1:10000) {
  RIN.length.median.random[i]=median(sample(STA.paper.exon$predicted.length,28))
}

RIN.median.exp.random=vector(length = 10000)
for (i in 1:10000) {
  RIN.median.exp.random[i]=median(sample(CIRI.JUNCTION.PASS.SCZ.voom.log.rowMeans,28))
}
median(CIRI.JUNCTION.PASS.SCZ.voom.log.rowMeans[RIN.DEG.circID])
range(RIN.median.exp.random)

######################plot RIN random sampling###############
median(STA.RIN.DEG$predicted.length)
{
  pdf("./pdf_C/RIN_randSample.pdf",height=5)
  hist(RIN.length.median.random,xaxt="n",main="",
       xlim=c(min(RIN.length.median.random),median(STA.RIN.DEG$predicted.length)))
  axis(side=1,at=c(seq(400,2000,400),median(STA.RIN.DEG$predicted.length)),
       labels =c(seq(400,2000,400),median(STA.RIN.DEG$predicted.length)))
  par(new=T)
  plot(2279,0,pch=16,axes=F,ylim=c(0,2600),xlab="",ylab="",
       xlim=c(min(RIN.length.median.random),median(STA.RIN.DEG$predicted.length)))
  dev.off()
}


#
RIN.DEG.color=rep(rgb(40,40,40,200,maxColorValue = 255),nrow(RIN.DEG))
RIN.DEG.color[which(RIN.DEG$adj.P.Val<0.05)]="red"
# plot(RIN.DEG$logFC~CIRI.JUNCTION.PASS.SCZ.voom.log.rowMeans[rownames(RIN.DEG)],
#      xlab="Mean log2(CPM)",ylab="Log2FC",bty="l",
#      pch=20,cex=0.6,col=RIN.DEG.color,las=1)
CIRI.JUNCTION.PASS.SCZ.voom.log.keepRIN=calcResiduals(geneBySampleValues = CIRI.JUNCTION.PASS.SCZ.voom.log,
                                                      samplesByCovariates = Circ.design.SCZwithControl,
                                                      varsToAddBackIn = c("RIN","RIN2"),
                                                      sampleWeights = CIRI.JUNCTION.PASS.SCZ.voom.weight)
for (i in 1:28) {
  RIN.df=data.frame(Exp=as.numeric(CIRI.JUNCTION.PASS.SCZ.voom.log[RIN.DEG.circID[i],]),
                    RIN=covariates.SCZ$RIN)
  tmp=ggplot(RIN.df,aes(x=RIN,y = Exp))+
    geom_point(size = 2,colour=rgb(215,48,31,maxColorValue = 255,alpha = 200)) + 
    stat_smooth(method="loess",colour=rgb(66,146,198,maxColorValue = 255,alpha = 250))+
    labs(title =paste(PGC.scz2.min.QTL[i,"geneSymbol"] ," ",RIN.DEG.circID[i]),
         x="RIN",y="CircRNA expression (log2 CPM)")+
    theme_pubr()
  print(tmp)
}
{
  pdf("./pdf_C/RIN_example.pdf",height = 5)
  tmpname="chr10:94653106|94733989"
  RIN.df=data.frame(Exp=as.numeric(CIRI.JUNCTION.PASS.SCZ.voom.log[tmpname,]),
                    RIN=covariates.SCZ$RIN)
  tmpcor=cor(RIN.df$Exp,RIN.df$RIN)
  tmp=ggplot(RIN.df,aes(x=RIN,y = Exp))+
    geom_point(size = 2,colour=rgb(215,48,31,maxColorValue = 255,alpha = 200)) + 
    stat_smooth(method="loess",colour=rgb(66,146,198,maxColorValue = 255,alpha = 250))+
    labs(title =paste(STA.paper[tmpname,"geneSymbol"] ," ",tmpname),
         x="RIN",y="CircRNA expression (log2 CPM)")+ theme_bw(base_size = 12)
  print(tmp)
  dev.off()
}
STA.paper[tmpname,]


abline(lm(as.numeric(CIRI.JUNCTION.PASS.SCZ.voom.log["chr1:240370099|240421332",])~covariates.SCZ$RIN))
plot(as.numeric(CIRI.JUNCTION.PASS.SCZ.voom.log["chr18:10759480|10773627",])~covariates.SCZ$RIN)
abline(lm(as.numeric(CIRI.JUNCTION.PASS.SCZ.voom.log["chr18:10759480|10773627",])~covariates.SCZ$RIN))

#########################Gender###########################
Gender.DEG=toptable(CIRI.JUNCTION.PASS.SCZ.voom.eBayes,
                 coef="SexMale",n=Inf,confint=T,sort.by = "p")
Gender.DEG.sig005=subset(Gender.DEG,adj.P.Val<0.05)
idx=union(grep("chrX",rownames(Gender.DEG)),grep("chrY",rownames(Gender.DEG)))
Gender.DEG.noSexChr=Gender.DEG[-idx,]
Gender.DEG.noSexChr$adj.P.Val=p.adjust(Gender.DEG.noSexChr$P.Value,"fdr")
plot(as.numeric(CIRI.JUNCTION.PASS.SCZ.voom.log["chr9:44106654|44114921",])~as.factor(covariates.SCZ$Sex))
STA.paper["chr9:44106654|44114921",]
gender.df=data.frame(value=as.numeric(CIRI.JUNCTION.PASS.SCZ.voom.log["chr9:44106654|44114921",]),
                     group=as.factor(covariates.SCZ$Sex))
circname="chr9:44106654|44114921"
genename=STA.paper[circname,'geneSymbol']
draw.ggviolin=ggviolin(gender.df,x="group",y="value",fill="group",
                       title =paste(genename,circname),cex=0.8,
                       xlab = "",
                       ylab="",
                       panel.labs = F,
                       size=0.1,
                       palette = c("#FC4E07","#00AFBB"),width = 0.8,
                       add="mean",add.params = list(color="black",size=0.1))+
  stat_compare_means(label.y=max(violin.df$value)+2*sqrt(var(violin.df$value)),label.x = 0.8)+
  theme(plot.title = element_text(hjust = 0.5))+
  font("title",size=9)+ theme(legend.position="none")+
  theme(legend.position="none")+ 
  theme_bw(base_size = 12)
print(draw.ggviolin)
{
  pdf("pdf_C/Gender_example.pdf",width=6)
  draw.ggbox=ggplot(gender.df,aes(x=group,y=value))+geom_boxplot(aes(fill=group))+
    theme_bw(base_size = 12)+theme(legend.position="none")+labs(title =paste(STA.paper[circname,"geneSymbol"] ," ",circname))
  print(draw.ggbox)
  dev.off()
}

##########################Residuals after remove covariates##############
CIRI.JUNCTION.PASS.voom.log.Res=calcResiduals(CIRI.JUNCTION.PASS.voom.log,Circ.design.all,sampleWeights = CIRI.JUNCTION.PASS.voom.weight)
write.table(CIRI.JUNCTION.PASS.voom.log.Res[,CaucasianID],"./qtl_C/CIRI.JUNCTION.PASS.voom.log.Res.CaucasianID.txt",row.names = T,col.names = T,sep="\t" ,quote=F)
#######################Focus2 gene################
########################Relation between Junctions Reads and Linear Reads###########
CIRI.JUNCTION.AND.Linear.counts.P=vector(length = nrow(STA.paper))
CIRI.JUNCTION.AND.Linear.counts.R=vector(length = nrow(STA.paper))
for (i in 1:nrow(STA.paper)) {
  idx=unique(c(which(is.na(CIRI.JUNCTION.PASS.mat[i,])),which(is.na(CIRI.Linear.mat[i,]))))
  if(length(idx)!=0){
    x=as.numeric(CIRI.JUNCTION.PASS.mat[i,-idx])
    y=as.numeric(CIRI.Linear.mat[i,-idx])
  }else{
    x=as.numeric(CIRI.JUNCTION.PASS.mat[i,])
    y=as.numeric(CIRI.Linear.mat[i,])
  }
  
  n=length(x)
  r=sum(scale(x)*scale(y))/(n-1)
  tmpSta=summary(lm(x~y))
  CIRI.JUNCTION.AND.Linear.counts.P[i]=tmpSta$coefficients[2,4]
  CIRI.JUNCTION.AND.Linear.counts.R[i]=r
}

names(CIRI.JUNCTION.AND.Linear.counts.P)=STA.paper$circRNA.ID
names(CIRI.JUNCTION.AND.Linear.counts.R)=STA.paper$circRNA.ID
CIRI.JUNCTION.AND.Linear.counts.FDR=p.adjust(CIRI.JUNCTION.AND.Linear.counts.P,"fdr")
sum(CIRI.JUNCTION.AND.Linear.counts.FDR<0.05)
idx=which(CIRI.JUNCTION.AND.Linear.counts.FDR<0.05)
Isoform.Cor.STA.Sig005=STA.paper[names(idx),]
CIRI.JUNCTION.AND.Linear.counts.sig.R=CIRI.JUNCTION.AND.Linear.counts.R[idx]
idx.neg=names(which(CIRI.JUNCTION.AND.Linear.counts.sig.R<0))
CIRI.JUNCTION.AND.Linear.counts.FDR[idx.neg]
for (i in 1:length(idx.neg)) {
  i=idx.neg[i]
  idx=unique(c(which(is.na(CIRI.JUNCTION.PASS.mat[i,])),which(is.na(CIRI.Linear.mat[i,]))))
  x=as.numeric(CIRI.JUNCTION.PASS.mat[i,-idx])
  y=as.numeric(CIRI.Linear.mat[i,-idx])
  tmp= "CircRNA:"%&%median(x)%&%" LinearRNA:"%&%median(y)
  plot(y~x,pch=20,main=tmp)
  abline(lm(y~x))
}
Cor.df=data.frame(R=CIRI.JUNCTION.AND.Linear.counts.R,type=STA.paper$region)
{
  pdf("./pdf_C/Isoforms_correlation.pdf",height = 5)
  tmp=ggplot(Cor.df)+
    geom_histogram(aes(R,fill=type,color=type),binwidth = 0.1)+
    theme_bw(base_size = 12)+
    scale_fill_manual(values = ggplot2::alpha(c("#F8766D","#00BA38","#619CFF"),0.7))
  print(tmp)
  dev.off()
}
idxnoSig005=which(CIRI.JUNCTION.AND.Linear.counts.FDR>=0.05)
Isoform.Cor.STA.noSig005=STA.paper[names(idxnoSig005),]
tmp=matrix(c(table(Isoform.Cor.STA.noSig005$region),
table(Isoform.Cor.STA.Sig005$region)),nrow=2,byrow = T)
type.chisq=chisq.test(tmp)

#####################################Plot Fig2D#####################
#choice chr15:65440557-65477680 for chr15:65471272|65472542
#PENN_RNA_PFC_45
#PENN_RNA_PFC_59
#PITT_RNA_PFC_1270
#PENN_RNA_PFC_28
#PITT_RNA_PFC_10003
#samtools depth -r chr15:65440557-65477680 -a /media/data3/CMC/data/BAM/aligned/PENN_RNA_PFC_45.accepted_hits.sort.coord.bam >/media/data3/circCMC/data/PENN_RNA_PFC_45.neg.pos
#samtools depth -r chr15:65440557-65477680 -a /media/data3/CMC/data/BAM/aligned/PENN_RNA_PFC_59.accepted_hits.sort.coord.bam >/media/data3/circCMC/data/PENN_RNA_PFC_59.neg.pos
###########################################################################
# choice.plot.sample=c("PENN_RNA_PFC_45","PENN_RNA_PFC_59","PITT_RNA_PFC_1270","PENN_RNA_PFC_28","PITT_RNA_PFC_10003")
# 
# neg.depth1=read.delim("./PENN_RNA_PFC_45.neg.pos",header=F)
# neg.depth2=read.delim("./PENN_RNA_PFC_59.neg.pos",header=F)
# neg.depth1.total=sum(neg.depth1$V3)
# neg.depth2.total=sum(neg.depth2$V3)
# neg.depth1$density=neg.depth1$V3/neg.depth1.total
# neg.depth2$density=neg.depth2$V3/neg.depth2.total
# par(mfrow=c(2,1))
# plot(neg.depth1$density~neg.depth1$V2,type="l",ylim=c(0,8*10^(-4)),col="red",xlim=c(65471272,65472542))
# par(new=T)
# plot(neg.depth2$density~neg.depth2$V2,type="l",ylim=c(0,8*10^(-4)),col="blue",xlim=c(65471272,65472542))
# Gene.sig=subset(Covarites.EVAL.FDR.mat,gene<0.05)
# Gene.sig=Gene.sig[order(Gene.sig$gene),]
# Gene.sig.STA=STA.paper[match(rownames(Gene.sig),STA.paper$circRNA.ID),]
# plot(as.numeric(CIRI.JUNCTION.PASS.SCZ.NA.mat[Gene.sig.STA$circRNA.ID[1],]),
#      as.numeric(Tophat2.SCZ.DGEList.sub.cpm[Gene.sig.STA$gene[1],]))
# Gene.sig.STA.exon=subset(Gene.sig.STA,region=="exon")
# Gene.nosig05=subset(Covarites.EVAL.P,gene>0.5)
# Gene.nosig05=Gene.nosig05[order(Gene.nosig05$gene),]
# Gene.nosig05.STA=STA.paper[match(rownames(Gene.nosig05),STA.paper$circRNA.ID),]
# Gene.nosig05.STA.exon=subset(Gene.nosig05.STA,region=="exon")
# sum(Covarites.EVAL.effect[,'gene']>0,na.rm = T)
# sum(Covarites.EVAL.effect[,'gene']<0,na.rm = T)
# #effect is negative
# Covarites.EVAL.effect.sig=Covarites.EVAL.effect[Gene.sig.STA$circRNA.ID,'gene']
# which(Covarites.EVAL.effect.sig<0)
# Gene.sig.neg=Gene.sig[which(Covarites.EVAL.effect.sig<0),]
# Gene.sig.neg.STA=Gene.sig.STA[which(Covarites.EVAL.effect.sig<0),]
# Gene.sig.pos.STA=Gene.sig.STA[which(Covarites.EVAL.effect.sig>0),]
# rownames(Gene.sig.STA)=Gene.sig.STA$circRNA.ID

###########################Cal Cor for Multiple Isoforms in a same gene###############
calCor=function(tmpmat){
  all_cor=corr.test(tmpmat, use='pairwise.complete.obs', method="pearson", adjust="none")
  all_cor_vals = all_cor$r
  all_cor_p = pt(-abs(all_cor$t),all_cor$n-1)*2
  all_cor_p[lower.tri(all_cor_p)]=0
  cor_mat = melt(all_cor_p, varnames=c("COMPARE", "COVAR"))
  colnames(cor_mat)[colnames(cor_mat) == "value"] = "pvalue"
  cor_mat$r = melt(all_cor_vals)$value
  return(subset(cor_mat,pvalue>0))
}

Isoform.Cor.list=list()

STA.paper.clearGene=unique(STA.paper.clear[,'gene'])
for (i in 1:length(STA.paper.clearGene)) {
  tmpgene=STA.paper.clearGene[i]
  tmplinear=intersect(gtfV19.transcript.detail[which(gtfV19.transcript.detail$geneID==tmpgene),'transID'],rownames(Isoform.log.cpm))
  if(length(tmplinear)>0){
    tmpcirc=STA.paper.clear[which(STA.paper.clear$gene==tmpgene),'circRNA.ID']
    tmpmat=t(as.matrix(rbind(Isoform.log.cpm[tmplinear,,drop=F],CIRI.JUNCTION.PASS.voom.log.Res[tmpcirc,,drop=F])))
    Isoform.Cor.list[[tmpgene]]=calCor(tmpmat)
  }
}
Isoform.Cor.df=as.data.frame(do.call("rbind",Isoform.Cor.list))
Isoform.Cor.df$fdr=p.adjust(Isoform.Cor.df$pvalue,"fdr")
#type classification 0: circRNA-linear 1: circRNA-circRNA 2: linear-linear
Isoform.Cor.df$type=0
idx1=which(Isoform.Cor.df$COMPARE%in%STA.paper.clear$circRNA.ID)
idx2=which(Isoform.Cor.df$COVAR%in%STA.paper.clear$circRNA.ID)
Isoform.Cor.df$type[intersect(idx1,idx2)]=1
Isoform.Cor.df$type[-union(idx1,idx2)]=2
Isoform.Cor.df$gene=substr(rownames(Isoform.Cor.df),1,15)
Isoform.Cor.df$COMPARE=as.character(Isoform.Cor.df$COMPARE)
Isoform.Cor.df$COVAR=as.character(Isoform.Cor.df$COVAR)

Isoform.Cor.df.type0=subset(Isoform.Cor.df,type==0)
Isoform.Cor.df.type0.sig005=subset(Isoform.Cor.df.type0,fdr<0.05)
nrow(Isoform.Cor.df.type0)
nrow(Isoform.Cor.df.type0.sig005)

Isoform.Cor.df.sig05=subset(Isoform.Cor.df,fdr<0.05)

#type 1
Isoform.Cor.df.type1=subset(Isoform.Cor.df,type==1)
Isoform.Cor.df.type1.sig005=subset(Isoform.Cor.df.type1,fdr<0.05)
nrow(Isoform.Cor.df.type1)
nrow(Isoform.Cor.df.type1.sig005)
Isoform.Cor.df.type1.pos=cbind(STA.pos[Isoform.Cor.df.type1$COMPARE,c("left","right")],
                               STA.pos[Isoform.Cor.df.type1$COVAR,c("left","right")])
colnames(Isoform.Cor.df.type1.pos)=c("l1","r1","l2","r2")
Isoform.Cor.df.type1.pos=data.frame(Isoform.Cor.df.type1.pos)
Isoform.Cor.df.type1.pos=cbind(Isoform.Cor.df.type1,Isoform.Cor.df.type1.pos)
Isoform.Cor.df.type1.pos$dl=Isoform.Cor.df.type1.pos$l1-Isoform.Cor.df.type1.pos$l2
Isoform.Cor.df.type1.pos$dr=Isoform.Cor.df.type1.pos$r1-Isoform.Cor.df.type1.pos$r2
#0: no overlap 1: one same site
Isoform.Cor.df.type1.pos$d=apply(Isoform.Cor.df.type1.pos[,c("dl","dr")],1,function(x){ifelse(min(abs(x))==0,1,0)})
Isoform.Cor.df.type1.pos$length1=STA.paper[Isoform.Cor.df.type1.pos$COMPARE,'predicted.length']
Isoform.Cor.df.type1.pos$length2=STA.paper[Isoform.Cor.df.type1.pos$COVAR,'predicted.length']
Isoform.Cor.df.type1.pos$exp1=CIRI.JUNCTION.PASS.SCZ.voom.log.rowMeans[Isoform.Cor.df.type1.pos$COMPARE]
Isoform.Cor.df.type1.pos$exp2=CIRI.JUNCTION.PASS.SCZ.voom.log.rowMeans[Isoform.Cor.df.type1.pos$COVAR]
Isoform.Cor.df.type1.pos$LEsig=sign(Isoform.Cor.df.type1.pos$length1-Isoform.Cor.df.type1.pos$length2)*sign(Isoform.Cor.df.type1.pos$exp1-Isoform.Cor.df.type1.pos$exp2)

Isoform.Cor.df.type1.pos.sig005=subset(Isoform.Cor.df.type1.pos,fdr<0.05)
Isoform.Cor.df.type1.pos.sig005.pos=subset(Isoform.Cor.df.type1.pos.sig005,r>0)
Isoform.Cor.df.type1.pos.sig005.neg=subset(Isoform.Cor.df.type1.pos.sig005,r<0)

table(Isoform.Cor.df.type1.pos.sig005$d)
table(Isoform.Cor.df.type1.pos.sig005.pos$d)
table(Isoform.Cor.df.type1.pos.sig005.neg$d)

hist(Isoform.Cor.df.type1.pos.sig005$r)
Isoform.Cor.df.type1.pos.nosig005=subset(Isoform.Cor.df.type1.pos,fdr>0.05)
table(Isoform.Cor.df.type1.pos.nosig005$d)

{
  pdf("./pdf_C/Isoform_circ_circ.pdf",width = 5)
  temp=cbind(table(Isoform.Cor.df.type1.pos.sig005.pos$d),
             table(Isoform.Cor.df.type1.pos.sig005.neg$d),
             table(Isoform.Cor.df.type1.pos.nosig005$d)
  )
  colnames(temp)=c("Positive","Negative","Non")
  bar=barplot(prop.table(temp,margin = 2)*100,ylab="Percentage",ylim=c(0,100),cex.axis=1.8,cex.names = 1.7)
  text(bar,80,paste(apply(temp,2,sum),sep=" "),cex=1.6)
  dev.off()
}

#type: 0
Isoform.Cor.df.type0=subset(Isoform.Cor.df,type==0)
hist(Isoform.Cor.df.type0$r)

idx=names(which(table(Isoform.Cor.df.type0$gene)>1))
Isoform.Cor.df.type0.single=subset(Isoform.Cor.df.type0,!(gene%in%idx))

plot(Isoform.Cor.df.type0.single$r~CIRI.JUNCTION.AND.Linear.counts.R[Isoform.Cor.df.type0.single$COVAR])


sum(table(Isoform.Cor.df.type0.sig005$gene)>1)

idx=names(which(table(Isoform.Cor.df.type0.sig005$gene)>1))
Isoform.Cor.df.type0.sig005.multi=subset(Isoform.Cor.df.type0.sig005,gene%in%idx)
Isoform.Cor.df.type0.sig005.single=subset(Isoform.Cor.df.type0.sig005,!(gene%in%idx))

singel.df=data.frame(isoR=Isoform.Cor.df.type0.single$r,
                     readR=CIRI.JUNCTION.AND.Linear.counts.R[Isoform.Cor.df.type0.single$COVAR])
{
  pdf("./pdf_C/Isoform_pearson_R.pdf")
  tmp=ggplot(singel.df,aes(x=isoR,y = readR))+
    geom_point(size = 3,colour=rgb(215,48,31,maxColorValue = 255,alpha = 200)) + 
    stat_smooth(method="lm",colour=rgb(66,146,198,maxColorValue = 255,alpha = 250))+
    labs(x="Pearson's r of circRNA and linRNA expression",y="Pearson's r of circular and linear junction reads")+theme_bw(base_size = 12)
  print(tmp)
  dev.off()
}
Isoform.Cor.df.type0.sig005.pos=subset(Isoform.Cor.df.type0.sig005,r>0)
Isoform.Cor.df.type0.sig005.neg=subset(Isoform.Cor.df.type0.sig005,r<0)

idx=which(CIRI.JUNCTION.AND.Linear.counts.R[Isoform.Cor.df.type0.sig005.neg$COVAR]<0)
Isoform.Cor.df.type0.sig005.neg.choice=Isoform.Cor.df.type0.sig005.neg[idx,]
Isoform.Cor.df.type0.sig005.neg.choice.multi=subset(Isoform.Cor.df.type0.sig005.neg.choice,gene%in%unique(Isoform.Cor.df.type0.sig005.multi$gene))
Isoform.Cor.df.type0.sig005.neg.choice.multi$readR=CIRI.JUNCTION.AND.Linear.counts.R[Isoform.Cor.df.type0.sig005.neg.choice.multi$COVAR]
Isoform.Cor.df.type0.sig005.neg.choice.multi$length=STA.paper[Isoform.Cor.df.type0.sig005.neg.choice.multi$COVAR,
                                                              "predicted.length"]
subset(Isoform.Cor.df.type0,COVAR=="chr14:36121011|36194359")
STA.paper["chr14:36121011|36194359",]
unGene=Isoform.Cor.df.type0.sig005.neg.choice.multi$gene #first time
unGene=unique(Isoform.Cor.df.type0.sig005.neg$gene)
for (i in 1:length(unGene)) {
  geneID=unGene[i]
  geneID.transcript=gtfV19.transcript.detail[which(gtfV19.transcript.detail$geneID==geneID),]
  geneID.detail=subset(gtfV19.gene,gene==geneID)
  region=geneID.detail[,1]%&%":"%&%geneID.detail[,4]%&%"-"%&%geneID.detail[,5]
  samplePool=c("PENN_RNA_PFC_45","PENN_RNA_PFC_59","PITT_RNA_PFC_1270","PENN_RNA_PFC_28","PITT_RNA_PFC_10003")
  for (j in samplePool) {
    baseDir="/media/data3/CMC/data/fastq/"
    sample=baseDir%&%j%&%".bam"
    out="/media/data3/circCMC/data/isoform/"%&%j%&%"."%&%geneID%&%".depth"
    command=paste("samtools depth -a -r",region,sample,">",out,sep=" ")
    try(system(command))
  }
}
genePlot=function(gene,circID,sample=samplePool,CircRegion=F){
  mycol=colorRampPalette(rev(brewer.pal(n = 7, name = "Set2")))(length(sample))
  tmp=unlist(strsplit(circID,":"))
  tmp2=as.numeric(unlist(strsplit(tmp[2],"\\|")))
  for (i in 1:length(sample)) {
    depth=read.delim("./isoform/"%&%sample[i]%&%"."%&%gene%&%".depth",header=F)
    depth$density=depth[,3]/sum(depth[,3])
    if(CircRegion==T){
      idx=seq((tmp2[1]-1000-depth[1,2]+1),(tmp2[2]+1000-depth[1,2])+1)
      depth=depth[idx,]
    }
    if(i==1){
      plot(density~V2,depth,xlab=gene,ylab="density",main=circID,las=1,cex.axis=0.8,col=mycol[i],
           type="l")
      abline(v=tmp2[1],col="red")
      abline(v=tmp2[2],col="red")
    }else{
      par(new=T)
      plot(density~V2,depth,axes=FALSE,xlab="",ylab="",las=1,cex.axis=0.8,col=mycol[i],
           type="l")
    }
  }
}  
which(Isoform.Cor.df.type0.sig005.neg$COVAR=="chr15:65471272|65472542")

gene=Isoform.Cor.df.type0.sig005.neg.choice.multi[1,"gene"]
circID=Isoform.Cor.df.type0.sig005.neg.choice.multi[1,2]

gene=Isoform.Cor.df.type0.sig005.neg[726,'gene']
circID=Isoform.Cor.df.type0.sig005.neg[726,2]

circGroup=Isoform.Cor.df.type0.sig005.neg.choice.multi[
  !duplicated(Isoform.Cor.df.type0.sig005.neg.choice.multi$COVAR),]
genePlot(gene,
         circID,CircRegion = T)
for (i in 1:nrow(circGroup)) {
  genePlot(circGroup[i,"gene"],
           circGroup[i,2],CircRegion = T)
}



##########################Remove gene and eQTL covariates###############
#read Gaucasian sample ID
CaucasianID=readLines("/media/data3/CMC/data/genotype/Caucasian465/Caucasian.id")
length(intersect(subset(covariates,Dx=="Control")[,'id'],CaucasianID))
length(intersect(subset(covariates,Dx=="SCZ")[,'id'],CaucasianID))
length(intersect(subset(covariates,Dx=="AFF")[,'id'],CaucasianID))


#############################DAVID analysis for parental gene of circRNA####################
#use DAVID on web
#combine it all
#read DAVID results
DAVID.RE=read.delim("./enrichment_C/DAVID.txt",header=F)
tmp=read.delim("./enrichment_C/COG_ONTOLOGY.txt")
colnames(DAVID.RE)=colnames(tmp)
rm(tmp)
DAVID.RE.sig=subset(DAVID.RE,Benjamini<0.05)
DAVID.RE.sig=DAVID.RE.sig[order(DAVID.RE.sig$Benjamini),]
DAVID.RE.sig.sub=DAVID.RE.sig[,c("Category","Term","Count","PValue","Fold.Enrichment","Benjamini")]
write.table(DAVID.RE.sig.sub,"./enrichment_C/DAVID.RE.sig.txt",row.names = F,col.names = T,sep="\t",quote=F)
###############################Imputation circ-to-linear ratio##########
CIRI.JUNCTION.PASS.mat[1,]
CIRI.Linear.mat[1,]
CIRI.RATIO.mat=2*CIRI.JUNCTION.PASS.mat/(2*CIRI.JUNCTION.PASS.mat+CIRI.Linear.mat)
# CIRI.RATIO.mat.t=t(CIRI.RATIO.mat)
# cl <- makeCluster(30)
# registerDoParallel(cl)
# set.seed(2018)
# CIRI.RATIO.IM.t=missForest(CIRI.RATIO.mat.t,verbose = T,parallelize = "variables")
# CIRI.RATIO.IM=t(CIRI.RATIO.IM.t$ximp)
# write.table(CIRI.RATIO.IM,"./qtl_C/CIRI.RATIO.IM.txt",col.names = T,row.names = T,sep="\t",quote=F)
# CIRI.RATIO.IM.Res=calcResiduals(CIRI.RATIO.IM,samplesByCovariates = Circ.design.all)
# CIRI.RATIO.IM.Res.Caucasian465=CIRI.RATIO.IM.Res[,CaucasianID]
# write.table(CIRI.RATIO.IM,"./qtl_C/CIRI.RATIO.IM.txt",col.names = T,row.names = T,sep="\t",quote=F)
# write.table(CIRI.RATIO.IM.Res.Caucasian465,
#             "./qtl_C/CIRI.RATIO.IM.Res.Cau465.txt",col.names = T,row.names = T,sep="\t",quote=F)
############################no imputation circRNA ratio##############
CIRI.RATIO.mat.zero=CIRI.RATIO.mat
CIRI.RATIO.mat.zero[is.na(CIRI.RATIO.mat.zero)]=0
CIRI.RATIO.mat.zero.Res=calcResiduals(CIRI.RATIO.mat.zero,samplesByCovariates = Circ.design.all)
CIRI.RATIO.mat.zero.Res.Caucasian465=CIRI.RATIO.mat.zero.Res[,CaucasianID]
write.table(CIRI.RATIO.mat.zero.Res.Caucasian465,"./qtl_C/CIRI.RATIO.mat.zero.Res.Caucasian465.txt",
            col.names = T,row.names = T,sep="\t",quote=F)

#############################QTL part############################
library(data.table)
library(qvalue)
library(ggplot2)
library(qqman)
############################################################
#nohup ./circCMCQTL_GeneKeep.sh &
#cd /media/data3/circCMC/data/qtl_C/GeneKeep/resultCis/
#for i in {1..22};do tail -n +2 cis.chr${i}.txt;done >cis.chrALL.txt
###########################################################
QTL.CIS.GeneKeep=read.delim("./qtl_C/GeneKeep/resultCis/cis.chrALL.txt",header=F)
colnames(QTL.CIS.GeneKeep)=c("SNP",     "gene",    "beta",    "tstat",  "pvalue", "FDR")
QTL.CIS.GeneKeep.pos=data.frame(do.call("rbind",strsplit(QTL.CIS.GeneKeep$SNP,"_")))
QTL.CIS.GeneKeep.key=paste(QTL.CIS.GeneKeep$X1,QTL.CIS.GeneKeep.pos$X2,QTL.CIS.GeneKeep.pos$X3,sep="_")
write.table(QTL.CIS.GeneKeep.key,"./qtl_C/QTL.CIS.GeneKeep.key",row.names = F,col.names = F,quote = F)

####################################################################
#nohup ./circCMCQTL_RATIO_Zero.sh &
#cd /media/data3/circCMC/data/qtl_C/RatioZero/resultCis/
#for i in {1..22};do tail -n +2 cis.chr${i}.txt;done >cis.chrALL.txt
###################################################################
QTL.CIS.Ratio=read.delim("./qtl_C/RatioZero/resultCis/cis.chrALL.txt",header=F)
colnames(QTL.CIS.Ratio)=c("SNP",     "gene",    "beta",    "tstat",  "pvalue", "FDR")
QTL.CIS.Ratio.pos=data.frame(do.call("rbind",strsplit(QTL.CIS.Ratio$SNP,"_")))
QTL.CIS.Ratio.key=paste(QTL.CIS.Ratio$X1,QTL.CIS.Ratio.pos$X2,QTL.CIS.Ratio.pos$X3,sep="_")
write.table(QTL.CIS.Ratio.key,"./qtl_C/QTL.CIS.Ratio.key",row.names = F,col.names = F,quote = F)

#############get minimum P-value##################
#cd /media/data3/circCMC/data/qtl_C/GeneKeep/minP/
# for i in {1..22};do cat minpvgenesnp.chr$i.txt ;done >minpvgenesnp.chrALL.txt
#cd /media/data3/circCMC/data/qtl_C/RatioZero/minP/
# for i in {1..22};do cat minpvgenesnp.chr$i.txt;done >minpvgenesnp.chrALL.txt
#read minimum P-value
QTL.MINP.GeneKeep=read.delim("./qtl_C/GeneKeep/minP/minpvgenesnp.chrALL.txt",header = F)
colnames(QTL.MINP.GeneKeep)=c("circID","pvalue")
QTL.MINP.GeneKeep=na.omit(QTL.MINP.GeneKeep)
QTL.MINP.GeneKeep=subset(QTL.MINP.GeneKeep,pvalue<1)
rownames(QTL.MINP.GeneKeep)=QTL.MINP.GeneKeep$circID
#############
QTL.MINP.Ratio=read.delim("./qtl_C/RatioZero/minP/minpvgenesnp.chrALL.txt",header = F)
colnames(QTL.MINP.Ratio)=c("circID","pvalue")
QTL.MINP.Ratio=na.omit(QTL.MINP.Ratio)
QTL.MINP.Ratio=subset(QTL.MINP.Ratio,pvalue<1)
rownames(QTL.MINP.Ratio)=QTL.MINP.Ratio$circID

#############permutation results###############
QTL.PERM.GeneKeep=read.delim("./qtl_C/GeneKeep/permutation/last/combineALL.txt",header=F,row.names = 1)
QTL.PERM.GeneKeep=QTL.PERM.GeneKeep[match(QTL.MINP.GeneKeep$circID,rownames(QTL.PERM.GeneKeep)),]
QTL.PERM.GeneKeep.Result=rowSums(QTL.PERM.GeneKeep<QTL.MINP.GeneKeep$pvalue,na.rm = T)
QTL.PERM.GeneKeep.Empirical=(QTL.PERM.GeneKeep.Result+1)/(ncol(QTL.PERM.GeneKeep)+1)
QTL.PERM.GeneKeep.Empirical.qvalue=qvalue(QTL.PERM.GeneKeep.Empirical,fdr.level = 0.05,pfdr = T)
{
  png(filename="./qtl_C/GeneKeep/QTL.PERM.png",width=480,height = 540,units="px")
  par(family="Arial",font.lab=1,mfrow=c(1,2),mar=c(5,5,1,0.5))
  hist(QTL.PERM.GeneKeep.Empirical,col = "grey",xlab="Empirical P-values",main="",axes = F)
  axis(side=1,at=seq(0,1,0.2),labels = seq(0,1,0.2),pos=0,las=1)
  axis(side=2,at=seq(0,3000,500),labels = seq(0,3000,500),pos=0,las=1)
  hist(QTL.PERM.GeneKeep.Empirical.qvalue$qvalues,col = "grey",xlab="Q-values",main="",axes=F)
  axis(side=1,at=seq(0,0.5,0.1),labels = seq(0,0.5,0.1),pos=0,las=1)
  axis(side=2,at=seq(0,3000,500),labels = seq(0,3000,500),pos=0,las=1)
  #abline(v = 0.05,col="blue",lty=5)
  dev.off()
}
############################
QTL.PERM.Ratio=read.delim("./qtl_C/RatioZero/permutation/last/combineALL.txt",header=F,row.names = 1)
QTL.PERM.Ratio=QTL.PERM.Ratio[match(QTL.MINP.Ratio$circID,rownames(QTL.PERM.Ratio)),]
QTL.PERM.Ratio.Result=rowSums(QTL.PERM.Ratio<QTL.MINP.Ratio$pvalue,na.rm = T)
QTL.PERM.Ratio.Empirical=(QTL.PERM.Ratio.Result+1)/(ncol(QTL.PERM.Ratio)+1)
QTL.PERM.Ratio.Empirical.qvalue=qvalue(QTL.PERM.Ratio.Empirical,fdr.level = 0.05,pfdr = T)
{
  png(filename="./qtl_C/RatioZero/QTL.PERM.png",width=480,height = 540,units="px")
  par(family="Arial",font.lab=1,mfrow=c(1,2),mar=c(5,5,1,0.5))
  hist(QTL.PERM.Ratio.Empirical,col = "grey",xlab="Empirical P-values",main="",axes = F)
  axis(side=1,at=seq(0,1,0.2),labels = seq(0,1,0.2),pos=0,las=1)
  axis(side=2,at=seq(0,3000,500),labels = seq(0,3000,500),pos=0,las=1)
  hist(QTL.PERM.Ratio.Empirical.qvalue$qvalues,col = "grey",xlab="Q-values",main="",axes=F)
  axis(side=1,at=seq(0,0.5,0.1),labels = seq(0,0.5,0.1),pos=0,las=1)
  axis(side=2,at=seq(0,3000,500),labels = seq(0,3000,500),pos=0,las=1)
  #abline(v = 0.05,col="blue",lty=5)
  dev.off()
}

####################get ecircQTL##################
QTL.CIS.GeneKeep.Adj=subset(QTL.CIS.GeneKeep,gene %in% names(which(QTL.PERM.GeneKeep.Empirical.qvalue$qvalues<0.05)))
QTL.EmpiricalThreshold.GeneKeep=unique(
  as.numeric(QTL.PERM.GeneKeep.Empirical.qvalue$pvalues[
    which(QTL.PERM.GeneKeep.Empirical.qvalue$qvalues==
            max(QTL.PERM.GeneKeep.Empirical.qvalue$qvalues[
              which(QTL.PERM.GeneKeep.Empirical.qvalue$qvalues<0.05)]))]))
##The QTL.CIS.PERM is too large, so using block to calculate the matrix 
each=10000
QTL.CIS.PERM.GeneKeep.Result=vector(length = nrow(QTL.CIS.GeneKeep.Adj))
for (i in seq(1,nrow(QTL.CIS.GeneKeep.Adj),by = each)) {
  print(i)
  QTL.CIS.PERM.GeneKeep.BLOCK=QTL.PERM.GeneKeep[
    match(QTL.CIS.GeneKeep.Adj$gene[i:min(i+each-1,nrow(QTL.CIS.GeneKeep.Adj))],rownames(QTL.PERM.GeneKeep)),]
  QTL.CIS.PERM.GeneKeep.Result[i:min(i+each-1,nrow(QTL.CIS.GeneKeep.Adj))]=
    rowSums(QTL.CIS.PERM.GeneKeep.BLOCK<QTL.CIS.GeneKeep.Adj[i:min(i+each-1,nrow(QTL.CIS.GeneKeep.Adj)),'pvalue'],na.rm = T)
}
QTL.CIS.PERM.GeneKeep.Empirical=(QTL.CIS.PERM.GeneKeep.Result+1)/(ncol(QTL.PERM.GeneKeep)+1)
QTL.CIS.GeneKeep.Adj=data.frame(QTL.CIS.GeneKeep.Adj,empiricalP=QTL.CIS.PERM.GeneKeep.Empirical)
QTL.CIS.GeneKeep.PASS=subset(QTL.CIS.GeneKeep.Adj,empiricalP<=QTL.EmpiricalThreshold.GeneKeep)
QTL.CIS.GeneKeep.PASS$key=paste(QTL.CIS.GeneKeep.PASS$SNP,QTL.CIS.GeneKeep.PASS$gene,sep="-")
length(unique(QTL.CIS.GeneKeep.PASS$SNP))
length(unique(QTL.CIS.GeneKeep.PASS$gene))
length(unlist(strsplit(unique(STA.paper[match(unique(QTL.CIS.GeneKeep.PASS$gene),STA.paper$circRNA.ID),"gene"]),",")))
#unique QTL 185210
#unique circRNA 2790
#unique gene 1623
QTL.CIS.Ratio.Adj=subset(QTL.CIS.Ratio,gene %in% names(which(QTL.PERM.Ratio.Empirical.qvalue$qvalues<0.05)))
QTL.EmpiricalThreshold.Ratio=unique(
  as.numeric(QTL.PERM.Ratio.Empirical.qvalue$pvalues[
    which(QTL.PERM.Ratio.Empirical.qvalue$qvalues==
            max(QTL.PERM.Ratio.Empirical.qvalue$qvalues[
              which(QTL.PERM.Ratio.Empirical.qvalue$qvalues<0.05)]))]))
each=10000
QTL.CIS.PERM.Ratio.Result=vector(length = nrow(QTL.CIS.Ratio.Adj))
for (i in seq(1,nrow(QTL.CIS.Ratio.Adj),by = each)) {
  print(i)
  QTL.CIS.PERM.Ratio.BLOCK=QTL.PERM.Ratio[
    match(QTL.CIS.Ratio.Adj$gene[i:min(i+each-1,nrow(QTL.CIS.Ratio.Adj))],rownames(QTL.PERM.Ratio)),]
  QTL.CIS.PERM.Ratio.Result[i:min(i+each-1,nrow(QTL.CIS.Ratio.Adj))]=
    rowSums(QTL.CIS.PERM.Ratio.BLOCK<QTL.CIS.Ratio.Adj[i:min(i+each-1,nrow(QTL.CIS.Ratio.Adj)),'pvalue'],na.rm = T)
}
QTL.CIS.PERM.Ratio.Empirical=(QTL.CIS.PERM.Ratio.Result+1)/(ncol(QTL.PERM.Ratio)+1)
QTL.CIS.Ratio.Adj=data.frame(QTL.CIS.Ratio.Adj,empiricalP=QTL.CIS.PERM.Ratio.Empirical)
QTL.CIS.Ratio.PASS=subset(QTL.CIS.Ratio.Adj,empiricalP<=QTL.EmpiricalThreshold.Ratio)
QTL.CIS.Ratio.PASS$key=paste(QTL.CIS.Ratio.PASS$SNP,QTL.CIS.Ratio.PASS$gene,sep="-")
length(unique(QTL.CIS.Ratio.PASS$SNP))
length(unique(QTL.CIS.Ratio.PASS$gene))
length(unlist(strsplit(unique(STA.paper[match(unique(QTL.CIS.Ratio.PASS$gene),STA.paper$circRNA.ID),"gene"]),",")))
#unique QTL 179410
#unique circRNA 2795
#unique gene 1652

# QTL.CIS.GeneKeepARatio=intersect(unique(QTL.CIS.Ratio.PASS$gene),unique(QTL.CIS.GeneKeep.PASS$gene))
# QTL.CIS.GeneKeepNoRatio=setdiff(unique(QTL.CIS.GeneKeep.PASS$gene),QTL.CIS.GeneKeepARatio)
# QTL.CIS.RatioNoGeneKeep=setdiff(unique(QTL.CIS.Ratio.PASS$gene),QTL.CIS.GeneKeepARatio)
# QTL.CIS.GeneKeepNoRatio.Intergenic=intersect(QTL.CIS.GeneKeepNoRatio,STA.paper.intergenic$circRNA.ID)
# length(QTL.CIS.GeneKeepARatio)
# length(QTL.CIS.GeneKeepNoRatio)
# length(QTL.CIS.RatioNoGeneKeep)
# length(QTL.CIS.GeneKeepNoRatio.Intergenic)
# length(setdiff(QTL.MINP.Ratio$circID,union(unique(QTL.CIS.Ratio.PASS$gene),unique(QTL.CIS.GeneKeep.PASS$gene))))
#################################
#  #           GeneKeep        #
#R #############################
#a #      Sig     # Non-sig   #
#t ############################
#i #    1514     #    844    #
#o ###########################
#  #1276 (177)  #    6645   #
#############################

#classify the QTL dependend exonic intronic and intergenic
QTL.CIS.GeneKeep.PASS.exon=subset(QTL.CIS.GeneKeep.PASS,gene%in%STA.paper.exon$circRNA.ID)
QTL.CIS.GeneKeep.PASS.intron=subset(QTL.CIS.GeneKeep.PASS,gene%in%STA.paper.intron$circRNA.ID)
QTL.CIS.GeneKeep.PASS.intergenic=subset(QTL.CIS.GeneKeep.PASS,gene%in%STA.paper.intergenic$circRNA.ID)

QTL.CIS.Ratio.PASS.exon=subset(QTL.CIS.Ratio.PASS,gene%in%STA.paper.exon$circRNA.ID)
QTL.CIS.Ratio.PASS.intron=subset(QTL.CIS.Ratio.PASS,gene%in%STA.paper.intron$circRNA.ID)
QTL.CIS.Ratio.PASS.intergenic=subset(QTL.CIS.Ratio.PASS,gene%in%STA.paper.intergenic$circRNA.ID)

length(intersect(QTL.CIS.GeneKeep.PASS.exon$key,QTL.CIS.Ratio.PASS.exon$key))
length(intersect(QTL.CIS.GeneKeep.PASS.intron$key,QTL.CIS.Ratio.PASS.intron$key))
length(intersect(QTL.CIS.GeneKeep.PASS.intergenic$key,QTL.CIS.Ratio.PASS.intergenic$key))

# length(intersect(QTL.CIS.GeneKeep.PASS.exon$key,QTL.CIS.Ratio.PASS.exon$key))/
#   length(union(QTL.CIS.GeneKeep.PASS.exon$key,QTL.CIS.Ratio.PASS.exon$key))
# length(intersect(QTL.CIS.GeneKeep.PASS.intron$key,QTL.CIS.Ratio.PASS.intron$key))/
#   length(union(QTL.CIS.GeneKeep.PASS.intron$key,QTL.CIS.Ratio.PASS.intron$key))
# length(intersect(QTL.CIS.GeneKeep.PASS.intergenic$key,QTL.CIS.Ratio.PASS.intergenic$key))/
#   length(union(QTL.CIS.GeneKeep.PASS.intergenic$key,QTL.CIS.Ratio.PASS.intergenic$key))

length(intersect(QTL.CIS.GeneKeep.PASS.exon$key,QTL.CIS.Ratio.PASS.exon$key))/
  length(QTL.CIS.GeneKeep.PASS.exon$key)
length(intersect(QTL.CIS.GeneKeep.PASS.intron$key,QTL.CIS.Ratio.PASS.intron$key))/
  length(QTL.CIS.GeneKeep.PASS.intron$key)
length(intersect(QTL.CIS.GeneKeep.PASS.intergenic$key,QTL.CIS.Ratio.PASS.intergenic$key))/
  length(QTL.CIS.GeneKeep.PASS.intergenic$key)


QTL.CIS.PASS.merge=merge(QTL.CIS.GeneKeep.PASS,QTL.CIS.Ratio.PASS,by="key")
QTL.CIS.PASS.merge.GeneBody=subset(QTL.CIS.PASS.merge,STA.paper[gene.x,'region']!="intergenic_region")
{
  pdf("./pdf_C/Beta_expression_ratio.pdf")
  par(mar=c(1,1,2,2))
  plot(QTL.CIS.PASS.merge.GeneBody$beta.y~
         QTL.CIS.PASS.merge.GeneBody$beta.x,pch=20,las=1,col=rgb(30,30,30,maxColorValue = 255,alpha = 100),
       xlab="",ylab="")
  abline(v=0,lty=2,lwd=2,col=rgb(30,30,30,maxColorValue = 255))
  abline(h=0,lty=2,lwd=2,col=rgb(30,30,30,maxColorValue = 255))
  dev.off()
}

QTL.CIS.PASS.merge.GeneBody$direct=sign(QTL.CIS.PASS.merge.GeneBody$beta.x)*sign(QTL.CIS.PASS.merge.GeneBody$beta.y)
length(unique(QTL.CIS.PASS.merge$gene.x))
#1991

QTL.CIS.PASS.common=QTL.CIS.PASS.merge[,c('SNP.x','gene.x','beta.x','tstat.x','pvalue.x','FDR.x','empiricalP.x','key')]
colnames(QTL.CIS.PASS.common)=c("SNP","gene","beta","tstat","pvalue","FDR","empiricalP","key")
QTL.CIS.PASS=data.frame(rbind(QTL.CIS.PASS.common,
                              QTL.CIS.GeneKeep.PASS.intergenic))
QTL.CIS.PASS=QTL.CIS.PASS[!duplicated(QTL.CIS.PASS$key),]
nrow(QTL.CIS.PASS)
length(unique(QTL.CIS.PASS$SNP))
length(unique(QTL.CIS.PASS$gene))
length(unlist(strsplit(unique(STA.paper[match(unique(QTL.CIS.PASS$gene),STA.paper$circRNA.ID),"gene"]),",")))
#total circQTL: 196,255
#unique QTL 148,202
#unique circRNA 2,086
#unique gene 1,269
write.table(QTL.CIS.PASS,"./qtl_C/QTL.CIS.PASS.txt",col.names = F,row.names = F,quote=F,sep="\t")
################Cal distance QTL to 5' and 3' site####################
QTL.CIS.PASS.QTL=data.frame(do.call("rbind",strsplit(QTL.CIS.PASS$SNP,split = "_")))
QTL.CIS.PASS.gene=data.frame(do.call("rbind",strsplit(QTL.CIS.PASS$gene,split = "\\|")))
tmp=data.frame(do.call("rbind",strsplit(QTL.CIS.PASS.gene$X1,split = ":")))
QTL.CIS.PASS.gene=data.frame(tmp,x3=QTL.CIS.PASS.gene$X2)
QTL.CIS.PASS.distance=data.frame(matrix(nrow=nrow(QTL.CIS.PASS),ncol=2))
QTL.CIS.PASS.QTL[,3]=as.numeric(QTL.CIS.PASS.QTL[,3])
QTL.CIS.PASS.gene[,2]=as.numeric(QTL.CIS.PASS.gene[,2])
QTL.CIS.PASS.gene[,3]=as.numeric(QTL.CIS.PASS.gene[,3])
QTL.CIS.PASS.distance$X1=QTL.CIS.PASS.QTL[,3]-QTL.CIS.PASS.gene[,2]
QTL.CIS.PASS.distance$X2=QTL.CIS.PASS.QTL[,3]-QTL.CIS.PASS.gene[,3]
##################read ANNOTATION file#########################
CIRI.ANNOTATION=fread("/media/data3/circCMC/data/combineCIRI.annotation.txt",header = F)
setkey(CIRI.ANNOTATION,"V1")
CIRI.ANNOTATION.PASS=CIRI.ANNOTATION[QTL.MINP.GeneKeep$circID,]
rm(CIRI.ANNOTATION)
#################Annotation##################################
QTL.CIS.PASS.Annotation=CIRI.ANNOTATION.PASS[match(QTL.CIS.PASS$gene,CIRI.ANNOTATION.PASS$V1),]
table(QTL.CIS.PASS.Annotation$V4)
QTL.CIS.PASS.distance=data.frame(QTL.CIS.PASS.distance,strand=QTL.CIS.PASS.Annotation$V4)
QTL.CIS.PASS.distance.Detail=data.frame(matrix(nrow=nrow(QTL.CIS.PASS),ncol=2))
for (i in 1:nrow(QTL.CIS.PASS.distance.Detail)){
  m=QTL.CIS.PASS.distance[i,1]*QTL.CIS.PASS.distance[i,2]
  if(i %% 10000 == 0){
    print(i)
  }
  if(QTL.CIS.PASS.distance[i,'strand']=="+"){
    if(m<0){
      QTL.CIS.PASS.distance.Detail[i,1]="Internal"
      QTL.CIS.PASS.distance.Detail[i,2]=abs(QTL.CIS.PASS.distance[i,1])/(abs(QTL.CIS.PASS.distance[i,1])+abs(QTL.CIS.PASS.distance[i,2]))
    }else{
      if(QTL.CIS.PASS.distance[i,1]<=0){
        QTL.CIS.PASS.distance.Detail[i,1]="splice5"
        QTL.CIS.PASS.distance.Detail[i,2]=QTL.CIS.PASS.distance[i,1]
      }else{
        QTL.CIS.PASS.distance.Detail[i,1]="splice3"
        QTL.CIS.PASS.distance.Detail[i,2]=QTL.CIS.PASS.distance[i,2]
      }
    }
  }else{
    if(m<0){
      QTL.CIS.PASS.distance.Detail[i,1]="Internal"
      QTL.CIS.PASS.distance.Detail[i,2]=abs(QTL.CIS.PASS.distance[i,2])/(abs(QTL.CIS.PASS.distance[i,1])+abs(QTL.CIS.PASS.distance[i,2]))
    }else{
      if(QTL.CIS.PASS.distance[i,1]<=0){
        QTL.CIS.PASS.distance.Detail[i,1]="splice3"
        QTL.CIS.PASS.distance.Detail[i,2]=-QTL.CIS.PASS.distance[i,1]
      }else{
        QTL.CIS.PASS.distance.Detail[i,1]="splice5"
        QTL.CIS.PASS.distance.Detail[i,2]=-QTL.CIS.PASS.distance[i,2]
      }
    }
  }
}
QTL.CIS.PASS.distance=data.frame(QTL.CIS.PASS.distance,QTL.CIS.PASS.distance.Detail)
colnames(QTL.CIS.PASS.distance)=c("d1","d2","strand","type","distance")

#private QTL (max p-value per gene)
QTL.private.index=as.data.frame(matrix(1,nrow=length(unique(QTL.CIS.PASS$gene)),ncol=2))
rownames(QTL.private.index)=unique(QTL.CIS.PASS$gene)
for(i in 1:nrow(QTL.CIS.PASS)){
  index=QTL.CIS.PASS[i,'gene']
  if(QTL.private.index[index,1]>=QTL.CIS.PASS[i,'pvalue']){
    QTL.private.index[index,1]=QTL.CIS.PASS[i,'pvalue']
    QTL.private.index[index,2]=i
  }
}
QTL.private.CIS=QTL.CIS.PASS[QTL.private.index$V2,]
write.table(QTL.private.CIS,"./qtl_C/QTL.private.CIS.txt",row.names = F,col.names = T,sep="\t",quote = F)
QTL.private.CIS.distance=QTL.CIS.PASS.distance[QTL.private.index$V2,]
plot(-log(QTL.private.CIS$pvalue[which(QTL.private.CIS.distance$type=="splice5")],
          10)~QTL.private.CIS.distance$distance[which(QTL.private.CIS.distance$type=="splice5")],col = "blue",pch = 19)
plot(-log(QTL.private.CIS$pvalue[which(QTL.private.CIS.distance$type=="Internal")],
          10)~QTL.private.CIS.distance$distance[which(QTL.private.CIS.distance$type=="Internal")],col = "blue",pch = 19)
plot(-log(QTL.private.CIS$pvalue[which(QTL.private.CIS.distance$type=="splice3")],
          10)~QTL.private.CIS.distance$distance[which(QTL.private.CIS.distance$type=="splice3")],col = "blue",pch = 19)

QTL.private.CIS.Detail=QTL.private.CIS
QTL.private.CIS.Detail$type=QTL.private.CIS.distance$type
QTL.private.CIS.Detail$distance=QTL.private.CIS.distance$distance
QTL.private.CIS.Detail.splice5=subset(QTL.private.CIS.Detail,type=="splice5")
QTL.private.CIS.Detail.splice5=QTL.private.CIS.Detail.splice5[
  order(QTL.private.CIS.Detail.splice5$distance,decreasing = T),]

QTL.private.CIS.Detail.splice3=subset(QTL.private.CIS.Detail,type=="splice3")
QTL.private.CIS.Detail.splice3=QTL.private.CIS.Detail.splice3[
  order(QTL.private.CIS.Detail.splice3$distance,decreasing = F),]

QTL.private.CIS.Detail.Internal=subset(QTL.private.CIS.Detail,type=="Internal")
QTL.private.CIS.Detail.Internal=QTL.private.CIS.Detail.Internal[
  order(QTL.private.CIS.Detail.Internal$distance,decreasing = F),]
##############draw the distribution###############

# plot(-log(QTL.private.CIS$pvalue,
#           10)~QTL.private.CIS.distance$distance,col = "blue",pch = 20,cex=0.4,
#      ylab=expression(paste(-log[10],italic(P),"-value")),xlab="Distance to the nearst splice site",main="")
# par(new=TRUE)
# plot(density(QTL.private.CIS.distance$distance,bw=5),ylab="",xlab="",main="",col="red",axes=F,lwd=2)
{
  png(filename="./qtl_C/QTL.private.distance.png",width=480,height = 380,units="px")
  QTL.private.CIS.distance.nonInternal=which(QTL.private.CIS.distance$type!="Internal")
  plot(-log(QTL.private.CIS$pvalue[QTL.private.CIS.distance.nonInternal],
            10)~QTL.private.CIS.distance$distance[QTL.private.CIS.distance.nonInternal],col = "blue",pch = 20,cex=0.4,
       ylab=expression(paste(-log[10]," ",italic(P),"-value")),xlab="Distance to the 5' or 3' splice site (bp)",main="")
  par(new=TRUE)
  plot(density(QTL.private.CIS.distance$distance[QTL.private.CIS.distance.nonInternal],bw=700),ylab="",xlab="",main="",col="red",axes=F,lwd=2)
  dev.off()
}


table(QTL.private.CIS.distance$type)
table(QTL.CIS.PASS.distance$type)
table(QTL.private.CIS.distance$type)/nrow(QTL.private.CIS.distance)

QTL.private.CIS.distance.boxplot=QTL.private.CIS.distance
splitBlock=c(seq(-100000,0,5000),seq(0,1,0.1),seq(1,100001,5000))
QTL.private.CIS.distance.boxplot$block=splitBlock[findInterval(QTL.private.CIS.distance.boxplot$distance,splitBlock)]

{
  pdf("./pdf_C/QTL_Distribution.pdf",width = 10)
  par(mar=c(6,4,1,1))
  tmp=table(QTL.private.CIS.distance.boxplot$block)
  names(tmp)=c("-100k~-95k","-95k~-90k","-90k~-85k","-85k~-80k","-80k~-75k","-75k~-70k",
               "-70k~-65k","-65k~-60k","-60k~-55k","-55k~-50k","-50k~-45k","-45k~-40k","-40k~-35k","-35k~-30k",
               "-30k~-25k","-25k~-20k","-20k~-15k","-15k~-10k","-10k~-5k","-5k~0","0~0.1","0.1~0.2","0.2~0.3",
               "0.3~0.4","0.4~0.5","0.5~0.6","0.6~0.7","0.7~0.8","0.8~0.9","0.9~1","1~5k","5k~10k","10k~15k",
               "15k~20k","20k~25k","25k~30k","30k~35k","35k~40k","40k~45k","45k~50k","50k~55k","55k~60k",
               "60k~65k","65k~70k","70k~75k","75k~80k","80k~85k","85k~90k","90k~95k","95k~100k")
  barplot(tmp,ylab = "Frequence",las=2,horiz=F,ylim = c(0,160),
          names.arg =c("-100k~-95k","","-90k~-85k","","-80k~-75k","",
                       "-70k~-65k","","-60k~-55k","","-50k~-45k","","-40k~-35k","",
                       "-30k~-25k","","-20k~-15k","","-10k~-5k"," ","0~0.1","","0.2~0.3",
                       "","0.4~0.5","","0.6~0.7","","0.8~0.9","","1~5k","","10k~15k",
                       "","20k~25k","","30k~35k","","40k~45k","","50k~55k","",
                       "60k~65k","","70k~75k","","80k~85k","","90k~95k",""),
          col = c(rep("#d53e4f",20),rep("#e6f598",10),rep("#3288bd",20)),
          cex.names = 1)
  legend(-2,150,legend = c("5' (N = 911)",
                           "Internal (N = 423)",
                           "3' (N = 752)"),
         fill=c("#d53e4f","#e6f598","#3288bd"),
         col=c("#d53e4f","#e6f598","#3288bd"),
         y.intersp =1,x.intersp = 0.5,text.width=2,
         bty="n")
  dev.off()
}


################Supplementary Data 4################
QTL.CIS.PASS.paper=QTL.CIS.PASS[,-6]
colnames(QTL.CIS.PASS.paper)[which(colnames(QTL.CIS.PASS.paper)=="gene")]="circRNA ID"
QTL.CIS.PASS.paper$gene=STA.paper[match(QTL.CIS.PASS.paper$`circRNA ID`,STA.paper$circRNA.ID),"gene"]
QTL.CIS.PASS.paper=QTL.CIS.PASS.paper[,c("SNP","circRNA ID","gene","beta","tstat","pvalue","empiricalP")]
colnames(QTL.CIS.PASS.paper)=c("SNP","circRNA ID","circRNA gene","beta","tstat","pvalue","empiricalP")
QTL.CIS.PASS.paper.pos=data.frame(do.call("rbind",strsplit(QTL.CIS.PASS.paper$SNP,"_")))
write.table(QTL.CIS.PASS.paper,"./qtl_C/QTL.CIS.PASS.paper.txt",col.names =T,row.names = F,sep="\t",quote = F)
QTL.CIS.PASS.paper.detail=data.frame(cbind(QTL.CIS.PASS.paper,QTL.CIS.PASS.paper.pos))
QTL.private.CIS.SNP=as.data.frame(do.call("rbind",strsplit(QTL.private.CIS$SNP,"_")))

#########################QTL colocalation##################
QTL.coloc=cbind(QTL.private.CIS,QTL.private.CIS.SNP)
QTL.coloc$parentalGene=STA.paper[QTL.coloc$gene,'gene']
QTL.coloc.sub=QTL.coloc[nchar(QTL.coloc$parentalGene)==15,]
QTL.coloc.sub$rsidGene=with(QTL.coloc.sub,V2%&%"_"%&%parentalGene)

###############Compare with eQTL in CMC####################
CMC.eQTL=read.table("/media/data3/CMC/data/qtl/CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedSVA.txt.gz",header=T)
CMC.eQTL.cis=subset(CMC.eQTL,eQTL_type=="cis")
CMC.eQTL.cis$rsidGene=with(CMC.eQTL.cis,SNP%&%"_"%&%Gene)

length(QTL.coloc.sub$rsidGene)
sum(QTL.coloc.sub$rsidGene%in%intersect(CMC.eQTL.cis$rsidGene,QTL.coloc.sub$rsidGene))
851/1750

CMC.eQTL.cis.overlap=subset(CMC.eQTL.cis,rsidGene%in%QTL.coloc.sub$rsidGene)
CMC.eQTL.cis.merge=merge(CMC.eQTL.cis.overlap,QTL.coloc.sub,by="rsidGene")
CMC.eQTL.cis.merge$MAF=SNP.MAF5[CMC.eQTL.cis.merge$SNP.y,'V2']/930
CMC.eQTL.cis.merge$direction=sign(CMC.eQTL.cis.merge$Beta)*sign(CMC.eQTL.cis.merge$beta)
CMC.eQTL.cis.merge.neg=subset(CMC.eQTL.cis.merge,direction=="-1")
CMC.eQTL.cis.merge.pos=subset(CMC.eQTL.cis.merge,direction=="1")


{
  pdf("./pdf_C/eQTL_circQTL.pdf")
  par(mfrow=c(1,1))
  ecQTL.Beta=data.frame(gene=-CMC.eQTL.cis.merge$Beta,circ=CMC.eQTL.cis.merge$beta)
  ecQTL.Beta$ID=CMC.eQTL.cis.merge$gene
  
  
  # ggplot(ecQTL.Beta,aes(x=gene,y = circ))+
  #   geom_point(color="red")+
  #   geom_text(aes(label=ID),nudge_y=0,check_overlap = T,nudge_x = -0.5)+theme_pubr()
  
  
  plot(-Beta~beta,CMC.eQTL.cis.merge,pch=20,col=rgb(30,30,30,maxColorValue = 255,alpha = 200),
       xlab="Beta from circQTL",ylab="Beta from eQTL",las=1)
  abline(v=0,lty=2,lwd=2,col=rgb(30,30,30,maxColorValue = 255))
  abline(h=0,lty=2,lwd=2,col=rgb(30,30,30,maxColorValue = 255))
  dev.off()
}
summary(lm(Beta~beta,CMC.eQTL.cis.merge))
table(CMC.eQTL.cis.merge$direction)
table(CMC.eQTL.cis.merge$direction)/851
#detail in the disconcordant QTLs
CMC.SVA.gene=read.delim("/media/data3/CMC/data/expression/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_gene-adjustedSVA-dataNormalization-includeAncestry-adjustedL.tsv.gz")
rownames(CMC.SVA.gene)=CMC.SVA.gene[,1]
CMC.SVA.gene=CMC.SVA.gene[,-1]
CMC.eQTL.cis.merge.neg.P=vector(length = nrow(CMC.eQTL.cis.merge.neg))
CMC.eQTL.cis.merge.neg.Beta=vector(length = nrow(CMC.eQTL.cis.merge.neg))

for (i in 1:nrow(CMC.eQTL.cis.merge.neg)) {
  tmpGene=CMC.eQTL.cis.merge.neg[i,'Gene']
  tmpCirc=CMC.eQTL.cis.merge.neg[i,'gene']
  tmp=summary(lm(as.numeric(CMC.SVA.gene[tmpGene,CaucasianID])~as.numeric(CIRI.JUNCTION.PASS.voom.log.Res[tmpCirc,CaucasianID])))
  CMC.eQTL.cis.merge.neg.P[i]=tmp$coefficients[2,4]
  CMC.eQTL.cis.merge.neg.Beta[i]=tmp$coefficients[2,1]
}
i=138
plot(as.numeric(CMC.SVA.gene[tmpGene,CaucasianID])~as.numeric(CIRI.JUNCTION.PASS.voom.log.Res[tmpCirc,CaucasianID]))
for (i in 1:nrow(CMC.eQTL.cis.merge)) {
  target.SNP=CMC.eQTL.cis.merge[i,'SNP.y']
  chr=CMC.eQTL.cis.merge[i,'V1']
  out="./qtl_C/genotype/"%&%target.SNP%&%".genotype"
  cmd="grep \'"%&%target.SNP%&%"\' /media/data3/CMC/data/genotype/Caucasian465/SNPmatrix/"%&%chr%&%".genotype.MAF5.mat " %&%">"%&%out 
  if(!file.exists(out)){
    try(system(cmd))
  }
}
tmpGene="ENSG00000165181"
tmpCirc="chr9:114549073|114550227"
tmpgenoytpe=read.table("./qtl_C/genotype/chr9_rs1853224_114457917_A_C.genotype",header=F)
tmpgenoytpe=as.numeric(tmpgenoytpe[1,-1])
plot(as.numeric(CMC.SVA.gene[tmpGene,CaucasianID])~tmpgenoytpe)
plot(as.numeric(CIRI.JUNCTION.PASS.voom.log.Res[tmpCirc,CaucasianID])~tmpgenoytpe)
CMC.eQTL.cis.merge.pos=CMC.eQTL.cis.merge.pos[order(CMC.eQTL.cis.merge.pos$pvalue.y,decreasing = F),]

CMC.eQTL.cis.merge.pos.choice=subset(CMC.eQTL.cis.merge.pos,pvalue.x<10^(-10) &pvalue.y<10^(-10))
CMC.eQTL.cis.merge.pos.choice$distance=QTL.private.CIS.Detail$distance[match(CMC.eQTL.cis.merge.pos.choice$gene,
                                                                             QTL.private.CIS.Detail$gene)]
par(mfrow=c(4,4),mar=c(1,2,2,1))
for (i in 1:nrow(CMC.eQTL.cis.merge.pos.choice)) {
  tmpGene=CMC.eQTL.cis.merge.pos.choice[i,'Gene']
  tmpCirc=CMC.eQTL.cis.merge.pos.choice[i,'gene']
  tmpgenoytpe=read.table("./qtl_C/genotype/"%&%CMC.eQTL.cis.merge.pos.choice[i,'SNP.y']%&%".genotype",header=F)
  tmpgenoytpe=as.numeric(tmpgenoytpe[1,-1])
  expGene=as.numeric(Tophat2.DGEList.voom.log[tmpGene,CaucasianID])
  expCirc=as.numeric(CIRI.JUNCTION.PASS.voom.log[tmpCirc,CaucasianID])
  tmpP=summary(lm(expGene~c(tmpgenoytpe)))$coefficients[2,4]
  a1=CMC.eQTL.cis.merge.pos.choice[i,'V4']
  a2=CMC.eQTL.cis.merge.pos.choice[i,'V5']
  tmpLabel=c(a1%&%a1,a1%&%a2,a2%&%a2)
  tmpgenoytpe=tmpLabel[fenRA(tmpgenoytpe)+1]
  #expGene=as.numeric(CMC.SVA.gene[tmpGene,CaucasianID])
  
   plot(expGene~factor(tmpgenoytpe),main=i%&%" "%&%tmpP%&%"\n"%&%
          CMC.eQTL.cis.merge.pos.choice[i,'pvalue.y'],pch=20,col="blue",ylim= range(c(expCirc,expGene)),
       xaxt="n",ylab="",xlab="")
   par(new=T)
   plot(expCirc~factor(tmpgenoytpe),pch=20,col="red",axes=F)
  # abline(lm(expGene~c(tmpgenoytpe-0.1)),col="blue")
  # abline(lm(expCirc~c(tmpgenoytpe+0.1)),col="red")
  
  plot(expGene~expCirc,pch=20,type="n",main=i)
  par(new=T)
  idx1=which(tmpgenoytpe==tmpLabel[1])
  plot(expGene[idx1]~expCirc[idx1],pch=20,col="red",axes=F)
  abline(lm(expGene[idx1]~expCirc[idx1]),col="red")
  par(new=T)
  idx2=which(tmpgenoytpe==tmpLabel[2])
  plot(expGene[idx2]~expCirc[idx2],pch=20,col="blue",axes=F)
  abline(lm(expGene[idx2]~expCirc[idx2]),col="blue")
  par(new=T)
  idx3=which(tmpgenoytpe==tmpLabel[3])
  plot(expGene[idx3]~expCirc[idx3],pch=20,col="green",axes=F)
  abline(lm(expGene[idx3]~expCirc[idx3]),col="green")
  abline(lm(expGene~expCirc),col="black")
  
  # expGene=as.numeric(CMC.SVA.gene[tmpGene,CaucasianID])
  # points(expGene~c(tmpgenoytpe-0.1),main=i,pch=20,col="green",xlim=c(-0.5,2.5),ylim= range(c(expCirc,expGene)),
  #      xaxt="n",ylab="",xlab="")
 
}
{
  i=5
  tmpGene=CMC.eQTL.cis.merge.pos.choice[i,'Gene']
  tmpCirc=CMC.eQTL.cis.merge.pos.choice[i,'gene']
  tmpgenoytpe=read.table("./qtl_C/genotype/"%&%CMC.eQTL.cis.merge.pos.choice[i,'SNP.y']%&%".genotype",header=F)
  tmpgenoytpe=as.numeric(tmpgenoytpe[1,-1])
  expGene=as.numeric(Tophat2.DGEList.voom.log[tmpGene,CaucasianID])
  expCirc=as.numeric(CIRI.JUNCTION.PASS.voom.log[tmpCirc,CaucasianID])
  tmpP=summary(lm(expGene~c(tmpgenoytpe)))$coefficients[2,4]
  a1=CMC.eQTL.cis.merge.pos.choice[i,'V4']
  a2=CMC.eQTL.cis.merge.pos.choice[i,'V5']
  tmpLabel=c(a1%&%a1,a1%&%a2,a2%&%a2)
  tmpgenoytpe=tmpLabel[fenRA(tmpgenoytpe)+1]
  tmp.df=data.frame(type=rep(c("Gene","CircRNA"),
                             each=length(expGene)),exp=c(expGene,expCirc),genotype=c(tmpgenoytpe,tmpgenoytpe))
  p=ggplot(tmp.df,aes(x=genotype,y=exp))+
    geom_boxplot(aes(fill=genotype))+ 
    theme_bw(base_size = 12)+
    theme(legend.position="top")+
    facet_wrap(~ type, scales="free")+labs(y = "Adjusted expression level (log2 CPM)")
  pdf("./pdf_C/boxplot_circQTL_eQTL.pdf")
  print(p)
  dev.off()
  tmp2.df=data.frame(gene=expGene,circ=expCirc,genotype=tmpgenoytpe)
  p2=ggplot(tmp2.df,aes(x=expGene,y=expCirc))+
    stat_smooth(method="lm",se=F,aes(colour=genotype))+
    geom_point(aes(col=genotype))+ 
    theme_bw(base_size = 12)+
    theme(legend.position="none")
  pdf("./pdf_C/expression_circQTL_eQTL.pdf")
  print(p2)
  dev.off()
}

#################annotation circQTL and non-circQTL###################
#read SNP and indel pos
cisDist=10^5
SNP.POS=read.delim("/media/data3/CMC/data/genotype/Caucasian465/chrALL.genotype.MAF5.info80.pos",header=F)

chrMax = max( STA.pos[, 4], SNP.POS[, 3], na.rm = TRUE) + cisDist
chrName=unique(c(STA.pos$chr,SNP.POS$V2))
chrName.Num=as.numeric(as.factor(chrName))
names(chrName.Num)=chrName
SNP.POS.Transform=(chrName.Num[SNP.POS$V2]-1)*chrMax+SNP.POS$V3
names(SNP.POS.Transform)=SNP.POS$V1
SNP.POS.Transform.order=sort.list(SNP.POS.Transform)
SNP.POS.Transform=SNP.POS.Transform[SNP.POS.Transform.order]

CIRI.POS.Transform.left=(chrName.Num[STA.pos$chr]-1)*chrMax+STA.pos$left
CIRI.POS.Transform.right=(chrName.Num[STA.pos$chr]-1)*chrMax+STA.pos$right
CIRI.POS.Transform.order=sort.list(CIRI.POS.Transform.left+CIRI.POS.Transform.right)
CIRI.POS.Transform.left=CIRI.POS.Transform.left[CIRI.POS.Transform.order]
CIRI.POS.Transform.right=CIRI.POS.Transform.right[CIRI.POS.Transform.order]

sn.l = findInterval(CIRI.POS.Transform.left - cisDist-1, SNP.POS.Transform);
sn.r = findInterval(CIRI.POS.Transform.right+ cisDist-1, SNP.POS.Transform);
SNPInter=unique(unlist(lapply(which(sn.r>sn.l),FUN=function(x){c(sn.l[x]+1):sn.r[x]})))
SNP.Inter=names(SNP.POS.Transform[SNPInter])
length(SNP.Inter)
length(intersect(SNP.Inter,QTL.CIS.PASS$SNP))
#in 100 kb of circRNA 1,611,445 SNPs

###################circQTL and non-circQTL comparison######################
#LD-based pruning of these max-QTL SNPs (exclude INDEL)
QTL.private.QTL=QTL.CIS.PASS.QTL[QTL.private.index$V2,]
QTL.private.QTL.SNP=subset(QTL.private.QTL,X4%in%baseAlp & X5%in%baseAlp)
QTL.private.QTL.SNP$key=with(QTL.private.QTL.SNP,X1%&%"_"%&%X3)
QTL.private.QTL.SNP.nonDuplicated=QTL.private.QTL.SNP[!duplicated(QTL.private.QTL.SNP$key),]

#LD-based pruning of these non-QTL SNPs (exclude INDEL)
SNP.Inter.POS=data.frame(do.call("rbind",strsplit(setdiff(SNP.Inter,unique(c(QTL.CIS.GeneKeep$SNP,QTL.CIS.Ratio$SNP))),"_")))#attention
SNP.Inter.POS.SNP=subset(SNP.Inter.POS,X4%in%baseAlp & X5%in%baseAlp)
SNP.Inter.POS.SNP$key=with(SNP.Inter.POS.SNP,X1%&%"_"%&%X3)
SNP.Inter.POS.SNP.nonDuplicated=SNP.Inter.POS.SNP[!duplicated(SNP.Inter.POS.SNP$key),]
dim(SNP.Inter.POS.SNP.nonDuplicated)
#uncorrected P > 0.05 951,208

#read MAF file
SNP.MAF5=fread("/media/data3/CMC/data/genotype/Caucasian465/chrALL.genotype.SNP.MAF5.info80.MAF",header = F)
setkey(SNP.MAF5,"V1")
SNP.MAF5$V3=SNP.MAF5$V2/930
SNP.MAF5$V3=apply(SNP.MAF5,1,function(x){ifelse(as.numeric(x[3])>0.5,1-as.numeric(x[3]),as.numeric(x[3]))})
QTL.private.QTL.SNP.nonDuplicated$allkey=with(QTL.private.QTL.SNP.nonDuplicated,X1%&%"_"%&%X2%&%"_"%&%X3%&%"_"%&%X4%&%"_"%&%X5)
SNP.Inter.POS.SNP.nonDuplicated$allkey=with(SNP.Inter.POS.SNP.nonDuplicated,X1%&%"_"%&%X2%&%"_"%&%X3%&%"_"%&%X4%&%"_"%&%X5)
#SNP.Inter.POS.SNP.nonDuplicated=SNP.Inter.POS.SNP.nonDuplicated[-which(SNP.Inter.POS.SNP.nonDuplicated$X2%in%QTL.private.QTL.SNP.nonDuplicated$X2),]
QTL.private.QTL.SNP.nonDuplicated$MAF=SNP.MAF5[QTL.private.QTL.SNP.nonDuplicated$allkey,'V2']/(465*2)
QTL.private.QTL.SNP.nonDuplicated$MAF=
  unlist(lapply(unlist(QTL.private.QTL.SNP.nonDuplicated$MAF), function(x){min(x,1-x)}))
SNP.Inter.POS.SNP.nonDuplicated$MAF=SNP.MAF5[SNP.Inter.POS.SNP.nonDuplicated$allkey,"V2"]/(465*2)
SNP.Inter.POS.SNP.nonDuplicated$MAF=
  unlist(lapply(unlist(SNP.Inter.POS.SNP.nonDuplicated$MAF), function(x){min(x,1-x)}))
hist(as.numeric(QTL.private.QTL.SNP.nonDuplicated$MAF),breaks = 100)
hist(as.numeric(SNP.Inter.POS.SNP.nonDuplicated$MAF),breaks = 100)

###bin the MAF
MAF.BIN=seq(0.05,0.5,0.025)
QTL.private.QTL.SNP.nonDuplicated$bin=MAF.BIN[findInterval(QTL.private.QTL.SNP.nonDuplicated$MAF,
                                                           MAF.BIN)]
SNP.Inter.POS.SNP.nonDuplicated$bin=MAF.BIN[findInterval(SNP.Inter.POS.SNP.nonDuplicated$MAF,
                                                         MAF.BIN)]
hist(QTL.private.QTL.SNP.nonDuplicated$bin,breaks = 20)

##PLINK extract SNPs files
QTL.private.QTL.SNP.nonDuplicated.PLINK=with(QTL.private.QTL.SNP.nonDuplicated,
                                             data.frame(chr=gsub("chr","",X1),
                                                        start=X3,
                                                        end=X3,label=X1%&%"_"%&%X3))
write.table(QTL.private.QTL.SNP.nonDuplicated.PLINK,"./qtl_C/PLINK/circQTL/QTL.private.QTL.SNP.nonDuplicated.PLINK",
            row.names = F,col.names = F,sep=" ",quote=F)
SNP.Inter.POS.SNP.nonDuplicated.PLINK=with(SNP.Inter.POS.SNP.nonDuplicated,
                                           data.frame(chr=gsub("chr","",X1),
                                                      start=X3,
                                                      end=X3,label=X1%&%"_"%&%X3))
write.table(SNP.Inter.POS.SNP.nonDuplicated.PLINK,"./qtl_C/PLINK/noncircQTL/SNP.Inter.POS.SNP.nonDuplicated.PLINK",
            row.names = F,col.names = F,sep=" ",quote=F)
##########################PLINK code##########################################
#  for i in {1..22};do plink --file /media/data3/1000Genome/eurPLINK_MAF5/$i \
#  --extract  range QTL.private.QTL.SNP.nonDuplicated.PLINK \
#  --make-bed --out $i;done
# for i in {1..22};do plink --bfile $i  --indep-pairwise 50 5 0.5 -out $i;done
# for i in {1..22};do plink --bfile $i --extract $i.prune.in --make-bed -out ${i}pd ;done
# for i in {1..22};do cat ${i}pd.bim ;done>ALLpd.bim
##############################################################################
# for i in {1..22};do plink --file /media/data3/1000Genome/eurPLINK_MAF5/$i \
#  --extract  range SNP.Inter.POS.SNP.nonDuplicated.PLINK \
#  --make-bed --out $i;done
# for i in {1..22};do plink --bfile $i  --indep-pairwise 50 5 0.5 -out $i;done
# for i in {1..22};do plink --bfile $i --extract $i.prune.in --make-bed -out ${i}pd ;done
# for i in {1..22};do cat ${i}pd.bim ;done>ALLpd.bim
########################################################
##read bim file
QTL.private.bim=read.delim("./qtl_C/PLINK/circQTL/ALLpd.bim",header = F)
QTL.private.bim$key=with(QTL.private.bim,"chr"%&%V1%&%"_"%&%V4)
SNP.Inter.bim=read.delim("./qtl_C/PLINK/noncircQTL/ALLpd.bim",header=F)
SNP.Inter.bim$key=with(SNP.Inter.bim,"chr"%&%V1%&%"_"%&%V4)

#extract pruning SNPs from nonDuplicated
QTL.private.QTL.SNP.nonDuplicated.pruning=
  QTL.private.QTL.SNP.nonDuplicated[match(QTL.private.bim$key,QTL.private.QTL.SNP.nonDuplicated$key),]
QTL.private.QTL.SNP.nonDuplicated.pruning$bin[which(QTL.private.QTL.SNP.nonDuplicated.pruning$bin==0.5)]=0.475
SNP.Inter.POS.SNP.nonDuplicated.pruning=
  SNP.Inter.POS.SNP.nonDuplicated[match(SNP.Inter.bim$key,SNP.Inter.POS.SNP.nonDuplicated$key),]
SNP.Inter.POS.SNP.nonDuplicated.pruning$bin[which(SNP.Inter.POS.SNP.nonDuplicated.pruning$bin==0.5)]=0.475
QTL.pruning.Cumulative.percentage=cumsum(table(QTL.private.QTL.SNP.nonDuplicated.pruning$bin))/nrow(QTL.private.QTL.SNP.nonDuplicated.pruning)
SNP.pruning.Cumulative.percentage=cumsum(table(SNP.Inter.POS.SNP.nonDuplicated.pruning$bin))/nrow(SNP.Inter.POS.SNP.nonDuplicated.pruning)

plot(QTL.pruning.Cumulative.percentage~seq(0.0625,0.4875,0.025),type="b",
     pch=20,col="red",xlim=c(0,0.51),ylim=c(0,1),
     xlab="Minor allel frequence",ylab="Cumulative percentage")
par(new=T)
plot(SNP.pruning.Cumulative.percentage~seq(0.0625,0.4875,0.025),type="b",ylab="",xlab="",
     pch=20,col="blue",axes=F,xlim=c(0,0.51),ylim=c(0,1))

#MAF match SNP pruning set
QTL.pruning.table.percentage=table(QTL.private.QTL.SNP.nonDuplicated.pruning$bin)/nrow(QTL.private.QTL.SNP.nonDuplicated.pruning)
SNP.pruning.table.percentage=table(SNP.Inter.POS.SNP.nonDuplicated.pruning$bin)/nrow(SNP.Inter.POS.SNP.nonDuplicated.pruning)
SNP.pruning.table=table(SNP.Inter.POS.SNP.nonDuplicated.pruning$bin)
SNP.pruning.ratio.QTL=SNP.pruning.table/QTL.pruning.table.percentage
SNP.pruning.matrix.QTL=SNP.pruning.ratio.QTL %*% t(QTL.pruning.table.percentage)
SNP.pruning.Diff.QTL=rowSums(t(t(SNP.pruning.matrix.QTL)-as.numeric(SNP.pruning.table))<=0)
SNP.pruning.choiceNum=floor(SNP.pruning.matrix.QTL[which(SNP.pruning.Diff.QTL==(length(MAF.BIN)-1)),])
SNP.pruning.choiceIndex=vector()
for(i in 1:length(SNP.pruning.choiceNum)){
  SNP.pruning.choiceIndex=c(SNP.pruning.choiceIndex,
                            sample(which(SNP.Inter.POS.SNP.nonDuplicated.pruning$bin==MAF.BIN[i]),size=SNP.pruning.choiceNum[i]))
}
SNP.pruning.choiceDf=SNP.Inter.POS.SNP.nonDuplicated.pruning[SNP.pruning.choiceIndex,]
SNP.pruning.choiceDf.Cumulative.percentage=cumsum(table(SNP.pruning.choiceDf$bin))/nrow(SNP.pruning.choiceDf)
{
  png("./qtl_C/MAF.cumulative.png",width=530,height=380)
  plot(QTL.pruning.Cumulative.percentage*100~seq(0.0625,0.4875,0.025),type="b",
       pch=16,col="red",xlim=c(0,0.51),ylim=c(0,100),lwd=3,bty="l",las=1,
       xlab="Minor allel frequence",ylab="Cumulative percentage")
  par(new=T)
  plot(SNP.pruning.Cumulative.percentage*100~seq(0.0625,0.4875,0.025),type="b",ylab="",xlab="",
       pch=16,col="blue",axes=F,lwd=3,xlim=c(0,0.51),ylim=c(0,100))
  par(new=T)
  plot(SNP.pruning.choiceDf.Cumulative.percentage*100+rnorm(18,sd=0.15)~seq(0.0625,0.4875,0.025),type="b",ylab="",xlab="",
       pch=16,col=rgb(0,1,0,0.5),axes=F,lwd=3,xlim=c(0,0.51),ylim=c(0,100))
  nrow(SNP.Inter.POS.SNP.nonDuplicated.pruning)
  nrow(SNP.pruning.choiceDf)
  nrow(QTL.private.QTL.SNP.nonDuplicated.pruning)
  
  #126,502
  #47,472
  #1,284
  legend(x=0,y=105,legend = c("All non-circQTL SNPs\n(N = 126,502)",
                              "MAF-matched non-circQTL SNPs\n(N = 47,472)",
                              "circQTL SNPs\n(N = 1,284)"),
         fill=c("blue","red",rgb(0,1,0,0.5)),
         col=c("blue","red",rgb(0,1,0,0.5)),
         y.intersp =1.35,x.intersp = 0.5,text.width=0.9,
         bty="n")
  dev.off()
}

##########################Enrichment analysis########################
#########################Format transformation#######################
########################transform to VCF format#####################
SNP.pruning.VCF=with(SNP.pruning.choiceDf,
                     data.frame(chr=X1,pos=X3,rsid=X2,ref=X4,alt=X5))
QTL.pruning.VCF=with(QTL.private.QTL.SNP.nonDuplicated.pruning,
                     data.frame(chr=X1,pos=X3,rsid=X2,ref=X4,alt=X5))
length(intersect(SNP.pruning.VCF$rsid,QTL.pruning.VCF$rsid))
write.table(SNP.pruning.VCF,"./qtl_C/SNP.pruning.VCF",sep="\t",row.names = F,col.names = F,quote = F)
write.table(QTL.pruning.VCF,"./qtl_C/QTL.pruning.VCF",sep="\t",row.names = F,col.names = F,quote = F)
#####################transform to BED format###################
SNP.pruning.BED=with(SNP.pruning.choiceDf,
                     data.frame(chr=X1,start=as.character(as.numeric(X3)-1),end=X3))
QTL.pruning.BED=with(QTL.private.QTL.SNP.nonDuplicated.pruning,
                     data.frame(chr=X1,start=as.character(as.numeric(X3)-1),end=X3))
write.table(SNP.pruning.BED,"./qtl_C/SNP.pruning.bed",sep="\t",row.names = F,col.names = F,quote = F)
write.table(QTL.pruning.BED,"./qtl_C/QTL.pruning.bed",sep="\t",row.names = F,col.names = F,quote = F)
###############################VEP annotation#########################
# vep -i /media/data3/circCMC/data/qtl_C/QTL.pruning.VCF \
# -o /media/data3/circCMC/data/qtl_C/VEP/QTL.pruning.most_severe.VEP \
# -gff /media/data3/annotation/GENCODE/gencode.v19.annotation.gff3.gz \
# -fasta /media/data3/genome/hg19/hg19.fa --most_severe --force_overwrite  2>1.log &
# 
# sort  -k 1,1d -k 2,2n /media/data3/circCMC/data/qtl_C/SNP.pruning.VCF >/media/data3/circCMC/data/qtl_C/SNP.pruning.sort.VCF
# vep -i /media/data3/circCMC/data/qtl_C/SNP.pruning.sort.VCF \
#  -o /media/data3/circCMC/data/qtl_C/VEP/SNP.pruning.sort.most_severe.VEP \
#  -gff /media/data3/annotation/GENCODE/gencode.v19.annotation.gff3.gz \
#  -fasta /media/data3/genome/hg19/hg19.fa --most_severe --force_overwrite 2>2.log &
#######################################################################
#annotation type
VEP.type=c("3_prime_UTR_variant","5_prime_UTR_variant","downstream_gene_variant",
           "intergenic_variant","intron_variant","missense_variant",
           "non_coding_transcript_exon_variant","splice_acceptor_variant",
           "splice_donor_variant","splice_region_variant","start_lost",
           "stop_lost","synonymous_variant",
           "upstream_gene_variant")
VEP.type.change=c("3'-UTR","5'-UTR","downstream gene","intergenic","intron","missense",
                  "non-coding transcript exon","canonical splice site","canonical splice site",
                  "splice region variant","loss-of-function","loss-of-function","synonymous",
                  "upstream gene")
names(VEP.type.change)=VEP.type
#read VEP file
QTL.pruning.VEP=read.delim("./qtl_C/VEP/QTL.pruning.most_severe.VEP",header = F,comment.char = "#")
QTL.pruning.VEP.noDuplicated=QTL.pruning.VEP[!duplicated(QTL.pruning.VEP$V1),]
QTL.pruning.VEP.noDuplicated$type=VEP.type.change[QTL.pruning.VEP.noDuplicated$V7]
QTL.pruning.VEP.noDuplicated.table=table(QTL.pruning.VEP.noDuplicated$type)
QTL.pruning.VEP.noDuplicated.percentage=QTL.pruning.VEP.noDuplicated.table/nrow(QTL.pruning.VEP.noDuplicated)*100

SNP.pruning.VEP=read.delim("./qtl_C/VEP/SNP.pruning.sort.most_severe.VEP",header = F,comment.char = "#")
SNP.pruning.VEP.noDuplicated=SNP.pruning.VEP[!duplicated(SNP.pruning.VEP$V1),]
SNP.pruning.VEP.noDuplicated$type=VEP.type.change[SNP.pruning.VEP.noDuplicated$V7]
SNP.pruning.VEP.noDuplicated.table=table(SNP.pruning.VEP.noDuplicated$type)
SNP.pruning.VEP.noDuplicated.percentage=SNP.pruning.VEP.noDuplicated.table/nrow(SNP.pruning.VEP.noDuplicated)*100

VEP.feature=intersect(names(QTL.pruning.VEP.noDuplicated.table),
                      names(SNP.pruning.VEP.noDuplicated.table))
VEP.fisher.test=data.frame(matrix(nrow=length(VEP.feature),ncol=4))
rownames(VEP.fisher.test)=VEP.feature
for (i in 1:length(VEP.feature)) {
  m=QTL.pruning.VEP.noDuplicated.table[VEP.feature[i]]
  n=SNP.pruning.VEP.noDuplicated.table[VEP.feature[i]]
  f=fisher.test(matrix(c(m,nrow(QTL.pruning.VEP.noDuplicated)-m,
                         n,nrow(SNP.pruning.VEP.noDuplicated)-n),nrow=2))
  VEP.fisher.test[i,]=as.numeric(c(f$p.value,f$estimate,m,n))
}
colnames(VEP.fisher.test)=c("p","OR","m","n")
VEP.fisher.test$fdr=p.adjust(VEP.fisher.test$p,"BH")

{
  pdf("./pdf_C/QTL_annotation_distribution.pdf",width = 7,height=6)
  
  tmp=rownames(subset(VEP.fisher.test,fdr<2 & m>-1))
  Fig3DB=data.frame(anno=rep(tmp,times=2),
                    Frequency=c(as.numeric(QTL.pruning.VEP.noDuplicated.percentage[tmp]/100),
                                as.numeric(SNP.pruning.VEP.noDuplicated.percentage[tmp]/100)),
                    type=rep(c("zcircQTL","non-circQTL"),each=length(tmp)))
  Fig3DB$per=round(100*Fig3DB$Frequency,1)
  p=ggplot(Fig3DB,aes(anno,per,fill=type))+
    geom_bar(stat="identity",position = "dodge")+theme_pubr()+
    #theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),legend.title = element_text(size = 0))+labs(x = "")+
    theme(legend.title = element_text(size = 0))+labs(y = "Percentage",x="Annotation")+
    scale_fill_manual(values=c('#377eb8','#e41a1c'))+
    geom_text(aes(x=anno,y=per+5,label=per%&%"%"),
              
              position=position_dodge(width=0.9),
              show.legend = F)+
    coord_flip()
  print(p)
  dev.off()
}

{
  pdf("./pdf_C/bubbleDistribution.pdf",width=8)
  symbols(VEP.fisher.test$OR,-log(VEP.fisher.test$fdr,10),
          #circle=sqrt(VEP.fisher.test$m/2),
          circle = rep(1,nrow(VEP.fisher.test)),
          las=1,bty="l",
          inches = 0.17,fg="white",bg="lightblue",
          main="",xlab="odds ratio",ylab=expression(paste(-log[10]," ","FDR")))
  #plot(VEP.fisher.test$OR,-log(VEP.fisher.test$fdr,10),pch=19,col="lightblue")
  text(VEP.fisher.test$OR,-log(VEP.fisher.test$fdr,10),
       rownames(VEP.fisher.test),cex =1)
  abline(h=2,lty=2,col="red")
  abline(v=1,lty=2,col="red")
  dev.off()
}
QTL.pruning.VEP.noDuplicated.CIS=as.data.frame(matrix(n=ncol(QTL.pruning.VEP.noDuplicated),
                                                      nrow = nrow(QTL.private.CIS.Detail)))
for (i in 1:nrow(QTL.pruning.VEP.noDuplicated.CIS)) {
  tmp=QTL.pruning.VEP.noDuplicated[which(QTL.pruning.VEP.noDuplicated$V1==
                                          QTL.private.CIS.SNP$V2[i]),]
  if(nrow(tmp)==1){
    QTL.pruning.VEP.noDuplicated.CIS[i,]=as.vector(tmp)
  }else{
    print(nrow(tmp))
  }
}
QTL.pruning.VEP.noDuplicated.CIS=as.data.frame(cbind(QTL.pruning.VEP.noDuplicated.CIS,
                                                     QTL.private.CIS.Detail))
QTL.pruning.VEP.noDuplicated.splice.region=subset(QTL.pruning.VEP.noDuplicated.CIS,V15=="splice region variant")
QTL.pruning.VEP.noDuplicated.canonical.splice.site=subset(QTL.pruning.VEP.noDuplicated.CIS,V15=="canonical splice site")

################genotype chr2_rs563601_45928630_G_A###########
genotype2=read.delim("/media/data3/CMC/data/genotype/Caucasian465/chr2_rs563601_45928630_G_A.mat",header=F,sep=" ")
genotype2=as.numeric(genotype2[,-1])
genotype2=apply(abs(genotype2-matrix(0:2,nrow = length(genotype2),ncol=3,byrow = T)),1,which.min)-1
boxplot(CIRI.JUNCTION.PASS.voom.log.Res["chr2:45925080|45928630",CaucasianID]~genotype2)
tmpRatio=as.numeric(CIRI.RATIO.mat.zero.Res["chr2:45925080|45928630",CaucasianID])
boxplot(tmpRatio~genotype2)


########################Repeat masker annotation########################
# bedtools intersect -a ./QTL.pruning.bed \
#  -b /media/data3/annotation/RepeatMasker/rmsk.bed.gz -wao >repeat/QTL.pruning.rmsk &
# ##############################################################################
# sort -k 1,1 -k 2,2n ./SNP.pruning.bed >./SNP.pruning.sort.bed
# bedtools intersect -a ./SNP.pruning.sort.bed \
#  -b /media/data3/annotation/RepeatMasker/rmsk.bed.gz -wao >repeat/SNP.pruning.rmsk &
###############################################################################
RMSK.QTL=read.delim("./qtl_C/repeat/QTL.pruning.rmsk",header=F)
RMSK.SNP=read.delim("./qtl_C/repeat/SNP.pruning.rmsk",header = F)
RMSK.QTL$key=with(RMSK.QTL,V1%&%":"%&%V3)
RMSK.SNP$key=with(RMSK.SNP,V1%&%":"%&%V3)
RMSK.QTL=RMSK.QTL[match(QTL.pruning.VEP$V2,RMSK.QTL$key),]
RMSK.SNP=RMSK.SNP[match(SNP.pruning.VEP.noDuplicated$V2,RMSK.SNP$key),]
RMSK.QTL.table=table(RMSK.QTL$V7)
RMSK.SNP.table=table(RMSK.SNP$V7)
RMSK.feature=intersect(names(RMSK.QTL.table),names(RMSK.SNP.table))

RMSK.fisher.test=data.frame(matrix(nrow = length(RMSK.feature),ncol=4))
for (i in 1:length(RMSK.feature)) {
  m=RMSK.QTL.table[RMSK.feature[i]]
  n=RMSK.SNP.table[RMSK.feature[i]]
  tmp=fisher.test(matrix(c(m,nrow(QTL.private.QTL.SNP.nonDuplicated.pruning)-m,
                           n,nrow(SNP.pruning.choiceDf)-n),nrow=2))
  RMSK.fisher.test[i,]=c(tmp$p.value,tmp$estimate,m,n)
}
colnames(RMSK.fisher.test)=c("p","OR","m","n")
rownames(RMSK.fisher.test)=RMSK.feature
RMSK.fisher.test$fdr=p.adjust(RMSK.fisher.test$p,"bonf")
#ALU enrichment 
#
##############################Intron repeat annotation###############################
#only consider SNPs located in intronic region
###########################################################################
RMSK.QTL.Intron=RMSK.QTL[which(QTL.pruning.VEP$V7=="intron_variant"),]
RMSK.SNP.Intron=RMSK.SNP[which(SNP.pruning.VEP.noDuplicated$V7=="intron_variant"),]
RMSK.QTL.Intron.table=table(RMSK.QTL.Intron$V7)
RMSK.SNP.Intron.table=table(RMSK.SNP.Intron$V7)
RMSK.Intron.feature=intersect(names(RMSK.QTL.Intron.table),names(RMSK.SNP.Intron.table))

RMSK.Intron.fisher.test=data.frame(matrix(nrow = length(RMSK.Intron.feature),ncol=4))
for (i in 1:length(RMSK.Intron.feature)) {
  m=RMSK.QTL.Intron.table[RMSK.Intron.feature[i]]
  n=RMSK.SNP.Intron.table[RMSK.Intron.feature[i]]
  tmp=fisher.test(matrix(c(m,nrow(RMSK.QTL.Intron)-m,
                           n,nrow(RMSK.SNP.Intron)-n),nrow=2))
  RMSK.Intron.fisher.test[i,]=c(tmp$p.value,tmp$estimate,m,n)
}
colnames(RMSK.Intron.fisher.test)=c("p","OR","m","n")
rownames(RMSK.Intron.fisher.test)=RMSK.Intron.feature
RMSK.Intron.fisher.test$fdr=p.adjust(RMSK.Intron.fisher.test$p,"BH")
#############################Different ALU region######################
RMSK.QTL.ALU=subset(RMSK.QTL,V7=="Alu")
RMSK.SNP.ALU=subset(RMSK.SNP,V7=="Alu")
RMSK.QTL.ALU.table=table(RMSK.QTL.ALU$V10)
RMSK.SNP.ALU.table=table(RMSK.SNP.ALU$V10)
RMSK.ALU.feature=intersect(names(RMSK.QTL.ALU.table),names(RMSK.SNP.ALU.table))
RMSK.ALU.fisher.test=data.frame(matrix(nrow = length(RMSK.ALU.feature),ncol=4))
for (i in 1:length(RMSK.ALU.feature)) {
  m=RMSK.QTL.ALU.table[RMSK.ALU.feature[i]]
  n=RMSK.SNP.ALU.table[RMSK.ALU.feature[i]]
  tmp=fisher.test(matrix(c(m,nrow(RMSK.QTL.ALU)-m,
                           n,nrow(RMSK.SNP.ALU)-n),nrow=2))
  RMSK.ALU.fisher.test[i,]=c(tmp$p.value,tmp$estimate,m,n)
}
colnames(RMSK.ALU.fisher.test)=c("p","OR","m","n")
rownames(RMSK.ALU.fisher.test)=RMSK.ALU.feature
RMSK.ALU.fisher.test$fdr=p.adjust(RMSK.ALU.fisher.test$p,"fdr")

#############################Intron base pairing enrichment#########
circFlankIntron.pos=read.delim("/media/data3/circCMC/data/FlankingIntron/circFlankIntron.pos",header=F)
table(STA.paper[circFlankIntron.pos$V1,'region'])
circFlankIntron.Blast=read.delim("/media/data3/circCMC/data/FlankingIntron/BlastResult.txt",header = F)
circFlankIntron.Blast=subset(circFlankIntron.Blast,V1%in%STA.paper.exon$circRNA.ID)
circFlankIntron.Blast.bed=data.frame(chr=STA.pos[circFlankIntron.Blast[,1],'chr'],
                                     start=circFlankIntron.Blast[,2]-1,
                                     end=circFlankIntron.Blast[,3])

circFlankIntron.Blast.bed=rbind(circFlankIntron.Blast.bed,data.frame(chr=STA.pos[circFlankIntron.Blast[,1],'chr'],
                                start=circFlankIntron.Blast[,5]-1,
                                end=circFlankIntron.Blast[,4]))
circFlankIntron.Blast.bed[,2]=as.character(circFlankIntron.Blast.bed[,2])
circFlankIntron.Blast.bed[,3]=as.character(circFlankIntron.Blast.bed[,3])
write.table(circFlankIntron.Blast.bed,"./FlankingIntron/circFlankIntron.Blast.bed",row.names = F,col.names = F,sep="\t",quote = F)
######################
#cd /media/data3/circCMC/data/FlankIntron/
#sort -k 1,1 -k 2,2g  -k  3,3g  circFlankIntron.Blast.bed|uniq >circFlankIntron.Blast.sort.bed
#cd /media/data3/circCMC/data/qtl_C/
# bedtools intersect -a ./QTL.pruning.bed 
#  -b ../FlankingIntron/circFlankIntron.Blast.sort.bed -wao >./FlankIntron/QTL.pruning.FlankIntron &
# bedtools intersect -a ./SNP.pruning.sort.bed \
#  -b ../FlankingIntron/circFlankIntron.Blast.sort.bed -wao >./FlankIntron/SNP.pruning.FlankIntron &
#########################
QTL.pruning.FlankIntron=read.delim("./qtl_C/FlankIntron/QTL.pruning.FlankIntron",header = F)
QTL.pruning.FlankIntron$SNPID=with(QTL.pruning.FlankIntron,V1%&%"_"%&%V3)
QTL.pruning.FlankIntron=subset(QTL.pruning.FlankIntron,V7!=0)
SNP.pruning.FlankIntron=read.delim("./qtl_C/FlankIntron/SNP.pruning.FlankIntron",header = F)
SNP.pruning.FlankIntron$SNPID=with(SNP.pruning.FlankIntron,V1%&%"_"%&%V3)
SNP.pruning.FlankIntron=subset(SNP.pruning.FlankIntron,V7!=0)
m=length(unique(QTL.pruning.FlankIntron$SNPID))
n=length(unique(SNP.pruning.FlankIntron$SNPID))
fisher.test(matrix(c(m,n,nrow(QTL.pruning.BED)-m,nrow(SNP.pruning.BED)-n),nrow=2))
QTL.FI.CIS=data.frame(QTL.private.CIS[match(QTL.pruning.FlankIntron$SNPID,QTL.private.QTL.SNP$key),],QTL.pruning.FlankIntron[,c(4,5,6)])

#only intron variant
QTL.pruning.FlankIntron$SNPID=with(QTL.pruning.FlankIntron,V1%&%":"%&%V3)
SNP.pruning.FlankIntron$SNPID=with(SNP.pruning.FlankIntron,V1%&%":"%&%V3)
QTL.pruning.VEP.intron=QTL.pruning.VEP[QTL.pruning.VEP$V7=="intron_variant",'V2']
SNP.pruning.VEP.noDuplicated.intron=SNP.pruning.VEP.noDuplicated[which(SNP.pruning.VEP.noDuplicated$V7=="intron_variant"),'V2']

m=length(intersect(unique(QTL.pruning.FlankIntron$SNPID),QTL.pruning.VEP.intron))
n=length(intersect(unique(SNP.pruning.FlankIntron$SNPID),SNP.pruning.VEP.noDuplicated.intron))

fisher.test(matrix(c(m,n,length(QTL.pruning.VEP.intron)-m,length(SNP.pruning.VEP.noDuplicated.intron)-n),nrow=2))


###############################GWAS enrichment#########################
#read GWAS file
GWAS=read.delim("/media/data3/GWAS/NIHcatalogy/GWAS_my_2_25_2018_hg19.txt")
GWAS=subset(GWAS,P.VALUE< 5* 10^(-8))
GWAS.SCZ=read.delim("/media/data3/GWAS/PGC/SCZ2/scz2.annel.108.position",header=F)
GWAS.key=unique(with(GWAS,seqnames%&%"_"%&%start))
GWAS.SCZ.key=unique(with(GWAS.SCZ,"chr"%&%V1%&%"_"%&%V2))

GWAS.list=list()
n=1
for (i in 1:nrow(GWAS)) {
  uri=unlist(strsplit(GWAS[i,"MAPPED_TRAIT_URI"],","))
  trait=unlist(strsplit(GWAS[i,"MAPPED_TRAIT"],","))
  num=length(uri)
  for (j in 1:num) {
    eachTrack=GWAS[i,]
    eachTrack["MAPPED_TRAIT_URI"]=uri[j]
    eachTrack["MAPPED_TRAIT"]=trait[j]
    GWAS.list[[n]]=eachTrack
    n=n+1
  }
}
GWAS.list=do.call("rbind",GWAS.list)
GWAS.list=GWAS.list[-which(is.na(GWAS.list$MAPPED_TRAIT_URI)),]
#EFO ID tags
GWAS.EFO.index=regexpr("EFO_",GWAS.list$MAPPED_TRAIT_URI)
GWAS.EFO=vector(length = nrow(GWAS.list))
for (i in 1:nrow(GWAS.list)) {
  if(GWAS.EFO.index[i]>0){
    GWAS.EFO[i]=substring(GWAS.list$MAPPED_TRAIT_URI[i],GWAS.EFO.index[i],GWAS.EFO.index[i]+10)
  }
  
}
GWAS.list$EFO=GWAS.EFO
writeLines(unique(GWAS.EFO),"GWAS.EFO.unique")
#######solve use ontoCAT in my local computer##############
#read the file
GWAS.EFO.isDisease=read.delim("efoDisese.df")
GWAS.EFO.disease=subset(GWAS.EFO.isDisease,isDisease==1)
GWAS.disease=subset(GWAS.list,EFO%in%GWAS.EFO.disease$GWAS.EFO)
GWAS.disease$efo.SNP=GWAS.disease$SNPS%&%"_"%&%GWAS.disease$EFO
GWAS.disease=GWAS.disease[!duplicated(GWAS.disease$efo.SNP),]
GWAS.disease$Trait=GWAS.EFO.isDisease[match(GWAS.disease$EFO,GWAS.EFO.isDisease$GWAS.EFO),3]
GWAS.rsID=unique(GWAS$SNPS)
writeLines(GWAS.rsID,"./qtl_C/GWAS.rsID")

###############################PLINK LD tag this GWAS SNPs############################
#  for i in {1..22};do plink \
# --file /media/data3/1000Genome/eurPLINK_MAF5/$i \
# --show-tags ../GWAS.rsID  --list-all -out $i  --tag-r2 0.6;done
# cat *.tags >ALL.tags
# cat *.tags.list >ALL.tags.list
# #########################################
# cd ../GWAS2
# for i in {1..22};do plink \
# --file /media/data3/1000Genome/eurPLINK_MAF5/$i \
# --show-tags ../GWAS.rsID  --list-all -out $i  --tag-r2 0.8;done
# cat *.tags >ALL.tags
# cat *.tags.list >ALL.tags.list
##########################################################################
#manually inspect max-circQTL SNPs with high LD GWAS SNPs
#######################
GWAS.tag.HLD.list=read.table("./qtl_C/GWAS2/ALL.tags.list",header=T)
GWAS.tag.HLD.list.disease=subset(GWAS.tag.HLD.list,SNP%in%unique(GWAS.disease$SNPS))
GWAS.tag.HLD.list.disease$TAGS=GWAS.tag.HLD.list.disease$TAGS%&%"|"%&%GWAS.tag.HLD.list.disease$SNP
GWAS.HLD.disease.SNPs=unique(unlist(strsplit(GWAS.tag.HLD.list.disease$TAGS,"\\|")))
GWAS.HLD.disease.SNPs.intersect=unique(intersect(QTL.CIS.PASS.paper.detail$X2,GWAS.HLD.disease.SNPs))
GWAS.HLD.inCircRNA.df=subset(QTL.CIS.PASS.paper.detail,X2%in%GWAS.HLD.disease.SNPs.intersect)

idxSNP=vector(length = nrow(GWAS.HLD.inCircRNA.df))
gwasDf=vector(length = nrow(GWAS.HLD.inCircRNA.df))
disDf=vector(length = nrow(GWAS.HLD.inCircRNA.df))
for (i in 1:nrow(GWAS.HLD.inCircRNA.df)) {
  tmpSNP=GWAS.HLD.inCircRNA.df[i,'X2']
  tmpidx=GWAS.tag.HLD.list.disease[grep(tmpSNP,GWAS.tag.HLD.list.disease$TAGS),'SNP']
  idxSNP[i]=paste(tmpidx,collapse = ";")
  idxSNPDP=vector(length=length(tmpidx))
  tmpDis=vector(length = length(tmpidx))
  for (j in 1:length(tmpidx)) {
    tmpGWAS=subset(GWAS.disease,SNPS==tmpidx[j])
    tmpDetail=vector(length = nrow(tmpGWAS))
    for (h in 1:nrow(tmpGWAS)) {
      tmpDetail[h]= paste(tmpGWAS[h,c("MAPPED_GENE",'Trait',"PUBMEDID")],collapse = '|')
    }
    idxSNPDP[j]=paste(tmpDetail,collapse = "+")
    tmpDis[j]=unique(tmpGWAS$Trait)
  }
  gwasDf[i]=paste(idxSNPDP,collapse = ";")
  disDf[i]=paste(unique(tmpDis),collapse = ";")
}

GWAS.HLD.inCircRNA.df$idxSNP=idxSNP
GWAS.HLD.inCircRNA.df$gwas=gwasDf
GWAS.HLD.inCircRNA.df$disease=disDf
GWAS.HLD.inCircRNA.df.paper=GWAS.HLD.inCircRNA.df[,c(1,2,13,14,15)]
write.table(GWAS.HLD.inCircRNA.df.paper,"./qtl_C/GWAS.HLD.inCircRNA.df.paper.txt",row.names = F,col.names = T,sep="\t",quote=F)
GWAS.HLD.inCircRNA.df.collapse=tapply(GWAS.HLD.inCircRNA.df$disease,
                                      as.factor(GWAS.HLD.inCircRNA.df$circRNA.ID),function(x){unique(unlist(strsplit(x,";")))})
sort(table(unlist(GWAS.HLD.inCircRNA.df.collapse)))
GWAS.HLD.inCircRNA.df.paper.SCZ=GWAS.HLD.inCircRNA.df.paper[grep("schizophrenia",GWAS.HLD.inCircRNA.df.paper$disease),]
GWAS.HLD.inCircRNA.df.paper.IBD=GWAS.HLD.inCircRNA.df.paper[grep("inflammatory bowel disease",GWAS.HLD.inCircRNA.df.paper$disease),]

#SCZ HLD
SCZ.HLD=subset(GWAS.tag.HLD.list.disease,SNP%in%unique(subset(GWAS.disease,Trait=="schizophrenia")$SNPS))
SCZ.HLD.maxcircQTL=intersect(unique(unlist(strsplit(SCZ.HLD$TAGS,"\\|"))),QTL.private.QTL.SNP.nonDuplicated.pruning$X2)
SCZ.HLD.df=QTL.private.CIS[which(QTL.private.CIS.SNP$V2%in%SCZ.HLD.maxcircQTL),]
for (i in 1:nrow(SCZ.HLD.df)) {
  tmpSNP=SCZ.HLD.maxcircQTL[i]
  tmpidx=SCZ.HLD[grep(tmpSNP,SCZ.HLD$TAGS),'SNP']
  tmpGWAS=subset(GWAS.disease,SNPS%in%tmpidx)
  print(tmpGWAS$JOURNAL)
  
}

#read the tag file
GWAS.tag.LD=readLines("./qtl_C/GWAS/ALL.tags")
GWAS.tag.LD.list=read.table("./qtl_C/GWAS/ALL.tags.list",header=F)
GWAS.tag.LD.list.disease=subset(GWAS.tag.LD.list,V1%in%unique(GWAS.disease$SNPS))
GWAS.disease.SNPs=unique(unlist(strsplit(GWAS.tag.LD.list.disease$V8,"\\|")))
########################GWAS analysis############################
m3=length(intersect(GWAS.disease.SNPs,QTL.private.QTL.SNP.nonDuplicated.pruning$X2))
n3=length(intersect(GWAS.disease.SNPs,SNP.pruning.choiceDf$X2))

GWAS.disease.all=fisher.test(matrix(c(m3,nrow(QTL.private.QTL.SNP.nonDuplicated.pruning)-m3,
                                      n3,nrow(SNP.pruning.choiceDf)-n3),nrow=2),alternative = "greater")
GWAS.disease.all$p.value
GWAS.disease.all$estimate

#P=3.919629e-17
#OR= 3.526221

GWAS.disease$Trait[which(GWAS.disease$Trait=="autism spectrum disorder")]="autism"
GWAS.disease$Trait[which(GWAS.disease$Trait=="Crohn's disease")]="inflammatory bowel disease"
GWAS.disease$Trait[which(GWAS.disease$Trait=="ulcerative colitis")]="inflammatory bowel disease"
GWAS.disease.table=sort(table(GWAS.disease$Trait),decreasing = T)
#gwas fisher test function
gwasFisherTest=function(disease){
  if(disease %in% names(GWAS.disease.table)){
    GWAS.tag.LD.list.disease=subset(GWAS.tag.LD.list,V1%in%unique(subset(GWAS.disease,Trait==disease)$SNPS))
    GWAS.disease.SNPs=unique(unlist(strsplit(GWAS.tag.LD.list.disease$V8,"\\|")))
    m3=length(intersect(GWAS.disease.SNPs,QTL.private.QTL.SNP.nonDuplicated.pruning$X2))
    n3=length(intersect(GWAS.disease.SNPs,SNP.pruning.choiceDf$X2))
    tmp=fisher.test(matrix(c(m3,nrow(QTL.private.QTL.SNP.nonDuplicated.pruning)-m3,
                             n3,nrow(SNP.pruning.choiceDf)-n3),nrow=2),alternative = "greater")
    
    return(c(tmp$estimate,tmp$p.value,m3,n3))
  }else{
    print("no this disease")
  }
}
#get SNPs > 100
GWAS.threshold=90
GWAS.disease.table.threshold=GWAS.disease.table[which(GWAS.disease.table>GWAS.threshold)]
  # GWAS.BRAIN=c("Alzheimer's disease","bipolar disorder","autism",
  #              "Parkinson's disease")
 GWAS.disease.table.threshold=GWAS.disease.table[unique(c(names(GWAS.disease.table.threshold),GWAS.BRAIN))]
GWAS.disease.ohterFT=data.frame(matrix(nrow=length(GWAS.disease.table.threshold),ncol = 6))
GWAS.disease.ohterFT$X1=names(GWAS.disease.table.threshold)
GWAS.disease.ohterFT$X2=GWAS.disease.table.threshold

for (i in 1:length(GWAS.disease.table.threshold)) {
  tmp=gwasFisherTest(GWAS.disease.ohterFT[i,1])
  GWAS.disease.ohterFT[i,3:6]=tmp
  
}
colnames(GWAS.disease.ohterFT)=c("disease","riskNum","OR","p","m","n")
GWAS.disease.ohterFT$fdr=p.adjust(GWAS.disease.ohterFT$p,"fdr")


#############################GWAS_disease_distribution##############################
{
  pdf("./pdf_C/GWAS_disease_distribution.pdf",height=6)
  layout(matrix(c(1,1,1,1,1,2), nrow = 1, ncol = 6, byrow = TRUE))
  par(mar=c(4,4,1,1))
  symbols(GWAS.disease.ohterFT$OR,-log(GWAS.disease.ohterFT$fdr,10),
          circle=(GWAS.disease.ohterFT$m+1)^(1/2),las=1,bty="l",
          inches = 0.2,fg=rgb(240,48,39,maxColorValue = 255),
          bg=rgb(240,48,39,maxColorValue = 255,alpha = 150),
          main="",xlab="odds ratio",ylab=expression(paste(-log[10]," ",FDR)))
  indx=which(GWAS.disease.ohterFT$OR>1 & GWAS.disease.ohterFT$fdr<0.05)
  text(GWAS.disease.ohterFT$OR[indx],-log(GWAS.disease.ohterFT$fdr[indx],10),
       GWAS.disease.ohterFT$disease[indx],cex =1)
  abline(h=-log(0.05,10),lty=1,col="#0005FF",lwd=2)
  abline(h=-log(0.01,10),lty=5,col="#0005FF",lwd=2)
  abline(v=1,lty=5,col="#0005FF",lwd=2)
  par(mar=c(2,0,0,1))
  #plot(2,axes=F,type="n")
  #points(rep(1,4),seq(2.15,2.35,length.out = 4),cex=seq(1,4,1),col=rgb(240,48,39,maxColorValue = 255,alpha = 150),pch=19)
  symbols(rep(1,4),seq(1.5,2.2,length.out = 4),circle=seq(1,4,1),ylim = c(0,3),inches=0.2,axes=F,
          fg=rgb(240,48,39,maxColorValue = 255),bg=rgb(240,48,39,maxColorValue = 255,alpha = 150))
  text(rep(1.25,4),seq(1.5,2.2,length.out = 4),labels = c(0,3,8,15),cex=1.5)
  text(1,2.41,label="Count",cex=1.5)
  dev.off() 
}
ggplot(GWAS.disease.ohterFT,aes(x=OR,y = -log(fdr,10)))+
  geom_point(size=sqrt(GWAS.disease.ohterFT$m/2),color="red")+
  geom_text(aes(label=disease),nudge_y=-0.25,check_overlap = T)+theme_pubr()

#########SCZ GWAS disease##########
disease="schizophrenia"
GWAS.tag.LD.list.disease=subset(GWAS.tag.LD.list,V1%in%unique(subset(GWAS.disease,Trait==disease)$SNPS))
GWAS.disease.SNPs=unique(unlist(strsplit(GWAS.tag.LD.list.disease$V8,"\\|")))
SCZ.GWAS.inQTL=intersect(GWAS.disease.SNPs,QTL.private.QTL.SNP.nonDuplicated.pruning$X2)
SCZ.risk.loci=list()
n=0
for(i in 1:nrow(GWAS.tag.LD.list.disease)){
  print(i)
  tmp=unlist(strsplit(GWAS.tag.LD.list.disease[i,8],"\\|"))
  for(j in 1:length(tmp)){
    n=n+1
    SCZ.risk.loci[[n]]=c(tmp[j],as.character(GWAS.tag.LD.list.disease[i,1:3]))
  }
}
SCZ.risk.loci=as.data.frame(do.call("rbind",SCZ.risk.loci))
colnames(SCZ.risk.loci)=c("LD_SNPs","GWAS_SNPs","GWAS_chr","GWAS_position")
SCZ.risk.loci.inQTL=merge(QTL.private.QTL.SNP.nonDuplicated.pruning,SCZ.risk.loci,by.x="X2",by.y="LD_SNPs")
SCZ.risk.loci.inQTL.detail=merge(SCZ.risk.loci.inQTL,GWAS,by.x="GWAS_SNPs",by.y="SNPS")
colnames(SCZ.risk.loci.inQTL.detail)[2]="LD_SNPs"
SCZ.risk.loci.inQTL.detail=merge(SCZ.risk.loci.inQTL.detail,QTL.pruning.VEP.noDuplicated[,c(1,2,7,15)],by.x="LD_SNPs",by.y="V1")
SCZ.risk.loci.inQTL.detail$distance=as.numeric(SCZ.risk.loci.inQTL.detail$X3)-as.numeric(SCZ.risk.loci.inQTL.detail$GWAS_position)
SCZ.risk.loci.inQTL.detail=merge(SCZ.risk.loci.inQTL.detail,QTL.CIS.PASS.paper,by.x="allkey",by.y = "SNP")
CMC.gene.DE.result=read.delim("/media/data3/CMC/data/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_gene-adjustedNoSVA-differentialExpression-includeAncestry-DxSCZ-DE.tsv")
rownames(CMC.gene.DE.result)=CMC.gene.DE.result$genes
SCZ.risk.loci.inQTL.detail$symbal=CMC.gene.DE.result[match(SCZ.risk.loci.inQTL.detail$`circRNA gene`,CMC.gene.DE.result$genes),"MAPPED_genes"]
SCZ.risk.loci.inQTL.detail.pgc=subset(SCZ.risk.loci.inQTL.detail,STUDY=="Biological insights from 108 schizophrenia-associated genetic loci.")

##############################################
#########find PGC causal variants#############
#############################################
PGC.scz2.credible.snps=read.delim("/media/data3/GWAS/PGC/SCZ2/pgc.scz2.credible.snps.txt")
#Pay attention:
#the rs115329265 chromosome is chr6 not chr5 because the wrong happend in pgc.scz2.credible.snps.txt
PGC.scz2.credible.snps[which(PGC.scz2.credible.snps$indexSNP=="rs115329265"),2]=6

PGC.scz2.credible.snps=subset(PGC.scz2.credible.snps,crediblePval <= (5*10^(-8)))
PGC.inQTL=unique(intersect(PGC.scz2.credible.snps$credibleSNP,QTL.CIS.PASS.paper.pos$X2))
PGC.scz2.credible.snps.inQTL=subset(PGC.scz2.credible.snps,credibleSNP %in% PGC.inQTL)
length(unique(PGC.scz2.credible.snps.inQTL$indexSNP))
##information of LD in the 1000 Genomes March 2012 EUR data set was not available for the index SNPs
##rs58120505 for chr7_2025096_I
tmp1=subset(PGC.scz2.credible.snps,indexSNP=="chr7_2025096_I")
# ##rs7951870 for chr11_46350213_D
# tmp2=subset(PGC.scz2.credible.snps,indexSNP=="chr11_46350213_D")
PGC.scz2.credible.snps.inQTL[which(PGC.scz2.credible.snps.inQTL$indexSNP=="chr7_2025096_I"),
                             'indexSNP']="rs58120505"
# PGC.scz2.credible.snps.inQTL[which(PGC.scz2.credible.snps.inQTL$indexSNP=="chr11_46350213_D"),
#                              'indexSNP']="rs7951870"
unique(PGC.scz2.credible.snps.inQTL$indexSNP)
write.table(unique(PGC.scz2.credible.snps.inQTL$indexSNP),"./qtl_C/PGC.scz2.credible.snps.inQTL.rsID",col.names = F,row.names = F,quote = F,sep="")

#############################coloc test#####################
# library(coloc)
# QTL.CIS.PASS.paper.detail.inPGC=subset(QTL.CIS.PASS.paper.detail,X2%in%PGC.scz2.credible.snps.inQTL$credibleSNP)
# choice.region=unique(QTL.CIS.PASS.paper.detail.inPGC$circRNA.ID)
# coloc.pgc.df=data.frame(matrix(nrow=length(choice.region),ncol=7))
# coloc.eQTL.df=data.frame(matrix(nrow=length(choice.region),ncol=7))
# 
# for (i in 1:length(choice.region)) {
#   tmpRegion=choice.region[i]
#   tmpQTL=subset(QTL.CIS.PASS.paper.detail.inPGC,circRNA.ID==tmpRegion)
#   tmpQTL$rsidGene=with(tmpQTL,paste(X2,circRNA.gene,sep = "_"))
#   #tmpQTL=tmpQTL[which.max(tmpQTL$pvalue),]
#   tmpPGC=PGC.scz2.credible.snps.inQTL[match(tmpQTL$X2,PGC.scz2.credible.snps.inQTL$credibleSNP),]
#   tmpres=coloc.abf(dataset1=list(pvalues=tmpQTL$pvalue,N=465,type="quant"),
#                    dataset2=list(pvalues=tmpPGC$crediblePval,N=150064,type="cc",s=0.246),
#                    MAF=unlist(SNP.MAF5[tmpQTL$SNP,'V3']),p12=10^(-5))
#   coloc.pgc.df[i,]=c(tmpRegion,tmpres[[1]])
#   tmpeQTL_circ=merge(CMC.eQTL.cis,tmpQTL,by="rsidGene")
#   if(nrow(tmpeQTL_circ)>0){
#     print("eQTL"%&%i)
#   }
# }
# coloc.pgc.df$geneSymbol=STA.paper[coloc.pgc.df$X1,'geneSymbol']
# coloc.pgc.df$gene=STA.paper[coloc.pgc.df$X1,'gene']

#############PLINK LD cal#######################
# cd /media/data3/circCMC/data/qtl_C/PLINK/LocusZoom
# for i in {1..22};do plink --file /media/data3/1000Genome/eurPLINK_MAF5/${i}\
#  --r2 --ld-snp-list ../../PGC.scz2.credible.snps.inQTL.rsID \
#  --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 -out ${i};done
#  for i in {1..22};do tail -n +2 ${i}.ld ;done> ALL.ld
###############################################
PGC.scz2.credible.snps.LD=read.table("./qtl_C/PLINK/LocusZoom/ALL.ld",header=F)
colnames(PGC.scz2.credible.snps.LD)=c( "CHR_A","BP_A","SNP_A","CHR_B","BP_B","SNP_B","R2")
PGC.scz2.credible.snps.LD.PASS=subset(PGC.scz2.credible.snps.LD,R2>0.25)
PGC.scz2.credible.snps.inQTL.PASS=subset(PGC.scz2.credible.snps.inQTL,credibleSNP%in%PGC.scz2.credible.snps.LD.PASS$SNP_B)
unique(PGC.scz2.credible.snps.inQTL.PASS$indexSNP)

QTL.CIS.PASS.info=data.frame(QTL.CIS.PASS.paper,QTL.CIS.PASS.paper.pos)
QTL.CIS.PASS.info[which(QTL.CIS.PASS.info$X2=="."),"X2"]=QTL.CIS.PASS.info[which(QTL.CIS.PASS.info$X2=="."),"X1"]%&%"_"%&%QTL.CIS.PASS.info[which(QTL.CIS.PASS.info$X2=="."),"X3"]

PGC.scz2.credible.snps.inQTL.PASS.detail=merge(PGC.scz2.credible.snps.inQTL.PASS,QTL.CIS.PASS.info,by.x="credibleSNP",by.y="X2")
PGC.scz2.credible.snps.inQTL.PASS.detail.SNP.pos=as.data.frame(do.call("rbind",strsplit(PGC.scz2.credible.snps.inQTL.PASS.detail$SNP,"_")))
PGC.scz2.credible.snps.inQTL.PASS.detail$DisLeft=as.numeric(PGC.scz2.credible.snps.inQTL.PASS.detail.SNP.pos$V3)-as.numeric(STA.pos[PGC.scz2.credible.snps.inQTL.PASS.detail$circRNA.ID,'left'])
PGC.scz2.credible.snps.inQTL.PASS.detail$DisRight=as.numeric(PGC.scz2.credible.snps.inQTL.PASS.detail.SNP.pos$V3)-as.numeric(STA.pos[PGC.scz2.credible.snps.inQTL.PASS.detail$circRNA.ID,'right'])
PGC.scz2.credible.snps.inQTL.PASS.detail$minDist=apply(PGC.scz2.credible.snps.inQTL.PASS.detail,1,function(x){min(abs(as.numeric(x[c("DisLeft","DisRight")])))})
index=tapply(PGC.scz2.credible.snps.inQTL.PASS$crediblePval,as.factor(PGC.scz2.credible.snps.inQTL.PASS$indexSNP),which.min)
PGC.scz2.credible.snps.inQTL.PASS.min=tapply(1:nrow(PGC.scz2.credible.snps.inQTL.PASS),
                                             PGC.scz2.credible.snps.inQTL.PASS$indexSNP,function(x){PGC.scz2.credible.snps.inQTL.PASS[x,][index[unique(PGC.scz2.credible.snps.inQTL.PASS[x,'indexSNP'])],]})
PGC.scz2.credible.snps.inQTL.PASS.min=do.call("rbind",PGC.scz2.credible.snps.inQTL.PASS.min)
#rs58120505 indexSNP_POS is wrong (PGC's data bug, not caused by my analysis)
PGC.scz2.credible.snps.inQTL.PASS.min[which(PGC.scz2.credible.snps.inQTL.PASS.min$indexSNP=="rs58120505"),"indexSNP_POS"]=2029867
##make QTL.CIS.PASS with QTL.CIS.PASS.paper.pos
PGC.scz2.min.QTL=merge(PGC.scz2.credible.snps.inQTL.PASS.min,QTL.CIS.PASS.info,by.x="credibleSNP",by.y="X2")
colnames(PGC.scz2.min.QTL)[which(colnames(PGC.scz2.min.QTL)=="X3")]="credibleSNP_POS"
PGC.scz2.min.QTL$distance=as.numeric(PGC.scz2.min.QTL$indexSNP_POS)-as.numeric(PGC.scz2.min.QTL$credibleSNP_POS)
PGC.scz2.min.QTL$credibleSNP_POS=as.numeric(PGC.scz2.min.QTL$credibleSNP_POS)
PGC.scz2.min.QTL.detail=data.frame(PGC.scz2.min.QTL,Pgene=QTL.CIS.PASS.paper[match(PGC.scz2.min.QTL$circRNA.ID,QTL.CIS.PASS.paper$`circRNA ID`),3])
PGC.scz2.min.QTL.detail$symbol=CMC.gene.DE.result[match(PGC.scz2.min.QTL.detail$Pgene,CMC.gene.DE.result$genes),2]

PGC.scz2.min.QTL.detail.circRNA=do.call("rbind",strsplit(PGC.scz2.min.QTL.detail$circRNA.ID,":"))
PGC.scz2.min.QTL.detail.circRNA=data.frame(circChr=PGC.scz2.min.QTL.detail.circRNA[,1],do.call("rbind",strsplit(PGC.scz2.min.QTL.detail.circRNA[,2],"\\|")))
colnames(PGC.scz2.min.QTL.detail.circRNA)=c("circChr","circStart","circEnd")
PGC.scz2.min.QTL.detail=data.frame(PGC.scz2.min.QTL.detail,PGC.scz2.min.QTL.detail.circRNA)
PGC.scz2.min.QTL.detail$credibleSNP_POS=as.numeric(PGC.scz2.min.QTL.detail$credibleSNP_POS)
PGC.scz2.min.QTL.detail$circStart=as.numeric(PGC.scz2.min.QTL.detail$circStart)
PGC.scz2.min.QTL.detail$circEnd=as.numeric(PGC.scz2.min.QTL.detail$circEnd)
PGC.scz2.min.QTL.detail$CSNP2Start=PGC.scz2.min.QTL.detail$credibleSNP_POS-PGC.scz2.min.QTL.detail$circStart
PGC.scz2.min.QTL.detail$CSNP2End=PGC.scz2.min.QTL.detail$credibleSNP_POS-PGC.scz2.min.QTL.detail$circEnd
PGC.scz2.min.QTL.detail$CSNP2Near=apply(PGC.scz2.min.QTL.detail[,c("CSNP2Start","CSNP2End")],1,function(x){min(abs(x))})
PGC.scz2.min.QTL$ensembl=STA.paper[match(PGC.scz2.min.QTL$circRNA.ID,STA.paper$circRNA.ID),"gene"]
PGC.scz2.min.QTL$geneSymbol=STA.paper[match(PGC.scz2.min.QTL$circRNA.ID,STA.paper$circRNA.ID),"geneSymbol"]
SNPLabel=with(PGC.scz2.min.QTL.detail,data.frame(V1=X4%&%"/"%&%X4,V2=X4%&%"/"%&%X5,V3=X5%&%"/"%&%X5))

########################LocusZoom prepare#############################
#tail -n +2  /media/data3/GWAS/PGC/SCZ2/ckqny.scz2snpres \ 
# |perl -alne '{print "$F[0]\t$F[4]\t$F[1]\t$F[2]\t$F[3]\t$F[8]"}' \
# >/media/data3/GWAS/PGC/SCZ2/ckqny.scz2snpres.vcf
#bgzip /media/data3/GWAS/PGC/SCZ2/ckqny.scz2snpres.vcf
#tabix  -p vcf /media/data3/GWAS/PGC/SCZ2/ckqny.scz2snpres.vcf.gz
##################################################################
LocusZoom.pos=data.frame(chr="chr"%&%PGC.scz2.min.QTL$CHR,
                         start=PGC.scz2.min.QTL$credibleSNP_POS-250000,
                         end=PGC.scz2.min.QTL$credibleSNP_POS+250000,
                         SNP=PGC.scz2.min.QTL$credibleSNP,
                         P=PGC.scz2.min.QTL$crediblePval,
                         geneStart=CIRI.POS[PGC.scz2.min.QTL$circRNA.ID,"left"],
                         geneEnd=CIRI.POS[PGC.scz2.min.QTL$circRNA.ID,"right"])
write.table(LocusZoom.pos,"./qtl_C/LocusZoom/LocusZoom.pos",sep="\t",row.names = F,col.names=F,quote=F)
################################tabix get VCF######################
# cd /media/data3/circCMC/data/qtl_C/LocusZoom/
# sort -k 1.4,1n -k 2,2n ./LocusZoom.pos >./LocusZoom.sort.pos
#  tabix   /media/data3/GWAS/PGC/SCZ2/ckqny.scz2snpres.vcf.gz -R \
#  ./LocusZoom.sort.pos >./LocusZoom.scz2snpres.vcf
##################################################################
# read LocusZoom vcf
LocusZoom.vcf=read.delim("./qtl_C/LocusZoom/LocusZoom.scz2snpres.vcf",header=F)
colnames(LocusZoom.vcf)=c("chr","pos","MarkerName","ref","alt","P.value")
write.table(LocusZoom.vcf,"./qtl_C/LocusZoom/LocusZoom.scz2snpres.vcf",sep="\t",row.names = F,col.names=T,quote=F)

#############################read T2D GWAS data###################
T2D.GWAS=read.delim("/media/data3/GWAS/DIAGRAM/1000GnoAdjBMI_2017/METAANALYSIS_DIAGRAM_SE1.txt",header=T)
T2D.GWAS.pos=do.call("rbind",strsplit(T2D.GWAS$Chr.Position,":"))
T2D.GWAS.pos=data.frame(T2D.GWAS.pos)
T2D.GWAS.pos$X1="chr"%&%T2D.GWAS.pos$X1
T2D.GWAS.pos$rsid=T2D.GWAS[,1]
T2D.GWAS.pos$ref=T2D.GWAS[,2]
T2D.GWAS.pos$alt=T2D.GWAS[,3]
T2D.GWAS.pos$P=T2D.GWAS[,6]
write.table(T2D.GWAS.pos,"/media/data3/GWAS/DIAGRAM/1000GnoAdjBMI_2017/METAANALYSIS_DIAGRAM_SE1.myZL.vcf",
            col.names =F,row.names=F,sep="\t",quote=F)
rm(T2D.GWAS)
###############################sort vcf########################
#sort -k 1,1 -k 2,2n METAANALYSIS_DIAGRAM_SE1.myZL.vcf >METAANALYSIS_DIAGRAM_SE1.myZL.sort.vcf
#bgzip METAANALYSIS_DIAGRAM_SE1.myZL.sort.vcf
###############################################################
GWAS.info.T2D=GWAS.HLD.inCircRNA.df[grep("22885922",GWAS.HLD.inCircRNA.df$gwas),]
##############################
#chr15_rs952471_77776498_C_G as index SNP beacuse the minimum P-value
# tabix  METAANALYSIS_DIAGRAM_SE1.myZL.sort.vcf.gz chr15:77526498-78026498 \
#  >/media/data3/circCMC/data/qtl_C/LocusZoom/LocusZoom.T2Dsnpres.vcf
################################
LocusZoom.T2Dsnpres=read.delim("/media/data3/circCMC/data/qtl_C/LocusZoom/LocusZoom.T2Dsnpres.vcf",header=F)
LocusZoom.T2Dsnpres[which(LocusZoom.T2Dsnpres$V3=="15:77776498"),'V3']="rs952471"
colnames(LocusZoom.T2Dsnpres)=c("chr","pos","MarkerName","ref","alt","P.value")
write.table(LocusZoom.T2Dsnpres,"/media/data3/circCMC/data/qtl_C/LocusZoom/LocusZoom.T2Dsnpres.vcf",
            col.names = T,row.names = F,sep="\t",quote=F)

###########################inflammatory bowel disease disese##########
GWAS.info.IBD=GWAS.HLD.inCircRNA.df[grep("inflammatory bowel disease",GWAS.HLD.inCircRNA.df$disease),]
IBD.GWAS=read.delim("/media/data3/GWAS/IBD/28067908/ibd_build37_59957_20161107.txt",header=T)
tmp=data.frame(do.call("rbind",strsplit(IBD.GWAS[,1],":")))
tmp1=data.frame(do.call("rbind",strsplit(tmp[,2],"_")))
IBD.GWAS.pos=data.frame(chr=tmp[,1],pos=tmp1[,1],SNP=IBD.GWAS[,1],ref=tmp1[,2],alt=tmp1[,3],p=IBD.GWAS[,6])
IBD.GWAS.pos[which(IBD.GWAS.pos$pos=="197701279"),3]="rs2488397"
colnames(IBD.GWAS.pos)=c("chr","pos","MarkerName","ref","alt","P.value")
IBD.GWAS.pos$chr="chr"%&%IBD.GWAS.pos$chr
write.table(IBD.GWAS.pos,"/media/data3/GWAS/IBD/28067908/ibd_build37_59957_20161107.vcf",
            col.names =F,row.names=F,sep="\t",quote=F)
#choice chr1_rs2224873_197645449_T_A to show
##########################################################
# sort -k 1,1 -k 2,2n ibd_build37_59957_20161107.vcf > ibd_build37_59957_20161107.sort.vcf
# bgzip ibd_build37_59957_20161107.sort.vcf
# tabix -p vcf ibd_build37_59957_20161107.sort.vcf.gz
# tabix ibd_build37_59957_20161107.sort.vcf.gz  chr1:197395449-197895449 \
#  > /media/data3/circCMC/data/qtl_C/LocusZoom/LocusZoom.IBDsnpres.vcf
##########################################################
LocusZoom.IBDsnpres=read.delim("/media/data3/circCMC/data/qtl_C/LocusZoom/LocusZoom.IBDsnpres.vcf",header=F)
colnames(LocusZoom.IBDsnpres)=c("chr","pos","MarkerName","ref","alt","P.value")
LocusZoom.IBDsnpres$MarkerName=LocusZoom.IBDsnpres$chr%&%":"%&%LocusZoom.IBDsnpres$pos
LocusZoom.IBDsnpres[which(LocusZoom.IBDsnpres$pos==197645449),'MarkerName']="rs2224873"

write.table(LocusZoom.IBDsnpres,"/media/data3/circCMC/data/qtl_C/LocusZoom/LocusZoom.IBDsnpres.vcf",
            col.names = T,row.names = F,sep="\t",quote=F)
#######################


###################################create LocusZoom genotype###############
# cd /media/data3/CMC/data/genotype/Caucasian465/SNPmatrix
# grep "chr3_rs2710313_52969529_T_A" chr3.genotype.MAF5.mat >/media/data3/circCMC/data/qtl_C/LocusZoom/LocusZoom.genoytpe.MAF5.mat
target.circRNA="chr3:52952500|52962357"
# read LocusZoom genotype
LocusZoom.genotype=read.table("./qtl_C/LocusZoom/LocusZoom.genoytpe.MAF5.mat",header = F)
LocusZoom.genotype.mat=LocusZoom.genotype[,-1]

fenRA=function(x){
  tmp=apply(abs(matrix(c(x,x,x),nrow = 3,byrow = T)-c(0,1,2)),2,which.min)
  return(tmp-1)
}
tmp1=as.numeric(CIRI.JUNCTION.PASS.voom.log.Res[target.circRNA,CaucasianID])
tmp2=as.numeric(LocusZoom.genotype.mat[1,])
tmp2=fenRA(tmp2)
plot(tmp1~tmp2)
abline(lm(tmp1~tmp2))
tmp3=as.numeric(CIRI.RATIO.IM.Res.Caucasian465[target.circRNA,])
plot(tmp3~tmp2)
abline(lm(tmp3~tmp2))
#violin
{
  idx=which(CIRI.JUNCTION.PASS.mat.Cau465[target.circRNA,]!=0)
  tmp1=as.numeric(CIRI.JUNCTION.PASS.voom.log.Res[target.circRNA,CaucasianID])[idx]
  tmplabel=c("TT","TA","AA")
  tmp2=as.character(tmplabel[fenRA(as.numeric(LocusZoom.genotype.mat[1,]))+1])[idx]
  pdf(file="./pdf_C/QTL_GWAS_example.pdf")
  violin.df=data.frame(value=tmp1,group=tmp2)
  violin.p <- ggplot(violin.df, aes(x=group,y=value),color=group) +
    geom_violin(aes(fill=factor(group)),width=0.5) +
    geom_boxplot(fill="white",width=.2) + theme_pubr()+
    theme(axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5)) +
    theme(legend.position="none")+
    labs(title =paste(STA.paper[target.circRNA,"geneSymbol"] ," ",target.circRNA))+
    xlab(label="rs2710313")+
    ylab(label="Adjusted circRNA expression")+
    theme(plot.title = element_text(hjust = 0.5))
  print(violin.p)
  dev.off()
}
{
  i=2
  idx=which(CIRI.JUNCTION.PASS.mat.Cau465[target.circRNA,]!=0)
  tmp1=as.numeric(CIRI.RATIO.IM.Res.Caucasian465[target.circRNA,])[idx]
  tmp1=as.numeric(CIRI.RATIO.IM[target.circRNA,CaucasianID])[idx]
  
  tmp2=as.character(SNPLabel[i,fenRA(as.numeric(LocusZoom.genotype.mat[1,]))+1])[idx]
  pdf(filename = paste0("./pdf_C/QTL_GWAS_example.pdf"),width=480,height = 360,units="px",bg="transparent")
  violin.df=data.frame(value=tmp1,group=tmp2)
  tmpname=PGC.scz2.min.QTL[i,"circRNA.ID"]
  violin.p <- ggplot(violin.df, aes(x=group,y=value),color=group) +
    geom_violin(aes(fill=factor(group)),width=0.5) +
    geom_boxplot(fill="white",width=.2) + theme_pubr()+
    theme(axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5)) +
    theme(legend.position="none")+
    labs(title =paste(STA.paper[target.circRNA,"geneSymbol"] ," ",tmpname))+
    xlab(label="rs2710313")+
    ylab(label="Adjusted circRNA expression")+
    theme(plot.title = element_text(hjust = 0.5))
  print(violin.p)
  dev.off()
}
#Attention:
#RP5-874C20.3's circQTL SNP frequence is C/C 2  T/C 56  T/T 477, we dont' consider it
#############################VeP all QTL#####################
QTL.CIS.PASS.VCF=with(QTL.CIS.PASS.info,data.frame(chr=X1,pos=X3,rsID=X2,ref=X4,alt=X5))
QTL.CIS.PASS.VCF[which(QTL.CIS.PASS.VCF$alt=="<DEL>"),"alt"]=""
write.table(QTL.CIS.PASS.VCF,"./qtl_C/QTL.CIS.PASS.VCF",col.names=F,row.names=F,sep="\t",quote=F)

#output PGC QTL info
write.table(STA.paper[match(PGC.scz2.min.QTL.detail$circRNA.ID,STA.paper$circRNA.ID),],
            "GWASinfo.paper.txt",sep="\t",col.names = T,row.names=F,quote = F)
#######################QTL.CIS.PASS VEP annotation###########################
# sort -k 1,1 -k 2,2n /media/data3/circCMC/data/qtl_C/QTL.CIS.PASS.VCF > /media/data3/circCMC/data/qtl_C/QTL.CIS.PASS.sort.VCF
# vep -i /media/data3/circCMC/data/qtl_C/QTL.CIS.PASS.sort.VCF \
# -o /media/data3/circCMC/data/qtl_C/VEP/QTL.CIS.PASS.sort.most_severe.VEP \
# -gff /media/data3/annotation/GENCODE/gencode.v19.annotation.gff3.gz \
# -fasta /media/data3/genome/hg19/hg19.fa --most_severe --force_overwrite 2>1.log &
# #######################################################################
# vep -i /media/data3/circCMC/data/qtl_C/QTL.CIS.PASS.sort.VCF \
# -o /media/data3/circCMC/data/qtl_C/VEP/QTL.CIS.PASS.sort.VEP \
# -gff /media/data3/annotation/GENCODE/gencode.v19.annotation.gff3.gz \
# -fasta /media/data3/genome/hg19/hg19.fa --force_overwrite 2>2.log &
# cd ./VEP
# perl -alne '{next  if  /^#/;$gene=substr($F[3],0,15);print "$F[0]\t$gene"}'\
#  QTL.CIS.PASS.sort.VEP |sort -k 1,1 -k 2,2n |uniq >QTL.CIS.PASS.rsid_gene
###########################################################################
QTL.CIS.PASS.sort.most_severe.VEP=read.table("./qtl_C/VEP/QTL.CIS.PASS.sort.most_severe.VEP",header=F)
QTL.CIS.PASS.sort.most_severe.VEP=QTL.CIS.PASS.sort.most_severe.VEP[!duplicated(QTL.CIS.PASS.sort.most_severe.VEP$V1),]
QTL.CIS.PASS.sort.most_severe.VEP=QTL.CIS.PASS.sort.most_severe.VEP[,c(1,7)]
QTL.CIS.PASS.most_severe=merge(QTL.CIS.PASS.sort.most_severe.VEP,cbind(QTL.CIS.PASS.info,QTL.CIS.PASS.distance),by.x="V1",by.y="X2")
QTL.CIS.PASS.most_severe$vepType=VEP.type.change[QTL.CIS.PASS.most_severe$V7]
QTL.CIS.PASS.most_severe.canonical.splice.site=
  subset(QTL.CIS.PASS.most_severe,vepType=="canonical splice site")
QTL.CIS.PASS.most_severe.splice.region.variant=
  subset(QTL.CIS.PASS.most_severe,vepType=="splice region variant")
length(unique(QTL.CIS.PASS.most_severe.canonical.splice.site$V1))
genotypeTmp=read.delim("/media/data3/CMC/data/genotype/Caucasian465/chr17_rs11077396_67175150_T_C.mat",header=F,sep=" ")
genotypeTmp=fenRA(as.numeric(genotypeTmp[,-1]))
boxplot(CIRI.JUNCTION.PASS.voom.log.Res["chr17:67151160|67175149",CaucasianID]~genotypeTmp)
boxplot(as.numeric(CIRI.RATIO.mat.zero.Res["chr17:67151160|67175149",CaucasianID])~genotypeTmp)


QTL.CIS.PASS.most_severe=QTL.CIS.PASS.most_severe[!duplicated(QTL.CIS.PASS.most_severe$SNP),]
QTL.CIS.PASS.rsid_gene=read.table("./qtl_C/VEP/QTL.CIS.PASS.rsid_gene",sep="\t",header = F)
#rsid_gene 's gene separate by
QTL.CIS.PASS.rsid_gene.rs=rep("",length = length(unique(QTL.CIS.PASS.rsid_gene$V1)))
names(QTL.CIS.PASS.rsid_gene.rs)=unique(QTL.CIS.PASS.rsid_gene$V1)
for (i in 1:nrow(QTL.CIS.PASS.rsid_gene)) {
  if(QTL.CIS.PASS.rsid_gene.rs[QTL.CIS.PASS.rsid_gene[i,1]]==""){
    QTL.CIS.PASS.rsid_gene.rs[QTL.CIS.PASS.rsid_gene[i,1]]=QTL.CIS.PASS.rsid_gene[i,2]
  }else{
    QTL.CIS.PASS.rsid_gene.rs[QTL.CIS.PASS.rsid_gene[i,1]]=QTL.CIS.PASS.rsid_gene.rs[QTL.CIS.PASS.rsid_gene[i,1]]%&%
      ","%&%QTL.CIS.PASS.rsid_gene[i,2]
  }
}
QTL.CIS.PASS.SNP.info=data.frame(QTL.CIS.PASS.most_severe[,c("SNP","X4","X5","V7")],gene=QTL.CIS.PASS.rsid_gene.rs[QTL.CIS.PASS.most_severe$V1])
colnames(QTL.CIS.PASS.SNP.info)=c("SNP","Ref","Alt","region","SNPgene")
write.table(QTL.CIS.PASS.SNP.info,"./QTL.CIS.PASS.SNP.info",col.names = T,row.names = F,sep="\t",quote=F)
QTL.CIS.PASS.SNP.info.private=merge(QTL.coloc,QTL.CIS.PASS.most_severe,by.x="V2",by.y="V1",all.x = T)
QCPSIP.canonical_splice_variant=subset(QTL.CIS.PASS.SNP.info.private,V7=="splice_donor_variant" | V7=="splice_acceptor_variant")

# CMC.sQTL=read.table("/media/data3/CMC/data/qtl/CMC_MSSM-Penn-Pitt_DLPFC_Analysis_isoQTL-adjustedSVA.txt.gz",header=T)
# CMC.sQTL.cis=subset(CMC.sQTL,eQTL_type=="cis")
# 
# QCPSIP.canonical_splice_variant.sQTL=merge(QCPSIP.canonical_splice_variant,CMC.sQTL.cis,by.x="V2",by.y="SNP")

canonical.splice.site.detail=merge(QTL.CIS.PASS.most_severe.canonical.splice.site,
                                                 QTL.CIS.PASS.SNP.info,by="SNP")
canonical.splice.site.detail.pass=subset(canonical.splice.site.detail,circRNA.gene==SNPgene&nchar(circRNA.gene)==15)
canonical.splice.site.detail.pass=canonical.splice.site.detail.pass[!duplicated(canonical.splice.site.detail.pass$SNP),]
write.table(canonical.splice.site.detail.pass,"./qtl_C/splicingSite/canonical.splice.site.detail.pass.txt",row.names = F,col.names = F,sep="\t",quote = F)
canonical.splice.site.detail.pass.choice=subset(canonical.splice.site.detail.pass,d2==1)

#####Read all genotype pos########
genotype.pos=read.delim("/media/data3/CMC/data/genotype/Caucasian465/chrALL.genotype.MAF5.info80.pos",header = F)
genotype.pos.key=genotype.pos$V2%&%":"%&%genotype.pos$V3
STA.pos.key=data.frame(matrix(nrow = nrow(STA.pos),ncol=4))
STA.pos.key$X1=STA.pos[,'chr']%&%":"%&%c(STA.pos[,'left']-1)
STA.pos.key$X2=STA.pos[,'chr']%&%":"%&%c(STA.pos[,'right']+1)
STA.pos.key$X3=STA.pos[,'chr']%&%":"%&%c(STA.pos[,'left']-2)
STA.pos.key$X4=STA.pos[,'chr']%&%":"%&%c(STA.pos[,'right']+2)

QTL.CIS.PASS.key=QTL.CIS.PASS.VCF$chr%&%":"%&%QTL.CIS.PASS.VCF$pos
genotype.BS=intersect(c(STA.pos.key$X1,STA.pos.key$X2),genotype.pos.key)
genotype.BSinQTL=intersect(genotype.BS,QTL.CIS.PASS.key)
QTL.cISS.PASS.inBS=QTL.CIS.PASS[which(QTL.CIS.PASS.key%in%genotype.BSinQTL),]
QTL.cISS.PASS.inBS$SNPpos=as.numeric(QTL.CIS.PASS.VCF[rownames(QTL.cISS.PASS.inBS),'pos'])
QTL.cISS.PASS.inBS$circLeft=STA.pos[QTL.cISS.PASS.inBS$gene,'left']
QTL.cISS.PASS.inBS$circRight=STA.pos[QTL.cISS.PASS.inBS$gene,'right']
QTL.cISS.PASS.inBS$disMin=apply(data.frame(left=abs(c(QTL.cISS.PASS.inBS$SNPpos-QTL.cISS.PASS.inBS$circLeft)),
           right=abs(c(QTL.cISS.PASS.inBS$SNPpos-QTL.cISS.PASS.inBS$circRight))),1,min)

QTL.cISS.PASS.inBS=subset(QTL.cISS.PASS.inBS,disMin==1)
write.table(QTL.cISS.PASS.inBS,"./qtl_C/splicingSite/QTL.cISS.PASS.inBS.txt",row.names = F,col.names = F,sep="\t",quote = F)

############figures for canonical splice site##########
#violin
{
  target.circRNA="chr1:65455118|65460771"
  target.genotype=as.numeric(read.table("/media/data3/CMC/data/genotype/Caucasian465/chr1_rs17127197_65460772_C_G.mat",header = F)[,-1])
  target.SNPLabel=c("CC","CG","GG")
  #idx=which(CIRI.JUNCTION.PASS.mat.Cau465[target.circRNA,]!=0)
  tmp1=as.numeric(CIRI.JUNCTION.PASS.voom.log.Res[target.circRNA,CaucasianID])
  tmp2=as.character(target.SNPLabel[fenRA(as.numeric(target.genotype))+1])
  pdf(file =  "./pdf_C/canonical_QTL_1.pdf")
  violin.df=data.frame(value=tmp1,group=tmp2)
  violin.p <- ggplot(violin.df, aes(x=group,y=value),color=group) +
    geom_violin(aes(fill=factor(group)),width=0.6) +
    geom_boxplot(fill="white",width=.2) + theme_pubr()+
    theme(axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5)) +
    theme(legend.position="none")+
    labs(title =paste(STA.paper[target.circRNA,"geneSymbol"] ," ",target.circRNA))+
    xlab(label="rs17127197")+
    ylab(label="Adjusted circRNA expression")+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_fill_manual(values = ggplot2::alpha(c("#F8766D","#00BA38","#619CFF"),1))
  print(violin.p)
  dev.off()
}
{
  target.circRNA="chr17:67151160|67175149"
  target.genotype=as.numeric(read.table("/media/data3/CMC/data/genotype/Caucasian465/chr17_rs11077396_67175150_T_C.mat",header = F)[,-1])
  target.SNPLabel=c("TT","TC","CC")
  tmp1=as.numeric(CIRI.JUNCTION.PASS.voom.log.Res[target.circRNA,CaucasianID])
  tmp2=as.character(target.SNPLabel[fenRA(as.numeric(target.genotype))+1])
  pdf(file =  "./pdf_C/canonical_QTL_2.pdf")
  violin.df=data.frame(value=tmp1,group=tmp2)
  violin.p <- ggplot(violin.df, aes(x=group,y=value),color=group) +
    geom_violin(aes(fill=factor(group)),width=0.8) +
    geom_boxplot(fill="white",width=.2) + theme_pubr()+
    theme(axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5)) +
    theme(legend.position="none")+
    labs(title =paste(STA.paper[target.circRNA,"geneSymbol"] ," ",target.circRNA))+
    xlab(label="rs11077396")+
    ylab(label="Adjusted circRNA expression")+
    theme(plot.title = element_text(hjust = 0.5))
  print(violin.p)
  dev.off()
}
{
  target.circRNA="chr8:124089351|124154696"
  target.genotype=as.numeric(read.table("/media/data3/CMC/data/genotype/Caucasian465/chr8_rs10101626_124154697_G_T.mat",header = F)[,-1])
  target.SNPLabel=c("GG","GT","TT")
  tmp1=as.numeric(CIRI.JUNCTION.PASS.voom.log.Res[target.circRNA,CaucasianID])
  tmp2=as.character(target.SNPLabel[fenRA(as.numeric(target.genotype))+1])
  pdf(file =  "./pdf_C/canonical_QTL_3.pdf")
  violin.df=data.frame(value=tmp1,group=tmp2)
  violin.p <- ggplot(violin.df, aes(x=group,y=value),color=group) +
    geom_violin(aes(fill=factor(group)),width=0.8) +
    geom_boxplot(fill="white",width=.2) + theme_pubr()+
    theme(axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5)) +
    theme(legend.position="none")+
    labs(title =paste(STA.paper[target.circRNA,"geneSymbol"] ," ",target.circRNA))+
    xlab(label="rs10101626")+
    ylab(label="Adjusted circRNA expression")+
    theme(plot.title = element_text(hjust = 0.5))
  print(violin.p)
  dev.off()
}
{
  target.circRNA="chr3:27215930|27326452"
  target.genotype=as.numeric(read.table("/media/data3/CMC/data/genotype/Caucasian465/chr3_rs1550769_27296796_C_G.mat",header = F)[,-1])
  target.SNPLabel=c("CC","CG","GG")
  tmp1=as.numeric(CIRI.JUNCTION.PASS.voom.log.Res[target.circRNA,CaucasianID])
  tmp2=as.character(target.SNPLabel[fenRA(as.numeric(target.genotype))+1])
  pdf(file =  "./pdf_C/internal_canonical_QTL_1.pdf")
  violin.df=data.frame(value=tmp1,group=tmp2)
  violin.p <- ggplot(violin.df, aes(x=group,y=value),color=group) +
    geom_violin(aes(fill=factor(group)),width=0.8) +
    geom_boxplot(fill="white",width=.2) + theme_pubr()+
    theme(axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5)) +
    theme(legend.position="none")+
    labs(title =paste(STA.paper[target.circRNA,"geneSymbol"] ," ",target.circRNA))+
    xlab(label="rs1550769")+
    ylab(label="Adjusted circRNA expression")+
    theme(plot.title = element_text(hjust = 0.5))
  print(violin.p)
  dev.off()
}
{
  target.circRNA="chr14:88883055|88904247"
  target.genotype=as.numeric(read.table("/media/data3/CMC/data/genotype/Caucasian465/chr14_rs3179969_88862529_G_A.mat",header = F)[,-1])
  target.SNPLabel=c("1GG","2GA","3AA")
  tmp1=as.numeric(CIRI.JUNCTION.PASS.voom.log.Res[target.circRNA,CaucasianID])
  tmp2=as.character(target.SNPLabel[fenRA(as.numeric(target.genotype))+1])
  pdf(file =  "./pdf_C/distant_canonical_QTL_1.pdf")
  violin.df=data.frame(value=tmp1,group=tmp2)
  violin.p <- ggplot(violin.df, aes(x=group,y=value),color=group) +
    geom_violin(aes(fill=factor(group)),width=0.8) +
    geom_boxplot(fill="white",width=.2) + theme_pubr()+
    theme(axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5)) +
    theme(legend.position="none")+
    labs(title =paste(STA.paper[target.circRNA,"geneSymbol"] ," ",target.circRNA))+
    xlab(label="rs3179969")+
    ylab(label="Adjusted circRNA expression")+
    theme(plot.title = element_text(hjust = 0.5))
  print(violin.p)
  dev.off()
}
#################figures for splice region variant###################
#QTL.CIS.PASS.SNP.info.private.detail=merge(QTL.private.CIS.Detail,QTL.CIS.PASS.SNP.info,by="SNP")
QCPSIP.splice_region_variant=subset(QTL.CIS.PASS.SNP.info.private,V7=="splice_region_variant")
# {
#   target.circRNA="chr2:45925080|45928630"
#   target.genotype=as.numeric(read.table("/media/data3/CMC/data/genotype/Caucasian465/chr2_rs563601_45928630_G_A.mat",header = F)[,-1])
#   target.SNPLabel=c("GG","GA","AA")
#   tmp1=as.numeric(CIRI.JUNCTION.PASS.voom.log.Res[target.circRNA,CaucasianID])
#   tmp2=as.character(target.SNPLabel[fenRA(as.numeric(target.genotype))+1])
#   pdf(file =  "./pdf_C/splice_region_QTL_1.pdf")
#   violin.df=data.frame(value=tmp1,group=tmp2)
#   violin.p <- ggplot(violin.df, aes(x=group,y=value),color=group) +
#     geom_violin(aes(fill=factor(group)),width=0.6) +
#     geom_boxplot(fill="white",width=.2) + theme_pubr()+
#     theme(axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5)) +
#     theme(legend.position="none")+
#     labs(title =paste(STA.paper[target.circRNA,"geneSymbol"] ," ",target.circRNA))+
#     xlab(label="rs11077396")+
#     ylab(label="Adjusted circRNA expression")+
#     theme(plot.title = element_text(hjust = 0.5))
#   print(violin.p)
#   dev.off()
# }
# "PITT_RNA_PFC_1433"#GG max
# "MSSM_RNA_PFC_316"#GA
# "MSSM_RNA_PFC_119"#GA
{
  target.circRNA="chr8:18656805|18662408"
  target.genotype=as.numeric(read.table("/media/data3/CMC/data/genotype/Caucasian465/chr8_rs2306825_18656805_G_A.mat",header = F)[,-1])
  target.SNPLabel=c("GG","GA","AA")
  tmp1=as.numeric(CIRI.JUNCTION.PASS.voom.log.Res[target.circRNA,CaucasianID])
  tmp2=as.character(target.SNPLabel[fenRA(as.numeric(target.genotype))+1])
  pdf(file =  "./pdf_C/splice_region_QTL_2.pdf")
  violin.df=data.frame(value=tmp1,group=tmp2)
  violin.p <- ggplot(violin.df, aes(x=group,y=value),color=group) +
    geom_violin(aes(fill=factor(group)),width=0.8) +
    geom_boxplot(fill="white",width=.2) + theme_pubr()+
    theme(axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5)) +
    theme(legend.position="none")+
    labs(title =paste(STA.paper[target.circRNA,"geneSymbol"] ," ",target.circRNA))+
    xlab(label="rs2306825")+
    ylab(label="Adjusted circRNA expression")+
    theme(plot.title = element_text(hjust = 0.5))
  print(violin.p)
  dev.off()
}
##########################plot chr2:45925080|45928630 density#################
target.SNPLabel=c("GG","GA","AA")
target.circRNA="chr8:18656805|18662408"
myBAM=c("PENN_RNA_PFC_28","PENN_RNA_PFC_45","PENN_RNA_PFC_59","PITT_RNA_PFC_10003","PITT_RNA_PFC_1270")
target.genotype=as.numeric(read.table("/media/data3/CMC/data/genotype/Caucasian465/chr2_rs563601_45928630_G_A.mat",header = F)[,-1])
target.genotype=as.character(target.SNPLabel[fenRA(as.numeric(target.genotype))+1])
names(target.genotype)=CaucasianID
target.genotype[myBAM]
CIRI.JUNCTION.PASS.voom.log.Res[target.circRNA,myBAM]
CIRI.JUNCTION.PASS.mat[target.circRNA,myBAM]
target.genotype[myBAM]
#"PITT_RNA_PFC_1454" #GA
#"MSSM_RNA_PFC_127" #GA
geneID="ENSG00000156011"
geneID.transcript=gtfV19.transcript.detail[which(gtfV19.transcript.detail$geneID==geneID),]
geneID.detail=subset(gtfV19.gene,gene==geneID)
region=geneID.detail[,1]%&%":"%&%geneID.detail[,4]%&%"-"%&%geneID.detail[,5]
samplePool=c("PENN_RNA_PFC_45","PITT_RNA_PFC_1454","PITT_RNA_PFC_1270")
for (j in samplePool) {
  baseDir="/media/data3/CMC/data/fastq/"
  sample=baseDir%&%j%&%".bam"
  out="/media/data3/circCMC/data/isoform/"%&%j%&%"."%&%geneID%&%".depth"
  command=paste("samtools depth -a -r",region,sample,">",out,sep=" ")
  try(system(command))
}
samplePool=c("PENN_RNA_PFC_45","PITT_RNA_PFC_1454","PITT_RNA_PFC_1270")
genePlot(gene = geneID,sample ="PENN_RNA_PFC_45" ,circID = target.circRNA,CircRegion = F)
genePlot(gene = geneID,sample ="PITT_RNA_PFC_1454" ,circID = target.circRNA,CircRegion = T)
genePlot(gene = geneID,sample ="PITT_RNA_PFC_1270" ,circID = target.circRNA,CircRegion = F)


######################manhattan plots for min-circQTL#######################
tmp=do.call("rbind",strsplit(QTL.private.CIS$SNP,"\\_"))
QTL.manhattan.CIS=data.frame(SNP=tmp[,2],CHR=as.numeric(substr(tmp[,1],4,100)),
                             BP=as.numeric(tmp[,3]),P=QTL.private.CIS$pvalue)
#manhattan(QTL.manhattan.CIS,col=colorRampPalette(c("navy", "lightblue","lightgreen","yellow", "orange","firebrick3"),bias=2)(22))
{
  pdf("./pdf_C/SupFig3d.pdf")
  manhattan(QTL.manhattan.CIS,suggestiveline = F,genomewideline = F)
  dev.off()
}
# {
#   pdf("./pdf_C/SupFig3a.pdf")
#   par(mar=c(5,5,1,1))
#   hist(QTL.PERM.Empirical,col = "grey",xlab="Empirical P-values",main="")
#   dev.off()
#   pdf("./pdf/SupFig3b.pdf")
#   par(mar=c(5,5,1,1))
#   hist(QTL.PERM.Empirical.qvalue$qvalues,col = "grey",xlab="Q-values",main="")
#   text(0.025,2790,"2,735",cex=1)
#   dev.off()
# }

#############################SupFig3c###########################
color3=c("#d53e4f","#7fcdbb","#3288bd")
names(color3)=c("splice5","Internal","splice3")
QTL.private.CIS.distance.Fig3=QTL.private.CIS.distance
QTL.private.CIS.distance.Fig3$distance[which(QTL.private.CIS.distance.Fig3$type=="Internal")]=
  50000*QTL.private.CIS.distance.Fig3$distance[which(QTL.private.CIS.distance.Fig3$type=="Internal")]
QTL.private.CIS.distance.Fig3$distance[which(QTL.private.CIS.distance.Fig3$type=="splice3")]=
  50000+QTL.private.CIS.distance.Fig3$distance[which(QTL.private.CIS.distance.Fig3$type=="splice3")]
{
  pdf("./pdf_C/SupFig3c.pdf")
  plot(-log(QTL.private.CIS$pvalue,10)~QTL.private.CIS.distance.Fig3$distance,pch = 19,cex=0.4,las=1,
       col = color3[QTL.private.CIS.distance.Fig3$type],
       ylab=expression(paste(-log[10]," ",italic(P),"-value")),xlab="Distance to the 5' or 3' splice site (bp)",main="",xaxt="n")
  axis(side=1,at=c(-100000,-50000,0,50000,100000,150000),labels = c(-100,-50,0,1,50,100),las=1)
  abline(v = 0,col="navy",lty=5)
  abline(v = 50000,col="navy",lty=5)
  par(new=TRUE)
  plot(density(QTL.private.CIS.distance.Fig3$distance,bw=700),ylab="",xlab="",main="",col="firebrick3",axes=F,lwd=2)
  legend("topleft",legend = c("5' (N = 1,410)",
                              "Internal (N = 548)",
                              "3' (N = 1,159)"),
         fill=c("#d53e4f","#7fcdbb","#3288bd"),
         col=c("#d53e4f","#7fcdbb","#3288bd"),
         y.intersp =1,x.intersp = 0.5,text.width=1,
         bty="n")
  dev.off()
}

##########Find DE circRNA with genetic risk gene related with QTL#####
DEC.gwas.qtl=intersect(DE.annotation[intersect(unique(QTL.CIS.PASS$gene),DE.annotation$`circRNA ID`),'geneSymbol'],intersect(DEgene,unlist(CMC.geneSet.list[["GWAS:PGC2 SCZ"]])))
DEC.gq.ID=DE.annotation[match(DEC.gwas.qtl,DE.annotation$geneSymbol),'circRNA ID']
DEC.qtl=QTL.CIS.PASS[which(QTL.CIS.PASS$gene==DEC.gq.ID),]
DEC.qtl.SNP=do.call("rbind",strsplit(unique(DEC.qtl$SNP),"_"))
DEC.qtl$rsID=DEC.qtl.SNP[,2]
pgc.scz2.credible.snps.old=read.delim("/media/data3/GWAS/PGC/SCZ2/pgc.scz2.credible.snps.txt")
DEC.qtl=merge(DEC.qtl,pgc.scz2.credible.snps.old,by.x="rsID",by.y="credibleSNP")

#############################
# date: 11/1/2018
# author: Zelin
# decription: revised paper
###########################
############Experimental validataion###############
circBase=read.table("/media/data3/circRNA/script/hsa_hg19_circRNA.bed",header=F)
rownames(circBase)=with(circBase,V1%&%":"%&%(V2+1)%&%"|"%&%V3)
#read circBase detail
circBase.detial=read.delim("/media/data3/annotation/circRNA/hsa_hg19_circRNA.txt",header=T)

MolCell=read.delim("/media/data3/annotation/circRNA/molecularCell.txt")
rownames(MolCell)=with(MolCell,chr%&%":"%&%(start+1)%&%"|"%&%end)

CIRIinDatabase=intersect(rownames(STA.paper),unique(c(rownames(circBase),rownames(MolCell))))
length(CIRIinDatabase)/10559
CIRI.not.Database=setdiff(rownames(STA.paper),unique(c(rownames(circBase),rownames(MolCell))))
STA.paper.notDB=STA.paper[CIRI.not.Database,]
#STA.notDB.Rand20=STA.paper.notDB[sample(1:nrow(STA.paper.notDB),20),]
CIRI.JUNCTION.noDB.mat=CIRI.JUNCTION.mat[CIRI.not.Database,]
CIRI.JUNCTION.noDB.nonZero=rowSums(CIRI.JUNCTION.noDB.mat>0)
CIRI.JUNCTION.noDB.Exp=rowMeans(CIRI.JUNCTION.noDB.mat)
length(which(CIRI.JUNCTION.noDB.Exp>30 & CIRI.JUNCTION.noDB.nonZero>500))
STA.paper.notDB$Exp=CIRI.JUNCTION.noDB.Exp
STA.paper.notDB$nonZero=CIRI.JUNCTION.noDB.nonZero
#
STA.paper.notDB.Genenot=subset(STA.paper.notDB,!(geneSymbol%in%circBase.detial$gene.symbol))
STA.paper.notDB.Genenot=STA.paper.notDB.Genenot[!is.na(STA.paper.notDB.Genenot$geneSymbol),]
STA.paper.notDB.Genenot=subset(STA.paper.notDB.Genenot,geneSymbol%in%names(which(table(STA.paper.notDB.Genenot$geneSymbol)==1)))
STA.paper.notDB.Genenot=subset(STA.paper.notDB.Genenot,region=="exon")
STA.paper.notDB.Genenot.cb=cbind(STA.paper.notDB.Genenot,circRNA.existinCellorTissue)
STA.paper.notDB.Genenot.cb$oneNum=apply(circRNA.existinCellorTissue,1,function(x)sum(x==1,na.rm = T))
STA.paper.notDB.Genenot.cb=STA.paper.notDB.Genenot.cb[order(STA.paper.notDB.Genenot.cb$Exp,decreasing = T),]
STA.paper.notDB.Genenot.cb$rseqLen=0
STA.paper.notDB.Genenot.cb$lseqLen=0

for (i in 1:nrow(STA.paper.notDB.Genenot.cb)) {
  id=STA.paper.notDB.Genenot.cb[i,1]
  chr=unlist(strsplit(id,":"))[1]
  exon.Start=unlist(strsplit(STA.paper.notDB.Genenot.cb[i,'exon.start.site'],","))
  exon.End=unlist(strsplit(STA.paper.notDB.Genenot.cb[i,'exon.end.site'],","))
  this.circ.strand=STA.paper.notDB.Genenot.cb[i,2]
  exon.Num=STA.paper.notDB.Genenot.cb[i,'exon.number']
  if(exon.Num>1){
    STA.paper.notDB.Genenot.cb$lseqLen[i]=as.numeric(exon.End[1])-as.numeric(exon.Start[1])
    STA.paper.notDB.Genenot.cb$rseqLen[i]=as.numeric(exon.End[exon.Num])-as.numeric(exon.Start[exon.Num])
  }
}
STA.paper.notDB.choice=subset(STA.paper.notDB.Genenot.cb,rseqLen>60&lseqLen>60)

write.table(STA.paper.notDB.choice,"Primer/STA.paper.notDB.choice.txt",row.names = F,col.names = T,sep="\t",quote=F,na = "")

STA.notDB.Choice20=STA.paper.notDB[which(CIRI.JUNCTION.noDB.Exp>30 & CIRI.JUNCTION.noDB.nonZero>500),]
STA.notDB.Choice20.nonDup=STA.notDB.Choice20[!duplicated(STA.notDB.Choice20$gene),]
STA.notDB.Choice20.final=rbind(STA.notDB.Choice20.nonDup,STA.notDB.Choice20["chr12:91069070|91073722",])
STA.notDB.Choice20.final=STA.notDB.Choice20.final[order(STA.notDB.Choice20.final$predicted.length),]
STA.notDB.Choice20.final$exp=CIRI.JUNCTION.noDB.Exp[STA.notDB.Choice20.final$circRNA.ID]
CIRI.JUNCTION.noDB.Exp.name=names(sort(CIRI.JUNCTION.noDB.Exp,decreasing = T))
validate13.list=c("chr6:117010483|117013555",
                  "chr12:1372200|1481143",
                  "chr14:47389215|47504507",
                  "chr2:162728803|162757521",
                  "chr7:14517786|14712652",
                  "chr14:71880665|71996087",
                  "chr13:101997652|102051516",
                  "chr10:116879949|117001514",
                  "chr7:137148235|137294355",
                  "chr1:243708812|243859018",
                  "chr8:42761316|42873637",
                  "chr2:50699436|50847321",
                  "chr6:16586011|16754322"
)
validate13.list.Exp=CIRI.JUNCTION.noDB.Exp[validate13.list]
order(validate13.list.Exp,decreasing = T)
for (i in 1:length(validate13.list)) {
  which(CIRI.JUNCTION.noDB.Exp.name==validate13.list[i])
}
write.table(STA.notDB.Choice20.final,"Primer/STA.notDB.Choice20.txt",row.names = F,col.names = F,sep="\t",quote=F)
#intergenic
STA.paper.notDB.Intergenic=subset(STA.paper.notDB,region=="intergenic_region")
STA.paper.notDB.Intergenic$nonZero=CIRI.JUNCTION.noDB.nonZero[rownames(STA.paper.notDB.Intergenic)]
STA.paper.notDB.Intergenic$Exp=CIRI.JUNCTION.noDB.Exp[rownames(STA.paper.notDB.Intergenic)]
#intron
STA.paper.notDB.Intron=subset(STA.paper.notDB,region=="intron")
STA.paper.notDB.Intron$nonZero=CIRI.JUNCTION.noDB.nonZero[rownames(STA.paper.notDB.Intron)]
STA.paper.notDB.Intron$Exp=CIRI.JUNCTION.noDB.Exp[rownames(STA.paper.notDB.Intron)]

# after manually check, we filter out circRNA with too short exon junction
circRNAexperVal=c("chr6:117010483|117013555",
                  "chr12:1372200|1481143","chr14:47389215|47504507",
                  "chr2:162728803|162757521","chr7:14517786|14712652",
                  "chr14:71880665|71996087","chr13:101997652|102051516",
                  "chr10:116879949|117001514","chr7:137148235|137294355",
                  "chr1:243708812|243859018","chr8:42761316|42873637",
                  "chr2:50699436|50847321","chr6:16586011|16754322")

# write.table(STA.notDB.Rand20[order(STA.notDB.Rand20$predicted.length),],"Primer/STA.notDB.Rand20.txt",row.names = F,col.names = F,sep="\t",quote=F)
getSeq=function(chr,start,end,strand){
  baseAlp=c("A","T","G","C")
  names(baseAlp)=c("T","A","C","G")
  cmd="samtools faidx /media/data3/genome/hg19/hg19.fa "%&%chr%&%":"%&%start%&%"-"%&%end
  each=toupper(paste(strsplit(as.character(system(cmd,intern=T)),"\n")[-1],collapse=""))
  if(strand=="-"){
    each=paste(baseAlp[rev(unlist(strsplit(each,"")))],collapse="")
  }
  return(each)
}

cat("",file="/media/data3/circCMC/data/Primer/design2.exon.Seq")
for (i in 1:nrow(STA.paper.notDB.choice)) {
  id=STA.paper.notDB.choice[i,1]
  chr=unlist(strsplit(id,":"))[1]
  exon.Start=unlist(strsplit(STA.paper.notDB.choice[i,'exon.start.site'],","))
  exon.End=unlist(strsplit(STA.paper.notDB.choice[i,'exon.end.site'],","))
  this.circ.strand=STA.paper.notDB.choice[i,2]
  if(STA.paper.notDB.choice[i,'region']=="exon"){
    exon.Num=STA.paper.notDB.choice[i,'exon.number']
    if(exon.Num>1){
      lseq=getSeq(chr,exon.Start[1],exon.End[1],this.circ.strand)
      lseq.len=nchar(lseq)
      rseq=getSeq(chr,exon.Start[exon.Num],exon.End[exon.Num],this.circ.strand)
      rseq.len=nchar(rseq)
      if(this.circ.strand=="-"){
        junctionSeq=substr(lseq,nchar(lseq)-9,nchar(lseq))%&%"|"%&%substr(rseq,1,10)
        targetSeq=">"%&%id%&%" length:"%&%lseq.len%&%"|"%&%rseq.len%&%" "%&%junctionSeq%&%"\n"%&%lseq%&%"|"%&%rseq%&%"\n"
      }else{
        junctionSeq=substr(rseq,nchar(rseq)-9,nchar(rseq))%&%"|"%&%substr(lseq,1,10)
        targetSeq=">"%&%id%&%" length:"%&%rseq.len%&%"|"%&%lseq.len%&%" "%&%junctionSeq%&%"\n"%&%rseq%&%"|"%&%lseq%&%"\n"
      }
      cat(targetSeq,file="/media/data3/circCMC/data/Primer/design2.exon.Seq",append=T)
    }
  }
}
for (i in 1:nrow(STA.uniqueCIRI)) {
  id=STA.uniqueCIRI[i,1]
  chr=unlist(strsplit(id,":"))[1]
  exon.Start=unlist(strsplit(STA.uniqueCIRI[i,'exon.start.site'],","))
  exon.End=unlist(strsplit(STA.uniqueCIRI[i,'exon.end.site'],","))
  this.circ.strand=STA.uniqueCIRI[i,2]
  if(STA.uniqueCIRI[i,'region']=="exon"){
    exon.Num=STA.uniqueCIRI[i,'exon.number']
    if(exon.Num>1){
      lseq=getSeq(chr,exon.Start[1],exon.End[1],this.circ.strand)
      lseq.len=nchar(lseq)
      rseq=getSeq(chr,exon.Start[exon.Num],exon.End[exon.Num],this.circ.strand)
      rseq.len=nchar(rseq)
      if(this.circ.strand=="-"){
        junctionSeq=substr(lseq,nchar(lseq)-19,nchar(lseq))%&%substr(rseq,1,20)
        print(id%&%" 1|"%&%(lseq.len-20)%&%"|"%&%(lseq.len+21)%&%"|"%&%(lseq.len+rseq.len)%&%" "%&%junctionSeq)
      }else{
        junctionSeq=substr(rseq,nchar(rseq)-19,nchar(rseq))%&%substr(lseq,1,20)
        print(id%&%" 1|"%&%(rseq.len-20)%&%"|"%&%(rseq.len+21)%&%"|"%&%(lseq.len+rseq.len)%&%" "%&%junctionSeq)
        
      }
    }
  }
  
}
tmp=CIRI.JUNCTION.mat[STA.uniqueCIRI$circRNA.ID,]
#####################Different method comparing############
#read detected circRNAs from CIRI
EACH.SAMPLE.circRNA.NUM=read.table("./compareMethods/CIRI.circRNAinEach.txt",header=F)
hist(EACH.SAMPLE.circRNA.NUM$V1)

#read findCIRI detected circRNA statistic
findCIRI.count=read.table("./compareMethods/findCIRI.circID.count.txt",header = F)
rownames(findCIRI.count)=findCIRI.count[,2]
findCIRI.inCIRI.count=findCIRI.count[STA.paper$circRNA.ID,]
rownames(findCIRI.inCIRI.count)=STA.paper$circRNA.ID
findCIRI.inCIRI.count$V1[is.na(findCIRI.inCIRI.count$V1)]=0
{
  pdf("./pdf_C/find_circDist.pdf")
  hist(findCIRI.inCIRI.count$V1,breaks = 25,col = "grey",ylab="Number of circRNAs",xlab="Number of samples")
  dev.off()
  
}
sum(findCIRI.inCIRI.count$V1>=100)
sum(findCIRI.inCIRI.count$V1>=100)/10559

CIRI.JUNCTION.SampleNum.inSTA=CIRI.JUNCTION.SampleNum[STA.paper$circRNA.ID]
plot(findCIRI.inCIRI.count$V1~CIRI.JUNCTION.SampleNum.inSTA)
findCIRI.inCIRI.nonDB=findCIRI.inCIRI.count[CIRI.not.Database,]
findCIRI.inCIRI.nonDB.STA=STA.paper[rownames(findCIRI.inCIRI.nonDB[which(findCIRI.inCIRI.nonDB$V1==0),]),]
findCirc.zero=rownames(findCIRI.inCIRI.nonDB)[which(findCIRI.inCIRI.nonDB$V1==0)]



#read CIRCexplorer2 detected circRNA statistic
CIRC.count=read.table("./compareMethods/CIRCexplorer2.circID.parse.count.10599.txt",header = F)
rownames(CIRC.count)=CIRC.count[,2]
CIRC.count=CIRC.count[STA.paper$circRNA.ID,]
CIRC.count$V1[is.na(CIRC.count$V1)]=0
{
  pdf("./pdf_C/CIRCDist.pdf")
  hist(CIRC.count$V1,breaks = 25,col = "grey",ylab="Number of circRNAs",xlab="Number of samples")
  dev.off()
}
sum(CIRC.count$V1>=100)
sum(CIRC.count$V1>=100)/10559

rownames(CIRC.count)=STA.paper$circRNA.ID
CIRC.inCIRI.nonDB=CIRC.count[CIRI.not.Database,]
CIRC.zero=rownames(CIRC.inCIRI.nonDB)[which(CIRC.inCIRI.nonDB$V1==0)]

STA.uniqueCIRI=STA.paper[intersect(findCirc.zero,CIRC.zero),]
(10559-length(intersect(rownames(CIRC.count)[which(CIRC.count$V1==0)],
                 rownames(findCIRI.inCIRI.count)[which(findCIRI.inCIRI.count$V1==0)])))/10559
combfindcircACIRC=data.frame(findcirc=findCIRI.inCIRI.count$V1,CIRC=CIRC.count$V1)
write.table(combfindcircACIRC,"./combfindcircACIRC.txt",
            row.names = F,col.names = F,quote = F,sep="\t")

#figure S1D
findCirc.pass.ID=rownames(subset(findCIRI.inCIRI.count,V1>=100))
CIRC.pass.ID=rownames(subset(CIRC.count,V1>=100))
findCircUCIRC=unique(c(findCirc.pass.ID,CIRC.pass.ID))
findCirc.CIRC=intersect(findCirc.pass.ID,CIRC.pass.ID)

(length(findCirc.pass.ID)-length(findCirc.CIRC))/10559
(length(CIRC.pass.ID)-length(findCirc.CIRC))/10559
length(findCirc.CIRC)/10559
length(findCircUCIRC)/10559

########################cal circRNA/total RNA fraction#######################
#cal gene length
gencodeV19exon=read.delim("/media/data3/annotation/GENCODE/gencode.v19.annotation.exon.sort.txt",header = F)
gencodeV19exon.geneID=unique(gencodeV19exon$V4)
gencodeV19exon.geneLength=vector(length =length(gencodeV19exon.geneID))
for(i in 1:length(gencodeV19exon.geneID)){
  print(i)
  tmpexon=subset(gencodeV19exon,V4==gencodeV19exon.geneID[i])
  tmpPos=c()
  for(j in 1:nrow(tmpexon)){
    tmpPos=c(tmpPos,tmpexon[j,2]:tmpexon[j,3])
  }
  gencodeV19exon.geneLength[i]=length(unique(tmpPos))
}
names(gencodeV19exon.geneLength)=gencodeV19exon.geneID
tmp=intersect(gencodeV19exon.geneID,names(CMC.exon2gene.length))
plot(gencodeV19exon.geneLength[tmp]~CMC.exon2gene.length[tmp])
#read the exon info
CMC.exon.length=read.delim("/media/data3/CMC/data/expression/exonResult.txt",header = F)
CMC.exon2gene.length=tapply(CMC.exon.length$V6, factor(CMC.exon.length$V1), sum)
geneIngtfID=intersect(names(CMC.exon2gene.length),rownames(Tophat2.raw.counts))
CMC.exon2gene.sub.length=CMC.exon2gene.length[geneIngtfID]
gtfV19Comp=intersect(setdiff(gencodeV19exon.geneID,
                   intersect(intersect(gencodeV19exon.geneID,rownames(Tophat2.raw.counts)),geneIngtfID)),
                   rownames(Tophat2.raw.counts))
geneIngtfID=c(geneIngtfID,gtfV19Comp)
CMC.exon2gene.sub.length=c(CMC.exon2gene.sub.length,gencodeV19exon.geneLength[gtfV19Comp])
#read raw gene counts
Tophat2.raw.counts.Inget=Tophat2.raw.counts[geneIngtfID,]
geneIngtf.depth=colSums(Tophat2.raw.counts.Inget)
CMC.Gene.RPKM=10^9*Tophat2.raw.counts.Inget/CMC.exon2gene.sub.length/matrix(rep(Tophat2.raw.counts.LIB_SIZE+CIRI.JUNCTION.depth,length(geneIngtfID)),nrow=length(geneIngtfID),byrow=T)
CMC.Gene.RPKM.rowMean=rowMeans(CMC.Gene.RPKM,na.rm = T)
CMC.Gene.RPKM.colSum=colSums(CMC.Gene.RPKM,na.rm = T)
sum(CMC.Gene.RPKM.rowMean>0.1,na.rm = T)
hist(CMC.Gene.RPKM.colSum)
# sum of circRNA RPKM
# circRNA RPKM = junction reads/(circRNA length * total mapped reads)
# circRNA length = (length of reads-20bp) * 2
CIRI.JUNCTION.sumRPKM=10^9*CIRI.JUNCTION.depth/160/(Tophat2.raw.counts.LIB_SIZE+CIRI.JUNCTION.depth)
circ2linear=CIRI.JUNCTION.sumRPKM/(CMC.Gene.RPKM.colSum+CIRI.JUNCTION.sumRPKM)
hist(circ2linear[which(covariates$Institution=="MSSM"&covariates$libBatch=="B")],breaks=20)
hist(circ2linear,breaks = 25)
mean(circ2linear)
which.max(circ2linear)

# another method to cal circ2linear
plot(CIRI.Linear.depth~CIRI.JUNCTION.depth)
circ2linear.ratio=CIRI.JUNCTION.depth/CIRI.Linear.depth
plot(sort(circ2linear.ratio))
{
  pdf("./pdf_C/ratioDepthC2L.pdf")
  hist(circ2linear.ratio,breaks = 40,col="grey",main="",ylab = "Number of samples",
       xlab="Ratio of circular to linear junction reads")
  abline(v = mean(circ2linear.ratio),col="blue",lty=5,cex=2)
  dev.off()
}
range(circ2linear.ratio)
median(circ2linear.ratio)
summary(lm(CIRI.eachSample.NumCirc~circ2linear.ratio))

# circular Ratio * gene counts
CIRI.ANNOTATION.nonNA=subset(CIRI.ANNOTATION,V3!="n/a")
CIRI.ANNOTATION.nonNA$V3=substr(CIRI.ANNOTATION.nonNA$V3,1,15)
CIRI.ANNOTATION.nonNA.sub=subset(CIRI.ANNOTATION.nonNA,V3%in%geneIngtfID)
CIRI.ANNOTATION.nonNA.sub=data.frame(CIRI.ANNOTATION.nonNA.sub)
rownames(CIRI.ANNOTATION.nonNA.sub)=CIRI.ANNOTATION.nonNA.sub[,1]
CIRI.RATIO.inGeneId=which(rownames(CIRI.RATIO)%in%CIRI.ANNOTATION.nonNA.sub$V1)
CIRI.RATIO.inGene.mat=CIRI.RATIO[CIRI.RATIO.inGeneId,-1]
geneID=CIRI.ANNOTATION.nonNA.sub[rownames(CIRI.RATIO.inGene.mat),'V3']
CMC.Gene.RPKM.gene2circLong=CMC.Gene.RPKM[geneID,]
CIRI.RATIO2GeneRPKM.inGene.mat=CIRI.RATIO.inGene.mat*CMC.Gene.RPKM.gene2circLong
CIRI.RATIO2GeneRPKM.inGene.colsum=colSums(CIRI.RATIO2GeneRPKM.inGene.mat,na.rm = T)
CIRI.GeneRPKM.ratio=CIRI.RATIO2GeneRPKM.inGene.colsum/CMC.Gene.RPKM.colSum
plot(CIRI.GeneRPKM.ratio)
range(CIRI.GeneRPKM.ratio)
median(CIRI.GeneRPKM.ratio)
mean(CIRI.GeneRPKM.ratio)
hist(CIRI.GeneRPKM.ratio)

#CIRI.JUNCTION.depth/Tophat2.raw.counts.LIB_SIZE
CIRI.JUN2Gene.Ratio=CIRI.JUNCTION.depth/Tophat2.raw.counts.LIB_SIZE
CIRI.JUN2Gene.Ratio.M1=CIRI.JUNCTION.depth.M1/Tophat2.raw.counts.LIB_SIZE
CIRI.JUN2Gene.Ratio.M2=CIRI.JUNCTION.depth.M2/Tophat2.raw.counts.LIB_SIZE
CIRI.JUN2Gene.Ratio.M3=CIRI.JUNCTION.depth.M3/Tophat2.raw.counts.LIB_SIZE
CIRI.JUN2Gene.Ratio.M4=CIRI.JUNCTION.depth.M4/Tophat2.raw.counts.LIB_SIZE
CIRI.JUN2Gene.Ratio.M5=CIRI.JUNCTION.depth.M5/Tophat2.raw.counts.LIB_SIZE
plot(sort(CIRI.JUN2Gene.Ratio))
boxplot(CIRI.JUN2Gene.Ratio)
median(CIRI.JUN2Gene.Ratio)
range(CIRI.JUN2Gene.Ratio)
hist(CIRI.JUN2Gene.Ratio,breaks = 100)
mean(CIRI.JUNCTION.depth/(Tophat2.raw.counts.LIB_SIZE+CIRI.JUNCTION.depth))

######################circRNA in different cell lines###########
for(i in  list.files("./encode/CIRI/")){
  assign(i,read.delim("./encode/CIRI/"%&%i,header = T))
  assign(i%&%".circ2linear",sum(get(i)$X.junction_reads)/sum(get(i)$X.non_junction_reads))
}
for(i in  list.files("./encode/CIRI/")){
  assign(i%&%".innonDBCIRI",intersect(get(i)$circRNA_ID,STA.paper.notDB$circRNA.ID))
  print(i)
  print(length(get(i%&%".innonDBCIRI")))
}
circRNAexperVal=STA.paper.notDB.Genenot$circRNA.ID

circRNA.existinCellorTissue=data.frame(matrix(nrow=length(circRNAexperVal),ncol=length(list.files("./encode/CIRI/"))))
rownames(circRNA.existinCellorTissue)=circRNAexperVal
n=0
for(i in list.files("./encode/CIRI/")){
  n=n+1
  circRNA.existinCellorTissue[which(circRNAexperVal%in%get(i)$circRNA_ID),n]=1
}
colnames(circRNA.existinCellorTissue)=tissueID[do.call("rbind",strsplit(list.files("./encode/CIRI/"),"\\."))[,1]]
#circRNA.existinCellorTissue$circExp=CIRI.JUNCTION.noDB.Exp[circRNAexperVal]
write.table(circRNA.existinCellorTissue,"circRNA.existinCellorTissue.txt",na = "",row.names = T,col.names = T,sep="\t",quote = F)

cellLine.circ2linear=vector(length = length(list.files("./encode/")))
names(cellLine.circ2linear)=list.files("./encode/")
for(i in list.files("./encode/")){
  cellLine.circ2linear[i]=get(i%&%".circ2linear")
}
for (i in list.files("./encode/")){
  print(i)
  print(intersect(STA.paper.notDB$circRNA.ID,get(i)[,1]))
}
######################circRNA in different tisses###########
# for i in `ls *`;do echo -n ${i%%\.*}" ";tail -n +5 $i |cut -f 2 |perl -alne '{$a=$a+$F[0]}END{print $a}';done  >../all.genecounts.txt
# for i in `ls *`;do echo -n ${i%%\.*}" ";tail -n +2 $i |cut -f 5 |perl -alne '{$a=$a+$F[0]}END{print $a}';done >../CJ.counts.txt
# for i in `ls *`;do echo -n ${i%%\.*}" ";tail -n +2 $i |cut -f 7 | perl -alne '{$a=$a+$F[0]}END{print $a}';done >../LJ.counts.txt
# for i in `ls *`;do echo -n ${i%%\.*}" ";tail -n +2 $i |cut -f 5 |perl -alne '{if($F[0]>1){$a=$a+$F[0]}}END{print $a}';done >../CJ.counts.more1.txt
# for i in `ls *`;do echo -n ${i%%\.*}" ";tail -n +2 $i |cut -f 5 |perl -alne '{if($F[0]>2){$a=$a+$F[0]}}END{print $a}';done >../CJ.counts.more2.txt
# for i in `ls *`;do echo -n ${i%%\.*}" ";tail -n +2 $i |cut -f 5 |perl -alne '{if($F[0]>3){$a=$a+$F[0]}}END{print $a}';done >../CJ.counts.more3.txt
# for i in `ls *`;do echo -n ${i%%\.*}" ";tail -n +2 $i |cut -f 5 |perl -alne '{if($F[0]>4){$a=$a+$F[0]}}END{print $a}';done >../CJ.counts.more4.txt
# for i in `ls *`;do echo -n ${i%%\.*}" ";tail -n +2 $i |cut -f 5 |perl -alne '{if($F[0]>5){$a=$a+$F[0]}}END{print $a}';done >../CJ.counts.more5.txt
tissueID=c("stomach", "thyroid gland female adult" ,"esophagus muscularis mucosa female adult" ,
           "gastrocnemius medialis","sigmoid colon", "spleen","uterus","pancreas" ,
           "heart left ventricle" , "ovary", "right atrium auricular region",
           "right lobe of liver" , "ascending aorta" , "vagina" ,
           "cerebellum","frontal cortex", "K562","dendritic cell" ,
           "Peyer's patch","adrenal gland","esophagus squamous epithelium",
           "gastroesophageal sphincter","lower leg skin female","omental fat pad",
           "subcutaneous adipose","suprapubic skin","tibial nerve","transverse colon",
           "upper lobe of left lung","breast epithelium","diencephalon","K562_2",
           "cd34mobilizedCell","Haoaf", "Haoec", "Hch00113082","Hfdpc",
           "Hmepc", "Hmncpb","Hmscat010260412","Hpiepc","Hsavec", "Hvmf", "Imr90"
)
names(tissueID)=c("SRR4422373","SRR4421529","SRR4421314","SRR4422107","SRR4421313", "SRR4422445",
                  "SRR4421350","SRR4422592","SRR4421747","SRR4422347","SRR4421966",
                  "SRR4421506","SRR4422192","SRR4421394","SRR3192427","SRR3192425",
                  "SRR4422162","SRR5210444","SRR4421545","SRR4422339","SRR4421678",
                  "SRR4422217","SRR4422571","SRR4422098","SRR4421875","SRR4421845",
                  "SRR4421566","SRR4422057","SRR4421631","SRR4422692","SRR3192432", "SRR3192410",
                  "cd34mobilizedCell","Haoaf", "Haoec", "Hch00113082","Hfdpc",
                  "Hmepc", "Hmncpb","Hmscat010260412","Hpiepc","Hsavec", "Hvmf", "Imr90"
)
cellLineID=c("K562","dendritic cell" ,"cd34mobilizedCell","Haoaf", "K562_2",
             "Haoec", "Hch00113082","Hfdpc", "Hmepc", "Hmncpb",
             "Hmscat010260412","Hpiepc","Hsavec", "Hvmf", "Imr90")

tissue.CJ=read.table("./encode/CJ.counts.txt",header = F)
tissue.m1.CJ=read.table("./encode/CJ.counts.more1.txt",header = F)
tissue.m2.CJ=read.table("./encode/CJ.counts.more2.txt",header = F)
tissue.m3.CJ=read.table("./encode/CJ.counts.more3.txt",header = F)
tissue.m4.CJ=read.table("./encode/CJ.counts.more4.txt",header = F)
tissue.m5.CJ=read.table("./encode/CJ.counts.more5.txt",header = F)
tissue.m9.CJ=read.table("./encode/CJ.counts.more9.txt",header = F)
tissue.mat.CJ=matrix(c(tissue.CJ$V2,tissue.m1.CJ$V2,tissue.m2.CJ$V2,tissue.m3.CJ$V2,
                       tissue.m4.CJ$V2,tissue.m5.CJ$V2,tissue.m9.CJ$V2),nrow=nrow(tissue.CJ),ncol=7)

rownames(tissue.CJ)=tissue.CJ$V1
#tissue.LJ=read.table("./encode/LJ.counts.txt",header = F)
#rownames(tissue.LJ)=tissue.LJ$V1
tissue.GeneCounts=read.table("./encode/all.genecounts.txt",header = F)
rownames(tissue.GeneCounts)=tissue.GeneCounts$V1
tissue.GeneCounts=tissue.GeneCounts[rownames(tissue.CJ),]
#tissue.CJ2LJ=tissue.CJ$V2/tissue.LJ$V2
tissue.CJ2Gene=tissue.CJ$V2/tissue.GeneCounts$V2

tissueID=tissueID[rownames(tissue.GeneCounts)]
tissue.CJ$tissue=tissueID
tissue.GeneCounts$tissue=tissueID
#tissue.CJ2LJ.df=data.frame(ratio=tissue.CJ2LJ,tissue=tissueID)
tissue.CJ2Gene.df=data.frame(ratio=tissue.CJ2Gene,tissue=tissueID)
tissue.all.df=data.frame(CJ=tissue.CJ$V2,
                         geneCounts=tissue.GeneCounts$V2,
                         CJ2Gene=tissue.CJ2Gene.df$ratio)
tissue.all.mat.df=cbind(tissue.all.df,tissue.mat.CJ/tissue.all.df$geneCounts)
tissue.all.mat.df=tissue.all.mat.df[-which(tissue.CJ$tissue=="K562"),]
tissue.all.mat.df=tissue.all.mat.df[order(tissue.all.mat.df$`1`),]
tissue.all.mat.df$Num=1:nrow(tissue.all.mat.df)

tissue.all.df$tissue=tissueID
tissue.all.df$type="tissue"
tissue.all.df[which(tissue.all.df$tissue%in%cellLineID),"type"]="cell line"
tissue.all.df=tissue.all.df[order(tissue.all.df$CJ2Gene),]
tissue.all.df=tissue.all.df[-which(tissue.all.df$tissue=="K562"),]
tissue.all.mat.df$tissue=tissue.all.df$tissue
tissue.all.mat.df$type=tissue.all.df$type
write.table(tissue.all.df,"encode/tissue.all.df.txt",
            col.names = T,row.names = F,sep="\t",quote = F)

write.table(tissue.all.mat.df,"encode/tissue.all.mat.df.txt",
            col.names = T,row.names = F,sep="\t",quote = F)
tmpColor=c("#f0027f","#66c2a5")
CIRI.JUN2Gene.Ratio.vec.Q2=100*c(as.numeric(quantile(CIRI.JUN2Gene.Ratio)[3]),as.numeric(quantile(CIRI.JUN2Gene.Ratio.M1)[3]),
                                 as.numeric(quantile(CIRI.JUN2Gene.Ratio.M2)[3]),as.numeric(quantile(CIRI.JUN2Gene.Ratio.M3)[3]),
                                 as.numeric(quantile(CIRI.JUN2Gene.Ratio.M4)[3]),as.numeric(quantile(CIRI.JUN2Gene.Ratio.M5)[3]))
CIRI.JUN2Gene.Ratio.vec.Q1=100*c(as.numeric(quantile(CIRI.JUN2Gene.Ratio)[2]),as.numeric(quantile(CIRI.JUN2Gene.Ratio.M1)[2]),
                             as.numeric(quantile(CIRI.JUN2Gene.Ratio.M2)[2]),as.numeric(quantile(CIRI.JUN2Gene.Ratio.M3)[2]),
                             as.numeric(quantile(CIRI.JUN2Gene.Ratio.M4)[2]),as.numeric(quantile(CIRI.JUN2Gene.Ratio.M5)[2]))
CIRI.JUN2Gene.Ratio.vec.Q3=100*c(as.numeric(quantile(CIRI.JUN2Gene.Ratio)[4]),as.numeric(quantile(CIRI.JUN2Gene.Ratio.M1)[4]),
                             as.numeric(quantile(CIRI.JUN2Gene.Ratio.M2)[4]),as.numeric(quantile(CIRI.JUN2Gene.Ratio.M3)[4]),
                             as.numeric(quantile(CIRI.JUN2Gene.Ratio.M4)[4]),as.numeric(quantile(CIRI.JUN2Gene.Ratio.M5)[4]))

for (i in 1:6) {
  {
    pdf("./pdf_C/revised_RelativeCircRank."%&%i%&%".pdf")
    tmp1=c(1:43,44)
    tmp2=c(tissue.all.mat.df[,i+5]*100,CIRI.JUN2Gene.Ratio.vec.Q2[i])
    plot(tmp2~tmp1,pch=c(c(15,16)[as.factor(tissue.all.df$type)],17),ylim=c(0,1.25),
         ylab="Relative circRNA junction counts (%)",las=1,
         col=c(tmpColor[as.factor(tissue.all.df$type)],"#619CFF"),xaxt="n",xlab="")
    lines(c(43.5,44.5),c(CIRI.JUN2Gene.Ratio.vec.Q1[i],CIRI.JUN2Gene.Ratio.vec.Q1[i]))
    lines(c(43.5,44.5),c(CIRI.JUN2Gene.Ratio.vec.Q3[i],CIRI.JUN2Gene.Ratio.vec.Q3[i]))
    lines(c(44,44),c(CIRI.JUN2Gene.Ratio.vec.Q3[i],CIRI.JUN2Gene.Ratio.vec.Q1[i]))
    legend(1,1.25,legend = c("14 cell lines",
                             "29 tissues","DLPFC"),
  
           col=c("#f0027f","#66c2a5","#619CFF"),pch = c(15,16,17),
           y.intersp =1,x.intersp = 0.5,text.width=2,
           bty="n")
    dev.off()
  }
}

##########################Encode frontal cortex##########
CIRI.Frontal.Cortex=read.delim("./encode/CIRI/SRR3192425.txt")
length(intersect(CIRI.Frontal.Cortex$circRNA_ID,STA.paper$circRNA.ID))
length(intersect(CIRI.Frontal.Cortex$circRNA_ID,STA.paper.notDB$circRNA.ID))

#how many circRNA higher than their host genes?
CIRI.RATIO.MoreThanlinear.NUM=colSums(CIRI.RATIO.mat>0.5)
range(CIRI.RATIO.MoreThanlinear.NUM)
median(CIRI.RATIO.MoreThanlinear.NUM)
range(CIRI.RATIO.MoreThanlinear.NUM/CIRI.eachSample.NumCirc)
median(CIRI.RATIO.MoreThanlinear.NUM/CIRI.eachSample.NumCirc)
CIRI.RATIO.rowMean=rowMeans(CIRI.RATIO.mat,na.rm = T)
range(CIRI.RATIO.rowMean)
sum(CIRI.RATIO.rowMean>0.5)
sum(CIRI.RATIO.rowMean>0.5)/10559
STA.ratio05=STA.paper[names(which(CIRI.RATIO.rowMean[rownames(STA.paper.ingene)]>0.5)),]
write.table(unlist(strsplit(unique(STA.ratio05$gene),",")),"./DAVID//ratio05.gene.txt",row.names = F,col.names = F,quote = F)
{
  pdf("./pdf_C/circRatioFreq.pdf")
  hist(CIRI.RATIO.rowMean,breaks = 40,col="grey",main="",ylab = "Number of circRNAs",xaxt="n",las=2,
       xlab="Circular ratio")
  axis(1,at =seq(0,1,by = 0.1),labels = seq(0,1,by = 0.1),pos = 0)
  abline(v = 0.5,col="blue",lty=5,lwd=1.5)
  dev.off()
}

# number of circRNAs for each samples
CIRI.eachSample.NumCirc=colSums(CIRI.JUNCTION.mat>0)
{
  pdf("./pdf_C/eachSampleNumCirc.pdf")
  hist(CIRI.eachSample.NumCirc,breaks = 20,col="grey",xlab="Number of circRNAs",ylab="Number of samples",main="")
  dev.off()
}
median(CIRI.eachSample.NumCirc)
mean(CIRI.eachSample.NumCirc)


# transform number to fraction in figure 1H
tmp=table(table(STA.paper$gene[which(STA.paper$gene!="n/a")]))*100/sum(table(table(STA.paper$gene[which(STA.paper$gene!="n/a")])))
sum(tmp[7:length(tmp)])
#############SCZ-Control differential analysis############
SAMPLE.SCZ.ID=rownames(subset(covariates.SCZ,Dx=="SCZ"))
SAMPLE.Control.ID=rownames(subset(covariates.SCZ,Dx=="Control"))
rand100.DE.top.adjP=vector(length = 1000)
rand100.DE.Pass.circ=vector()
set.seed(2018)
for (i in 1:1000) {
  print(i)
  choice.SCZ.ID=sample(SAMPLE.SCZ.ID,100)
  choice.Control.ID=sample(SAMPLE.Control.ID,100)
  choiceID=c(choice.SCZ.ID,choice.Control.ID)
  CIRI.JUNCTION.PASS.SCZ.rand100.voom=voom(CIRI.JUNCTION.PASS.DEG[,choiceID,keep.lib.sizes=T],design=Circ.design.SCZwithControl[choiceID,],plot=F)
  CIRI.JUNCTION.PASS.SCZ.rand100.voom.lmFit=lmFit(CIRI.JUNCTION.PASS.SCZ.rand100.voom,design=Circ.design.SCZwithControl[choiceID,])
  CIRI.JUNCTION.PASS.SCZ.rand100.voom.eBayes=eBayes(CIRI.JUNCTION.PASS.SCZ.rand100.voom.lmFit)
  CIRI.JUNCTION.PASS.SCZ.rand100.voom.DGE=subset(toptable(CIRI.JUNCTION.PASS.SCZ.rand100.voom.eBayes,
                                                   coef="DxSCZ",n=10559,confint=T,sort.by = "p"),adj.P.Val<0.05)
  if(nrow(CIRI.JUNCTION.PASS.SCZ.rand100.voom.DGE)>0){
    rand100.DE.Pass.circ=c(rand100.DE.Pass.circ,rownames(CIRI.JUNCTION.PASS.SCZ.rand100.voom.DGE))
  }
  rand100.DE.top.adjP[i]=CIRI.JUNCTION.PASS.SCZ.rand100.voom.DGE[1,6]
}
{
  pdf("./pdf_C/FDR1000RS.pdf")
  hist(rand100.DE.top.adjP,breaks = 50,main="Histogram of Top FDR in 1000 times random sampling",
       xlab="FDR",ylab="Frequency",cex.lab=1.5,cex.axis=1.5,las=1,cex.main=1.5)
  dev.off()
}
{
  pdf("./pdf_C/rand100DEPasstable.pdf")
  hist(table(rand100.DE.Pass.circ),breaks = 10,main="",
       xlab="Numbers of recurrence for potential DE circRNA",ylab="Frequency",cex.lab=1.5,cex.axis=1.5,las=1,cex.main=1.5)
  dev.off()
}
{
  png("./pic/rand100DEPasstable.png")
  hist(table(rand100.DE.Pass.circ),breaks = 10,main="",
       xlab="Numbers of recurrence for potential DE circRNA",ylab="Frequency",cex.lab=1.5,cex.axis=1.5,las=1,cex.main=1.5)
  dev.off()
}
rand50.DE.top.adjP=vector(length = 1000)
rand50.DE.Pass.circ=vector()

set.seed(2018)
for (i in 1:1000) {
  print(i)
  choice.SCZ.ID=sample(SAMPLE.SCZ.ID,50)
  choice.Control.ID=sample(SAMPLE.Control.ID,50)
  choiceID=c(choice.SCZ.ID,choice.Control.ID)
  CIRI.JUNCTION.PASS.SCZ.rand50.voom=voom(CIRI.JUNCTION.PASS.DEG[,choiceID,keep.lib.sizes=T],design=Circ.design.SCZwithControl[choiceID,],plot=F)
  CIRI.JUNCTION.PASS.SCZ.rand50.voom.lmFit=lmFit(CIRI.JUNCTION.PASS.SCZ.rand50.voom,design=Circ.design.SCZwithControl[choiceID,])
  CIRI.JUNCTION.PASS.SCZ.rand50.voom.eBayes=eBayes(CIRI.JUNCTION.PASS.SCZ.rand50.voom.lmFit)
  CIRI.JUNCTION.PASS.SCZ.rand50.voom.DGE=toptable(CIRI.JUNCTION.PASS.SCZ.rand50.voom.eBayes,
                                                   coef="DxSCZ",n=10559,confint=T,sort.by = "p")
  if(nrow(CIRI.JUNCTION.PASS.SCZ.rand50.voom.DGE)>0){
    rand50.DE.Pass.circ=c(rand50.DE.Pass.circ,rownames(CIRI.JUNCTION.PASS.SCZ.rand50.voom.DGE))
  }
  rand50.DE.top.adjP[i]=CIRI.JUNCTION.PASS.SCZ.rand50.voom.DGE[1,6]
}
sum(rand50.DE.top.adjP<0.05)
# DE analysis according to gender
maleId=covariates.SCZ$id[!covariates.SCZ$Sex=="Female"]
CIRI.JUNCTION.PASS.SCZ.Accmale.voom=voom(CIRI.JUNCTION.PASS.DEG[,maleId,keep.lib.sizes=T],design=Circ.design.SCZwithControl[maleId,],plot=F)
CIRI.JUNCTION.PASS.SCZ.Accmale.voom.lmFit=lmFit(CIRI.JUNCTION.PASS.SCZ.Accmale.voom,design=Circ.design.SCZwithControl[maleId,])
CIRI.JUNCTION.PASS.SCZ.Accmale.voom.eBayes=eBayes(CIRI.JUNCTION.PASS.SCZ.Accmale.voom.lmFit)
CIRI.JUNCTION.PASS.SCZ.Accmale.voom.DGE=toptable(CIRI.JUNCTION.PASS.SCZ.Accmale.voom.eBayes,
                                                coef="DxSCZ",n=Inf,confint=T,sort.by = "p")
femaleId=covariates.SCZ$id[covariates.SCZ$Sex=="Female"]
CIRI.JUNCTION.PASS.SCZ.AccFemale.voom=voom(CIRI.JUNCTION.PASS.DEG[,femaleId,keep.lib.sizes=T],design=Circ.design.SCZwithControl[femaleId,],plot=F)
CIRI.JUNCTION.PASS.SCZ.AccFemale.voom.lmFit=lmFit(CIRI.JUNCTION.PASS.SCZ.AccFemale.voom,design=Circ.design.SCZwithControl[femaleId,])
CIRI.JUNCTION.PASS.SCZ.AccFemale.voom.eBayes=eBayes(CIRI.JUNCTION.PASS.SCZ.AccFemale.voom.lmFit)
CIRI.JUNCTION.PASS.SCZ.AccFemale.voom.DGE=toptable(CIRI.JUNCTION.PASS.SCZ.AccFemale.voom.eBayes,
                                                   coef="DxSCZ",n=Inf,confint=T,sort.by = "p")

# DE analysis according to LIB
#LIB B
libBid=covariates.SCZ$id[covariates.SCZ$libBatch=="B"]
CIRI.JUNCTION.PASS.SCZ.libBid.voom=voom(CIRI.JUNCTION.PASS.DEG[,libBid,keep.lib.sizes=T],design=Circ.design.SCZwithControl[libBid,],plot=F)
CIRI.JUNCTION.PASS.SCZ.libBid.voom.lmFit=lmFit(CIRI.JUNCTION.PASS.SCZ.libBid.voom,design=Circ.design.SCZwithControl[libBid,])
CIRI.JUNCTION.PASS.SCZ.libBid.voom.eBayes=eBayes(CIRI.JUNCTION.PASS.SCZ.libBid.voom.lmFit)
CIRI.JUNCTION.PASS.SCZ.libBid.voom.DGE=toptable(CIRI.JUNCTION.PASS.SCZ.libBid.voom.eBayes,
                                                 coef="DxSCZ",n=Inf,confint=T,sort.by = "p")
#LIB D
libDid=covariates.SCZ$id[covariates.SCZ$libBatch=="D"]
CIRI.JUNCTION.PASS.SCZ.libDid.voom=voom(CIRI.JUNCTION.PASS.DEG[,libDid,keep.lib.sizes=T],design=Circ.design.SCZwithControl[libDid,],plot=F)
CIRI.JUNCTION.PASS.SCZ.libDid.voom.lmFit=lmFit(CIRI.JUNCTION.PASS.SCZ.libDid.voom,design=Circ.design.SCZwithControl[libDid,])
CIRI.JUNCTION.PASS.SCZ.libDid.voom.eBayes=eBayes(CIRI.JUNCTION.PASS.SCZ.libDid.voom.lmFit)
CIRI.JUNCTION.PASS.SCZ.libDid.voom.DGE=toptable(CIRI.JUNCTION.PASS.SCZ.libDid.voom.eBayes,
                                                coef="DxSCZ",n=Inf,confint=T,sort.by = "p")
#LIB E
#LIB D
libEid=covariates.SCZ$id[covariates.SCZ$libBatch=="E"]
CIRI.JUNCTION.PASS.SCZ.libEid.voom=voom(CIRI.JUNCTION.PASS.DEG[,libEid,keep.lib.sizes=T],design=Circ.design.SCZwithControl[libEid,],plot=F)
CIRI.JUNCTION.PASS.SCZ.libEid.voom.lmFit=lmFit(CIRI.JUNCTION.PASS.SCZ.libEid.voom,design=Circ.design.SCZwithControl[libEid,])
CIRI.JUNCTION.PASS.SCZ.libEid.voom.eBayes=eBayes(CIRI.JUNCTION.PASS.SCZ.libEid.voom.lmFit)
CIRI.JUNCTION.PASS.SCZ.libEid.voom.DGE=toptable(CIRI.JUNCTION.PASS.SCZ.libEid.voom.eBayes,
                                                coef="DxSCZ",n=Inf,confint=T,sort.by = "p")



# DE ratio
CIRI.RATIO.10559.SCZAControl.mat=CIRI.RATIO.10559.mat[,covariates.SCZ$id]
Ratio.DE.P=vector(length = 10559)
for (i in 1:10559) {
  tmpRatio=as.numeric(CIRI.RATIO.10559.SCZAControl.mat[i,])
  index=which(tmpRatio==0)
  if(length(index)>0){
    tmpRatio=tmpRatio[-index]
    tmpDesign=Circ.design.SCZwithControl[-index,]
  }else{
    tmpDesign=Circ.design.SCZwithControl
    
  }
  tmplm=summary(lm(tmpRatio~-1+ tmpDesign))
  Ratio.DE.P[i]=tmplm$coefficients[2,4]
}
Ratio.DE.P.adj=p.adjust(Ratio.DE.P,"fdr")
Circ.design.SCZwithControl

# try edgeR
Circ.design.SCZwithControl
CIRI.JUNCTION.PASS.DEG.edgeR=DGEList(CIRI.JUNCTION.PASS.mat[,covariates.SCZ$id],
                               lib.size = CIRI.JUNCTION.depth[covariates.SCZ$id],
                               genes=rownames(CIRI.JUNCTION.PASS.mat),
                               group = covariates.SCZ$Dx)
CIRI.JUNCTION.PASS.DEG.edgeR=calcNormFactors(CIRI.JUNCTION.PASS.DEG.edgeR)
# classical method
CIRI.JUNCTION.PASS.DEG.edgeR.eCD=estimateCommonDisp(CIRI.JUNCTION.PASS.DEG.edgeR,verbose = T)
CIRI.JUNCTION.PASS.DEG.edgeR.eCD=estimateTagwiseDisp(CIRI.JUNCTION.PASS.DEG.edgeR.eCD)
plotBCV(CIRI.JUNCTION.PASS.DEG.edgeR.eCD)
CIRI.JUNCTION.PASS.DEG.edgeR.et=exactTest(CIRI.JUNCTION.PASS.DEG.edgeR.eCD,pair = c("SCZ","Control"))
CIRI.JUNCTION.PASS.DEG.edgeR.top=topTags(CIRI.JUNCTION.PASS.DEG.edgeR.et,n=10559)
# glm method
CIRI.JUNCTION.PASS.DEG.edgeR.eGLMCD=estimateGLMCommonDisp(CIRI.JUNCTION.PASS.DEG.edgeR,Circ.design.SCZwithControl,verbose = T)
CIRI.JUNCTION.PASS.DEG.edgeR.eGLMCD <- estimateGLMTrendedDisp(CIRI.JUNCTION.PASS.DEG.edgeR.eGLMCD, Circ.design.SCZwithControl);  
CIRI.JUNCTION.PASS.DEG.edgeR.eGLMCD <- estimateGLMTagwiseDisp(CIRI.JUNCTION.PASS.DEG.edgeR.eGLMCD,Circ.design.SCZwithControl);
CIRI.JUNCTION.PASS.DEG.edgeR.eGLMCD.fit <- glmFit(CIRI.JUNCTION.PASS.DEG.edgeR.eGLMCD, Circ.design.SCZwithControl);
lrt <- glmLRT(CIRI.JUNCTION.PASS.DEG.edgeR.eGLMCD.fit,coef="DxSCZ"); 
CIRI.JUNCTION.PASS.DEG.edgeR.eGLMCD.top=topTags(lrt,n = 10559); 
CIRI.JUNCTION.PASS.DEG.edgeR.eGLMCD.top.table=CIRI.JUNCTION.PASS.DEG.edgeR.eGLMCD.top$table

######################Flank Intron Position##################
circFlankIntron.pos$len1=circFlankIntron.pos$V4-circFlankIntron.pos$V3
circFlankIntron.pos$len2=circFlankIntron.pos$V6-circFlankIntron.pos$V5
circFlankIntron.pos.exon=subset(circFlankIntron.pos,V1%in%STA.paper.exon$circRNA.ID)
range(circFlankIntron.pos.exon$len1)
range(circFlankIntron.pos.exon$len2)
sum(circFlankIntron.pos.exon$len1<100000&sum(circFlankIntron.pos.exon$len2<100000))

# in back-splicing, internal, and external sites
# circRNA located in back-splcing of the parental gene
QTL.inBS=QTL.cISS.PASS.inBS[,c(1,2,3,5,7,13)]
QTL.inBS$type="back-splicing"
QTL.inBS$circRNA.gene=STA.paper[QTL.inBS$gene,'gene']
QTL.inBS=QTL.inBS[,c(1,2,8,6,7,3,4,5)]
QTL.inBS$qtlPos=c("splice_acceptor","splice_acceptor","splice_acceptor",
                  "splice_acceptor","splice_acceptor","splice_donor","splice_acceptor","splice_donor")
#in internal or external
QTL.inInEx=canonical.splice.site.detail.pass[which(canonical.splice.site.detail.pass$d2!=1),
                                             c(1,3,4,5,6,8,9,16,17)]
QTL.inInEx=QTL.inInEx[,c(1,3,4,8,9,5,6,7,2)]
colnames(QTL.inInEx)=colnames(QTL.inBS)
#combine above
QTL.SplicingSite=data.frame(rbind(QTL.inBS,QTL.inInEx))
QTL.SplicingSite$rsid=do.call("rbind",strsplit(QTL.SplicingSite$SNP,"_"))[,2]
QTL.SplicingSite$key=QTL.SplicingSite$rsid%&%":"%&%QTL.SplicingSite$circRNA.gene
SplicingSitein.sQTL.cis=subset(CMC.sQTL.cis,SNP%in%unique(QTL.SplicingSite$rsid))
SplicingSitein.sQTL.cis$key=SplicingSitein.sQTL.cis$SNP%&%":"%&%SplicingSitein.sQTL.cis$Gene
SplicingSitein.sQTL.cis$Beta=-SplicingSitein.sQTL.cis$Beta
SplicingSitein.eQTL.cis=subset(CMC.eQTLs,SNP%in%unique(QTL.SplicingSite$rsid))
SplicingSitein.eQTL.cis$key=SplicingSitein.eQTL.cis$SNP%&%":"%&%SplicingSitein.eQTL.cis$Gene
SplicingSitein.eQTL.cis$Beta=-SplicingSitein.eQTL.cis$Beta

QTL.SplicingSite=merge(QTL.SplicingSite,SplicingSitein.sQTL.cis[,c(2,3,5,12,16)],by="key",all.x=T)
QTL.SplicingSite=merge(QTL.SplicingSite,SplicingSitein.eQTL.cis[,c(1,2,4,6,16)],by="key",all.x=T)

write.table(QTL.SplicingSite,"./qtl_C/splicingSite/QTL.SplicingSite.txt",row.names = F,col.names = T,quote = F,sep="\t")

#############################make detailed QTL info##################
QTL.CIS.PASS.info.withVEP=merge(QTL.CIS.PASS.info,QTL.CIS.PASS.sort.most_severe.VEP,by.x="X2",by.y="V1")
QTL.CIS.PASS.info.withVEP$rsidGene=QTL.CIS.PASS.info.withVEP$X2%&%"_"%&%QTL.CIS.PASS.info.withVEP$circRNA.gene
QTL.CIS.PASS.info.withVEP.eQTL=merge(QTL.CIS.PASS.info.withVEP,CMC.eQTL.cis[,c("Beta","pvalue","rsidGene")],by="rsidGene",all.x=T)
QTL.CIS.PASS.info.withVEP.eQTL$circLeft=STA.pos[QTL.CIS.PASS.info.withVEP.eQTL$circRNA.ID,'left']
QTL.CIS.PASS.info.withVEP.eQTL$circRight=STA.pos[QTL.CIS.PASS.info.withVEP.eQTL$circRNA.ID,'right']
QTL.CIS.PASS.info.withVEP.eQTL$X3=as.numeric(QTL.CIS.PASS.info.withVEP.eQTL$X3)
QTL.CIS.PASS.info.withVEP.eQTL$region=sign((QTL.CIS.PASS.info.withVEP.eQTL$circLeft-QTL.CIS.PASS.info.withVEP.eQTL$X3))*
  sign((QTL.CIS.PASS.info.withVEP.eQTL$circRight-QTL.CIS.PASS.info.withVEP.eQTL$X3))
QTL.CIS.PASS.info.withVEP.eQTL=merge(QTL.CIS.PASS.info.withVEP.eQTL,circFlankIntron.pos,by.x="circRNA.ID",by.y="V1",all.x=T)
QTL.CIS.PASS.info.withVEP.eQTL$intronL=sign((QTL.CIS.PASS.info.withVEP.eQTL$V3-QTL.CIS.PASS.info.withVEP.eQTL$X3))*
  sign((QTL.CIS.PASS.info.withVEP.eQTL$V4-QTL.CIS.PASS.info.withVEP.eQTL$X3))
QTL.CIS.PASS.info.withVEP.eQTL$intronR=sign((QTL.CIS.PASS.info.withVEP.eQTL$V5-QTL.CIS.PASS.info.withVEP.eQTL$X3))*
  sign((QTL.CIS.PASS.info.withVEP.eQTL$V6-QTL.CIS.PASS.info.withVEP.eQTL$X3))
QTL.CIS.PASS.info.withVEP.eQTL$intronRegion=""
QTL.CIS.PASS.info.withVEP.eQTL$intronRegion[which(QTL.CIS.PASS.info.withVEP.eQTL$intronL==-1)]=-1
QTL.CIS.PASS.info.withVEP.eQTL$intronRegion[which(QTL.CIS.PASS.info.withVEP.eQTL$intronR==-1)]=-1
QTL.CIS.PASS.info.withVEP.eQTL=QTL.CIS.PASS.info.withVEP.eQTL[,c(1,4,5,6,7,8,9,14,15,16,19,27)]
QTL.CIS.PASS.info.withVEP.eQTL[which(QTL.CIS.PASS.info.withVEP.eQTL$region==-1),'region']="circRNA region"
QTL.CIS.PASS.info.withVEP.eQTL[which(QTL.CIS.PASS.info.withVEP.eQTL$intronRegion==-1),'region']="flanking intron"
QTL.CIS.PASS.info.withVEP.eQTL[which(QTL.CIS.PASS.info.withVEP.eQTL$region==0),'region']="circRNA region"
QTL.CIS.PASS.info.withVEP.eQTL[which(QTL.CIS.PASS.info.withVEP.eQTL$region==1),'region']=""
QTL.CIS.PASS.info.withVEP.eQTL=QTL.CIS.PASS.info.withVEP.eQTL[,c(1,2,3,8,11,4,5,6,7,9,10)]
QTL.CIS.PASS.info.withVEP.eQTL$Beta=-QTL.CIS.PASS.info.withVEP.eQTL$Beta
write.table(QTL.CIS.PASS.info.withVEP.eQTL,"QTL.CIS.PASS.info.withVEP.eQTL.txt",
            row.names = F,col.names = T,sep="\t",quote=F,na = "")
QTL.CIS.PASS.info.withVEP.eQTL$key=QTL.CIS.PASS.info.withVEP.eQTL$SNP%&%"-"%&%QTL.CIS.PASS.info.withVEP.eQTL$circRNA.ID
QTL.CIS.PASS.info.withVEP.eQTL=QTL.CIS.PASS.info.withVEP.eQTL[match(QTL.CIS.PASS$key,QTL.CIS.PASS.info.withVEP.eQTL$key),]
QTL.CIS.PASS.info.withVEP.eQTL.private=QTL.CIS.PASS.info.withVEP.eQTL[QTL.private.index$V2,]
##############Where are the QTL SNPs located ##########
GWAS.HLD.inCircRNA.df.paper.Trans=GWAS.HLD.inCircRNA.df.paper
GWAS.HLD.inCircRNA.df.paper.Trans$key=GWAS.HLD.inCircRNA.df.paper.Trans$circRNA.ID%&%":"%&%
  GWAS.HLD.inCircRNA.df.paper.Trans$SNP
tmp=QTL.CIS.PASS.info.withVEP.eQTL
tmp$key=tmp$circRNA.ID%&%":"%&%tmp$SNP
GWAS.HLD.inCircRNA.df.paper.Trans=merge(GWAS.HLD.inCircRNA.df.paper.Trans,
                                        tmp,by="key",all.x = T)
index1=tapply(GWAS.HLD.inCircRNA.df.paper.Trans$pvalue.x,
              as.factor(GWAS.HLD.inCircRNA.df.paper.Trans$circRNA.ID.x),function(x){which.min(x)})
index2=tapply(1:nrow(GWAS.HLD.inCircRNA.df.paper.Trans),
              as.factor(GWAS.HLD.inCircRNA.df.paper.Trans$circRNA.ID.x),function(x){x[1]-1})
index3=index1+index2
GWAS.HLD.inCircRNA.df.paper.Trans.minP=GWAS.HLD.inCircRNA.df.paper.Trans[index3,]
GWAS.HLD.inCircRNA.df.paper.Trans.nonDup=GWAS.HLD.inCircRNA.df.paper.Trans[!duplicated(GWAS.HLD.inCircRNA.df.paper.Trans$SNP),]
table(GWAS.HLD.inCircRNA.df.paper.Trans.nonDup$V7)
table(GWAS.HLD.inCircRNA.df.paper.Trans.nonDup$region)
GWAS.circRNA.inRegion=tapply(GWAS.HLD.inCircRNA.df.paper.Trans$region, as.factor(GWAS.HLD.inCircRNA.df.paper.Trans$circRNA.ID), function(x){sum(x%in%c("flanking intron","circRNA region"))})

##################Splicing#########
{
  target.circRNA="chr8:124089351|124154696"
  target.genotype=as.numeric(read.table("/media/data3/CMC/data/genotype/Caucasian465/chr8_rs10101626_124154697_G_T.mat",header = F)[,-1])
  target.SNPLabel=c("GG","GT","TT")
  expCirc=as.numeric(CIRI.JUNCTION.PASS.voom.log.Res[target.circRNA,CaucasianID])
  tmpgenoytpe=as.character(target.SNPLabel[fenRA(as.numeric(target.genotype))+1])
  tmpGene=STA.paper[target.circRNA,"geneSymbol"]
  tmpEnsembl=STA.paper[target.circRNA,"gene"]
  #expGene=as.numeric(Tophat2.DGEList.voom.log[tmpEnsembl,CaucasianID])
  expGene=as.numeric(CMC.SVA.gene[tmpEnsembl,CaucasianID])
  tmp.df=data.frame(type=rep(c("gene","circRNA"),
                             each=length(expGene)),exp=c(expGene,expCirc),genotype=c(tmpgenoytpe,tmpgenoytpe))
  pdf(file =  "./pdf_C/Re_canonical_QTL_3.pdf",width = 8)
  violin.p <- ggplot(tmp.df,aes(x=genotype,y=exp)) +
    geom_violin(aes(fill=factor(genotype)),width=0.8) +
    geom_boxplot(fill="white",width=.2) + 
    theme_bw(base_size = 12)+
    theme(legend.position="none")+
    labs(title =paste(tmpGene," ",target.circRNA))+
    xlab(label="rs10101626")+
    ylab(label="Adjusted circRNA expression")+
    theme(legend.title = element_text(size = 0))+
    theme(axis.text.y = element_text(size = 24),axis.text.x = element_text(size = 24))+
    facet_wrap(~ type, scales="free")+labs(y = "Adjusted expression level (log2 CPM)")
  print(violin.p)
  dev.off()
}
{
  target.circRNA="chr3:27215930|27326452"
  target.genotype=as.numeric(read.table("/media/data3/CMC/data/genotype/Caucasian465/chr3_rs1550769_27296796_C_G.mat",header = F)[,-1])
  target.SNPLabel=c("CC","CG","GG")
  expCirc=as.numeric(CIRI.JUNCTION.PASS.voom.log.Res[target.circRNA,CaucasianID])
  tmpgenoytpe=as.character(target.SNPLabel[fenRA(as.numeric(target.genotype))+1])
  tmpGene=STA.paper[target.circRNA,"geneSymbol"]
  tmpEnsembl=STA.paper[target.circRNA,"gene"]
  #expGene=as.numeric(Tophat2.DGEList.voom.log[tmpEnsembl,CaucasianID])
  expGene=as.numeric(CMC.SVA.gene[tmpEnsembl,CaucasianID])
  tmp.df=data.frame(type=rep(c("gene","circRNA"),
                             each=length(expGene)),exp=c(expGene,expCirc),genotype=c(tmpgenoytpe,tmpgenoytpe))
  pdf(file =  "./pdf_C/internal_canonical_QTL_1.pdf",width = 8)
  violin.p <- ggplot(tmp.df,aes(x=genotype,y=exp)) +
    geom_violin(aes(fill=factor(genotype)),width=0.8) +
    geom_boxplot(fill="white",width=.2) + 
    theme_bw(base_size = 12)+
    theme(legend.position="none")+
    labs(title =paste(tmpGene," ",target.circRNA))+
    xlab(label="rs1550769")+
    ylab(label="Adjusted circRNA expression")+
    theme(legend.title = element_text(size = 0))+
    theme(axis.text.y = element_text(size = 24),axis.text.x = element_text(size = 24))+
    facet_wrap(~ type, scales="free")+labs(y = "Adjusted expression level (log2 CPM)")
  print(violin.p)
  dev.off()
}
{
  target.circRNA="chr14:88883055|88904247"
  target.genotype=as.numeric(read.table("/media/data3/CMC/data/genotype/Caucasian465/chr14_rs3179969_88862529_G_A.mat",header = F)[,-1])
  target.SNPLabel=c("1GG","2GA","3AA")
  expCirc=as.numeric(CIRI.JUNCTION.PASS.voom.log.Res[target.circRNA,CaucasianID])
  tmpgenoytpe=as.character(target.SNPLabel[fenRA(as.numeric(target.genotype))+1])
  tmpGene=STA.paper[target.circRNA,"geneSymbol"]
  tmpEnsembl=STA.paper[target.circRNA,"gene"]
  #expGene=as.numeric(Tophat2.DGEList.voom.log[tmpEnsembl,CaucasianID])
  expGene=as.numeric(CMC.SVA.gene[tmpEnsembl,CaucasianID])
  tmp.df=data.frame(type=rep(c("gene","circRNA"),
                             each=length(expGene)),exp=c(expGene,expCirc),genotype=c(tmpgenoytpe,tmpgenoytpe))
  
  pdf(file =  "./pdf_C/distant_canonical_QTL_1.pdf",width = 8)
  violin.p <- ggplot(tmp.df,aes(x=genotype,y=exp)) +
    geom_violin(aes(fill=factor(genotype)),width=0.8) +
    geom_boxplot(fill="white",width=.2) + 
    theme_bw(base_size = 12)+
    theme(legend.position="none")+
    labs(title =paste(tmpGene," ",target.circRNA))+
    xlab(label="rs3179969")+
    ylab(label="Adjusted circRNA expression")+
    theme(legend.title = element_text(size = 0))+
    theme(axis.text.y = element_text(size = 24),axis.text.x = element_text(size = 24))+
    facet_wrap(~ type, scales="free")+labs(y = "Adjusted expression level (log2 CPM)")
  print(violin.p)
  dev.off()
}

################revise Fig 1c##############
Depth.merge=data.frame(Gene=Tophat2.raw.counts.LIB_SIZE+CIRI.JUNCTION.depth,
                       Circ=CIRI.JUNCTION.depth,
                       Linear=CIRI.Linear.depth)

Depth.merge=Depth.merge[order(Depth.merge$Gene,decreasing = T),]

Depth.merge.long=data.frame(id=rep(1:nrow(Depth.merge),3),depth=c(Depth.merge$Gene,Depth.merge$Linear,Depth.merge$Circ),
                            type=rep(c("a","b","c"),each=nrow(Depth.merge)))
Depth.merge.long$col=c("#F8766D","#00BA38","#619CFF")[as.numeric(as.factor(Depth.merge.long$type))]
Depth.merge.long=subset(Depth.merge.long,type!="b")
{
  pdf("./pdf_C/sample_depth_distribution.pdf")
  plot(log(depth,10)~id,Depth.merge.long,type="h",col=Depth.merge.long$col,
       las=1,bty="l",xlab="samples",ylab="reads",xaxt="n")
  dev.off()
}
####################revised geneSymbol##########
library(biomaRt)
listMarts()
ensembl=useMart(biomart = "ENSEMBL_MART_ENSEMBL")
head(listDatasets(ensembl))
ensembleDB=useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembleFilter=listFilters(ensembleDB)
ensembl2gene=getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
      filters="ensembl_gene_id",values=unique(unlist(strsplit(STA.paper$gene,","))),mart=ensembleDB)
rownames(ensembl2gene)=ensembl2gene$ensembl_gene_id
ensembl2gene=ensembl2gene[unique(unlist(strsplit(STA.paper$gene,","))),]
rownames(ensembl2gene)=unique(unlist(strsplit(STA.paper$gene,",")))
ensembl2gene[which(ensembl2gene$hgnc_symbol==""),]="n.a."
ensembl2gene[which(is.na(ensembl2gene$hgnc_symbol)),]="n.a."

ensembl2gene["ENSG00000139915",]
STA.paper.new=STA.paper
for (i in 1:nrow(STA.paper.new)) {
  if(STA.paper.new[i,"region"]!="intergenic_region"){
    STA.paper.new[i,"geneSymbol"]=paste(ensembl2gene[unlist(strsplit(STA.paper.new[i,'gene'],",")),2],collapse = ",")
  }
}
write.table(STA.paper.new,"./STA.paper.new.txt",row.names = F,col.names = F,sep="\t",quote=F,na = "n.a.")
###################
# data: 1/31/2019
# author: Zelin Liu
# description: the second time revised paper
#################
library(pwr)

addFile3=read.delim('AnnotationCirc/AddFile3.txt',header=T)
DEwithAnno=cbind(SCZ.DGE[addFile3$circRNA.ID,],addFile3)
plot(-log(DEwithAnno$P.Value,base = 10)~DEwithAnno$Mean.expression..log2.circCPM.)
summary(lm(-log(DEwithAnno$P.Value,base = 10)~DEwithAnno$Mean.expression..log2.circCPM.))
subDEwithAnno=subset(DEwithAnno,Mean.circRNA.ratio>0.5)
tmp=p.adjust(subDEwithAnno$P.Value,'fdr')

hist(DEwithAnno$logFC)
hist(Gene.DEG$logFC)

#######CV############
#circRNA
matCircCount=CIRI.JUNCTION.PASS.SCZ.voom$E
meanCircCount=rowMeans(matCircCount)
sdCircCount=apply(matCircCount, 1, sd)
cvCircCount=sdCircCount/meanCircCount
hist(sdCircCount)
hist(cvCircCount,xlim = c(0,1))

matCircCount.SCZ=matCircCount[,which(covariates.SCZ$Dx=='SCZ')]
meanCircCount.SCZ=rowMeans(matCircCount.SCZ)
sdCircCount.SCZ=apply(matCircCount.SCZ, 1, sd)
cvCircCount.SCZ=sdCircCount.SCZ/meanCircCount.SCZ

matCircCount.CTL=matCircCount[,which(covariates.SCZ$Dx=='Control')]
meanCircCount.CTL=rowMeans(matCircCount.CTL)
sdCircCount.CTL=apply(matCircCount.CTL, 1, sd)
cvCircCount.CTL=sdCircCount.CTL/meanCircCount.CTL



#linearRNA
matLinearCount=Tophat2.DGEList.SCZ.voom$E
meanLinearCount=rowMeans(matLinearCount)
sdLinearCount=apply(matLinearCount, 1, sd)
cvLinearCount=sdLinearCount/meanLinearCount

idx=which(meanLinearCount>=min(meanCircCount) & meanLinearCount<=max(meanCircCount))
hist(sdLinearCount[idx])
hist(cvLinearCount)
hist(cvLinearCount[which(cvLinearCount>0&cvLinearCount<1)],xlim = c(0,1))
hist(cvLinearCount[rownames(Gene.DEG)[which(Gene.DEG$adj.P.Val<0.05)]])
hist(cvLinearCount[which(cvLinearCount>0&cvLinearCount<1)][rownames(Gene.DEG)[which(Gene.DEG$adj.P.Val<0.05)]],xlim = c(0,1))
plot(meanLinearCount[idx]~cvLinearCount[idx])

matLinearCount.SCZ=matLinearCount[,which(covariates.SCZ$Dx=='SCZ')]
meanLinearCount.SCZ=rowMeans(matLinearCount.SCZ)
sdLinearCount.SCZ=apply(matLinearCount.SCZ, 1, sd)
cvLinearCount.SCZ=sdLinearCount.SCZ/meanLinearCount.SCZ

matLinearCount.CTL=matLinearCount[,which(covariates.SCZ$Dx=='Control')]
meanLinearCount.CTL=rowMeans(matLinearCount.CTL)
sdLinearCount.CTL=apply(matLinearCount.CTL, 1, sd)
cvLinearCount.CTL=sdLinearCount.CTL/meanLinearCount.CTL


############evaluate the sample size needed to detect DE########
effectSize=DEwithAnno$logFC/sdCircCount
range(effectSize)
hist(effectSize)
pwr.t.test(d=max(abs(effectSize)),sig.level=0.05/10559,power=0.8,type="two.sample",alternative="two.sided")$n


effectLinear=Gene.DEG$logFC/sdLinearCount[rownames(Gene.DEG)]
range(effectLinear)
pwr.t.test(d=max(abs(effectLinear)),sig.level=0.05/length(effectLinear),power=0.8,type="two.sample",alternative="two.sided")
{
  pdf('./pdf_C/effect_size_hist.gene.pdf')
  par(mfrow=c(2,1),mar=c(2,5,1,1))
  hist(effectLinear,col=rgb(5,5,250,alpha = 180,maxColorValue = 255),las=1,
       ylim = c(0,2500),xlim=c(-0.5,0.5),main='',xlab='')
  par(new=T)
  hist(effectLinear[rownames(Gene.DEG)[which(Gene.DEG$adj.P.Val<0.05)]],
       col=rgb(5,250,5,alpha = 180,maxColorValue = 255),ylim = c(0,2500),ylab = '',xaxt='n',
       breaks = 20,xlim=c(-0.5,0.5),main='',xlab='',border=NA,yaxt='n')
  dev.off()
}
#hist(effectSize,col=rgb(230,13,13,alpha = 180,maxColorValue = 255),ylim = c(0,2500),xlim=c(-0.5,0.5),main='',xlab='')
{
  pdf('./pdf_C/effect_size_density.pdf')
  plot(density(effectLinear),col=rgb(5,5,250,alpha = 180,maxColorValue = 255),lwd=2,ylim=c(0,4),main='',las=1)
  lines(density(effectSize),col=rgb(230,13,13,alpha = 180,maxColorValue = 255),lwd=2)
  dev.off()
}

############

tmp=SCZ.DGE[names(cvCircCount)[which(cvCircCount<0.2)],]
tmp$newFDR=p.adjust(tmp$P.Value,'fdr')

#####group-variance####
matCircScale=abs(t(scale(t(matCircCount))))
matCircScale=(t(scale(t(matCircCount),scale = F)))^2

varP=vector(length = nrow(matCircScale))
idSCZ=which(covariates.SCZ$Dx=='SCZ')
idCTL=which(covariates.SCZ$Dx=='Control')
for (i in 1:nrow(matCircScale)) {
  varP[i]=t.test(matCircScale[i,idSCZ],matCircScale[i,idCTL])$p.value
}
varFDR=p.adjust(varP,'fdr')
