#!/usr/bin/Rscript
options(warn=-1)
#########################################################################################################################################
#
#                  This R script allows to replicate the trait analysis proposed in Fourneir-Level el al. (2018, submitted)
#
# - The input data can be sourced at the following URL:  
# - The input data folder should be located in the R working directory without altering the folder names or structure
#
#########################################################################################################################################

#########################################################################################################################################
#
#                                     Distribution of serquencing coverage by windows of 500bp
# - the parameter q can be tuned to exclude the extreme range of coverage.
# - by default interquarile windows of coverage for each pool were kept
#
#########################################################################################################################################

rm(Data)
setwd("Sequencing_Summary")

q=0.25
Graph=F

CovFiles=list.files()
CovFiles=CovFiles[grep(".coverage",CovFiles)]

for (i in CovFiles) {
	Data.sub=read.table(i,sep="\t",header=T)
	if (!exists("Data")) {
		Data=Data.sub
	} else {Data=cbind(Data,Data.sub[,5])}
}

Data=Data[,-c(3,4)]
colnames(Data)=c("Chr","Position",gsub(".coverage","",CovFiles))

CHR.std=c("2L","2R","3L","3R","4","X")
Data.std=matrix(NA,ncol=ncol(Data),nrow=0)
for (x in CHR.std) Data.std=rbind(Data.std,Data[Data$Chr==x,])

i=1
Trans.std=vector()
while (Data.std$Chr[i] %in% CHR.std) {
	i=i+1
	if (i==nrow(Data.std)) break
	if (Data.std$Chr[i]!=Data.std$Chr[i+1]) Trans.std=c(Trans.std,i)
}

CovFrac=function(k,X) {sum(Data.std[,k]>rep(X,nrow(Data.std)))/(2*nrow(Data.std))}

CovMat=matrix(NA,ncol=ncol(Data)-2,nrow=6)
colnames(CovMat)=names(Data.std)[-c(1,2)]
rownames(CovMat)=c(5,10,15,20,30,40)
for (k in names(Data.std)[-c(1,2)]) {
	for (X in c(5,10,15,20,30,40)) {
		CovMat[paste(X),k]=CovFrac(k,X)
	}
}

if (Graph==T) {
	dev.new()
	par(mfrow=c(3,5))
	for (i in 3:ncol(Data.std)) {
		hist(log1p(Data.std[,i]),main=names(Data.std)[i])
		abline(v=quantile(log1p(Data.std[,i]),.05),lty=2,col=2)
		abline(v=quantile(log1p(Data.std[,i]),.25),lty=1,col=2)
		abline(v=quantile(log1p(Data.std[,i]),.75),lty=1,col=2)
		abline(v=quantile(log1p(Data.std[,i]),.95),lty=2,col=2)
		if ((i-2) %% 15 == 0) {
			print(i)
			dev.new()
			par(mfrow=c(3,5))
		}
	}
}

MedCov_min=apply(Data.std[,-c(1,2)]/2,2,function(x) {quantile(x,q)})
MedCov_max=apply(Data.std[,-c(1,2)]/2,2,function(x) {quantile(x,1-q)})

Lib.names=matrix(unlist(strsplit(names(MedCov_min),"-")),ncol=3,byrow=T)
Lib.names=unique(paste(Lib.names[,1],Lib.names[,2],sep="-"))

setwd("../")

cat("LibMedCov={\n",file=paste("LibMedCov_q",q*100,".txt",sep=""))
for (x in Lib.names) {
	cat(paste("\'",x,"_min\':[",paste(MedCov_min[grep(x,names(MedCov_min))],collapse=","),"],\n",sep=""),file=paste("LibMedCov_q",q*100,".txt",sep=""),append=T)
	cat(paste("\'",x,"_max\':[",paste(MedCov_max[grep(x,names(MedCov_max))],collapse=","),"],\n",sep=""),file=paste("LibMedCov_q",q*100,".txt",sep=""),append=T)
}
cat("}\n",file=paste("LibMedCov_q",q*100,".txt",sep=""),append=T)



