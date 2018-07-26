#!/usr/bin/Rscript
options(warn=-1)
#########################################################################################################################################
#
#                  This R script allows to replicate the trait analysis proposed in Fourneir-Level el al. (2018, submitted)
#
# - It requires the vegan, fields and raster packages to be installed to perform the mixed-linear modelling
# - The input data can be sourced at the following URL:  
# - The input data folder should be located in the R working directory without altering the folder names or structure
#
#########################################################################################################################################

##################################################################################################
#                                        BioGeo/BioClim inferences
##################################################################################################
library(vegan)
library(raster)

Pheno.name="STX-20C_All"
TAJ=read.table(Sys.glob(paste("PopGen/",unlist(strsplit(Pheno.name,"-"))[1],".npstats",sep="")),header=T)
FST=read.table(Sys.glob("PopGen/FST.txt"),header=F) #FST were computed from the mpileup files using the Popoolation package
FST[,2]=FST[,2]+FST[1,2]
FST=FST[paste(FST[,1],FST[,2]) %in% paste(TAJ[,1],TAJ[,2]),]
FST.colnames=vector(length=ncol(FST))
FST.colnames[1:5]=c("CHR","Window","nSNP","percCOV","minCOV")
for (i in 6:ncol(FST)) {
	POPs=strsplit(as.character(FST[1,i]),"=")[[1]][1]
	FST.colnames[i]=POPs
	FST[,i]=as.numeric(gsub(paste(POPs,"=",sep=""),"",FST[,i]))
}
colnames(FST)=FST.colnames

PopFX=read.csv(Sys.glob("PopGen/PopFX.csv"),header=T) #Population effects are the same as the ones computed in model 3 of the PhenotypeAnalysis.r script
rownames(PopFX)=PopFX[,1]

POP=unique(unlist(strsplit(names(FST)[-(1:5)],":")))

FST.mean=apply(FST[,6:ncol(FST)],2,mean)
FST.dist=IM.dist=TMP.dist=matrix(NA,length(POP),length(POP))
colnames(FST.dist)=colnames(IM.dist)=colnames(TMP.dist)=POP
rownames(FST.dist)=rownames(IM.dist)=rownames(TMP.dist)=POP
for (i in 1:length(POP)) for (j in 1:length(POP)) if (j>i) {
	FST.dist[i,j]=as.numeric(FST.mean[names(FST.mean)==paste(POP[i],POP[j],sep=":")])
	IM.dist[i,j]=sqrt((PopFX$IM_CTIM[rownames(PopFX)==POP[i]]-PopFX$IM_CTIM[rownames(PopFX)==POP[j]])^2)
	TMP.dist[i,j]=sqrt((PopFX$Temp30C[rownames(PopFX)==POP[i]]-PopFX$Temp30C[rownames(PopFX)==POP[j]])^2)
}

CLIM.dist=matrix(c(NA,0.5,0.5,0.5,0.5,0,0,1,0.5,0.5,0.5,1,1,0,1,0.5,
					NA,NA,1,1,1,0.5,0.5,0.5,0,0,1,0.5,0.5,0.5,0.5,0,
					NA,NA,NA,0,0,0.5,0.5,0.5,1,1,0,0.5,0.5,0.5,0.5,1,
					NA,NA,NA,NA,0,0.5,0.5,0.5,1,1,0,0.5,0.5,1,0.5,0.5,
					NA,NA,NA,NA,NA,0.5,0.5,0.5,1,1,0,0.5,0.5,1,0.5,0.5,
					NA,NA,NA,NA,NA,NA,0,1,0.5,0.5,0.5,1,1,0,1,0.5,
					NA,NA,NA,NA,NA,NA,NA,1,0.5,0.5,0.5,1,1,0,1,0.5,
					NA,NA,NA,NA,NA,NA,NA,NA,0.5,0.5,0.5,0,0,1,0,0.5,
					NA,NA,NA,NA,NA,NA,NA,NA,NA,0,1,0.5,0.5,0.5,0.5,0,
					NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1,0.5,0.5,0.5,0.5,0,
					NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.5,0.5,0.5,0.5,1,
					NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0,1,0,0.5,
					NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1,0,0.5,
					NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1,0.5,
					NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.5,
					NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),length(POP),length(POP),byrow=T)
colnames(CLIM.dist)=POP
rownames(CLIM.dist)=POP

Dist.OrthoDrom<-function(Lat1,long1,Lat2,long2) {
Lat1=Lat1*(3.1416/180)
long1=long1*(3.1416/180)
Lat2=Lat2*(3.1416/180)
long2=long2*(3.1416/180)
result<-6371*2*asin(sqrt(sin((Lat2-Lat1)/2)^2+cos(Lat2)*cos(Lat1)*sin((long2-long1)/2)^2))
return(result)
}

LOC=matrix(c("SLA",	30.43,	-91.12,	0,
           	"NLA",	32.5,	-92.14,0,
           	"NTX",	33.57,	-100.9,0,
           	"STX",	30.23,	-99.66,0,
           	"NCA",	39.15,	-123.65,0,
           	"SCA",	39,		-123.21,0,
           	"NB",	45.25,	-66.07,0,
           	"NS",	45.08,	-64.37,0,
           	"NQ1",	-17.54,	145,1,
           	"NQ2",	-17.50,	145,1,
           	"SQ1",	-25.17,	149.93,1,
           	"SQ2",	-23.13,	150.73,1,
           	"SA",	-37.82,	140.77,1,
           	"VIC",	-37.43,	142.02,1,
           	"NTS",	-41.23,	146.98,1,
           	"STS",	-43,	146.9,1
           	),ncol=4,byrow=T)

xy=data.frame(as.numeric(LOC[,3]),as.numeric(LOC[,2]))

GEO.dist=matrix(NA,length(POP),length(POP))
colnames(GEO.dist)=POP
rownames(GEO.dist)=POP
for (i in 1:length(POP)) for (j in 1:length(POP)) if (j>i) GEO.dist[i,j]=Dist.OrthoDrom(Lat1=as.numeric(LOC[LOC[,1]==POP[i],2]),
																						long1=as.numeric(LOC[LOC[,1]==POP[i],3]),
																						Lat2=as.numeric(LOC[LOC[,1]==POP[j],2]),
																						long2=as.numeric(LOC[LOC[,1]==POP[j],3]))

LAT.dist=matrix(NA,length(POP),length(POP))
colnames(LAT.dist)=POP
rownames(LAT.dist)=POP
for (i in 1:length(POP)) for (j in 1:length(POP)) if (j>i) LAT.dist[i,j]=abs(as.numeric(LOC[LOC[,1]==POP[i],2]))-abs(as.numeric(LOC[LOC[,1]==POP[j],2]))

tmax <- raster(Sys.glob("GIS/Bioclim/bio10.bil"))
pmin <- raster(Sys.glob("GIS/Bioclim/bio17.bil"))
tmax_LOC=as.matrix(extract(tmax,xy),ncol=1)
pmin_LOC=as.matrix(extract(pmin,xy),ncol=1)
rownames(tmax_LOC) = rownames(pmin_LOC) = LOC[,1]
tmax_LOC=apply(tmax_LOC,1,function(x) sqrt((x-tmax_LOC)^2))
pmin_LOC=apply(pmin_LOC,1,function(x) sqrt((x-pmin_LOC)^2))
rownames(tmax_LOC) = rownames(pmin_LOC) = LOC[,1]
tmax_LOC=tmax_LOC[order(rownames(tmax_LOC)),]
tmax_LOC=tmax_LOC[,order(colnames(tmax_LOC))]
pmin_LOC=pmin_LOC[order(rownames(pmin_LOC)),]
pmin_LOC=pmin_LOC[,order(colnames(pmin_LOC))]

print("Mantel statistic for the correlation of pairwise Fst between populations and geographic distance")
mantel(as.dist(t(FST.dist)),as.dist(t(GEO.dist)))$statistic
print("Partial Mantel statistic for the correlation of pairwise Fst between populations and climate controlling for geographic distance")
mantel.partial(as.dist(t(FST.dist)),as.dist(t(CLIM.dist)),as.dist(t(GEO.dist)))$statistic
print("Mantel statistic for the correlation of pairwise Fst between populations and temperature difference controlling for geographic distance")
mantel.partial(as.dist(t(FST.dist)),as.dist(tmax_LOC),as.dist(t(GEO.dist)))$statistic
print("Mantel statistic for the correlation of pairwise Fst between populations and precipitation difference controlling for geographic distance")
mantel.partial(as.dist(t(FST.dist)),as.dist(pmin_LOC),as.dist(t(GEO.dist)))$statistic

if (file.exists("PopGen/Mantel.RData")) {
	load(Sys.glob("PopGen/Mantel.RData"))
} else {
Mantel=matrix(NA,nrow=0,ncol=4)
colnames(Mantel)=c("MantGeo","pMantClim","pMantmaxtemp","pMantminprec")
for (z in 1:nrow(FST)) {
	#At the whole genome level
	#FST.mean=apply(FST[,6:ncol(FST)],2,mean)
	#or for each window
	FST.mean=FST[z,6:ncol(FST)]
	FST.dist=matrix(NA,length(POP),length(POP))
	colnames(FST.dist)=POP
	rownames(FST.dist)=POP
	for (i in 1:length(POP)) for (j in 1:length(POP)) if (j>i) FST.dist[i,j]=as.numeric(FST.mean[names(FST.mean)==paste(POP[i],POP[j],sep=":")])
	test=c(mantel.partial(as.dist(t(FST.dist)),as.dist(t(GEO.dist)),as.dist(t(CLIM.dist)),permutations=0)$statistic,
           mantel.partial(as.dist(t(FST.dist)),as.dist(t(CLIM.dist)),as.dist(t(GEO.dist)),permutations=0)$statistic,
           mantel.partial(as.dist(t(FST.dist)),as.dist(t(CLIM.dist)),as.dist(t(GEO.dist)),permutations=0)$statistic,
           mantel.partial(as.dist(t(FST.dist)),as.dist(tmax_LOC),as.dist(t(GEO.dist)),permutations=0)$statistic,
           mantel.partial(as.dist(t(FST.dist)),as.dist(pmin_LOC),as.dist(t(GEO.dist)),permutations=0)$statistic)
	Mantel=rbind(Mantel,test)
	}
	save(Mantel,file="PopGen/Mantel.RData")
}
FST=cbind(FST,Mantel)

#dev.new()
#plot('n',xlim=c(0,nrow(FST)),ylim=c(0,.07))
#for (i in 6:(ncol(FST)-2)) { #lines(FST[,i],col=i)
#	lines(predict(smooth.spline(FST[,i]^2,spar=0.05))$y,col=i)
#}
#Chr.break=which(c("2L",FST[,1])!=c(FST[,1],"X"))
#abline(v=Chr.break,lty=2)
#lines(FST[,i],col=i)

Mantel[is.na(Mantel)]=0
#dev.new()
Clim.mantel=predict(smooth.spline(Mantel[,2],spar=0.05))$y
Geo.mantel=predict(smooth.spline(Mantel[,1],spar=0.05))$y
FST.avg=predict(smooth.spline(apply(FST[,grep(":",names(FST))],1,mean),spar=0.05))$y
png("Fig4b_MantelClimvsGeo.png",width=1200,height=300)
par(mar=c(3.5,9,2.5,5),mgp=c(1.8,0.6,0),xaxs="i")
plot(Clim.mantel,ylim=c(min(c(Clim.mantel,Geo.mantel)),max(c(Clim.mantel,Geo.mantel))),type="line",lwd=2,lty=3,col=grey(.5),xaxt='n',ylab="",xlab="Genomic Position")
mtext(text="Geographic and Climate Differentiation", side = 3, cex=1.3, font=2, adj=0)
lines(Geo.mantel,col=grey(.7),lwd=1.5)
lines(FST.avg,lwd=1.5)
X=c(1,which(c(FST[1,1],FST[,1]) != c(FST[,1],FST[nrow(FST),1])))
BK=which(c(FST[1,1],FST[,1]) != c(FST[,1],FST[nrow(FST),1]))-1
CHR=FST[c(1,which(c(FST[1,1],FST[,1]) != c(FST[,1],FST[nrow(FST),1]))),1]
axis(side=1,at=X+(c(BK,nrow(FST))-c(1,BK))/2,par("usr")[3]+.1*abs(par("usr")[3]),labels=CHR,tick=F,las=1)
abline(v=which(c(FST[1,1],FST[,1]) != c(FST[,1],FST[nrow(FST),1])),lty=2, col="blue")
legend(.82*par("usr")[2],.9*par("usr")[4],legend=c("r(Climate Distance)","r(Geographic Distance)","avg. Fst"),lwd=1.5,col=c(grey(.5),grey(.8),1),lty=c(3,1,1),bty='n')
par(xpd=T)
mtext(text="B", side = 3, cex=1.5, font=2, at=-.1*par("usr")[2])
dev.off()

##################################################################################################
#                                        Testing reduced Fst at Candidate loci
##################################################################################################

PermList=list.files(path="PopGen",pattern="FstMantel_permutation")
PermFull=list()
avg_FST=c()
PermFST=matrix(NA,nrow=0,ncol=2)
for (f in PermList) {
	load(paste("PopGen/",f,sep=""))
	PermFull=c(PermFull,PermAnal)
	PermFST=rbind(PermFST,cbind(names(PermAnal),PermAnal[[1]][["avg_FST_perm"]]))
	avg_FST=c(avg_FST,PermAnal[[1]][["avg_FST_obs"]])
}

#dev.new()
png("Fig4d_FstAmongCandidates.png",width=400,height=300)
par(mar=c(4,3,3,1),mgp=c(2,.7,0))
boxplot(as.numeric(PermFST[,2])~factor(PermFST[,1],levels=1:length(PermList),labels=as.character(length(PermList):1)),border=grey(.5),ylab=expression('average F'[st]),xlab="Number of GWAlpha tests with gene overlap",main="Genetic differention for resistance loci")
lines(avg_FST[order(as.numeric(names(PermFull)),decreasing=F)],lwd=1.5,col=2)
par(xpd=T)
mtext(text="D", side = 3, cex=1.5, font=2, at=-.05*par("usr")[2],line=1)
dev.off()

##################################################################################################
#
#                       Analysis of the level of resistance for Ancestral alleles
#
##################################################################################################

PHENO=c("NLA-20C_All","NLA-30C_All","SLA-20C_All","SLA-30C_All","NTX-20C_All","STX-20C_All",
        "NTX-30C_All","STX-30C_All","NB-20C_All","NB-30C_All","NS-20C_All","NS-30C_All",
        "SCA-20C_All","SCA-30C_All","NCA-20C_All","NCA-30C_All","NQ1-30C_All","NQ2-30C_All",
        "SQ1-30C_All","SQ2-30C_All","NTS-30C_All","STS-30C_All","SA-30C_All","VIC-30C_All",
        "NQ1-20C_All","NQ2-20C_All","SQ1-20C_All","SQ2-20C_All","NTS-20C_All","STS-20C_All",
        "SA-20C_All","VIC-20C_All")

Pheno.name="NQ1-30C_All"

for (Pheno.name in PHENO) {

print(Pheno.name)
GWAlpha=read.csv(paste("GWAS/GWAlpha_",Pheno.name,"_out.csv",sep=""),header=T)
GWAlpha=GWAlpha[GWAlpha[,1]!=4,] # Population genetics analysis excluded the chromosome 4
GWAlpha=GWAlpha[GWAlpha[,3]!="INV"&GWAlpha[,3]!="TE",] # Population genetics analysis excluded the inversions and the transposable elements
TAJ=read.table(Sys.glob(paste("PopGen/",unlist(strsplit(Pheno.name,"-"))[1],".npstats",sep="")),header=T)
ANC=read.table(paste("PopGen/GWAlpha_",Pheno.name,"_out_anc.bed",sep=""),header=F)

#CHR.labs=c("9"="4","2"="2LHet","4"="2RHet","6"="3LHet","8"="3RHet","11"="XHet","1"="2L","3"="2R","5"="3L","7"="3R","10"="X","12"="YHet")
if ("2L" %in% GWAlpha[,1]) {
CHR.labs=c("5"="4","1"="2L","2"="2R","3"="3L","4"="3R","6"="X")
GWAlpha=GWAlpha[order(as.character(GWAlpha[,1]),as.numeric(as.character(GWAlpha[,2]))),]
} else {
CHR.labs=unique(as.character(GWAlpha[,1]))
names(CHR.labs)=unique(as.character(GWAlpha[,1]))
GWAlpha=GWAlpha[order(as.character(GWAlpha[,1]),as.numeric(as.character(GWAlpha[,2]))),]
}

OUT=GWAlpha[GWAlpha[,1] %in% CHR.labs,]

#Compute the Threshold
GWAlpha=abs(as.numeric(as.character(OUT$Alpha)))
GWAlpha.order=GWAlpha[order(GWAlpha)]
GWAlpha.QQ=predict(smooth.spline(GWAlpha.order,spar=.6),deriv=2)$y
Thresh=GWAlpha.order[which(GWAlpha.QQ==max(GWAlpha.QQ))]

GWAlpha=as.numeric(as.character(OUT$Alpha))
LL=function(x,y) sum(dnorm(y,x[1],x[2],log=T))
EST=optim(par=c(mu="0",sd=sd(GWAlpha)),fn=LL,y=GWAlpha,control=list(fnscale=-1,reltol=1e-10))$par
pval=ifelse(pnorm(GWAlpha,EST[1],EST[2]) > 0.5, 2*(1-pnorm(GWAlpha,EST[1],EST[2])), 2*pnorm(GWAlpha,EST[1],EST[2]))
p_score=-log10(pval)
Padj=p.adjust(pval,method="BH")
q=0.05
if (sum(Padj<q)!=0) {
	pThresh=max(pval[Padj<q])
} else pThresh=min(pval)/2

#Effect of Ancestrality of GWAlpha loci
ANC=ANC[paste(ANC[,1],ANC[,2]) %in% paste(OUT[,1],OUT[,2]),]
OUT=OUT[paste(OUT[,1],OUT[,2]) %in% paste(ANC[,1],ANC[,2]),]
if (Pheno.name=="SLA-30C_All") OUT=OUT[-62400,]
sum(paste(ANC[,1],ANC[,2])==paste(OUT[,1],OUT[,2]))==dim(OUT)[1]

Stater=function(x) {
	x=table(x[x!='N'])
	x=x[order(x,decreasing=T)]
	if (length(x)>1) {
		if (x[1]==x[2]) return(paste(names(x)[1],names(x)[2],sep="/"))
	} else return(names(x)[1])
}

Anc=apply(ANC[,-c(1:2)],1,Stater)
Anc[sapply(Anc, is.null)]=NA
OUT$Anc=toupper(unlist(Anc))

OUT$is.anc=apply(OUT,1,function(x) grepl(x[3],x[7]))
ANC=OUT[!is.na(OUT$Anc),]
Anc.table=cbind(table(ANC$is.anc[as.numeric(as.character(ANC$Alpha))<(-Thresh)]),table(ANC$is.anc[as.numeric(as.character(ANC$Alpha))>Thresh]))
print(paste("Chi-sq for Ancestral enrichment:",chisq.test(Anc.table)$p.value))
Anc.ratio=Anc.table[2,]/Anc.table[1,]
if (Anc.ratio[1]>Anc.ratio[2]) {
	print("with ancestral resistance ")
} else print("with derived resistance")
}

##################################################################################################
#                     Selective sweeps from Pool-HMM
##################################################################################################

#compiling the data + graphs
#for (k in c("0.001","0.0001","0.00001","0.000001","0.0000001")) {
k="0.000001"
ChrPos=read.table(paste("PopGen/",list.files(path="PopGen",pattern="npstats")[1],sep=""),header=T)[,1:2]
Sweep.list=paste("PopGen/k_",k,"/",list.files(path=paste("PopGen/k_",k,sep=""),pattern="\\.stat"),sep="")

Pop=c(North_Queensland1="NQ1",North_Queensland2="NQ2",
      South_Queensland1="SQ1",South_Queensland2="SQ2",
      South_Australia="SA",Victoria="VIC",
      North_Tasmania="NTS",South_Tasmania="STS",
      North_Lousiana="NLA",South_Lousiana="SLA",
      North_Texas="NTX",South_Texas="STX",
      North_California1="NCA",North_California2="SCA",
      New_Brunswick="NB",Nova_Scotia="NS")

Pop=Pop[c(5,13,8,7,16,15,4,3,11,2,6,14,12,1,10,9)]

Chr=c("2L","2R","3L","3R","X")

Sweep.matrix=matrix(0,nrow=nrow(ChrPos),ncol=length(Pop))
colnames(Sweep.matrix)=Pop
for (p in Pop) {
	Sweep.file=Sweep.list[grep(p,Sweep.list)]
	Sweep.block=matrix(NA,nrow=0,ncol=3)
	for (x in Chr) {
		Sweep.chr=Sweep.file[grep(paste("down_",x,sep=""),Sweep.file)]
		if (length(Sweep.chr)==0) next
		test=read.table(Sweep.chr,header=F,nrows=1)
		if (test[1,1]==0) next
		Sweep=read.table(Sweep.chr,header=F,skip=1)
		for (i in 1:nrow(Sweep)) Sweep.block=rbind(Sweep.block,cbind(x,seq(round(Sweep$V1[i] / ChrPos[1,2])*ChrPos[1,2],round(Sweep$V2[i]/ChrPos[1,2])*ChrPos[1,2],by=ChrPos[1,2]),Sweep$V3[i]))
		}
		dbl=names(table(paste(Sweep.block[,1],as.integer(Sweep.block[,2]))))[table(paste(Sweep.block[,1],as.integer(Sweep.block[,2])))>1]
		idx.rm=c()
		if (length(dbl)>0) {
			for (i in dbl) {
				idx=which(paste(Sweep.block[,1],as.integer(Sweep.block[,2]))==i)
				Sweep.block[paste(Sweep.block[,1],as.integer(Sweep.block[,2]))==i,3]=mean(as.numeric(Sweep.block[paste(Sweep.block[,1],as.integer(Sweep.block[,2]))==i,3]))
				idx.rm=c(idx.rm,idx[-1])
			}
			Sweep.block=Sweep.block[-idx.rm,]
		}
	Sweep.matrix[paste(ChrPos[,1],ChrPos[,2]) %in% paste(Sweep.block[,1],as.integer(Sweep.block[,2])),p]=as.numeric(Sweep.block[,3])
}

Sweep=cbind(ChrPos,Sweep.matrix)
save(Pop,Sweep,file=paste("PopGen/Hmm_SelSweep_k",paste(k),".RData",sep=""))
#}

#Getting the main candidate sweeps+ producing the Figure

load(paste("PopGen/Hmm_SelSweep_k",paste(k),".RData",sep=""))
load("GWAS/Candidate_Genes.RData")
Dm3Genes=read.table("Ref_Genome/Dmel_v5.57.genes", header=F)
Pop=Pop[]
Sweep.matrix=Sweep[,-c(1:2)]
ChrPos=Sweep[,c(1:2)]
Int.size=2500
Sweep_Thresh=0.99
Sweep.mean=apply(Sweep.matrix,1,mean)
Loci.list=ChrPos[Sweep.mean>quantile(Sweep.mean,Sweep_Thresh),]
Loci.list=Loci.list[order(Sweep.mean[Sweep.mean>quantile(Sweep.mean,Sweep_Thresh)],decreasing=T),]
Loci.list[,2]=Loci.list[,2]-Int.size
Gene.list=rep("Intergenic",nrow(Loci.list))
for (i in 1:nrow(Loci.list)) {
	gene=Dm3Genes[Dm3Genes[,1] == as.character(Loci.list[i,1]) &
	              Dm3Genes[,4]-Int.size < as.numeric(as.character(Loci.list[i,2])) & 
	              Dm3Genes[,5]+Int.size > as.numeric(as.character(Loci.list[i,2])), 10]
	if (length(gene)>0) {
	if (length(gene)>1) Gene.list[i]=paste(gene,collapse=", ")
	else Gene.list[i]=as.character(gene)
	}
}

Sweep_table=cbind(Loci.list,Gene.list,Sweep.mean[Sweep.mean>quantile(Sweep.mean,Sweep_Thresh)][order(Sweep.mean[Sweep.mean>quantile(Sweep.mean,Sweep_Thresh)],decreasing=T)])
Sweep_genes=unlist(strsplit(gsub(" ","",Sweep_table[,3]),","))
sum(unique(unlist(All_Genes))[table(unlist(All_Genes))>1] %in% unique(Sweep_genes))/sum(table(unlist(All_Genes))>1)
unique(unlist(All_Genes))[table(unlist(All_Genes))>1][unique(unlist(All_Genes))[table(unlist(All_Genes))>1] %in% unique(Sweep_genes)]
#sum(unique(Sweep_genes) %in% unique(unlist(All_Genes))[table(unlist(All_Genes))>1])/length(unique(Sweep_genes))
#unique(Sweep_genes)[unique(Sweep_genes) %in% unique(unlist(All_Genes))[table(unlist(All_Genes))>1]]
unique(Sweep_table[,1])
min(Sweep_table[Sweep_table[,1]=="2R",2])
max(Sweep_table[Sweep_table[,1]=="2R",2])

#dev.new()
png("Fig4c_pHMM.png",width=1200,height=400)
par(mar=c(3,9,2,5),xaxs="i")
COL=colorRampPalette(c("#ffffff","#f7b7ca","#fc6a95","#ed0449","#7a0125","#020000"))(100)
image(as.matrix(100*Sweep.matrix/max(Sweep.matrix)),col=COL,xaxt='n',yaxt='n')
mtext(text="Selective Sweeps", side = 3, cex=1.3, font=2, adj=0)
X=c(1,which(c(ChrPos[1,1],ChrPos[,1]) != c(ChrPos[,1],ChrPos[nrow(ChrPos),1])))
BK=which(c(ChrPos[1,1],ChrPos[,1]) != c(ChrPos[,1],ChrPos[nrow(ChrPos),1]))-1
CHR=ChrPos[c(1,which(c(ChrPos[1,1],ChrPos[,1]) != c(ChrPos[,1],ChrPos[nrow(ChrPos),1]))),1]
axis(side=1,at=(X+(c(BK,nrow(ChrPos))-c(1,BK))/2)/nrow(ChrPos),(par("usr")[3]+.1*abs(par("usr")[3]))/nrow(ChrPos),labels=CHR,tick=F,las=1)
LABS=gsub("2"," 2",gsub("1"," 1",gsub("_"," ",names(Pop))))
axis(side=2,at=seq(0,1,1/(length(Pop)-1)),LABS,las=2)
rect(0,mean(seq(0,1,1/(length(Pop)-1))[10:11]),1,par("usr")[4],border="#41d894",lwd=1.3)
legend("topright","Populations with derived susceptibility",text.col="#41d894",text.font=3,bty='n')
abline(v=(which(c(ChrPos[1,1],ChrPos[,1]) != c(ChrPos[,1],ChrPos[nrow(ChrPos),1])))/nrow(ChrPos),lty=2, col="blue")
library(fields)
image.plot(as.matrix(Sweep.matrix), col=COL,legend.only=T,legend.lab=expression(-log[10](1-italic(p[Selection]))))
par(xpd=T)
mtext(text="C", side = 3, cex=1.5, font=2, at=-.1*par("usr")[2])
dev.off()

#plotting the outcome of the MonteCarlo analysis that can be produced using the pHMM_MonteCarloAnalysis.r
PermList=list.files(path="PopGen/",pattern="SweepAll")
PermFull=list()
avg_pHMM=c()
PermpHMM=matrix(NA,nrow=0,ncol=2)
for (f in PermList) {
	load(paste("PopGen/",f,sep=""))
	PermFull=c(PermFull,PermAnal)
	PermpHMM=rbind(PermpHMM,cbind(names(PermAnal),PermAnal[[1]][["avg_pHMM_perm"]]))
	avg_pHMM=c(avg_pHMM,PermAnal[[1]][["avg_pHMM_obs"]])
}

png("Fig4e_pHMMAmongCandidates.png",width=400,height=300)
par(mar=c(4,3,3,1),mgp=c(2,.7,0))
boxplot(as.numeric(PermpHMM[,2])~factor(PermpHMM[,1],levels=1:length(PermFull),labels=as.character(length(PermFull):1)),border=grey(.5),ylim=c(0,1),ylab='average pHMM',xlab="Number of GWAlpha tests intersected",main="Selective sweeps")#paste("All populations (",unlist(strsplit(PermList[1],"_"))[2],")",sep=""))
mtext("All populations")
lines(avg_pHMM[order(as.numeric(names(PermFull)),decreasing=T)],lwd=1.5,col=2)
par(xpd=T)
mtext(text="E", side = 3, cex=1.5, font=2, at=-.05*par("usr")[2],line=1)
dev.off()

#plotting the outcome of the MonteCarlo analysis that can be produced using the pHMM_Ancestral_MonteCarloAnalysis.r
PermList=list.files(path="PopGen/",pattern="SweepAnc")
PermFull=list()
avg_pAncestral=c()
PermpAncestral=matrix(NA,nrow=0,ncol=2)
for (f in PermList) {
	load(paste("PopGen/",f,sep=""))
	PermFull=c(PermFull,PermAnal)
	PermpAncestral=rbind(PermpAncestral,cbind(names(PermAnal),PermAnal[[1]][["avg_pAncestral_perm"]]))
	avg_pAncestral=c(avg_pAncestral,PermAnal[[1]][["avg_pAncestral_obs"]])
}

png("Fig4f_pAncestralAmongCandidates.png",width=400,height=300)
par(mar=c(4,3,3,1),mgp=c(2,.7,0))
boxplot(as.numeric(PermpAncestral[,2])~factor(PermpAncestral[,1],levels=1:length(PermFull),labels=as.character(length(PermFull):1)),border=grey(.5),ylim=c(0,1),ylab='average pHMM',xlab="Number of GWAlpha tests intersected",main="Selective sweeps")#paste("Populations with ancestral resistance (",unlist(strsplit(PermList[1],"_"))[2],")",sep=""))
box(col="#41d894")
mtext("Populations with derived susceptibility")
lines(avg_pAncestral[order(as.numeric(names(PermFull)),decreasing=T)],lwd=1.5,col=2)
par(xpd=T)
mtext(text="F", side = 3, cex=1.5, font=2, at=-.05*par("usr")[2],line=1)
dev.off()


