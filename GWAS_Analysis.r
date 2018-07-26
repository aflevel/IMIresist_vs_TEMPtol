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

##################################################################################################
#        GWAlpha test
##################################################################################################

GRAPH=T
GRAPH_dia=F
Highlight_Candies=F
Int.size=2500
SNPinKB=1*Int.size/1000
All_Genes=list()
All_OUTPUT=matrix(NA,nrow=0,ncol=8)
colnames(All_OUTPUT)=c("Assay","Chromosome","Position","Mutation","GWAlpha","MAF","COV","Gene_Name")
Pheno.name="STX-20C_All"

Pheno.lab=list('STX-20C_All'=expression(paste("South Texas (USA, hot and dry) at ",20*degree*C)),
               'STX-30C_All'=expression(paste("South Texas (USA, hot and dry) at ",30*degree*C)),
               'NTX-20C_All'=expression(paste("North Texas (USA, hot and dry) at ",20*degree*C)),
               'NTX-30C_All'=expression(paste("North Texas (USA, hot and dry) at ",30*degree*C)),
               'NLA-20C_All'=expression(paste("North Louisiana (USA, hot and wet) at ",20*degree*C)),
               'NLA-30C_All'=expression(paste("North Louisiana (USA, hot and wet) at ",30*degree*C)),
               'SLA-20C_All'=expression(paste("South Louisiana (USA, hot and wet) at ",20*degree*C)),
               'SLA-30C_All'=expression(paste("South Louisiana (USA, hot and wet) at ",30*degree*C)),
               'NCA-20C_All'=expression(paste("North California 1 (USA, cold and dry) at ",20*degree*C)),
               'NCA-30C_All'=expression(paste("North California 1 (USA, cold and dry) at ",30*degree*C)),
               'SCA-20C_All'=expression(paste("North California 2 (USA, cold and dry) at ",20*degree*C)),
               'SCA-30C_All'=expression(paste("North California 2 (USA, cold and dry) at ",30*degree*C)),
               'NB-20C_All'=expression(paste("New Brunswick (CAN, cold and wet) at ",20*degree*C)),
               'NB-30C_All'=expression(paste("New Brunswick (CAN, cold and wet) at ",30*degree*C)),
               'NS-20C_All'=expression(paste("Nova Scotia (CAN, cold and wet) at ",20*degree*C)),
               'NS-30C_All'=expression(paste("Nova Scotia (CAN, cold and wet) at ",30*degree*C)),
               'STS-20C_All'=expression(paste("South Tasmania (AUS, cold and wet) at ",20*degree*C)),
               'STS-30C_All'=expression(paste("South Tasmania (AUS, cold and wet) at ",30*degree*C)),
               'NTS-20C_All'=expression(paste("North Tasmania (AUS, cold and wet) at ",20*degree*C)),
               'NTS-30C_All'=expression(paste("North Tasmania (AUS, cold and wet) at ",30*degree*C)),
               'NQ1-20C_All'=expression(paste("North Queensland 1 (AUS, hot and wet) at ",20*degree*C)),
               'NQ1-30C_All'=expression(paste("North Queensland 1 (AUS, hot and wet) at ",30*degree*C)),
               'NQ2-20C_All'=expression(paste("North Queensland 2 (AUS, hot and wet) at ",20*degree*C)),
               'NQ2-30C_All'=expression(paste("North Queensland 2 (AUS, hot and wet) at ",30*degree*C)),
               'SQ1-20C_All'=expression(paste("South Queensland 1 (AUS, hot and wet) at ",20*degree*C)),
               'SQ1-30C_All'=expression(paste("South Queensland 1 (AUS, hot and wet) at ",30*degree*C)),
               'SQ2-20C_All'=expression(paste("South Queensland 2 (AUS, hot and wet) at ",20*degree*C)),
               'SQ2-30C_All'=expression(paste("South Queensland 2 (AUS, hot and wet) at ",30*degree*C)),
               'SA-20C_All'=expression(paste("South Australia (AUS, cold and dry) at ",20*degree*C)),
               'SA-30C_All'=expression(paste("South Australia (AUS, cold and dry) at ",30*degree*C)),
               'VIC-20C_All'=expression(paste("Victoria (AUS, cold and dry) at ",20*degree*C)),
               'VIC-30C_All'=expression(paste("Victoria (AUS, cold and dry) at ",30*degree*C))
               )

PHENO=c("NLA-20C_All","NLA-30C_All","SLA-20C_All","SLA-30C_All","NTX-20C_All","STX-20C_All",
        "NTX-30C_All","STX-30C_All","NB-20C_All","NB-30C_All","NS-20C_All","NS-30C_All",
        "SCA-20C_All","SCA-30C_All","NCA-20C_All","NCA-30C_All","NQ1-30C_All","NQ2-30C_All",
        "SQ1-30C_All","SQ2-30C_All","NTS-30C_All","STS-30C_All","SA-30C_All","VIC-30C_All",
        "NQ1-20C_All","NQ2-20C_All","SQ1-20C_All","SQ2-20C_All","NTS-20C_All","STS-20C_All",
        "SA-20C_All","VIC-20C_All")

for (Pheno.name in PHENO) {

GWAlpha=read.csv(paste("GWAS/GWAlpha_",Pheno.name,"_out.csv",sep=""),header=T)
GWAlpha=GWAlpha[GWAlpha[,1]!=4,] # The chromosome 4 was removed from the manhattan plots but was retained for defining the candidate genes set

if ("2L" %in% GWAlpha[,1]) {
CHR.labs=c("5"="4","1"="2L","2"="2R","3"="3L","4"="3R","6"="X")
GWAlpha=GWAlpha[order(as.character(GWAlpha[,1]),as.numeric(as.character(GWAlpha[,2]))),]
} else {
CHR.labs=unique(as.character(GWAlpha[,1]))
names(CHR.labs)=unique(as.character(GWAlpha[,1]))
GWAlpha=GWAlpha[order(as.character(GWAlpha[,1]),as.numeric(as.character(GWAlpha[,2]))),]
}

OUT=GWAlpha[GWAlpha[,1] %in% CHR.labs,]
CHR=OUT[,1]

for (i in 1:length(CHR.labs)) CHR=gsub(CHR.labs[i],names(CHR.labs)[i],CHR)
CHR.labs=CHR.labs[names(CHR.labs) %in% unique(CHR)]
CHR.labs=CHR.labs[order(names(CHR.labs))]

Position=as.numeric(as.character(OUT[,2]))
CHR.start=vector()
CHR.end=vector()
for (i in unique(CHR)) {
	ID=which(CHR==i)
	CHR.start=c(CHR.start,ID[1])
	CHR.end=c(CHR.end,ID[length(ID)])
}
CHR.summary=cbind(Position[CHR.start],Position[CHR.end])
rownames(CHR.summary)=unique(CHR)

Position.cum=CHR.summary[,2]-CHR.summary[,1]
CHR.cum=0
for (i in 2:length(unique(CHR))) CHR.cum=c(CHR.cum,sum(Position.cum[1:(i-1)]))
names(CHR.cum)=unique(CHR)

Position.cum=vector()
for (i in unique(CHR)) Position.cum=c(Position.cum,Position[which(CHR==i)]+as.numeric(CHR.cum[which(names(CHR.cum)==i)])-1)

CHR.mid=CHR.cum+(c(CHR.cum[-1],max(CHR.cum)+Position[length(Position)])-CHR.cum)/2

CHR=as.numeric(as.factor(CHR))
GWAlpha=abs(as.numeric(as.character(OUT$Alpha)))

GWAlpha.order=GWAlpha[order(GWAlpha)]
GWAlpha.QQ=predict(smooth.spline(GWAlpha.order,spar=.6),deriv=2)$y
Thresh=GWAlpha.order[which(GWAlpha.QQ==max(GWAlpha.QQ))]
p_alpha=sum(GWAlpha>Thresh)/length(GWAlpha)

if (GRAPH_dia==T) {
	dev.new()
	plot(GWAlpha[order(GWAlpha)])
	abline(h=Thresh,col=3)

	dev.new()
	par(mfrow=c(2,1),mgp=c(2,1,0))
}
	SNP_density=hist(Position.cum,breaks=round(max(Position.cum)/Int.size)+1,plot=GRAPH_dia)
	Peak=hist(Position.cum[GWAlpha>Thresh],breaks=round(max(Position.cum)/Int.size)+1,plot=GRAPH_dia)
	Peak.Inter=SNP_density
	Peak.Inter$counts=rep(0,length(Peak.Inter$counts))
	Peak.Inter$counts[SNP_density$mids %in% Peak$mids]=Peak$counts
	Peak=Peak.Inter
	Peak.position=Peak$mids[(Peak$counts>=SNPinKB)&Peak$counts/SNP_density$counts>=p_alpha]
	SNP_Peak=SNP_density$counts[Peak$count>=SNPinKB]

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

if (Highlight_Candies) {
	load(Sys.glob("Ref_Genome/Candidate_Genes.RData"))
	Dm3Genes=read.table("Ref_Genome/Dmel_v5.57.genes", header=F)
	Dm3Genes=Dm3Genes[Dm3Genes$V1 %in% unique(OUT[,1]),]

	Candies=names(table(unlist(All_Genes[[Pheno.name]]))) #Where you can focus on a single Geneset
#OR
	#Candies=names(table(unlist(All_Genes)))[table(unlist(All_Genes))>8] #Where you increase the stringency of the Geneset
	names(Candies)=Candies
	Candies_idx=which(Dm3Genes$V10 %in% Candies)
	Candidates=Dm3Genes[Candies_idx,c(1,4,5,10)]
} else {
	Dm3Genes=read.table("Ref_Genome/Dmel_v5.57.genes", header=F)
	Candies=c("Prm","nAChRalpha3")#,"Cyp6g1","Pde9","Snap25")
	names(Candies)=c("Prm","nAChRalpha3")#,"Cyp6g1","Pde9","Snap25")
	Candidates=Dm3Genes[Dm3Genes$V10 %in% Candies,c(1,4,5,10)]
}

p_draw=p_score
p_draw[(p_draw<1)&(Position.cum %%2 == 0)]=NA
if (GRAPH==T) {
	png(file=paste("GWAlpha_",Pheno.name,"_Highlight.png",sep=""),width=1200,height=350)
	par(mar=c(3.5,3.7,2.5,2),mgp=c(1.8,0.6,0),cex.axis=1.5)
	plot(0,pch='',xlim=c(0,max(Position.cum)),xlab="",xaxt='n',ylab="",ylim=c(0,max(p_score,na.rm=T)))
	mtext(text=Pheno.lab[[Pheno.name]], side = 3, cex=1.6, font=2, adj=0)
	title(xlab="Genomic position", ylab=expression(-log[10](italic(p))),cex.lab=1.5)
	for (i in 1:length(Candies)) {
		CandiPos=Position.cum[(as.character(OUT[,1]) == Candidates[i,1]) & (as.numeric(as.character(OUT[,2])) < Candidates[i,2])]
		abline(v=CandiPos[length(CandiPos)],col=grey(.9),lwd=10)
		if (Candidates[i,4]=="nAChRalpha3") text(CandiPos[length(CandiPos)],.85*par("usr")[4],expression(paste("nAChR",alpha,3)),srt=45)
		else text(CandiPos[length(CandiPos)],.85*par("usr")[4],as.expression(paste(names(Candies)[Candies==Candidates[i,4]])),srt=45)
	}
	points(p_draw[CHR %% 2 !=0]~Position.cum[CHR %% 2 !=0],pch=20, cex=1.5, col="#0e58a0")
	points(p_draw[CHR %% 2 ==0]~Position.cum[CHR %% 2 ==0],pch=20, cex=1.5, col="#ADD8E6")
	points(p_draw[OUT$Allele=="TE"]~Position.cum[OUT$Allele=="TE"],pch=20, cex=1.5, col="#b642f4")
	points(p_draw[OUT$Allele=="INV"]~Position.cum[OUT$Allele=="INV"],pch=20, cex=1.5, col="#2eb279")
	abline(v=CHR.cum[-1],lty=2)
	abline(h=-log10(pThresh),col=2)
	axis(side=1,at=CHR.mid,labels=CHR.labs,tick=F)
	box()
	dev.off()
}

#Using the peak calling
Loci.list=matrix(NA,ncol=ncol(OUT),nrow=0)
for (p in Peak.position) Loci.list=rbind(Loci.list,OUT[Position.cum>(p-Int.size/2) & Position.cum<(p+Int.size/2) & GWAlpha>Thresh,])
#or using BH-corrected p-values
#Loci.list=OUT[pval<pThresh,]
#if (nrow(Loci.list)==0) next

#exclude Prm and nACHRalpha3 regions
#Loci.list=Loci.list[!(Loci.list$Chromosome=="3L" & Loci.list$Position>8700000 & Loci.list$Position<8750000),]
#Loci.list=Loci.list[!(Loci.list$Chromosome=="X" & Loci.list$Position>8350000 & Loci.list$Position<8400000),]


##################################################################################################
#        Gene List inferences
##################################################################################################

chr.gene=which(as.character(OUT[,1]) == Candidates[1,1])
pos.gene=which(OUT[,2] < Candidates[1,2])
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

OUTPUT=cbind(Pheno.name,Loci.list,Gene.list)
colnames(OUTPUT)=colnames(All_OUTPUT)
All_OUTPUT=rbind(All_OUTPUT,OUTPUT)
print(Pheno.name)
print(table(gsub(" ","",unlist(strsplit(as.character(OUTPUT[,"Gene_Name"]),",")))))
#print("PRM")
#PRM=OUT[OUT[,1]=="3L" & Position>8725958 & Position<8738410 & abs(GWAlpha)>Thresh,]
#if (nrow(PRM)>0) print(PRM)
#print("nAChRalpha3")
#nAChRalpha3=OUT[OUT[,1]=="X" & Position>8170233 & Position<8335000 & abs(GWAlpha)>Thresh,]
#if (nrow(nAChRalpha3)>0) print(nAChRalpha3)
#write.csv(OUTPUT,paste("GWAList_",Pheno.name,".csv",sep=""),row.names=F)
#print("TEs")
#TE=OUT[OUT$Allele=="TE" & GWAlpha>Thresh,]
#if (nrow(TE)>0) print(TE)

All_Genes[[Pheno.name]]=unique(gsub(" ","",unlist(strsplit(as.character(OUTPUT[OUTPUT[,"Gene_Name"]!="Intergenic","Gene_Name"]),","))))
}

#Create a candidate gene list 
save(All_Genes,file="Candidate_Genes.RData")

#export gene list, with or without singletons 
#the candidate gene list was generated including Chr 4 (comment out line 70) and excluding singletons
#write.table(unique(unlist(All_Genes)),"GWAlist_dm5.57.txt",row.names=F,quote=F,col.names=F)
#OR
write.table(names(table(unlist(All_Genes)))[table(unlist(All_Genes))>1],"GWAlist_dm5.57.txt",row.names=F,quote=F,col.names=F)


############################################################################################################################################
#
#                                                  Effect of CNV on resistance
#
############################################################################################################################################

GRAPH=F
CNV_Genes=list()
Pheno.name="STX-20C_All"

Pheno.lab=list('STX-20C_All'=expression(paste("South Texas (USA, hot and dry) at ",20*degree*C)),
               'STX-30C_All'=expression(paste("South Texas (USA, hot and dry) at ",30*degree*C)),
               'NTX-20C_All'=expression(paste("North Texas (USA, hot and dry) at ",20*degree*C)),
               'NTX-30C_All'=expression(paste("North Texas (USA, hot and dry) at ",30*degree*C)),
               'NLA-20C_All'=expression(paste("North Louisiana (USA, hot and wet) at ",20*degree*C)),
               'NLA-30C_All'=expression(paste("North Louisiana (USA, hot and wet) at ",30*degree*C)),
               'SLA-20C_All'=expression(paste("South Louisiana (USA, hot and wet) at ",20*degree*C)),
               'SLA-30C_All'=expression(paste("South Louisiana (USA, hot and wet) at ",30*degree*C)),
               'NCA-20C_All'=expression(paste("North California 1 (USA, cold and dry) at ",20*degree*C)),
               'NCA-30C_All'=expression(paste("North California 1 (USA, cold and dry) at ",30*degree*C)),
               'SCA-20C_All'=expression(paste("North California 2 (USA, cold and dry) at ",20*degree*C)),
               'SCA-30C_All'=expression(paste("North California 2 (USA, cold and dry) at ",30*degree*C)),
               'NB-20C_All'=expression(paste("New Brunswick (CAN, cold and wet) at ",20*degree*C)),
               'NB-30C_All'=expression(paste("New Brunswick (CAN, cold and wet) at ",30*degree*C)),
               'NS-20C_All'=expression(paste("Nova Scotia (CAN, cold and wet) at ",20*degree*C)),
               'NS-30C_All'=expression(paste("Nova Scotia (CAN, cold and wet) at ",30*degree*C)),
               'STS-20C_All'=expression(paste("South Tasmania (AUS, cold and wet) at ",20*degree*C)),
               'STS-30C_All'=expression(paste("South Tasmania (AUS, cold and wet) at ",30*degree*C)),
               'NTS-20C_All'=expression(paste("North Tasmania (AUS, cold and wet) at ",20*degree*C)),
               'NTS-30C_All'=expression(paste("North Tasmania (AUS, cold and wet) at ",30*degree*C)),
               'NQ1-20C_All'=expression(paste("North Queensland 1 (AUS, hot and wet) at ",20*degree*C)),
               'NQ1-30C_All'=expression(paste("North Queensland 1 (AUS, hot and wet) at ",30*degree*C)),
               'NQ2-20C_All'=expression(paste("North Queensland 2 (AUS, hot and wet) at ",20*degree*C)),
               'NQ2-30C_All'=expression(paste("North Queensland 2 (AUS, hot and wet) at ",30*degree*C)),
               'SQ1-20C_All'=expression(paste("South Queensland 1 (AUS, hot and wet) at ",20*degree*C)),
               'SQ1-30C_All'=expression(paste("South Queensland 1 (AUS, hot and wet) at ",30*degree*C)),
               'SQ2-20C_All'=expression(paste("South Queensland 2 (AUS, hot and wet) at ",20*degree*C)),
               'SQ2-30C_All'=expression(paste("South Queensland 2 (AUS, hot and wet) at ",30*degree*C)),
               'SA-20C_All'=expression(paste("South Australia (AUS, cold and dry) at ",20*degree*C)),
               'SA-30C_All'=expression(paste("South Australia (AUS, cold and dry) at ",30*degree*C)),
               'VIC-20C_All'=expression(paste("Victoria (AUS, cold and dry) at ",20*degree*C)),
               'VIC-30C_All'=expression(paste("Victoria (AUS, cold and dry) at ",30*degree*C))
               )

PHENO=c("NLA-20C_All","NLA-30C_All","SLA-20C_All","SLA-30C_All","NTX-20C_All","STX-20C_All",
        "NTX-30C_All","STX-30C_All","NB-20C_All","NB-30C_All","NS-20C_All","NS-30C_All",
        "SCA-20C_All","SCA-30C_All","NCA-20C_All","NCA-30C_All","NQ1-30C_All","NQ2-30C_All",
        "SQ1-30C_All","SQ2-30C_All","NTS-30C_All","STS-30C_All","SA-30C_All","VIC-30C_All",
        "NQ1-20C_All","NQ2-20C_All","SQ1-20C_All","SQ2-20C_All","NTS-20C_All","STS-20C_All",
        "SA-20C_All","VIC-20C_All")

for (Pheno.name in PHENO) {


#CNV summary statistics were computed from the pileup file using the CNVkit package 

if (file.exists(paste("CNVtest/",gsub("All","CNV.txt",Pheno.name),sep=""))) {
	CNV=read.table(paste("CNVtest/",gsub("All","CNV.txt",Pheno.name),sep=""),header=T)
} else {
CNV.files=list.files(paste("CNVtest/",unlist(strsplit(Pheno.name,"-"))[1],sep=""),pattern="trim.cnr")
CNV.files=CNV.files[grepl(substr(unlist(strsplit(Pheno.name,"-"))[2],1,3),CNV.files)]

CNV=read.table(paste("CNVtest/",unlist(strsplit(Pheno.name,"-"))[1],"/",CNV.files[1],sep=""),header=T)
CNV=CNV[,c(1:3,6:7)]
for (f in CNV.files[-1]) CNV=cbind(CNV,read.table(paste("CNVtest/",unlist(strsplit(Pheno.name,"-"))[1],"/",f,sep=""),header=T)[,c(6:7)])
Reg=apply(CNV[,names(CNV)=="log2"],1,function(x) {summary(lm(x~c(4,20,62.5,105,122)))$coefficients[2,1]})
CNV=cbind(CNV,Reg)
CNV=CNV[order(CNV[,1]),]
write.table(CNV,paste("CNVtest/",gsub("All","CNV.txt",Pheno.name),sep=""),row.names=F,quote=F,sep="\t")
}
#}

LL=function(x,y) sum(dnorm(y,x[1],x[2],log=T))
EST=optim(par=c(mu="0",sd=sd(CNV$Reg)),fn=LL,y=CNV$Reg,control=list(fnscale=-1,reltol=1e-10))$par
pval=ifelse(pnorm(CNV$Reg,EST[1],EST[2]) > 0.5, 2*(1-pnorm(CNV$Reg,EST[1],EST[2])), 2*pnorm(CNV$Reg,EST[1],EST[2]))
p_score=-log10(pval)
Padj=p.adjust(pval,method="BH")
q=0.05
if (sum(Padj<q)!=0) {
	pThresh=max(pval[Padj<q])
} else pThresh=min(pval)/2

if (GRAPH==T) {
	#png(file=paste("GWA_CNV_",Pheno.name,"_Highlight.png",sep=""),width=1200,height=350)
	par(mar=c(3.5,3.7,2.5,2),mgp=c(1.8,0.6,0))
	plot('n',ylim=c(0,max(p_score)),xlim=c(1,nrow(CNV)),xaxt='n',xlab="",ylab="")
	points(p_score[as.numeric(CNV[,1]) %% 2 !=0]~c(1:nrow(CNV))[as.numeric(CNV[,1]) %% 2 !=0],pch=20, col="#104E8B")
	points(p_score[as.numeric(CNV[,1]) %% 2 ==0]~c(1:nrow(CNV))[as.numeric(CNV[,1]) %% 2 ==0],pch=20, col="#ADD8E6")
	mtext(text=Pheno.lab[[Pheno.name]], side = 3, cex=1.6, font=2, adj=0)
	title(xlab="Genomic position", ylab="-log(p) log-copy number",cex.lab=1.5)
	X=c(1,which(c(CNV[1,1],CNV[,1]) != c(CNV[,1],CNV[nrow(CNV),1])))
	BK=which(c(CNV[1,1],CNV[,1]) != c(CNV[,1],CNV[nrow(CNV),1]))-1
	CHR=CNV[c(1,which(c(CNV[1,1],CNV[,1]) != c(CNV[,1],CNV[nrow(CNV),1]))),1]
	axis(side=1,at=X+(c(BK,nrow(CNV))-c(1,BK))/2,par("usr")[3]+.1*abs(par("usr")[3]),labels=CHR,las=1,tick=F)
	abline(v=which(c(CNV[1,1],CNV[,1]) != c(CNV[,1],CNV[nrow(CNV),1])),lty=2)
	abline(h=-log10(pThresh),col=2)
	#dev.off()
}

OUT=CNV[p_score>-log10(pThresh),c(1:3,14)]
Dm3Genes=read.table("Ref_Genome/Dmel_v5.57.genes", header=F)
Gene.list=c()
for (i in 1:nrow(OUT)) {
	Dm3.sub=Dm3Genes[as.character(Dm3Genes$V1)==as.character(OUT[i,1]) & ((OUT[i,2]<Dm3Genes$V4 & OUT[i,3]>Dm3Genes$V4) | 
		                                     (OUT[i,2]<Dm3Genes$V5 & OUT[i,2]>Dm3Genes$V5) | 
		                                     (OUT[i,2]>Dm3Genes$V4 & OUT[i,3]<Dm3Genes$V5)),10]
	if (length(Dm3.sub)>0) Gene.list=c(Gene.list,as.character(Dm3.sub))
}

CNV_Genes[[Pheno.name]]=unique(Gene.list)

#load(Sys.glob("Candidate_Genes.RData"))
#for (i in names(CNV_Genes)) {
#	CNV_over=CNV_Genes[[i]][CNV_Genes[[i]] %in% All_Genes[[gsub("All","trim2",i)]]]
#	print(paste(i,":",length(CNV_over),"/",length(All_Genes[[gsub("All","trim2",i)]])))
#}
Dm3Genes=read.table("Ref_Genome/Dmel_v5.57.genes", header=F)
Candies=c("Prm","nAChRalpha3","Ir41a","kirre","rg")#,"Cyp6g1","Pde9","Snap25")
names(Candies)=c("Prm","nAChRalpha3","Ir41a","kirre","rg")#,"Cyp6g1","Pde9","Snap25")
Candidates=Dm3Genes[Dm3Genes$V10 %in% Candies,c(1,4,5,10)]

#dev.new()
png(paste("GWA_CNV_",unlist(strsplit(Pheno.name,"_"))[1],".png",sep=""),width=1200,height=400)
par(mar=c(3.5,3.7,2.5,2),mgp=c(1.8,0.6,0))
plot('n',xlim=c(0,nrow(CNV)),ylim=c(min(CNV$Reg),max(CNV$Reg)),col=grey(.6),xaxt='n',xlab="",ylab="")
mtext(text=Pheno.lab[[Pheno.name]], side = 3, cex=1.6, font=2, adj=0)
title(xlab="Genomic position", ylab="effect of log-copy number",cex.lab=1.5)
X=c(1,which(c(CNV[1,1],CNV[,1]) != c(CNV[,1],CNV[nrow(CNV),1])))
BK=which(c(CNV[1,1],CNV[,1]) != c(CNV[,1],CNV[nrow(CNV),1]))-1
CHR=CNV[c(1,which(c(CNV[1,1],CNV[,1]) != c(CNV[,1],CNV[nrow(CNV),1]))),1]
axis(side=1,at=X+(c(BK,nrow(CNV))-c(1,BK))/2,par("usr")[3]+.1*abs(par("usr")[3]),labels=CHR,las=1,tick=F,cex=1.5)
for (i in 1:nrow(Candidates)) {
	CandiPos=round(mean(which((as.character(CNV[,1]) == Candidates[i,1]) & (as.numeric(as.character(CNV[,3])) > Candidates[i,2]-5000) & (as.numeric(as.character(CNV[,2])) < Candidates[i,3]+5000))))
	abline(v=CandiPos,col=grey(.9),lwd=10)
	if (Candidates[i,4]=="nAChRalpha3") text(CandiPos,.85*par("usr")[4],expression(paste("nAChR",alpha,3)),srt=45)
	else text(CandiPos,.85*par("usr")[4],as.expression(paste(names(Candies)[Candies==Candidates[i,4]])),srt=45)
}
points(CNV$Reg,col=grey(.6))
abline(v=which(c(CNV[1,1],CNV[,1]) != c(CNV[,1],CNV[nrow(CNV),1])),lty=2)
smoo <- smooth.spline(CNV$Reg,spar=0.03)
lines(3*predict(smoo)$y,col="#0e58a0")
points(CNV$Reg[pval<pThresh]~which(pval<pThresh),pch=19)
abline(h=c(min(CNV$Reg[pval<pThresh&CNV$Reg>0]),max(CNV$Reg[pval<pThresh&CNV$Reg<0])),col=2)
box()
dev.off()
if (sum(pval==pThresh)>0) { print(Pheno.name)
} else print(paste("no exact pval match for",Pheno.name))

}

#Update the Cnadidate Genes list with the CNV association candidates
load("Candidate_Genes.RData")
for (i in PHENO) All_Genes[[i]]=c(All_Genes[[i]],CNV_Genes[[i]])
save(All_Genes,file="Candidate_Genes.RData")


############################################################################################################################################
#
#                                                  Effect of Endo on resistance
#
############################################################################################################################################

Endo=read.table(Sys.glob(paste("ENDOtest/endosymb_presence.txt")))
Genome=Endo[Endo$V4=="genome",]

#Since it is needed Coverage summary statistic
AvCov=c()
for (i in unique(Genome$V1)) {
	Genome.sub=Genome[Genome$V1==i,]
	AvCov=c(AvCov,mean(Genome.sub$V5))
}
names(AvCov)=unique(Genome$V1)
as.matrix(AvCov,ncol=1)

Endo=Endo[Endo$V4!="genome",]

out=matrix(NA,nrow=0,ncol=6)
colnames(out)=c("Pop","Temp","Endo","dose","coef","p-val")
for (i in unique(paste(Endo$V1,Endo$V2,Endo$V4))) {
	Endo.sub=Endo[paste(Endo$V1,Endo$V2,Endo$V4)==i,]
	Genome.sub=Genome[paste(Genome$V1,Genome$V2)==paste(unlist(strsplit(i," "))[1:2],collapse=" "),]
	dose=mean(Endo.sub$V5/Genome.sub$V5*c(8,24,61,24,8)/125)
	test=summary(lm(Endo.sub$V5/Genome.sub$V5~c(4,20,62.5,105,122)))
	out=rbind(out,c(unlist(strsplit(i," ")),dose,test$coefficients[2,1],test$coefficients[2,4]))
}

out[out[,6]<0.01,]

Endo=matrix(NA,nrow=length(unique(out[,1])),ncol=length(unique(out[,3])))
rownames(Endo)=unique(out[,1])
colnames(Endo)=unique(out[,3])
for (i in unique(out[,1])) for (j in unique(out[,3])) Endo[i,j]=mean(as.numeric(out[,4][out[,1]==i & out[,3]==j]))
print(Endo)



##################################################################################################
#                                        Density map of the Candidate Genes
##################################################################################################

ChrPos=read.table(Sys.glob(paste("PopGen/",list.files(path=Sys.glob("PopGen"),pattern="npstats"),sep=""))[1],header=T)[,1:2]
load(Sys.glob("PopGen/Fst_Mantel.RData"))
load(Sys.glob("Candidate_Genes.RData"))
Dm3Genes=read.table("Ref_Genome/Dmel_v5.57.genes", header=F)
Dm3Genes=Dm3Genes[Dm3Genes$V1 %in% unique(Pop[,1]),]

png(file="Fig4a_CandidatesDensity.png",width=1200,height=100)
par(mar=c(3,9,2,5),mgp=c(1.8,0.6,0),xaxs="i")
plot('n',xlim=c(0,nrow(ChrPos)),ylim=c(0,1),yaxt='n',ylab="",xaxt='n',xlab="Genomic Position")
mtext(text="Density of Candidate Genes", side = 3, cex=1.3, font=2, adj=0)

Couleur=colorRampPalette(c("#f2eded", "pink","red","#751919","black"))(15)


for (i in 2:16) {
Global_Geneset=unlist(All_Genes)
Global_Geneset=names(table(Global_Geneset))[table(Global_Geneset)>(i-1)] #Where you increase the stringency of the Geneset
Global_idx=which(Dm3Genes$V10 %in% Global_Geneset)
Global_Candy=Dm3Genes[Global_idx,]

Candy.matrix=matrix(0,nrow=nrow(ChrPos),ncol=1)
Chr_Candy=Global_Candy$V1
Pos_Candy=round(Global_Candy[,5]/ChrPos[1,2])*ChrPos[1,2]
Candy.matrix[ChrPos[,1] %in% Chr_Candy & ChrPos[,2] %in% Pos_Candy,1]=1

abline(v=which(Candy.matrix[,1]==1),col=Couleur[i-1])
}

X=c(1,which(c(ChrPos[1,1],ChrPos[,1]) != c(ChrPos[,1],ChrPos[nrow(ChrPos),1])))
BK=which(c(ChrPos[1,1],ChrPos[,1]) != c(ChrPos[,1],ChrPos[nrow(ChrPos),1]))-1
abline(v=BK,col="blue")
CHR=ChrPos[c(1,which(c(ChrPos[1,1],ChrPos[,1]) != c(ChrPos[,1],ChrPos[nrow(ChrPos),1]))),1]
axis(side=1,at=(X+(c(BK,nrow(ChrPos))-c(1,BK))/2),(par("usr")[3]+.1*abs(par("usr")[3])),labels=CHR,tick=F,las=1)
axis(side=1,at=BK,(par("usr")[3]+.1*abs(par("usr")[3])),labels=F,col="blue",las=1)
box()
par(xpd=T)
mtext(text="A", side = 3, cex=1.5, font=2, at=-.1*par("usr")[2])
dev.off()

