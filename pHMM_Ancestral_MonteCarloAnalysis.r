#!/usr/bin/Rscript

#########################################################################################################################################
#
#                  This R script allows to replicate the trait analysis proposed in Fourneir-Level el al. (2018, submitted)
#
# - The script was intended to be run in parallel using the following command parallel -j 12 'Rscript pHMM_Ancestral_MomnteCarloAnalysis.r 0.000001 {}' ::: `echo {1..7}`
# - The input data can be sourced at the following URL:  
# - The input data folder should be located in the R working directory without altering the folder names or structure
#
#########################################################################################################################################

options(scipen=999)
args=commandArgs(TRUE)
k=as.numeric(args[1])
i=as.numeric(args[2])

load(Sys.glob(paste("PopGen/Hmm_SelSweep_k",k,".RData",sep="")))
load(Sys.glob("GWAS/Candidate_Genes.RData"))
Sweep$avg_pHMM=apply(Sweep[,-c(1:2)],1,mean)
Sweep$max_pHMM=apply(Sweep[,-c(1:2)],1,max)
Sweep$avg_pTropical=apply(Sweep[,-c(1:2)][,Pop[grepl("LA|NQ",Pop)]],1,mean)
Sweep$max_pTropical=apply(Sweep[,-c(1:2)][,Pop[grepl("LA|NQ",Pop)]],1,max)
Sweep$avg_pAncestral=apply(Sweep[,-c(1:2)][,Pop[grepl("LA|NQ1|SCA|VIC|STX",Pop)]],1,mean)
Sweep$max_pAncestral=apply(Sweep[,-c(1:2)][,Pop[grepl("LA|NQ1|SCA|VIC|STX",Pop)]],1,max)
Dm3Genes=read.table("Ref_Genome/Dmel_v5.57.genes", header=F)
Dm3Genes=Dm3Genes[Dm3Genes$V1 %in% unique(Sweep[,1]),]

#For a global set
#Global_Geneset=unlist(All_Genes)
#Global_Geneset=names(table(Global_Geneset))[table(Global_Geneset)>(i-1)] #Where you increase the stringency of the Geneset
#Global_idx=which(Dm3Genes$V10 %in% Global_Geneset)
#Global_Candy=Dm3Genes[Global_idx,]

#Sweep.Global=matrix(NA,ncol=ncol(Sweep),nrow=0)
#for (j in 1:nrow(Global_Candy)) Sweep.Global=rbind(Sweep.Global,Sweep[as.character(Sweep[,1])==Global_Candy[j,1]&Sweep[,2]>Global_Candy[j,5],][1,])
# OR for a local set
Local_Geneset=unlist(All_Genes[grepl("LA|NQ1|SCA|VIC|STX",names(All_Genes))])
Local_Geneset=names(table(Local_Geneset))[table(Local_Geneset)>(i-1)] #Where you increase the stringency of the Geneset
Local_idx=which(Dm3Genes$V10 %in% Local_Geneset)
Local_Candy=Dm3Genes[Local_idx,]

Sweep.Global=matrix(NA,ncol=ncol(Sweep),nrow=0)
for (j in 1:nrow(Local_Candy)) Sweep.Global=rbind(Sweep.Global,Sweep[as.character(Sweep[,1])==Local_Candy[j,1]&Sweep[,2]>Local_Candy[j,5],][1,])

perm.avg_pHMM=c()
perm.max_pHMM=c()
perm.avg_pTropical=c()
perm.max_pTropical=c()
perm.avg_pAncestral=c()
perm.max_pAncestral=c()
z=0
while (length(perm.avg_pHMM)<1000) {
	z=z+1
	Kut=sample(2:(nrow(Dm3Genes)-1),1)
	Dm.perm=Dm3Genes[c((Kut+1):nrow(Dm3Genes),1:Kut),]
	Dm.perm=Dm.perm[Local_idx,]
	Sweep.perm=matrix(NA,ncol=ncol(Sweep),nrow=0)
	for (j in 1:nrow(Dm.perm)) Sweep.perm=rbind(Sweep.perm,Sweep[as.character(Sweep[,1])==Dm.perm[j,1]&Sweep[,2]>Dm.perm[j,5],][1,])
	if (nrow(Sweep.perm[!is.na(Sweep.perm$avg_pHMM),])==0) next
	perm.avg_pHMM=c(perm.avg_pHMM,mean(Sweep.perm$avg_pHMM,na.rm=T))
	perm.max_pHMM=c(perm.max_pHMM,mean(Sweep.perm$max_pHMM,na.rm=T))
	perm.avg_pTropical=c(perm.avg_pTropical,mean(Sweep.perm$avg_pTropical,na.rm=T))
	perm.max_pTropical=c(perm.max_pTropical,mean(Sweep.perm$max_pTropical,na.rm=T))
	perm.avg_pAncestral=c(perm.avg_pAncestral,mean(Sweep.perm$avg_pAncestral,na.rm=T))
	perm.max_pAncestral=c(perm.max_pAncestral,mean(Sweep.perm$max_pAncestral,na.rm=T))
	if (length(perm.avg_pHMM) %% 100 ==0) print(paste(length(perm.avg_pHMM),"permutations done for Candidates shared across at least",i,"GWAlphas"))
}

if (!sum(mean(Sweep.Global$avg_pHMM,na.rm=T)>perm.avg_pHMM)) { p_avg_pHMM=0
} else p_avg_pHMM=uniroot(function(x) mean(Sweep.Global$avg_pHMM,na.rm=T)-quantile(perm.avg_pHMM,x),c(0,1))$root
if (!sum(mean(Sweep.Global$max_pHMM,na.rm=T)>perm.max_pHMM)) { p_max_pHMM=0
} else p_max_pHMM=uniroot(function(x) mean(Sweep.Global$max_pHMM,na.rm=T)-quantile(perm.max_pHMM,x),c(0,1))$root
if (!sum(mean(Sweep.Global$avg_pTropical,na.rm=T)>perm.avg_pTropical)) { p_avg_pTropical=0
} else p_avg_pTropical=uniroot(function(x) mean(Sweep.Global$avg_pTropical,na.rm=T)-quantile(perm.avg_pTropical,x),c(0,1))$root
if (!sum(mean(Sweep.Global$max_pTropical,na.rm=T)>perm.max_pTropical)) { p_max_pTropical=0
} else p_max_pTropical=uniroot(function(x) mean(Sweep.Global$max_pTropical,na.rm=T)-quantile(perm.max_pTropical,x),c(0,1))$root
if (!sum(mean(Sweep.Global$avg_pAncestral,na.rm=T)>perm.avg_pAncestral)) { p_avg_pAncestral=0
} else p_avg_pAncestral=uniroot(function(x) mean(Sweep.Global$avg_pAncestral,na.rm=T)-quantile(perm.avg_pAncestral,x),c(0,1))$root
if (!sum(mean(Sweep.Global$max_pAncestral,na.rm=T)>perm.max_pAncestral)) { p_max_pAncestral=0
} else p_max_pAncestral=uniroot(function(x) mean(Sweep.Global$max_pAncestral,na.rm=T)-quantile(perm.max_pAncestral,x),c(0,1))$root

PermAnal=list()
PermAnal[[paste(i)]]=list(avg_pHMM_obs=mean(Sweep.Global$avg_pHMM,na.rm=T),
	                      avg_pHMM_perm=perm.avg_pHMM,
	                      p_avg_pHMM=p_avg_pHMM,
	                      max_pHMM_obs=mean(Sweep.Global$max_pHMM,na.rm=T),
	                      max_pHMM_perm=perm.max_pHMM,
	                      p_max_pHMM=p_max_pHMM,
	                      avg_pTropical_obs=mean(Sweep.Global$avg_pTropical,na.rm=T),
	                      avg_pTropical_perm=perm.avg_pTropical,
	                      p_avg_pTropical=p_avg_pTropical,
	                      max_pTropical_obs=mean(Sweep.Global$max_pTropical,na.rm=T),
	                      max_pTropical_perm=perm.max_pTropical,
	                      p_max_pTropical=p_max_pTropical,
	                      avg_pAncestral_obs=mean(Sweep.Global$avg_pAncestral,na.rm=T),
	                      avg_pAncestral_perm=perm.avg_pAncestral,
	                      p_avg_pAncestral=p_avg_pAncestral,
	                      max_pAncestral_obs=mean(Sweep.Global$max_pAncestral,na.rm=T),
	                      max_pAncestral_perm=perm.max_pAncestral,
	                      p_max_pAncestral=p_max_pAncestral
	                      )

save(PermAnal,file=paste("PopGen/SweepAnc_k",paste(k),"_permutation_x",i,".RData",sep=""))

