#!/usr/bin/Rscript
options(warn=-1)
#########################################################################################################################################
#
#                  This R script allows to replicate the trait analysis proposed in Fourneir-Level el al. (2018, submitted)
#
# - It requires the lme4, FactoMineR and plotrix packages to be installed to perform the mixed-linear modelling
# - The input data can be sourced at the following URL:  
# - The input data folder should be located in the R working directory without altering the folder names or structure
# -lines 46 64-66 and 72 can be commented out to perform a maximum likelihood estimation of mean, std dev. minimun and maximum for the insecticide data only (no controls) that is required for the GWAlpha association analysis
#########################################################################################################################################


rm(list=ls(all=T))
options(help_type = "html")
library(lme4)
library(compiler)
setwd("Pheno_data")

File.List=list.files(pattern=".csv")
File.List=File.List[!grepl("Asso|CMH|MU|IMI|SeqPoly|mapped|CG_Synonyms|GWAlpha",File.List)]

par(mfrow=c(3,4))

Data.full=data.frame(Quantile=NA,Time=NA,Temp=NA,Pop=NA,Geno=NA)[-1,]
MLest=list()
LogLik=list()

#########################################################################################################################################
#
#                                     Estimates of mean/sd using normal distribution
#
#########################################################################################################################################

k=0
for (f in File.List) {
	k=k+1
	if ((k-1) %% 12 ==0) {
		dev.new()
		par(mfrow=c(3,4))
	}
	Data=read.csv(f,header=T,na.strings="NA",sep=";")
	Data$Geno=gsub("i","I",Data$Geno)
	t.start=Data$Hour[1]+Data$Minute[1]/60
	#lines 46 64-66 and 72 can be commented out to perform a maximum likelihood estimation of mean, std dev. minimun and maximum for the insecticide data only (no controls) that is required for the GWAlpha association analysis
	#Data=Data[which(Data$Num!="CT"),]
	Data=Data[which(!is.na(Data$Day)),]
	Time=24*Data$Day+Data$Hour+Data$Minute/60-t.start
	Data=data.frame(Quantile=Data$Percentile,Time,Temp=unlist(strsplit(f,"_"))[2],Pop=Data$Pop,Geno=paste(Data$Pop,Data$Geno,sep="_"))
	liveCT=unique(Data$Geno)[-1][!(unique(Data$Geno)[-1] %in% unique(Data$Geno[Data$Quantile=="CT"])[-1])]
	if (Data$Pop[2]=="STJ") {
		liveCT=c(rep(unique(Data$Geno)[-1][!(unique(Data$Geno)[-1] %in% unique(Data$Geno[Data$Quantile=="CT"])[-1])],2),
	             names(table(Data$Geno[Data$Quantile=="CT"]))[table(Data$Geno[Data$Quantile=="CT"])<2])
	}
	extraCT=data.frame(Quantile=rep("CT",length(liveCT)),
		               Time=as.numeric(rep(max(Data$Time)+12,length(liveCT))),
		               Temp=as.character(rep(Data$Temp[1],length(liveCT))),
		               Pop=as.character(rep(Data$Pop[10],length(liveCT))),
		               Geno=as.character(liveCT))
	Data.full=rbind(Data.full,Data,extraCT)
	t.max=Time[order(Time)]
	t.min=0
	for (time in t.max[-1]) t.min=c(t.min,unique(c(0,t.max))[which(unique(c(0,t.max))==time)-1])
	#lambda=function(mu, sigma, left, right) {sum(log(pnorm(right, mu, sigma) - pnorm(left, mu, sigma)))}
	#MLest[[f]]=optim(c(mu=100,sig=50),fn=function(theta) -lambda(theta[1], theta[2], t.min, t.max))$par
	#LogLik[[f]]=lambda(MLest[[f]][1], MLest[[f]][2], t.min, t.max)
	Table.Mort=table(t.max)
	Cumul.Mort=vector()
	for (i in 1:length(Table.Mort)) Cumul.Mort=c(Cumul.Mort,sum(Table.Mort[1:i])/sum(Table.Mort))
	plot(Cumul.Mort~as.numeric(names(table(t.max))), ylab="Cumulative mortality",xlab="time in hours",main=f,xlim=c(0,400))
	x=seq(from=t.min[1],to=t.max[length(t.max)],by=1)
	#curve(pnorm(x, MLest[[f]][1], MLest[[f]][2]), from=min(x), to=max(x), col=2, lwd=2,add=TRUE)
}

Data.full$Geno=gsub("i","I",Data.full$Geno)
Data.full=Data.full[Data.full$Geno!="4"&Data.full$Geno!="8",]

Pop=Data.full$Pop
Pop.dict = setNames(c("NLA","SLA","SLA","NTX","STX","NB","NS","NS","NCA","NCA","SCA","SCA",
                      "NQ1","NQ2","SQ1","SQ2","SQ2","NTS","STS","SA","VIC"),
                    c("LDY","BTR","OLS","LBK","BEC","STJ","BBR","GAS","HUS","HDL","FST","ELK",
                      "CHI","CHE","ISG","ADD","SLV","CFY","CGC","MTG","HMT"))

for (i in 1:length(Pop.dict)) Pop=gsub(names(Pop.dict)[i],Pop.dict[i],Pop)
Data.full$Pop=Pop

Clim=Pop
Clim.dict = setNames(c("HW","HW","HD","HD","CW","CW","CD","CD",
                      "HW","HW","HD","HD","CW","CW","CD","CD"),
                     c("NLA","SLA","NTX","STX","NB","NS","NCA","SCA",
                      "NQ1","NQ2","SQ1","SQ2","NTS","STS","SA","VIC"))
for (i in 1:length(Clim.dict)) Clim=gsub(names(Clim.dict)[i],Clim.dict[i],Clim)

Cont=Pop
Cont.dict = setNames(c("USA","USA","USA","USA","USA","USA","USA","USA",
                      "AUS","AUS","AUS","AUS","AUS","AUS","AUS","AUS"),
                     c("NQ1","NQ2","SQ1","SQ2","NTS","STS","SA","VIC",
                      "NLA","SLA","NTX","STX","NB","NS","NCA","SCA"))
for (i in 1:length(Cont.dict)) Cont=gsub(names(Cont.dict)[i],Cont.dict[i],Cont)
Cont=gsub("UUSA","USA",Cont)

for (i in dev.list()) dev.off(i)
setwd("../")

#########################################################################################################################################
#
#                                         Variance partionning in the observed data
# - To reproduce the complete trait analysis insecticide AND control, lines 46 64-66 and 72 need to be commented out
#
#########################################################################################################################################

if (length(MLest)==0) {

Data=data.frame(Quantile=Data.full$Quantile,Death_Time=Data.full$Time,Temp=Data.full$Temp,Clim,Cont,Pop,Geno=Data.full$Geno)
Data=Data[!is.na(Data[,1]),]
Data$IM_CT=ifelse(Data[,1]=="CT","CT","IM")
Data$TP=ifelse(grepl("C",Data$Clim,fixed=T),"C","H")
Data$PR=ifelse(grepl("D",Data$Clim,fixed=T),"D","W")

summary(lm(Death_Time~TP,data=Data[Data$Temp=="30C",]))

library(plotrix)
dev.new()
multhist(list(Data[Data$IM_CT=="CT"&Data$Temp=="20C",2],
              Data[Data$IM_CT!="CT"&Data$Temp=="20C",2],
              Data[Data$IM_CT=="CT"&Data$Temp=="30C",2],
              Data[Data$IM_CT!="CT"&Data$Temp=="30C",2]),
         freq=F)

model_0=lm(Death_Time~Temp*IM_CT,data=Data)
model_1=lmer(Death_Time~Temp*IM_CT+(1|Pop),data=Data)
model_2=lmer(Death_Time~Temp*IM_CT+(0+Temp+IM_CT|Pop),data=Data)
model_3=lmer(Death_Time~Temp*IM_CT+(0+Temp+IM_CT|Pop)+(0+Temp+IM_CT|Geno:Pop),data=Data)
model_4=lmer(Death_Time~Temp*IM_CT+(0+Temp*IM_CT|Pop)+(0+Temp*IM_CT|Geno:Pop),data=Data)

library(broom)
model.names <- c("1 Null", "2 Population", "3 Population + Genotype", "4 Population + Genotype w inter")
summ.table <- do.call(rbind, lapply(list(model_0, model_2, model_3, model_4), function(x) broom::glance(x)[c("sigma","logLik","AIC","BIC","deviance","df.residual")]))
table.cols <- c("df.residual", "deviance", "AIC")
reported.table <- summ.table[table.cols]
names(reported.table) <- c("Resid. Df", "Resid. Dev", "AIC")
reported.table[['dAIC']] <-  with(reported.table, AIC - min(AIC))
reported.table[['weight']] <- with(reported.table, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
reported.table$AIC <- NULL
reported.table$weight <- round(reported.table$weight, 2)
reported.table$dAIC <- round(reported.table$dAIC, 1)
row.names(reported.table) <- model.names

reported.table

Y_hat=getME(model_4,'X') %*% fixef(model_4) + getME(model_4,'Z') %*% unlist(ranef(model_4))
Y_obs=getME(model_4,'y')
Z=getME(model_4,'Z')
X=getME(model_4,'X')

for (i in c("_1","_2","_3","_4")) {
	VCV=as.data.frame(VarCorr(get(paste("model",i,sep=""))))
	VCV=VCV[is.na(VCV[,3]),]
	VARPOP=sum(VCV[VCV$grp=="Pop",4])/sum(VCV[1:nrow(VCV),4])
	VARGENO=sum(VCV[grepl("Geno",VCV$grp),4])/sum(VCV[1:nrow(VCV),4])
	VAREXP=sum(VCV[1:nrow(VCV)-1,4])/sum(VCV[1:nrow(VCV),4])
	print(i)
	print(paste(VARPOP,VARGENO,VAREXP))
}

VCV=as.data.frame(VarCorr(model_2))
VCV=VCV[is.na(VCV[,3]),]
VAREXP=sum(VCV[1:nrow(VCV)-1,4])/sum(VCV[1:nrow(VCV),4])

PopFX=as.data.frame(ranef(model_3)[[2]])
GenoFX=as.data.frame(ranef(model_3)[[1]])
Plot.dict = matrix(c("NLA","HW","USA",23,"#0ad15d",
                     "SLA","HW","USA",23,"#0ad15d",
                     "NTX","HD","USA",23,"#d10a2e",
                     "STX","HD","USA",23,"#d10a2e",
                     "NB","CW","USA",23,"#0a7bd2",
                     "NS","CW","USA",23,"#0a7bd2",
                     "NCA","CD","USA",23,"#c98806",
                     "SCA","CD","USA",23,"#c98806",
                     "NQ1","HW","AUS",21,"#0ad15d",
                     "NQ2","HW","AUS",21,"#0ad15d",
                     "SQ1","HD","AUS",21,"#d10a2e",
                     "SQ2","HD","AUS",21,"#d10a2e",
                     "NTS","CW","AUS",21,"#0a7bd2",
                     "STS","CW","AUS",21,"#0a7bd2",
                     "SA","CD","AUS",21,"#c98806",
                     "VIC","CD","AUS",21,"#c98806"),byrow=T,ncol=5)

png(file="SuppMat_PredvsObs_Resid.png",width=680,height=480)
par(mfrow=c(1,2))
plot('n',xlim=c(min(Y_hat),max(Y_hat)),ylim=c(min(Y_obs),max(Y_obs)),main="Goodness-of-fit of full model for longevity",ylab="Observed longevity (in hours)",xlab="Predicted longevity (in hours)")
points(Y_hat[X[,"Temp30C"]==1&X[,"IM_CTIM"]==1],Y_obs[X[,"Temp30C"]==1&X[,"IM_CTIM"]==1],pch=1,cex=.75)
points(Y_hat[X[,"Temp30C"]==1&X[,"IM_CTIM"]==0],Y_obs[X[,"Temp30C"]==1&X[,"IM_CTIM"]==0],pch=4,cex=.75)
points(Y_hat[X[,"Temp30C"]==0&X[,"IM_CTIM"]==1],Y_obs[X[,"Temp30C"]==0&X[,"IM_CTIM"]==1],pch=1,cex=.75,col=grey(.5))
points(Y_hat[X[,"Temp30C"]==0&X[,"IM_CTIM"]==0],Y_obs[X[,"Temp30C"]==0&X[,"IM_CTIM"]==0],pch=4,cex=.75,col=grey(.5))
abline(0,1,col=2)
legend(par("usr")[1]+.0*(par("usr")[2]-par("usr")[1]),
	   par("usr")[3]+.97*(par("usr")[4]-par("usr")[3]),
	   col=c(grey(.5),1),
	   legend=c(expression(20*degree*C),expression(30*degree*C)),
	   pch=15,bty='n',pt.cex=1,x.intersp=.5)
legend(par("usr")[1]+.15*(par("usr")[2]-par("usr")[1]),
	   par("usr")[3]+.97*(par("usr")[4]-par("usr")[3]),
	   legend=c("Control","Imidacloprid"),
	   pch=c(4,1),bty='n',pt.cex=.7,x.intersp=.4)
hist(residuals(model_4),main="Distribution of residuals for full model",breaks=20,xlim=c(-250,250))
dev.off()

png(file="Fig1c_DistMort.png",width=600,height=800)
plot(density(Data[Data$IM_CT!="CT"&Data$Temp=="20C",2],bw=10),lwd=2,ylim=c(0,0.02),col=4,xlim=c(min(Data$Death_Time,na.rm=T),max(Data$Death_Time,na.rm=T)),cex.lab=1.5,cex.main=2,main="Density distribution of time to death",ylab="Frequency",xlab="Time to death in hours")
lines(density(Data[Data$IM_CT=="CT"&Data$Temp=="20C",2],bw=10),lwd=2,col=4,lty=3)
lines(density(Data[Data$IM_CT!="CT"&Data$Temp=="30C",2],bw=10),lwd=2,ylim=c(0,0.02),col=2)
lines(density(Data[Data$IM_CT=="CT"&Data$Temp=="30C",2],bw=10),lwd=2,col=2,lty=3)
legend("topright",
       col=c("blue","blue","red","red"),
       lty=c(1,3,1,3),
       pch='',bty='n',lwd=2,
       expression(paste(20*degree*C," with Imidacloprid"),paste(20*degree*C," no Imidacloprid"),paste(30*degree*C," with Imidacloprid"),paste(30*degree*C," no Imidacloprid")))
dev.off()

TRADE=data.frame(Pop=rownames(PopFX),
                 IM=as.numeric(PopFX[,3]),
                 Temp30C=as.numeric(PopFX[,2]),
                 Cont=c("USA","USA","USA","AUS","AUS","USA","AUS","USA","AUS","USA","USA","AUS","AUS","AUS","USA","AUS"),
                 ANC=grepl("LA|NQ1|SCA|VIC|STX",rownames(PopFX)))
TRADE=TRADE[order(TRADE$Pop),]

summary(lm(Temp30C/IM~0+ANC,data=TRADE))

png("Fig1d-e_GFX.png",width=400,height=800)
par(mfrow=c(2,1))
par(mar=c(3,3,3,1),mgp=c(1.7,.4,0))
plot('n',xlim=c(min(PopFX[,2]),max(PopFX[,2])),ylim=c(min(PopFX[,3]),max(PopFX[,3])),ylab="Response to imidacloprid (hours)",xlab=expression(paste("Response to ",30*degree*C," (hours)")),cex.lab=1.5,cex.main=2,main="Population effects")
abline(0,1,lty=1,col=grey(.5))
abline(0,-1,lty=2,col=grey(.5))
abline(h=0,lty=3,col=grey(.7))
abline(v=0,lty=3,col=grey(.7))
for (i in Plot.dict[,1]) points(PopFX[grep(i,rownames(PopFX)),3]~PopFX[grep(i,rownames(PopFX)),2],pch=as.numeric(Plot.dict[grep(i,Plot.dict[,1]),4]),bg=Plot.dict[grep(i,Plot.dict[,1]),5],cex=2)
legend(par("usr")[1]+.57*(par("usr")[2]-par("usr")[1]),
	   par("usr")[3]+(par("usr")[4]-par("usr")[3]),
	   legend=c("USA","AUS"),
	   pch=c(23,21),bty='n',pt.cex=2,x.intersp=.8)
legend(par("usr")[1]+.7*(par("usr")[2]-par("usr")[1]),
	   par("usr")[3]+(par("usr")[4]-par("usr")[3]),
	   col=c("#0ad15d","#d10a2e","#c98806","#0a7bd2"),
	   legend=c("Hot and Wet","Hot and Dry","Cold and Dry","Cold and Wet"),pch=15,bty='n',pt.cex=2,x.intersp=.8)
par(mar=c(3,3,3,1),mgp=c(1.7,.4,0))
plot('n',xlim=c(min(GenoFX[,2]),max(GenoFX[,2])),ylim=c(min(GenoFX[,3]),max(GenoFX[,3])),ylab="Response to imidacloprid (hours)",xlab=expression(paste("Response to ",30*degree*C," (hours)")),cex.lab=1.5,cex.main=2,main="Line within population effects")
for (i in Plot.dict[,1]) points(GenoFX[grep(i,rownames(GenoFX)),3]~GenoFX[grep(i,rownames(GenoFX)),2],pch=as.numeric(Plot.dict[grep(i,Plot.dict[,1]),4]),bg=Plot.dict[grep(i,Plot.dict[,1]),5],cex=1.5)
abline(0,1,lty=1,col=grey(.5))
abline(0,-1,lty=2,col=grey(.5))
abline(h=0,lty=3,col=grey(.7))
abline(v=0,lty=3,col=grey(.7))
dev.off()

#########################################################################################################################################
#                                         Climate Origin effect on Polymorphism
#########################################################################################################################################

Data=read.csv("Sequencing_Summary/mapped_reads_summary.csv",header=T)
for (i in unique(Data[,1])) print(paste(i,sum(Data[Data[,1]==i,4])))

Data=read.csv("Sequencing_Summary/SeqPoly_descr.csv",header=T)
t.test(Data$Pi.Gono,Data$Theta.Gono,paired=T)
summary(aov(In.2L.t ~ Cont + Temp * Prec, data = Data))
summary(aov(In.3R.Payne ~ Cont + Temp * Prec, data = Data))
summary(aov(In.3R.Mo ~ Cont + Temp * Prec, data = Data))

for (i in names(Data)[-c(1:4)]) {
	print(i)
	print(anova(lm(get(i)~Cont+Temp*Prec,data=Data)))
}

for (i in colnames(PopFX)[-1]) {
	for (j in names(Data)[-c(1:4)]) {
		print(paste("Effect of",j,"on",i))
		model=paste("PopFX$",i,"~Data$",j,sep="")
		#model=paste("Data$",j,"~Data$AvCov",sep="") # an alternative model is needed to test the effect of the average coverage on polymorphism, requires to comment line 325
		print(anova(lm(model)))
	}
}

q()
}


#########################################################################################################################################
#
#                                         Generates the pheno.py files to be used for GWAlpha.py
# - This section of the script only focuses on flies exposed to insecticide and NOT the control that have not been sequenced
# - should be run after line 1-104 commenting out lines 46 64-66 and 72
#
#########################################################################################################################################

library(FactoMineR)

for (pop in unique(Data.full$Pop)) {
	for (temp in c("20C","30C")) {
		Data.sub=Data.full[(Data.full$Pop==pop)&(Data.full$Temp==temp),]
		MLest.sub=MLest[grep(paste(pop,temp,sep="_"),names(MLest))]
		par=apply(matrix(unlist(MLest.sub),nrow=length(MLest.sub)),1,mean)
		print("__________________________")
		print(paste(pop,temp,sep="-"))
		print(paste("sig=",sd(Data.sub$Time), "- min=",min(Data.sub$Time), "- max=",max(Data.sub$Time)))
		quant=quantile(Data.sub$Time,c(0.064, 0.256, 0.744, 0.936))
		print(quant)
		write.infile(paste("Pheno_name=",paste("\"",pop,"-",temp,"\"",sep="")),file=paste(pop,"-",temp,"_pheno.py",sep=""),append=F)
		write.infile(paste("sig=",sd(Data.sub$Time)),file=paste(pop,"-",temp,"_pheno.py",sep=""),append=T)
		write.infile(paste("MIN=",min(Data.sub$Time)),file=paste(pop,"-",temp,"_pheno.py",sep=""),append=T)
		write.infile(paste("MAX=",max(Data.sub$Time)),file=paste(pop,"-",temp,"_pheno.py",sep=""),append=T)
		write.infile(paste("perc=[ 0.064, 0.256, 0.744, 0.936]"),file=paste(pop,"-",temp,"_pheno.py",sep=""),append=T)
		write.infile(paste("q=[",paste(quant,collapse=","),"]"),file=paste(pop,"-",temp,"_pheno.py",sep=""),append=T)
	}
}

