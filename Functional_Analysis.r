#!/usr/bin/Rscript
options(warn=-1)
#########################################################################################################################################
#
#                  This R script allows to replicate the trait analysis proposed in Fourneir-Level el al. (2018, submitted)
#
# - It requires the lme4 and ggplot2 packages to be installed to perform the mixed-linear modelling
# - The input data can be sourced at the following URL:  
# - The input data folder should be located in the R working directory without altering the folder names or structure
#
#########################################################################################################################################

library(lme4)
library(ggplot2)

setwd("Functional_tests/")

for (Temp in c("20C","30C")) {

dev.new()

File.list=list.files(pattern=paste("MU_",Temp,"_.\\.csv",sep=""))

Data_all=as.data.frame(matrix(NA,ncol=8,nrow=0))

for (file in File.list) {

Data=read.csv(file)

Start=Data$hour[1]+Data$Minute[1]/60
Data=Data[-1,]
Data$TTD=24*Data$Day+Data$hour+Data$Minute/60-Start
IM=Data$Treat=="IM"

unique(Data$Geno)

CTmod=lmer(TTD~1+(1|Geno),data=Data[!IM,])
CTmort=fixef(CTmod)+ranef(CTmod)[[1]][,1]
names(CTmort)=rownames(ranef(CTmod)[[1]])

Data$TTD_CT=Data$TTD
for (g in names(CTmort)) Data$TTD_CT[Data$Geno==g]=Data$TTD[Data$Geno==g]/CTmort[g]
Data=Data[Data$Geno %in% names(CTmort),]
IM=Data$Treat=="IM"

Data_all=rbind(Data_all,Data)
}

names(Data_all)=names(Data)
Data=Data_all
IM=Data$Treat=="IM"


alpha<-0.05 # for a (1.00-alpha) confidence interval
Data.summary=data.frame(Genotype=levels(Data[IM,]$Geno),
                        mean=tapply(Data[IM,]$TTD_CT, Data[IM,]$Geno, mean),
                        n=tapply(Data[IM,]$TTD_CT, Data[IM,]$Geno, length),
                        sd=tapply(Data[IM,]$TTD_CT, Data[IM,]$Geno, sd)
                        )
Data.summary$sem=Data.summary$sd/sqrt(Data.summary$n)
Data.summary$me=qt(1-alpha/2, df=Data.summary$n)*Data.summary$sem

p=ggplot(Data.summary, aes(x = Genotype, y = mean))+ylab("")
p=p+geom_bar(position = position_dodge(), stat="identity")
p=p+geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem, width=.5))
p=p+ggtitle(paste("Longevity exposed to imidacloprid relative to control at ",gsub("C","",Temp),"°C",sep=""))
p=p+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)


print("Prm")
#dev.new()
png(paste("../Fig4a-b_Prm_",Temp,".png",sep=""),width=450,height=450)

PRMs=c("CT","KO","KO","KO","KO","KO","KO","KO","KO","CT","KO","CT","KO")
names(PRMs)=c("yw",
              "ywx11752",
              "ywx11752W",
              "11752xyw",
              "11752xywW",
              "ywx56623",
              "ywx56623W",
              "56623xyw",
              "56623xywW",
              "W1118GDxTub",
              "33615GDxTub",
              "TubxW1118GD",
              "Tubx33615GD")
#Selecting only MiMIC and P-ele
PRMs=c("CT","KO","KO")#,"KO","KO")
names(PRMs)=c("yw",
              #"ywx11752W",
              #"11752xywW",
              "ywx56623W",
              "56623xywW")

Data.sum.PRMs=Data.summary[rownames(Data.summary) %in% names(PRMs),]
Data.sum.PRMs=Data.sum.PRMs[order(match(Data.sum.PRMs$Genotype,names(PRMs))),]
Data.sum.PRMs$Status=PRMs[match(Data.sum.PRMs$Genotype,names(PRMs))]

Data.PRMs=Data[IM & Data$Geno %in% names(PRMs),]
Data.PRMs$Status=PRMs[match(Data.PRMs$Geno,names(PRMs))]
print(summary(lm(TTD_CT~Status,data=Data.PRMs)))

par(mar=c(9.5,3,2,1))
BAR=barplot(Data.sum.PRMs$mean,ylim=c(0,max(c(.75,1.1*max(Data.sum.PRMs$mean+Data.sum.PRMs$sem,na.rm=T)))),col=ifelse(Data.sum.PRMs$Status=="CT","#4bf2c2","#d198ea"),las=2)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=grey(.9))
abline(h=seq(0.2,2,.2),col="white",lty=2)
BAR=barplot(Data.sum.PRMs$mean,ylim=c(0,max(c(.75,1.1*max(Data.sum.PRMs$mean+Data.sum.PRMs$sem,na.rm=T)))),col=ifelse(Data.sum.PRMs$Status=="CT","#4bf2c2","#d198ea"),las=2,add=T)
segments(BAR, Data.sum.PRMs$mean - Data.sum.PRMs$sem * 2, BAR,Data.sum.PRMs$mean + Data.sum.PRMs$sem * 2, lwd = 1.5)
arrows(BAR, Data.sum.PRMs$mean - Data.sum.PRMs$sem * 2, BAR,Data.sum.PRMs$mean + Data.sum.PRMs$sem * 2, lwd = 1.5, angle = 90,code = 3, length = 0.1)
SamSize=table(Data$Geno)[names(table(Data$Geno)) %in% names(PRMs)]
SamSize=SamSize[order(names(SamSize))]
names(SamSize)=c("yw[PRM]","yw[prm/+]","yw[+/prm]")
labels=c(expression(paste("yw; +/+")),expression(paste("yw; prm/+")),expression(paste("yw; +/prm")))
axis(1,at=BAR,labels = NA, padj = 1,las=2)
axis(1,at=BAR-.2,labels = labels, padj = 1,las=2,tick=F)
axis(1,at=BAR,labels = paste("\n(n=",SamSize,")",sep=""),mgp=c(3, 2, 0), padj = 1,las=2,tick=F)


mtext("Time-to-death relative to control",cex=1.5)

u=par("usr")
v=c(grconvertX(u[1:2], "user", "ndc"),grconvertY(u[3:4], "user", "ndc"))
v=c( (v[1]+v[2])*1/2, v[2], (v[3]+v[4])*3/5, v[4] )
par(fig=v, new=TRUE, mar=c(0,0,1,0.5) ,mgp=c(.5,.5,0),las=2)
boxplot(TTD_CT~Status,data=Data.PRMs,main=paste("Effect of Prm mutations at ",gsub("C","",Temp),"°C",sep=""),tck=-.1,cex.main=.8,cex.axis=.6,col.main=ifelse(Temp=="20C","blue","red"))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="white")
boxplot(TTD_CT~Status,data=Data.PRMs,main=paste("Effect of Prm mutations at ",gsub("C","",Temp),"°C",sep=""),tck=-.1,cex.main=.8,cex.axis=.6,col.main=ifelse(Temp=="20C","blue","red"),add=T)
box(col=ifelse(Temp=="20C","blue","red"))
dev.off()


print("Dalpha3s")
#dev.new()
png(paste("../Fig4c-d_Dalpha3_",Temp,".png",sep=""),width=450,height=450)

Dalphas=c("CT","KO","KO","KO","KO")
names(Dalphas)=c("CAS9","Dalpha3_1.7","Dalpha3_5.3","Dalpha3_11.4","Dalpha3_15.1")

Data.sum.Dalphas=Data.summary[rownames(Data.summary) %in% names(Dalphas),]
Data.sum.Dalphas=Data.sum.Dalphas[order(match(Data.sum.Dalphas$Genotype,names(Dalphas))),]
Data.sum.Dalphas$Status=Dalphas[match(Data.sum.Dalphas$Genotype,names(Dalphas))]

Data.Dalphas=Data[IM & Data$Geno %in% names(Dalphas),]
Data.Dalphas$Status=Dalphas[match(Data.Dalphas$Geno,names(Dalphas))]
print(summary(lm(TTD_CT~Status,data=Data.Dalphas)))

par(mar=c(9.5,3,2,1))
BAR=barplot(Data.sum.Dalphas$mean,ylim=c(0,max(c(.75,1.1*max(Data.sum.Dalphas$mean+Data.sum.Dalphas$sem,na.rm=T)))),col=ifelse(Data.sum.Dalphas$Status=="CT","#4bf2c2","#d198ea"),las=2)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=grey(.9))
abline(h=seq(0.2,2,.2),col="white",lty=2)
BAR=barplot(Data.sum.Dalphas$mean,ylim=c(0,max(c(.75,1.1*max(Data.sum.Dalphas$mean+Data.sum.Dalphas$sem,na.rm=T)))),col=ifelse(Data.sum.Dalphas$Status=="CT","#4bf2c2","#d198ea"),las=2,add=T)
segments(BAR, Data.sum.Dalphas$mean - Data.sum.Dalphas$sem * 2, BAR,Data.sum.Dalphas$mean + Data.sum.Dalphas$sem * 2, lwd = 1.5)
arrows(BAR, Data.sum.Dalphas$mean - Data.sum.Dalphas$sem * 2, BAR,Data.sum.Dalphas$mean + Data.sum.Dalphas$sem * 2, lwd = 1.5, angle = 90,code = 3, length = 0.1)
SamSize=table(Data$Geno)[names(table(Data$Geno)) %in% names(Dalphas)]
SamSize=SamSize[order(names(SamSize))]

labels=c("yw;{CAS9};+",expression(paste("yw;{CAS9};nachrα3"^"1.7")),expression(paste("yw;{CAS9};nachrα3"^"5.3")),expression(paste("yw;{CAS9};nachrα3"^"11.4")),expression(paste("yw{CAS9};nachrα3"^"15.1")))
axis(1,at=BAR,labels = NA, padj = 1,las=2)
axis(1,at=BAR-.2,labels = labels, padj = 1,las=2,tick=F)
axis(1,at=BAR,labels = paste("\n(n=",SamSize,")",sep=""),mgp=c(3, 2, 0), padj = 1,las=2,tick=F)
mtext("Time-to-death relative to control",cex=1.5)

u=par("usr")
v=c(grconvertX(u[1:2], "user", "ndc"),grconvertY(u[3:4], "user", "ndc"))
v=c( (v[1]+v[2])*1/2, v[2], (v[3]+v[4])*3/5, v[4] )
par(fig=v, new=TRUE, mar=c(0,0,1,0.5) ,mgp=c(.5,.5,0),las=2)
boxplot(TTD_CT~Status,data=Data.Dalphas,main=paste("Effect of nACHRα3 mutations at ",gsub("C","",Temp),"°C",sep=""),tck=-.1,cex.main=.8,cex.axis=.6,col.main=ifelse(Temp=="20C","blue","red"))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="white")
boxplot(TTD_CT~Status,data=Data.Dalphas,main=paste("Effect of nACHRα3 mutations at ",gsub("C","",Temp),"°C",sep=""),tck=-.1,cex.main=.8,cex.axis=.6,col.main=ifelse(Temp=="20C","blue","red"),add=T)
box(col=ifelse(Temp=="20C","blue","red"))
dev.off()

print("Pde9s")
#dev.new()
png(paste("../FigS3c-d_Pde9_",Temp,".png",sep=""),width=450,height=450)

Pde9s=c("CT","KO","CT","KO","CT","KO","CT","KO","CT","KO","CT","KO")
names(Pde9s)=c("W1118GDxTub",
               "40468GDxTub",
               "TubxW1118GD",
               "Tubx40468GD",
               "KKxTub",
               "106959KKxTub",
               "TubxKK",
               "Tubx106959KK",
               "545x851bar",
               "545x851W",
               "851x545bar",
               "851x545W")

Pde9s=c("CT","KO","CT","KO")
names(Pde9s)=c("545x851bar",
               "545x851W",
               "851x545bar",
               "851x545W")

Data.sum.Pde9s=Data.summary[rownames(Data.summary) %in% names(Pde9s),]
Data.sum.Pde9s=Data.sum.Pde9s[order(match(Data.sum.Pde9s$Genotype,names(Pde9s))),]
Data.sum.Pde9s$Status=Pde9s[match(Data.sum.Pde9s$Genotype,names(Pde9s))]

Data.Pde9s=Data[IM & Data$Geno %in% names(Pde9s),]
Data.Pde9s$Status=Pde9s[match(Data.Pde9s$Geno,names(Pde9s))]
print(summary(lm(TTD_CT~Status,data=Data.Pde9s)))

par(mar=c(10.5,3,2,1))
BAR=barplot(Data.sum.Pde9s$mean,ylim=c(0,max(c(.75,1.1*max(Data.sum.Pde9s$mean+Data.sum.Pde9s$sem,na.rm=T)))),col=ifelse(Data.sum.Pde9s$Status=="CT","#4bf2c2","#d198ea"),las=2)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=grey(.9))
abline(h=seq(0.2,2,.2),col="white",lty=2)
BAR=barplot(Data.sum.Pde9s$mean,ylim=c(0,max(c(.75,1.1*max(Data.sum.Pde9s$mean+Data.sum.Pde9s$sem,na.rm=T)))),col=ifelse(Data.sum.Pde9s$Status=="CT","#4bf2c2","#d198ea"),las=2,add=T)
segments(BAR, Data.sum.Pde9s$mean - Data.sum.Pde9s$sem * 2, BAR,Data.sum.Pde9s$mean + Data.sum.Pde9s$sem * 2, lwd = 1.5)
arrows(BAR, Data.sum.Pde9s$mean - Data.sum.Pde9s$sem * 2, BAR,Data.sum.Pde9s$mean + Data.sum.Pde9s$sem * 2, lwd = 1.5, angle = 90,code = 3, length = 0.1)
SamSize=table(Data$Geno)[names(table(Data$Geno)) %in% names(Pde9s)]
SamSize=SamSize[order(names(SamSize))]

labels=c(expression(paste("w;+/ΔBSC545")),expression(paste("w;ΔBSC851/ΔBSC545")),expression(paste("w;ΔBSC851/+")),expression(paste("w;ΔBSC545/ΔBSC851")))
axis(1,at=BAR,labels = NA, padj = 1,las=2)
axis(1,at=BAR-.2,labels = labels, padj = 1,las=2,tick=F)
axis(1,at=BAR,labels = paste("\n(n=",SamSize,")",sep=""),mgp=c(3, 2, 0), padj = 1,las=2,tick=F)
mtext("Time-to-death relative to control",cex=1.5)

u=par("usr")
v=c(grconvertX(u[1:2], "user", "ndc"),grconvertY(u[3:4], "user", "ndc"))
v=c( (v[1]+v[2])*1/2, v[2], (v[3]+v[4])*3/5, v[4] )
par(fig=v, new=TRUE, mar=c(0,0,1,0.5) ,mgp=c(.5,.5,0),las=2)
boxplot(TTD_CT~Status,data=Data.Pde9s,main="",tck=-.1,cex.main=.8,cex.axis=.6,col.main=ifelse(Temp=="20C","blue","red"))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="white")
boxplot(TTD_CT~Status,data=Data.Pde9s,main=paste("Effect of Pde9 mutations at ",gsub("C","",Temp),"°C",sep=""),tck=-.1,cex.main=.8,cex.axis=.6,col.main=ifelse(Temp=="20C","blue","red"),add=T)
box(col=ifelse(Temp=="20C","blue","red"))
dev.off()


print("UGT86Ds")
png(paste("../FigS3e-f_UGT86Ds_",Temp,".png",sep=""),width=450,height=450)


UGT86D=c("CT","KO","CT","KO")
names(UGT86D)=c("17306X5506sb",
                "17306X5506W",
                "5506x17306sb",
                "5506x17306W")

Data.sum.UGT86D=Data.summary[rownames(Data.summary) %in% names(UGT86D),]
Data.sum.UGT86D=Data.sum.UGT86D[order(match(Data.sum.UGT86D$Genotype,names(UGT86D))),]
Data.sum.UGT86D$Status=UGT86D[match(Data.sum.UGT86D$Genotype,names(UGT86D))]

Data.UGT86D=Data[IM & Data$Geno %in% names(UGT86D),]
Data.UGT86D$Status=UGT86D[match(Data.UGT86D$Geno,names(UGT86D))]
print(summary(lm(TTD_CT~Status,data=Data.UGT86D)))

par(mar=c(10.5,3,2,1))
BAR=barplot(Data.sum.UGT86D$mean,ylim=c(0,max(c(.75,1.1*max(Data.sum.UGT86D$mean+Data.sum.UGT86D$sem,na.rm=T)))),col=ifelse(Data.sum.UGT86D$Status=="CT","#4bf2c2","#d198ea"),las=2)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=grey(.9))
abline(h=seq(0.2,2,.2),col="white",lty=2)
BAR=barplot(Data.sum.UGT86D$mean,ylim=c(0,max(c(.75,1.1*max(Data.sum.UGT86D$mean+Data.sum.UGT86D$sem,na.rm=T)))),col=ifelse(Data.sum.UGT86D$Status=="CT","#4bf2c2","#d198ea"),las=2,add=T)
segments(BAR, Data.sum.UGT86D$mean - Data.sum.UGT86D$sem * 2, BAR,Data.sum.UGT86D$mean + Data.sum.UGT86D$sem * 2, lwd = 1.5)
arrows(BAR, Data.sum.UGT86D$mean - Data.sum.UGT86D$sem * 2, BAR,Data.sum.UGT86D$mean + Data.sum.UGT86D$sem * 2, lwd = 1.5, angle = 90,code = 3, length = 0.1)
SamSize=table(Data$Geno)[names(table(Data$Geno)) %in% names(UGT86D)]
SamSize=SamSize[order(names(SamSize))]
labels=c(expression(paste("w;+/ΔED5506")),expression(paste("w;ΔExel17306/ΔED5506")),expression(paste("w;+/ΔExel17306")),expression(paste("w;ΔED5506/ΔExel17306")))
axis(1,at=BAR,labels = NA, padj = 1,las=2)
axis(1,at=BAR-.2,labels = labels, padj = 1,las=2,tick=F)
axis(1,at=BAR,labels = paste("\n(n=",SamSize,")",sep=""),mgp=c(3, 2, 0), padj = 1,las=2,tick=F)
mtext("Time-to-death relative to control",cex=1.5)

u=par("usr")
v=c(grconvertX(u[1:2], "user", "ndc"),grconvertY(u[3:4], "user", "ndc"))
v=c( (v[1]+v[2])*1/2, v[2], (v[3]+v[4])*3/5, v[4] )
par(fig=v, new=TRUE, mar=c(0,0,1,0.5) ,mgp=c(.5,.5,0),las=2)
boxplot(TTD_CT~Status,data=Data.UGT86D,main="",tck=-.1,cex.main=.75,cex.axis=.6,col.main=ifelse(Temp=="20C","blue","red"))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="white")
boxplot(TTD_CT~Status,data=Data.UGT86D,main=paste("Effect of UGT86Ds mutations at ",gsub("C","",Temp),"°C",sep=""),tck=-.1,cex.main=.75,cex.axis=.6,col.main=ifelse(Temp=="20C","blue","red"),add=T)
box(col=ifelse(Temp=="20C","blue","red"))
dev.off()


print("Anks")
png(paste("../FigS3a-b_Ank_",Temp,".png",sep=""),width=450,height=450)
#dev.new()

Ank=c("CT","KO","CT","CT","KO","KO")
names(Ank)=c("yw",
               "Ank")#,
               #"KKxTub",
               #"TubxKK",
               #"107369KKxTub",
               #"Tubx107369KK")

#Ank=c("CT","CT","KO","KO")
#names(Ank)=c( "KKxTub",
#               "TubxKK",
#               "107369KKxTub",
#               "Tubx107369KK")

Data.sum.Ank=Data.summary[rownames(Data.summary) %in% names(Ank),]
Data.sum.Ank=Data.sum.Ank[order(match(Data.sum.Ank$Genotype,names(Ank))),]
Data.sum.Ank$Status=Ank[match(Data.sum.Ank$Genotype,names(Ank))]

Data.Ank=Data[IM & Data$Geno %in% names(Ank),]
Data.Ank$Status=Ank[match(Data.Ank$Geno,names(Ank))]
print(summary(lm(TTD_CT~Status,data=Data.Ank)))

par(mar=c(10.5,3,2,1))
BAR=barplot(Data.sum.Ank$mean,ylim=c(0,max(c(.75,1.1*max(Data.sum.Ank$mean+Data.sum.Ank$sem,na.rm=T)))),col=ifelse(Data.sum.Ank$Status=="CT","#4bf2c2","#d198ea"),las=2)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=grey(.9))
abline(h=seq(0.2,2,.2),col="white",lty=2)
BAR=barplot(Data.sum.Ank$mean,ylim=c(0,max(c(.75,1.1*max(Data.sum.Ank$mean+Data.sum.Ank$sem,na.rm=T)))),col=ifelse(Data.sum.Ank$Status=="CT","#4bf2c2","#d198ea"),las=2,add=T)
segments(BAR, Data.sum.Ank$mean - Data.sum.Ank$sem * 2, BAR,Data.sum.Ank$mean + Data.sum.Ank$sem * 2, lwd = 1.5)
arrows(BAR, Data.sum.Ank$mean - Data.sum.Ank$sem * 2, BAR,Data.sum.Ank$mean + Data.sum.Ank$sem * 2, lwd = 1.5, angle = 90,code = 3, length = 0.1)
SamSize=table(Data$Geno)[names(table(Data$Geno)) %in% names(Ank)]
SamSize=SamSize[order(names(SamSize))]
labels=c(expression(paste("yw;+")),expression(paste("yw;ank")))
axis(1,at=BAR,labels = NA, padj = 1,las=2)
axis(1,at=BAR-.2,labels = labels, padj = 1,las=2,tick=F)
axis(1,at=BAR,labels = paste("\n(n=",SamSize,")",sep=""),mgp=c(3, 2, 0), padj = 1,las=2,tick=F)
mtext("Time-to-death relative to control",cex=1.5)

u=par("usr")
v=c(grconvertX(u[1:2], "user", "ndc"),grconvertY(u[3:4], "user", "ndc"))
v=c( (v[1]+v[2])*1/2, v[2], (v[3]+v[4])*3/5, v[4] )
par(fig=v, new=TRUE, mar=c(0,0,1,0.5) ,mgp=c(.5,.5,0),las=2)
boxplot(TTD_CT~Status,data=Data.Ank,main="",tck=-.1,cex.main=.8,cex.axis=.6,col.main=ifelse(Temp=="20C","blue","red"))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="white")
boxplot(TTD_CT~Status,data=Data.Ank,main=paste("Effect of Anks mutations at ",gsub("C","",Temp),"°C",sep=""),tck=-.1,cex.main=.8,cex.axis=.6,col.main=ifelse(Temp=="20C","blue","red"),add=T)
box(col=ifelse(Temp=="20C","blue","red"))
dev.off()

}

#Rogues
png(paste("../Fig4k_Rogues_",Temp,".png",sep=""),width=450,height=450)

Misc=c("CT","KO","KO","KO","KO")
names(Misc)=c( "yw",
               "Or65a",
               "Sona",
               "AlphaEST-7",
               "DGRP354")

Data.sum.Misc=Data.summary[rownames(Data.summary) %in% names(Misc),]
Data.sum.Misc=Data.sum.Misc[order(match(Data.sum.Misc$Genotype,names(Misc))),]
Data.sum.Misc$Status=Misc[match(Data.sum.Misc$Genotype,names(Misc))]

Data.Misc=Data[IM & Data$Geno %in% names(Misc),]
Data.Misc$Status=Misc[match(Data.Misc$Geno,names(Misc))]
print(summary(lm(TTD_CT~Status,data=Data.Misc)))

par(mar=c(10,3,2,1))
BAR=barplot(Data.sum.Misc$mean,ylim=c(0,max(c(.75,1.1*max(Data.sum.Misc$mean+Data.sum.Misc$sem,na.rm=T)))),col=ifelse(Data.sum.Misc$Status=="CT","#4bf2c2","#d198ea"),las=2)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=grey(.9))
abline(h=seq(0.2,2,.2),col="white",lty=2)
BAR=barplot(Data.sum.Misc$mean,ylim=c(0,max(c(.75,1.1*max(Data.sum.Misc$mean+Data.sum.Misc$sem,na.rm=T)))),col=ifelse(Data.sum.Misc$Status=="CT","#4bf2c2","#d198ea"),las=2,add=T)
segments(BAR, Data.sum.Misc$mean - Data.sum.Misc$sem * 2, BAR,Data.sum.Misc$mean + Data.sum.Misc$sem * 2, lwd = 1.5)
arrows(BAR, Data.sum.Misc$mean - Data.sum.Misc$sem * 2, BAR,Data.sum.Misc$mean + Data.sum.Misc$sem * 2, lwd = 1.5, angle = 90,code = 3, length = 0.1)
SamSize=table(Data$Geno)[names(table(Data$Geno)) %in% names(Misc)]
SamSize=SamSize[order(names(SamSize))]
axis(1,at=BAR,labels = paste("\n(n=",SamSize,")",sep=""), padj = 1,las=2)
mtext("Time-to-death relative to control",cex=1.5)

u=par("usr")
v=c(grconvertX(u[1:2], "user", "ndc"),grconvertY(u[3:4], "user", "ndc"))
v=c( (v[1]+v[2])*1/2, v[2], (v[3]+v[4])*3/5, v[4] )
par(fig=v, new=TRUE, mar=c(0,0,1,0.5) ,mgp=c(.5,.5,0),las=2)
boxplot(TTD_CT~Status,data=Data.Misc,main=paste("Effect of Miscs mutations at ",gsub("C","",Temp),"°C",sep=""),tck=-.1,cex.main=.6,cex.axis=.6,col.main=ifelse(Temp=="20C","blue","red"))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="white")
boxplot(TTD_CT~Status,data=Data.Misc,main=paste("Effect of Miscs mutations at ",gsub("C","",Temp),"°C",sep=""),tck=-.1,cex.main=.6,cex.axis=.6,col.main=ifelse(Temp=="20C","blue","red"),add=T)
box(col=ifelse(Temp=="20C","blue","red"))
dev.off()
