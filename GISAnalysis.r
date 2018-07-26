#!usr/bin/Rscript
options(warn=-1)
#########################################################################################################################################
#
#                  This R script allows to replicate the GIS analysis proposed in Fourneir-Level el al. (2018, submitted)
#
# - It requires the raster package to be installed to perform the mixed-linear modelling
# - The input data can be sourced at the following URL:  
# - The input data folder should be located in the R working directory without altering the folder names or structure
#
#########################################################################################################################################


library(raster)

########################################################################################
#
#                                IMIGlobal: sampling zones
#
########################################################################################

#Loading a raster layer
tmax <- raster(paste(getwd(), "/GIS/Bioclim/bio10.bil", sep = ""))
pmax <- raster(paste(getwd(), "/GIS/Bioclim/bio16.bil", sep = ""))
pmin <- raster(paste(getwd(), "/GIS/Bioclim/bio17.bil", sep = ""))
#tmin <- tmin/10  # Worldclim temperature data come in decimal degrees 
tmax=tmax/10  # look at the info
alt <- getData('worldclim', var='alt',res=2.5)

USA <- c(-125, -60, 9, 61)
tmax_USA <- crop(tmax, USA)
#tmin_USA <- crop(tmin, USA)
#pmax_USA <- crop(pmax, USA)
pmin_USA <- crop(pmin, USA)
alt_USA = raster::crop(alt,USA)
Zones_USA=tmax_USA

Zones_USA[tmax_USA!=5553.7]=0.1
#Hot and Wet
Zones_USA[tmax_USA>26&pmin_USA>250]=1.1
#Cold and Wet
Zones_USA[tmax_USA<18&pmin_USA>250]=2.1
#Cold and Dry
Zones_USA[tmax_USA<18&pmin_USA<100]=3.1
#Hot and Dry
Zones_USA[tmax_USA>26&pmin_USA<100]=4.1
#And remove high altitude
Zones_USA[alt_USA>700]=0.1

Zones_USA[tmax_USA==5553.7]=-0.9

breakpoints <- c(-1,0,1,2,3,4,5)
colors <- c(grey(.2),grey(.6),"#0ad15d","#0a7bd2","#c98806","#d10a2e")

png("Fig1a_SamplingScheme.png",width=420,height=800)
par(oma = c(5, 1, 1, 0))

par(mfrow=c(2,1),mar=c(2,3,2,0))

plot(Zones_USA,breaks=breakpoints,col=colors,legend=F)

Loc_USA=matrix(c(	"Baton Rouge",	30.43,-91.12,
                	"West Monroe",	32.5,-92.14,
                 	"Lubbock",		33.57,-100.9,
                 	"Stonewall",	30.23,-99.66,
                	"Anderson Valley Nth",	39.15,-123.65,
                	"Anderson Valley Sth",39,-123.21,
                	"St John",		45.25,-66.07,
                	"Gaspereau",	45.08,-64.37
                	),ncol=3,byrow=T)
points(as.numeric(Loc_USA[,2])~as.numeric(Loc_USA[,3]),col="white",pch=19)
mtext("North America", outer=T, line=-1.7,cex=1.5)

AUS <- c(110, 155, -45, -5)
tmax_AUS <- crop(tmax, AUS)
#tmin_AUS <- crop(tmin, AUS)
pmax_AUS <- crop(pmax, AUS)
pmin_AUS <- crop(pmin, AUS)
alt_AUS = crop(alt,AUS)

Zones_AUS=tmax_AUS
Zones_AUS[tmax_AUS!=5553.7]=0.1

#Hot and Wet
Zones_AUS[tmax_AUS>26&pmin_AUS>250]=1.1
#Cold and Wet
Zones_AUS[tmax_AUS<18&pmax_AUS>300]=2.1
#Cold and Dry
Zones_AUS[tmax_AUS<18&pmin_AUS<100]=3.1
#Hot and Dry
Zones_AUS[tmax_AUS>26&pmin_AUS<100]=4.1

Zones_AUS[alt_AUS>700]=0.1

Zones_AUS[tmax_AUS==5553.7]=-0.9

breakpoints <- c(-1,0,1,2,3,4,5)
colors <- c(grey(.2),grey(.6),"#0ad15d","#0a7bd2","#c98806","#d10a2e")
plot(Zones_AUS,breaks=breakpoints,col=colors,legend=F)
mtext("Australia", outer=T, line=-26.5,cex=1.5)

arrows(147,-17,146.2,-17.48,col="#0ad15d",length=0.05,angle=20)
arrows(147.1,-18,146.3,-17.617,col="#0ad15d",length=0.05,angle=20)

Loc_AUS=matrix(c(	"Innisfail1",	-17.48,146,
                	"Innisfail2",	-17.51,145.98,
                 	"Isla Gorge",	-25.17,149.93,
                 	"Yeppoon",		-23.13,150.73,
                	"Mt Gambier",	-37.82,140.77,
                	"Hamilton",-37.43,142.02,
                	"Lucaston",-41.23,146.98,
                	"Hillwood",-43,146.9
                	),ncol=3,byrow=T)
points(as.numeric(Loc_AUS[,2])~as.numeric(Loc_AUS[,3]),col="white",pch=19)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("bottom",
       col=c("#0ad15d","#d10a2e","#c98806","#0a7bd2"),
       legend=c("Hot (>26° in warmest quarter) and Wet (>250mm in driest quarter)",
                "Hot (>26° in warmest quarter) and Dry (<100mm in driest quarter)",
                "Cold (<18° in warmest quarter) and Dry (<100mm in driest quarter)",
                "Cold (<18° in warmest quarter) and Wet (>250mm  in driest quarter)"),
       inset=c(0,0),xpd = TRUE,ncol=1,pch=15,bty='n',pt.cex=1.5,x.intersp=.8)

dev.off()

########################################################################################
#
#                Pesticide usage ( et al (2010) Nature)
#
########################################################################################

Pest=raster(paste(getwd(), "/GIS/Pesticide_Loading.asc", sep = ""))
png(file="Pesticide_usage.png")
plot(Pest)
points(as.numeric(Loc_USA[,2])~as.numeric(Loc_USA[,3]),col="purple",pch=19)
points(as.numeric(Loc_AUS[,2])~as.numeric(Loc_AUS[,3]),col="purple",pch=19)

Pest_USA=extract(Pest,matrix(as.numeric(Loc_USA[,c(3,2)]),ncol=2))
Pest_AUS=extract(Pest,matrix(as.numeric(Loc_AUS[,c(3,2)]),ncol=2))
Pest_use=rbind(cbind(Loc_USA,Pest_USA),cbind(Loc_AUS,Pest_AUS))
rownames(Pest_use)=c("SLA","NLA","NTX","STX","NCA","SCA","NB","NS","NQ1","NQ2","SQ1","SQ2","SA","VIC","NTS","STS")
Pest_use=Pest_use[order(rownames(Pest_use)),]
print("Pesticide use at locatio of origin of populations")
colnames(Pest_use)=c("Population","lat","Long","Pesticide concentration ")
print(Pest_use)
