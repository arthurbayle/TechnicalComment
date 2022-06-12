# Technical comment of "From white to green: Snow cover loss and increased vegetation productivity in the European Alps#
# *correspond author : arthur.bayle.env@gmail.com

# Load package
library(terra)

# Set working directory
setwd("F:/Recherche/Publications/2022/Science/")

# -----------------------------------#
# ---------- PREPARE DATA ---------- #
# -----------------------------------#

# Proj : +proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs
# ETRS89-extended / LAEA Europe (EPSG:3035)
ETRS89 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"

# Prepare standard grid by reprojecting to ETRS89 and aggregate to 300 m
# to limit processing time.
GRID = rast("DATArumpf/RAW/PROB/alps_Prob_2021.tif")
GRID = project(GRID,ETRS89)
GRID = aggregate(GRID,fact=13.3333333)
GRID[] = 1
writeRaster(GRID,"DATArumpf/GRID.tif")

GRID = rast("DATArumpf/GRID.tif")

########## COUNT
LIST = list.files("F:/Recherche/Publications/2022/Science/DATArumpf/RAW/COUNT/",full.names=T)

YEAR = 1984:2021
for(i in 1:length(LIST)){
  print(YEAR[i])
  NDVI = rast(LIST[i])
  NDVI = project(NDVI,GRID,method="near")
  writeRaster(NDVI,paste0("DATArumpf/PROJ/COUNT/alps_imageCount_",YEAR[i],".tif"),overwrite=T)
}

######## NDVI
LIST = list.files("F:/Recherche/Publications/2022/Science/DATArumpf/RAW/NDVI/",full.names=T)

YEAR = 1984:2021
for(i in 1:length(LIST)){
  print(YEAR[i])
  NDVI = rast(LIST[i])
  NDVI = project(NDVI,GRID)
  writeRaster(NDVI,paste0("DATArumpf/PROJ/NDVI/alps_NDVI_",YEAR[i],".tif"))
}

######### PROB
LIST = list.files("F:/Recherche/Publications/2022/Science/DATArumpf/RAW/PROB/",full.names=T)

YEAR = 1984:2021
for(i in 1:length(LIST)){
  print(YEAR[i])
  NDVI = rast(LIST[i])
  NDVI = project(NDVI,GRID)
  writeRaster(NDVI,paste0("DATArumpf/PROJ/PROB/alps_Prob_",YEAR[i],".tif"))
}

######### PERM
LIST = list.files("F:/Recherche/Publications/2022/Science/DATArumpf/RAW/PERM/",full.names=T)

YEAR = 1984:2021
for(i in 1:length(LIST)){
  print(YEAR[i])
  NDVI = rast(LIST[i])
  NDVI = project(NDVI,GRID,method = "near")
  writeRaster(NDVI,paste0("DATArumpf/PROJ/PERM/alps_Perm_",YEAR[i],".tif"),overwrite=T)
}

# -------------------------------#
# ---------- FIGURE 1 ---------- #
# -------------------------------#
GRID = rast("DATArumpf/GRID.tif")

# Proj : +proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs
# ETRS89-extended / LAEA Europe (EPSG:3035)
ETRS89 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"

# DEM 
DEM = rast("F:/Recherche/Travail/01_Choro/BIS/DATA/DEM/DEM_EU-ALPS.tif")
DEM = project(DEM,GRID)

# -----------------------#
# ----- FIGURE 1 B ----- #
# -----------------------#

# We randomly selected 100.000 pixels around the Alps to conduct the analysis
n = sample(which(!is.na(GRID[])),100000)

# Load 'Image Count'
LIST.IC = list.files("DATArumpf/PROJ/COUNT/",full.names=T)
NUMB = sds(LIST.IC)
DF.N = data.frame(ID = 1:100000)
YEAR = 1984:2021
for(i in 1:length(nlyr(NUMB))){
  print(YEAR[i])
  CUR = NUMB[i]
  CUR[which(CUR[] == 0)] = NA
  DF.N[,i] = CUR[n][,1]
  colnames(DF.N)[i] = paste0("NUMB_",YEAR[i])
}

# Plot
png("FIGURES/FIGURE1A.png",width=1500,height=2500,res = 300)
# 1984 to 1996
X = apply(DF.N[,1:13],MARGIN = 1,FUN=function(x){mean(x,na.rm=T)})
plot(y=DEM[n][,1], x=X,pch=".",cex=1.2,col=adjustcolor("red",alpha.f = 0.7),xlim=c(1,10),xaxs="i",xlab="Mean number of images per year",
     ylab="Altitude (m)",ylim=c(1600,3200),yaxs="i")
abline(h=seq(1000,4000,500),v=seq(0,10,2),lty=2,col="grey")
LOESS = loess.smooth(X,DEM[n][,1], span = 2/3, degree = 2)
lines(LOESS$x,LOESS$y,col="red",lwd=4,lty=3)

# 1997 to 2009
X = apply(DF.N[,14:26],MARGIN = 1,FUN=function(x){mean(x,na.rm=T)})
points(y=DEM[n][,1], x=X,pch=".",cex=1.2,col=adjustcolor("darkgreen",alpha.f = 0.7),xlim=c(1,8),xaxs="i")
LOESS = loess.smooth(X,DEM[n][,1], span = 2/3, degree = 2)
lines(LOESS$x,LOESS$y,col="darkgreen",lwd=4,lty=3)

# 2010 to 2021
X = apply(DF.N[,27:38],MARGIN = 1,FUN=function(x){mean(x,na.rm=T)})
points(y=DEM[n][,1], x=X,pch=".",cex=1.2,col=adjustcolor("royalblue",alpha.f = 0.7),xlim=c(1,8),xaxs="i")
LOESS = loess.smooth(X,DEM[n][,1], span = 2/3, degree = 2)
lines(LOESS$x,LOESS$y,col="royalblue",lwd=4,lty=3)

legend("topright",legend=c("1984 - 1996","1997 - 2009","2010 - 2021"),col=c("red","darkgreen","blue"),lwd=4,lty=3,bty="n",cex=1.5)
abline(h=1700,lwd=2)
text(3.5,1740,"Rumpf et al. (1) study limit")
polygon(x=c(1,10,10,1),y=c(1600,1600,1700,1700),col=adjustcolor("grey",alpha.f = 0.7))
dev.off()

# -----------------------#
# ----- FIGURE 1 A ----- #
# -----------------------#

png("FIGURES/FIGURE1B.png",width=1500,height=1400,res = 300)

Y = apply(DF.N,MARGIN=2,function(x){mean(x,na.rm=T)})
Y.se = apply(DF.N,MARGIN=2,function(x){sd(x,na.rm=T)})
X = 1984:2021

plot(X,Y,type="l",lwd=4,ylim=c(0,12),yaxs="i",xaxs="i",
     xlab="Time",ylab="",lty=3)
polygon(x=c(X,rev(X)),y=c(Y-2*Y.se,rev(Y+2*Y.se)),border = NA,col = "grey90")
polygon(x=c(X,rev(X)),y=c(Y-Y.se,rev(Y+Y.se)),border = NA,col = "grey")
lines(X,Y,lwd=4,lty=3)

arrows(x0=1984,x1=2011,y0=8,y1=8,col = "red",lwd=2,angle = 90,length = 0.1)
arrows(x0=2011,x1=1984,y0=8,y1=8,col = "red",lwd=2,angle = 90,length = 0.1)
text(x=1992,y=9,"L5 TM",col="red",cex=1.5)

arrows(x0=2021,x1=1999,y0=9,y1=9,col = "black",lwd=2,angle = 90,length = 0.1)
text(x=2006,y=10,"L7 ETM+",col="black",cex=1.5)

arrows(x0=2021,x1=2013,y0=10,y1=10,col = "orange",lwd=2,angle = 90,length = 0.1)
text(x=2017,y=11,"L8 OLI",col="orange",cex=1.5)
box()
dev.off()

# ------------------------#
# ----- FIGURE BIAS ----- #
# ------------------------#


# SELECT 100.000 SAMPLES FROM VEGETATION PRODUCTIVITY
NDVIL =  list.files("F:/Recherche/Publications/2022/Science/DATArumpf/PROJ/",recursive = T,full.names=T,pattern="NDVI")
YEAR = 1984:2021

NDVI = sds(NDVIL)
DF.NDVI = data.frame(ID = 1:100000)
YEAR = 1984:2021
for(i in 1:length(nlyr(NUMB))){
  print(YEAR[i])
  DF.NDVI[,i] = NDVI[i][n][,1]
  colnames(DF.NDVI)[i] = paste0("NDVI_",YEAR[i])
}

# SELECT 100.000 SAMPLES FROM SUMMER SNOW
SSL =  list.files("F:/Recherche/Publications/2022/Science/DATArumpf/PROJ/",recursive = T,full.names=T,pattern="Prob")
YEAR = 1984:2021

SS = sds(SSL)
DF.SS = data.frame(ID = 1:100000)
YEAR = 1984:2021
for(i in 1:length(nlyr(NUMB))){
  print(YEAR[i])
  DF.SS[,i] = SS[i][n][,1]
  colnames(DF.SS)[i] = paste0("SS_",YEAR[i])
}

# SELECT 100.000 SAMPLES FROM PERMANENT SNOW
PERML =  list.files("F:/Recherche/Publications/2022/Science/DATArumpf/PROJ/",recursive = T,full.names=T,pattern="Perm")
YEAR = 1984:2021

PERM = sds(PERML)
DF.PERM = data.frame(ID = 1:100000)
YEAR = 1984:2021
for(i in 1:length(nlyr(NUMB))){
  print(YEAR[i])
  DF.PERM[,i] = PERM[i][n][,1]
  colnames(DF.PERM)[i] = paste0("PERM_",YEAR[i])
}


# ------------------------#
# ----- FIGURE 1 C ------ #
# ------------------------#

png("FIGURES/FIGURE1c_ndvi.png",width=1500,height=1400,res = 300)

DF.N2 = round(data.frame(newcol = c(t(DF.N)), stringsAsFactors=FALSE))
DF.NDVI2 = data.frame(newcol = c(t(DF.NDVI)), stringsAsFactors=FALSE)
DF.ALL = data.frame(NUMB = DF.N2, NDVI = DF.NDVI2, YEAR = rep(1984:2021, each=10000))
colnames(DF.ALL) = c("NUMB","NDVI","YEAR")
DF.ALL[which(DF.ALL$NUMB == 0),]=NA
DF.ALL[which(DF.ALL$NUMB > 10),]=NA

boxplot(DF.ALL$NDVI ~ DF.ALL$NUMB,xaxs="i",col=adjustcolor("darkgreen",alpha.f = 0.5),
        xlab="Mean number of images per year",outline=F,ylim=c(0,1),
        ylab="Vegetation productivity")
LOESS = loess.smooth(DF.ALL$NUMB,DF.ALL$NDVI, span = 1, degree = 2)
lines(LOESS$x,LOESS$y, col="red",lwd=4,lty=1)

dev.off()

# ------------------------#
# ----- FIGURE 1 D ------ #
# ------------------------#

png("FIGURES/FIGURE1d_prob.png",width=1500,height=1400,res = 300)

DF.N2 = round(data.frame(newcol = c(t(DF.N)), stringsAsFactors=FALSE))
DF.SS2 = data.frame(newcol = c(t(DF.SS)), stringsAsFactors=FALSE)
DF.ALL = data.frame(NUMB = DF.N2, SS = DF.SS2)
colnames(DF.ALL) = c("NUMB","SS")
DF.ALL[which(DF.ALL$NUMB == 0),]=NA
DF.ALL[which(DF.ALL$NUMB > 10),]=NA

# Remove pixels that are almost never snowy
DF.ALL$SS[which(DF.ALL$SS < 0.1)] = NA

boxplot(DF.ALL$SS ~ DF.ALL$NUMB,xaxs="i",col=adjustcolor("royalblue",alpha.f = 0.5),
        xlab="Mean number of images per year",outline=F,ylim=c(0,1),
        ylab="Summer snow proportion")
LOESS = loess.smooth(DF.ALL$NUMB,DF.ALL$SS, span = 1, degree = 2)
lines(LOESS$x,LOESS$y, col="red",lwd=4)

dev.off()

# ------------------------#
# ----- FIGURE 1 E ------ #
# ------------------------#

png("FIGURE_PERMobs.png",width=1500,height=1400,res = 300)

# Prepare data
DF.N2 = round(data.frame(newcol = c(t(DF.N)), stringsAsFactors=FALSE))
DF.PERM2 = data.frame(newcol = c(t(DF.PERM)), stringsAsFactors=FALSE)
DF.ALL = data.frame(NUMB = DF.N2, PERM = DF.PERM2, YEAR = rep(1984:2021, each=10000))
colnames(DF.ALL) = c("NUMB","PERM","YEAR")

# Range between 1 to 10 images
DF.ALL[which(DF.ALL$NUMB == 0),]=NA
DF.ALL[which(DF.ALL$NUMB > 10),]=NA

# Compute percentage of permanent snow cover
AGG = aggregate(DF.ALL$PERM,by=list(DF.ALL$NUMB),function(x){table(x)})
AGG = data.frame(NUMB = AGG$Group.1,NOSNOW=AGG$x[,1],SNOW=AGG$x[,2])
AGG$PERC = (AGG$SNOW/(AGG$NOSNOW+AGG$SNOW))*100

BP=barplot(AGG$PERC ~ AGG$NUMB,col=adjustcolor("darkblue",alpha.f = 1),
           xlab="Mean number of images per year",ylim=c(0,20),
           ylab="Proportion of permanent snow",xaxt="n")
axis(1,at = BP[,1],labels = 1:10)
box()
LOESS = loess.smooth(AGG$NUMB,AGG$PERC, span = 1, degree = 2)
lines(LOESS$x,LOESS$y, col="red",lwd=4)

dev.off()