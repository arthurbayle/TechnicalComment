# Technical comment of "From white to green: Snow cover loss and increased vegetation productivity in the European Alps#
# *correspond author : arthur.bayle.env@gmail.com

# Load package
library(terra)
library(DescTools)
library(mblm)
library(Kendall)
library(zoo)

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
# ----- FIGURE 1 A ----- #
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
# ----- FIGURE 1 B ----- #
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
DF.ALL = data.frame(NUMB = DF.N2, NDVI = DF.NDVI2, YEAR = rep(1984:2021, each=100000))
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
DF.ALL = data.frame(NUMB = DF.N2, PERM = DF.PERM2, YEAR = rep(1984:2021, each=100000))
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

# --------------------------#
# ----- SELECT PIXELS ----- #
# --------------------------#

setwd("F:/Recherche/Publications/2022/Science/")
LIST = list.files("F:/Recherche/Publications/2022/Science/DATArumpf/SCIENCE-20220614T050323Z-001/SCIENCE/",full.names=T,pattern="NDVI")
NAMES = list.files("F:/Recherche/Publications/2022/Science/DATArumpf/SCIENCE-20220614T050323Z-001/SCIENCE/",pattern="NDVI")
NAMES = StrRight(NAMES,29)

# SELECT IMAGE TO REPROJECT ALL IMAGES ON
GRID = rast(LIST[503])

GRID[which(GRID[] == 0)] = NA
SAMPLE = 10000
n = sample(x = which(!is.na(GRID[])),size = SAMPLE)

# -------------------------------------------#
# ----- COMPUTE TRENDS WITH ALL IMAGES ----- #
# -------------------------------------------#

# LIST IMAGES
LIST.NDVI = list.files("F:/Recherche/Publications/2022/Science/DATArumpf/SCIENCE-20220614T050323Z-001/SCIENCE/",full.names=T,pattern="NDVI")
NAMES = list.files("F:/Recherche/Publications/2022/Science/DATArumpf/SCIENCE-20220614T050323Z-001/SCIENCE/",pattern="NDVI")
NAMES = StrRight(NAMES,29)
PR = substr(NAMES,6,11)
LIST.NDVI = LIST.NDVI[which(PR == "195028")]

LISTYEAR = as.numeric(substr(StrRight(LIST.NDVI,17),1,4))

YEARS = 1984:2021

DF.SLP.NDVI.BIASED = data.frame(ID = 1:SAMPLE)
DF.PVAL.NDVI.BIASED = data.frame(ID = 1:SAMPLE)

DF.NUMBTOT.BIASED = data.frame(YEAR = 1984:2021)
z=2

DF.NUMB.BIASED = data.frame(ID = 1:SAMPLE)
DF.NDVI.BIASED = data.frame(ID = 1:SAMPLE)

for(i in 1:length(YEARS)){
  print(paste0(YEARS[i]))
  year = YEARS[i]
  
  # SELECT PER YEAR
  LIST.NDVIcur = LIST.NDVI[which(LISTYEAR == year)]
  
  ############################################## NDVI 0.75 QUANTILE ----
  
  # LOAD IMAGES AND KEEP ONLY SAMPLES
  DF.Nndvi = data.frame(ID = 1:SAMPLE)
  STACK = sds(LIST.NDVIcur)
  for(v in 1:length(nlyr(STACK))){
    CUR = STACK[v]
    CUR[which(CUR[] == 0)] = NA
    DF.Nndvi[,v] = CUR[n][,1]
    colnames(DF.Nndvi)[v] = paste0("NDVI")
  }
  
  # COMPUTE 0.75 NDVI QUANTILE
  QUANT = apply(DF.Nndvi,1,function(x){quantile(x,0.75,na.rm=T)})
  DF.NDVI.BIASED[,i+1] = QUANT
  colnames(DF.NDVI.BIASED)[i+1] = as.character(year)
  
  # COMPUTE NUMBER OF VALUES
  for(k in 1:ncol(DF.Nndvi)){DF.Nndvi[which(!is.na(DF.Nndvi[,k])),k] = 1}
  DF.NUMB.BIASED[,i+1] = apply(DF.Nndvi,1,function(x){sum(x,na.rm=T)})
  colnames(DF.NUMB.BIASED)[i+1] = as.character(year)
  
}

for(k in 1:ncol(DF.NUMB.BIASED)){DF.NUMB.BIASED[which(DF.NUMB.BIASED[,k]==0),k] = NA}
DF.NUMBTOT.BIASED[,z] = as.numeric(apply(DF.NUMB.BIASED[,-1],2,function(x){mean(x,na.rm=T)}))

# ----- COMPUTE TRENDS FOR THE THREE VARIABLES ----- #

for(i in 1:nrow(DF.NDVI.BIASED)){
  print(paste0(i))
  
  dat.ndvi = DF.NDVI.BIASED[i,-1]
  dat.ndvi = data.frame(YEAR = as.numeric(colnames(dat.ndvi)), NDVI = as.numeric(dat.ndvi))
  dat.ndvi$NDVI[which(dat.ndvi$NDVI == 0)] = NA
  dat.ndvi = as.data.frame(na.omit(dat.ndvi))
  
  if(length(dat.ndvi$YEAR) > 12){
    
    MOD.NDVI = mblm(NDVI ~ YEAR,dataframe = dat.ndvi)
    DF.SLP.NDVI.BIASED[i,z] = MOD.NDVI$coefficients[[2]]
    dat.ts = as.ts(read.zoo(dat.ndvi))
    DF.PVAL.NDVI.BIASED[i,z] = MannKendall(dat.ts)[[2]][1]
    
  } else {
    DF.SLP.NDVI.BIASED[i,z] = NA; DF.PVAL.NDVI.BIASED[i,z] = NA}
}


# ----------------------------------------------#
# ----- COMPUTE TRENDS WITH RANDOMIZATION ----- #
# ----------------------------------------------#


# LIST IMAGES
LIST.NDVI = list.files("F:/Recherche/Publications/2022/Science/DATArumpf/SCIENCE-20220614T050323Z-001/SCIENCE/",full.names=T,pattern="NDVI")
NAMES = list.files("F:/Recherche/Publications/2022/Science/DATArumpf/SCIENCE-20220614T050323Z-001/SCIENCE/",pattern="NDVI")
NAMES = StrRight(NAMES,29)
PR = substr(NAMES,6,11)
LIST.NDVI = LIST.NDVI[which(PR == "195028")]

LISTYEAR = as.numeric(substr(StrRight(LIST.NDVI,17),1,4))

YEARS = 1984:2021

DF.SLP.NDVI = data.frame(ID = 1:SAMPLE)
DF.PVAL.NDVI = data.frame(ID = 1:SAMPLE)

DF.NUMBTOT = data.frame(YEAR = 1984:2021)
DF.NUMBTOTsd = data.frame(YEAR = 1984:2021)
DF.NDVITOT = data.frame(YEAR = 1984:2021)

DF.NUMB = data.frame(ID = 1:SAMPLE)
DF.NDVI = data.frame(ID = 1:SAMPLE)

for(z in 2:101){
  
  # ----- COMPUTE ANNUAL METRICS BY RANDOMLY SELECTING 3 IMAGES PER YEAR ----- #
  
  DF.NUMB = data.frame(ID = 1:SAMPLE)
  DF.NDVI = data.frame(ID = 1:SAMPLE)
  
  for(i in 1:length(YEARS)){
    print(paste0(z," ----- ",YEARS[i]))
    year = YEARS[i]
    
    # SELECT PER YEAR
    LIST.NDVIcur = LIST.NDVI[which(LISTYEAR == year)]
    
    # RANDOMLY SELECT 3 IMAGES
    if(length(LIST.NDVIcur) < 4){RANDOM = 1:length(LIST.NDVIcur)}else{
      RANDOM = sample(1:length(LIST.NDVIcur),3)}
    LIST.NDVIcurRANDOM = LIST.NDVIcur[RANDOM]
    
    DATE = StrRight(LIST.NDVIcurRANDOM,17)
    DATE = as.Date(paste(substr(DATE,1,4),substr(DATE,5,6),substr(DATE,7,8),sep="-"))
    
    ############################################## NDVI 0.75 QUANTILE ----
    
    # LOAD IMAGES AND KEEP ONLY SAMPLES
    DF.N = data.frame(ID = 1:SAMPLE)
    if(length(LIST.NDVIcurRANDOM) > 1){
      
      STACK = sds(LIST.NDVIcurRANDOM)
      for(v in 1:length(nlyr(STACK))){
        CUR = STACK[v]
        CUR[which(CUR[] == 0)] = NA
        DF.N[,v] = CUR[n][,1]
        colnames(DF.N)[v] = paste0("NDVI")
      }
      
    } else {
      STACK = rast(LIST.NDVIcurRANDOM)
      CUR = STACK
      CUR[which(CUR[] == 0)] = NA
      DF.N[,i] = CUR[n][,1]
      colnames(DF.N)[i] = paste0("NDVI")
    }
    
    # COMPUTE 0.75 NDVI QUANTILE
    QUANT = apply(DF.N,1,function(x){quantile(x,0.75,na.rm=T)})
    DF.NDVI[i+1] = QUANT
    colnames(DF.NDVI)[i+1] = as.character(year)
    
    # COMPUTE NUMBER OF VALUES
    for(k in 1:ncol(DF.N)){DF.N[which(!is.na(DF.N[,k])),k] = 1}
    DF.NUMB[i+1] = apply(DF.N,1,function(x){sum(x,na.rm=T)})
    colnames(DF.NUMB)[i+1] = as.character(year)
    
  }
  
  for(k in 1:ncol(DF.NUMB)){DF.NUMB[which(DF.NUMB[,k]==0),k] = NA}
  DF.NUMBTOT[,z] = apply(DF.NUMB[,-1],2,function(x){mean(x,na.rm=T)})
  DF.NUMBTOTsd[,z] = apply(DF.NUMB[,-1],2,function(x){sd(x,na.rm=T)})
  
  DF.NDVITOT[,z] = apply(DF.NDVI[,-1],2,function(x){mean(x,na.rm=T)})
  
  
  # ----- COMPUTE TRENDS FOR THE THREE VARIABLES ----- #
  
  for(i in 1:nrow(DF.NDVI)){
    print(paste0(i))
    
    dat.ndvi = DF.NDVI[i,-1]
    dat.ndvi = data.frame(YEAR = as.numeric(colnames(dat.ndvi)), NDVI = as.numeric(dat.ndvi))
    dat.ndvi$NDVI[which(dat.ndvi$NDVI == 0)] = NA
    dat.ndvi = as.data.frame(na.omit(dat.ndvi))
    
    if(length(dat.ndvi$YEAR) > 12){
      
      MOD.NDVI = mblm(NDVI ~ YEAR,dataframe = dat.ndvi)
      DF.SLP.NDVI[i,z] = MOD.NDVI$coefficients[[2]]
      dat.ts = as.ts(read.zoo(dat.ndvi))
      DF.PVAL.NDVI[i,z] = MannKendall(dat.ts)[[2]][1]
      
    } else {
      DF.SLP.NDVI[i,z] = NA; DF.PVAL.NDVI[i,z] = NA}
  }
}



# -----------------------#
# ----- FIGURE 1 F ----- #
# -----------------------#

png("FIGURES/FIGURE1F.png",width=1500,height=1400,res = 300)

plot(DF.NUMBTOT.BIASED,ylim=c(0,8),type="l",lwd=2,lty=1,xaxs="i",xlab="Time",ylab="",col="black")
abline(lm(DF.NUMBTOT.BIASED$V2~DF.NUMBTOT.BIASED$YEAR),col="black",lwd=2)
lines(DF.NUMBTOT,lwd=2,lty=1,col="royalblue")
abline(lm(DF.NUMBTOT$V2~DF.NUMBTOT$YEAR),col="royalblue",lwd=2)
polygon(x=c(DF.NUMBTOT$V2-DF.NUMBTOTsd$V2,rev(DF.NUMBTOT$V2+DF.NUMBTOTsd$V2)),
        y=c(DF.NUMBTOT$YEAR,rev(DF.NUMBTOT$YEAR)),col=adjustcolor("red",alpha.f = 0.2),border = NA)
legend("topleft",legend=c("Using all images available (ALL)","Degrading to 3 images per \nyear maximum (3MAX)"),
       col=c("black","royalblue"),lty=c(2,1),bty="n")

dev.off()


# -----------------------#
# ----- FIGURE 1 G ----- #
# -----------------------#

PERC.PVAL.BIASED = na.omit(DF.PVAL.NDVI.BIASED)
PERC.SLP.BIASED = na.omit(DF.SLP.NDVI.BIASED)
PERC.BIASED.POS0.01 = (length(which(PERC.PVAL.BIASED$V2 < 0.01 & PERC.SLP.BIASED$V2 > 0))/length(PERC.PVAL.BIASED$V2))*100
PERC.BIASED.POS0.05 = (length(which(PERC.PVAL.BIASED$V2 > 0.01 & PERC.PVAL.BIASED$V2 < 0.05 & PERC.SLP.BIASED$V2 > 0))/length(PERC.PVAL.BIASED$V2))*100
PERC.BIASED.NEG0.01 = (length(which(PERC.PVAL.BIASED$V2 < 0.01 & PERC.SLP.BIASED$V2 < 0))/length(PERC.PVAL.BIASED$V2))*100
PERC.BIASED.NEG0.05 = (length(which(PERC.PVAL.BIASED$V2 > 0.01 & PERC.PVAL.BIASED$V2 < 0.05 & PERC.SLP.BIASED$V2 < 0))/length(PERC.PVAL.BIASED$V2))*100

PERC.POS0.01 = c()
PERC.POS0.05 = c()
PERC.NEG0.01 = c()
PERC.NEG0.05 = c()

for(i in 1:10){
  PERC.PVAL.NDVI = na.omit(DF.PVAL.NDVI)
  PERC.SLP.NDVI = na.omit(DF.SLP.NDVI)
  PERC.POS0.01[i] = (length(which(PERC.PVAL.NDVI[,i+1] < 0.01 & PERC.SLP.NDVI[,i+1] > 0))/length(PERC.PVAL.NDVI[,i+1]))*100
  PERC.POS0.05[i] = (length(which(PERC.PVAL.NDVI[,i+1] > 0.01 & PERC.PVAL.NDVI[,i+1] < 0.05 & PERC.SLP.NDVI[,i+1] > 0))/length(PERC.PVAL.NDVI[,i+1]))*100
  PERC.NEG0.01[i] = (length(which(PERC.PVAL.NDVI[,i+1] < 0.01 & PERC.SLP.NDVI[,i+1] < 0))/length(PERC.PVAL.NDVI[,i+1]))*100
  PERC.NEG0.05[i] = (length(which(PERC.PVAL.NDVI[,i+1] > 0.01 & PERC.PVAL.NDVI[,i+1] < 0.05 & PERC.SLP.NDVI[,i+1] < 0))/length(PERC.PVAL.NDVI[,i+1]))*100
  
}

barplot(c(PERC.BIASED,mean(PERC.NDVI)),ylim=c(0,100))
sd(PERC.NDVI)
MAT = matrix(c(PERC.BIASED.POS0.01 , PERC.BIASED.POS0.05,
               mean(PERC.POS0.01) , mean(PERC.POS0.05)),ncol=2,nrow=2)
colnames(MAT) = c("ALL","3MAX")
rownames(MAT) = c("pos001","pos005")


png("FIGURES/FIGURE1G.png",width=1500,height=1400,res = 300)

bp = barplot(MAT,ylim=c(0,100),col=c(adjustcolor("darkgreen",alpha.f = 0.5),adjustcolor("palegreen",alpha.f = 0.5)),ylab="% of pixels",las=1,
             width=4,xlim=c(0,10))
text(bp[1],40,paste0(round(PERC.BIASED.POS0.01,2),"%"),cex=2)
text(bp[1],30,"(P-val < 0.01)",cex=0.7)

text(bp[2],30,paste0(round(mean(PERC.POS0.01),2),"%"),cex=2)
text(bp[2],20,"(P-val < 0.01)",cex=0.7)

box(which="plot", bty="]")
dev.off()














