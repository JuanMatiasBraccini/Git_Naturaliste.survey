####   SCRIPT FOR ANALYSING NATURALISTE SURVEY DATA ####

# INDEX :
# -- DATA SECTION --
# -- CONTROL SECTION --
# -- PROCEDURE SECTION --
#1. CATCH RATES FROM FISHERY INDEPENDENT SURVEYS   
#2. Multivariate analysis and Ecosystem indicators


# -- REPORT SECTION --


#NOTES:   shot is defined by SHEET_NO
#         Periods: prior to 2001, drumlines; 2001 onwards, longlines
#         fixed stations: 2001:2008, 2015: onwards
#         Station Definition within 10 nm
#         Timeline: 
#             First period: targeted at sandbar sharks for tagging
#             Then WAMSI objectivies
#             Then FRDC tagging of dusky sharks
#             Then standard survey
#         Zero inflated and Hurdle model don't produce SE for predictions,
#         hence, did Monte Carlo or bootstrapping



# DATA SECTION ------------------------------------------------------------

#Define user
User="Matias"

#Sharks data base 
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R")

library(lubridate)  #for dates manipulation
library(geosphere)
library(ggplot2)
library(vcd)    #correlation between factors
library(pscl)  #zero-inflated models
library(MCMCglmm) #zero-inflated models
library(glmmADMB) #zero-inflated models
library(tweedie)
library(statmod) # Provides  tweedie  family functions
library(bbmle)  #AIC table
library(lmtest) #likelihood ratio tests
library(mgcv)
library(zigam)
library(PBSmapping)
data(worldLLhigh)
library(plotrix) #multihistogram
library(gplots)  #annotated boxplot
#install.packages("C:/Matias/R/COZIGAM_2.0-3.tar.gz",type="source")
#library(COZIGAM)
library(mvtnorm)      #for multivariate normal pdf
library(caret)
library(ReporteRs)
library(gridExtra)
#library(cede)
library(tidyr)
library(dplyr)
library(mpath)
library(glmulti)  #model selection
library(emmeans)
library(foreach)
library(doParallel)
library(geosphere)
library(imputeTS)


#install.packages("countreg", repos = "http://R-Forge.R-project.org")
#see great vignette: https://cran.r-project.org/web/packages/pscl/vignettes/countreg.pdf
library("countreg")



source("C:/Matias/R/HighstatLibV6.R")  #for corvif ()
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Plot.Map.R")
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Compare.error.structure.R")  
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/fn.fig.R")
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Smart_par.R")



#Get info for tiger sharks tissue samples
Tiger.gen.samples=read.csv("C:/Matias/Data/Shark_bio/Tiger.gen.samples.csv")
Tiger.gen=subset(DATA,BAG_NO%in%Tiger.gen.samples$Bag_NO,select=c(BAG_NO,SPECIES,SHEET_NO,
                                                                  Month,year,BOTDEPTH,Mid.Lat,Mid.Long))
not.in.list=as.character(Tiger.gen.samples$Bag_NO[which(!Tiger.gen.samples$Bag_NO%in%Tiger.gen$BAG_NO)])



#sampling stations                   
Fixed.Stations=read.csv("C:/Matias/Data/Fixed station sites.csv")
Mid.Point=with(Fixed.Stations,gcIntermediate(cbind(Long.1,Lat.1), cbind(Long.2,Lat.2), n=1, addStartEnd=F))
Mid.Point=as.data.frame(do.call(rbind,Mid.Point))
names(Mid.Point)=c("Fix.St.mid.lon","Fix.St.mid.lat")
Fixed.Stations=cbind(Fixed.Stations,Mid.Point)
South.lat.bound=with(Fixed.Stations,min(c(Lat.1,Lat.2)))
North.lat.bound=with(Fixed.Stations,max(c(Lat.1,Lat.2)))
West.lon.bound=with(Fixed.Stations,min(c(Long.1,Long.2)))
East.lon.bound=with(Fixed.Stations,max(c(Long.1,Long.2)))

#Southern Oscillation Index
SOI=read.csv("C:/Matias/Data/Oceanography/SOI.csv")



#Mean Freo sea level
Freo=read.csv("C:/Matias/Data/Oceanography/Freo_mean_sea_level.csv")
names(Freo)[c(1,3)]=c("Year","Freo")

#Depth
Bathymetry_120=read.table("C:/Matias/Data/Mapping/get_data112_120.cgi")
Bathymetry_138=read.table("C:/Matias/Data/Mapping/get_data120.05_138.cgi") 
Bathymetry=rbind(Bathymetry_120,Bathymetry_138)

#Genetic tissue in stock
#library(xlsx)
#Gen.tissue=read.xlsx("M:/Fisheries Research/FinFish/Shark/Dani/Samples/Stored samples - STOCKTAKE.xlsx", 1)

#Shapefiles
PerthIs=read.table("C:/Matias/Data/Mapping/WAislandsPointsNew.txt", header=T) #function for reading txt file
Rottnest.Is=subset(PerthIs,ID%in%c("ROTT1"))
Garden.Is=subset(PerthIs,ID%in%c("ROTT3"))



# CONTROL SECTION ---------------------------------------------------------

export.dat.hard.drive="NO"  #control if exporting copy of data to external hard drive
bySEX="YES"  #control if size frequency separated by sex or not

#control if doing abundance analysis
Do.abundance="YES"
do.GAM="YES"  #define if using GAM or GLM approach
do.GLM="NO"
First.LL.yr=2001
These.Months=5:8
explore.cpue="NO"
Select.term="NO"  #Run if need to select terms, otherwise skip to next point
Select.error="NO"  #Run if selecting best error structure
Select.error.AIC="NO"
Do.K.fold.test="NO"  #control if doing k fold test
check.GAM="NO"
do.like.ratio.test="NO"  #see significance of each term (anova doesn't work)
do.MC="NO"  #get confidence bounds thru MC
niter=1000  #Number of iterations
do.boot="NO"  #get confidence bounds thru bootstrapping
n.boot=2e3   #Number of bootstrap iterations (for some species boot dat is rejected so bumped up iterations)
Do.tiff="NO"    #select figure extension
Do.jpeg="YES"

#control if doing sex ratio anlysis
do.sex.ratio="NO"

#control if doing multivariate analysis
Do.multivariate="NO"
Do.multivariate.glm="NO"

#control if doing ecosystem indicators
Do.ecosystems="NO"
if(Do.multivariate=="YES") Do.ecosystems="YES"


#species selection
MIN.per.yr=5     #minimum number of records per year
MIN.yrs=5        #minimum years with minimum number of records per year
MIN.RECORDS=100   #minimum total number of observations 

MaxDepth=450

#maximum distance and depth from station to classify as "fixed station"
Max.Dis=10   #(in km) 
Max.depth=3  #(in metres)



# PROCEDURE SECTION -------------------------------------------------------

#Check size distribution by gear
fn.check.size.gr=function(SPEC)
{
  a=subset(DATA,Mid.Lat<=(-26) & SPECIES==SPEC & !is.na(FL))
  a$Method=with(a,ifelse(Method=="LLL","LL",Method))
  
  fn.fig(paste("C:/Matias/Analyses/Surveys/Naturaliste_longline/outputs/Size.GN.LL_",SPEC,sep=""),2000,2400)
  par(mfcol=c(2,1),mai=c(.3,.55,.5,.1),oma=c(3,1.25,2,.1),las=1,mgp=c(.04,.6,0))
  
  with(subset(a,Method=="GN"),hist(FL,main="Gillnet",xlim=c(0,200),ylab="",xlab=""))
  with(subset(a,Method=="LL"),hist(FL,main="Longline",xlim=c(0,200),ylab="",xlab=""))
  mtext("Fork length (cm)",1,outer=T)
  dev.off()
}
fn.check.size.gr(SPEC='GM')
fn.check.size.gr(SPEC='WH')

#Output scalefish size data for Jef Norris
Jeff=DATA%>%filter(!BOAT%in%c("FLIN","HAM","HOU","NAT","RV BREAKSEA","RV Gannet","RV GANNET","RV SNIPE 2") &
                     Mid.Lat<=(-31) & Mid.Long>115.5 & Method=="GN" & !is.na(TL) & 
                     SPECIES%in%c("PS.T","QS.T","RS.T","BG.T"))%>%
            select(c("SPECIES","SCIENTIFIC_NAME","COMMON_NAME","TYPE","SHEET_NO","date","year",
                      "Month","TL","Method","BOAT","BLOCK","Mid.Lat","Mid.Long","BOTDEPTH",
                      "MESH_SIZE","MESH_DROP","NET_LENGTH","SOAK.TIME"))
colnames(Jeff)=tolower(colnames(Jeff))
write.csv(Jeff,"C:/Matias/Analyses/Catch and effort/Data_Resquests/Jeff N/Scalefish_size.csv",row.names = F)

#Depth distribution of sandbar off Perth
fn.chck.dep.range=function(SP,LAT.RANGE,LONG.RANGE,GEAR,YLIM,XLIM)
{
  a=subset(DATA,SPECIES==SP & Mid.Lat>=LAT.RANGE[1] & Mid.Lat<=LAT.RANGE[2] &Method%in%GEAR)
  Bath=subset(Bathymetry,V2>=LAT.RANGE[1] & V2<=LAT.RANGE[2] & V1>=LONG.RANGE[1] & V1<=LONG.RANGE[2] )
  Bath=Bath[order(Bath$V1,Bath$V2),]
  xbat=sort(unique(Bath$V1))
  ybat=sort(unique(Bath$V2))
  reshaped=as.matrix(reshape(Bath,idvar="V1",timevar="V2",v.names="V3", direction="wide"))
  
  smoothScatter(a[,match(c("Mid.Long","Mid.Lat"),names(a))], nrpoints = 0,ylab="Lat",xlab="Long",
                ylim=YLIM,xlim=XLIM,cex.axis=1.25,cex.lab=1.9)
  points(a$Mid.Long,a$Mid.Lat,cex=1.25,col="white")
  contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=Ylim,xlim=Xlim, zlim=c(-1,-300),
          nlevels = 3,labcex=2,lty = c(1,2,3,4),col=1,add=T)
  legend("topright",paste("Species=",SP),bty="n",cex=1.5)
  legend("topleft",c(paste("Number of shots=",length(unique(a$SHEET_NO))),
                     paste("Number of sharks=",nrow(a))),bty="n",cex=1.5)
  polygon(x=Rottnest.Is$Longitude,y=Rottnest.Is$Latitude,col="dark grey")  #add missing islands
  polygon(x=Garden.Is$Longitude,y=Garden.Is$Latitude,col="dark grey")
  
  #add location of best shots
  a$Mid.Lat_Mid.Long=with(a,paste(Mid.Lat,Mid.Long))
  TABL=table(a$Mid.Lat_Mid.Long)
  best=subset(a,Mid.Lat_Mid.Long%in%names(TABL[which(TABL>4)]))
  best.agg=aggregate(Number~Mid.Long+Mid.Lat,best,sum)
  with(best.agg,points(Mid.Long,Mid.Lat,cex=Number/2,col=2))
  legend("bottomright","highest catches",col=2,pch=1,cex=2,bty="n")
  return(best.agg)
}

# fn.chck.dep.range(SP="TK",LAT.RANGE=c(-33,-31),LONG.RANGE=c(114.5,116),GEAR=c("LL","GN"),
#                   XLIM=c(114.7,115.7),YLIM=c(-33,-30.99))
hnd.TK.Perth="C:/Matias/Analyses/Surveys/Naturaliste_longline/outputs/"
fn.fig(paste(hnd.TK.Perth,"Sandbar_off_Perth_LL_survey",sep=''),2400,2400)
aa=fn.chck.dep.range(SP="TK",LAT.RANGE=c(-33,-31),LONG.RANGE=c(114.5,116),GEAR=c("LL"),
                     XLIM=c(114.7,115.7),YLIM=c(-33,-30.99))
dev.off()
write.csv(aa,paste(hnd.TK.Perth,"Sandbar_off_Perth_LL_survey.csv",sep=''),row.names=F)


fun.dummy=function(YLIM,XLIM,scaler)
{
  a=subset(DATA,SPECIES=="TK" &Method%in%"LL")
  Last.trip=subset(DATA,year==2016&Method%in%"LL"&BOAT%in%"NAT")
  Last.trip=Last.trip[!duplicated(Last.trip$SHEET_NO),]
  
  KTCH=aggregate(Number~SHEET_NO+Mid.Long+Mid.Lat,a,sum)
  EFF=aggregate(cbind(SOAK.TIME,N.hooks)~SHEET_NO,a,max,na.rm=T)
  KTCH.EFF=merge(KTCH,EFF,by="SHEET_NO")
  KTCH.EFF$CPUE=with(KTCH.EFF,Number/(SOAK.TIME* N.hooks))
  KTCH.EFF=subset(KTCH.EFF,!CPUE=="Inf")
  
  
  plotMap(worldLLhigh, ylim=YLIM,xlim=XLIM,
          plt = c(.1, 1, 0.075, 1),
          col="grey80",tck = 0.025, 
          tckMinor = 0.0125, xlab="",ylab="",axes=F)
  box()
  with(KTCH.EFF,points(Mid.Long,Mid.Lat,cex=CPUE*scaler,col="grey30"))
  with(Fixed.Stations[1:20,],points(Fix.St.mid.lon,Fix.St.mid.lat,
                                    pch=19,col=2,cex=1.25))
  
  points(113.661,-24.884,cex=1.25,pch=17,col=3)
  text(113.661,-24.884,"Carnarvon",pos=4)
  points(114.162,-27.720,cex=1.25,pch=17,col=3)
  text(114.162,-27.720,"Kalbarri",pos=4)
  points(122.2359,-17.9614,cex=1.25,pch=17,col=3)
  text(122.3,-18.2,"Broome",pos=2, srt=90)
  points(115.86,-31.95,cex=1.25,pch=17,col=3)
  text(115.86,-31.95,"Perth",pos=4)
  
  #Add last years shots
  points(Last.trip$Mid.Long,Last.trip$Mid.Lat,col="blue",pch=19)
  legend("bottomright",c("Fixed stations","2016 trip"),
         pch=c(19,19),col=c("red","blue"),bty="n")
  
  with(Fixed.Stations[1:20,],text(Fix.St.mid.lon,Fix.St.mid.lat,Station.no.,
                                  pos=3,cex=1.1,col="firebrick"))
  if(YLIM[1]<(-29) & YLIM[2] <(-25))
  {
    s=subset(KTCH.EFF,Mid.Lat>=YLIM[1] & Mid.Lat<=YLIM[2],
             select=c(Mid.Long,Mid.Lat,CPUE))
    s$CPUE=round(s$CPUE,4)
    s=s[order(-s$CPUE),]
    return(s)
  }
}
fn.fig(paste(hnd.TK.Perth,"Sandbar_cpue",sep=''),2400,2400)
fun.dummy(YLIM=c(-34.5,-15.5),XLIM=c(112.25,123.25),scaler=20)
dev.off()

fn.fig(paste(hnd.TK.Perth,"Sandbar_cpue_Kalbarri",sep=''),2400,2400)
aa=fun.dummy(YLIM=c(-30.5,-26),XLIM=c(112.25,115.9),scaler=100)
dev.off()
write.csv(aa,paste(hnd.TK.Perth,"Sandbar_cpue_Kalbarri.csv",sep=''),row.names=F)

fn.fig(paste(hnd.TK.Perth,"Sandbar_cpue_First_leg",sep=''),2400,2400)
fun.dummy(YLIM=c(-26.5,-19),XLIM=c(112.25,116.85),scaler=20)
dev.off()



#Spatial plot of Hammerheads for species listing issue

hammers=c("HZ","HS","HG")
Hamm=subset(DATA,SPECIES%in%hammers)
Hamm$COMMON_NAME=as.character(Hamm$COMMON_NAME)
Hamm$COMMON_NAME=with(Hamm,ifelse(COMMON_NAME=="HAMMERHEAD (GREAT)","Great Hammerhead",COMMON_NAME))

#display spatially
fn.hamm=function(sP)
{
  d=subset(Hamm,SPECIES==sP)
  plot(d$Mid.Long,d$Mid.Lat,main="",ylim=c(-36,-12),xlim=c(113.5,128.5),ylab="",
       xlab="",cex.axis=1.5)
  polygon(x=c(113,129,129,113),y=c(-36,-36,-26,-26),col=rgb(.1,.1,.2,alpha=.3))
  mtext(unique(d$COMMON_NAME),3,-2,cex=1.65)
}
fn.fig("C:/Matias/Analyses/Surveys/Naturaliste_longline/outputs/Hammerheads",1400,2400)
par(mfcol=c(3,1),mai=c(.5,.6,.1,.1),las=1,mgp=c(2.5,0.65,0))
for(i in 1:length(hammers)) fn.hamm(sP=hammers[i])
mtext("Latitude",2,-2,outer=T,cex=1.65,las=3)
mtext("Longitude",1,-1.5,outer=T,cex=1.65)
dev.off()


#Prop by species north (for Catch-MSY)
Hamm.North=subset(Hamm,Mid.Lat>(-26))
Tab.hh.sp=table(Hamm.North$COMMON_NAME)
Tab.hh.sp=round(Tab.hh.sp/sum(Tab.hh.sp),2)

#drop unnecessary variables
drop=c("MESH_SIZE","MESH_DROP","NET_LENGTH")
DATA=DATA[,-match(drop,names(DATA))]


#Select boats (scientific survey only) within Fixed stations area
Boats=c("NAT","HOU","HAM","FLIN")   

#Drop non-longline methods
Drop.methods=c("GN","Mandy J","TW","DL")

#Select data for survevys in northwestern WA 
DATA=subset(DATA,BOAT%in%Boats & N.hooks>=0  &
              END1LATD>(South.lat.bound-.5) & END1LATD<(North.lat.bound+.5) &
              END1LNGD>(West.lon.bound-.5) & END1LNGD<(East.lon.bound+.5) )   


#Add Freo and SOI
DATA=merge(DATA,SOI, by.x=c("Month","year"),by.y=c("Month","Year"),all.x=T)
DATA=merge(DATA,Freo, by.x=c("Month","year"),by.y=c("Month","Year"),all.x=T)

#remove Gummy, missindentification
DATA=subset(DATA,!SPECIES%in%c("GM","CO"))

#remove unidentified species
DATA=subset(DATA,!SPECIES%in%c("XX","XX.T"))


#dropline data base
DATA.DL=subset(DATA, Method=="DL" & year<2001)   

#longline data base
DATA=subset(DATA,Method=="LL")   

Yrs.DL=sort(unique(DATA.DL$year))
Yrs.LL=sort(unique(DATA$year))

#Table 1. Number of species  by gear   
Sks=sort(unique(subset(DATA,Taxa=="Elasmobranch")$SPECIES))
Scalies=sort(unique(subset(DATA,Taxa=="Teleost")$SPECIES))
shark.other="XX"
unknown.species=c("SX")

DATA=subset(DATA,!SPECIES%in%unknown.species)
all.species=unique(c(DATA$SPECIES,DATA.DL$SPECIES))

N.species=function(DAT,what,SP)
{
  Sks.only=subset(DAT,SPECIES%in%what)
  Table.species=table(as.character(Sks.only$SPECIES))
  N=sum(Table.species)
  N.species=length(Table.species)
  Cum=cumsum(rev(sort(Table.species)))
  Cum=round(100*Cum/N,1)
  Table.species=data.frame(SPECIES=names(Table.species),Numbers=c(Table.species))
  Cum=data.frame(SPECIES=names(Cum),Cum.Percent=Cum)
  MERGED=merge(Table.species,Cum,by="SPECIES")
  
  if(SP=="shark")
  {
    #add sex ratio
    SEX=aggregate(Number~SPECIES+SEX,Sks.only,sum)
    SEX=reshape(SEX,v.names = "Number", idvar = "SPECIES",timevar = "SEX", direction = "wide")
    #SEX$Fem.to.male.ratio=paste("1:",round(SEX$Number.M/SEX$Number.F,3))
    MERGED=merge(MERGED,SEX,by="SPECIES")
    if(is.na(match("Number.M",names(MERGED)))) MERGED$Number.M=NA
    if(is.na(match("Number.F",names(MERGED)))) MERGED$Number.F=NA
    MERGED$Number.U=with(MERGED,Numbers-(Number.F+ Number.M))
    MERGED[is.na(MERGED)]=""
    MERGED$Number.U=with(MERGED,ifelse(Number.U==0,"",Number.U))
  }
  return(list(N=N,N.species=N.species,Table.species=MERGED))
}
Size.species=function(DAT,what,SP)
{
  Sks.only=subset(DAT,SPECIES%in%what)
  sp.unik=as.character(unique(Sks.only$SPECIES))
  DatA=data.frame(SPECIES=sp.unik)
  DatA$N.measured=DatA$Max=DatA$Mean=DatA$Min=NA
  for (p in 1:length(sp.unik))
  {
    a=subset(Sks.only,SPECIES==sp.unik[p])
    if(SP=="shark" &  !is.na(sum(a$FL,na.rm=T))& sum(a$FL,na.rm=T)>0)
    {
      DatA$Min[p]=round(min(a$FL[a$FL>0],na.rm=T),0)
      DatA$Mean[p]=round(mean(a$FL[a$FL>0],na.rm=T),0)
      DatA$Max[p]=round(max(a$FL[a$FL>0],na.rm=T),0)
      DatA$N.measured[p]=length(a$FL[!is.na(a$FL)])
    }
    if(SP=="scalie" &  !is.na(sum(a$TL,na.rm=T))& sum(a$TL,na.rm=T)>0)
    {
      DatA$Min[p]=round(min(a$TL[a$TL>0],na.rm=T),0)
      DatA$Mean[p]=round(mean(a$TL[a$TL>0],na.rm=T),0)
      DatA$Max[p]=round(max(a$TL[a$TL>0],na.rm=T),0)
      DatA$N.measured[p]=length(a$TL[!is.na(a$TL)])
    }
  }
  
  return(DatA)
}

#LL
LL.species=N.species(DATA,Sks,SP="shark")
LL.size=Size.species(DATA,Sks,SP="shark")

#DL
DL.species=N.species(DATA.DL,Sks,SP="shark")
DL.size=Size.species(DATA.DL,Sks,SP="shark")

#LL scalies
LL.species.scalies=N.species(DATA,Scalies,SP="scalie")

DATA$TL=ifelse(DATA$TL<5,NA,DATA$TL)
DATA$TL=with(DATA,ifelse(SPECIES=="NW.T" & TL> 100,NA,
             ifelse(SPECIES=="BD.T" & TL> 200,NA,
             ifelse(SPECIES=="TV.T" & TL> 130,NA,TL))))
LL.size.scalies=Size.species(DATA,Scalies,SP="scalie")

#DL scalies
#DL.species.scalies=N.species(DATA.DL,Scalies,SP="scalie")
#DL.size.scalies=Size.species(DATA.DL,Scalies,SP="scalie")
#no scalies caught in drum lines

Table.1.LL=merge(LL.species$Table.species,LL.size,by="SPECIES")
ID=match("SPECIES",names(Table.1.LL))
names(Table.1.LL)[-ID]=paste("LL.",names(Table.1.LL)[-ID],sep="")

Table.1.DL=merge(DL.species$Table.species,DL.size,by="SPECIES")
names(Table.1.DL)[-ID]=paste("DL.",names(Table.1.DL)[-ID],sep="")

Table.1=merge(Table.1.LL,Table.1.DL,by="SPECIES",all=T)
Table.1=merge(SPECIES.names,Table.1,by.x="Species",by.y="SPECIES",all.y=T)
Table.1=Table.1[order(Table.1$LL.Cum.Percent,-Table.1$LL.Numbers),]
Table.1[is.na(Table.1)]=""
names(Table.1)[match("LL.Numbers",names(Table.1))]=paste("LL.Numbers_",Yrs.LL[1],"-",Yrs.LL[length(Yrs.LL)],sep="")
names(Table.1)[match("DL.Numbers",names(Table.1))]=paste("DL.Numbers_",Yrs.DL[1],"-",Yrs.DL[length(Yrs.DL)],sep="")

Req.shrk=vapply(strsplit(Table.1$SCIENTIFIC_NAME, " ", fixed = TRUE), "[", "", 1)
Req.shrk=c(Table.1$Species[which(Req.shrk=="Carcharhinus")],"MI","TG","SE","LE","TA","GB")
Prop.requiem.sharks=Table.1%>%mutate(Requiem=ifelse(Species%in%Req.shrk,"YES","NO")) %>%
                              rename(N='LL.Numbers_2001-2017')%>%
                              group_by(Requiem) %>%
                              summarise(N = sum(N)) %>%
                              mutate(prop=N/sum(N)) %>%
                              as.data.frame()
                              

Table.1.scalies=merge(LL.species.scalies$Table.species,LL.size.scalies,by="SPECIES")
names(Table.1.scalies)[-ID]=paste("LL.",names(Table.1.scalies)[-ID],sep="")
Table.1.scalies=merge(SPECIES.names,Table.1.scalies,by.x="Species",by.y="SPECIES",all.y=T)
Table.1.scalies=Table.1.scalies[order(Table.1.scalies$LL.Cum.Percent,-Table.1.scalies$LL.Numbers),]
Table.1.scalies[is.na(Table.1.scalies)]=""


#number of droplines and longlines done
N.droplines=length(unique(DATA.DL$SHEET_NO))
N.longlines=length(unique(DATA$SHEET_NO))


#Create useful vars
YEAR=sort(unique(DATA$year))
N.yrs=length(YEAR) 
id.n=match(paste("LL.Numbers_",Yrs.LL[1],"-",Yrs.LL[length(Yrs.LL)],sep=""),names(Table.1))

SHOTS=as.character(unique(DATA$SHEET_NO))

#Combine a few common blacktips with black tip complex
Blk.tps=c("Carcharhinus limbatus","Carcharhinus limbatus/tilstoni","Carcharhinus tilstoni")

DATA$SPECIES.old=DATA$SPECIES
DATA.DL$SPECIES.old=DATA.DL$SPECIES
DATA$SPECIES=with(DATA,ifelse(SPECIES=="BL","BT",SPECIES))
DATA.DL$SPECIES=with(DATA.DL,ifelse(SPECIES=="BL","BT",SPECIES))

DATA$SPECIES=with(DATA,ifelse(SPECIES=="BL","BT",SPECIES))
DATA.DL$SPECIES=with(DATA.DL,ifelse(SPECIES=="BL","BT",SPECIES))

DATA$COMMON_NAME=with(DATA,ifelse(COMMON_NAME=="Common blacktip shark","Blacktip sharks",COMMON_NAME))
DATA.DL$COMMON_NAME=with(DATA.DL,ifelse(COMMON_NAME=="Common blacktip shark","Blacktip sharks",COMMON_NAME))

DATA$SCIENTIFIC_NAME=with(DATA,ifelse(SCIENTIFIC_NAME=="Carcharhinus limbatus/tilstoni","Carcharhinus limbatus",SCIENTIFIC_NAME))
DATA.DL$SCIENTIFIC_NAME=with(DATA.DL,ifelse(SCIENTIFIC_NAME=="Carcharhinus limbatus/tilstoni","Carcharhinus limbatus",SCIENTIFIC_NAME))


#Fix N.hook numbers
Hand.line.Shot=c("J00938","J01199","N00271","N00341","N00409","N00435","N00445","N00469",
                 "N00500","N00518","N00546","Y00047")
Beach.net=c("R00941","R00948","R00958","R00966","R00972","R00980","R00986","R00990","R00994","R00998",
            "J01129","J01137","J01150","J01152","J01165","J01168")
DATA$N.hooks.Fixed=with(DATA,
                     ifelse(SHEET_NO=="R00418",45,
                     ifelse(SHEET_NO=="R00454",39,
                     ifelse(SHEET_NO=="J00735",40,
                     ifelse(SHEET_NO=="J00754",48,
                     ifelse(SHEET_NO=="J00789",42,
                     ifelse(SHEET_NO%in%Hand.line.Shot,1,
                     ifelse(SHEET_NO%in%Beach.net,0,N.hooks))))))))

DATA=subset(DATA,!is.na(Mid.Lat))   

D.map=DATA

#Define fixed stations
Check.Station=Fixed.Stations[,match(c("Station.no.","Fix.St.mid.lon","Fix.St.mid.lat","Min.depth","max.depth","Point.depth"),names(Fixed.Stations))]
N.rowCheck.St=nrow(Check.Station)
Depth.dispersal=rep(NA,N.rowCheck.St)
for(i in 1:N.rowCheck.St)Depth.dispersal[i]=(Check.Station$max.depth[i]-Check.Station$Min.depth[i])
Depth.dispersal=round(mean(Depth.dispersal,na.rm=T))/2
Check.Station$Min.depth=with(Check.Station,ifelse(is.na(Min.depth),Point.depth-Depth.dispersal,Min.depth))
Check.Station$max.depth=with(Check.Station,ifelse(is.na(max.depth),Point.depth+Depth.dispersal,max.depth))


Shots.loc=DATA[,match(c("SHEET_NO","Mid.Lat","Mid.Long","BOTDEPTH"),names(DATA))]
Shots.loc=Shots.loc[!duplicated(Shots.loc$SHEET_NO),]
N.rowShots.loc=nrow(Shots.loc)

#replicate data frames to have all combos of sheet# and fixed stations
Shots.loc=do.call("rbind", replicate(N.rowCheck.St,Shots.loc,simplify = FALSE))
Shots.loc=Shots.loc[order(Shots.loc$SHEET_NO),]

Check.Station$Station.no.=as.character(Check.Station$Station.no.)
Check.Station$Station.no.=with(Check.Station,ifelse(Station.no.=="additional",paste(Station.no.,rownames(Check.Station)),Station.no.))
Check.Station=do.call("rbind", replicate(N.rowShots.loc,Check.Station,simplify = FALSE))

Shots.loc=cbind(Shots.loc,Check.Station)

Shots.loc$Dist.Station.km=distCosine(cbind(Shots.loc$Mid.Long,Shots.loc$Mid.Lat),cbind(Shots.loc$Fix.St.mid.lon,Shots.loc$Fix.St.mid.lat))/1000

Shots.loc$FixedStation=with(Shots.loc,ifelse(Dist.Station.km<=Max.Dis,"YES","NO"))
Shots.loc$FixedStation=with(Shots.loc,ifelse(FixedStation =="YES" 
        & (BOTDEPTH>=(Min.depth-Max.depth) & BOTDEPTH<=(max.depth+Max.depth)),"YES",FixedStation))


Yes.Fix.ST=subset(Shots.loc,FixedStation=="YES",select=c("SHEET_NO","FixedStation","Station.no."))
Yes.Fix.ST=Yes.Fix.ST[!duplicated(Yes.Fix.ST$SHEET_NO),]
NO.Fix.ST=subset(Shots.loc,FixedStation=="NO",select=c("SHEET_NO","FixedStation","Station.no."))
NO.Fix.ST$Station.no.=NA
NO.Fix.ST=NO.Fix.ST[!duplicated(NO.Fix.ST$SHEET_NO),]
NO.Fix.ST=subset(NO.Fix.ST,!(SHEET_NO%in%Yes.Fix.ST$SHEET_NO))
Yes.No.Fix.ST=rbind(Yes.Fix.ST,NO.Fix.ST)

DATA=merge(DATA,Yes.No.Fix.ST,by="SHEET_NO",all.x=T)


#set all lines within a station to Fix or not
UniC=DATA[!duplicated(DATA$SHEET_NO),]
Fixed.only=subset(UniC,FixedStation=="YES")
Fixed.only=unique(Fixed.only$date)
DATA$FixedStation=with(DATA,ifelse(date%in%Fixed.only,"YES",FixedStation))

#Set additional stations to non-fixed as only a few years revisited
DATA$FixedStation=with(DATA,ifelse(!Station.no.%in%as.character(1:20),"NO",FixedStation))


setwd('C:/Matias/Analyses/Surveys/Naturaliste_longline/outputs/Abundance')

#see what records were set to fixed stations
hndl.expl=paste(getwd(),"Exploratory",sep="/")
pdf(paste(hndl.expl,"Kept_Stations.pdf",sep="/"))
with(subset(DATA,FixedStation=="YES"),plot(Mid.Long,Mid.Lat))
points(Fixed.Stations$Fix.St.mid.lon,Fixed.Stations$Fix.St.mid.lat,pch=19,cex=.5,col=2)
with(subset(DATA,FixedStation=="NO"),points(Mid.Long,Mid.Lat,pch=19,cex=.5,col=3))
legend("topleft",c("used","station","dropped"),pch=c(1,19,19),col=c(1,2,3),bty='n')
dev.off()

Table.dropped=table(DATA$year,DATA$FixedStation)
PROPR.dropped=Table.dropped[,1]/rowSums(Table.dropped)
write.csv(PROPR.dropped,"prop.not.fixed.station.csv")

#Check stations visited each year (should be ~ 5 shots per station per year)
TAB=with(DATA[!duplicated(paste(DATA$SHEET_NO,DATA$Station.no.)),],table(year,Station.no.))
pdf(paste(hndl.expl,"Stations_by_year.pdf",sep="/"))
tmp <- TAB
mytheme <- gridExtra::ttheme_default(
  core = list(padding=unit(c(1, 1), "mm"),fg_params=list(cex = .5)),
  colhead = list(fg_params=list(cex = .7)),
  rowhead = list(fg_params=list(cex = .7)))
myt <- gridExtra::tableGrob(tmp[,1:19], theme = mytheme)
grid.draw(myt)
plot.new()
myt <- gridExtra::tableGrob(tmp[,20:ncol(tmp)], theme = mytheme)
grid.draw(myt)
dev.off()

#Select species for analysis
Species.by.yr=with(subset(DATA,FixedStation=="YES" & TYPE=="Elasmo"),table(as.character(SPECIES),year))
Species.by.yr.Min.per.year=Species.by.yr
Species.by.yr.Min.per.year[Species.by.yr.Min.per.year<MIN.per.yr]=0
Species.by.yr.Min.per.year[Species.by.yr.Min.per.year>=MIN.per.yr]=1
Tot.records.yr=rowSums(Species.by.yr.Min.per.year)
Tot.records.yr=rev(sort(Tot.records.yr))
TARGETS=names(Tot.records.yr[Tot.records.yr>=MIN.yrs])
TARGETS.name=subset(Table.1,Species%in%TARGETS,select=c(Species,COMMON_NAME))


Tot.records=rowSums(Species.by.yr)
Tot.records=rev(sort(Tot.records))

Tot.yrs.pos.ktch=Species.by.yr
Tot.yrs.pos.ktch[Tot.yrs.pos.ktch>0]=1
Tot.yrs.pos.ktch=rowSums(Tot.yrs.pos.ktch)
Tot.yrs.pos.ktch=sort(Tot.yrs.pos.ktch)

IID=match(names(Tot.records[Tot.records>=MIN.RECORDS]),TARGETS)
TARGETS=TARGETS[IID]
TARGETS.name=subset(TARGETS.name,Species%in%TARGETS)
TARGETS=TARGETS.name$Species
TARGETS.name=TARGETS.name$COMMON_NAME
names(TARGETS)=TARGETS.name
N.species=length(TARGETS)


#Number of species per stations
dummy=subset(DATA,FixedStation=="YES",select=c(year,SPECIES,Mid.Long,Mid.Lat))
dummy$dummy=with(dummy,paste(substr(Mid.Lat,1,5),substr(Mid.Long,1,5)))  
TAB=table(dummy$dummy)
#plot(as.numeric(substr(names(TAB),7,11)),as.numeric(substr(names(TAB),1,5)),cex=5*TAB/max(TAB),col="yellow")
rm(dummy)


#table Fixed Stations per year
write.csv(table(UniC$year,UniC$FixedStation),"Table.fixed Vs NotFixed.csv")
          
if(export.dat.hard.drive=="YES")write.csv(DATA,"E:/Matias_WA_Fisheries/Matias/Analyses/Surveys/DATA.csv",row.names=F)


DATA$Moon=as.character(DATA$Moon)


#Export nice table
Export.tbl=function(WD,Tbl,Doc.nm,caption,paragph,HdR.col,HdR.bg,Hdr.fnt.sze,Hdr.bld,
                    body.fnt.sze,Zebra,Zebra.col,Grid.col,Fnt.hdr,Fnt.body,
                    HDR.names,HDR.span,HDR.2nd)
{
  mydoc = docx(Doc.nm)  #create r object
  mydoc = addSection( mydoc, landscape = T )   #landscape table
  # add title
  if(!is.na(caption))mydoc = addParagraph(mydoc, caption, stylename = "TitleDoc" )
  
  # add a paragraph
  if(!is.na(paragph))mydoc = addParagraph(mydoc , paragph, stylename="Citationintense")
  
  #add table
  MyFTable=FlexTable(Tbl,header.column=F,add.rownames =F,
                     header.cell.props = cellProperties(background.color=HdR.bg), 
                     header.text.props = textProperties(color=HdR.col,font.size=Hdr.fnt.sze,
                                                        font.weight="bold",font.family =Fnt.hdr), 
                     body.text.props = textProperties(font.size=body.fnt.sze,font.family =Fnt.body))
  
  #Add header
  MyFTable = addHeaderRow(MyFTable,text.properties=textBold(),value=HDR.names,colspan=HDR.span)
  
  #Add second header
  MyFTable = addHeaderRow(MyFTable, text.properties = textBold(),value =HDR.2nd)
  
  
  # zebra stripes - alternate colored backgrounds on table rows
  if(Zebra=="YES") MyFTable = setZebraStyle(MyFTable, odd = Zebra.col, even = "white" )
  
  # table borders
  MyFTable = setFlexTableBorders(MyFTable,
                                 inner.vertical = borderNone(),inner.horizontal = borderNone(),
                                 outer.vertical = borderNone(),
                                 outer.horizontal = borderProperties(color=Grid.col, style="solid", width=4))
  
  # set columns widths (in inches)
  #MyFTable = setFlexTableWidths( MyFTable, widths = Col.width)
  
  mydoc = addFlexTable( mydoc, MyFTable)   
  mydoc = addSection( mydoc, landscape = F ) 
  
  # write the doc 
  writeDoc( mydoc, file = paste(Doc.nm,".docx",sep=''))
}


#interpolate missing Temperature (linear interpolation)
Unk.Temp=subset(DATA,!is.na(TEMP),select=c(SHEET_NO,TEMP))
Unk.Temp=Unk.Temp[!duplicated(Unk.Temp$SHEET_NO),]
add.Sheet.no=unique(DATA$SHEET_NO)
add.Sheet.no=add.Sheet.no[which(!add.Sheet.no%in%Unk.Temp$SHEET_NO)]
add.Sheet.no=data.frame(SHEET_NO=add.Sheet.no,TEMP=NA)
add.Sheet.no$SHEET_NO=as.character(add.Sheet.no$SHEET_NO)
Unk.Temp=rbind(Unk.Temp,add.Sheet.no)
Unk.Temp=Unk.Temp[order(Unk.Temp$SHEET_NO),]


Unk.Temp$TEMP1=na.interpolation(ts(Unk.Temp$TEMP))
DATA=merge(DATA,subset(Unk.Temp,select=c(SHEET_NO,TEMP1)),by='SHEET_NO')


#Set some predictors as factors
DATA$depth.bin=10*round(DATA$BOTDEPTH/10)
DATA$TEMP1.bin=floor(DATA$TEMP1)
DATA$Mid.Lat.bin=floor(abs(DATA$Mid.Lat))


# 1. CATCH RATES FROM FISHERY INDEPENDENT SURVEYS ----------------------------------------------------------------------
#References: Zuur et al. Zero inflated models and generalized linear mixed models with R
if(Do.abundance=="YES")
{
  #1.1.  Figure 2. Size freq of each species to show its adult abundance 
  SizeFreq.fn=function(SPEC,Species.names,MAT,BIN)
  {
    datos=subset(DATA,SPECIES==SPEC)
    
    if(bySEX=="YES")
    {
      Males=subset(datos,SEX=="M")
      Females=subset(datos,SEX=="F")
      Rango=range(c(Males$FL,Females$FL),na.rm=T)
      Rango[1]=floor(Rango[1]/BIN)*BIN
      Rango[2]=ceiling(Rango[2]/BIN)*BIN
      Bks=seq(Rango[1],Rango[2],BIN)
      Mal=hist(Males$FL,breaks=Bks,plot=F)
      Fem=hist(Females$FL,breaks=Bks,plot=F)
      YMAX=max(Fem$counts,Mal$counts)
      xvals=barplot(rbind(Fem$counts,Mal$counts),beside=T,ylim=c(0,YMAX*1.05),
                    names.arg=Mal$mids,col=c("white","grey70"),cex.axis=1.25,cex.names=1.25)
      axis(1,xvals[1,]+0.5,F)
      
      YARROW=round(YMAX*.2)
      id=min(which(abs(Fem$mids-MAT)==min(abs(Fem$mids-MAT))))#find closest value
      arrows(xvals[1,id],YARROW,xvals[1,id],0,lwd=2,length = 0.10)
      
    }
    
    if(bySEX=="NO")
    {
      a=hist(datos$FL,breaks=seq(20,340,BIN),main="",ylab="",xlab="",col="grey",xlim=c(0,XLIM),cex.axis=1.25)
      #a=hist(datos$FL,breaks=BKS,main="",ylab="",xlab="",col="grey",xaxt="n",cex.axis=1.25)
      #axis(1,at=BKS,labels=F,tck=-0.035,las=1,cex.axis=1.5)
    }
    
    
    box()
  }
  
  Tar.names=TARGETS.name
  
  #add 50% maturity
  Fem.Size.Mat=list()
  Fem.Size.Mat$"Sandbar shark"=135
  Fem.Size.Mat$"Milk shark"=66
  Fem.Size.Mat$"Spot-tail shark"=76
  Fem.Size.Mat$"Tiger shark"=260 
  Fem.Size.Mat$"Blacktip sharks"=99
  Fem.Size.Mat$"Scalloped hammerhead"=210  
  Fem.Size.Mat$"Dusky shark"=254
  Fem.Size.Mat$"Sliteye shark"=60
  Fem.Size.Mat$"Great hammerhead"=220
  Fem.Size.Mat$"Silvertip shark"=195
  Fem.Size.Mat$"Pigeye shark"=215
  Fem.Size.Mat$"Whitespot shovelnose"=155
  Fem.Size.Mat$"Grey reef shark"=135
  Fem.Size.Mat$"Lemon shark"=220
  
  iDs=match(TARGETS.name,names(Fem.Size.Mat))
  Fem.Size.Mat=Fem.Size.Mat[iDs]
  
 
  #1.2.  Temporal variability of size by year (including Kruskal Wallis test)
  #Lat.ranges=list(c(-20,-24),c(-20,-24),c(-20,-24),c(-20,-23),c(-20,-22),c(-20,-23),c(-20,-24))
  SizeTemp.fn=function(SPEC,Species.names,LATS)
  {
    datos=subset(DATA,SPECIES==SPEC)
    #  datos=subset(DATA,SPECIES==SPEC & Mid.Lat>=LATS[2] & Mid.Lat<=LATS[1])
    a=subset(datos,!is.na(FL),select=c(year,FL))
    Id=unique(a$year)
    ID=YEAR[-which(YEAR%in%Id)]
    if(length(ID)>0)a=rbind(a,data.frame(year=ID,FL=NA))
    a$year=as.factor(a$year)
    boxplot(a$FL~a$year,main="",ylab="",xlab="",varwidth=T,col="grey80",
            ylim=c(min(a$FL,na.rm=T)*0.95,max(a$FL,na.rm=T)))
    Ns=tapply(a$FL,a$year,length)
    Ns=ifelse(Ns==1,0,Ns)
    text(1:N.yrs,rep(min(a$FL,na.rm=T)*0.95,N.yrs),Ns,cex=0.65)
    KW=kruskal.test(FL ~ year, data = datos)
    return(KW)
  }
  Kruskal=vector('list',N.species)
  names(Kruskal)=TARGETS
  
  
  #1.3   Sex ratios
  Chi.squar=function(SPEC)
  {
    datos=subset(DATA,SPECIES==SPEC & SEX%in%c("M","F"))
    TABLA=table(datos$SEX)
    X2=chisq.test(TABLA)
    return(list(numbers=TABLA,X2=X2))
  }
  StoreX2=Kruskal
  for(i in 1:N.species)StoreX2[[i]]=Chi.squar(TARGETS[i])
  
  #1.4   Table of month by years
  Table.yr.Mn=with(DATA,table(year,Month))
  
  
  #1.5    Put data as list and get 0 shots
  DATA.list=vector("list",length=N.species)
  names(DATA.list)=TARGETS.name
  Prop.Ktch.list=Nom.dummy=DATA.list
  
  # #species ranges                             
  # DISTRIB$"Dusky shark"=c(-25.7,-17.35)
  # DISTRIB$"Sandbar shark"=c(-27,-16.86)
  # DISTRIB$"Tiger shark"=c(-27,-15.6)
  
  
  #Depth distributions
  Z.dist=function(SPEC)
  {
    dat=subset(DATA,SPECIES==SPEC)
    a=density(dat$BOTDEPTH,na.rm=T,adjust = 2,from=0, to=MaxDepth)  
  }
  
  
  #Table sheet number vs number of lines per sheet
  Tabla.Sheet.Lines=table(as.character(DATA$SHEET_NO),as.numeric(DATA$LINE_NO))
  ColSum.Sheet.Lines=rowSums(Tabla.Sheet.Lines)
  hist(ColSum.Sheet.Lines,breaks=seq(1,37,1),xlab="Number of lines per sheet")
  
  #Methods used by year
  DATA$Method=as.character(DATA$Method)
  Yr.Method=with(DATA,table(Method,year,useNA="ifany"))
  
  #add Temperature residuals to remove correlation between Temp and Month
  DATA=DATA%>%group_by(Month,year)%>%
    mutate(Temp.res=TEMP1/mean(TEMP1,na.rm=T))%>%
    as.data.frame()
  
  # Select relevant variables
  THESEVARS=c("SHEET_NO","LINE_NO","date","Month","year","BOAT","BOTDEPTH","depth.bin",
              "Moon","SOI","Freo","BLOCK","Mid.Lat","Mid.Long","N.hooks","TEMP1",
              "TEMP1.bin","Temp.res","Mid.Lat.bin",
              "N.hooks.Fixed","SOAK.TIME","Method","FixedStation","Station.no.")
  
  #1.6    Construct wide database for analysis
  Effort.data.fun=function(target,SUBSET)
  {
    #remove record if no effort data or outside distribution
    if(SUBSET=="NO")DAT=subset(DATA,N.hooks.Fixed>20)
    if(SUBSET=="YES")DAT=subset(DATA,N.hooks.Fixed>20 & BOTDEPTH <MaxDepth &
                                  Month%in%These.Months & Set.time <"08:00" & !(year==2004))
    DAT=subset(DAT,!is.na(SPECIES))
    
    #target species catch 
    DAT$Catch.Target=with(DAT,ifelse(SPECIES==target,Number,0))
    
    #other species catch
    DAT$Catch.Dusky=with(DAT,ifelse(SPECIES=="BW",Number,0))
    DAT$Catch.Sandbar=with(DAT,ifelse(SPECIES=="TK",Number,0))
    DAT$Catch.Tiger=with(DAT,ifelse(SPECIES=="TG",Number,0))
    DAT$Catch.Scalloped=with(DAT,ifelse(SPECIES=="HS",Number,0))
    DAT$Catch.Blacktip=with(DAT,ifelse(SPECIES=="BT",Number,0))
    DAT$Catch.SpotTail=with(DAT,ifelse(SPECIES=="SO",Number,0))
    DAT$Catch.Milk=with(DAT,ifelse(SPECIES=="MI",Number,0))
    
    #reshape catch data
    TABLE=aggregate(cbind(Catch.Target,Catch.Dusky,Catch.Sandbar,Catch.Tiger,Catch.Scalloped,
                          Catch.Blacktip,Catch.SpotTail,Catch.Milk)~SHEET_NO,data=DAT,sum,na.rm=T)
    
    DAT1=DAT[!duplicated(DAT$SHEET_NO),match(THESEVARS,names(DAT))]
    TABLE=merge(TABLE,DAT1,by="SHEET_NO",all.x=TRUE)
    
    #proportion of records with target catch
    prop.with.catch=round(100*sum(TABLE$Catch.Target>0)/length(TABLE$Catch.Target),0)
    return(list(dat=TABLE,prop.ktch=prop.with.catch))
  }
  for ( i in 1:N.species)
  {
    dum=Effort.data.fun(TARGETS[i],SUBSET="NO")
    DATA.list[[i]]=dum$dat
    Prop.Ktch.list[[i]]=dum$prop.ktch
  }
  
  
  #1.7    Plot location of positive catches
  
  #By species and year
  REFS=list(c(1,5,10),c(10,20,40),c(1,5,10),c(1,5,10),c(5,10,20),c(1,5,10),c(1,5,10))
  SCALES=c(5,1,5,5,5,5,5)
  pos.catch.pos=function(Dat,REF,SCALER)
  {
    Dat1=Dat
    Dat$Catch.Target=with(Dat,ifelse(Catch.Target>0,log(Catch.Target+SCALER),0))
    
    Xlim=range(Dat$Mid.Long,na.rm=T)
    Ylim=range(Dat$Mid.Lat,na.rm=T)
    LEG=c(as.character(REF[1]),"",paste("",as.character(REF[2])),"",paste("",as.character(REF[3])))
    
    Yr=sort(unique(Dat$year))
    n.Yr=length(Yr) 
    fn.fig(paste(TARGETS.name[i],"dist.by.yr",sep=""),2400,2400)
    par(mfcol=c(4,4),las=1,mai=c(.3,.35,.01,.01),oma=c(3,2,1,1),mgp=c(.05,.7,0))
    for(p in 1:n.Yr)
    {
      a=subset(Dat,year==Yr[p])
      a1=subset(Dat1,year==Yr[p])
      
      plot(a$Mid.Long,a$Mid.Lat,cex=a$Catch.Target,xlim=Xlim,ylim=Ylim,ylab="",xlab="",
           xaxt='n',yaxt='n',pch=21,bg="grey80")
      legend("topleft",paste(Yr[p]),bty='n',cex=1.25) 
      legend("bottomright",paste("N=",sum(a1$Catch.Target)),bty='n',cex=1.25)
      axis(1,seq(round(Xlim[1]),round(Xlim[2]),2),F)
      axis(2,seq(round(Ylim[1]),round(Ylim[2]),2),F)
      if(p %in%c(4,8,12))axis(1,seq(round(Xlim[1]),round(Xlim[2]),2),seq(round(Xlim[1]),round(Xlim[2]),2))
      if(p %in%1:4)axis(2,seq(round(Ylim[1]),round(Ylim[2]),2),-seq(round(Ylim[1]),round(Ylim[2]),2))
      
    }
    plot(1:10,1:10,bty='n',pch='',xaxt='n',yaxt='n',ann=F)
    legend("center",TARGETS.name[i],text.col=1,bty='n',cex=2)
    
    plot(1:10,1:10,bty='n',pch='',xaxt='n',yaxt='n',ann=F)
    legend("center",LEG,pch=21,col=c("black","white","black","white","black"),
           pt.cex=c(log(REF[1]+SCALER),1,log(REF[2]+SCALER),1,log(REF[3]+SCALER)),bty='n',cex=1.75)
    mtext("Longitude",side=1,line=1,font=1,las=0,cex=1.25,outer=T)
    mtext("Latitude",side=2,line=0.5,font=1,las=0,cex=1.25,outer=T)
    dev.off()
  }
  
  #Species combined
  pos.catch.pos=function(Dat,REF,SCALER)
  {
    Dat1=Dat
    Dat$Catch.Target=with(Dat,ifelse(Catch.Target>0,log(Catch.Target+SCALER),0))
    
    Xlim=range(Dat$Mid.Long,na.rm=T)
    Ylim=range(Dat$Mid.Lat,na.rm=T)
    LEG=c(as.character(REF[1]),"",paste("",as.character(REF[2])),"",paste("",as.character(REF[3])))
    
    plotMap(worldLLhigh, xlim=c(Xlim[1]*0.995,Xlim[2]*1.0025),ylim=c(Ylim[1]*1.015,Ylim[2]*0.995),
            plt = c(.1, .99, 0.1, 1),col="grey65",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
    points(Dat$Mid.Long,Dat$Mid.Lat,cex=Dat$Catch.Target,pch=21,bg="grey95")
    axis(1,seq(round(Xlim[1]),round(Xlim[2]),2),F)
    axis(2,seq(round(Ylim[1]),round(Ylim[2]),2),F)
    if(i%in%c(3,6,7))axis(1,seq(round(Xlim[1]),round(Xlim[2]),2),seq(round(Xlim[1]),round(Xlim[2]),2),cex.axis=1.25)
    if(i%in%1:3)axis(2,seq(round(Ylim[1]),round(Ylim[2]),2),-seq(round(Ylim[1]),round(Ylim[2]),2),cex.axis=1.25)
    box()
    mtext("Longitude",side=1,line=1,font=1,las=0,cex=1.85,outer=T)
    mtext("Latitude",side=2,line=0,font=1,las=0,cex=1.85,outer=T)
    mtext(Tar.names[i],side=3,line=0.05,font=1,las=0,cex=1.15,outer=F)
    legend("topleft",paste("N=",sum(Dat1$Catch.Target)),bty='n',cex=1.5)
    legend("bottomright",LEG,pch=21,col=c("black","transparent","black","transparent","black"),
           pt.cex=c(log(REF[1]+SCALER),1,log(REF[2]+SCALER),1,log(REF[3]+SCALER)),bty='n',cex=1.25)
    
  }
  
  
  #1.8   Figure 1. Map of Catch rates species combined
  pos.cpue=function(Dat,SCALER,ADD.zero)
  {
    Dat$CPUE=Perhooks*Dat$Catch.Target/(Dat$N.hooks.Fixed*Dat$SOAK.TIME)
    ZeroCatch=subset(Dat,Catch.Target==0)
    Dat=subset(Dat,CPUE>0)
    x=ceiling(quantile(Dat$CPUE,probs=c(.5,1)))
    Dat$CPUE=fn.scale(Dat$CPUE,max(Dat$CPUE,na.rm=T),SCALER)
    
    
    a=c(Xlim[1]*0.995,Xlim[2]*1.0001)
    b=c(Ylim[1]*1.015,Ylim[2]*0.995)
    PLATE=c(.01,.9,.075,.9)
    plotmap(a,b,PLATE,Col.land,a,b)
    
    #Add shots
    points(Dat$Mid.Long,Dat$Mid.Lat,cex=Dat$CPUE,pch=21,bg=sht.col,col="grey60")
    if(ADD.zero=="YES")points(ZeroCatch$Mid.Long,ZeroCatch$Mid.Lat,cex=0.75,pch=4,col="grey35")
    axis(1,seq(round(Xlim[1]),round(Xlim[2]),2),F)
    axis(2,seq(round(Ylim[1]),round(Ylim[2]),2),F)
    box()
    
    if(add.depth=="YES")contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=Ylim,xlim=Xlim, zlim=c(-1,-300),
                                nlevels = 3,labcex=.1,lty = c(1,2,3),col=c("gray20","gray20","gray20","transparent"),add=T)
 
    legend("bottomright",paste(x),pch=21,col="grey60",
            pt.cex=mapply(fn.scale,x,max(x),SCALER),bty='n',cex=1,
            pt.bg=sht.col,title="catch rate")
    
  }
  Bathymetry=Bathymetry[order(Bathymetry$V1,Bathymetry$V2),]
  xbat=sort(unique(Bathymetry$V1))
  ybat=sort(unique(Bathymetry$V2))
  
  
  #1.9   Figure 4. Proportion with and without catch for each species
  plot0=function(DAT,Ymax,SPECIES)
  {
    DATA.no0=subset(DAT,Catch.Target>0)
    Tab=table(DATA.no0$Catch.Target)
    Ns=as.numeric(names(Tab))
    id=1:Ns[length(Ns)]
    id1=which(!id%in%Ns)
    id2=rep(0,length(id1))
    names(id2)=id1
    Tab=c(Tab,id2)
    Tab=Tab[match(id,names(Tab))]
    Prop=Tab/nrow(DAT)
    a=barplot(Prop,ylim=c(0,Ymax),las=1,cex.axis=1.25,cex.names=1.25)
    axis(1,at=a[1:(length(Prop)+1)],labels=F,tck=-0.032,cex.axis=1.25)
    legend("top",SPECIES,bty='n',cex=1.25)
    Catch.0=round((nrow(DAT)-nrow(DATA.no0))/nrow(DAT),2)
    legend("bottomright",paste("Prop. 0=",round(Catch.0,2)),bty='n',cex=1.1)
    box()
  }
  
  
  
  #1.10   ---- Nominal cpue ----#
  
  #create data sets for cpue standardisation
  for ( i in 1:N.species)
  {
    dum=Effort.data.fun(TARGETS[i],SUBSET="YES")
    DATA.list[[i]]=dum$dat
    Prop.Ktch.list[[i]]=dum$prop.ktch
  }
  
  YEAR=sort(unique(DATA.list[[i]]$year))
  N.yrs=length(YEAR) 
  
  #1.10.1. function for calculating cpue (all records)
  All=function(DAT)
  {
    DAT$CPUE=DAT$Catch.Target/(DAT$N.hooks.Fixed*DAT$SOAK.TIME)
    return(DAT[,match(c("year","CPUE"),names(DAT))])
  }
  
  All.Nom.dummy.FS=All.Nom.dummy=Nom.dummy
  for ( i in 1:N.species)
  {
    All.Nom.dummy[[i]]=All(DATA.list[[i]])
    All.Nom.dummy.FS[[i]]=All(subset(DATA.list[[i]],FixedStation=="YES"))
  }
  
  #1.10.2. Calculate annual mean CPUE (aggregated by year)
  Nom.dummy.FS=Nom.dummy
  
  Nominal=function(DAT)CPUE=aggregate(Catch.Target/(N.hooks.Fixed*SOAK.TIME)~year,DAT,mean)
  
  for ( i in 1:N.species)
  {
    Nom.dummy[[i]]=Nominal(DATA.list[[i]])
    Nom.dummy.FS[[i]]=Nominal(subset(DATA.list[[i]],FixedStation=="YES"))
  }
  
  #show all data and mean
  if(explore.cpue=="YES")
  {
    plot.nom.cpue=function(A,B,LEGE,ALL)
    {
      if(ALL=="YES")YLIM=c(0,max(A$CPUE))
      if(ALL=="NO")YLIM=c(0,max(B[,2]))
      plot(A$year,A$CPUE,ylab="",xlab="",ylim=YLIM)
      lines(B$year,B[,2],lwd=2,col=2)
      legend("topright",LEGE,bty='n',cex=1.15)
    }
    par(mfcol=c(4,3),mai=c(.3,.4,.01,.1),oma=c(3,2,2,.1),las=1,mgp=c(.04,.6,0))
    for ( i in 1:N.species) plot.nom.cpue(All.Nom.dummy[[i]],Nom.dummy[[i]],Tar.names[i],"NO")
    mtext("Year",side=1,line=1,font=1,las=0,cex=1.85,outer=T)
    mtext("CPUE",side=2,line=0,font=1,las=0,cex=1.85,outer=T)
  }
  
  
  #Mean All stations Vs Fixed stations
  hndLs="C:/Matias/Analyses/Surveys/Naturaliste_longline/outputs/Abundance/"
  
  if(explore.cpue=="YES")
  {
    fn.fig(paste(hndLs,"Nominal.cpue.All_vs_fixed.stations",sep=""),2000,2400)
    par(mfcol=c(5,3),mai=c(.3,.4,.1,.1),oma=c(3,2,2,.1),las=1,mgp=c(.04,.6,0))
    for ( i in 1:N.species)
    {
      YMAX=max(Nom.dummy[[i]][,2],Nom.dummy.FS[[i]][,2],na.rm=T)
      plot(YEAR,Nom.dummy[[i]][,2],type='l',ylab="",xlab="",main="",ylim=c(0,YMAX),xaxt='n')
      lines(YEAR,Nom.dummy.FS[[i]][,2],col=2)
      axis(1,YEAR,YEAR)
      mtext(Tar.names[i],side=3,line=.25,font=1,las=0,cex=.9)
    }
    plot(1:10,1:10,bty='n',pch='',xaxt='n',yaxt='n',ann=F)
    legend("center",c("All stations","Fixed stations only"),lty=1,col=1:2,bty='n',cex=1.5)
    mtext("Year",side=1,line=1,font=1,las=0,cex=1.85,outer=T)
    mtext("CPUE",side=2,line=0,font=1,las=0,cex=1.85,outer=T)
    dev.off()
  }

  
  #Effort series
  Get.effort=function(DAT)EFF=aggregate(N.hooks.Fixed*SOAK.TIME~year,DAT,sum)
  Eff.All.stations=Get.effort(DATA.list[[2]])
  Eff.Fixed.stations=Get.effort(subset(DATA.list[[2]],FixedStation=="YES"))
  
  
  
  #1.11   ---- Standardised cpue ----#
  cfac=function(x,breaks=NULL)
  {
    x=round(x,2)  
    if(is.null(breaks)) breaks=unique(quantile(x))
    x=cut(x,breaks,include.lowest=T,right=F)
    levels(x)=paste(breaks[-length(breaks)],ifelse(diff(breaks)>1,
                    c(paste('-',breaks[-c(1,length(breaks))]-1,sep=''),'+'),''),sep='')
    return(x)
  }
  fn.nmnl=function(dat,REL)
  {
    dat$cpue=with(dat,Catch.Target/(N.hooks.Fixed*SOAK.TIME)) 
    TT=with(subset(dat,cpue>0),table(year))
    dat=subset(dat,year%in%names(TT))
    
    d = dat %>%
      group_by(year) %>%
      summarise(n = length(cpue),
                m = length(cpue[cpue>0]),
                mean.lognz = mean(log(cpue[cpue>0])),
                sd.lognz = sd(log(cpue[cpue>0]))) %>%
      as.data.frame
    d$sd.lognz=ifelse(is.na(d$sd.lognz),0,d$sd.lognz)
    
    d=d%>%
      mutate(p.nz = m/n,
             theta = log(p.nz)+mean.lognz+sd.lognz^2/2,
             c = (1-p.nz)^(n-1),
             d = 1+(n-1)*p.nz,
             vartheta = ((d-c)*(1-c*d)-m*(1-c)^2)/(m*(1-c*d)^2)+
               sd.lognz^2/m+sd.lognz^4/(2*(m+1)),
             mean = exp(theta),
             lowCL = exp(theta - 1.96*sqrt(vartheta)),
             uppCL = exp(theta + 1.96*sqrt(vartheta))) 
    
    if(REL=="YES")
    {
      Mn=mean(d$mean,na.rm=T)
      d$mean=d$mean/Mn
      d$low95=d$lowCL/Mn
      d$up95=d$uppCL/Mn
    }
    all.yrs=seq(YEAR[1],d$year[length(d$year)])
    msn.yr=all.yrs[which(!all.yrs%in%d$year)]  
    msn.yr.dummy=d[1:length(msn.yr),]
    msn.yr.dummy[,]=NA
    msn.yr.dummy$year=msn.yr
    d=rbind(d,msn.yr.dummy)
    d=d[order(d$year),]
    return(d)
  }
  
  #1.11.1. Exploratory analysis 
  if(explore.cpue=="YES")
  {
    library(grid)
    # Explore 2004 high soak time and different set time
    a=subset(DATA,year==2004)
    a$Z=cfac(a$BOTDEPTH)
    a$LaT=cfac(a$Mid.Lat)
    b=table(a$Z,a$LaT)
    bb=aggregate(SOAK.TIME~Z+LaT,a,mean)
    bb
    table(a$Set.time)
    
    
    fun.prelim=function(DAT,NAMES)
    {
      par(mfcol=c(3,2),mai=c(.6,.65,.3,.1),oma=c(.2,.2,.1,.1),mgp=c(2.5, .5, 0))
      
      #positive catch
      TAb1.all=table(DAT$year)
      Pos=subset(DAT,Catch.Target>0)
      TAb1.pos=table(Pos$year)
      id=names(TAb1.all)[which(!names(TAb1.all)%in%names(TAb1.pos))]
      b=names(TAb1.pos)
      TAb1.pos=as.data.frame(matrix(TAb1.pos,ncol=length(TAb1.pos)))
      names(TAb1.pos)=b
      
      ss=names(TAb1.all)
      TAb1.all=as.data.frame(matrix(TAb1.all,ncol=length(TAb1.all)))
      names(TAb1.all)=ss
      
      if(length(id)>0)
      {
        a=as.data.frame(matrix(0,ncol=length(id)))
        names(a)=id
        TAb1.pos=cbind(TAb1.pos,a)
        TAb1.pos=TAb1.pos[,match(YEAR,names(TAb1.pos))]
      }
      TAB=rbind(TAb1.pos,TAb1.all)
      barplot(rbind(as.matrix(TAB[1,]),as.matrix(TAB[2,]-TAB[1,])))
      legend("top",legend = c("pos. catch","0 catch"), 
             fill = c("grey30", "grey80"),bty='n')
      box()
      
      plot(table(DAT$Catch.Target),type='h',xlab="Catch",ylab="Count",
           main="Zero infl. & right tail")
      
      boxcox(Catch.Target+0.0001 ~ as.factor(year), data = DAT,lambda = seq(0, 1, length = 10))
      legend("topright","Box Cox (should be small)",bty='n')
      
      # Cook distance to see outliers or overdisperse data (if distance >1)
      M1.1=glm(Catch.Target~as.factor(year)+offset(log(N.hooks.Fixed*SOAK.TIME)),family=poisson,data=DAT)
      plot(M1.1,which=4)
      legend("topright","outliers or overdispersion if distance >1",bty='n',cex=.85, text.col=2)
      
      #Outliers response var
      boxplot(DAT$Catch.Target~as.factor(DAT$year),main="Outliers in response var?",
              ylab="Numbers caught")
      
      plot(1:10,1:10,bty='n',pch='',xaxt='n',yaxt='n',ann=F)
      legend("center",NAMES,bty='n',cex=1.75)
    }
    pdf(paste(hndl.expl,"overdispersion.pdf",sep="/"))
    for (i in 1:N.species) fun.prelim(DATA.list[[i]],Tar.names[i])
    dev.off()
    
    #Number of positive records per species
    TAb.n.species=aggregate(Number~COMMON_NAME,subset(DATA,FixedStation=="YES" & SPECIES%in%TARGETS),sum)
    TAb.n.species=TAb.n.species[order(TAb.n.species$Number),]
    pdf(paste(hndl.expl,"Number.positive.records.per.species.pdf",sep="/"))
    grid.table(TAb.n.species)
    dev.off()
    
    PREDS=c("depth.bin","Month","Moon","BLOCK","SOI","Freo","Mid.Lat","Mid.Long")
    
    #Colinearity
    fn.colinearity=function(dat)
    {
      Covars=dat[,match(c("year","Month","BOTDEPTH","SOI","TEMP1","Freo","Mid.Lat","Mid.Long"),names(dat))]
      panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
      {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r <- abs(cor(x, y))
        txt <- format(c(r, 0.123456789), digits=digits)[1]
        txt <- paste(prefix, txt, sep="")
        if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
        text(0.5, 0.5, txt, cex = cex.cor * r)
      }
      pairs(Covars, lower.panel=panel.smooth, upper.panel=panel.cor)
      Variance.Inflation.Factor=corvif(Covars)
      
      # #catch correlation
      # par(mfcol=c(3,3),mai=c(.3,.4,.01,.1),oma=c(3,2,1,.1),las=1,mgp=c(1.5,.6,0))
      # for(j in 1:7)
      # {
      #   plot(dat[,2],dat[,j+2],ylab=TARGETS.name[j],xlab="")
      #   legend("topright",paste("correl.=",round(cor(dat[,2],dat[,j+2]),2)),bty='n')
      # }
      # mtext(TARGETS.name[2],1,1,outer=T)
       grid.newpage()
      grid.table(as.data.frame.matrix(Variance.Inflation.Factor))
    }
    pdf(paste(hndl.expl,"/colinearity.pdf",sep=""))
    List.explor=fn.colinearity(DATA.list[[1]]) 
    dev.off()
    
    
    #box plots of covariates (outliers in covariates?)
    fn.plt=function(a,y,TITL)
    {
      plot(1:nrow(a),ylim=c(0,max(a,na.rm=T)),col="transparent",ann=F,axes=F)
      CL=rainbow(ncol(a))
      for(pp in 1:ncol(a)) lines(1:nrow(a),a[,pp],col=CL[pp],lwd=4)
      axis(1,1:nrow(a),rownames(a))
      nn=seq(0,max(a,na.rm=T),length.out=5)
      axis(2,nn,round(nn))
      mtext(y,2,3,las=3,cex=1.5)
      legend("topright",colnames(a),text.col=CL,bty='n',title=TITL)
    }
    Yrs=length(unique(DATA.list[[1]]$year))
    div=1
    names(PREDS)=c(rep("Cat",4),rep("Cont",4))
    for (i in 1:N.species)
    {
      pdf(paste(hndl.expl,"/boxplot_",TARGETS.name[i],".pdf",sep=""))
      d=subset(DATA.list[[i]],FixedStation=="YES")
      #fn.plt(tapsum(d,"Catch.Target","year","Month",div=div),"Catch","month")
      par(mfcol=c(2,2),mar=c(2,2,2,.1))
      barplot(table(trunc(d$BOTDEPTH/2) * 2),main="2 m bin")
      barplot(table(trunc(d$BOTDEPTH/5) * 5),main="5 m bin")
      barplot(table(trunc(d$BOTDEPTH/10) * 10),main="10 m bin")
      barplot(table(trunc(d$BOTDEPTH/25) * 25),main="25 m bin")
      mtext("Depth categories",3,-2,outer=T,col=2)
      smart.par(length(PREDS),c(2.5,2.5,.1,.1),c(1,1,1,1),c(1.5,.5,0))
      for(pp in 1:length(PREDS))
      {
        x=d[,match(c("Catch.Target",PREDS[pp],"N.hooks.Fixed","SOAK.TIME","year"),names(d))]
        x$cpue=x$Catch.Target/(x$N.hooks.Fixed*x$SOAK.TIME)
        if(names(PREDS[pp])=="Cont") x[,2]=cut(x[,2],breaks=quantile(x[,2]))
        if(names(PREDS[pp])=="Cat") x[,2]=as.factor(x[,2])
        boxplot(x$cpue~x[,2],ylab="cpue",xlab=PREDS[pp],notch=F,varwidth=T)
      }
      dev.off()
    }

    #records with Positive catch and variation of cpue by factor level
    fn.explore.catch=function(dat)
    {
      names(dat)[match("Catch.Target",names(dat))]="Ktch"
      taB=aggregate(Ktch~year,dat,sum)
      grid.table(taB)
      
      for(pp in 1:length(PREDS))
      {
        x=dat[,match(c(PREDS[pp],"year","Ktch"),names(dat))]
        if(names(PREDS[pp])=="Cont") x[,1]=cut(x[,1],breaks=quantile(x[,1]))
        if(names(PREDS[pp])=="Cat") x[,1]=as.factor(x[,1])
        Form=as.formula(paste("Ktch~year",PREDS[pp],sep="+"))
        d=aggregate(Form,dat,sum)
        d.w=d%>%spread(key = PREDS[pp], value = Ktch)
        table=tableGrob(d.w)
        grid.newpage()
        h <- grobHeight(table)
        w <- grobWidth(table)
        title <- textGrob(PREDS[pp], y=unit(0.5,"npc") + 0.5*h, 
                          vjust=0, gp=gpar(fontsize=20))
         gt <- gTree(children=gList(table, title))
        grid.draw(gt)
      }
      
      # Main Effect of factors
      par(mfcol=c(1,1),mai=c(.3,1.1,.1,.1),oma=c(1,.1,.1,.1),las=1)
      dd=dat[,match(c("Ktch","N.hooks.Fixed","SOAK.TIME","year",PREDS),names(dat))]
      dd$year=as.factor(dd$year)
        for(pp in 1:length(PREDS))
        {
          id=match(PREDS[pp],names(dd))
          if(names(PREDS[pp])=="Cont") dd[,id]=cut(dd[,id],breaks=quantile(dd[,id]))
          if(names(PREDS[pp])=="Cat") dd[,id]=as.factor(dd[,id])
        }
      plot.design(Ktch/(N.hooks.Fixed*SOAK.TIME)~.,data=dd,cex.lab=1,
                  cex.axis=.75,ylab="cpue")
    }
    for (i in 1:N.species)
    {
      pdf(paste(hndl.expl,"/Pos_catch_",TARGETS.name[i],".pdf",sep=""))
      fn.explore.catch(dat=subset(DATA.list[[i]],FixedStation=="YES"))  
      dev.off()
    }
    
    #Correlation between factors
    pdf(paste(hndl.expl,"/Factor.correlation.pdf",sep=""))
    dd=subset(DATA.list[[1]],FixedStation=="YES")
    dd$year=as.factor(dd$year)
    dd$Moon=as.factor(dd$Moon)
    dd$Month=as.factor(dd$Month)
    
    ggplot(aes(y = Catch.Target, x = year, fill =Moon), data = dd) +
      geom_boxplot()
    grid.newpage()
    grid.table(with(dd,table(year,Moon)))
    ggplot(aes(y = Catch.Target, x = year, fill =Month), data = dd) +
      geom_boxplot()
    grid.newpage()
    grid.table(with(dd,table(year,Month)))
    
      # calculate Cramer's V varies from 0 to 1, with a 1 indicting a perfect association
    catcorrm <- function(vars, dat) sapply(vars, function(y) sapply(vars, function(x) assocstats(table(dat[,x], dat[,y]))$cramer))
    TAB=catcorrm(vars=c("year","Moon","Month"), dat=dd)
    grid.newpage()
    grid.table(TAB)
    dev.off()
    
    #no linear relation between predictors and response
    fn.expl.non.linear=function(d,VARS,NamE)
    {
      d=d %>% filter(BOTDEPTH<210 & FixedStation=="YES" & Catch.Target>0 ) %>%
        select(VARS) %>%
        mutate(Mid.Lat=abs(Mid.Lat),
               log.Ef=log(SOAK.TIME*N.hooks.Fixed),
               year=factor(year,levels=unique(sort(year))))
      par(mfcol=c(3,1),mar=c(2.5,2.2,1,.1),oma=c(1,1,1,1),mgp=c(1.5,.5,0),cex.lab=1.25)
      with(d,
           {
             boxplot(Catch.Target~as.factor(Month),col="grey80",xlab="Month")
             boxplot(Catch.Target~as.factor(10*round(BOTDEPTH/10)),col="grey80",xlab="Depth")
             boxplot(Catch.Target~as.factor(round(abs(Mid.Lat))),col="grey80",xlab="Lat")
           })
      mtext(NamE,3,outer=T,line=-1)
    }
    pdf(paste(hndl.expl,"Non.linear.pattern.pdf",sep="/"))
    for(i in Species.cpue)fn.expl.non.linear(d=DATA.list[[i]],VARS=c(VARIABLES,"Month"),NamE=names(DATA.list)[i])
    dev.off()
   }
  
  Species.cpue=1:N.species #define species vector
  Res.var="Catch.Target"   #define response var and offset
  Offset='offset(log.Ef)'
  
  #1.11.2. Select terms and error structure
  
  if(do.GAM=="YES" & Select.term=="YES")
  {
    LOgLik=function(MoD)  #get model Loglikelihood
    {
      if(!is.null(MoD))
      {
        if(class(MoD)[1]%in%c("zinbgam","zipgam")) LOgL=MoD$logL[length(MoD$logL)] else
          LOgL=as.numeric(logLik(MoD))
      }
      
      if(is.null(MoD))LOgL=NA
      return(LOgL)
    }
    fn.AIC=function(MoD)  #get AIC
    {
      return(MoD$aic)
    }
    fn.AICc=function(MoD,LoGLike)  #get AICc (not applicable if all models have same number of parameter)
    {
      if(!is.null(MoD))
      {
        Sample.size=dim(model.matrix(MoD))[1]
        Num.pars=dim(model.matrix(MoD))[2]
        AIC.c=-2*LoGLike+(2*Num.pars*(Sample.size/(Sample.size-Num.pars-1)))
      }
      return(AIC.c)
    }
    fn.AIC.ratio=function(DAT)
    {
      MIN=min(DAT,na.rm=T)
      Delta=DAT-MIN
      Like.model.give.dat=exp(-Delta/2)
      Weight=Like.model.give.dat/sum(Like.model.give.dat,na.rm=T)
      id=which(Weight==max(Weight,na.rm=T))
      names(Weight)=NULL
      Evidence.ratio=c(outer(Weight[id],Weight, "/"))
      return(list(Best.Mod=id,Delta=Delta,Like=Like.model.give.dat,Weight=Weight,Evidence.ratio_how.much.better=Evidence.ratio))
    }
    
    cl <- makeCluster(detectCores()-1)
    registerDoParallel(cl)
    getDoParWorkers()
    system.time({Select.term_error=foreach(i=Species.cpue,.packages=c('dplyr','doParallel')) %dopar%
      {
        if(names(DATA.list)[i]=="Sandbar shark")
        {
          Terms=c("year","Month","s(Mid.Lat,k=3)","s(BOTDEPTH,k=3)","s(Temp.res,k=3)","Moon")
        }else
        {
          Terms=c("year","s(Mid.Lat,k=3)","s(BOTDEPTH,k=3)","s(Temp.res,k=3)","Moon")
        }
        d=DATA.list[[i]] %>% filter(BOTDEPTH<210 & FixedStation=="YES") %>%
          mutate(Mid.Lat=abs(Mid.Lat),
                 log.Ef=log(SOAK.TIME*N.hooks.Fixed),
                 year=factor(year,levels=unique(sort(year))),
                 Month=factor(Month,levels=unique(sort(Month))),
                 Moon=factor(Moon,levels=c("Full","Waning","New","Waxing")))
        tested.mods=foreach(t =1:length(Terms),.errorhandling='remove',.packages=c('zigam','mgcv','doParallel'))%dopar%
          {
            FRMLA=formula(paste(Res.var,paste(c(Terms[1:t],Offset),collapse="+"),sep="~"))
            Gam.Pois=gam(FRMLA,data=d,method = "REML",family=poisson)
            Gam.Nb=gam(FRMLA,data=d,method = "REML",family = nb)
            Gam.ZIP=zipgam(lambda.formula=FRMLA,pi.formula=FRMLA,data=d)
            Gam.ZINB=zinbgam(mu.formula=FRMLA,pi.formula=FRMLA,data=d)
            dummy=list(Pois=Gam.Pois,NB=Gam.Nb,ZIP=Gam.ZIP,ZINB=Gam.ZINB)
            return(dummy)
          }
        for(t in 1:length(Terms))names(tested.mods)[t]=paste(Terms[1:t],collapse="+")
        return(tested.mods)
      }})    #takes 8 minutes
    names(Select.term_error)=names(DATA.list)
    stopCluster(cl) 
    
    #  Select best model terms (i.e. model structure) based on AIC and deviance explained
    Select.best.modl.str=function(d)
    {
      #set up table
      TeRms=names(d)
      Dev.Exp=as.data.frame(matrix(NA,nrow=length(TeRms),ncol=1+length(names(d[[1]]))))
      colnames(Dev.Exp)=c("Term",names(d[[1]]))
      for(i in 1:nrow(Dev.Exp))
      {
        n=which(strsplit(TeRms[i], "")[[1]]=="+")
        if(length(n)==0)
        {
          Nme=TeRms[i]
        }else
          Nme=substr(TeRms[i],n[length(n)]+1,1000)
        Dev.Exp$Term[i]=Nme
      }
      AiC=Dev.Exp
      
      #fill in table
      for(t in 1:nrow(Dev.Exp))
      {
        id=match(names(d[[t]]),names(Dev.Exp))
        for(x in 1:length(id))
        {
          if(class(d[[t]][[x]])[1]%in%c("zinbgam","zipgam"))
          {
            Dev.Exp[t,id[x]]=with(d[[t]][[x]][[1]],(null.deviance-deviance)/null.deviance)*100
            AiC[t,id[x]]=d[[t]][[x]][[1]]$aic
          }else
          {
            Dev.Exp[t,id[x]]=with(d[[t]][[x]],(null.deviance-deviance)/null.deviance)*100
            AiC[t,id[x]]=d[[t]][[x]]$aic
          }
        }
      }
      
      #Get formulas
      teRmS=rep(NA,length(d))
      for(t in 1:length(d)) teRmS[t]=paste(Res.var,paste(c(names(d)[t],Offset),collapse="+"),sep="~")
      
      #Get AIC derives
      AIC.der=vector('list',length=length(2:ncol(AiC)))
      names(AIC.der)=colnames(AiC)[2:ncol(AiC)]
      Best=AIC.der
      for(t in 2:ncol(AiC))
      {
        AIC.w.r=fn.AIC.ratio(AiC[,t])
        AIC.der[[t-1]]=data.frame(AIC=round(AiC[,t],1),
                     deltaAIC=round(AIC.w.r$Delta,1),
                     AIC.weight=formatC(AIC.w.r$Weight, format = "e", digits = 1))
        Best[[t-1]]=teRmS[AIC.w.r$Best.Mod]
      }
      AiC.tabl=cbind(Model=teRmS,do.call(cbind,AIC.der))
      return(list(AIC=AiC.tabl,Deviance.exp=Dev.Exp,Best=unlist(Best)))
    }
    SLCT.TRM=vector('list',length(Species.cpue))
    names(SLCT.TRM)=names(DATA.list)
    for(i in Species.cpue) SLCT.TRM[[i]]=Select.best.modl.str(d=Select.term_error[[i]])
    AIC.terms=do.call(rbind,lapply(SLCT.TRM, `[[`, 1))
    Resid.deviance.terms=do.call(rbind,lapply(SLCT.TRM, `[[`, 2))
    write.csv(AIC.terms,paste(getwd(),"/Model selection/AIC.terms.csv",sep=''))
    Resid.deviance.terms=Resid.deviance.terms%>%rename(Model=Term)
    Resid.deviance.terms$Model=AIC.terms$Model
    write.csv(Resid.deviance.terms,paste(getwd(),"/Model selection/Deviance.explained.by.terms.csv",sep=''))
    Best.SLCT.TRM=SLCT.TRM  
    for(i in Species.cpue)
    {
      d=Best.SLCT.TRM[[i]]$Deviance.exp
      Best.terms=vector('list',ncol(d)-1)
      names(Best.terms)=colnames(d)[-1]
      for(p in 1:length(Best.terms))
      {
        d1=d[,c(1,p+1)]%>%mutate(Prev=lag(d[,p+1],1),
                                 Exp.per=d[,p+1]-Prev,
                                 Exp.per=ifelse(is.na(Exp.per),10,Exp.per))
        Best.terms[[p]]=d1$Term[which(d1$Exp.per>1)]
      }
      Best.SLCT.TRM[[i]]=Best.terms
    }
    
    #  Select best error 
    cl <- makeCluster(detectCores()-1)
    registerDoParallel(cl)
    getDoParWorkers()
    system.time({Select.error=foreach(i=Species.cpue,.packages=c('zigam','mgcv','dplyr','doParallel')) %dopar%
      {
        Terms=Best.SLCT.TRM[[i]]
        
        TT=with(subset(DATA.list[[i]],Catch.Target>0),table(year))
        d=subset(DATA.list[[i]],year%in%names(TT[TT>0]))
        
        d=d %>% filter(BOTDEPTH<210 & FixedStation=="YES") %>%
                mutate(Mid.Lat=abs(Mid.Lat),
                       log.Ef=log(SOAK.TIME*N.hooks.Fixed),
                       year=factor(year,levels=unique(sort(year))),
                       Month=factor(Month,levels=unique(sort(Month))),
                       Moon=factor(Moon,levels=c("Full","Waning","New","Waxing")))
        
      
        Gam.Pois=gam(formula(paste(Res.var,paste(c(Terms$Pois,Offset),collapse="+"),sep="~")),
                     data=d,method = "REML",family=poisson)
        Gam.Nb=gam(formula(paste(Res.var,paste(c(Terms$NB,Offset),collapse="+"),sep="~")),
                   data=d,method = "REML",family = nb)
        Gam.ZIP=zipgam(lambda.formula=formula(paste(Res.var,paste(c(Terms$ZIP,Offset),collapse="+"),sep="~")),
                       pi.formula=formula(paste(Res.var,paste(c(Terms$ZIP,Offset),collapse="+"),sep="~")),
                       data=d)
        Gam.ZINB=zinbgam(mu.formula=formula(paste(Res.var,paste(c(Terms$ZINB,Offset),collapse="+"),sep="~")),
                         pi.formula=formula(paste(Res.var,paste(c(Terms$ZINB,Offset),collapse="+"),sep="~")),
                         data=d)
        
        dummy=list(Pois=Gam.Pois,NB=Gam.Nb,ZIP=Gam.ZIP,ZINB=Gam.ZINB)
        return(dummy)
      }})    #takes 20 seconds
    names(Select.error)=names(DATA.list)
    stopCluster(cl) 
    
    #output pdf  #ACA
    CL=c("grey50","grey80","blue","cornflowerblue")
    fn.pred=function(term,DAT,MOD,Nml)
    {
      TT=with(subset(DAT,Catch.Target>0),table(year))
      DAT=subset(DAT,year%in%names(TT[TT>0]))
      
      newdata=data.frame(year=factor(sort(unique(DAT$year)),levels=unique(sort(DAT$year))),
                         Month=mx.lev(DAT$Month),
                         Mid.Lat=mean(DAT$Mid.Lat),
                         BOTDEPTH=mean(DAT$BOTDEPTH),
                         Temp.res=mean(DAT$Temp.res),
                         Moon=mx.lev(DAT$Moon),
                         log.Ef=log(mean(DAT$SOAK.TIME*DAT$N.hooks.Fixed))) 
      
      PRD=vector('list',length(MOD))
      for(m in 1:length(MOD))
      {
        d1=data.frame(year=newdata$year)
        pp=c(predict(MOD[[m]],type='response',newdata=newdata))
        d1$year=as.numeric(as.character(d1$year))
        d1$mean=pp/mean(pp)
        misn=which(!Nml$year%in%d1$year)
        if(length(misn)>0)
        {
          d1=rbind(d1,Nml[misn,c('year','mean')])
          d1=d1[order(d1$year),]
        }
        PRD[[m]]=d1
      }
      
      plot(Nml$year,Nml$mean,pch=19,cex=2,ylab="Relative cpue",xlab="year",
           ylim=c(0,max(Nml$mean,na.rm=T)))
      for(m in 1:length(MOD)) 
      {
        lines(PRD[[m]]$year,PRD[[m]]$mean,col=CL[m],lwd=2)
        points(PRD[[m]]$year,PRD[[m]]$mean,col=CL[m],pch=19,cex=1.5)
      }
        
      legend('topleft',paste(c("nominal",names(MOD))),col=c('black',CL),bty='n',
             pch=c(19,rep(NA,length(MOD))),lwd=c(1,rep(2,length(MOD))),lty=c(NA,rep(1,length(MOD))))
      
    }
    mx.lev=function(x) factor(names(which(table(x)==max(table(x)))) ,levels=levels(x))
    Gam.fit.diag=function(MOD,SP)
    {
      par(mfcol=c(2,2))
      for(m in 1:length(MOD))
      {
        Nm=ifelse(names(MOD)[m]=="Pois","Poisson",
                  ifelse(names(MOD)[m]=="NB","Negative binomial",
                         ifelse(names(MOD)[m]=="ZIP","Zero-inflated Poisson",
                                "Zero-inflated negative binomial")))
        if(class(MOD[[m]])[1]%in%c("zipgam","zinbgam"))
        {
          gam.check(MOD[[m]][[1]])
          mtext(paste("Figure ",i,".",m,". ",SP," fit diagnostics. ",Nm,
                      " distribution, counts part",sep=''),3,outer=T,line=-1.5,cex=.8)
          gam.check(MOD[[m]][[2]])
          mtext(paste("Figure ",i,".",m,". ",SP," fit diagnostics. ",Nm,
                      " distribution, binomial part",sep=''),3,outer=T,line=-1.5,cex=.8)
        }else
        {
          gam.check(MOD[[m]])
          mtext(paste("Figure ",i,".",m,". ",SP," fit diagnostics. ",Nm,
                      " distribution",sep=''),3,outer=T,line=-1.5,cex=1)
        }
      }
    }
    
    pdf(paste(hndLs,"/Model selection/Compare error.pdf",sep=""))
    
    #table of deviance explained for each term
    # plot.new()
    # mytheme <- gridExtra::ttheme_default(
    #   core = list(padding=unit(c(1, 1), "mm"),fg_params=list(cex = .75)),
    #   colhead = list(fg_params=list(cex = .8)),
    #   rowhead = list(fg_params=list(cex = .8)))
    # myt <- gridExtra::tableGrob(Resid.deviance.terms[1:21,], theme = mytheme)
    # grid.draw(myt)
    # mtext("Deviance explained",3,line=2,cex=1.5)
    # 
    # plot.new()
    # mytheme <- gridExtra::ttheme_default(
    #   core = list(padding=unit(c(1, 1), "mm"),fg_params=list(cex = .75)),
    #   colhead = list(fg_params=list(cex = .8)),
    #   rowhead = list(fg_params=list(cex = .8)))
    # myt <- gridExtra::tableGrob(Resid.deviance.terms[22:41,], theme = mytheme)
    # grid.draw(myt)

    for(i in Species.cpue)
    {
        #Table of information criteria
      Title=c("Log.Like","AIC.c","AIC.delta","AIC.w","AIC.Ev.ratio")
      Compare=as.data.frame(matrix(ncol=length(Title),nrow=length(Select.error[[i]])))
      colnames(Compare)=Title
      rownames(Compare)=names(Select.error[[i]])
      Compare$Log.Like=round(unlist(sapply(Select.error[[i]],LOgLik)))
      store.aicc=unlist(sapply(Select.error[[i]],fn.AIC))
      Compare$AIC.c=round(store.aicc)
      AIC.w.r=fn.AIC.ratio(store.aicc)
      Compare$AIC.delta=round(AIC.w.r$Delta)
      Compare$AIC.w=formatC(AIC.w.r$Weight, format = "e", digits = 2)
      Compare$AIC.Ev.ratio=formatC(AIC.w.r$Evidence.ratio_how.much.better, format = "e", digits = 2)
      id.best=which(AIC.w.r$Weight==max(AIC.w.r$Weight,na.rm=T))
      Best.error.AIC.w=rownames(Compare)[id.best]
      Best.formula=formula(paste(Res.var,paste(c(Terms[id.best],Offset),collapse="+"),sep="~"))
      Compare=Compare%>%rename('Log likelihood'=Log.Like,
                               AIC=AIC.c,
                               'delta AIC'=AIC.delta,
                               'AIC weight'=AIC.w,
                               'Evidence ratio'=AIC.Ev.ratio)
      plot.new()
      mytheme <- gridExtra::ttheme_default(
        core = list(padding=unit(c(1, 1), "mm"),fg_params=list(cex = 1)),
        colhead = list(fg_params=list(cex = 1.25)),
        rowhead = list(fg_params=list(cex = 1.25)))
      myt <- gridExtra::tableGrob(Compare, theme = mytheme)
      grid.draw(myt)
    
      mtext(paste('Table ',i,'. Summary of error selection for the best model structure of ',
                  Tar.names[i],sep=""),3,cex=.9,adj=0)
      mtext(c(" (Pois, Poisson; NB, Negative binomial;"),3,-1,cex=.9,adj=0)
      mtext(c("ZIP, zero-inflated Poisson; ZINB, zero-inflated negative binomial)"),3,-2,cex=.9,adj=0)
      
      #model fit to data
      Gam.fit.diag(MOD=Select.error[[i]],SP=names(Select.error)[i])
      
      par(mfcol=c(1,1))
      fn.pred(term="year",
              DAT=DATA.list[[i]] %>% filter(BOTDEPTH<210 & FixedStation=="YES") %>%
                mutate(Mid.Lat=abs(Mid.Lat),
                       log.Ef=log(SOAK.TIME*N.hooks.Fixed),
                       Month=factor(Month,levels=unique(sort(Month))),
                       Moon=factor(Moon,levels=c("Full","Waning","New","Waxing"))),
              MOD=Select.error[[i]],
              Nml=fn.nmnl(dat=subset(DATA.list[[i]],FixedStation=="YES"),REL="YES"))
      mtext(paste("Figure ",i,".",5,". ",names(Select.error)[i]," Predicted annual relative cpue",
                  sep=''),3,line=2.5,cex=1,adj=0)
      mtext(paste("(Pois, Poisson; NB, negative binomial;",sep=''),3,line=1.5,cex=1,adj=0)
      mtext(paste("ZIP, zero-inflated Poisson; ZINB, zero-inflated negative binomial)",sep=''),3,
            line=0.5,cex=1,adj=0)
      
    }
    dev.off()
    
  }
  
  if(do.GLM=="YES" & Select.term=="YES")
  {
    # Set up different distribution function
    #Poisson and NB
    Count.Pois.fn=function(DAT,FORMULA) Pois=glm(FORMULA,data=DAT,family="poisson")
    Count.NB.fn=function(DAT,FORMULA) NB=glm.nb(FORMULA, data = DAT)
    
    #Zero Inflated
    ZIP.fn=function(DAT,FORMULA)  ZIP <- zeroinfl(FORMULA, data = DAT, dist = "poisson")
    ZINB.fn=function(DAT,FORMULA) if(!is.null(FORMULA))ZINB <- zeroinfl(FORMULA, data = DAT, dist = "negbin")
    
    #Tweedie
    Tweedie.fn=function(DAT,FORMULA)
    {
      #Tweedie
      # 1.Find p value thru likelihood profiling
      #note: expect  xi  about 1.5
      P.eval=seq(1.1, 1.9, length=9)
      out <- tweedie.profile(FORMULA,xi.vec=P.eval,data=DAT, do.plot=F)
      P=out$xi.max
      
      # 2.Fit the glm
      Tweed=glm( FORMULA,data=DAT, family=tweedie(var.power=P, link.power=0) )
      
      return(list(Tweed=Tweed))
    }
    
    #Hurdle
    Hurdle.Pois.fn=function(DAT,FORMULA) Hurdle.Pois = hurdle(FORMULA, data = DAT , dist="poisson")
    Hurdle.NB.fn=function(DAT,FORMULA) if(!is.null(FORMULA))Hurdle.NB = hurdle(FORMULA, data = DAT , dist="negbin")
    
    #Fit all model structures to all species
    Preds=c("year","Mid.Lat","BOTDEPTH")
    VARIABLES=c("Catch.Target",Preds,"N.hooks.Fixed","SOAK.TIME")
    names(Preds)=c("Cat","Cont","Cont")  
    
    # Select terms for each error structure
    if(Select.term=="YES")
    {
      fitFun=function(FORMULA, data,...) zeroinfl(FORMULA, data,dist = ErroR)
      choose.term=function(DAT,FORMULA,ErroR)
      {
        id=match(Preds[which(names(Preds)=="Cat")],names(DAT))
        for(pp in 1:length(id))DAT[,id[pp]]=as.factor(DAT[,id[pp]])
        DAT$log.Effort=log(DAT$N.hooks.Fixed*DAT$SOAK.TIME)
        DAT$Mid.Lat=abs(DAT$Mid.Lat)
        #m1 = zeroinfl(FORMULA, data = DAT, dist = ErroR)
        #fitbe = be.zeroinfl(m1, data = DAT, dist = ErroR,alpha = 0.01, trace = FALSE)
        res <- glmulti(FORMULA,data=DAT, level=1,method="h", fitfunction=fitFun,
                       crit="aicc",plotty = F, report = T)
        BEST=res@formulas[[1]]
        
        return(BEST)
      }
      choose.term(DAT=subset(DATA.list[[i]],FixedStation=="YES",select=VARIABLES),
                  FORMULA=as.formula(paste("Catch.Target",paste(c(Preds,'offset(log.Effort)'),collapse="+"),sep="~")),
                  ErroR="poisson")
      
      #Formulas count
      Formulas=list(
        as.formula(Catch.Target ~ year + offset(log.Effort)),
        as.formula(Catch.Target ~ year + Month + offset(log.Effort)),
        as.formula(Catch.Target ~ year + Month + BOTDEPTH + offset(log.Effort)),
        as.formula(Catch.Target ~ year + Month + BOTDEPTH + Moon + offset(log.Effort)),
        as.formula(Catch.Target ~ year + Month + BOTDEPTH + Moon + Mid.Lat +  offset(log.Effort)))
      
      #Formulas zero inflated
      Formulas.ZI=list(
        as.formula(Catch.Target ~ year + offset(log.Effort)|year + offset(log.Effort)),
        as.formula(Catch.Target ~ year + Month + offset(log.Effort)|year + Month + offset(log.Effort)),
        as.formula(Catch.Target ~ year + Month + BOTDEPTH + offset(log.Effort)|year + Month + BOTDEPTH + offset(log.Effort)),
        as.formula(Catch.Target ~ year + Month + BOTDEPTH + Moon + offset(log.Effort)|year + Month + BOTDEPTH + Moon+ offset(log.Effort)),
        as.formula(Catch.Target ~ year + Month + BOTDEPTH + Moon + Mid.Lat +  offset(log.Effort)|year + Month + BOTDEPTH + Moon + Mid.Lat+ offset(log.Effort)))
      
      # Formulas.ZI.Du=Formulas.ZI.SH=NULL
      # Formulas.ZI.SB=list(
      #   as.formula(Catch.Target ~ year + BOTDEPTH + offset(log.Effort)|year + offset(log.Effort)),
      #   as.formula(Catch.Target ~ year + BOTDEPTH + Moon+ offset(log.Effort)|year + offset(log.Effort)),
      #   as.formula(Catch.Target ~ year + BOTDEPTH + Moon+ Mid.Lat+ offset(log.Effort)|year + Mid.Lat+ offset(log.Effort)))
      # 
      # Formulas.ZI.Ti=list(
      #   as.formula(Catch.Target ~ year + Moon+ offset(log.Effort)|year + offset(log.Effort)),
      #   as.formula(Catch.Target ~ year + Moon+ Mid.Lat+ offset(log.Effort)|year + offset(log.Effort)),
      #   as.formula(Catch.Target ~ year + Moon+ Mid.Lat+ BOTDEPTH+ offset(log.Effort)|year + offset(log.Effort)))
      # 
      # Formulas.ZI.BT=list(
      #   as.formula(Catch.Target ~ year + BOTDEPTH + Moon+ offset(log.Effort)|year + Moon+ offset(log.Effort)),
      #   as.formula(Catch.Target ~ year + BOTDEPTH + Moon+ Mid.Lat+ offset(log.Effort)|year + Moon+ offset(log.Effort)))
      # 
      # Formulas.ZI.ST=list(
      #   as.formula(Catch.Target ~ year + BOTDEPTH + offset(log.Effort)|year + offset(log.Effort)),
      #   as.formula(Catch.Target ~ year + BOTDEPTH + Moon+ Mid.Lat+ offset(log.Effort)|year + Mid.Lat+ offset(log.Effort)))
      # 
      # Formulas.ZI.Mi=Formulas.ZI.SB
      # 
      # Formulas.ZI=list(Formulas.ZI.Du,Formulas.ZI.SB,Formulas.ZI.Ti,Formulas.ZI.SH,
      #                  Formulas.ZI.BT,Formulas.ZI.ST,Formulas.ZI.Mi)
      
      
      #Formulas Tweedie
      Formulas.Twe=Formulas
      
      #Formulas Hurdle
      Formulas.Hu=Formulas.ZI
      # Formulas.Hu.Du=Formulas.Hu.SH=NULL
      # Formulas.Hu.SB=Formulas
      # Formulas.Hu.Ti=list(
      #   as.formula(Catch.Target ~ BOTDEPTH + offset(log.Effort) | year + BOTDEPTH+ offset(log.Effort)),
      #   as.formula(Catch.Target ~ BOTDEPTH +Moon+ offset(log.Effort) | year + BOTDEPTH+Moon+ offset(log.Effort)),
      #   as.formula(Catch.Target ~ BOTDEPTH +Moon+Mid.Lat+SOI+ offset(log.Effort) | year + BOTDEPTH+Moon+Mid.Lat+SOI+ offset(log.Effort)))
      # 
      # Formulas.Hu.BT=Formulas 
      # Formulas.Hu.ST=Formulas 
      # Formulas.Hu.Mi=Formulas 
      # 
      # Formulas.Hu=list(Formulas.Hu.Du,Formulas.Hu.SB,Formulas.Hu.Ti,Formulas.Hu.SH,
      #                  Formulas.Hu.BT,Formulas.Hu.ST,Formulas.Hu.Mi)
      
      
      
      Store.count=vector('list',N.species)
      names(Store.count)=names(DATA.list)
      Store.count.ZI=Store.count.Twee=Store.count.Hur=Store.count
      system.time(for(i in Species.cpue)
      {
        Dat=DATA.list[[i]]
        
        #count glms
        S.count=vector('list',length(Formulas))
        for(p in 1:length(Formulas)) S.count[[p]]=Count.fn(Dat,Formulas[[p]])
        Store.count[[i]]=S.count
        
        #Zero inflated glms
        Form.ZI=Formulas.ZI
        #Form.ZI=Formulas.ZI[[i]]
        Z.count=vector('list',length(Form.ZI))
        for(p in 1:length(Form.ZI)) Z.count[[p]]=ZI.fn(Dat,Form.ZI[[p]])
        Store.count.ZI[[i]]=Z.count
        
        #Tweedie glms
        Tw.count=vector('list',length(Formulas.Twe))
        for(p in 1:length(Formulas.Twe)) Tw.count[[p]]=Tweedie.fn(Dat,Formulas.Twe[[p]])
        Store.count.Twee[[i]]=Tw.count
        
        #Hurdle glms
        Form.Hu=Formulas.Hu[[i]]
        H.count=vector('list',length(Form.Hu))
        for(p in 1:length(Form.Hu)) H.count[[p]]=Hurdle.fn(Dat,Form.Hu[[p]])
        Store.count.Hur[[i]]=H.count
      })
    }
    
    #Predicting function
    Pred.Model=function(MOD,Dta)
    {
      Dta$log.Effort=log(Dta$N.hooks.Fixed*Dta$SOAK.TIME)
      PREDS=NULL
      if(!is.null(MOD))PREDS=predict(MOD,Dta, type='response')
      return(list(PREDS=PREDS))
    }
    
    # Select best error structure   
    if(Select.error=="YES")
    {
      #Formulas
      #count
      Formulas=as.formula(Catch.Target ~ year+log.Mid.Lat+log.BOTDEPTH+offset(log.Effort))
      
      #zero inflated and Hurdle
      Formulas.ZI=Formulas.Hu=DATA.list
      for(i in Species.cpue)
      {
        ll=list(P=Formulas,NB=Formulas)
        Formulas.ZI[[i]]=ll
        Formulas.Hu[[i]]=ll
      }
      
      #Specify submodels for some species
      Formulas.ZI$"Sandbar shark"$NB=
        as.formula(Catch.Target ~ year + log.Mid.Lat + offset(log.Effort) |
                     year + log.Mid.Lat + log.BOTDEPTH + offset(log.Effort))
      
      Formulas.ZI$"Spot-tail shark"$NB=
        as.formula(Catch.Target ~ year + log.Mid.Lat + log.BOTDEPTH + offset(log.Effort) | 
                     year + offset(log.Effort))
      Formulas.Hu$"Spot-tail shark"$P=Formulas.ZI$"Sandbar shark"$NB
      
      Formulas.ZI$"Tiger shark"$NB=
        as.formula(Catch.Target ~ year + log.Mid.Lat + offset(log.Effort) | 
                     year +log.Mid.Lat + log.BOTDEPTH + offset(log.Effort))
      
      Formulas.Hu$"Tiger shark"$P=Formulas.Hu$"Tiger shark"$NB=
        as.formula(Catch.Target ~ log.Mid.Lat + log.BOTDEPTH + offset(log.Effort) | 
                     year + log.Mid.Lat + log.BOTDEPTH + offset(log.Effort))
      
      Formulas.ZI$"Blacktip sharks"$NB=
        as.formula(Catch.Target ~ year + log.BOTDEPTH + offset(log.Effort) | 
                     year + log.Mid.Lat + log.BOTDEPTH + offset(log.Effort))
      Formulas.Hu$"Blacktip sharks"$NB=
        as.formula(Catch.Target ~ log.Mid.Lat + log.BOTDEPTH + offset(log.Effort) | 
                     year + log.Mid.Lat + log.BOTDEPTH + offset(log.Effort))
      
      Formulas.Hu$"Scalloped hammerhead"$P=
        as.formula(Catch.Target ~ year + log.BOTDEPTH + offset(log.Effort) | 
                     year + log.Mid.Lat + log.BOTDEPTH + offset(log.Effort))
      Formulas.Hu$"Scalloped hammerhead"$NB=
        as.formula(Catch.Target ~ log.Mid.Lat + log.BOTDEPTH + offset(log.Effort) | 
                     year + log.Mid.Lat + log.BOTDEPTH + offset(log.Effort))
      
      Formulas.Hu$"Dusky shark"$P=
        as.formula(Catch.Target ~ year + log.BOTDEPTH + offset(log.Effort) | 
                     year + log.Mid.Lat + log.BOTDEPTH + offset(log.Effort))
      Formulas.ZI$"Dusky shark"$NB=
        as.formula(Catch.Target ~ year + offset(log.Effort) | 
                     year + log.Mid.Lat + log.BOTDEPTH + offset(log.Effort))
      
      Formulas.ZI$"Dusky shark"$P=Formulas.Hu$"Dusky shark"$NB=
        as.formula(Catch.Target ~ log.Mid.Lat + log.BOTDEPTH + offset(log.Effort) | 
                     year + log.Mid.Lat + log.BOTDEPTH + offset(log.Effort))
      
      Formulas.ZI$"Sliteye shark"$P=
        as.formula(Catch.Target ~ year + log.Mid.Lat + offset(log.Effort) | 
                     year + log.Mid.Lat + log.BOTDEPTH + offset(log.Effort))
      Formulas.ZI$"Sliteye shark"$NB=Formulas.ZI$"Dusky shark"$NB
      Formulas.Hu$"Sliteye shark"$P=Formulas.Hu$"Sliteye shark"$NB=
        as.formula(Catch.Target ~ log.Mid.Lat + log.BOTDEPTH + offset(log.Effort) | 
                     year + log.Mid.Lat + log.BOTDEPTH + offset(log.Effort))
      
      #Fit all error structures
      Store=vector('list',N.species)
      names(Store)=names(DATA.list)
      
      fn.chk.error=function(d)
      {
        TT=with(subset(d,Catch.Target>0),table(year))
        d=subset(d,year%in%names(TT))
        d$year=factor(d$year,levels=sort(unique(d$year)))
        d$log.Effort=log(d$N.hooks.Fixed*d$SOAK.TIME)
        d$Mid.Lat=abs(d$Mid.Lat)
        d$log.BOTDEPTH=log(d$BOTDEPTH)
        d$log.Mid.Lat=log(d$Mid.Lat)
        
        Pois=NB=ZINB=ZIP=Hurdle.Pois=Hurdle.NB=NULL
        
        Pois=Count.Pois.fn(DAT=d,FORMULA=Formulas)
        NB=Count.NB.fn(DAT=d,FORMULA=Formulas)
        
        tryCatch({
          ZIP=ZIP.fn(DAT=d,FORMULA=Formulas.ZI[[i]]$P)
        }, error = function(e) 
        {
        })
        tryCatch({
          ZINB=ZINB.fn(DAT=d,FORMULA=Formulas.ZI[[i]]$NB)
        }, error = function(e) 
        {
        })
        
        tryCatch({
          Hurdle.Pois=Hurdle.Pois.fn(DAT=d,FORMULA=Formulas.Hu[[i]]$P)
        }, error = function(e) 
        {
        })
        
        tryCatch({
          Hurdle.NB=Hurdle.NB.fn(DAT=d,FORMULA=Formulas.Hu[[i]]$NB)
        }, error = function(e) 
        {
        })
        
        return=list(MODS=list(Pois=Pois,NB=NB,ZIP=ZIP,ZINB=ZINB,
                              Hurdle.Pois=Hurdle.Pois,Hurdle.NB=Hurdle.NB),
                    DAT=d,
                    Forms=list(Pois=Formulas,NB=Formulas,
                               ZIP=Formulas.ZI[[i]]$P,ZINB=Formulas.ZI[[i]]$NB,
                               Hurdle.Pois=Formulas.Hu[[i]]$P,Hurdle.NB=Formulas.Hu[[i]]$NB))
      }
      system.time(for(i in Species.cpue)  #takes 13 seconds
      {
        Store[[i]]=fn.chk.error(d=subset(DATA.list[[i]],FixedStation=="YES",select=VARIABLES)) 
      })
      
      #output pdf
      kmpr.error.pred=function(MOD)
      {
        lsm=MOD
        for(m in 1:length(lsm))
        {
          if(!is.null(MOD[[m]]))
          {
            lsm[[m]]=summary(emmeans(MOD[[m]],'year', type="response"))
            names(lsm[[m]])=c("year","mean","SE","df","Low.CI","UP.CI")
          }
        }
        Yrs=as.numeric(as.character(levels(DAT$year)))
        COL=c("black","grey50","blue","cyan","grey80","cornflowerblue")
        Dlt=seq(-0.2,.2,length.out = 6)
        
        plot(Yrs,Yrs,main="",xlab="",ylab="CPUE",ylim=c(0,min(c(10,max(do.call(rbind,lsm)$UP.CI,na.rm = T)))))
        for(m in 1:length(lsm))
        {
          if(!is.null(MOD[[m]]))
          {
            with(lsm[[m]],points(Yrs+Dlt[m],mean,cex=1.25,"o", pch=16, lty=2,col=COL[m]))
            suppressWarnings(with(lsm[[m]],arrows(x0=Yrs+Dlt[m], y0=Low.CI, x1=Yrs+Dlt[m],
                                                  y1=UP.CI,code = 3,angle=90,length=.025,col=COL[m])))
            
          }
        }
        not.converged=names(MOD)[unlist(lapply(MOD,is.null))]
        ids=match(not.converged,names(MOD))
        LGn=c("Pois","NB","ZIP","ZINB","HuP","HuNB")
        LGn[ids]=paste(LGn[ids],"(no convergence)")
        legend("top",LGn,col=COL,pch=19,bty='n')
      }
      kmpr.error.BIC=function(MOD)
      {
        id=which(!names(MOD)%in%names(MOD)[unlist(lapply(MOD,is.null))])
        par(las=3)
        barplot(unlist(lapply(MOD[id],BIC)),ylab="BIC")
      }
      
      kmpr.error.fit.diag=function(MOD)
      {
        #Show coefficients
        for(m in 1:length(MOD))
        {
          if(!is.null(MOD[[m]]))
          {
            plot.new()
            legend('topleft',names(MOD)[m],bty='n')
            tmp <- round(as.data.frame.matrix(coeftest(MOD[[m]])),3)
            mytheme <- gridExtra::ttheme_default(
              core = list(padding=unit(c(1, 1), "mm"),fg_params=list(cex = .5)),
              colhead = list(fg_params=list(cex = .7)),
              rowhead = list(fg_params=list(cex = .7)))
            myt <- gridExtra::tableGrob(tmp, theme = mytheme)
            grid.draw(myt)
          }
        }
        
        #Quantiles
        par(mfcol=c(3,2),las=1,mar=c(3,3,1,1),mgp=c(1.5,.7,0))
        for(m in 1:length(MOD))
        {
          if(!is.null(MOD[[m]])) qqrplot(MOD[[m]], main = names(MOD)[m])
        }
        
        #Rotogram
        #note: if wave like pattern, then overdispersion in data
        par(mfcol=c(3,2),las=1,mar=c(3,3,1,1),mgp=c(1.5,.7,0))
        for(m in 1:length(MOD))
        {
          if(!is.null(MOD[[m]])) rootogram(MOD[[m]], main = names(MOD)[m])
        }
        
        
        #Pearson Residuals (more appropriate for count distributions)
        par(mfcol=c(3,2),las=1,mar=c(3,3,1,1),mgp=c(1.5,.7,0))
        for(m in 1:length(MOD))
        {
          if(!is.null(MOD[[m]]))
          {
            E1=resid(MOD[[m]],type='pearson')
            hist(E1,xlab="Pearson residuals",main=names(MOD[m]))
          }
          if(is.null(MOD[[m]]))
          {
            plot(1,col="transparent",ann=F,axes=F)
            legend("center",paste(names(MOD[m]),"(no convergence)"),bty='n')
          }
        }
      }
      
      OverDisp=function(MoD)  #If Overdispersion >> 1, then Poisson...
      {
        E1=resid(MoD,type='pearson')
        Dispersion=sum(E1^2)/MoD$df.resid 
        if(is.null(MoD))Dispersion=NA
        return(list(Dispersion=Dispersion))
      }
      
      LOgLik=function(MoD)  #get model Loglikelihood
      {
        if(!is.null(MoD))LOgL=logLik(MoD)
        if(is.null(MoD))LOgL=NA
        return(LOgL)
      }
      
      fn.AIC=function(MoD)  #get AIC
      {
        AIc=NA
        if(!is.null(MoD)) AIc=AIC(MoD)
        return(AIc)
      }
      fn.AICc=function(MoD)  #get AICc
      {
        AIC.c=NA
        if(!is.null(MoD))
        {
          Sample.size=dim(model.matrix(MoD))[1]
          Num.pars=dim(model.matrix(MoD))[2]
          AIC.c=AIC(MoD)+((2*Num.pars^2+2*Num.pars)/(Sample.size+Num.pars+1))
        }
        #if(names(DAT)=="Tweed" & !is.null(DAT[[1]])) AIC.c=AICtweedie(DAT[[1]])
        return(AIC.c)
      }
      fn.AIC.ratio=function(DAT)
      {
        MIN=min(DAT,na.rm=T)
        Delta=DAT-MIN
        Like.model.give.dat=exp(-Delta/2)
        Weight=Like.model.give.dat/sum(Like.model.give.dat,na.rm=T)
        id=which(Weight==max(Weight,na.rm=T))
        names(Weight)=NULL
        Evidence.ratio=outer(Weight[id],Weight, "/")
        return(list(Best.Mod=id,Delta=Delta,Like=Like.model.give.dat,Weight=Weight,Evidence.ratio_how.much.better=Evidence.ratio))
      }
      
      Deg.Free=function(MoD)   #get degrees of freedom
      {
        if(is.null(MoD$df.residual))return(NA)
        if(!is.null(MoD$df.residual)) return(MoD$df.residual)  
      }
      
      Get.dev=function(MoD,DAT)
      {
        Null.2times.like=Res.2times.like=NA
        if(!is.null(MoD))
        {
          MOD.null=update(MoD,.~1)
          Null.2times.like=as.numeric(-2*logLik(MOD.null))
          Res.2times.like=as.numeric(-2*logLik(MoD)) 
        }
        return(list(Null.2times.like=Null.2times.like,Res.2times.like=Res.2times.like))
      }
      
      Dsquared <- function(model, adjust = F)
      {
        d1=NA
        d2=NA
        d3=NA
        
        if(!is.null(model$deviance))
        {
          d2 <- (model$null.deviance - model$deviance) / model$null.deviance
          if (adjust)
          {
            n <- length(model$fitted.values)
            p <- length(model$coefficients)
            d2 <- 1 - ((n - 1) / (n - p)) * (1 - d2)
          }
          d3=d2*100
          d1=model$deviance
          d2=model$null.deviance
        }
        
        return(list(d1=d1,d2=d2,d3=d3))
      }
      
      Inc.Dev.Exp=function(A,B)round((100*A/B)-100)  #increase in deviance explained from one model to next
      
      for(i in Species.cpue)
      {
        pdf(paste(hndLs,"/Model selection/Comapre error.",Tar.names[i],".pdf",sep=""))
        DAT=Store[[i]]$DAT
        kmpr.error.pred(MOD=Store[[i]]$MODS)
        rm(DAT)
        kmpr.error.BIC(MOD=Store[[i]]$MODS)
        kmpr.error.fit.diag(MOD=Store[[i]]$MODS)
        
        
        #Put model comparisons in table
        MoDs=Store[[i]]$MODS
        Compare=matrix(nrow=9,ncol=length(MoDs))
        colnames(Compare)=names(MoDs)
        rownames(Compare)=c("Dispersion","Deg.free","Log Like","AIC.c","AIC.delta","AIC.w",
                            "AIC.Ev.ratio","Null.2times.like","Res.2times.like")
        Compare[1,]=unlist(sapply(MoDs,OverDisp))
        Compare[2,]=unlist(sapply(MoDs,Deg.Free))
        Compare[3,]=unlist(sapply(MoDs,LOgLik))
        store.aicc=rep(NA,length(MoDs))
        for(p in 1:length(MoDs))store.aicc[p]=fn.AIC(MoDs[[p]])
        #for(p in 1:length(MoDs))store.aicc[p]=fn.AICc(MoDs[p])
        Compare[4,]=store.aicc
        AIC.w.r=fn.AIC.ratio(store.aicc)
        Compare[5,]=AIC.w.r$Delta
        Compare[6,]=AIC.w.r$Weight
        Compare[7,]=AIC.w.r$Evidence.ratio_how.much.better
        store.deviance=vector("list",length(MoDs))
        #for(p in 1:length(MoDs))store.deviance[[p]]=Dsquared(MoDs[[p]],adjust=T)
        for(p in 1:length(MoDs)) store.deviance[[p]]=Get.dev(MoD=MoDs[[p]],DAT=Store[[i]]$DAT)
        store.deviance=do.call(cbind,store.deviance)
        Compare[8,]=unlist(store.deviance[1,])
        Compare[9,]=unlist(store.deviance[2,])
        Best.error.AIC.w=names(which(Compare[6,]==max(Compare[6,],na.rm=T)))
        Best.formula=Store[[i]]$Forms[[match(Best.error.AIC.w,names(Store[[i]]$Forms))]]
        plot.new()
        mytheme <- gridExtra::ttheme_default(
          core = list(padding=unit(c(1, 1), "mm"),fg_params=list(cex = .5)),
          colhead = list(fg_params=list(cex = .7)),
          rowhead = list(fg_params=list(cex = .7)))
        myt <- gridExtra::tableGrob(Compare, theme = mytheme)
        grid.draw(myt)
        legend('bottom',c(paste("Best.AIC.w=",Best.error.AIC.w),paste("Formula:",Reduce(paste, deparse(Best.formula)))),bty='n',cex=.7)
        dev.off()
      }
      
      
      # Vuong statistic
      #note: usefull for comparing non-nested models (e.g. Poisson Vs zero inflated Poisson)
      #      A large, positive test statistic provides evidence of the superiority of model 1 over model 2,
      #     while a large, negative test statistic is evidence of the superiority of model 2 over model 1.
      #     Had to remove Tweedie, not supported by vuong()
      VUong=function(m1, m2, digits = getOption("digits")) 
      {
        if(!is.null(m1) & !is.null(m2))
        {
          m1y <- m1$y
          m2y <- m2$y
          m1n <- length(m1y)
          m2n <- length(m2y)
          if (m1n == 0 | m2n == 0) 
            stop("Could not extract dependent variables from models.")
          if (m1n != m2n) 
            stop(paste("Models appear to have different numbers of observations.\n", 
                       "Model 1 has ", m1n, " observations.\n", "Model 2 has ", 
                       m2n, " observations.\n", sep = ""))
          if (any(m1y != m2y)) {
            stop(paste("Models appear to have different values on dependent variables.\n"))
          }
          p1 <- predprob(m1)
          p2 <- predprob(m2)
          if (!all(colnames(p1) == colnames(p2))) {
            stop("Models appear to have different values on dependent variables.\n")
          }
          whichCol <- match(m1y, colnames(p1))
          whichCol2 <- match(m2y, colnames(p2))
          if (!all(whichCol == whichCol2)) {
            stop("Models appear to have different values on dependent variables.\n")
          }
          m1p <- rep(NA, m1n)
          m2p <- rep(NA, m2n)
          for (i in 1:m1n) {
            m1p[i] <- p1[i, whichCol[i]]
            m2p[i] <- p2[i, whichCol[i]]
          }
          k1 <- length(coef(m1))
          k2 <- length(coef(m2))
          lm1p <- log(m1p)
          lm2p <- log(m2p)
          m <- lm1p - lm2p
          bad1 <- is.na(lm1p) | is.nan(lm1p) | is.infinite(lm1p)
          bad2 <- is.na(lm2p) | is.nan(lm2p) | is.infinite(lm2p)
          bad3 <- is.na(m) | is.nan(m) | is.infinite(m)
          bad <- bad1 | bad2 | bad3
          neff <- sum(!bad)
          if (any(bad)) {
            cat("NA or numerical zeros or ones encountered in fitted probabilities\n")
            cat(paste("dropping these", sum(bad), "cases, but proceed with caution\n"))
          }
          aic.factor <- (k1 - k2)/neff
          bic.factor <- (k1 - k2)/(2 * neff) * log(neff)
          v <- rep(NA, 3)
          arg1 <- matrix(m[!bad], nrow = neff, ncol = 3, byrow = FALSE)
          arg2 <- matrix(c(0, aic.factor, bic.factor), nrow = neff, 
                         ncol = 3, byrow = TRUE)
          num <- arg1 - arg2
          s <- apply(num, 2, sd)
          numsum <- apply(num, 2, sum)
          v <- numsum/(s * sqrt(neff))
          names(v) <- c("Raw", "AIC-corrected", "BIC-corrected")
          pval <- rep(NA, 3)
          msg <- rep("", 3)
          for (j in 1:3) {
            if (v[j] > 0) {
              pval[j] <- 1 - pnorm(v[j])
              msg[j] <- "model1 > model2"
            }
            else {
              pval[j] <- pnorm(v[j])
              msg[j] <- "model2 > model1"
            }
          }
          out <- data.frame(v, msg, format.pval(pval))
          names(out) <- c("Vuong.z.statistic", "H_A", "p.value")
          
          return(out)
        }
        if(is.null(m1) | is.null(m2))
        {
          out=data.frame('Vuong z-statistic'=NA,'H_A'=NA,'p-value'=NA)
          return(out)
        }
        
      }
      Store.vuong=vector('list',length(Species.cpue))
      for(n in Species.cpue)
      {
        MOD=Store[[n]]$MODS
        
        if(!is.null(MOD))
        {
          if("Tweed"%in%names(names(MOD)))MOD=MOD[-match("Tweed",names(MOD))]
          PAirs=combn(1:length(MOD), 2)
          Names.pairs=combn(names(MOD), 2)
          NaMes=rep(NA,ncol(PAirs));for(x in  1:ncol(PAirs))NaMes[x]=paste(Names.pairs[1,x],"vs",Names.pairs[2,x])
          dummy=vector('list',ncol(PAirs))
          for(p in 1:ncol(PAirs))
          {
            dummy[[p]]=cbind(Comparison=NaMes[p],VUong(MOD[[PAirs[,p][1]]],MOD[[PAirs[,p][2]]]))
          }
          Store.vuong[[n]]=cbind(Species=names(Store)[n],do.call(rbind,dummy))
        }
      }
      Vuong.table=do.call(rbind,Store.vuong)
      write.csv(Vuong.table,paste(hndLs,"Model selection/Voung.compare.models.csv",sep=""))
      Compas=unique(Vuong.table$Comparison)
      fn.compr.vuong=function(SP,XX)
      {
        dd=subset(Vuong.table, Species==SP )
        dd$dummy=substr(row.names(dd),1,3)
        dd=subset(dd,dummy==XX)
        dd$p.value=as.numeric(as.character(dd$p.value))
        
        Misn=Compas[which(!Compas%in%dd$Comparison)]
        if(length(Misn)>0)
        {
          Add=dd[1:length(Misn),]
          Add[,]=NA
          Add$Comparison=Misn
          dd=rbind(dd,Add)
          dd=dd[match(Compas,dd$Comparison),]
        }
        
        dd$shape=with(dd,ifelse(H_A=="model2 > model1",21,ifelse(H_A=="model1 > model2",24,NA)))
        dd$col=with(dd,ifelse(p.value>=0.05,"grey95",
                              ifelse(p.value<0.05 & p.value>=0.01,"grey60",
                                     ifelse(p.value<0.01 & p.value>=0.001,"grey30",
                                            ifelse(p.value<0.001,"black",NA)))))  
        points(1:15,rep(i,15),pch=dd$shape,bg=dd$col,cex=3)
      }
      VuOn.sp=as.character(unique(Vuong.table$Species))
      Voung.metrics=c("Raw","AIC","BIC")
      for(v in 1:length(Voung.metrics))
      {
        fn.fig(paste(hndLs,"Model selection/Voung.compare_",Voung.metrics[v],sep=""),2400,2400)
        par(mar=c(2,2,3,1),oma=c(6,5.5,.1,.1),las=1,mgp=c(2.5,.7,0),xpd=TRUE)
        plot(1:15,1:15,col="transparent",ylab="",xlab="",xaxt='n',yaxt='n',ylim=c(1,length(VuOn.sp)))
        axis(1,1:15,Compas,las=3,cex.axis=.75)
        axis(2,1:length(VuOn.sp),VuOn.sp,las=1,cex.axis=.75)
        for(i in 1:length(VuOn.sp)) fn.compr.vuong(SP=VuOn.sp[i],XX=Voung.metrics[v])
        legend("topright", inset=c(0,-0.1), legend=c("model2 > model1","model1 > model2"), 
               pch=c(21,24), bg="white", bty='n',pt.cex = 1.5)
        legend("topleft", inset=c(0,-0.1), legend=c("p>=0.05 ","p<0.05 ","p<0.01","p<0.001"), 
               pch=22, pt.bg=c("grey95","grey60","grey30","black"), bty='n',horiz=T,pt.cex = 1.5)
        dev.off()
      }
    }
    
    # Compare model terms for each error using AIC, etc
    if(Select.error.AIC=="YES")
    {
      # calculate AIC weights and ratios
      fn.AIC.ratio=function(DAT)
      {
        MIN=min(DAT,na.rm=T)
        Delta=DAT-MIN
        Like.model.give.dat=exp(-Delta/2)
        Weight=Like.model.give.dat/sum(Like.model.give.dat,na.rm=T)
        id=which(Weight==max(Weight,na.rm=T))
        names(Weight)=NULL
        Evidence.ratio=outer(Weight[id],Weight, "/")
        return(list(Best.Mod=id,Delta=Delta,Like=Like.model.give.dat,Weight=Weight,Evidence.ratio_how.much.better=Evidence.ratio))
      }
      
      #Deviance explained
      #note: this doesn't work for zero inflated models. Zero-inflated models are not associated 
      #       with a deviance in the GLM sense. They are a mixture of two models rather 
      #       than a single model from an exponential family. Hence, zeroinfl/hurdle in "pscl" 
      #       reports the log-likelihood but not the deviance. Possibly -2 times the likelihood is
      
      # adjust: logical, whether or not to use the adjusted deviance taking into 
      #acount the nr of observations and parameters (Weisberg 1980; Guisan & Zimmermann 2000)
      
      
      Get.dev=function(MOD,MOD.null) as.numeric(-2*(logLik(MOD.null)-logLik(MOD)))
      
      Dsquared <- function(model, adjust = F)
      {
        if(!is.null(model$deviance))
        {
          d2 <- (model$null.deviance - model$deviance) / model$null.deviance
          if (adjust)
          {
            n <- length(model$fitted.values)
            p <- length(model$coefficients)
            d2 <- 1 - ((n - 1) / (n - p)) * (1 - d2)
          }
          d3=d2*100
          d1=model$deviance
          d2=model$null.deviance
        }
        
        if(is.null(model$deviance))
        {
          d1=NA
          d2=NA
          d3=NA
        }
        return(list(d1=d1,d2=d2,d3=d3))
      }
      
      #increase in deviance explained from one model to next
      Inc.Dev.Exp=function(A,B)round((100*A/B)-100)
      
      
      #LikeRatio test
      fun.LRT=function(MOD)drop1(MOD,test="Chi")
      #fn.LRT1=function(MOD) anova(POIS[[1]],POIS[[2]],POIS[[3]],POIS[[4]], test="Chisq")
      #fn.LRT2=function(MOD) lrtest(POIS[[1]],POIS[[2]],POIS[[3]],POIS[[4]])
      # Like.Ratio1=Like.Ratio2=Dev.explained
      # for(i in 1:N.species)Like.Ratio1[[i]]=fn.LRT1(MODELS[[i]])
      # for(i in 1:N.species)Like.Ratio2[[i]]=fn.LRT2(MODELS[[i]])
      
      
      if(Select.term=="YES")
      {
        #fit nested models for mixture models         Missing: see Zuu et al 2012 page 194
        
        for (i in Species.cpue)
        {
          Full.ZI=Formulas.ZI[[i]][[length(Formulas.ZI[[i]])]]
          Full.Hu=Formulas.Hu[[i]][[length(Formulas.Hu[[i]])]]
        }
        
        
        
        #2.6.4 Store comparisons
        Store.compare.terms=vector('list',N.species)
        names(Store.compare.terms)=names(DATA.list)
        for(i in Species.cpue)
        {
          #1.Count
          COUNT=Store.count[[i]]
          NC=length(COUNT)
          POIS=NB=vector("list",NC)
          for(z in 1:NC)
          {
            POIS[[z]]=COUNT[[z]]$Pois
            NB[[z]]=COUNT[[z]]$NB
          }
          
          #get AICs
          Pois.AIC.w.r=fn.AIC.ratio(unlist(sapply(POIS,AIC)))
          NB.AIC.w.r=fn.AIC.ratio(unlist(sapply(NB,AIC)))
          
          #LikeRatio test
          LR.POIS=fun.LRT(POIS[[NC]])
          LR.NB=fun.LRT(NB[[NC]])
          
          #Deviance explained
          store.deviance=vector("list",NC)
          for(z in 1:NC)store.deviance[[z]]=Dsquared(POIS[[z]],adjust=T)$d3
          Pois.dev.expl=do.call(cbind,store.deviance)
          for(z in 1:NC)store.deviance[[z]]=Dsquared(NB[[z]],adjust=T)$d3
          NB.dev.expl=do.call(cbind,store.deviance)
          
          
          #2.Zero inflated
          ZERO=Store.count.ZI[[i]]
          NC=length(ZERO)
          POIS=NB=vector("list",NC)
          for(z in 1:NC)
          {
            POIS[[z]]=ZERO[[z]]$ZIP
            NB[[z]]=ZERO[[z]]$ZINB
          }
          
          #get AICs
          ZIP.AIC.w.r=fn.AIC.ratio(unlist(sapply(POIS,AIC)))
          ZINB.AIC.w.r=fn.AIC.ratio(unlist(sapply(NB,AIC)))
          
          #LikeRatio test     #MISSING: here goes ZUU page 194
          LR.ZIP=NULL
          LR.ZINB=NULL
          
          #Deviance explained   NA for Zero inflated models
          ZIP.dev.expl=NULL
          ZINB.dev.expl=NULL
          
          
          #3.Tweedie
          TWEEDIE=Store.count.Twee[[i]]
          NC=length(TWEEDIE)
          
          Tweedie.AIC.w.r=NA
          LR.TWEEDIE=NA
          TWEEDIE.dev.expl=NA
          
          if(!is.null(TWEEDIE[[1]]$Tweed))
          {
            #get AICs
            aic.twee=rep(NA,NC)
            for(z in 1:NC) aic.twee[z]=AICtweedie(TWEEDIE[[z]][[1]])
            Tweedie.AIC.w.r=fn.AIC.ratio(aic.twee)
            
            #LikeRatio test
            LR.TWEEDIE=fun.LRT(TWEEDIE[[NC]][[1]])
            
            #Deviance explained
            store.deviance=vector("list",NC)
            for(z in 1:NC)store.deviance[[z]]=Dsquared(TWEEDIE[[z]][[1]],adjust=T)$d3
            TWEEDIE.dev.expl=do.call(cbind,store.deviance)
            
          }
          
          
          #4.Hurdle
          HURDLE=Store.count.Hur[[i]]
          NC=length(HURDLE)
          POIS=NB=vector("list",NC)
          for(z in 1:NC)
          {
            POIS[[z]]=HURDLE[[z]]$Hurdle.Pois
            NB[[z]]=HURDLE[[z]]$Hurdle.NB
          }
          
          #get AICs
          Hur.P.AIC.w.r=fn.AIC.ratio(unlist(sapply(POIS,AIC)))
          Hur.NB.AIC.w.r=fn.AIC.ratio(unlist(sapply(NB,AIC)))
          
          #LikeRatio test     #MISSING: here goes ZUU page 194
          LR.Hur.P=NULL
          LR.Hur.NB=NULL
          
          #Deviance explained   NA for Zero inflated models
          Hur.P.dev.expl=NULL
          Hur.NB.dev.expl=NULL
          
          
          #Store everything
          Store.compare.terms[[i]]=list(Pois.AIC.w.r=Pois.AIC.w.r,NB.AIC.w.r=NB.AIC.w.r,LR.POIS=LR.POIS,LR.NB=LR.NB,
                                        Pois.dev.expl=Pois.dev.expl,NB.dev.expl=NB.dev.expl,ZIP.AIC.w.r=ZIP.AIC.w.r,ZINB.AIC.w.r=ZINB.AIC.w.r,
                                        LR.ZIP=LR.ZIP,LR.ZINB=LR.ZINB,ZIP.dev.expl=ZIP.dev.expl,ZINB.dev.expl=ZINB.dev.expl,
                                        Tweedie.AIC.w.r=Tweedie.AIC.w.r,LR.TWEEDIE=LR.TWEEDIE,TWEEDIE.dev.expl=TWEEDIE.dev.expl,
                                        Hur.P.AIC.w.r=Hur.P.AIC.w.r,Hur.NB.AIC.w.r=Hur.NB.AIC.w.r,LR.Hur.P=LR.Hur.P,
                                        LR.Hur.NB=LR.Hur.NB,Hur.P.dev.expl=Hur.P.dev.expl,Hur.NB.dev.expl=Hur.NB.dev.expl)
        }
        
      }
      
      
      #Select error structure
      
      #FORMULAS
      
      #Poisson and NB
      Best.count.form=list(
        NULL,
        as.formula(Catch.Target ~ year + BOTDEPTH + Mid.Lat +  SOI+ offset(log.Effort)),
        as.formula(Catch.Target ~ year + BOTDEPTH + Mid.Lat +  offset(log.Effort)),
        NULL,
        as.formula(Catch.Target ~ year +  Mid.Lat +  offset(log.Effort)),
        as.formula(Catch.Target ~ year +  Mid.Lat +  Moon + offset(log.Effort)),
        as.formula(Catch.Target ~ year + Moon + Mid.Lat +  SOI+ offset(log.Effort)))
      names(Best.count.form)=names(DATA.list)
      
      #Zero Inflated
      Best.ZI.form=list(
        NULL,
        as.formula(Catch.Target ~ year + BOTDEPTH + Mid.Lat+ offset(log.Effort)|year + Mid.Lat+ offset(log.Effort)),
        as.formula(Catch.Target ~ year +  Moon+ Mid.Lat+ offset(log.Effort)|year +BOTDEPTH+ Mid.Lat+ offset(log.Effort)),
        NULL,
        as.formula(Catch.Target ~ year + BOTDEPTH + Moon+ Mid.Lat+ offset(log.Effort)|year + Moon+ offset(log.Effort)),
        as.formula(Catch.Target ~ year + BOTDEPTH + Moon+ Mid.Lat+ offset(log.Effort)|year + Mid.Lat+ offset(log.Effort)),
        as.formula(Catch.Target ~ year + BOTDEPTH + Moon+ Mid.Lat+ offset(log.Effort)|year + Mid.Lat+ offset(log.Effort)))
      
      #Tweedie
      Best.Twee.form=list(
        NULL,
        as.formula(Catch.Target ~ year + BOTDEPTH + Mid.Lat + SOI + offset(log.Effort)),
        NULL,
        NULL,
        NULL,
        NULL,
        as.formula(Catch.Target ~ year + Moon + Mid.Lat + SOI + offset(log.Effort)))
      
      
      #Hurdle
      Best.HU.Pois.form=list(
        NULL,
        as.formula(Catch.Target ~ year + BOTDEPTH + Moon + Mid.Lat + SOI + offset(log.Effort)),
        as.formula(Catch.Target ~ BOTDEPTH +Moon+Mid.Lat+SOI+ offset(log.Effort) | year + BOTDEPTH+Moon+Mid.Lat+SOI+ offset(log.Effort)),
        NULL,
        as.formula(Catch.Target ~ year + BOTDEPTH + Moon + Mid.Lat + SOI + offset(log.Effort)),
        as.formula(Catch.Target ~ year + BOTDEPTH + Moon + Mid.Lat + offset(log.Effort)),
        as.formula(Catch.Target ~ year + BOTDEPTH + Moon + Mid.Lat + SOI + offset(log.Effort)))
      
      Best.HU.NB.form=Best.HU.Pois.form
      Best.HU.NB.form[[5]]=Catch.Target ~ year + BOTDEPTH + Moon + Mid.Lat + offset(log.Effort)
      
      
    }
    
    # K-fold cross validation
    if(Do.K.fold.test=="YES")
    {
      #formulas of each species and each error structure
      THESE=Best.count.form
      
      THESE$Sandbar=list(Pois=c("Catch.Target","year","BOTDEPTH","Mid.Lat","SOI","N.hooks.Fixed","SOAK.TIME"),
                         NB=c("Catch.Target","year","BOTDEPTH","Mid.Lat","SOI","N.hooks.Fixed","SOAK.TIME"),
                         ZIP=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","N.hooks.Fixed","SOAK.TIME"),
                         ZINB=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","N.hooks.Fixed","SOAK.TIME"),                   
                         Twee=c("Catch.Target","year","BOTDEPTH","Mid.Lat","SOI","N.hooks.Fixed","SOAK.TIME"),
                         Hu.P=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","SOI","N.hooks.Fixed","SOAK.TIME"),
                         Hu.NB=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","SOI","N.hooks.Fixed","SOAK.TIME"))
      
      
      THESE$Tiger=list(Pois=c("Catch.Target","year","BOTDEPTH","Mid.Lat","N.hooks.Fixed","SOAK.TIME"),
                       NB=c("Catch.Target","year","BOTDEPTH","Mid.Lat","N.hooks.Fixed","SOAK.TIME"),
                       ZIP=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","N.hooks.Fixed","SOAK.TIME"),
                       ZINB=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","N.hooks.Fixed","SOAK.TIME"),                   
                       Twee=c("Catch.Target","year","BOTDEPTH","Mid.Lat","N.hooks.Fixed","SOAK.TIME"),
                       Hu.P=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","SOI","N.hooks.Fixed","SOAK.TIME"),
                       Hu.NB=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","SOI","N.hooks.Fixed","SOAK.TIME"))
      
      THESE$Blacktip=list(Pois=c("Catch.Target","year","Mid.Lat","N.hooks.Fixed","SOAK.TIME"),
                          NB=c("Catch.Target","year","Mid.Lat","N.hooks.Fixed","SOAK.TIME"),
                          ZIP=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","N.hooks.Fixed","SOAK.TIME"),
                          ZINB=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","N.hooks.Fixed","SOAK.TIME"),                   
                          Twee=c("Catch.Target","year","Mid.Lat","N.hooks.Fixed","SOAK.TIME"),
                          Hu.P=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","SOI","N.hooks.Fixed","SOAK.TIME"),
                          Hu.NB=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","N.hooks.Fixed","SOAK.TIME"))
      
      THESE$SpotTail=list(Pois=c("Catch.Target","year","Mid.Lat","Moon","N.hooks.Fixed","SOAK.TIME"),
                          NB=c("Catch.Target","year","Mid.Lat","Moon","N.hooks.Fixed","SOAK.TIME"),
                          ZIP=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","N.hooks.Fixed","SOAK.TIME"),
                          ZINB=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","N.hooks.Fixed","SOAK.TIME"),                   
                          Twee=c("Catch.Target","year","Mid.Lat","Moon","N.hooks.Fixed","SOAK.TIME"),
                          Hu.P=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","N.hooks.Fixed","SOAK.TIME"),
                          Hu.NB=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","N.hooks.Fixed","SOAK.TIME"))
      
      THESE$Milk=list(Pois=c("Catch.Target","year","Moon","Mid.Lat","SOI","N.hooks.Fixed","SOAK.TIME"),
                      NB=c("Catch.Target","year","Moon","Mid.Lat","SOI","N.hooks.Fixed","SOAK.TIME"),
                      ZIP=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","N.hooks.Fixed","SOAK.TIME"),
                      ZINB=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","N.hooks.Fixed","SOAK.TIME"),                   
                      Twee=c("Catch.Target","year","Moon","Mid.Lat","SOI","N.hooks.Fixed","SOAK.TIME"),
                      Hu.P=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","SOI","N.hooks.Fixed","SOAK.TIME"),
                      Hu.NB=c("Catch.Target","year","BOTDEPTH","Moon","Mid.Lat","SOI","N.hooks.Fixed","SOAK.TIME"))
      
      
      
      K=20
      training.data=0.95 #training subset (95% of data)
      
      
      
      fn.cor=function(OBS,PREDS)
      {
        CORE=NULL
        if(!is.null(FITS[[a]]))CORE=cor(OBS,PREDS,method='pearson')
        return(list(CORE=CORE))
      }
      fn.RMSE=function(OBS,PREDS)
      {
        CORE=NULL
        if(!is.null(FITS[[a]]))CORE=sqrt(mean((PREDS-OBS)^2))/Rng.obs
        return(list(CORE=CORE))
      }
      
      fn.best.error.k_fold=function(DAT,FORMULA,FORMULA.ZI,FORMULA.Twee,FORMULA.Hu.P,FORMULA.Hu.NB,These)
      {
        
        # 1. Put data in proper shape
        DataFile=DAT[,match(VARIABLES,names(DAT))]
        
        # 2. select training subsets
        id=round(nrow(DataFile)*training.data)
        list.ID=list()
        for(t in 1:K)list.ID[[t]]=sort(sample(1:nrow(DataFile),id,replace =F))
        
        
        #3. loop for k-fold validation
        Cor=RMSE=TEST=PREDS=vector('list',length=K)
        
        for(k in 1:K)
        {
          tryCatch({
            # 3.1. select data subsets
            All.data.train=DataFile[list.ID[[k]],]
            All.data.test=DataFile[-list.ID[[k]],]
            
            
            # 3.2. Fit different error structures to training subset
            DAT=All.data.train
            DAT$year=as.factor(DAT$year)
            DAT$log.Effort=log(DAT$N.hooks.Fixed*DAT$SOAK.TIME)
            
            #GLM
            #Count models
            Pois=glm(FORMULA,data=DAT,family="poisson")
            NB=glm.nb(FORMULA, data = DAT)
            
            #Zero inflated models
            #ZIP=NULL
            ZIP <- zeroinfl(FORMULA.ZI, data = DAT, dist = "poisson")
            #ZINB=NULL
            ZINB <- zeroinfl(FORMULA.ZI, data = DAT, dist = "negbin")
            
            Tweed=NULL
            if(i %in%c(2,7))
            {
              #Tweedie
              # 1.Find p value thru likelihood profiling
              #note: expect  xi  about 1.5
              P.eval=seq(1.1, 1.9, length=9)
              out <- tweedie.profile(FORMULA.Twee,xi.vec=P.eval,data=DAT, do.plot=F)
              P=out$xi.max
              
              # 2.Fit the glm
              Tweed=glm( FORMULA.Twee,data=DAT, family=tweedie(var.power=P, link.power=0) )
            }
            
            #Hurdle modles
            Hurdle.Pois = hurdle(FORMULA.Hu.P, data = DAT , dist="poisson")
            Hurdle.NB = hurdle(FORMULA.Hu.NB, data = DAT , dist="negbin")
            
            FITS=list(Pois=Pois,NB=NB,ZIP=ZIP, ZINB=ZINB,Tweed=Tweed,Hurdle.Pois=Hurdle.Pois,
                      Hurdle.NB=Hurdle.NB)  
            
            # 3.3. Predict test.data  
            Store.Preds=vector('list',length(FITS))
            names(Store.Preds)=names(FITS)
            for (a in 1:length(FITS))Store.Preds[[a]]=Pred.Model(FITS[[a]],All.data.test[,match(These[[a]],names(All.data.test))])
            
            #Store quantities
            TEST[[k]]=All.data.test
            PREDS[[k]]=Store.Preds
            
            
            # 3.4. Measure fit
            
            #3.4.1 Correlation
            Store.Corr=vector('list',length(FITS))
            names(Store.Corr)=names(FITS)      
            for (a in 1:length(FITS))Store.Corr[[a]]=fn.cor(All.data.test$Catch.Target,Store.Preds[[a]]$PREDS)
            Cor[[k]]=Store.Corr
            
            
            #3.4.2 Normalised RMSE
            Store.RMSE=vector('list',length(FITS))
            names(Store.RMSE)=names(FITS)
            Rng.obs=range(All.data.test$Catch.Target)[2]-range(All.data.test$Catch.Target)[1]       
            for (a in 1:length(FITS))Store.RMSE[[a]]=fn.RMSE(All.data.test$Catch.Target,Store.Preds[[a]]$PREDS)      
            RMSE[[k]]=Store.RMSE
            
            
          }, error = function(e) {
            
          })
          
        }
        
        return(list(Cor=Cor,RMSE=RMSE,TEST=TEST,PREDS=PREDS))
      }
      
      Store.CrossVal=vector('list',N.species)
      names(Store.CrossVal)=TARGETS.name
      system.time(for(i in Species.cpue)
      {
        Store.CrossVal[[i]]=fn.best.error.k_fold(DATA.list[[i]],Best.count.form[[i]],Best.ZI.form[[i]],
                                                 Best.Twee.form[[i]],Best.HU.Pois.form[[i]],Best.HU.NB.form[[i]],THESE[[i]])
      })
      
      
      #extract values for average cor and RMSE
      See.K.fold.fun.mean=function(StOrE)
      {
        CoR=unlist(StOrE$Cor)
        RmSe=unlist(StOrE$RMSE)
        Cor.names=unique(names(CoR))
        LISTA=vector("list",length(Cor.names))
        for(t in 1:length(LISTA))LISTA[[t]]=which(Cor.names[t]==names(CoR))
        Data.cor=Data.rmse=matrix(NA,nrow=length(LISTA[[1]]),ncol=length(Cor.names))
        colnames(Data.cor)=colnames(Data.rmse)=Cor.names
        for(x in 1:ncol(Data.cor)) Data.cor[,x]=CoR[LISTA[[x]]]
        for(x in 1:ncol(Data.cor)) Data.rmse[,x]=RmSe[LISTA[[x]]]
        
        Mean.cor=colMeans(Data.cor)
        Mean.rmse=colMeans(Data.rmse)
        
        return(list(Mean.cor=Mean.cor,Mean.rmse=Mean.rmse))
      }
      Store.k.fold.Cross.v=vector('list',N.species)
      names(Store.k.fold.Cross.v)=names(DATA.list)
      for(i in Species.cpue) Store.k.fold.Cross.v[[i]]=See.K.fold.fun.mean(Store.CrossVal[[i]])
      
      
      #extract values for overall cor 
      See.K.fold.fun.overall=function(StOrE)
      {
        Test.Dat=do.call(rbind,StOrE$TEST)
        Preds.Pois=Preds.NB=Preds.ZIP=Preds.ZINB=Preds.Tweed=Preds.Hurdle.Pois=Preds.Hurdle.NB=NULL
        for(q in 1:length(StOrE$PREDS))
        {
          Preds.Pois=c(Preds.Pois,StOrE$PREDS[[q]]$Pois$PREDS)
          Preds.NB=c(Preds.NB,StOrE$PREDS[[q]]$NB$PREDS)
          if(!is.null(StOrE$PREDS[[q]]$ZIP$PREDS))Preds.ZIP=c(Preds.ZIP,StOrE$PREDS[[q]]$ZIP$PREDS)
          if(!is.null(StOrE$PREDS[[q]]$ZINB$PREDS))Preds.ZINB=c(Preds.ZINB,StOrE$PREDS[[q]]$ZINB$PREDS)
          if(!is.null(StOrE$PREDS[[q]]$Tweed$PREDS))Preds.Tweed=c(Preds.Tweed,StOrE$PREDS[[q]]$Tweed$PREDS)
          Preds.Hurdle.Pois=c(Preds.Hurdle.Pois,StOrE$PREDS[[q]]$Hurdle.Pois$PREDS)
          Preds.Hurdle.NB=c(Preds.Hurdle.NB,StOrE$PREDS[[q]]$Hurdle.NB$PREDS)
        }
        
        FITS=names(StOrE$PREDS[[1]])
        Store.Corr=vector('list',length(FITS))
        names(Store.Corr)=FITS
        Store.Corr$Pois=cor(Test.Dat$Catch.Target[1:length(Preds.Pois)],Preds.Pois,method='pearson')
        Store.Corr$NB=cor(Test.Dat$Catch.Target[1:length(Preds.NB)],Preds.NB,method='pearson')
        Store.Corr$ZIP=cor(Test.Dat$Catch.Target[1:length(Preds.ZIP)],Preds.ZIP,method='pearson')
        Store.Corr$ZINB=cor(Test.Dat$Catch.Target[1:length(Preds.ZINB)],Preds.ZINB,method='pearson')
        if(!is.null(Preds.Tweed))Store.Corr$Tweed=cor(Test.Dat$Catch.Target[1:length(Preds.Tweed)],Preds.Tweed,method='pearson')
        Store.Corr$Hurdle.Pois=cor(Test.Dat$Catch.Target[1:length(Preds.Hurdle.Pois)],Preds.Hurdle.Pois,method='pearson')
        Store.Corr$Hurdle.NB=cor(Test.Dat$Catch.Target[1:length(Preds.Hurdle.NB)],Preds.Hurdle.NB,method='pearson')
        
        
        return(Store.Corr=Store.Corr)
        
      }
      Store.k.fold.Cross.v.overall=vector('list',N.species)
      names(Store.k.fold.Cross.v.overall)=names(DATA.list)
      for(i in Species.cpue) Store.k.fold.Cross.v.overall[[i]]=See.K.fold.fun.overall(Store.CrossVal[[i]])
      
    }
  }

  
  #1.11.3  Fit Best model and error structure 
  if(do.GAM=="YES")
  {
    BEST.model=vector('list',N.species)
    names(BEST.model)=names(DATA.list)
    BEST.model$'Milk shark'=BEST.model$'Tiger shark'=
      formula(paste(Res.var,paste(c("year","s(Mid.Lat,k=3)",
                                    "s(BOTDEPTH,k=3)",Offset),collapse="+"),sep="~"))
    BEST.model$'Sandbar shark'=formula(paste(Res.var,paste(c("year","Month","s(Mid.Lat,k=3)",
                                    "s(BOTDEPTH,k=3)","s(Temp.res,k=3)","Moon",Offset),
                                    collapse="+"),sep="~"))
    BEST.model$`Spot-tail shark`=BEST.model$`Blacktip sharks`=
      formula(paste(Res.var,paste(c("year","s(Mid.Lat,k=3)",
                                    "s(BOTDEPTH,k=3)","s(Temp.res,k=3)",Offset),collapse="+"),sep="~"))
    BEST.model$'Scalloped hammerhead'= formula(paste(Res.var,paste(c("year",
                                         "s(BOTDEPTH,k=3)",Offset),collapse="+"),sep="~"))
    BEST.model$"Dusky shark"=formula(paste(Res.var,paste(c("year","s(Mid.Lat,k=3)","Moon",
                                     Offset),collapse="+"),sep="~"))
    BEST.model$'Sliteye shark'=formula(paste(Res.var,paste(c("year",
                                "s(BOTDEPTH,k=3)",Offset),collapse="+"),sep="~"))   
    ERROR.st=BEST.model
    for(i in 1:N.species)ERROR.st[[i]]="NB"
  }
  
  if(do.GLM=="YES")
  {
    BEST.model=vector('list',N.species)
    names(BEST.model)=names(DATA.list)
    BEST.model$"Sandbar shark"=as.formula(Catch.Target ~ year+log.Mid.Lat+log.BOTDEPTH+offset(log.Effort))
    
    BEST.model$"Milk shark"=BEST.model$"Spot-tail shark"=
      BEST.model$"Scalloped hammerhead"=BEST.model$"Tiger shark"=
      BEST.model$"Dusky shark"=BEST.model$"Sandbar shark"
    
    BEST.model$"Blacktip sharks"=
      as.formula(Catch.Target ~ year + log.BOTDEPTH + offset(log.Effort) | 
                   year + log.Mid.Lat + log.BOTDEPTH + offset(log.Effort))
    BEST.model$"Sliteye shark"=
      as.formula(Catch.Target ~ year + log.Mid.Lat + offset(log.Effort) | 
                   year + log.Mid.Lat + log.BOTDEPTH + offset(log.Effort))
    
    ERROR.st=BEST.model
    ERROR.st$"Sandbar shark"="HU.NB"
    ERROR.st$"Milk shark"=ERROR.st$"Blacktip sharks"="ZINB"
    ERROR.st$"Tiger shark"=ERROR.st$"Dusky shark"="NB"
    ERROR.st$"Spot-tail shark"=ERROR.st$"Scalloped hammerhead"=ERROR.st$"Sliteye shark"="ZIP"
  }


  fit.best=function(d,FORMULA,ErroR)
  {
    if(do.GAM=="YES")
    {
      if(ErroR=="Pois") Fit=gam(FORMULA,data=d,method = "REML",family=poisson)
      if(ErroR=="NB")  Fit=gam(FORMULA, data=d,method = "REML",family = nb)
      if(ErroR=="ZIP") Fit=zipgam(lambda.formula=FORMULA,pi.formula=FORMULA,data=d)
      if(ErroR=="ZINB") Fit=zinbgam(mu.formula=FORMULA,pi.formula=FORMULA, data=d)
      
      #Predictions
      year.pred=summary(emmeans(Fit,"year", type="response"))
      
      Lat.pred=NULL
      used.term=grepl("Mid.Lat", attr(terms(FORMULA),'term.labels'))
      if(sum(used.term)>0)
      {
        Lata=with(subset(d,Catch.Target>0),range(abs(floor(d$Mid.Lat))))
        Lat.pred=summary(emmeans(Fit,"Mid.Lat", type="response",at=list(Mid.Lat=seq(Lata[1],Lata[2],.1))))
      }
      
      Depth.pred=NULL
      used.term=grepl("BOTDEPTH", attr(terms(FORMULA),'term.labels'))
      if(sum(used.term)>0)
      {
        ZZ=with(subset(d,Catch.Target>0),range(10*floor(BOTDEPTH/10)))
        Depth.pred=summary(emmeans(Fit,"BOTDEPTH", type="response",at=list(BOTDEPTH=seq(ZZ[1],ZZ[2],10))))
      }
    }
    
    if(do.GLM=="YES")
    {
      if(ErroR=="Pois") Fit=Count.Pois.fn(DAT=d,FORMULA=FORMULA)
      if(ErroR=="NB") Fit=Count.NB.fn(DAT=d,FORMULA=FORMULA)
      if(ErroR=="ZIP") Fit=zeroinfl(FORMULA, data = d, dist = "poisson")
      if(ErroR=="ZINB") Fit=zeroinfl(FORMULA, data = d, dist = "negbin")
      if(ErroR=="HU.P") Fit=hurdle(FORMULA, data = d , dist="poisson")
      if(ErroR=="HU.NB") Fit=hurdle(FORMULA, data = d , dist="negbin")
      
      #Predictions
      year.pred=summary(emmeans(Fit,"year", type="response"))
      
      Lat.pred=NULL
      used.term=grepl("log.Mid.Lat", attr(terms(FORMULA),'term.labels'))
      if(sum(used.term)>0)
      {
        Lata=with(subset(d,Catch.Target>0),range(abs(floor(d$Mid.Lat))))
        Lat.pred=summary(emmeans(Fit,"log.Mid.Lat", type="response",at=list(log.Mid.Lat=log(seq(Lata[1],Lata[2],.1)))))
      }
      
      Depth.pred=NULL
      used.term=grepl("log.BOTDEPTH", attr(terms(FORMULA),'term.labels'))
      if(sum(used.term)>0)
      {
        ZZ=with(subset(d,Catch.Target>0),range(10*floor(BOTDEPTH/10)))
        Depth.pred=summary(emmeans(Fit,"log.BOTDEPTH", type="response",at=list(log.BOTDEPTH=log(seq(ZZ[1],ZZ[2],10)))))
      }
    }
    
    #Null model
    Null=update(Fit,.~1)
    
    return(list(Fit=Fit,DAT=d,Null=Null,year.pred=year.pred,
                Lat.pred=Lat.pred,Depth.pred=Depth.pred))
  }

  Store=vector('list',N.species)   
  names(Store)=names(DATA.list)
  system.time(for(i in Species.cpue)     
  {
    DAT=subset(DATA.list[[i]],FixedStation=="YES" & BOTDEPTH<210)
    TT=with(subset(DAT,Catch.Target>0),table(year))
    d=DAT %>% filter(year%in%names(TT[TT>0])) %>%
      mutate(Mid.Lat=abs(Mid.Lat),
             log.Ef=log(SOAK.TIME*N.hooks.Fixed),
             year=factor(year,levels=unique(sort(year))),
             Month=factor(Month,levels=unique(sort(Month))),
             Moon=factor(Moon,levels=c("Full","Waning","New","Waxing")))
    
    Store[[i]]=fit.best(d=d,FORMULA=BEST.model[[i]],ErroR=ERROR.st[[i]])
    rm(DAT,TT,d)
  })
  
  #Is there a non-linear relation in residuals that requires GAM?    NOT applicable anymore as I am already using GAMs....
  #note: if Pearson residuals vs the covariate shows a clear pattern 
  #       (i.e. Significant Anova), then GAM is needed (i.e. covariates
  #       must be modelled thru a polynomial as it is not linearly related
  #       to the response variable)
  if(check.GAM=="YES")
  {
    fn.need.Gam=function(MOD,DAT,VAR,Var.name)
    {
      E1=resid(MOD,type='pearson')
      plot(VAR,E1,ylab='',xlab=Var.name)
      lo <- loess(E1~VAR)
      lines(predict(lo), col='red', lwd=2)
      
      tmp=gam(E1~s(VAR),data=DAT)
      return(anova(tmp))
    }
    
    #note: if anova is Significant, then there's a pattern in residuals (ie non-linear relation
    #       between covariate and response var)
    fn.fig(paste(hndl.expl,"/Need_gam_plot",sep=""),2400,1200)
    par(mfcol=c(2,length(Species.cpue)),las=1,mar=c(2.5,.5,.5,.75),oma=c(1,2.25,1,.1),mgp=c(1.5,0.6,0))
    STORE.GAM=vector('list',N.species)
    names(STORE.GAM)=names(Store)
    for(i in Species.cpue)
    {
      MoD=Store[[i]]$Fit
      DAt=Store[[i]]$DAT
      ths=c("BOTDEPTH","Mid.Lat")
      id=match(ths,names(DAt))
      Store.Gam=vector('list',length(ths))
      names(Store.Gam)=ths
      for(x in 1:length(id))
      {
        Store.Gam[[x]]=fn.need.Gam(MoD,DAt,DAt[,id[x]],ths[x])
        if(x==1)mtext(names(Store)[i],3,line=0,cex=.65)
        legend('topright',paste("F=",round(Store.Gam[[x]]$s.table[,3],2),
                                "; p-value=",round(Store.Gam[[x]]$s.table[,4],2)),bty='n',cex=.65)
      }
      STORE.GAM[[i]]=Store.Gam
    }
    mtext("Pearson residuals",2,outer=T,las=3,line=0.5,cex=1.25)
    dev.off()
    
    #1.11.14. Fitting GAM and ZIGAM
    # 
    #   #Formulas
    #   
    #   #Count model
    #   pois.GAM=gam(Catch.Target ~ year s(BOTDEPTH)+ offset(log.Effort), family = poisson, data=DAT)
    #   NB.GAM=gam(Catch.Target ~ year s(BOTDEPTH)+ offset(log.Effort), family = negbin(NB$theta), data=DAT)
    #   
    #   #Zero inflated (see page 84-86 Zuur, how to get deviance)
    #   library(VGAM) 
    #   detach("package:mgcv")
    #   ZI.pois.GAM=vgam(Catch.Target ~ year +s(BOTDEPTH)+ offset(log.Effort), family = zipoisson, data=DAT)
    #   ZI.NB.GAM=vgam(Catch.Target ~ year +s(BOTDEPTH)+ offset(log.Effort), family = zinegbinomial, data=DAT,
    #                  control=vgam.control(maxit=100,epsilon=1e-4))
    
  }
  
  ##compare fitted model to null model to test if it's an improvement
   Ch.sq=ERROR.st
   for(i in Species.cpue) Ch.sq[[i]]=pchisq(2 * (logLik(Store[[i]]$Fit) - logLik(Store[[i]]$Null)), df = 3, lower.tail = FALSE)

    #rootgram
  fn.fig(paste(getwd(),"/Model selection/rootgram_abundance",sep=""),2000,2400)
  smart.par(n.plots=length(Species.cpue),MAR=c(2,2,1,1),OMA=c(1,1.5,.1,.1),MGP=c(2.5,.7,0))
  for(i in Species.cpue) rootogram(Store[[i]]$Fit, main = TARGETS.name[i])
  dev.off()
  
     #Anova tables    
  Anova.tab=function(mod,Fcol)  
  {
    ANVA=anova.gam(mod)
    Param=ANVA$pTerms.table
    Non.param=ANVA$s.table
    Tab=rbind(Param[,match(c(Fcol,"p-value"),colnames(Param))],
              Non.param[,match(c(Fcol,"p-value"),colnames(Non.param))])
    row.names(Tab)=NULL
    Tab=cbind(data.frame(Term=tolower(c(row.names(Param), row.names(Non.param)))),
              as.data.frame(Tab))
    Tab=Tab%>%rename(p=`p-value`)%>%mutate(p=ifelse(p<0.001,"<0.001",round(p,3)))
    id=match(Fcol,names(Tab))
    Tab[,id]=round(Tab[,id],3)
    dummy=Tab[1,]
    dummy[,]=""
    SP=dummy
    SP$Term=TARGETS.name[i]
    Dev.Exp=dummy
    Dev.Exp$Term="Dev.exp"
    Dev.Exp[,id]=paste(round(ANVA$dev.expl,2)*100,"%",sep="")
    Tab=rbind(SP,Tab,Dev.Exp)
    Tab$Term[which(grepl("mid.lat",Tab$Term))]="latitude"
    Tab$Term[which(grepl("botdepth",Tab$Term))]="bottom depth"
    Tab$Term[which(grepl("Dev.exp",Tab$Term))]="Deviance explained"
    Tab$Join=with(Tab,paste(TARGETS.name[i],Tab$Term,sep="_"))

    return(Tab)
  }
 
  
    #Get terms significance and coefficients 
  Sig.terms=vector('list',N.species)
  names(Sig.terms)=TARGETS.name
  Sig.term.coeff=Sig.terms
  for (i in 1:N.species)
  {
    a=as.data.frame.matrix(coeftest(Store[[i]]$Fit))
    write.csv(a,paste("Paper/Anovas/Anova.abundance_",names(Store)[i],".csv",sep=""),row.names=T)
    
    id=match('Pr(>|z|)',names(a))
    if(length(id)>0) names(a)[id]='Pr(>|t|)'
    
    Sigi=rownames(a[which(a$`Pr(>|t|)`<0.05),])
    Sigi=Sigi[!grepl("\\bIntercept\\b",Sigi)]
    Sigi=Sigi[!grepl('year',Sigi)]
    Sigi.count=sapply(Sigi[grepl('count',Sigi)], function(x) substr(x,7,30))
    Sigi.zero=sapply(Sigi[grepl('zero',Sigi)], function(x) substr(x,6,30))
    if(ERROR.st[[i]]%in%c("Pois","NB")) Sigi.count=Sigi
    Sig.terms[[i]]=list(Sigi.count=Sigi.count,Sigi.zero=Sigi.zero)
    
    dd=a$`Pr(>|t|)`[match(Sigi,row.names(a))]
    names(dd)=Sigi
    if(!ERROR.st[[i]]%in%c("Pois","NB")) 
    {
      Sigi.count=dd[grepl('count',names(dd))]
      names(Sigi.count)=sapply( names(Sigi.count)[grepl('count', names(Sigi.count))], function(x) substr(x,7,30))
      Sigi.zero=dd[grepl('zero',names(dd))]
      names(Sigi.zero)=sapply( names(Sigi.zero)[grepl('zero', names(Sigi.zero))], function(x) substr(x,6,30))
      Sig.term.coeff[[i]]=list(Sigi.count=Sigi.count,Sigi.zero=Sigi.zero)
    }
    if(ERROR.st[[i]]%in%c("Pois","NB")) Sig.term.coeff[[i]]=list(Sigi.count=dd,Sigi.zero=NULL)
    
  }
  
    #Wald statistic
  #Fit is full model, Fit1 is full model - term
  #waldtest(Fit,Fit1)
  
    #stepwise likelihood ratios
  if(do.like.ratio.test=="YES")
  {
    lrtest(Full.mod,Nested.mod)  #use this
    
    fun.LRT.mixture=function(M1,M2,df)
    {
      #M1 is full model, M2 is the nested model without the term
      L1=logLik(M1)
      L2=logLik(M2)
      L=2*(L1-L2)
      L=abs(as.numeric(L))
      p=1-pchisq(L,df)
      output=round(c(L,p),digits=4)
      names(output)=c("LikeRatio","p value")
      return(output)
      
    }
    
    #Nested models
    Nested.mod=BEST.model
    Nested.mod$`Sandbar shark`=list(
      as.formula(Catch.Target ~ year+log.Mid.Lat+log.BOTDEPTH+offset(log.Effort)|
                   year+log.Mid.Lat+offset(log.Effort)),
      as.formula(Catch.Target ~ year+log.Mid.Lat+log.BOTDEPTH+offset(log.Effort)|
                   year+log.BOTDEPTH+offset(log.Effort)),
      as.formula(Catch.Target ~ year+log.BOTDEPTH+offset(log.Effort)|
                   year+log.Mid.Lat+log.Mid.Lat+offset(log.Effort)),
      as.formula(Catch.Target ~ year+log.Mid.Lat+offset(log.Effort)|
                   year+log.Mid.Lat+log.BOTDEPTH+offset(log.Effort))
    )
    Nested.mod$`Milk shark`=Nested.mod$`Spot-tail shark`=Nested.mod$`Sandbar shark`
    
    Nested.mod$`Blacktip sharks`=list(
      as.formula(Catch.Target ~ log.Mid.Lat + log.BOTDEPTH + offset(log.Effort) | 
                   year + log.BOTDEPTH + offset(log.Effort)),
      as.formula(Catch.Target ~ log.Mid.Lat + log.BOTDEPTH + offset(log.Effort) | 
                   year + log.Mid.Lat  + offset(log.Effort)),
      as.formula(Catch.Target ~  log.BOTDEPTH + offset(log.Effort) | 
                   year +log.Mid.Lat + log.BOTDEPTH + offset(log.Effort)),
      as.formula(Catch.Target ~ log.Mid.Lat  + offset(log.Effort) | 
                   year + log.Mid.Lat + log.BOTDEPTH + offset(log.Effort))
    )
    
    Nested.mod$`Tiger shark`=Nested.mod$`Scalloped hammerhead`=
    Nested.mod$`Sliteye shark`=Nested.mod$`Dusky shark`=Nested.mod$`Blacktip sharks`

    Store.Likes=vector('list',N.species)
    names(Store.Likes)=names(DATA.list)
    
    system.time(for(i in Species.cpue)
    {
      Forms=Nested.mod[[i]]
      Store.Fits=vector('list',length(Forms))
      for (j in 1:length(Forms))Store.Fits[[j]]=fit.best(subset(DATA.list[[i]],FixedStation=="YES"),
                                                         Forms[[j]],ERROR.st[[i]])
      Store.Likes[[i]]=Store.Fits
    })
    
    Store.LRT=vector('list',N.species)
    names(Store.LRT)=names(DATA.list)
    for(i in Species.cpue)
    {
      M1=Store[[i]]
      Nested=Store.Likes[[i]]
      Stores=vector('list',length(Nested))
      for(j in 1:length(Nested))
      {
        Stores[[j]]=NULL
        if(!is.null(Nested[[j]]))Stores[[j]]=fun.LRT.mixture(M1,Nested[[j]],df=Nested[[j]]$df.residual-M1$df.residual)
      }
      
      Store.LRT[[i]]=Stores
    }
    
  }
  
  
  #1.12.5. Predict year, lat and depth effect and get confidence intervals
  
    #1.12.5.1 lsmeans approach (asymptotic errors)
  PRED.CPUE=vector('list',length(Species.cpue))
  names(PRED.CPUE)=names(DATA.list)
  PRED.lat=PRED.z=PRED.CPUE
  for(i in Species.cpue)
  {
    PRED.CPUE[[i]]=Store[[i]]$year.pred
    if(!is.null(Store[[i]]$Lat.pred))PRED.lat[[i]]=Store[[i]]$Lat.pred
    if(!is.null(Store[[i]]$Depth.pred))PRED.z[[i]]=Store[[i]]$Depth.pred
  }

    #1.12.5.2 Bootstrapping approach   #takes 0.2 secs per n.boot iteration   
  system.time({
    if(do.boot=="YES")
    {
      library(sampling)
      fn.boot=function(dd)   # stratified sampling with replacement
      {
        s=strata(dd,c("year"),size=table(dd$year), method="srswr") 
        boot.d=getdata(dd,s)
        return(boot.d)
      }
      fit.best.boot=function(DAT,FORMULA,ErroR) # fit models to bootstrapped data
      {
        if(ErroR=="ZIP") Fit=zeroinfl(FORMULA, data = DAT, dist = "poisson")
        if(ErroR=="ZINB") Fit=zeroinfl(FORMULA, data = DAT, dist = "negbin")
        if(ErroR=="HU.P") Fit=hurdle(FORMULA, data = DAT , dist="poisson")
        if(ErroR=="HU.NB") Fit=hurdle(FORMULA, data = DAT , dist="negbin")
        return(list(Fit=Fit,DAT=DAT))
      }
      Pred.Model.boot=function(MOD,Dta)
      {
        if(!is.null(MOD))PREDS=predict(MOD,Dta, type='response')
        return(list(PREDS=PREDS))
      }
      Pred.Model.boot_count_zero=function(MOD,Dta,PRd.knt,PRd.zro)
      {
        PREDS.count=PREDS.zero=NULL
        if(!is.null(MOD) & PRd.knt=='YES')PREDS.count=predict(MOD,Dta, type='count')
        if(!is.null(MOD) & PRd.zro=='YES')PREDS.zero=predict(MOD,Dta, type='prob')[,1]
        return(list(PREDS.count=PREDS.count,PREDS.zero=PREDS.zero))
      }
      fn.seQ=function(Min,Max) seq(round(Min),round(Max)) 
      
      #fixed stations
      cores=detectCores()      #setup parallel backend to use many processors
      cl <- makeCluster(cores[1]-1) #leave 1 core not to overload your computer
      registerDoParallel(cl)
      STORE.BOOT=vector('list',N.species)
      names(STORE.BOOT)=names(DATA.list)
      for(i in Species.cpue)   
      {
        #parallel processing
        Store.boot=foreach(k=1:n.boot,.errorhandling='remove',.packages=c('sampling','pscl')) %dopar%
        {
          mod=fit.best.boot(fn.boot(dd=Store[[i]]$DAT),BEST.model[[i]],ERROR.st[[i]])
          return(mod)
          rm(mod)
        }
        STORE.BOOT[[i]]=Store.boot
      }
      stopCluster(cl)  #stop cluster
      
      
      #predict each boot run   
      
      #Vars=all.vars(BEST.model$`Sandbar shark`)[-1]
      
      #Year effect
        #standardised
      PRED.CPUE=vector('list',N.species)
      names(PRED.CPUE)=names(DATA.list)
      for(i in Species.cpue)
      {
        dat=subset(Store[[i]]$DAT,Catch.Target>0)  #within species range
        NEWDATA=data.frame(year=factor(levels(dat$year)),
                           log.BOTDEPTH=mean(dat$log.BOTDEPTH,na.rm=T),
                           log.Mid.Lat=mean(dat$log.Mid.Lat,na.rm=T),
                           log.Effort=mean(dat$log.Effort,na.rm=T))
        MOD=STORE.BOOT[[i]]
        An.Ktch=matrix(nrow=length(MOD),ncol=length(NEWDATA$year))
        colnames(An.Ktch)=NEWDATA$year
        for (x in 1:length(MOD))
        {
          a=Pred.Model.boot(MOD[[x]]$Fit,NEWDATA)$PREDS
          if(!is.null(a))An.Ktch[x,]=a
        }
        PRD=data.frame(Year=NEWDATA$year,
                       MEAN=apply(An.Ktch, 2, function(x) median(x, na.rm=T)),
                       SD=apply(An.Ktch,2,function(x) sd(x,na.rm=T)),
                       LOW=apply(An.Ktch, 2, function(x) quantile(x, 0.025,na.rm=T)),
                       UP=apply(An.Ktch, 2, function(x) quantile(x, 0.975,na.rm=T)))
        PRD$CV=(PRD$SD/PRD$MEAN)*100
        PRED.CPUE[[i]]=PRD
        rm(dat,MOD)
      }
        #nominal
      nominal.boot=function(DAT) 
      {
        DAT$cpue=DAT$Catch.Target/(DAT$N.hooks.Fixed*DAT$SOAK.TIME)
        Nmnl = DAT %>% 
          group_by(year) %>%
          summarise(n = length(cpue),
                    m = length(cpue[cpue>0]),
                    mean = mean(cpue[cpue>0])) %>%
          mutate(p.nz = m/n,
                 mean.delta = p.nz*mean)%>%
          as.data.frame
        Mean.delta=Nmnl$mean.delta
        names(Mean.delta)=Nmnl$year
        return(Mean.delta)
      }
      Nominal.CPUE=vector('list',N.species)
      names(Nominal.CPUE)=names(DATA.list)
      for(i in Species.cpue)
      {
        dummy=vector('list',length(STORE.BOOT[[i]]))
        for(b in 1:length(dummy))dummy[[b]]=nominal.boot(DAT=STORE.BOOT[[i]][[b]]$DAT)
        dummy=do.call(rbind,dummy)
        Nominal.CPUE[[i]]=data.frame(year=as.numeric(colnames(dummy)),
                                     mean=apply(dummy, 2, function(x) median(x, na.rm=T)),
                                     low95=apply(dummy, 2, function(x) quantile(x, 0.025,na.rm=T)),
                                     up95=apply(dummy, 2, function(x) quantile(x, 0.975,na.rm=T)))
      }
      
    
      #Latitude effect  
      PRED.lat=vector('list',length(Species.cpue))
      names(PRED.lat)=names(DATA.list)
      PRED.lat.count=PRED.lat.zero=PRED.lat
      for(i in Species.cpue)
      {
        dat=subset(Store[[i]]$DAT,Catch.Target>0)  #predict within species range
        Yr.ref=names(sort(table(dat$year)))
        Yr.ref=Yr.ref[length(Yr.ref)]
        NEWDATA=data.frame(year=factor(Yr.ref,levels(dat$year)),
                           log.BOTDEPTH=mean(dat$log.BOTDEPTH,na.rm=T),
                           log.Mid.Lat=log(fn.seQ(min(dat$Mid.Lat,na.rm=T),max(dat$Mid.Lat,na.rm=T))),
                           log.Effort=mean(dat$log.Effort,na.rm=T))
        MOD=STORE.BOOT[[i]]
        An.Ktch=matrix(nrow=length(MOD),ncol=length(NEWDATA$log.Mid.Lat))
        colnames(An.Ktch)=exp(NEWDATA$log.Mid.Lat)
        An.Ktch.count=An.Ktch.zero=An.Ktch
        PRd.knt=PRd.zro="NO"  #select species for which lat is significant
        if("log.Mid.Lat"%in%Sig.terms[[i]]$Sigi.count) PRd.knt="YES"
        if("log.Mid.Lat"%in%Sig.terms[[i]]$Sigi.zero) PRd.zro="YES"
        
        if(PRd.knt=="YES")
        {
          for (x in 1:length(MOD))An.Ktch.count[x,]=Pred.Model.boot_count_zero(MOD[[x]]$Fit,NEWDATA,PRd.knt,PRd.zro)$PREDS.count
          PRD=data.frame(Mid.Lat=exp(NEWDATA$log.Mid.Lat),
                         MEAN=apply(An.Ktch.count, 2, function(x) median(x, na.rm=T)),
                         SD=apply(An.Ktch.count,2,function(x) sd(x,na.rm=T)),
                         LOW=apply(An.Ktch.count, 2, function(x) quantile(x, 0.025,na.rm=T)),
                         UP=apply(An.Ktch.count, 2, function(x) quantile(x, 0.975,na.rm=T)))
          PRD$CV=PRD$SD/PRD$MEAN*100
          PRED.lat.count[[i]]=PRD
        }
        if(PRd.zro=="YES")
        {
          for (x in 1:length(MOD))An.Ktch.zero[x,]=Pred.Model.boot_count_zero(MOD[[x]]$Fit,NEWDATA,PRd.knt,PRd.zro)$PREDS.zero
          PRD=data.frame(Mid.Lat=exp(NEWDATA$log.Mid.Lat),
                         MEAN=apply(An.Ktch.zero, 2, function(x) median(x, na.rm=T)),
                         SD=apply(An.Ktch.zero,2,function(x) sd(x,na.rm=T)),
                         LOW=apply(An.Ktch.zero, 2, function(x) quantile(x, 0.025,na.rm=T)),
                         UP=apply(An.Ktch.zero, 2, function(x) quantile(x, 0.975,na.rm=T)))
          PRD$CV=PRD$SD/PRD$MEAN*100
          PRED.lat.zero[[i]]=PRD
        }
        rm(dat,MOD)
      }
      
      #Depth effect
      PRED.z=vector('list',length(Species.cpue))
      names(PRED.z)=names(DATA.list)
      PRED.z.count=PRED.z.zero=PRED.z
      for(i in Species.cpue)
      {
        dat=subset(Store[[i]]$DAT,Catch.Target>0)  #within species range
        Yr.ref=names(sort(table(dat$year)))
        Yr.ref=Yr.ref[length(Yr.ref)]
        Z.seq=10*round(seq(round(min(dat$BOTDEPTH,na.rm=T)),round(max(dat$BOTDEPTH,na.rm=T)),by=10)/10)
        NEWDATA=data.frame(year=factor(Yr.ref,levels(dat$year)),
                           log.BOTDEPTH=log(Z.seq),
                           log.Mid.Lat=mean(dat$log.Mid.Lat,na.rm=T),
                           log.Effort=mean(dat$log.Effort,na.rm=T))
        MOD=STORE.BOOT[[i]]
        An.Ktch=matrix(nrow=length(MOD),ncol=length(NEWDATA$log.BOTDEPTH))
        colnames(An.Ktch)=exp(NEWDATA$log.BOTDEPTH)
        An.Ktch.count=An.Ktch.zero=An.Ktch
        PRd.knt=PRd.zro="NO"  #select species for which lat is significant
        if("log.BOTDEPTH"%in%Sig.terms[[i]]$Sigi.count) PRd.knt="YES"
        if("log.BOTDEPTH"%in%Sig.terms[[i]]$Sigi.zero) PRd.zro="YES"
        
        if(PRd.knt=="YES")
        {
          for (x in 1:length(MOD))An.Ktch.count[x,]=Pred.Model.boot_count_zero(MOD[[x]]$Fit,NEWDATA,PRd.knt,PRd.zro)$PREDS.count
          PRD=data.frame(BOTDEPTH=exp(NEWDATA$log.BOTDEPTH),
                         MEAN=apply(An.Ktch.count, 2, function(x) median(x, na.rm=T)),
                         SD=apply(An.Ktch.count,2,function(x) sd(x,na.rm=T)),
                         LOW=apply(An.Ktch.count, 2, function(x) quantile(x, 0.025,na.rm=T)),
                         UP=apply(An.Ktch.count, 2, function(x) quantile(x, 0.975,na.rm=T)))
          PRD$CV=PRD$SD/PRD$MEAN*100
          PRED.z.count[[i]]=PRD
        }
        if(PRd.zro=="YES")
        {
          for (x in 1:length(MOD))An.Ktch.zero[x,]=Pred.Model.boot_count_zero(MOD[[x]]$Fit,NEWDATA,PRd.knt,PRd.zro)$PREDS.zero
          PRD=data.frame(BOTDEPTH=exp(NEWDATA$log.BOTDEPTH),
                         MEAN=apply(An.Ktch.zero, 2, function(x) median(x, na.rm=T)),
                         SD=apply(An.Ktch.zero,2,function(x) sd(x,na.rm=T)),
                         LOW=apply(An.Ktch.zero, 2, function(x) quantile(x, 0.025,na.rm=T)),
                         UP=apply(An.Ktch.zero, 2, function(x) quantile(x, 0.975,na.rm=T)))
          PRD$CV=PRD$SD/PRD$MEAN*100
          PRED.z.zero[[i]]=PRD
        }
        rm(dat,MOD)
      }
  }
  })
  
  
  #1.12.6. Plot standardised cpue
  fn.relative=function(D)
  {
    Mn=mean(D$MEAN)
    D$LOW=D$LOW/Mn
    D$UP=D$UP/Mn
    D$MEAN=D$MEAN/Mn
    return(D)
  }
  CI.fun=function(YR,LOW1,UP1,Colr,Colr2)
  {
    Year.Vec <- c(YR, tail(YR, 1), rev(YR), YR[1]) 
    Biom.Vec <- c(LOW1, tail(UP1, 1), rev(UP1), LOW1[1])
    polygon(Year.Vec, Biom.Vec, col = Colr, border = Colr2)
  }
  fun.plot.yr.pred=function(PRD,X,normalised,REV,n.seq,YLIM,XLIM,Type)
  {
    THIS=match(X,names(PRD))
    if(X=="BOTDEPTH")
    { 
      dd=subset(DATA.list[[i]],Catch.Target>0)
      that=match(X,names(dd))
      X.max=max(dd[,that],na.rm=T)
      PRD=PRD[PRD[,THIS]<=X.max,]
    }
    if(REV=="YES")
    {
      PRD[,THIS]=-PRD[,THIS]
      PRD=PRD[order(PRD[,THIS]),]
    }
    #standardise to a mean score of 1
    if(normalised=="YES") PRD=fn.relative(PRD)
    yr=as.numeric(as.character(PRD[,THIS]))
    ALL.yrs=seq(yr[1],yr[length(yr)])
    MeAn=PRD$MEAN
    UppCI=PRD$UP
    LowCI=PRD$LOW
    dat.plt=data.frame(yr=yr,MeAn=MeAn,UppCI=UppCI,LowCI=LowCI,CV=PRD$CV)
    #add missing years
    if(X=="year")
    {
      mis.yr=ALL.yrs[which(!ALL.yrs%in%yr)]
      if(length(mis.yr)>0)
      {
        dummy=dat.plt[1:length(mis.yr),]
        dummy[,]=NA
        dummy$yr=mis.yr
        dat.plt=rbind(dat.plt,dummy)
      }
    }
    if(X=="year")dat.plt=dat.plt[order(dat.plt$yr),]
    if(is.null(YLIM))YLIM=c(0,max(dat.plt$UppCI,na.rm=T))
    if(YLIM[2]>1e5) YLIM[2]=quantile(dat.plt$UppCI,0.5)
    if(is.null(XLIM)) XLIM=range(dat.plt[,1])
    if(Type=="points")
    {
      with(dat.plt,plot(yr,MeAn,pch=19,main="",xlab="",ylab="",
                        cex=1.25,xaxt="n",cex.axis=1.25,ylim=YLIM,xlim=XLIM))
      suppressWarnings(with(dat.plt,arrows(x0=yr, y0=LowCI, x1=yr, y1=UppCI,code = 3,angle=90,length=.025)))
    }
    if(Type=="polygon")
    {
      with(dat.plt,plot(yr,MeAn,type='l',main="",xlab="",ylab="",
                        lwd=2,xaxt="n",cex.axis=1.25,ylim=YLIM,xlim=XLIM))
      with(dat.plt,CI.fun(yr,UppCI,LowCI,"grey80","transparent"))
      with(dat.plt,lines(yr,MeAn,type='l',lwd=2))
    }
    axis(1,seq(XLIM[1],XLIM[2],n.seq),seq(XLIM[1],XLIM[2],n.seq),tck=-0.04,cex.axis=1.25)
    if(n.seq<10)axis(1,seq(XLIM[1],XLIM[2],1),F,tck=-0.02,cex.axis=1.25)
    #with(dat.plt,axis(1,yr,F,tck=-0.025))
    #with(dat.plt,axis(1,seq(yr[1],yr[length(yr)],n.seq),F,tck=-0.05))
    #with(dat.plt,axis(1,seq(yr[1],yr[length(yr)],n.seq),seq(yr[1],yr[length(yr)],n.seq),tck=-0.05,cex.axis=1.25))
    return(dat.plt)
  }
 
  
  #1.13.   Trends in Sex ratio           
  if(do.sex.ratio=="YES")
  {
    #fit binomial model
    SEX.fun=function(SPEC)
    {
      dat=subset(DATA,SPECIES==SPEC & SEX%in%c("F","M") & year%in%YEAR & BOTDEPTH <MaxDepth &
                   Month%in%These.Months & Set.time <"08:00" & !(year==2004))
      dat$Sex.bin=with(dat,ifelse(SEX=="M",1,0))
      dat$year=as.factor(dat$year)
      
      model<- glm(Sex.bin~year+BOTDEPTH+Mid.Lat, data=dat, family="binomial", maxit=500)
      Signifcance=anova(model,test="Chisq")
      
      #Predictions
      Yr=sort(unique(dat$year))
      #years
      NEWDATA=expand.grid(year=Yr,BOTDEPTH=mean(dat$BOTDEPTH,na.rm=T),Mid.Lat=mean(dat$Mid.Lat))
      YEAR.preds=predict(model,newdata=NEWDATA, type='response')
      
      #latitude
      LAT.seq=range(dat$Mid.Lat)
      LAT.seq=seq(round(LAT.seq[1]),round(LAT.seq[2]),.5)
      NEWDATA=expand.grid(year=factor(Yr[1],levels=levels(Yr)),BOTDEPTH=mean(dat$BOTDEPTH,na.rm=T),Mid.Lat=LAT.seq)
      LAT.preds=predict(model,newdata=NEWDATA, type='response',se.fit=T)
      
      #Depth
      Z.seq=range(dat$BOTDEPTH,na.rm=T)
      Z.seq=round(Z.seq/10)*10
      Z.seq=seq(round(Z.seq[1]),round(Z.seq[2]),10)
      NEWDATA=expand.grid(year=factor(Yr[1],levels=levels(Yr)),BOTDEPTH=Z.seq,Mid.Lat=mean(dat$Mid.Lat))
      Z.preds=predict(model,newdata=NEWDATA, type='response',se.fit=T)
      
      return(list(model=model,Signifcance=Signifcance,
                  YEAR=Yr,YEAR.preds=YEAR.preds,
                  LAT.seq=LAT.seq,LAT.preds=LAT.preds,
                  Z.seq=Z.seq,Z.preds=Z.preds,n=nrow(dat)))
    }
    #keep species with at least 10 males and 10 females
    TARGETS.sex=TARGETS[-match(c("WR","GR"),TARGETS)]
    N.species.sex=length(TARGETS.sex)
    
    Store.SEX=vector('list',N.species.sex)
    names(Store.SEX)=names(TARGETS.sex)
    for (i in 1:N.species.sex) Store.SEX[[i]]=SEX.fun(TARGETS.sex[i])
    
    #Predictions and CI thru Bootstrapping approach
    if(do.boot=="YES")
    {
      STORE.BOOT.sex=vector('list',N.species.sex)
      names(STORE.BOOT.sex)=names(TARGETS.sex)
      #STORE.BOOT.sex.all.stations=STORE.BOOT.sex
      
      #fit models to bootstrapped data
      system.time(for(i in 1:N.species.sex)  #takes 0.7 secs per iteration
      {
        Store.boot=Store.boot.all.stations=vector('list',n.boot)
        for(k in 1:n.boot)
        {
          Sex.fun.boot=function(SPEC)
          {
            dat=subset(DAT,SPECIES==SPEC & SEX%in%c("F","M") & year%in%YEAR & BOTDEPTH <MaxDepth &
                         Month%in%These.Months & Set.time <"08:00" & !(year==2004))
            dat$Sex.bin=with(dat,ifelse(SEX=="M",1,0))
            dat$year=as.factor(dat$year)
            model<- glm(Sex.bin~year+BOTDEPTH+Mid.Lat, data=dat, family="binomial", maxit=500)
            
            
            #Predictions
            Yr=factor(levels(dat$year),levels(dat$year))
            
            #years
            NEWDATA=data.frame(year=Yr,
                               SEX=factor("F",levels(dat$SEX)),
                               BOTDEPTH=mean(dat$BOTDEPTH,na.rm=T),
                               Mid.Lat=mean(dat$Mid.Lat,na.rm=T))
            YEAR.preds=predict(model,newdata=NEWDATA, type='response')
            
            #latitude
            LAT.seq=range(dat$Mid.Lat)
            LAT.seq=seq(round(LAT.seq[1]),round(LAT.seq[2]),.5)
            Yr.ref=names(sort(table(dat$year)))
            Yr.ref=Yr.ref[length(Yr.ref)]
            NEWDATA=data.frame(year=factor(Yr.ref,levels(dat$year)),
                               SEX=factor("F",levels(dat$SEX)),
                               BOTDEPTH=mean(dat$BOTDEPTH,na.rm=T),
                               Mid.Lat=LAT.seq)
            LAT.preds=predict(model,newdata=NEWDATA, type='response')
            
            #Depth
            Z.seq=range(dat$BOTDEPTH,na.rm=T)
            Z.seq=round(Z.seq/10)*10
            Z.seq=seq(round(Z.seq[1]),round(Z.seq[2]),10)
            NEWDATA=data.frame(year=factor(Yr.ref,levels(dat$year)),
                               SEX=factor("F",levels(dat$SEX)),
                               BOTDEPTH=Z.seq,
                               Mid.Lat=mean(dat$Mid.Lat,na.rm=T))
            Z.preds=predict(model,newdata=NEWDATA, type='response')
            
            return(list(YEAR=Yr,YEAR.preds=YEAR.preds,
                        LAT.seq=LAT.seq,LAT.preds=LAT.preds,
                        Z.seq=Z.seq,Z.preds=Z.preds))
          }
          
          #fixed stations
          DAT=fn.boot(subset(DATA,FixedStation=="YES"))
          mod=try(Sex.fun.boot(TARGETS.sex[i]),TRUE)  #deal with dodgy boot data
          if(isTRUE(class(mod)=="try-error")) { next } else { Store.boot[[k]]=mod } 
          
          
          #All stations
          # DAT=fn.boot(DATA) 
          # mod=try(Sex.fun.boot(TARGETS.sex[i]),TRUE)  #deal with dodgy boot data
          # if(isTRUE(class(mod)=="try-error")) { next } else { Store.boot.all.stations[[k]]=mod } 
          
        }
        STORE.BOOT.sex[[i]]=Store.boot
        #       STORE.BOOT.sex.all.stations[[i]]=Store.boot.all.stations
      })
    }
    
  }
  
  
  #1.14.   Trends in size
  #fit gaussian model to Fixed stations
  TARGETS.size=TARGETS
  N.species.size=length(TARGETS.size)
  BEST.model.size=BEST.model
  for(i in 1:N.species.size) BEST.model.size[[i]]=formula(paste("FL",paste(c("year",
                                  's(Mid.Lat,k=3)','s(BOTDEPTH,k=3)'),collapse="+"),sep="~"))
  BEST.model.size$`Milk shark`=BEST.model.size$`Spot-tail shark`=BEST.model.size$`Sliteye shark`=
      formula(paste("FL",paste(c("year",'s(Mid.Lat,k=3)'),collapse="+"),sep="~"))
  BEST.model.size$`Scalloped hammerhead`=formula(paste("FL",paste(c("year",
                                  's(BOTDEPTH,k=3)'),collapse="+"),sep="~"))
  Size.fun=function(dat,FORMULA)
  {
    if(do.GLM=="YES")
    {
      model<- glm(FL~year+log.BOTDEPTH+log.Mid.Lat, data=dat,family=gaussian,maxit=500)
      #Signifcance=anova(model,test="Chisq")
      
      #Predictions
      year.pred=summary(emmeans(model,"year", type="response"))
      Lata=range(floor(dat$Mid.Lat))
      Lat.pred=summary(emmeans(model,"log.Mid.Lat", type="response",at=list(log.Mid.Lat=log(seq(Lata[1],Lata[2],.1)))))
      ZZ=range(10*floor(dat$BOTDEPTH/10))
      Depth.pred=summary(emmeans(model,"log.BOTDEPTH", type="response",at=list(log.BOTDEPTH=log(seq(ZZ[1],ZZ[2],10)))))
    }
    if(do.GAM=="YES")
    {
      model<- gam(FORMULA, data=dat,method = "REML",family=gaussian)
      #Signifcance=anova(model)
      
      #Predictions
      year.pred=summary(emmeans(model,"year", type="response"))
      Lat.pred=NULL
      used.term=grepl("Mid.Lat", attr(terms(FORMULA),'term.labels'))
      if(sum(used.term)>0)
      {
        Lata=range(floor(dat$Mid.Lat))
        Lat.pred=summary(emmeans(model,"Mid.Lat", type="response",at=list(Mid.Lat=seq(Lata[1],Lata[2],.1))))
      }
      Depth.pred=NULL
      used.term=grepl("BOTDEPTH", attr(terms(FORMULA),'term.labels'))
      if(sum(used.term)>0)
      {
        ZZ=range(10*floor(dat$BOTDEPTH/10))
        Depth.pred=summary(emmeans(model,"BOTDEPTH", type="response",at=list(BOTDEPTH=seq(ZZ[1],ZZ[2],10))))
      }
    }
    
    return(list(Fit=model,DAT=dat,year.pred=year.pred,Lat.pred=Lat.pred,Depth.pred=Depth.pred))
  }
  
  Store.size=vector('list',N.species.size)
  names(Store.size)=names(TARGETS.size)
  for (i in 1:N.species.size)
  {
    dat=subset(DATA,SPECIES==TARGETS.size[i] & !is.na(FL) & year%in%YEAR 
               & BOTDEPTH <210 & Month%in%These.Months 
               & Set.time <"08:00" & FixedStation=="YES")
    dat=subset(dat,!year==2004)
    
    yr.tabl=table(dat$year)
    
    #at least 5 observations per year
    dat=subset(dat,year%in%as.numeric(names(yr.tabl[yr.tabl>=5])))
    
    dat$year=factor(dat$year,levels=sort(unique(dat$year)))
    dat$SEX=as.factor(dat$SEX)
    
    dat$Mid.Lat=abs(dat$Mid.Lat)
    dat$log.BOTDEPTH=log(dat$BOTDEPTH)
    dat$log.Mid.Lat=log(dat$Mid.Lat)
    
    Store.size[[i]]=Size.fun(dat,FORMULA=BEST.model.size[[i]])
  }
    
  #Goodnes of fit 
  Pos.Diag.fn=function(MODEL,M=1)   #function for fit diagnostics
  {
    RES=MODEL$residuals   #residuals
    Std.RES=RES/sd(RES)   #standardised residuals (res/SD(res))
    PRED=predict(MODEL)
    
    qqnorm(RES,main="",ylab="",xlab="")
    qqline(RES, col = 'grey40',lwd=1.5,lty=2)
    #mtext(SPECIES,3,outer=F,line=0.25,cex=1.3)
    if(i==4) mtext("Residuals",2,outer=F,line=1.3,las=3,cex=M)
    # mtext("                        Theoretical quantiles",1,outer=F,line=1.5,cex=M)
    
    hist(Std.RES,xlim=c(-5,5),ylab="",xlab="",main="",col="grey",breaks=50)
    box()
    if(i==4)  mtext("Frequency",2,outer=F,line=1.3,las=3,cex=M)
    #mtext("                      Standardised residuals",1,outer=F,line=1.5,cex=M)
    
    plot(PRED,Std.RES,ylab="",xlab="",ylim=c(-5,5))
    abline(0,0,lwd=1.5,lty=2,col='grey40')
    if(i==4)  mtext("Standardised residuals",2,outer=F,line=1.3,las=3,cex=M)
    #mtext("                         Fitted values",1,outer=F,line=1.5,cex=M)
    
  }
  hndl.fit.size="C:/Matias/Analyses/Surveys/Naturaliste_longline/outputs/Size/"
  fn.fig(paste(hndl.fit.size,"fit.diagnostics",sep=""),1600,2400)
  par(mfrow=c(8,3),mar=c(2,2,1.5,1),oma=c(1,1.5,.1,2),mgp=c(2,.7,0))
  for(i in 1:N.species.size)
  {
    Pos.Diag.fn(MODEL=Store.size[[i]]$Fit,M=1)
    mtext(names(TARGETS)[i],4,las=3,line=1,cex=.6)
  }
  dev.off()

  #Predict year, lat and depth effects
  PRED.size=vector('list',length(Species.cpue))
  names(PRED.size)=names(DATA.list)
  PRED.lat.size=PRED.z.size=PRED.size
  for(i in Species.cpue)
  {
    PRED.size[[i]]=Store.size[[i]]$year.pred
    if(!is.null(Store.size[[i]]$Lat.pred))PRED.lat.size[[i]]=Store.size[[i]]$Lat.pred
    if(!is.null(Store.size[[i]]$Depth.pred))PRED.z.size[[i]]=Store.size[[i]]$Depth.pred
  }
}

# 2. Multivariate analysis and Ecosystem indicators   MISSING:  USE Fixed Stations only! ----------------------------------------------------------------------
if(Do.multivariate=="YES")
{
  hndl.eco="C:/Matias/Analyses/Surveys/Naturaliste_longline/outputs/Ecosystems/"
  
  DATA.eco=subset(DATA,Month%in%These.Months & Taxa=="Elasmobranch" & BOTDEPTH<500 & FixedStation=="YES" )
  DATA.eco$SPECIES=as.factor(DATA.eco$SPECIES)
  
  # #for consistency, use only shots used for standardisation
  # if(exists("DATA.list"))      #DATA.list doesn't exist if abundance analysis not executed
  # {
  #   d=do.call(rbind,DATA.list)    
  #   d=subset(d,FixedStation=="YES")
  #   d=unique(d$SHEET_NO)
  #   DATA.eco=subset(DATA.eco,SHEET_NO%in%d)
  # }
  
  #some inputs
  ResVar="Number"
  MultiVar="SPECIES"
  IDVAR=c("year","BOTDEPTH.mean","Mid.Lat.mean","Effort")   #removed Long due to high correlation and Month as only 5-8 months used
  Predictors=IDVAR[-match(c("Effort"),IDVAR)]
  num.vars=c("Number","BOTDEPTH","Mid.Lat","Mid.Long","N.hooks","SOAK.TIME")
  Formula=formula("d.res.var~.")
  Formula.mvglm=formula("y~.")
  #Formula=formula("d.res.var~year*Mid.Lat.mean+BOTDEPTH.mean")
  Keep.vars=c("SHEET_NO","SPECIES","COMMON_NAME","SCIENTIFIC_NAME","Station.no.","date","Day",
              "Month","year","TL","FL","SEX","Number","BOTDEPTH","Mid.Lat","Mid.Long","N.hooks","N.hooks.Fixed",
              "SOAK.TIME","Moon","TEMP","Lat.round","Long.round")
  DATA.eco=DATA.eco[,match(Keep.vars,names(DATA.eco))]
  
  #remove additional stations as not sampled extensively
  DATA.eco=subset(DATA.eco,Station.no.%in%
                    as.character(Fixed.Stations$Station.no.[which(!Fixed.Stations$Station.no.=="additional")])) 
  DATA.eco$Station.no.=as.numeric(DATA.eco$Station.no.)
  DATA.eco$Statn_yr=with(DATA.eco,paste(Station.no.,year))
  
  DATA.eco$SPECIES=as.character(DATA.eco$SPECIES)
  DATA.eco$Month=as.factor(DATA.eco$Month)
  DATA.eco$year=as.factor(DATA.eco$year)
  
  for(n in 1:length(num.vars)) DATA.eco[,match(num.vars[n],names(DATA.eco))]=as.numeric(as.character(DATA.eco[,match(num.vars[n],names(DATA.eco))]))
  DATA.eco$Effort=with(DATA.eco,N.hooks*SOAK.TIME)
  
  #Aggreate shots by station
  #note: aggregate by Station-year
  Table.stations.per.year=with(DATA.eco[!duplicated(DATA.eco$SHEET_NO),],table(year,Station.no.))  #unbalanced design    
  
  #remove stations 2 and 3 as hardly visited
  DATA.eco=subset(DATA.eco,!Station.no.%in%c(2,3))
  
  #mean depth and lat by station-year
  Mean.Mid.Lat.Long.depth=aggregate(cbind(Mid.Lat,Mid.Long,BOTDEPTH)~Station.no.+year,DATA.eco,mean)   
  names(Mean.Mid.Lat.Long.depth)[match(c("Mid.Lat","Mid.Long","BOTDEPTH"),names(Mean.Mid.Lat.Long.depth))]=c("Mid.Lat.mean","Mid.Long.mean","BOTDEPTH.mean")
  DATA.eco=merge(DATA.eco,Mean.Mid.Lat.Long.depth,by=c("Station.no.","year"))
  
  #effort by station-year
  Effort.Station.year=aggregate(Effort~Station.no.+year,DATA.eco[!duplicated(DATA.eco$SHEET_NO),],sum)
  names(Effort.Station.year)[match("Effort",names(Effort.Station.year))]="Effort.statn_year"
  DATA.eco=merge(DATA.eco,Effort.Station.year,by=c("Station.no.","year"))  
  
  #round mean position and depth
  DATA.eco$Mid.Long.mean=round(DATA.eco$Mid.Long.mean,1)
  DATA.eco$Mid.Lat.mean=round(abs(DATA.eco$Mid.Lat.mean),1)
  DATA.eco$BOTDEPTH.mean=round(DATA.eco$BOTDEPTH.mean)
  
  #numbers by station-year
  Numbers.Station.year=aggregate(Number~Station.no.+year+SPECIES+COMMON_NAME+SCIENTIFIC_NAME+
                                   +Mid.Lat.mean+Mid.Long.mean+BOTDEPTH.mean+Effort.statn_year,DATA.eco,sum)
  
  # #keep species occurring in >0.1%
  # Tab.sp=table(Numbers.Station.year$COMMON_NAME)
  # Tab.sp=Tab.sp/sum(Tab.sp)
  # Tab.sp=round(100*Tab.sp,2)
  # Tab.sp=Tab.sp[Tab.sp>=0.1]
  # Numbers.Station.year=subset(Numbers.Station.year,COMMON_NAME%in%names(Tab.sp))
  
  #See spatial cpue of species ( in numbers per 1000 hooks)
  Numbers.Station.year$cpue=1000*with(Numbers.Station.year,Number/Effort.statn_year)
  
  Dist.sptl=aggregate(cpue~Mid.Lat.mean+Mid.Long.mean+COMMON_NAME,Numbers.Station.year,sum)
  Unik.sp=unique(Dist.sptl$COMMON_NAME)
  smart.par(n.plots=length(Unik.sp),MAR=c(2,2,1,1),OMA=c(1,1.5,.1,.1),MGP=c(2.5,.7,0))
  for(sp in 1:length(Unik.sp)) with(subset(Dist.sptl,COMMON_NAME==Unik.sp[sp]),plot(Mid.Long.mean,Mid.Lat.mean,
                                                                                    main=Unik.sp[sp],pch=19,col=2,cex.main=.85,cex=cpue/max(cpue)*3,ylab="",xlab="",
                                                                                    ylim=range(Dist.sptl$Mid.Lat.mean),xlim=range(Dist.sptl$Mid.Long.mean)))
  
  #2.1 Multivariate
  #2.1. Classic
  # review: unbalanced design (some years only  Broome, some only Ningaloo)
  #         grouping by station-year, should drop stations 2 and 3?
  #         can have only factors or also covariate as predictors? how to show lat effect then?
  #         should have blocks rathen than lat? if so, then depth nested in block?
  if(Do.multivariate=="YES")   
  {
    source("C:/Matias/Analyses/SOURCE_SCRIPTS/Multivariate_statistics.R")
    DataSets=c("catch","cpue","proportion")   #response variables
    IDVAR=c("year","BOTDEPTH.mean","Mid.Lat.mean","Effort.statn_year")
    STore.multi.var.trad=Multivar.fn(DATA=Numbers.Station.year,ResVar=ResVar,MultiVar=MultiVar,
                                     Predictors=Predictors,IDVAR=IDVAR,
                                     Formula=formula("d.res.var~."),DataSets=DataSets)
  }               
  
  #2.2 Multivariate glm
  if(Do.multivariate.glm=="YES")
  {
    #create data matrix
    fn.reshp=function(d,Y,TimeVar,IDVAR)
    {
      DATA.agg=aggregate(formula(paste(Y,paste(c(TimeVar,IDVAR),collapse="+"),sep="~")),d,sum)
      DATA.wide=reshape(DATA.agg,v.names=Y,idvar=IDVAR,timevar=TimeVar,direction="wide")
      DATA.wide[is.na(DATA.wide)]=0
      colnames(DATA.wide)=gsub(paste(Y,".",sep=""), "", names(DATA.wide))
      #props=DATA.wide
      #props[,-which(IDVAR%in%names(props))]=props[,-which(IDVAR%in%names(props))]/rowSums(props[,-which(IDVAR%in%names(props))])
      #return(list(catch=DATA.wide,proportion=props))
      return(list(DATA.wide))
    }
    Dat=fn.reshp(d=Numbers.Station.year,Y=ResVar,TimeVar=MultiVar,IDVAR=IDVAR)
    mod.fn=function(X) factor(names(rev(sort(table(X)))[1]),levels=levels(X))
    fun.percent.dev.exp=function(null,modl) 100*(abs(null-modl))/null
    fn.multivar=function(d,FORMULA.mvglm,OFFSET)
    {
      library(mvabund)
      d.res.var=d[,-which(IDVAR%in%names(d))]
      d.preds=d[,c(Predictors)]
      
      y=mvabund(d.res.var)
      
      null <- manyglm(y~1, family="negative.binomial")
      full <- manyglm(FORMULA.mvglm, family="negative.binomial", data=d.preds)
      #plot(full)
      
      ANOVA=anova(full, nBoot=199, test="wald")
      Percen.dev.exp.modl=fun.percent.dev.exp(null$deviance,full$deviance)
      
      fn.pred=function(what,dd)
      {
        prdktr=dd[,match(what,names(dd))]
        if(is.factor(prdktr)) prdktr=factor(unique(prdktr),levels(prdktr)) else
          prdktr=seq(min(prdktr),max(prdktr),length.out=20)
        
        others=dd[,-match(what,names(dd))]
        for(cl in 1:ncol(others))
        {
          
        }
        NEW=data.frame(
          year=factor("2009",levels(dd$year)),
          BOTDEPTH.mean=seq(min(dd$BOTDEPTH.mean),max(dd$BOTDEPTH.mean),length.out=20),
          Mid.Lat.mean=mean(dd$Mid.Lat.mean),
          Effort=mean(dd$Effort))
        
        
        PRED=predict(full,newdata=NEW,type='response', se.fit =T)
        
        Stand=apply(PRED$fit,2,max)
        PRED.dot=t(t(PRED$fit)/Stand)
        PRED.dot.min2SE=t(t(PRED$fit-2*PRED$se.fit)/Stand)
        PRED.dot.plus2SE=t(t(PRED$fit+2+PRED$se.fit)/Stand)
        plot(0,type='n',axes=FALSE,ann=FALSE,
             ylim=c(1,ncol(PRED.dot)),xlim=c(1,length(unique(NEW$year))))
        for(n in 1:nrow(NEW))
        {
          x=NEW[n,]
          y=PRED.dot[n,]
          #Mean
          points(rep(x$year,length(y)),1:length(y),cex=PRED.dot[n,]*3,col="black")
          #CI
          points(rep(x$year,length(y)),1:length(y),cex=PRED.dot.min2SE[n,]*3,col="grey60") 
          points(rep(x$year,length(y)),1:length(y),cex=PRED.dot.plus2SE[n,]*3,col="grey60") 
        }
        axis(2,at=1:ncol(PRED.dot),colnames(PRED.dot),las=1,cex.axis=.85)
        axis(1,at=1:length(unique(NEW$year)),unique(NEW$year),las=1,cex.axis=.85)
        box()
      }
      fn.pred(what=names(d.preds)[p],dd=d.preds)
      
      
      
      return(list(ANOVA=ANOVA,Percen.dev.exp.modl=Percen.dev.exp.modl))
      
    }
    STore.multi.var=Dat
    for(s in 1:length(STore.multi.var)) STore.multi.var[[s]]=fn.multivar(d=Dat[[s]],FORMULA.mvglm=Formula.mvglm)  
  }
  
  #2.2 Ecosystems indicators
  if(Do.ecosystems=="YES")
  {
    require(vegan)
    require(BEQI2)
    require(asbio)
    
    #Function for calculating Diversity indices and Ecosystem indicators
    source("C:/Matias/Analyses/Ecosystem indices/Shark-bycatch/Git_bycatch_TDGDLF/Ecosystem_functions.R")
    
    #add trophic level
    TL=read.csv("C:/Matias/Analyses/Ecosystem indices/Shark-bycatch/SPECIES+PCS+FATE.csv",stringsAsFactors=F)
    DATA.eco=merge(DATA.eco,subset(TL,select=c(SPECIES,TROPHIC_LEVEL)),by="SPECIES",all.x=T)
    
    #balance desing
    #plot shots by year     
    yr.dummy=sort(unique(DATA.eco$year))
    fn.fig(paste(hndl.eco,"Shots by year",sep=""),2000,2400)
    smart.par(n.plots=length(yr.dummy),MAR=c(2,2,1,1),OMA=c(1,1.5,.1,.1),MGP=c(2.5,.7,0))
    for(y in 1:length(yr.dummy))
    {
      a=subset(DATA.eco,year==yr.dummy[y])
      plot(a$Mid.Long,a$Mid.Lat,pch=19,ylab="",xlab="",main=yr.dummy[y],
           xlim=c(min(DATA.eco$Mid.Long),max(DATA.eco$Mid.Long)),ylim=c(min(DATA.eco$Mid.Lat),max(DATA.eco$Mid.Lat)))
      abline(v=117,lwd=2.5,col="grey70")
    }
    mtext("Latitude (S)",side=2,line=-0.35,outer=T,cex=1.5,las=3)
    mtext("Longitude (E)",side=1,line=1,outer=T,cex=1.5)
    dev.off()
    
    DATA.eco$Leg=with(DATA.eco,ifelse(Mid.Long<=117,1,ifelse(Mid.Long>117,2,NA)))
    Table.yr.leg=table(DATA.eco$Leg,DATA.eco$year,useNA='ifany')
    
    if(is.character(DATA.eco$SPECIES))DATA.eco$SPECIES=factor(DATA.eco$SPECIES)
    
    #calculate ecosystem indicators  on cpue      
    gp.by="Statn_yr"
    #gp.by="SHEET_NO"
    
    if(gp.by=="Statn_yr")
    {
      fn.reshp=function(d,Y,TimeVar,IdVAR)
      {
        DATA.wide.cpue=reshape(d[,match(c("cpue",IdVAR,TimeVar),names(d))],v.names="cpue",
                               idvar=IdVAR,timevar=TimeVar,direction="wide")
        DATA.wide.cpue[is.na(DATA.wide.cpue)]=0
        colnames(DATA.wide.cpue)=gsub(paste("cpue",".",sep=""), "", names(DATA.wide.cpue))
        return(DATA.wide.cpue)
      }
      Dat=fn.reshp(d=Numbers.Station.year,Y=ResVar,TimeVar=MultiVar,IdVAR=c("Station.no.",Predictors)) 
      DATA.shots.diversity=Dat[,match(c("Station.no.",Predictors),names(Dat))]
      Dat.y=Dat[,-match(c("Station.no.",Predictors),names(Dat))]
      
      #diversity indices
      DATA.shots.diversity$Shannon = diversity(Dat.y, index = "shannon")
      DATA.shots.diversity$Simpson = diversity(Dat.y, index = "simpson")
      sp_names=colnames(Dat.y)
      Evenness = function(H, S) {J = H / log(S)}  #where H = Shannon index for each community; S = Number of species
      DATA.shots.diversity$Pielou = Evenness(H = DATA.shots.diversity$Shannon, S = length(sp_names))
      DATA.shots.diversity$Pielou[DATA.shots.diversity$Pielou == 0] = NA
      
      #ecosystems indicators
      MTL=aggregate(TROPHIC_LEVEL~Station.no.+year,DATA.eco,mean,na.rm=T) 
      names(MTL)[match("TROPHIC_LEVEL",names(MTL))]="MTL"
      MML=aggregate(FL~Station.no.+year,DATA.eco,mean,na.rm=T)
      names(MML)[match("FL",names(MML))]="MML"
      DATA.shots.diversity=merge(DATA.shots.diversity,MTL,by=c("Station.no.","year"),all=T)
      DATA.shots.diversity=merge(DATA.shots.diversity,MML,by=c("Station.no.","year"),all=T)
      
    }
    if(gp.by=="SHEET_NO")
    {
      SHOTS=unique(DATA.eco$SHEET_NO)   #by Shot
      Store.shots=vector('list',length(SHOTS))
      names(Store.shots)=SHOTS
      system.time(for(i in 1:length(SHOTS))   #takes 12 seconds
      {
        dat=subset(DATA.eco,SHEET_NO==SHOTS[i])
        d=table(dat$SPECIES)
        
        #diversity indices
        Div.InDX=Div.ind.shot(data=d)
        Div.InDX=as.data.frame(do.call(cbind,Div.InDX))
        
        dat.indx=dat[!duplicated(dat$SHEET_NO),]
        dat.indx=cbind(dat.indx,Div.InDX)
        
        dat.indx$EFFORT=with(dat.indx,N.hooks.Fixed*SOAK.TIME)
        
        dat.indx=subset(dat.indx,select=c(SHEET_NO,Leg,year,Month,Mid.Lat,Mid.Long,BOTDEPTH,
                                          SOAK.TIME,N.hooks.Fixed,EFFORT,Shannon,Margalef,Pielou,Simpson))
        
        #ecosystem indicators 
        dd=d
        dd=subset(dd,dd>0)
        t_level=data.frame(SPECIES=character(),TROPHIC_LEVEL=numeric())
        ddat=subset(dat,!is.na(TROPHIC_LEVEL))
        if(nrow(ddat)>0)t_level=aggregate(TROPHIC_LEVEL ~ SPECIES, FUN=mean, data=ddat)
        dummy=subset(dat,!is.na(FL))
        dummy=subset(dummy,FL>0)
        fl_max=NULL
        if(nrow(dummy)>0)  fl_max=aggregate(FL~SPECIES,FUN=max,data=dummy,na.rm=T) 
        Eco.InDX=Eco.ind.shot(data=dd,TROPHIC.LEVEL=t_level,MAX.BODY.LENGTH=fl_max)
        Eco.InDX=as.data.frame(do.call(cbind,Eco.InDX))
        dat.indx=cbind(dat.indx,Eco.InDX)
        Store.shots[[i]]=dat.indx
      })
      DATA.shots.diversity=do.call(rbind,Store.shots)
      
    }
    
    #remove shots with no or single species as indices cannot be calculated
    DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(Shannon))  
    DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(Pielou))
    DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(Simpson))
    DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(MTL))
    DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(MML))
    DATA.shots.diversity=subset(DATA.shots.diversity,Simpson>0)
    DATA.shots.diversity=subset(DATA.shots.diversity,MML>0)
    #DATA.shots.diversity=subset(DATA.shots.diversity,!year==2004)
    
    resp.vars=c("Shannon","Simpson","MTL","MML")
    n.rv=length(resp.vars)
    
    #nominal indices
    fn.nminl=function(d,RESVAR,PRED)
    {
      # a=aggregate(DATA.shots.diversity[,match(RESVAR,names(DATA.shots.diversity))]~year,DATA.shots.diversity,mean)
      # a.sd=aggregate(DATA.shots.diversity[,match(RESVAR,names(DATA.shots.diversity))]~year,DATA.shots.diversity,sd)
      # plot(a[,1],a[,2],pch=19,main=RESVAR)
      # segments(a[,1],a[,2]-a.sd[,2],a[,1],a[,2]+a.sd[,2])
      
      d$dummy=DATA.shots.diversity[,match(PRED,names(DATA.shots.diversity))]
      if(!is.factor(d$dummy)) d$dummy=cut(d$dummy,20)
      boxplot(d[,match(RESVAR,names(d))]~d$dummy,main=RESVAR)
    }
    #year
    fn.fig(paste(hndl.eco,"bxplt_year",sep=""),2000,2400)
    smart.par(n.plots=length(resp.vars),MAR=c(2,2,1,1),OMA=c(1,1.5,.1,.1),MGP=c(2.5,.7,0))
    for(n in 1:n.rv) fn.nminl(d=DATA.shots.diversity,RESVAR=resp.vars[n],PRED="year")
    dev.off()
    
    #BOTDEPTH.mean
    fn.fig(paste(hndl.eco,"bxplt_depth",sep=""),2000,2400)
    smart.par(n.plots=length(resp.vars),MAR=c(2,2,1,1),OMA=c(1,1.5,.1,.1),MGP=c(2.5,.7,0))
    for(n in 1:n.rv) fn.nminl(d=DATA.shots.diversity,RESVAR=resp.vars[n],PRED="BOTDEPTH.mean")
    dev.off()
    
    #Mid.Lat.mean
    fn.fig(paste(hndl.eco,"bxplt_mid.lat",sep=""),2000,2400)
    smart.par(n.plots=length(resp.vars),MAR=c(2,2,1,1),OMA=c(1,1.5,.1,.1),MGP=c(2.5,.7,0))
    for(n in 1:n.rv) fn.nminl(d=DATA.shots.diversity,RESVAR=resp.vars[n],PRED="Mid.Lat.mean")
    dev.off()
    
    # #check if including effort as offset
    # fn.fig(paste(hndl.eco,"Diversity vs effort",sep=""),2000,2400)
    # par(mfcol=c(3,2),mai=c(.4,.45,.1,.1),mgp=c(2,.8,0))
    # for(n in 1:n.rv)
    # {
    #   x=DATA.shots.diversity$EFFORT
    #   y= DATA.shots.diversity[,match(resp.vars[n],names(DATA.shots.diversity))] 
    #   plot(x,y,ylab=resp.vars[n],xlab='Effort')
    #   MODL=lm(y~x)
    #   legend('topright',paste("slope=",round(coef(MODL)[2],3)),bty='n')
    # }
    # dev.off()
    
    # fn.fig(paste(hndl.eco,"Diversity vs hooks",sep=""),2000,2400)
    # par(mfcol=c(3,2),mai=c(.4,.45,.1,.1),mgp=c(2,.8,0))
    # for(n in 1:n.rv)
    # {
    #   x=DATA.shots.diversity$N.hooks.Fixed
    #   y= DATA.shots.diversity[,match(resp.vars[n],names(DATA.shots.diversity))] 
    #   plot(x,y,ylab=resp.vars[n],xlab='Effort')
    #   MODL=lm(y~x)
    #   legend('topright',paste("slope=",round(coef(MODL)[2],3)),bty='n')
    # }
    # dev.off()
    
    # fn.fig(paste(hndl.eco,"Diversity vs SOAK.TIME",sep=""),2000,2400)
    # par(mfcol=c(3,2),mai=c(.4,.45,.1,.1),mgp=c(2,.8,0))
    # for(n in 1:n.rv)
    # {
    #   x=DATA.shots.diversity$SOAK.TIME
    #   y= DATA.shots.diversity[,match(resp.vars[n],names(DATA.shots.diversity))] 
    #   plot(x,y,ylab=resp.vars[n],xlab='Effort')
    #   MODL=lm(y~x)
    #   legend('topright',paste("slope=",round(coef(MODL)[2],3)),bty='n')
    # }
    # dev.off()
    
    # fn.fig(paste(hndl.eco,"Effort vars vs time",sep=""),2000,2400)
    # par(mfcol=c(3,1),mai=c(.4,.45,.1,.1),mgp=c(2,.8,0))
    # y=DATA.shots.diversity$N.hooks.Fixed
    # x= DATA.shots.diversity$year 
    # boxplot(y~x,ylab="N.hooks.Fixed",xlab='year')
    # 
    # y=DATA.shots.diversity$SOAK.TIME
    # x= DATA.shots.diversity$year 
    # boxplot(y~x,ylab="SOAK.TIME",xlab='year')
    # 
    # y=DATA.shots.diversity$EFFORT
    # x= DATA.shots.diversity$year 
    # boxplot(y~x,ylab="EFFORT",xlab='year')
    # 
    # dev.off()
    
    
    #no linear increase in indices with effort, hence no effort offset. Also, gear snaps (the reason for which
    #   some shots have less hooks) shows no systematic change thru time so remove effort
    OFFSETT=NA
    #OFFSETT="offset(log.EFFORT)"
    
    
    #fit glm to ecosystem indicators     
    Mod.fn.glm=function(d,ResVar,Expl.vars,Predictrs,FactoRs,OFFSET,log.var)
    {
      d=d[,match(c(ResVar,Expl.vars),names(d))]
      #d=d[!is.na(d[,match(ResVar,names(d))]),]
      #d=d[d[,match(ResVar,names(d))]>0,]
      #d=subset(d,!is.na(EFFORT))
      # if(length(d$BOTDEPTH)>0)
      # {
      #   d=subset(d,!is.na(BOTDEPTH))
      #   d=subset(d,BOTDEPTH>=3)
      #   d$log.BOTDEPTH=log(d$BOTDEPTH)
      # }
      Y=d[,match(ResVar,names(d))]
      #d=subset(d,EFFORT>0)
      #d$log.EFFORT=log(d$EFFORT)
      #d$log.Mid.Lat=log(-d$Mid.Lat)
      
      #set factors
      #for(x in 1:length(FactoRs)) d[,match(FactoRs[x],names(d))]=as.factor(d[,match(FactoRs[x],names(d))])
      Pred.form=paste(c(Predictrs),collapse="+")
      if(!is.na(OFFSET))Pred.form=paste(c(Predictrs,OFFSET),collapse="+")
      if(log.var=="YES")Formula=as.formula(paste("log(",ResVar,")", "~", Pred.form,collapse=NULL))
      if(log.var=="NO")Formula=as.formula(paste(ResVar, "~", Pred.form,collapse=NULL))
      
      
      #Fit model 
      model <- glm(Formula, data=d)
      
      return(list(model=model, data=d))
    }
    
    
    #Predictors=c("year","Mid.Lat","BOTDEPTH")
    #Predictors=c("year","log.Mid.Lat","log.BOTDEPTH","log.SOAK.TIME","log.N.hooks.Fixed")
    Expl.varS=Predictors
    FactoRS="year"
    #Expl.varS=c("year","Mid.Lat","BOTDEPTH","EFFORT","SOAK.TIME","N.hooks.Fixed")
    #FactoRS=Expl.varS[-match(c("Mid.Lat","BOTDEPTH","EFFORT","SOAK.TIME","N.hooks.Fixed"),Expl.varS)]
    
    
    Res.var.in.log=rep("NO",length(resp.vars))    #fit response var in log space or not?
    
    Store.mod.out.observer=vector('list',length(resp.vars))
    names(Store.mod.out.observer)=resp.vars
    
    
    #Prelim analysis
    # fn.bxplt=function(D,var,FactR)
    # {
    #   
    #   for(mn in 1:n.rv)
    #   {
    #     idres=match(resp.vars[mn],names(D))
    #     idterm=match(var,names(D))
    #     X=D[,idterm]
    #     if(!FactR=="YES") X=cut(X,10)
    #     boxplot(D[,idres]~X,ylab=resp.vars[mn],col="grey75",cex.lab=1.5,cex.axis=1.25)
    #   }
    # }
    # par(mfcol=c(3,2),mai=c(.4,.45,.1,.1),mgp=c(2,.8,0))
    # fn.bxplt(D=DATA.shots.diversity,var="year",FactR="YES")
    # fn.bxplt(D=DATA.shots.diversity,var="Mid.Lat.mean",FactR="NO")
    # fn.bxplt(D=DATA.shots.diversity,var="BOTDEPTH.mean",FactR="NO")
    #fn.bxplt(D=DATA.shots.diversity,var="SOAK.TIME",FactR="NO")
    #fn.bxplt(D=DATA.shots.diversity,var="N.hooks.Fixed",FactR="NO")
    
    
    #respon var dist
    fn.fig(paste(hndl.eco,"dist_error_str",sep=""),2000,2400)
    smart.par(n.plots=4*length(resp.vars),MAR=c(3.5,2,1,1),OMA=c(1,1.5,.1,.1),MGP=c(2,.7,0))
    for(mn in 1:n.rv)
    {
      id=match(resp.vars[mn],names(DATA.shots.diversity))
      hist(DATA.shots.diversity[,id],main="",ylab="",xlab=resp.vars[mn])
      if(mn==1) mtext("normal scale",3)
      hist(log(DATA.shots.diversity[,id]),main='',ylab="",xlab=resp.vars[mn])
      if(mn==1) mtext("log",3)
      hist((DATA.shots.diversity[,id])^0.25,main='',ylab="",xlab=resp.vars[mn])
      if(mn==1) mtext("2root",3)
      hist(scale(DATA.shots.diversity[,id]),main='',ylab="",xlab=resp.vars[mn])
      if(mn==1) mtext("scaled",3)
    }
    dev.off()
    
    # #Standardise data to same number of shots per year   
    dummy=DATA.shots.diversity
    
    Standard.to.same.shots="NO"
    if(Standard.to.same.shots=="YES")
    {
      Min.shts=min(table(DATA.shots.diversity$Leg,DATA.shots.diversity$year))
      Min.shts=min(table(DATA.shots.diversity$year))
      Unic.yr=sort(unique(DATA.shots.diversity$year))
      dummy=vector('list',length(Unic.yr))
      for(i in 1:length(Unic.yr))
      {
        a=subset(DATA.shots.diversity,year==Unic.yr[i])
        id=sample(1:nrow(a),Min.shts)
        dummy[[i]]=a[id,]
      }
      dummy=do.call(rbind,dummy)
      
      #Show no year effect
      fn.violin=function(D)
      {
        D$year=factor(D$year)
        for(mn in 1:n.rv)
        {
          idres=match(resp.vars[mn],names(D))
          Lis=vector('list',length(levels(D$year)))
          for(l in 1:length(Lis)) Lis[[l]]=subset(D,year==levels(D$year)[l])[,idres]
          Lis=do.call(cbind,Lis)
          violin_plot(Lis,x_axis_labels=levels(D$year),main=resp.vars[mn],col="grey70")
        }
      }
      
      fn.fig(paste(hndl.eco,"paper/Violin",sep=""),2000,2400)
      par(mfcol=c(3,2),mai=c(.4,.45,.1,.1),mgp=c(2,.8,0))
      fn.violin(D=dummy)
      dev.off()
      
      fn.fig(paste(hndl.eco,"paper/Boxplot",sep=""),2000,2400)
      par(mfcol=c(3,2),mai=c(.4,.5,.1,.1),oma=c(1,1,.1,.1),mgp=c(2.5,.6,0),las=1)
      fn.bxplt(D=DATA.shots.diversity,var="year",FactR="YES")
      #fn.bxplt(D=dummy,var="year",FactR="YES")
      mtext("Year",1,-1,outer=T,cex=1.5)
      dev.off()
    }
    
    #GLM
    do.glm.ecos="YES"
    if(do.glm.ecos=="YES")
    {
      #glm approach on data subset
      for(i in 1:n.rv)
      {
        Store.mod.out.observer[[i]]=Mod.fn.glm(d=dummy,
                                               ResVar=resp.vars[i],
                                               Expl.vars=Expl.varS,
                                               Predictrs=Predictors,
                                               FactoRs=FactoRS,
                                               OFFSET=OFFSETT,
                                               log.var=Res.var.in.log[i])
      }
      
      #glm approach on data set as is
      # for(i in 1:n.rv)
      # {
      #   Store.mod.out.observer[[i]]=Mod.fn.glm(d=DATA.shots.diversity,
      #      ResVar=resp.vars[i],Expl.vars=Expl.varS,Predictrs=Predictors,
      #      FactoRs=FactoRS,OFFSET="offset(log.EFFORT)",log.var=Res.var.in.log[i])
      # }
      
      #Anova tables
      TABL=vector('list',length(resp.vars))
      names(TABL)=resp.vars
      
      for(i in 1:n.rv)
      {
        Modl=Store.mod.out.observer[[i]]$model
        
        #each term
        Anova.tab=anova(Modl, test = "Chisq")
        n=2:length(Anova.tab$Deviance)
        Term.dev.exp=100*(Anova.tab$Deviance[n]/Modl$null.deviance)
        names(Term.dev.exp)=rownames(Anova.tab)[n]
        
        #nice table
        Anov.tab=as.data.frame.matrix(Anova.tab)
        Term.tab=data.frame(Percent.dev.exp=Term.dev.exp)
        Anova.tab=Anova.tab[-1,match(c("Deviance","Pr(>Chi)"),names(Anova.tab))]
        Anova.tab=cbind(Anova.tab,Term.tab)
        Anova.tab=Anova.tab[,-match("Deviance",names(Anova.tab))]
        Anova.tab$"Pr(>Chi)"=ifelse(Anova.tab$"Pr(>Chi)"<0.001,"<0.001",round(Anova.tab$"Pr(>Chi)",3))
        Total=Anova.tab[1,]
        Total$"Pr(>Chi)"=""
        Total$Percent.dev.exp=sum(Anova.tab$Percent.dev.exp)
        rownames(Total)="Total"
        Anova.tab=rbind(Anova.tab,Total)
        Anova.tab$Percent.dev.exp=round(Anova.tab$Percent.dev.exp,2)
        TABL[[i]]=Anova.tab
      }
      
      TABL=do.call(cbind,TABL)
      TABL.ecosys=cbind(TERM=row.names(TABL),TABL)
      
    }
    
  }
  
}


# REPORT_CATCH RATES FROM FISHERY INDEPENDENT SURVEYS----------------------------------------------------------------------
if(Do.abundance=="YES")
{
  setwd("C:/Matias/Analyses/Surveys/Naturaliste_longline/outputs/Abundance")
  
  
  #Table 1
  write.csv(Table.1,"Paper/Table.1.csv",row.names=F)   
  write.csv(Table.1.scalies,"Paper/Table.1.scalies.csv",row.names=F)
  
  
  #Word cloud of shark species  
  library(tm)    #text mining
  library(SnowballC)
  library(wordcloud)
  #Longline
  fn.fig("Paper/Species cloud_longline",2400,2400)
  a=subset(Table.1,!is.na('LL.Numbers_2001-2017'))
  wordcloud(a$COMMON_NAME,a$'LL.Numbers_2001-2017',min.freq=1,
            scale=c(6,1),
            colors="black")
  #colors=rev(gray.colors(30, start = 0.1, end = 0.6, gamma = 2.2)))
  #colors=brewer.pal(12, "Paired"))
  dev.off()

  #dropline
  fn.fig("Paper/Species cloud_dropline",2400,2400)
  a=Table.1
  names(a)[match('DL.Numbers_1994-2000',names(a))]='dummy'
  a$dummy=as.numeric(a$dummy)
  a=subset(a,!is.na(dummy))
  wordcloud(a$COMMON_NAME,a$dummy,min.freq=1,
            scale=c(2,1),
            colors="black")
  dev.off()
  
  
  fn.fig("Paper/Species cloud",2400,2400)
  par(mfrow=c(2,1),mar=c(.1,.1,.1,.1),oma=c(.1,.1,.1,.1))
  #dropline
  a=Table.1
  names(a)[match('DL.Numbers_1994-2000',names(a))]='dummy'
  a$dummy=as.numeric(a$dummy)
  a=subset(a,!is.na(dummy))
  wordcloud(a$COMMON_NAME,a$dummy,min.freq=1,
            scale=c(3,.7),
            colors="black")
  legend('topleft',"Dropline",bty='n')
  #longline
  a=subset(Table.1,!is.na('LL.Numbers_2001-2017'))
  wordcloud(a$COMMON_NAME,a$'LL.Numbers_2001-2017',min.freq=1,
            scale=c(3,.7),
            colors="black")
  legend('topleft',"Longline",bty='n')
  #colors=rev(gray.colors(30, start = 0.1, end = 0.6, gamma = 2.2)))
  #colors=brewer.pal(12, "Paired"))
  dev.off()
  

  #stacked area plot of top shark species by time
  n=10
  fn.stack.area=function(tops,DAT,NME)
  {
    topS=names(tops[1:n])
    not.topS=names(tops[n+1:length(tops)])
    a=subset(DAT,TYPE=="Elasmo" & !SPECIES%in%not.topS)
    a$SPECIES=factor(a$SPECIES,levels=topS)
    d.stacked.plot=aggregate(Number~SPECIES+year,a,sum)
    d.stacked.plot=reshape(d.stacked.plot,v.names = "Number", idvar = "year",
                           timevar = "SPECIES", direction = "wide")
    d.stacked.plot[is.na(d.stacked.plot)]=0
    CL=(rainbow(nrow(d.stacked.plot),start=0.45,end=0.9,alpha=0.95))
    
    fn.fig(NME,2400,2400)
    par(mai=c(.8,.8,.5,.5),xpd=TRUE,mgp=c(2.5,.8,0),las=1)
    stackpoly(d.stacked.plot[,2:ncol(d.stacked.plot)],stack=T,xaxlab=d.stacked.plot$year,
              col=CL,border=1,ylab="Number of sharks",xlab="Year",cex.lab=2.25,cex.axis=1.35)
    legend("top",substr(names(d.stacked.plot)[2:ncol(d.stacked.plot)],8,9),bty='n',fill=CL,
           horiz=T,inset=c(0,-0.05),cex=1.125)
    dev.off()
    
  }
  
  #Longlines
  fn.stack.area(tops=rev(sort(table(subset(DATA,TYPE=="Elasmo")$SPECIES))),
                DAT=DATA,NME="Report/Stack.up.plot.report")
  
  #Droplines
  fn.stack.area(tops=rev(sort(table(subset(DATA.DL,TYPE=="Elasmo")$SPECIES))),
                DAT=DATA.DL,NME="Report/Stack.up.plot.DL.report")
  
  
 
  #Map of shots and stations
  fn.map.sht=function(d,NME)
  {
    fn.fig(NME,2400,2400)
    par(mgp=c(1,.8,0))
    plotMap(worldLLhigh, xlim= c(112,122),ylim=c(-26,-16.8),plt = c(.1, 1, 0.075, 1),
            col="firebrick",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
    box()
    text(118,-22,("Western"),col="white", cex=3)
    text(118,-23,("Australia"),col="white", cex=3)
    points(d$Mid.Long,d$Mid.Lat,cex=2.5,col=rgb(.5,.5,.5,alph=0.075),pch=19)
    with(Fixed.Stations[1:20,],points(Fix.St.mid.lon,Fix.St.mid.lat,pch=21,bg="white",cex=1.5))
    with(Fixed.Stations[1:20,],text(Fix.St.mid.lon*0.999,Fix.St.mid.lat,Station.no.,pos=2,cex=1.5,col="dodgerblue4"))
    mtext("Longitude (?E)",1,-1.5,outer=T,cex=2)
    mtext("Latitude (?S)",2,line=-1.85,outer=T,cex=2)
    axis(1,seq(112,120,2),seq(112,120,2),cex.axis=1.25)
    axis(2,seq(-16,-26,-2),seq(16,26,2),las=1,cex.axis=1.25)
    dev.off()
    
  }
  #Longlines
  fn.map.sht(d=D.map[!duplicated(D.map$SHEET_NO),],
             NME="Report/Map.stations")
  
  #Droplines
  fn.map.sht(d=DATA.DL[!duplicated(DATA.DL$SHEET_NO),],
             NME="Report/Map.stations.DL")
  
  #Scalefish size distribution
  Scale.size=a=subset(DATA,Taxa=="Teleost")
  Tab.scl=sort(table(Scale.size$SPECIES))
  Tab.scl=names(Tab.scl[Tab.scl>10])
  fun.his=function(d)
  {
    d=subset(d,TL>0 | !is.na(TL))
    NM=unique(d$COMMON_NAME)
    hist(d$TL,ylab="",xlab="",col="grey60",main=NM,cex.main=1.5,cex.axis=1.25)
    box()
  }
  fn.fig("Report/Scalefish.size.comp",2400,2400)
  par(mfcol=n2mfrow(length(Tab.scl)),mai=c(.275,.15,.3,.3),oma=c(2.5,3.5,.1,.1),las=1,mgp=c(1,.6,0))
  for(i in 1:length(Tab.scl))fun.his(d=subset(Scale.size,SPECIES==Tab.scl[i]))
  mtext("Total length(cm)",1,0,outer=T,cex=1.5)
  mtext("Frequency",2,1.5,outer=T,cex=1.5,las=3)
  dev.off()
  
  
  #Scatterplot of species by station  (only postive records plot)   
  plt.dis.sp=TARGETS
  for(i in 1:length(plt.dis.sp))
  {
    fn.fig(paste("density_main_sp/",plt.dis.sp[i],sep=""),2400,2400)
    par(mfcol=c(2,1),mar=c(2,2,1,1))
    DAT=subset(DATA.list[[i]],Catch.Target>0)
    DAT$cpue=250*DAT$Catch.Target/(DAT$SOAK.TIME*DAT$N.hooks.Fixed)
    smoothScatter(DAT[,match(c("Mid.Long","Mid.Lat"),names(DAT))],main=paste(names(DATA.list)[i],"catch"),
                  nrpoints = 0,ylab="Lat",xlab="Long",ylim=c(-25.98050,-16.67147),xlim=c(112.6228,122.2861))
    with(Fixed.Stations[1:20,],points(Fix.St.mid.lon,Fix.St.mid.lat,pch=19,col=2,cex=1.5))
    with(Fixed.Stations[1:20,],text(Fix.St.mid.lon,Fix.St.mid.lat,Station.no.,pos=3,cex=0.75,col=1))
    
    points(113.661,-24.884,cex=1.5,pch=17,col=3)
    text(113.661,-24.884,"Carnarvon",pos=4)
    points(115.5417,-20.45,cex=1.5,pch=17,col=3)
    text(115.5417,-20.8,"Montebellos",pos=3, srt=330)
    points(122.2359,-17.9614,cex=1.5,pch=17,col=3)
    text(122.3,-18.2,"Broome",pos=2, srt=90)
    
    
    Agg=aggregate(cpue~round(Mid.Long,1)+round((Mid.Lat),1),DAT,mean)
    names(Agg)=c("Long","Lat","cpue")
    wide <- reshape(Agg, v.names = "cpue", idvar = "Long",
                    timevar = "Lat", direction = "wide")
    wide=wide[order(wide[,1]),]
    
    numInt=20
    Breaks=quantile(Agg$cpue,probs=seq(0,1,1/numInt),na.rm=T)
    couleurs=rev(heat.colors(numInt))
    
    image(wide$Long,sort(unique(round((DAT$Mid.Lat),1))),as.matrix(wide[,-1]),
          col =couleurs,breaks=Breaks,ylab="Lat",xlab="Long",
          xlim=c(112,124),ylim=c(-26,-15))
    color.legend(123,-15,124,-25,
                 paste(round(Breaks[seq(1,numInt,2)],1)," sharks/250 hour hook"),rect.col=couleurs,gradient="y",cex=.75)
    
    
     dev.off()
  }
  
  #check classification of stations by plotting the YES and NO separately against the stations
  for (i in 1:N.yrs)
  {
    dum=subset(DATA,year==YEAR[i] & N.hooks.Fixed>20 & BOTDEPTH <MaxDepth)
    dum=dum[!duplicated(dum$SHEET_NO),]
    YLIM=c(min(dum$Mid.Lat)-0.5,max(dum$Mid.Lat)+0.5)
    XLIM=c(min(dum$Mid.Long)-0.5,max(dum$Mid.Long)+0.5)   
    dum.FixeS=subset(dum,FixedStation=="YES")
    dum.No.FixeS=subset(dum,FixedStation=="NO")  
    fn.fig(paste("Stations_by_year/",YEAR[i],".Station",sep=""),2400,2400)
    par(las=1)
    plot(Fixed.Stations$Long.1,Fixed.Stations$Lat.1,pch=19,cex=1,xlim=XLIM,ylim=YLIM,ylab="LAT",xlab="LONG")
    points(Fixed.Stations$Long.2,Fixed.Stations$Lat.2,pch=19,cex=1)
    points(Fixed.Stations$Fix.St.mid.lon,Fixed.Stations$Fix.St.mid.lat,pch=19,cex=1,col=2)
    points(dum.FixeS$Mid.Long,dum.FixeS$Mid.Lat,pch=19,cex=.5,col=4)
    points(dum.No.FixeS$Mid.Long,dum.No.FixeS$Mid.Lat,pch=19,cex=.5,col=3)
    legend("topleft",paste(YEAR[i]),bty='n',cex=2)
    legend("bottomright",c("Start station","End station","mid point Station","Classed as FS","Class as No FS"),
           bty='n',cex=1.5,pch=19,col=c(1,1,2,4,3))
    dev.off()
  }
  
  
  #Figure S2.  
  BIN=10
  fn.fig("Paper/Figure S2",2000,2400)
  par(mfcol=n2mfrow(N.species),mai=c(.3,.395,.1,.1),oma=c(3,1.5,2,.1),las=1,mgp=c(.04,.6,0))
  for(i in 1:N.species)
  {
    BINs=BIN
    if(Tar.names[i]=="Tiger shark") BINs=20
    SizeFreq.fn(TARGETS[i],Tar.names[i],
                Fem.Size.Mat[[match(names(TARGETS)[i],names(Fem.Size.Mat))]],BINs)
    mtext(Tar.names[i],3,line=.05,cex=1.1)
    if(bySEX=="YES" & i==1)
    {
      #plot(1:10,1:10,bty='n',pch='',xaxt='n',yaxt='n',ann=F)
      legend("topleft",c("female","male"),fill=c("white","grey70"),bty='n',cex=1.4)
    }
    
  }
  mtext("Frequency",side=2,line=-0.35,outer=T,cex=1.5,las=3)
  mtext("Fork length (cm)",side=1,line=1,outer=T,cex=1.5)
  dev.off()
  
  
  #Boxplot FL
  fn.fig("Paper/Boxplot_FL",2000,2400)
  par(mfcol=n2mfrow(N.species),mai=c(.3,.395,.1,.1),oma=c(3,1.5,2,.1),las=1,mgp=c(.04,.6,0))
  for(i in 1:N.species)
  {
    Kruskal[[i]]=SizeTemp.fn(TARGETS[i],Tar.names[i])
    mtext(Tar.names[i],3,line=.05,cex=1.1)
  }
  mtext("Fork length (cm)",side=2,line=-0.2,outer=T,cex=1.5,las=3)
  mtext("Year",side=1,line=1,outer=T,cex=1.5)
  dev.off()
  
  
  
  #Depth distributions
  Store.dens=vector('list',N.species)
  for ( i in 1:N.species) Store.dens[[i]]=Z.dist(TARGETS[i])
  COLS=rainbow(N.species)
  LTY=rep(c(1:5,1,3),2)
  LWD=rep(c(2,2,2,3,2,2,3),2)
  fn.fig("Depth.dist",2400,2400)
  par(mfcol=c(1,1),mai=c(.6,.6,.1,.1),oma=c(1,.1,.1,.1),las=1)
  plot(Store.dens[[4]],col="transparent",xlim=c(0,350),ylim=c(0,.021),yaxt='n',
       main="",xlab="",ylab="",cex.axis=1.5)
  for ( i in 1:N.species) lines(Store.dens[[i]],lwd=LWD[i],col=COLS[i],lty=LTY[i])
  legend("topright",Tar.names,lty=LTY,lwd=LWD,col=COLS,bty='n',cex=1.5)
  mtext("Depth (m)",side=1,line=3,font=1,las=0,cex=2.,outer=F)
  mtext("Density",side=2,line=0.5,font=1,las=0,cex=2.5,outer=F)
  dev.off()
  
  # Plot location of positive catches
  #by species and year
  #for ( i in 1:N.species)pos.catch.pos(DATA.list[[i]],REFS[[i]],SCALES[i])
  
  #species combined
  #par(mfcol=c(5,3),mai=c(.3,.395,.1,.1),oma=c(3,1.5,2,.1),las=1,mgp=c(.04,.6,0))
  #for ( i in 1:N.species)pos.catch.pos(DATA.list[[i]],REFS[[i]],SCALES[i]) 
  
  
  #Figure 1.
  do.map="NO"
  add.depth="NO"
  if(do.map=="YES")  
  {
    #North West WA
    library(rgdal)
    source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Plot.Map.R")
    
    JA_Northern_Shark=readOGR("C:/Matias/Data/Mapping/Shark_shape_files/JA_Northern_Shark.shp", layer="JA_Northern_Shark") 
    WA_Northern_Shark=readOGR("C:/Matias/Data/Mapping/Shark_shape_files/NorthCoastShark_s43.shp", layer="NorthCoastShark_s43") 
    WA_Northern_Shark_2=readOGR("C:/Matias/Data/Mapping/Shark_shape_files/NorthWestCoastShark_s43.shp", layer="NorthWestCoastShark_s43") 
    #Shark.zones=readOGR("C:/Matias/Data/Mapping/Shark_shape_files/FisheriesGuide_ConsolidatedNoticesandOrdersDPIRD_061.shp", layer="FisheriesGuide_ConsolidatedNoticesandOrdersDPIRD_061") 
    #Exception.zones <- readOGR("C:/Matias/Data/Mapping/Shark_shape_files/FisheriesGuide_InstrumentofExemptionDPIRD_051.shp", layer="FisheriesGuide_InstrumentofExemptionDPIRD_051")
    
    
    if(add.depth=="YES") if(!exists("reshaped"))reshaped=as.matrix(reshape(Bathymetry,idvar="V1",timevar="V2",v.names="V3", direction="wide"))
    
    LEG=c(as.character(5),"",as.character(10),"",as.character(20))
    
    fn.scale=function(x,max,scaler) ((x/max)^0.5)*scaler
    
    South.WA.long=c(109,129)
    South.WA.lat=c(-26,-12)
    Xlim=South.WA.long
    Ylim=South.WA.lat
    Perhooks=100  #present cpue for 100 hooks rather than 1 hook
    SCALER=3
    sht.col=rgb(.1, .1, .1, alpha=.15)
    
    #Col.Ning=rgb(.2,.2,.2,alpha=.4)
    Col.Ning="grey80"
    #Col.NSF=rgb(.4,.4,.4,alpha=.1)
    Col.NSF="grey50"
    Col.land="grey95"
    
    fn.fig("Paper/Figure 1",2400,2400)
    par(mfcol=n2mfrow(N.species+1),mar=c(1,1,.5,.5),oma=c(3,3,1,.3),las=1,mgp=c(.04,.6,0))
    
    #Add Australia and Closures
    plot(1,xlim=c(Xlim[1]*0.9995,Xlim[2]*0.99),ylim=Ylim,xlab="",ylab="",axes=F,main="")
    #plot(1,xlim=South.WA.long,ylim=South.WA.lat,xlab="",ylab="",axes=F,main="")
    
    #NSF
    plot(WA_Northern_Shark,add=T,col=Col.Ning)
    
    #WANCS open
    #mpts <- raster::geom(subset(Exception.zones,name=="The WA North Coast Shark Fishery"))
    #mpts=as.data.frame(mpts)
    #polygon((mpts$x/1e5)-14,(mpts$y/1e5)+1.5,col="transparent")
    polygon(c(123.75,123.75,122.0204,121.4045,120,
              120,122.3488,123.211, 122.9647),
            c(-16.33,-11.50386,-11.65960,-12.39937,-12.63298,
              -17.96708,-17.96708,-17.69453,-16.33),col=Col.NSF)
    
    #WANCS closure
    plot(JA_Northern_Shark,add=T,col=Col.NSF)
    text(116.8,-16.75,"Closed",srt=45,cex=1.1)
    text(117.5,-18,"since 2005",srt=45,cex=1.1)
    #Ningaloo closure
    plot(WA_Northern_Shark_2,add=T,col=Col.Ning)
    text(111.1,-21.75,"Closed",srt=65,cex=1.1)
    text(112,-23,"since 1993",srt=65,cex=1.1)
    
    polygon(WAcoast$Longitude,WAcoast$Latitude, col=Col.land)
    
    #axis(side = 1, at =South.WA.long[1]:South.WA.long[2], labels = F, tck = -.015)
    #axis(side = 2, at = South.WA.lat[2]:South.WA.lat[1], labels = F, tck = -.015)
    axis(side = 1, seq(South.WA.long[1],South.WA.long[2],2), labels =F, tck = -.035)
    axis(side = 2, seq(South.WA.lat[1],South.WA.lat[2],2), labels = F, tck = -.035)
    
    #Stations
    ddd=Fixed.Stations
    names(ddd)[match('Station.no.',names(ddd))]='STNum'
    with(subset(ddd,!STNum=='additional'),points(Long.1,Lat.1,pch=21,bg="black",col=1))
    axis(2,seq(round(Ylim[1]),round(Ylim[2]),2),-seq(round(Ylim[1]),round(Ylim[2]),2),cex.axis=1.25)
    box()
    
    #Each species   
    for ( i in 1:N.species)
    {
      pos.cpue(DATA.list[[i]],SCALER,ADD.zero="NO")    
      mtext(paste(TARGETS.name[i]),side=3,line=-1.5,font=1,las=0,cex=1)
      if(i%in%c(2,5,8))axis(1,seq(round(Xlim[1]),round(Xlim[2]),2),seq(round(Xlim[1]),round(Xlim[2]),2),cex.axis=1.25)
      if(i%in%1:2)axis(2,seq(round(Ylim[1]),round(Ylim[2]),2),-seq(round(Ylim[1]),round(Ylim[2]),2),cex.axis=1.25)
    }
    
    # #Australia
    #par(fig = c(.825, .975, 0.025, .125), mar=c(0,0,0,0), new=TRUE)
    par(fig = c(.175, .3, .675, .8), mar=c(0,0,0,0), new=TRUE)
    OZ.lat=c(-44.5,South.WA.lat[2]);OZ.long=c(South.WA.long[1],155)
    plotMap(worldLLhigh, xlim=OZ.long,ylim=OZ.lat,plt = c(.1, 1, 0.075, 1),
            col='black',tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
    #col=Col.land,tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
    box()
    polygon(x=c(rep(South.WA.long[2],2),rep(South.WA.long[1],2)),y=c(South.WA.lat,rev(South.WA.lat)),lwd=1.5,col=rgb(.4,.2,.2,alpha=.4))
    text(135,-25,("Australia"),col="white", cex=1.35)
    mtext(expression(paste("Longitude (",degree,"E)",sep="")),side=1,line=1.5,font=1,las=0,cex=1.35,outer=T)
    mtext(expression(paste("Latitude (",degree,"S)",sep="")),side=2,line=1,font=1,las=0,cex=1.35,outer=T)
    dev.off() 
  }
  
  do.GitHub.map="NO"
  if(do.GitHub.map=="YES")  
  {
    #North West WA
    library(rgdal)
    source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Plot.Map.R")
    
    JA_Northern_Shark=readOGR("C:/Matias/Data/Mapping/Shark_shape_files/JA_Northern_Shark.shp", layer="JA_Northern_Shark") 
    WA_Northern_Shark=readOGR("C:/Matias/Data/Mapping/Shark_shape_files/NorthCoastShark_s43.shp", layer="NorthCoastShark_s43") 
    WA_Northern_Shark_2=readOGR("C:/Matias/Data/Mapping/Shark_shape_files/NorthWestCoastShark_s43.shp", layer="NorthWestCoastShark_s43") 
    
    
    SCALER=3
    LEG=c(as.character(5),"",as.character(10),"",as.character(20))
    
    fn.scale=function(x,max,scaler) ((x/max)^0.5)*scaler
    
    South.WA.long=c(109,129)
    South.WA.lat=c(-26,-12)
    Xlim=South.WA.long
    Ylim=South.WA.lat
    Col.Ning="grey80"
    Col.NSF="grey50"
    Col.land="grey95"
    
    fn.fig("Paper/Figure 1_GitHub",2400,2400)
    par(mfcol=c(1,1),mar=c(1,1,.5,.5),oma=c(3,3,1,.1),las=1,mgp=c(.04,.6,0))
    
    #Add Australia and Closures
    plot(1,xlim=c(Xlim[1]*0.9995,Xlim[2]*.9975),ylim=Ylim,xlab="",ylab="",axes=F,main="")
    
    #NSF
    plot(WA_Northern_Shark,add=T,col=Col.Ning)
    
    #WANCS open
    polygon(c(123.75,123.75,122.0204,121.4045,120,
              120,122.3488,123.211, 122.9647),
            c(-16.33,-11.50386,-11.65960,-12.39937,-12.63298,
              -17.96708,-17.96708,-17.69453,-16.33),col=Col.NSF)
    
    #WANCS closure
    plot(JA_Northern_Shark,add=T,col=Col.NSF)
    text(116.8,-16.75,"Closed",srt=45,cex=1.1)
    text(117.5,-18,"since 2005",srt=45,cex=1.1)
    #Ningaloo closure
    plot(WA_Northern_Shark_2,add=T,col=Col.Ning)
    text(111.1,-21.75,"Closed",srt=65,cex=1.1)
    text(112,-23,"since 1993",srt=65,cex=1.1)
    
    polygon(WAcoast$Longitude,WAcoast$Latitude, col=Col.land)
    box()
    
    #Stations
    ddd=Fixed.Stations
    names(ddd)[match('Station.no.',names(ddd))]='STNum'
    with(subset(ddd,!STNum=='additional'),points(Long.1,Lat.1,pch=21,bg="black",col=1))
    axis(2,seq(round(Ylim[1]),round(Ylim[2]),2),-seq(round(Ylim[1]),round(Ylim[2]),2),cex.axis=1.25)
    axis(side = 1, seq(South.WA.long[1],South.WA.long[2],2), 
         labels =seq(South.WA.long[1],South.WA.long[2],2), tck = -.015,cex.axis=1.25)
    
    
    # #Australia
    par(fig = c(.6, .9, .1, .4), mar=c(0,0,0,0), new=TRUE)
    OZ.lat=c(-44.5,South.WA.lat[2]);OZ.long=c(South.WA.long[1],155)
    plotMap(worldLLhigh, xlim=OZ.long,ylim=OZ.lat,plt = c(.1, 1, 0.075, 1),
            col='black',tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
    box()
    polygon(x=c(rep(South.WA.long[2],2),rep(South.WA.long[1],2)),y=c(South.WA.lat,rev(South.WA.lat)),lwd=1.5,col=rgb(.4,.2,.2,alpha=.4))
    text(135,-25,("Australia"),col="white", cex=1.35)
    mtext(expression(paste("Longitude (",degree,"E)",sep="")),side=1,line=1.2,font=1,las=0,cex=1.35,outer=T)
    mtext(expression(paste("Latitude (",degree,"S)",sep="")),side=2,line=1,font=1,las=0,cex=1.35,outer=T)
    dev.off() 
  }

  #Figure S1.   
  Ymax.vec=rep(.3,N.species)
  fn.fig("Paper/Figure S1",2000,2400)
  par(mfcol=n2mfrow(N.species),mai=c(.3,.3,.01,.1),oma=c(3,3,2,.1),las=1,mgp=c(.04,.6,0))
  for ( i in 1:N.species) plot0(DATA.list[[i]],Ymax.vec[i],Tar.names[i])
  mtext("Number of animals per set",side=1,line=1,font=1,las=0,cex=1.25,outer=T)
  mtext("Proportion of sets with catch",side=2,line=.5,font=1,las=0,cex=1.25,outer=T)
  dev.off()
  
  

  #Effort all vs fixed stations
  fn.fig("Paper/Effort all vs fixed stations",2400,2400)
  par(mfcol=c(1,1),mai=c(.6,.8,.1,.1),oma=c(1,.1,.1,.1),las=1,mgp=c(2.15, .65, 0))
  plot(2002:YEAR[length(YEAR)],2002:YEAR[length(YEAR)],ylim=c(0,26),ylab="Effort (1000 hook hours)",xlab="Year",
       cex.lab=1.95, col="transparent",cex.axis=1.5)
  All.Eff=Eff.All.stations[,2]/1000
  lines(YEAR[1:2],All.Eff[1:2],lwd=2.5,col=1)
  lines(YEAR[3:length(YEAR)],All.Eff[3:length(YEAR)],lwd=2.5,col=1)
  
  Fixed.Eff=Eff.Fixed.stations[,2]/1000
  lines(YEAR[1:2],Fixed.Eff[1:2],lwd=2.5,col="grey40",lty=2)
  lines(YEAR[3:length(YEAR)],Fixed.Eff[3:length(YEAR)],lwd=2.5,col="grey40",lty=2)
  legend("topleft",c("All stations", "Fixed stations only"),col=c("black","grey40"),lty=1:2,bty='n',lwd=2,cex=1.75)
  axis(1,2002:YEAR[length(YEAR)],F)
  dev.off()
  
  # Anova tables for abundance and size trends
  ANOV.TAB=vector('list',N.species) 
  for (i in 1:N.species)  ANOV.TAB[[i]]=Anova.tab(mod=Store[[i]]$Fit,Fcol="Chi.sq")
  ANVA.abundance=do.call(rbind,ANOV.TAB)
  ANVA.abundance$counter=1:nrow(ANVA.abundance)
  
  ANOV.TAB.FL=vector('list',N.species) 
  for (i in 1:N.species)  ANOV.TAB.FL[[i]]=Anova.tab(mod=Store.size[[i]]$Fit,Fcol="F")
  ANVA.size=do.call(rbind,ANOV.TAB.FL)
  ANVA.size$counter=1:nrow(ANVA.size)
  
  ANoVA=full_join(ANVA.abundance,ANVA.size,by="Join",all=T)%>%
        arrange(counter.y)
  write.csv(ANoVA,paste(getwd(),"/Paper/Anovas/Table1_anova.csv",sep=""),row.names=F)
  
  
  #---Abundance trends---   
  ADD.P="NO"
  #2. Effect of latitude on abundance   
      #fixed stations
  fn.fig("Paper/Figure 3",2400,2400)
  par(mfcol=n2mfrow(N.species-sum(sapply(PRED.lat,is.null))),mai=c(.3,.38,.15,.1),oma=c(2,1.25,1,.1),las=1,mgp=c(.04,.6,0))
  for(i in 1:N.species)
  {
    dummy=PRED.lat[[i]]
    if(!is.null(dummy))
    {
      #names(dummy)[6:7]=c("lower.CL","upper.CL")
      dummy=cbind(dummy,CV=NA)
      dummy$yr=dummy$Mid.Lat
      id=match("response",names(dummy))
      if(length(id)>0) names(dummy)[id]="emmean"
      dummy$MEAN=dummy$emmean
      dummy$UP=dummy$upper.CL
      dummy$LOW=dummy$lower.CL
      a=fun.plot.yr.pred(dummy,X="Mid.Lat",normalised="YES",REV="NO",n.seq=1,YLIM=NULL,XLIM=NULL,Type="polygon")
      mtext(Tar.names[i],3,line=.05,cex=1.1)
    }
    
    #add p values
    if(ADD.P=="YES")
    {
          Count.p=Zero.p=NULL
    Nms=Sig.term.coeff[[i]]$Sigi.count   
    id=which(names(Nms)=='log.Mid.Lat')
    if(length(id)>0)
    {
      Count.p=Nms[id]
      dd1=Count.p
      if(dd1>=0.01 & dd1<0.05) Count.p=paste("Count","*",sep="")
      if(dd1>=0.001 & dd1<0.01) Count.p=paste("Count","**",sep="")
      if(dd1<0.001) Count.p=paste("Count","***",sep="")
    }
    Nms=Sig.term.coeff[[i]]$Sigi.zero
    id.z=which(names(Nms)=='log.Mid.Lat')
    if(length(id.z)>0) 
    {
      Zero.p=Nms[id.z]
      dd1=Zero.p
      if(dd1>=0.01 & dd1<0.05) Zero.p=paste("Zero","*",sep="")
      if(dd1>=0.001 & dd1<0.01) Zero.p=paste("Zero","**",sep="")
      if(dd1<0.001) Zero.p=paste("Zero","***",sep="")
    }
    LGN=c(Count.p,Zero.p)
    if(!is.null(LGN))legend("bottomleft",LGN,bty='n',cex=1.25)
    }
  }
  mtext("Relative CPUE",2,outer=T,line=-0.5,cex=1.5,las=3)
  mtext(expression(paste("Latitude (",degree,"S)",sep="")),1,outer=T,line=0.5,cex=1.5)
  dev.off()
  

  #3. Effect of depth on abundance
      #fixed stations 
  fn.fig("Paper/Figure 4",2400,2400)
  par(mfcol=n2mfrow(N.species-sum(sapply(PRED.z,is.null))),mai=c(.3,.38,.15,.1),oma=c(2,1.25,1,.1),las=1,mgp=c(.04,.6,0))
  limX=c(20,200)
  for(i in 1:N.species)
  {
    dummy=PRED.z[[i]]
    if(!is.null(dummy))
    {
      #names(dummy)[6:7]=c("lower.CL","upper.CL")
      dummy=cbind(dummy,CV=NA)
      dummy$yr=dummy$BOTDEPTH
      id=match("response",names(dummy))
      if(length(id)>0) names(dummy)[id]="emmean"
      dummy$MEAN=dummy$emmean
      dummy$UP=dummy$upper.CL
      dummy$LOW=dummy$lower.CL
      a=fun.plot.yr.pred(dummy,X="BOTDEPTH",normalised="YES",REV="NO",n.seq=10,YLIM=NULL,XLIM=limX,Type="polygon")
      mtext(Tar.names[i],3,line=.05,cex=1.1)
    }
    if(ADD.P=="YES")
    {
      #add p values
      Count.p=Zero.p=NULL
      Nms=Sig.term.coeff[[i]]$Sigi.count
      id=which(names(Nms)=='log.BOTDEPTH')
      if(length(id)>0)
      {
        Count.p=Nms[id]
        dd1=Count.p
        if(dd1>=0.01 & dd1<0.05) Count.p=paste("Count","*",sep="")
        if(dd1>=0.001 & dd1<0.01) Count.p=paste("Count","**",sep="")
        if(dd1<0.001) Count.p=paste("Count","***",sep="")
      }
      Nms=Sig.term.coeff[[i]]$Sigi.zero
      id.z=which(names(Nms)=='log.BOTDEPTH')
      if(length(id.z)>0) 
      {
        Zero.p=Nms[id.z]
        dd1=Zero.p
        if(dd1>=0.01 & dd1<0.05) Zero.p=paste("Zero","*",sep="")
        if(dd1>=0.001 & dd1<0.01) Zero.p=paste("Zero","**",sep="")
        if(dd1<0.001) Zero.p=paste("Zero","***",sep="")
      }
      LGN=c(Count.p,Zero.p)
      Whre="topleft"
      if(Tar.names[i]=="Milk shark")Whre="bottomright"
      if(!is.null(LGN))legend(Whre,LGN,bty='n',cex=1.25)
    }
  }
  mtext("Relative CPUE",2,outer=T,line=-0.5,cex=1.5,las=3)
  mtext("Depth (m)",1,outer=T,line=0.5,cex=1.5)
  dev.off()
 
   
  #4. Effect of time on abundance
    #fixed stations
  INDEX=PRED.CPUE
  Plus=0.3
  fn.fig("Paper/Figure 2",2400,2400)
  par(mfcol=n2mfrow(N.species),mai=c(.3,.38,.15,.1),oma=c(2,1.25,1,.1),las=1,mgp=c(.04,.6,0))
  LimX=c(min(DATA.list$`Sandbar shark`$year),max(DATA.list$`Sandbar shark`$year))
  for(i in Species.cpue)
  {
    Nml=fn.nmnl(dat=subset(DATA.list[[i]],FixedStation=="YES" & BOTDEPTH<210),REL="YES")
    dummy=PRED.CPUE[[i]]
    #names(dummy)[5:6]=c("lower.CL","upper.CL")
    dummy$yr=dummy$year
    id=match("response",names(dummy))
    if(length(id)>0) names(dummy)[id]="emmean"
    dummy$MEAN=dummy$emmean
    dummy$UP=dummy$upper.CL
    dummy$LOW=dummy$lower.CL
    dummy$CV=100*dummy$SE/dummy$MEAN
    
    MaX=max(c(dummy$UP/mean(dummy$MEAN),Nml$up95),na.rm=T)
    
    INDEX[[i]]=fun.plot.yr.pred(dummy,X="year",normalised="YES",REV="NO",n.seq=5,YLIM=c(0,MaX),XLIM=LimX,Type="points")
    
    #add nominal
    #with(Nml,points(year+Plus, mean, "o", pch=16, lty=2, col="grey50"))
    with(Nml,points(year+Plus, mean, pch=16, cex=1.25, col="grey70"))
    suppressWarnings(with(Nml,arrows(x0=year+Plus, y0=low95,x1=year+Plus, y1=up95, 
           code=3, angle=90, length=0.025, col="grey70")))
    mtext(Tar.names[i],3,line=.05,cex=1.1)
  }
  mtext("Relative CPUE",2,outer=T,line=-0.5,cex=1.5,las=3)
  mtext("Year",1,outer=T,line=0.5,cex=1.5)
  plot.new()
  legend('topright',c("standardised","nominal"),bty='n',cex=1.5,pch=19,col=c("black","grey70"))
  dev.off()
  
  
  #5. Export Sandbar and Dusky sharks index   
      #Fixed stations
  hnd.indx="C:/Matias/Analyses/Data_outs/"
  for(i in 1:length(INDEX)) write.csv(INDEX[[i]],paste(hnd.indx,names(INDEX)[i],'/',names(INDEX)[i],".Srvy.FixSt.csv",sep=""),row.names=F)
  
  
  
  #6. create dusky and sandbar figures for RAR
      #fixed stations
  San.dusky.species=match(c("Sandbar shark","Dusky shark"),names(INDEX))
  fn.fig("Paper/Figure 2_sandbar_dusky",1200,2400)
  par(mfcol=c(2,1),mai=c(.35,.55,.25,.1),oma=c(1.5,1.25,.1,.1),las=1,mgp=c(.04,.6,0))
  for(i in San.dusky.species)
  {
    Nml=fn.nmnl(dat=subset(DATA.list[[i]],FixedStation=="YES"),REL="YES")
    
    dummy=PRED.CPUE[[i]]
    #names(dummy)[5:6]=c("lower.CL","upper.CL")
    dummy$yr=dummy$year
    id=match("response",names(dummy))
    if(length(id)>0) names(dummy)[id]="emmean"
    dummy$MEAN=dummy$emmean
    dummy$UP=dummy$upper.CL
    dummy$LOW=dummy$lower.CL
    dummy$CV=100*dummy$SE/dummy$MEAN
    
    MaX=max(c(dummy$UP/mean(dummy$MEAN),Nml$up95),na.rm=T)
    
    INDEX[[i]]=fun.plot.yr.pred(dummy,X="year",normalised="YES",REV="NO",n.seq=5,
                                YLIM=c(0,MaX),XLIM=LimX,Type="points")
    #add nominal
    with(Nml,points(year+Plus, mean, pch=16, cex=1.25, col="grey70"))
    suppressWarnings(with(Nml,arrows(x0=year+Plus, y0=low95,x1=year+Plus, y1=up95, 
                                     code=3, angle=90, length=0.025, col="grey70")))
    mtext(Tar.names[i],3,line=.05,cex=1.1)
  }
  mtext("Relative CPUE",2,outer=T,line=-0.5,cex=1.5,las=3)
  mtext("Year",1,outer=T,line=.25,cex=1.5)
  legend('topleft',c("standardised","nominal"),bty='n',cex=1.25,pch=19,col=c("black","grey70"))
  dev.off()
  

  
  #---Sex ratio trends---
  if(do.sex.ratio=="YES")
  {
    #fixed stations
    fn.fig("Paper/Prop_male_year",2000,2400)
    par(mfcol=c(5,3),mai=c(.3,.38,.01,.1),oma=c(3,1.25,2,.1),las=1,mgp=c(.04,.6,0))
    for(i in 1:N.species.sex)
    {
      fun.plot.yr.pred.arrows.size(BOOTS=STORE.BOOT.sex[[i]],normalised="NO")
      mtext(names(STORE.BOOT.sex)[i],3,line=.05,cex=1.1)
    }
    mtext("Proportion of males",2,outer=T,line=-0.5,cex=1.75,las=3)
    mtext("Year",1,outer=T,line=1,cex=1.75)
    dev.off()
    

    #Anova tables
    TABL.sex=vector('list',length(Store.SEX))
    names(TABL.sex)=names(Store.SEX)
    for (i in 1:N.species.sex) TABL.sex[[i]]=fn.anova.table(Store.SEX[[i]])
    TABL.sex=do.call(cbind,TABL.sex)
    TABL.sex=cbind(Terms=rownames(TABL.sex),TABL.sex)
    write.csv(TABL.sex,"Paper/Anovas/TABL.sex.csv",row.names=F)
    
    #Depth and latitude trends
    #select species for which depth or latitude are significant
    depth.sp=c("Blacktip sharks","Scalloped hammerhead",
               "Great hammerhead","Pigeye shark","Lemon shark")	
    Lat.sp=c("Sandbar shark","Milk shark","Blacktip sharks","Pigeye shark")
    fn.fig("Paper/Depth_Lat_effect_prop_male",2000,2400)
    par(mfcol=c(2,1),mai=c(.3,.3,.3,.1),oma=c(2,2,.3,.1),mgp=c(2.5, .5, 0),las=1)
    
    #Proportion of males changes with latitude
    Store.SEX.show.lat=Store.SEX[match(Lat.sp,names(Store.SEX))]
    fn.plt.pred.lat.z.sx.size(D=Store.SEX.show.lat,varX="LAT.seq",
                              varY="LAT.preds",REVRT="YES",CL=1:15,XLIM=c(17,25),YLIM=c(0,1))
    mtext("Latitude (?S)",1,line=1.9,cex=1.5)
    
    #Proportion of males changes with depth
    Store.SEX.show.depth=Store.SEX[match(depth.sp,names(Store.SEX))]
    
    fn.plt.pred.lat.z.sx.size(D=Store.SEX.show.depth,varX="Z.seq",
                              varY="Z.preds",REVRT="NO",CL=1:15,XLIM=c(0,250),YLIM=c(0,1))
    mtext("Depth (m)",1,line=1.9,cex=1.5)
    mtext("Probability of catching a male",2,0.75,outer=T,las=3,cex=1.65)
    dev.off()
    
  }
  
  
  
  #---Size trends---
  
  #2. Effect of latitude on size   
    #fixed stations
  fn.fig("Paper/Figure 6",2400,2400)
  par(mfcol=n2mfrow(N.species.size-sum(sapply(PRED.lat.size,is.null))),mai=c(.3,.38,.15,.1),oma=c(2,1.25,1,.1),las=1,mgp=c(.04,.6,0))
  for(i in 1:N.species.size)
  {
    dummy=PRED.lat.size[[i]]
    if(!is.null(dummy))
    {
      #names(dummy)[6:7]=c("lower.CL","upper.CL")
      dummy=cbind(dummy,CV=NA)
      dummy$yr=dummy$Mid.Lat
      dummy$MEAN=dummy$emmean
      dummy$UP=dummy$upper.CL
      dummy$LOW=dummy$lower.CL
      a=fun.plot.yr.pred(dummy,X="Mid.Lat",normalised="YES",REV="NO",n.seq=1,YLIM=NULL,XLIM=NULL,Type="polygon")
      LGn=names(PRED.lat.size)[i]
    }
    
    if(ADD.P=="YES")
    {
      dd1=as.data.frame(Store.size[[i]]$Signifcance)    
      dd1=dd1[match("log.Mid.Lat",row.names(dd1)),]$'Pr(>Chi)'
      P=NULL
      if(dd1>=0.01 & dd1<0.05) P="*"
      if(dd1>=0.001 & dd1<0.01) P="**"
      if(dd1<0.001) P="***"
      if(!is.null(P)) LGn=paste(LGn,P)
    }
    mtext(LGn,3,line=.05,cex=1.1)
  }
  mtext("Relative size",2,outer=T,line=-0.5,cex=1.5,las=3)
  mtext(expression(paste("Latitude (",degree,"S)",sep="")),1,outer=T,line=0.5,cex=1.5)
  dev.off()
  
  
  #3. Effect of depth on size
    #fixed stations
  limX=c(20,200)
  fn.fig("Paper/Figure 7",2400,2400)
  par(mfcol=n2mfrow(N.species.size-sum(sapply(PRED.z.size,is.null))),mai=c(.3,.38,.15,.1),oma=c(2,1.25,1,.1),las=1,mgp=c(.04,.6,0))
  for(i in 1:N.species.size)
  {
    dummy=PRED.z.size[[i]]
    if(!is.null(dummy))
    {
      #names(dummy)[6:7]=c("lower.CL","upper.CL")
      dummy=cbind(dummy,CV=NA)
      dummy$yr=dummy$Mid.Lat
      dummy$MEAN=dummy$emmean
      dummy$UP=dummy$upper.CL
      dummy$LOW=dummy$lower.CL
      a=fun.plot.yr.pred(dummy,X="BOTDEPTH",normalised="YES",REV="NO",n.seq=10,YLIM=NULL,XLIM=limX,Type="polygon")
      LGn=names(PRED.z.size)[i]
    }
    
    if(ADD.P=="YES")
    {
      dd1=as.data.frame(Store.size[[i]]$Signifcance)    #add p values
      dd1=dd1[match("log.BOTDEPTH",row.names(dd1)),]$'Pr(>Chi)'
      P=NULL
      if(dd1>=0.01 & dd1<0.05) P="*"
      if(dd1>=0.001 & dd1<0.01) P="**"
      if(dd1<0.001) P="***"
      if(!is.null(P)) LGn=paste(LGn,P)
    }
    mtext(LGn,3,line=.05,cex=1.1)
  }
  mtext("Relative size",2,outer=T,line=-0.5,cex=1.5,las=3)
  mtext("Depth (m)",1,outer=T,line=0.5,cex=1.5)
  dev.off()
  
  
  #4. Effect of time on size
    #fixed stations 
  fn.nmnl.size=function(dat,REL)
  {
    d = dat %>%
      mutate(year=as.numeric(as.character(year)))%>%
      group_by(year) %>%
      summarise(n = length(FL),
                mean = mean(FL),
                sd = sd(FL),
                se = sd/(sqrt(n)),
                lowCL= mean-1.96*se,
                uppCL= mean+1.96*se) %>%
      as.data.frame
    
    if(REL=="YES")
    {
      Mn=mean(d$mean,na.rm=T)
      d$mean=d$mean/Mn
      d$low95=d$lowCL/Mn
      d$up95=d$uppCL/Mn
    }
    all.yrs=seq(YEAR[1],d$year[length(d$year)])
    msn.yr=all.yrs[which(!all.yrs%in%d$year)]  
    msn.yr.dummy=d[1:length(msn.yr),]
    msn.yr.dummy[,]=NA
    msn.yr.dummy$year=msn.yr
    d=rbind(d,msn.yr.dummy)
    d=d[order(d$year),]
    return(d)
  }
  INDEX.size=PRED.size
  LimX=c(min(DATA.list$`Sandbar shark`$year),max(DATA.list$`Sandbar shark`$year))
  fn.fig("Paper/Figure 5",2400,2400)
  par(mfcol=n2mfrow(N.species.size),mai=c(.3,.38,.15,.25),oma=c(2,1.25,1,2),las=1,mgp=c(.04,.6,0))
  for(i in 1:N.species.size)
  {
    dummy=PRED.size[[i]]
    #names(dummy)[5:6]=c("lower.CL","upper.CL")
    dummy$yr=dummy$Mid.Lat
    dummy$MEAN=dummy$emmean
    dummy$UP=dummy$upper.CL
    dummy$LOW=dummy$lower.CL
    dummy$CV=100*dummy$SE/dummy$MEAN
    INDEX.size[[i]]=fun.plot.yr.pred(dummy,X="year",normalised="YES",REV="NO",n.seq=2,YLIM=NULL,XLIM=LimX,Type="points")
    
    #add observed mean size and size at maturity
    par(new = T)
    Nml.size=fn.nmnl.size(dat=Store.size[[i]]$DAT,REL="NO")
    with(Nml.size,plot(year+Plus, mean, pch=16, cex=1.25, col="grey70",axes=F,
                       xlab=NA, ylab=NA,ylim=c(0,max(Nml.size$uppCL,na.rm=T)),xlim=LimX))
    # with(Nml.size,points(year+Plus, mean, pch=16, cex=1.25, col="grey70"))
    suppressWarnings(with(Nml.size,arrows(x0=year+Plus, y0=lowCL,x1=year+Plus, y1=uppCL, 
                                     code=3, angle=90, length=0.025, col="grey70")))
    axis(side = 4,cex.axis=1.25)
    abline(h=Fem.Size.Mat[[match(TARGETS.name[i],names(Fem.Size.Mat))]],col="grey40",lwd=1.5,lty=3)   
    
    
    mtext(names(PRED.z.size)[i],3,line=.05,cex=1.1)
  }
  mtext("Relative size",2,outer=T,line=-0.5,cex=1.5,las=3)
  mtext("Year",1,outer=T,line=0.5,cex=1.5)
  plot.new()
  legend('topright',c("standardised","nominal"),bty='n',cex=1.5,pch=19,col=c("black","grey70"))
  mtext("Size (cm)",4,outer=T,line=0.5,cex=1.5,las=3,col="grey70")
  dev.off()
  
  
  #5. Export Sandbar and Dusky sharks index   
    #Fixed stations
  for(i in 1:length(INDEX.size)) write.csv(INDEX.size[[i]],paste(hnd.indx,names(INDEX.size)[i],'/',names(INDEX)[i],".Srvy.FixSt_size.csv",sep=""),row.names=F)
  

  
  #6. create dusky and sandbar figures for RAR
    #fixed stations
  San.dusky.species=match(c("Sandbar shark","Dusky shark"),names(INDEX.size))
  fn.fig("Paper/Figure 5_sandbar_dusky",1200,2400)
  par(mfcol=c(2,1),mai=c(.35,.55,.25,.1),oma=c(1.5,1.25,.1,.1),las=1,mgp=c(.04,.6,0))
  for(i in San.dusky.species)
  {
    dummy=PRED.size[[i]]
    #names(dummy)[5:6]=c("lower.CL","upper.CL")
    dummy$yr=dummy$Mid.Lat
    dummy$MEAN=dummy$emmean
    dummy$UP=dummy$upper.CL
    dummy$LOW=dummy$lower.CL
    dummy$CV=100*dummy$SE/dummy$MEAN
    INDEX.size[[i]]=fun.plot.yr.pred(dummy,X="year",normalised="YES",REV="NO",n.seq=2,YLIM=NULL,XLIM=LimX,Type="points")
    mtext(names(PRED.z.size)[i],3,line=.05,cex=1.1)
  }
  mtext("Relative size",2,outer=T,line=-0.5,cex=1.5,las=3)
  mtext("Year",1,outer=T,line=.25,cex=1.5)
  dev.off()
}


# REPORT_Ecosystem indicators FROM FISHERY INDEPENDENT SURVEYS----------------------------------------------------------------------
if(Do.ecosystems=="YES")
{
  if(do.glm.ecos=="YES")    
  {
    setwd("C:/Matias/Analyses/Surveys/Naturaliste_longline/outputs/Ecosystems")
    
    #explore model fit
    fn.fig("Model.fits",2000,2400)
    smart.par(n.plots=4*length(resp.vars),MAR=c(1,3.5,2.7,.1),OMA=c(3,1,.1,3),MGP=c(1.8,.8,0))
    for(i in 1:length(Store.mod.out.observer)) 
    {
      plot(Store.mod.out.observer[[i]]$model)
      mtext(names(Store.mod.out.observer)[i],4,las=3,line=0.75)
    }
    dev.off()  
    
    
    #paper outputs
    setwd(paste(getwd(),"/paper",sep=""))
      #export anova tables as word doc
    Export.tbl(WD=getwd(),Tbl=TABL.ecosys,Doc.nm="Anova.table.Ecosystem",caption=NA,paragph=NA,
               HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
               Zebra='NO',Zebra.col='grey60',Grid.col='black',
               Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman",
               HDR.names=c('TERM', resp.vars),
               HDR.span=c(1,rep(2,length(resp.vars))),
               HDR.2nd=c("",rep(c("P","% dev. expl."),length(resp.vars))))
    
    
      #Plot diversity and ecosystem indices
    #note: present in relative terms as we are interested in the trend
    Mst.cmn=function(D)
    {
      D=table(D)
      return(names(sort(D)[length(D)]))
    }
    
    #Function for expressing index in relative terms
    fn.relative=function(D)
    {
      Mn=mean(D$MEAN)
      D$LOW=D$LOW/Mn
      D$UP=D$UP/Mn
      D$MEAN=D$MEAN/Mn
      return(D)
    }
    
    #function for making model predictions with CI form Monte Carlo simulations   
    fun.plt.pred=function(d,Show.pred,normalised,PredictorS,MAIN,log.var,Cx,YLIM,Cx.axs)
    {
      #create new data
      Dat=d$data
      Other.preds=c(PredictorS[-match(Show.pred,PredictorS)])
      if(!is.na(OFFSETT))Other.preds=c(PredictorS[-match("year",PredictorS)],"log.EFFORT")
      
      idPred=match(Show.pred,names(Dat))
      if(is.factor(Dat[,idPred]))
      {
        NEWDATA=data.frame(x=sort(unique(as.character(Dat[,idPred]))))
        NEWDATA[,1]=factor(NEWDATA[,1],levels=levels(Dat[,idPred]))
        names(NEWDATA)=Show.pred
      }else
      {
        NEWDATA=data.frame(x=sort(seq(min(Dat[,idPred]),max(Dat[,idPred]),by=BY)))
        names(NEWDATA)=Show.pred
      }
      
      d.frm=t(matrix(rep(NA,length(Other.preds))))
      colnames(d.frm)=Other.preds
      d.frm=as.data.frame(d.frm)
      for(q in 1:length(Other.preds))  
      {
        Vec=Dat[,match(Other.preds[q],names(Dat))]
        CLs=class(Vec)
        if(CLs=='factor') Value=factor(Mst.cmn(Vec),levels=levels(Vec))
        if(CLs=='numeric') Value=mean(Vec)
        d.frm[,q]=Value
      }                  
      NEWDATA=cbind(NEWDATA,d.frm)
      Covar.pos=as.matrix(vcov(d$model))
      set.seed(999);Pos.pars.rand=rmvnorm(niter,mean=coef(d$model),sigma=Covar.pos)    
      MC.preds=matrix(nrow=niter,ncol=length(NEWDATA[,match(Show.pred,names(NEWDATA))]))
      
      for(n in 1:niter)
      {
        model=d$model
        model$coefficients=Pos.pars.rand[n,]
        a=predict(model,newdata=NEWDATA,type='response',se.fit=T)
        if(log.var=="YES") Pred=exp(a$fit+(a$se.fit^2)/2)  #apply bias correction for log transf
        if(log.var=="NO") Pred=a$fit
        MC.preds[n,]=Pred
      }
      
      PRD=data.frame(X=NEWDATA[,match(Show.pred,names(NEWDATA))],
                     MEAN=colMeans(MC.preds,na.rm=T),
                     LOW=apply(MC.preds, 2, function(x) quantile(x, 0.025,na.rm=T)),
                     UP=apply(MC.preds, 2, function(x) quantile(x, 0.975,na.rm=T)))
      names(PRD)[1]=Show.pred
      
      #standardise to a mean score of 1
      if(normalised=="YES") PRD=fn.relative(PRD)
      XX=as.numeric(as.character(PRD[,match(Show.pred,names(PRD))]))
      MeAn=PRD$MEAN
      UppCI=PRD$UP
      LowCI=PRD$LOW
      dat.plt=data.frame(yr=XX,MeAn=MeAn,UppCI=UppCI,LowCI=LowCI)
      
      #add missing years
      if(Show.pred=="year")
      {
        mis.yr=seq(ALL.yrs[1],ALL.yrs[length(ALL.yrs)])
        mis.yr=mis.yr[which(!mis.yr%in%dat.plt$yr)]
        if(length(mis.yr)>0)
        {
          dummy=dat.plt[1:length(mis.yr),]
          dummy[,]=NA
          dummy$yr=mis.yr
          dat.plt=rbind(dat.plt,dummy)
        }
        dat.plt=dat.plt[order(dat.plt$yr),]
      }
      
      if(is.null(YLIM)) YLIM=c(min(dat.plt$LowCI,na.rm=T),max(dat.plt$UppCI,na.rm=T))
      
      with(dat.plt,plot(yr,MeAn,pch=19,main="",xlab="",ylab="",
                        cex=1.25,cex.axis=1.25,ylim=YLIM))
      #with(dat.plt,plot(yr,MeAn,pch=19,main=MAIN,xlab="",ylab="",
      #                  cex=1.25,xaxt="n",cex.axis=1.25,ylim=YLIM))
      with(dat.plt,arrows(x0=yr, y0=LowCI, x1=yr, y1=UppCI,code = 3,angle=90,length=.025))
      #axis(1,dat.plt$yr,F,tck=-0.02)
      #with(dat.plt,axis(1,seq(yr[1],yr[length(yr)],5),F,tck=-0.05))
      #if(i %in%c(3,6)) with(dat.plt,axis(1,seq(yr[1],
      #     yr[length(yr)],5),seq(yr[1],yr[length(yr)],5),tck=-0.05,cex.axis=1.25))
      mtext(MAIN,3,cex=1.25)
    }
    ALL.yrs=as.numeric(as.character(levels(DATA.eco$year)))
    Main.title=resp.vars
    
    fn.fig("Div & Eco indices_normalised_year",1800,2400)  
    smart.par(n.plots=length(resp.vars),MAR=c(1,3.5,3,.1),OMA=c(3,1,.1,2),MGP=c(1,.8,0))
    system.time(for(i in 1:length(Store.mod.out.observer))  #takes 0.02 secs per iteration
    {
      fun.plt.pred(d=Store.mod.out.observer[[i]],Show.pred="year",normalised="YES",
                      PredictorS=Predictors,MAIN=Main.title[i],log.var=Res.var.in.log[i],
                      Cx=1.5,YLIM=NULL,Cx.axs=1.5)
    })
    mtext("Year",1,outer=T,cex=1.5,line=1.5)
    mtext("Relative value",2,-.75,outer=T,cex=1.5,las=3)
    dev.off()
    
    fn.fig("Div & Eco indices_normalised_BOTDEPTH.mean",1800,2400)  
    BY=10
    smart.par(n.plots=length(resp.vars),MAR=c(1,3.5,3,.1),OMA=c(3,1,.1,2),MGP=c(1,.8,0))
    system.time(for(i in 1:length(Store.mod.out.observer))  #takes 0.02 secs per iteration
    {
      fun.plt.pred(d=Store.mod.out.observer[[i]],Show.pred="BOTDEPTH.mean",normalised="YES",
                   PredictorS=Predictors,MAIN=Main.title[i],log.var=Res.var.in.log[i],
                   Cx=1.5,YLIM=NULL,Cx.axs=1.5)
    })
    mtext("Depth (m)",1,outer=T,cex=1.5,line=1.5)
    mtext("Relative value",2,-.75,outer=T,cex=1.5,las=3)
    dev.off()
    
    fn.fig("Div & Eco indices_normalised_Mid.Lat.mean",1800,2400) 
    BY=.5
    smart.par(n.plots=length(resp.vars),MAR=c(1,3.5,3,.1),OMA=c(3,1,.1,2),MGP=c(1,.8,0))
    system.time(for(i in 1:length(Store.mod.out.observer))  #takes 0.02 secs per iteration
    {
      fun.plt.pred(d=Store.mod.out.observer[[i]],Show.pred="Mid.Lat.mean",normalised="YES",
                   PredictorS=Predictors,MAIN=Main.title[i],log.var=Res.var.in.log[i],
                   Cx=1.5,YLIM=NULL,Cx.axs=1.5)
    })
    mtext("Latitude (?S)",1,outer=T,cex=1.5,line=1.5)
    mtext("Relative value",2,-.75,outer=T,cex=1.5,las=3)
    dev.off()
  }
}


# REPORT_MULTIVARIATE ANALYSIS OF CATCH COMPOSITION FROM FISHERY INDEPENDENT SURVEYS----------------------------------------------------------------------
if(Do.multivariate=="YES")
{
  Hndl="C:/Matias/Analyses/Surveys/Naturaliste_longline/outputs/Multivariate/"
  for(s in 1:length(STore.multi.var.trad)) 
  {
    with(STore.multi.var.trad[[s]],
      {
      fn.display.multivar(d=d,IDVAR=IDVAR,Predictors=Predictors,MDS=MDS,
                          Permanova.table=Permanova.table,
                          permanova.pairwise=permanova.pairwise,
                          Simper=Simper,NM=names(STore.multi.var.trad)[s],hndl=Hndl,cexMDS=.75)
      })
  }

}



