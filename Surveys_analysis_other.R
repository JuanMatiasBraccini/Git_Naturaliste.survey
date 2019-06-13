####   SCRIPT FOR ANALYSING NATURALISTE SURVEY DATA ####

# INDEX :
          # -- DATA SECTION --
          # -- CONTROL SECTION --
          # -- PROCEDURE SECTION --
            #1. IMMEDIATE POST CAPTURE MORTALITY
            #2. REPRODUCTION 
            #3. DIET



#NOTES:   shot is defined by SHEET_NO
#         Periods: prior to 2001, drumlines; 2001 onwards, longlines
#         fixed stations: 2001:2008, 2015: onwards
#         Timeline: 
#             First period: targeted at sandbar sharks for tagging
#             Then WAMSI objectivies
#             Then FRDC tagging of dusky sharks
#             Then standard survey
#         PCS: use Release Condition: if 1:3, then Alive, if blank then DEAD. Cross check
#             no blank with Tag number; no diet/reprod data with Release Condition


rm(list=ls(all=TRUE))
require(lubridate)  #for dates manipulation
library(geosphere)
library(ggplot2)
library(pscl)  #zero-inflated models
library(MCMCglmm) #zero-inflated models
library(glmmADMB) #zero-inflated models
library(tweedie)
require(statmod) # Provides  tweedie  family functions
library(bbmle)  #AIC table
library(lmtest) #likelihood ratio tests
library(mgcv)
library(PBSmapping)
data(worldLLhigh)
require(plotrix) #multihistogram
library(gplots)  #annotated boxplot
#install.packages("C:/Matias/R/COZIGAM_2.0-3.tar.gz",type="source")
#library(COZIGAM)
library(dplyr)
library(data.table)
library(sizeMat)

#Define user
User="Matias"
#User="Dani"

#handy function for exporting figures
fn.fig=function(NAME,Width,Height)
{
  if(Do.tiff=="YES") tiff(file=paste(NAME,".tiff",sep=""),width=Width,height=Height,units="px",res=300,compression="lzw")
  if(Do.jpeg=="YES") jpeg(file=paste(NAME,".jpeg",sep=""),width=Width,height=Height,units="px",res=300)
}

#choose if doing .jpeg or .tiff figures
Do.jpeg="YES"
Do.tiff="NO"


# -- DATA SECTION --

#Source Sharks data base 
if(User=="Matias") source("C:/Matias/Analyses/SOURCE_SCRIPTS/Source_Shark_bio.R")
if(User=="Dani") source("C:/Users/ddw/Documents/Projects/Naturaliste_longline/Source codes/Source_Shark_bio.R")
rm(DATA)

#Comments data  
if(User=="Matias") 
{
  Comnt_Alive=read.csv("C:/Matias/Data/Naturaliste/Comments_Alive.csv")
  Comnt_Dead=read.csv("C:/Matias/Data/Naturaliste/Comments_Dead.csv")
}
if(User=="Dani") 
{
  Comnt_Alive=read.csv("C:/Users/ddw/Documents/Projects/Naturaliste_longline/Source codes/Comments_Alive.csv")
  Comnt_Dead=read.csv("C:/Users/ddw/Documents/Projects/Naturaliste_longline/Source codes/Comments_Dead.csv")  
}



#Diet
#note: Diet is a csv file that must be manually updated each year. Data is presence/absence
if(User=="Matias")
{
  Diet=read.csv("C:/Matias/Data/Naturaliste/Look.at.diet_1.csv")
  Predator_TL=read.csv("C:/Matias/Data/Naturaliste/Predator_TL.csv")
  Prey_TL=read.csv("C:/Matias/Data/Naturaliste/Prey_TL.csv")
  Prey_cnvrsn=read.csv("C:/Matias/Data/Naturaliste/Prey.csv")
  
}
if(User=="Dani")   Diet=read.csv("C:/Users/ddw/Documents/Projects/Naturaliste_longline/Diet/Look.at.diet_1.csv")



# -- CONTROL SECTION --
Current.Yr=2017


#control if doing survival analysis
#Do.survival="YES"
Do.survival="NO"

#control if doing reproduction and diet analysis
Do.rep="NO"
Do.diet="YES"

#control if doing ecosystem analysis
Do.ecosystem="NO"



# -- PROCEDURE SECTION --

#Create data sets for PCS and biological analysis 
#note: we assume that any shark that was dissected was not sacrified so it was dissected 
#       because it arrived dead on deck

  #is there biological data of any sort?
DATA.bio$Biol.dat=with(DATA.bio,ifelse(
  !is.na(CLASPLENTH)| !is.na(CLASP_CALC) | !is.na(GON_STAGE) | !is.na(RUN_SPERM)  | !is.na(MAXOVRYDIA)
  | !is.na(NO_YOLKOVA) | !is.na(UTERINESTG) | !is.na(NO_EMBRYOS) | !is.na(NO_UNDEVELOPED) 
  | !is.na(EMBLEN_1)   | !is.na(STMCH_FULL)| !is.na(STMCH_CONT),"YES","NO"))  

#was the shark tagged?
DATA.bio$Tagged=with(DATA.bio,ifelse(!is.na(ATAG_NO) | !is.na(DART_TAG_NO)  
                                     | !is.na(FINTAG_NO)| !is.na(FINTAG_2),"YES","NO")) 

#is there a sample?
DATA.bio$Sample=with(DATA.bio,ifelse(!is.na(BAG_NO) ,"YES","NO"))  

DATA.bio=subset(DATA.bio,!SPECIES=="XX")

  #was the shark dissected?
DATA.bio$Dissected=with(DATA.bio,ifelse(
  !is.na(GON_STAGE) | !is.na(MAXOVRYDIA)
  | !is.na(NO_YOLKOVA) | !is.na(UTERINESTG) | !is.na(NO_EMBRYOS) | !is.na(NO_UNDEVELOPED)
  | !is.na(EMBLEN_1)   | !is.na(STMCH_FULL)| !is.na(STMCH_CONT),"YES","NO"))  

DATA.bio$Method=with(DATA.bio,ifelse(Method=="LLL","LL",Method))


  #Create PCS data
Yr.start=2015  #First year for which we can confirm sharks were not sacrificed

DATA.PCS=DATA.bio
DATA.PCS$ALIVE=DATA.PCS$"RELEASE CONDITION"
DATA.PCS$ALIVE=with(DATA.PCS,ifelse(ALIVE%in%1:3,1,ifelse(ALIVE==4,0,ALIVE)))
DATA.PCS$ALIVE=with(DATA.PCS,ifelse(is.na(ALIVE) & Dissected=="YES",0,
              ifelse(Tagged=="YES" & is.na(ALIVE),1,ALIVE))) 

#check comments
#CMNTS=subset(DATA.PCS,is.na(ALIVE) & year>=Yr.start & Method=="LL",select=c(SHEET_NO,SPECIES,ALIVE,COMMENTS))
#CMNTS=as.character(sort(unique(CMNTS$COMMENTS)))
#write.csv(as.data.frame(CMNTS),"C:/Users/myb/Desktop/New folder (2)/cmnts.csv",row.names=F)

DATA.PCS$ALIVE=with(DATA.PCS,ifelse(COMMENTS%in%Comnt_Alive$COMMENTS,1,
                   ifelse(COMMENTS%in%Comnt_Dead$COMMENTS,0,ALIVE)))



Tab=table(DATA.PCS$ALIVE,DATA.PCS$SPECIES,useNA='ifany')
DATA.PCS=subset(DATA.PCS, year>=Yr.start & Method=="LL")  
Tab1=table(DATA.PCS$ALIVE,DATA.PCS$SPECIES,useNA='ifany')

DATA.PCS=subset(DATA.PCS,!is.na(ALIVE))
Tab1=table(DATA.PCS$ALIVE,DATA.PCS$SPECIES,useNA='ifany')


Prop.Dead=round(Tab1[1,]/colSums(Tab1),2)


DATA.bio=subset(DATA.bio,Biol.dat=="YES")
Tab.sp=sort(table(as.character(DATA.bio$COMMON_NAME)))



#1. IMMEDIATE POST CAPTURE MORTALITY
if(Do.survival=="YES")  
{
  #Combine blacktips
  DATA.PCS$SPECIES=with(DATA.PCS,ifelse(SPECIES=="BL","BT",SPECIES))
  DATA.PCS$COMMON_NAME=with(DATA.PCS,ifelse(SPECIES=="Common blacktip shark","Blacktip sharks",COMMON_NAME))
  DATA.PCS$SCIENTIFIC_NAME=with(DATA.PCS,ifelse(SPECIES=="Carcharhinus limbatus","Carcharhinus limbatus/tilstoni",SCIENTIFIC_NAME))
  
  N.min=5  #minimum number of records per species for PCS analysis
  
  if(User=="Matias") setwd("C:/Matias/Analyses/Surveys/Naturaliste_longline/outputs/Survival")
  if(User=="Dani") setwd("C:/Users/ddw/Documents/Projects/Naturaliste_longline/outputs/Survival")
  

  #HABITAT 
  Demersal<- c("BN","BU","CD","GG","GM","GN","LE","LP","MI","PE","PJ","PN","SD", ## Demersal living on or near the bottom
               "SE","SI","SO","SR","ST","TN","WE","WG","WH","ZE")
  Benthic <- c("FR","PZ","SH","WD","WR","WW") ## Demersal can be divided into benthic and bentho-pelagic
  Pelagic <- c("BL","BT","BW", "CP", "GR", "HG","HS","HZ","LG","MS","TG","TK","WP") ## Pelagic or oceanic sharks live in the open waters of the seas and oceans. 
  
  # SPIRACLE
  Large <- c("PZ","SD","SH","TG","WD","WG","WW","GM") # GM, WG - Gummy shark has moderate sized spiracles
  Small <- c("TN","PJ","ZE")
  Absent<- c("BL","BT","BW","GN","GR","HG","HS","HZ","LE","LG","MI","PE","SO","TK") 
  
  names(DATA.PCS)[match("RELEASE CONDITION",names(DATA.PCS))]="CONDITION"
  
  Shots=length(unique(DATA.PCS$SHEET_NO))
  
  #Figure 1
  do.map="NO"
  if(do.map=="YES")
  {
    library(PBSmapping)
    data(worldLLhigh)
    plt.dat=DATA.PCS[!duplicated(DATA.PCS$SHEET_NO),]
    
    fn.fig("Map",1600,2400)
    par(mai=c(.1,.1,.6,.1),oma=c(1,1,1,.1),mgp=c(.04,.6,0))
    plotMap(worldLLhigh, xlim= c(112,129),ylim=c(-36,-13),plt = c(.1, 1, 0.075, 1),
            col="grey90",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
    box()
    text(122,-23,("Western"),col="black", cex=3)
    text(122,-26,("Australia"),col="black", cex=3)
    points(plt.dat$Mid.Long,plt.dat$Mid.Lat,cex=1.5,
           col=rgb(.2,.2,.2,alpha=0.25),pch=19)
    mtext("Longitude (°E)",1,-.8,outer=T,cex=2)
    mtext("Latitude (°S)",2,line=-0.6,outer=T,cex=2)
    axis(1,seq(112,129,2),seq(112,129,2),cex.axis=1.25)
    axis(2,seq(-13,-36,-2),seq(13,36,2),las=1,cex.axis=1.25)
    dev.off()
    
  }
   
  

  #Exclude <- c("WR", "ST","SE")    
  Exclude=NA
  PCS <- DATA.PCS%>%
    select(year, Day, Month, year, SPECIES, SCIENTIFIC_NAME, COMMON_NAME,CONDITION, SEX, FL, TL, ALIVE, 
           SOAK.TIME, BOTDEPTH, TEMP, BLOCK, Mid.Lat, Mid.Long, Lat.round, Long.round, zone, SKIPPER)%>%
    mutate(HABITAT=ifelse(SPECIES %in% Demersal, "D",
                          ifelse(SPECIES %in% Pelagic, "P",
                                 ifelse(SPECIES %in% Benthic, "B", "NA")))) %>%
    mutate(SPIRACLE=ifelse(SPECIES %in% c(Large,Small), "1",
                             ifelse(SPECIES %in% Absent, "0", "NA"))) %>%
    filter(!(SPECIES %in% Exclude),
           !is.na(SEX))%>%
    mutate(SPECIES=factor(SPECIES) ,
           SEX=factor(SEX) , 
           Month=factor(Month),
           HABITAT=factor(HABITAT),
           SPIRACLE=factor(SPIRACLE,levels=c("1","0")))
  
  
  
  plot(PCS$SOAK.TIME,PCS$Lat.round)
  PCS$SKIPPER=as.factor(as.character(PCS$SKIPPER))
  plot(PCS$SKIPPER,PCS$Lat.round)
  
  plot(PCS$Month,PCS$Lat.round)
  boxplot(ALIVE~SEX,data=PCS)
  
  tab=with(PCS,table(ALIVE,SPIRACLE))
  tab=with(PCS,table(ALIVE,HABITAT))
  


#Table 1: Observed PCS by species
  SURVIVAL <- PCS%>%
    group_by(SPECIES)%>%
    summarise(Dead=length(ALIVE[ALIVE==0]),
              Alive=length(ALIVE[ALIVE==1]))%>%
    mutate(Total=Dead+Alive)%>%
    mutate(Prop.dead=Dead/Total)%>%
    mutate(Prop.dead.SE=sqrt(Prop.dead*(1-Prop.dead)/Total))%>%
    arrange(., Prop.dead)%>%
    as.data.frame
  Add.sp=PCS %>% select(SPECIES,SCIENTIFIC_NAME,COMMON_NAME) %>%
    distinct(SPECIES,SCIENTIFIC_NAME,COMMON_NAME)
  Add.sp=Add.sp[!duplicated(Add.sp$SPECIES),]
  SURVIVAL=merge(Add.sp,SURVIVAL,by="SPECIES") %>%
            filter(!is.na(SCIENTIFIC_NAME))
  
  #some manipulations
  PCS$SPECIES=as.character(PCS$SPECIES)
  PCS$FL=with(PCS,ifelse(is.na(FL) & SPECIES %in%c("SH","ZE"),TL,FL))  #explain in Manuscript that TL was used for these 2 species
  PCS.glm=PCS
  PCS.glm$Size.range=cut(PCS.glm$FL,5)
  PCS.glm=subset(PCS.glm,!SPECIES=="PF")

  #Select species with at least 5 records
  TAB=table(PCS.glm$SPECIES)
  This.sp=names(subset(TAB,TAB>N.min))
  
  #select species with contrast (i.e. dead and alive records)
  Tab.contrast=with(subset(PCS.glm,SPECIES%in%This.sp),table(ALIVE,SPECIES))
  Tab.contrast[Tab.contrast>0]=1
  This.sp=names(which(colSums(Tab.contrast)>1))
  
  PCS.glm=subset(PCS.glm,SPECIES%in%This.sp)
  PCS.glm$SPECIES=factor(as.character(PCS.glm$SPECIES))

  
  
# GLM model was used to test the effects of predictors on at vessel mortality.
#note: no contrast in habitat or spiracles so don't use in model
  
  
  #add missing FL for dusky
  Add.misn.FL.dusky=round(mean(subset(PCS.glm,
            year==2017 & Day==16 & Month==5 & SPECIES=="BW")$FL,na.rm=T))
  PCS.glm$FL=with(PCS.glm,ifelse(SPECIES=="BW" & is.na(FL) &
            year==2017 & Day==16 & Month==5,Add.misn.FL.dusky,FL))
  
  #put continuous vars in log space
  PCS.glm$log_FL=log(PCS.glm$FL)
  PCS.glm$log_Lat.round=log(-PCS.glm$Lat.round)
  PCS.glm$log_BOTDEPTH=log(PCS.glm$BOTDEPTH)
  PCS.glm$log_SOAK.TIME=log(PCS.glm$SOAK.TIME)
  
  #define response var
  PCS.glm$DEAD=1-PCS.glm$ALIVE
  
  
  #First step. Look at differences among species 
  model1 <- glm(DEAD~SPECIES, family='binomial', data=PCS.glm)
  
    #Export ANOVA table
  fn.anov=function(model)
  {
    Anova.tab=anova(model, test = "Chisq")
    n=2:length(Anova.tab$Deviance)
    Term.dev.exp=100*(Anova.tab$Deviance[n]/model$null.deviance)
    names(Term.dev.exp)=rownames(Anova.tab)[n]
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
    return(Anova.tab)
  }
  write.csv(fn.anov(model1),"Anova.tab_species.only.csv")
  
    #Predict species effect: no point for Species only as it's the same as the observed
  fn.pred=function(mod,NEWDATA)
  {
    SP.pred.avg=predict(mod,type="response",newdata=NEWDATA,se.fit=T)
    NEWDATA$"Predicted proportion dead"=SP.pred.avg$fit
    NEWDATA$SE=SP.pred.avg$se.fit
    return(NEWDATA)
  }
  
  Table1=SURVIVAL%>%arrange(., Prop.dead,desc(Total))%>%
    rename('Observed proportion dead'=Prop.dead)%>%
    rename('SE obs.'=Prop.dead.SE)%>%
    rename('Common name'=COMMON_NAME)%>%
    rename('Scientific name'=SCIENTIFIC_NAME)%>%
    mutate_if(is.numeric, round, 3)
  write.csv(Table1%>%select(-SPECIES),"Table1.csv",row.names=F)
  
  
  
  
  #Second, fit full model for selected species
  library(gridExtra)
  pdf("Predictor.range_Analysed.species.pdf")
  for(s in 1:length(This.sp))
  {
    d=subset(PCS.glm,SPECIES==This.sp[s])   #aca: issues with missiing length, depth, etc
    d$dead=as.factor(d$DEAD)
    table(d$SEX)
    p1=ggplot(d, aes(x=exp(log_FL), fill=dead)) + xlab("FL") +
      geom_histogram(binwidth=10, position="dodge")
    p2=ggplot(d, aes(x=exp(log_BOTDEPTH), fill=dead)) +xlab("Depth") +
      geom_histogram(binwidth=10, position="dodge") + ggtitle(d$COMMON_NAME[1])
    p3=ggplot(d, aes(x=exp(log_SOAK.TIME), fill=dead)) +xlab("Soak") +
      geom_histogram(binwidth=1, position="dodge")
    p4=ggplot(d, aes(x=exp(log_Lat.round), fill=dead)) +xlab("Lat") +
      geom_histogram(binwidth=, position="dodge")
    grid.arrange(p1, p2,p3,p4, nrow = 2,ncol=2)
  }
  dev.off()

  #keep only species with a good contrast in predictor
  Spi=c("TK","MI","SO")
  names(Spi)=c("Sandbar shark","Milk shark","Spot-tail shark")
  # This.sp=c("BT","TK","MI","SO")
  # names(This.sp)=c("Blacktip sharks","Sandbar shark","Milk shark","Spot-tail shark")
  
  
  MOD=vector('list',length(Spi))
  names(MOD)=Spi
  ANOV=MOD
  for(s in 1:length(Spi))
  {
    Dat.mod=subset(PCS.glm,SPECIES==Spi[s])
    Dat.mod$SPECIES=factor(Dat.mod$SPECIES)
    if(Spi[s]=="BT")MOD[[s]]=glm(DEAD~log_FL+log_BOTDEPTH+log_Lat.round,
                             family='binomial', data=Dat.mod)else
                    MOD[[s]]=glm(DEAD~log_FL+log_BOTDEPTH+log_SOAK.TIME+log_Lat.round+SEX,
                             family='binomial', data=Dat.mod) 
    ANOV[[s]]=fn.anov(MOD[[s]])
  }
  
  
  #Export anova table
   write.csv(do.call(rbind,ANOV),"Anova.tab_all_terms.csv")
  

  #Model predictions
  fn.rango=function(D)
  {
    a=subset(PCS.glm,SPECIES==D)
    FFL=sort(unique(round(a$FL/10)*10))
    FFL=seq(FFL[1],FFL[length(FFL)],10)
    return(data.frame(SPECIES=factor(D,levels=Spi),log_FL=log(FFL)))
  }
  fn.rango.depth=function(D)
  {
    a=subset(PCS.glm,SPECIES==D)
    FFL=sort(unique(round(a$BOTDEPTH/10)*10))
    FFL=seq(FFL[1],FFL[length(FFL)],10)
    return(data.frame(SPECIES=factor(D,levels=Spi),log_BOTDEPTH=log(FFL)))
  }
  fn.rango.soak=function(D)
  {
    a=subset(PCS.glm,SPECIES==D)
    FFL=sort(unique(round(a$SOAK.TIME)))
    FFL=seq(FFL[1],FFL[length(FFL)],.1)
    return(data.frame(SPECIES=factor(D,levels=Spi),log_SOAK.TIME=log(FFL)))
  }
  Rng=vector('list',length(Spi))
  names(Rng)=Spi
  FL.effect=Depth.effect=Soak.effect=Rng
  for(s in 1:length(Spi))
  {
    Dat.mod=subset(PCS.glm,SPECIES==Spi[s])
    
    #Predict for range of size
    NEWDATA=fn.rango(D=Spi[s])
    NEWDATA=cbind(NEWDATA,    
                  SEX=factor("F",levels=levels(Dat.mod$SEX)),
                  log_Lat.round=mean(log(-Dat.mod$Lat.round),na.rm=T),
                  log_BOTDEPTH=mean(log(Dat.mod$BOTDEPTH),na.rm=T),
                  log_SOAK.TIME=mean(log(Dat.mod$SOAK.TIME),na.rm=T))
    FL.effect[[s]]=fn.pred(mod=MOD[[s]],NEWDATA=NEWDATA)
    
    
    #Predict depth effect
    NEWDATA=fn.rango.depth(D=Spi[s])
    NEWDATA=cbind(NEWDATA,    
                  SEX=factor("F",levels=levels(Dat.mod$SEX)),
                  log_Lat.round=mean(log(-Dat.mod$Lat.round),na.rm=T),
                  log_FL=mean(log(Dat.mod$FL),na.rm=T),
                  log_SOAK.TIME=mean(log(Dat.mod$SOAK.TIME),na.rm=T))
    Depth.effect[[s]]=fn.pred(mod=MOD[[s]],NEWDATA=NEWDATA)
    
    
    #Predict soak time effect
    NEWDATA=fn.rango.soak(D=Spi[s])
    NEWDATA=cbind(NEWDATA,    
                  SEX=factor("F",levels=levels(Dat.mod$SEX)),
                  log_Lat.round=mean(log(-Dat.mod$Lat.round),na.rm=T),
                  log_FL=mean(log(Dat.mod$FL),na.rm=T),
                  log_BOTDEPTH=mean(log(Dat.mod$BOTDEPTH),na.rm=T))
    Soak.effect[[s]]=fn.pred(mod=MOD[[s]],NEWDATA=NEWDATA)
    
  }
  
  #Figure 1
  fn.plt.fig1=function(D,var,KL,delta)
  {
    id=match(var,names(D))
    a=D[order(D[,id]),]
    a[,id]=exp(a[,id])
    a$pred=a$'Predicted proportion dead'
    points(a[,id]+delta,a$pred,pch=19,cex=2,col=CL[s])
    segments(a[,id]+delta,a$pred,a[,id]+delta,a$pred+1.96*a$SE,lwd=2,col=KL)
    segments(a[,id]+delta,a$pred,a[,id]+delta,a$pred-1.96*a$SE,lwd=2,col=KL)
  }
  
  CL=c("black","grey50","grey70","grey85")
  
  tiff(file="Figure1.tiff",width = 1400, height = 2400,units = "px", res = 300, compression = "lzw")
  par(mfcol=c(3,1),mai=c(.3,.175,.4,.1),oma=c(3,3.5,.01,.1),las=1,mgp=c(1,.6,0))
  
  #Predict FL effect 
  dummy=do.call(rbind,FL.effect)
  Dlt=c(1,2,-2,0)
  plot(1:10,ylim=c(0,1),xlim=c(min(exp(dummy$log_FL)),max(1.025*exp(dummy$log_FL))),ylab="",xlab="",cex.axis=1.5)
  for(s in 1:length(Spi))fn.plt.fig1(D=FL.effect[[s]],var='log_FL',KL=CL[s],delta=Dlt[s]) 
  mtext("Body size (cm)",1,2.5,cex=1.5)
  legend("topright",names(Spi),bty='n',pch=19,col=CL,cex=1.25)
  
  #Predict depth effect
  dummy=do.call(rbind,Depth.effect)
  Dlt=c(1,2,-2,0)
  plot(1:10,ylim=c(0,1),xlim=c(min(exp(dummy$log_BOTDEPTH)),max(1.025*exp(dummy$log_BOTDEPTH))),ylab="",xlab="",cex.axis=1.5)
  for(s in 1:length(Spi))fn.plt.fig1(D=Depth.effect[[s]],var='log_BOTDEPTH',KL=CL[s],delta=Dlt[s]) 
  mtext("Bottom depth (m)",1,2.5,cex=1.5)
  
  #Predict soak time effect
  dummy=do.call(rbind,Soak.effect)
  Dlt=c(.05,.075,.1,.125)
  plot(1:10,ylim=c(0,1),xlim=c(min(exp(dummy$log_SOAK.TIME)),max(1.025*exp(dummy$log_SOAK.TIME))),ylab="",xlab="",cex.axis=1.5)
  for(s in 1:length(Spi))fn.plt.fig1(D=Soak.effect[[s]],var='log_SOAK.TIME',KL=CL[s],delta=Dlt[s]) 
  mtext("Soak time (hrs)",1,2.5,cex=1.5)
  
  mtext("Predicted proportion dead",2,1.65,outer=T,cex=1.5,las=3)
  dev.off()
  
  
  
  
  
  #Appendix 1
  fun.his=function(d)
  {
    barplot(table(d),cex.axis=1.25,cex.names=1.25)
    box()
  }
  tiff(file="Appendix1.tiff",width = 1200, height = 2400,units = "px", res = 300, compression = "lzw")
  par(mfcol=c(3,1),mai=c(.2,.15,.3,.2),oma=c(2,3.5,.01,.1),las=1,mgp=c(1,.6,0))
  fun.his(round(PCS.glm$SOAK.TIME,1))
  mtext("Soak time (hrs)",1,2,cex=1.35)
  fun.his(10*round(PCS.glm$BOTDEPTH/10))
  mtext("Bottom depth (m)",1,2,cex=1.35)
  fun.his(2*round(-PCS.glm$Lat.round/2))
  mtext(expression(paste("Latitude (",degree,"S)")) ,1,2.5,cex=1.35)
  mtext("Frequency",2,1.5,outer=T,cex=1.5,las=3)
  dev.off()
  
  #Appendix 2
  SPI=subset(Table1,Total>=10)$SPECIES
  tiff(file="Appendix2.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")
  par(mfcol=c(4,4),mai=c(.1,.15,.485,.2),oma=c(3,3.5,.01,.1),las=1,mgp=c(1,.6,0))
  for(s in 1:length(SPI))
  {
    a=subset(PCS,SPECIES==SPI[s])
    fun.his(10*round(a$FL/10))
    Lab=unique(a$COMMON_NAME)
    Lab=gsub(" shark", "", Lab)
    if(Lab[1]=="Common blacktip") Lab="Blacktips"
    if(SPI[s]%in%c("WR","SH","ZE"))Lab=paste(Lab,"*",sep='')
    mtext(Lab,3,0.5,cex=.95)
  }
  mtext("Body size (cm)",1,1.5,cex=1.35,outer=T)
  mtext("Frequency",2,1.5,cex=1.35,outer=T,las=3)
  dev.off()
  
  
  #ACA
  #Show condition of species used in GLM
  This.sp.cond=c(This.sp,"LE","SH","PE","WR","ZE")
  A=subset(PCS,SPECIES%in%This.sp.cond & !(is.na(CONDITION)| CONDITION==4))
  Cond.tab=as.matrix(table(A$CONDITION,A$SPECIES,useNA='ifany'))
  Ns=colSums(Cond.tab)
  Cond.tab=Cond.tab/matrix(rep(colSums(Cond.tab),3),nrow=nrow(Cond.tab),byrow=T)
 
  #Correlation Condition 1 and IM
  Corr=Table1[,match(c("Common name","Observed proportion dead"),names(Table1))]
  names(Corr)[2]="Prop.dead"
  Corr=subset(Corr,Prop.dead<1)
  SNAMES=subset(PCS,SPECIES%in%This.sp.cond, select=c(SPECIES,COMMON_NAME,SCIENTIFIC_NAME))
  SNAMES=SNAMES[!duplicated(SNAMES$SPECIES),]
  Corr=merge(Corr,SNAMES,by.x="Common name",by.y="COMMON_NAME",all.y=T)
  a=data.frame(COND1=Cond.tab[1,],SPECIES=colnames(Cond.tab))
  Corr=merge(Corr,a,by="SPECIES")
  Corr=Corr[!is.na(Corr$COND1),]
    
  Corr=subset(Corr,select=c(SPECIES,Prop.dead,COND1))
  Cor=round(cor(Corr$COND1,Corr$"Prop.dead"),2)
  Corr$Prop.dead=round(Corr$Prop.dead,2)
  Corr$COND1=round(Corr$COND1,2)
  Corr$COND1=with(Corr,ifelse(Prop.dead==0,jitter(COND1,10),COND1))
  Corr$Prop.dead=with(Corr,ifelse(Prop.dead==0,jitter(Prop.dead,10),Prop.dead))
  
  Corr=Corr[order(Corr$COND1),]
 
  tiff(file="Figure2.tiff",width = 2000, height = 2400,units = "px", res = 300, compression = "lzw")
  par(mfcol=c(2,1),mar=c(2,4,1.5,.1),oma=c(2.5,.7,.5,.1),las=1,mgp=c(1,.6,0),xpd=T)
  Cond.tab=Cond.tab[,match(Corr$SPECIES,colnames(Cond.tab))]
  barplot(Cond.tab,beside=F,cex.axis=1.2,cex.names=.8,legend.text=rownames(Cond.tab),
          args.legend=list(x=16,y=1.175,bty='n',cex=1.65,horiz=T))
  mtext("Condition",2,line=2.75,cex=1.5,las=3)
  mtext("Species",1,line=2,cex=1.5)
 # text(x=seq(0.7,21,length.out=ncol(Cond.tab)),y=0.05,Ns,col="white",cex=.8)
  box()
 
   plot(Corr$COND1,Corr$"Prop.dead",
       xlab="",ylab="",
       col="transparent",ylim=c(0,1))
  text(Corr$COND1,Corr$"Prop.dead",Corr$SPECIES,cex=.8,srt=45)
  legend('topright',paste("correlation= ",Cor,sep=""),cex=1.2,bty='n')
  mtext("Proportion in conditon 1",1,line=2,cex=1.5)
  mtext("Observed proportion dead",2,line=2.75,cex=1.5,las=3)
  dev.off()
  
  
}


if(Do.rep=="YES"| Do.diet=="YES")
{
  #Exclude all bony fish from data
  n.min=5
  Exclude <- c("CD", "BN", "FR", "GB", "GG", "SB","SE", "ST", "WR") #Ostheichthyes 
  
  DATA.rep <-   filter(DATA.bio,year>=1990)%>%
    dplyr::select(SHEET_NO,date, Day, Month, year, Lat.round, Long.round, SPECIES, SEX,TL, FL, COMMON_NAME, SCIENTIFIC_NAME, UMBIL_SCAR, 
                  CLASPLENTH, CLASP_CALC, RUN_SPERM,
                  UTERINESTG, NO_EMBRYOS,EMBLEN_1, NO_UNDEVELOPED,STMCH_FULL,STMCH_CONT, GON_STAGE, NO_YOLKOVA, MAXOVRYDIA,Biol.dat)%>%
    filter(!(SPECIES %in% Exclude))%>%
    mutate(SPECIES=factor(SPECIES) ,
           SEX=factor(SEX) , 
           Month=factor(Month))
  
  STORE.l=vector('list',length(Biol.vars))
  names(STORE.l)=Biol.vars
  for(q in 1:length(Biol.vars))
  {
    d=DATA.rep[,match(c("SPECIES",Biol.vars[q]),names(DATA.rep))]  
    id=match(Biol.vars[q],names(d))
    d=d[!is.na(d[,id]),]
    TABLA=table(d$SPECIES)
    STORE.l[[q]]=TABLA
  }
  
  STORE.l=do.call(rbind,STORE.l)
  
  Dummy=STORE.l
  Dummy[Dummy<n.min]=0
  Drop=colSums(Dummy)
  Drop=subset(Drop,Drop==0)
  ID=match(names(Drop),colnames(Dummy))
  Dummy=Dummy[,-ID]
  Tabl.rep.info.sp=t(Dummy)
  
  #Select species for analysis

  Dummy.1=Dummy[-match(c("STMCH_FULL","STMCH_CONT"),rownames(Dummy)),]  
  Reprod.vars=rownames(Dummy.1)
  Reprod.vars=Reprod.vars[-match("NO_UNDEVELOPED",Reprod.vars)]  #remove NO_UNDEVELOPED because no obs for these species
  RnK.sp=sort(colSums(Dummy.1))
  RnK.sp=RnK.sp[-match(c("GM","TK","BW","WH","ES","SC"),names(RnK.sp))]    #remove commercial species (already done)
  Reprod.sp=names(RnK.sp)
  Look.at.repro=subset(DATA.bio,SPECIES%in%Reprod.sp & (!is.na(UMBIL_SCAR) | !is.na(CLASPLENTH)
                  | !is.na(CLASP_CALC) | !is.na(GON_STAGE) | !is.na(RUN_SPERM) | !is.na(MAXOVRYDIA) 
                  | !is.na(NO_YOLKOVA) | !is.na(UTERINESTG) | !is.na(NO_EMBRYOS)
                  | !is.na(NO_UNDEVELOPED) | !is.na(EMBLEN_1) ))
  Look.at.repro=Look.at.repro[,match(c("SHEET_NO","Month","year","SPECIES","COMMON_NAME","SCIENTIFIC_NAME",
                            "TL","FL","SEX",Reprod.vars),names(Look.at.repro)),]
  
  #remove species with ID issues
  sp.ID.issues=c("Guitarfish & shovelnose ray","Wobbegong (general)","Western Wobbegong","Spotted wobbegong",
                 "Angel Shark (general)","Banded wobbegong","Blacktip sharks","Sawsharks","Cobbler Wobbegong")
  Look.at.repro=subset(Look.at.repro,!COMMON_NAME%in%sp.ID.issues)
  #Look.at.repro[Look.at.repro=="NA"]<-NA
  these.sp.rep=unique(Look.at.repro$SPECIES)
  n.sp=length(these.sp.rep)
  n.vars=length(Reprod.vars)
  
  if(User=="Matias") setwd("C:/Matias/Analyses/Surveys/Naturaliste_longline/outputs/Reproduction")
  if(User=="Dani") setwd("whatever")
  
  #see available data by species
  for(s in 1:n.sp)
  {
    d=subset(Look.at.repro,SPECIES==these.sp.rep[s])
    SPnme=unique(d$COMMON_NAME)
    fn.fig(paste("Prelim/Rep.data_",SPnme,sep=""),2400, 2400)    
    par(mfcol=c(4,3),mar=c(1,1,1,1),oma=c(.75,.75,1,.1),las=1,mgp=c(1,.6,0))
    for(v in 1:n.vars)
    {
      ID.v=match(Reprod.vars[v],names(d))
      if(is.numeric(d[,ID.v]))
      {
        if(sum(d[,ID.v],na.rm=T)>0)hist(d[,ID.v],ylab="",xlab="FL",main=Reprod.vars[v])else
          {
            plot(0,ann=F,col='transparent',axes=F)  
            legend("center","No data",bty='n',cex=1.5)
            legend("top",Reprod.vars[v],cex=1.25)
          }
      }else
      {
        if(Reprod.vars[v]=="GON_STAGE")
        {
          TAB=as.matrix(table(d[,ID.v],d$SEX))
          if(nrow(TAB)>0)barplot((TAB),beside=F,legend = rownames(TAB),main="") else
          {
            plot(0,ann=F,col='transparent',axes=F)
            legend("center","No data",bty='n',cex=1.5)
          }
        }else
        {
          TAB=as.matrix(table(d[,ID.v]))
          if(nrow(TAB)>0)barplot((TAB),beside=F,legend = rownames(TAB),main="") else
          {
            plot(0,ann=F,col='transparent',axes=F)
            legend("center","No data",bty='n',cex=1.5)
          }
        }
        legend("top",Reprod.vars[v],cex=1.25)
      }
      
    }
    mtext(SPnme,3,-.45,outer=T,cex=1.25)
    dev.off()
  }
}


#2. REPRODUCTION     
#MISSING: EMBLEN_ there's several so look at variability and report! (use mean..?)
#         TL-FL conversion, review outputs
#note: in MS say that 4 main commercial and species with ID issues were not included
if(Do.rep=="YES")
{
    #Analyses
    #1. Umbilical scars    #data only reliable for a few species
    bin=10  #size bins
    ID.v=match("UMBIL_SCAR",names(Look.at.repro))
    Select.sp=with(subset(Look.at.repro,!UMBIL_SCAR=="N"),as.matrix(table(SPECIES,UMBIL_SCAR)))
    Select.sp=Select.sp[rowSums(Select.sp)>1,]
    UMB_sca_sp=rownames(Select.sp)
    UMB_sca_sp=UMB_sca_sp[-match(c("HZ"),UMB_sca_sp)]
    fn.fig('Fig_Umb_scar',1200, 2400)
    par(mfcol=c(3,1),mar=c(1,4,3,1.5),oma=c(3,.1,.1,.1),las=1,mgp=c(.9,.6,0))
    for(s in 1:length(UMB_sca_sp))
    {
      d=subset(Look.at.repro,SPECIES==UMB_sca_sp[s]&!UMBIL_SCAR=="N")
      SPnme=unique(d$COMMON_NAME)
      d$FL=floor(d$FL/bin)*bin
      TAB=as.matrix(table(d[,ID.v],d$FL))
      if(ncol(TAB)==1)
      {
        TAB1=matrix(nrow=nrow(TAB),ncol=1)
        colnames(TAB1)=as.numeric(colnames(TAB))+bin
        TAB=cbind(TAB,TAB1)
      }
      
      MAX=as.numeric(colnames(TAB)[ncol(TAB)])
      MIN=as.numeric(colnames(TAB)[1])
      SEQ=seq(MIN,MAX,bin)
      Add1=SEQ[which(!SEQ%in%as.numeric(colnames(TAB)))]
      Add=matrix(ncol=length(Add1),nrow=nrow(TAB))
      colnames(Add)=Add1
      TAB=cbind(TAB,Add)
      TAB=TAB[,match(SEQ,colnames(TAB))]
      barplot(TAB,legend = rownames(TAB),cex.axis=1.25,cex.names=1.25,args.legend=list(bty='n',cex=1.25))
      mtext(paste(SPnme," (n=",nrow(d),")",sep=""),3,.5,cex=1.25)
      box()
    }
    mtext("Fork length (cm)",1,1.5,outer=T,cex=1.5)
    mtext("Frequency",2,-1.8,outer=T,las=3,cex=1.5)
    dev.off()
    
    
    #MALE ANALYSES
    
    #1. Clasper length, calcification and maturity
    #note: logistic ogive of the form Prop Mature=pmax/(1+exp((dat-a)/b)); a=inflex; b=slope
    
    Look.at.repro$CLASP_CALC=with(Look.at.repro,ifelse(SPECIES=="GN"&FL<140,"N",CLASP_CALC))
    Select.sp=with(subset(Look.at.repro,!is.na(CLASPLENTH) & SEX=="M"),table(SPECIES))
    Select.sp=Select.sp[Select.sp>=n.min]
    CLas_cal_sp=rownames(Select.sp)
    #remove some species for which there's not enough data
    CLas_cal_sp=CLas_cal_sp[-match(c("PJ","TN","ZE","HZ","PA","NS","TS"),CLas_cal_sp)]
    
    fn.fig('Fig_male_clas_length',1400, 2400)    
    par(mfcol=c(6,3),mar=c(1,1,2,.9),oma=c(3,3,.1,.1),las=1,mgp=c(.9,.6,0))
    for(s in 1:length(CLas_cal_sp))
    {
      d=subset(Look.at.repro,SPECIES==CLas_cal_sp[s] & CLASPLENTH>0)
      SPnme=unique(d$COMMON_NAME)
      d$CLASP_CALC=with(d,ifelse(is.na(CLASP_CALC),"NA",CLASP_CALC))
      d$CLASP_CALC_col=with(d,
                            ifelse(CLASP_CALC=="N","white",
                                   ifelse(CLASP_CALC=="P","grey60",
                                          ifelse(CLASP_CALC=="Y","black",
                                                 "grey90"))))
      d$CLASP_CALC_pt=with(d,ifelse(CLASP_CALC%in%c("N","P","Y"),21,3))
      
      
      plot(d$FL,d$CLASPLENTH,pch=d$CLASP_CALC_pt,bg=d$CLASP_CALC_col,cex=1.25,
           ylab="",xlab="")
      mtext(paste(SPnme," (n=",nrow(d),")",sep=""),3,0,cex=.72)
    }
    plot(1,ann=F,xaxt='n',yaxt='n',col="transparent")
    box(col="white")
    legend('top',c("Y","P","N","Unknown"),pch=c(21,21,21,3),pt.bg=c("black","grey60","white","black"),
           bty='n',cex=1.25,title="Calcification")
    mtext("Fork length (cm)",1,1.5,outer=T,cex=1.5)
    mtext("Clasper length (cm)",2,1.2,outer=T,las=3,cex=1.5)
    dev.off()
    
    
    #2. get L50
    dat.size.brth=data.frame(SPECIES=CLas_cal_sp,Birth=c(60,50,100,50,65,45,50,
                                                         60,35,60,60,35,52,50,50,130))
    
    fn.logis=function(dat,pmax,inflx,slop) pmax/(1+exp((dat-inflx)/slop))
    fn.logis_L50=function(dat,pmax,L50,L95) pmax/(1+exp(-log(19)*((dat-L50)/(L95-L50))))
    objfun=function(theta)
    {
      #d$Prob.n = fn.logis(d$FL,1,theta[1],theta[2])
      d$Prob.n = fn.logis_L50(d$FL,1,theta[1],theta[2])
      d$Like.n <-  with(d,ifelse(N==1,Prob.n,(1-Prob.n)))
      d$log.like=log(d$Like.n+1e-100)
      ll=sum(-d$log.like)
      return(ll)
    }
    show.classification="NO"
    #Fitting.approach="optim"
    Fitting.approach="nls"
    
    Store=vector('list',length(CLas_cal_sp))  
    names(Store)=CLas_cal_sp
    fn.fig('Fig_male_maturity_fit',1800, 2400)    
    par(mfcol=c(6,3),mar=c(1,1,2,.9),oma=c(3,3,.1,.1),las=1,mgp=c(.9,.6,0))
    for(s in 1:length(CLas_cal_sp))
    {
      d=subset(Look.at.repro,SPECIES==CLas_cal_sp[s] & CLASPLENTH>0)
      d$FL=with(d,ifelse(is.na(FL),.85*TL,FL))
      d=subset(d,!is.na(FL))
      d$N=d$CLASPLENTH/max(d$CLASPLENTH)
      
      #consider calcification
      d$N=with(d,ifelse(CLASP_CALC=="N",0,N))
      
      #add dummies to anchor model
      BirtH=dat.size.brth[s,]$Birth
      Extra.n=5
      d.optim=d[1:Extra.n,]
      d.optim[,]=NA
      d.optim$FL=seq(BirtH,BirtH*1.1,length.out=Extra.n)
      d.optim$N=rep(0,Extra.n)
      
      d=rbind(d,d.optim)
      d=subset(d,!is.na(N))
      
      #fit model
      #nls approach
      if(Fitting.approach=='nls')
      {
        mod <- nls(N~1/(1+exp((FL-inflex)/slope)), data=d,
                   start=list(inflex=quantile(d$FL,probs=.6), slope=-2),
                   nls.control(maxiter = 100, tol = 1e-06, minFactor = 1/1024,
                               printEval = FALSE, warnOnly = T))
        #L50 and L95 approach
        #note: L50 == inflex point
        # mod1 <- nls(N~1/(1+exp(-log(19)*((FL-L50)/(L95-L50)))), data=d, 
        #                 start=list(L50=quantile(d$FL,probs=.6), L95=quantile(d$FL,probs=.8)),
        #                 nls.control(maxiter = 100, tol = 1e-06, minFactor = 1/1024,
        #                             printEval = FALSE, warnOnly = T))
        
        NM=unique(d$COMMON_NAME)[1]
        newD=data.frame(FL=seq(min(d$FL),max(d$FL)))
        PRED=predict(mod,newdata=newD)
        plot(d$FL,d$N,pch=19,ylab="",xlab="",main="")
        mtext(paste(NM," (n=",nrow(d),")",sep=""),3,0,cex=.85)
        lines(newD$FL,PRED,lwd=2,col="grey50")
        
        TABL=as.data.frame(summary(mod)$coefficients[,c(1:2,4)])
        TABL$SPECIES=NM
        TABL=cbind(TABL[1,c(4,1:3)],TABL[2,1:3])
        names(TABL)[2:7]=c("inflex","inflex.SE","P","slope","slope.SE","P")
        TABL$"a (SE)"=paste(round(TABL$inflex,2)," (",round(TABL$inflex.SE,2),")",sep="")
        TABL$"b (SE)"=paste(round(TABL$slope,2)," (",round(TABL$slope.SE,2),")",sep="")
        Store[[s]]=TABL[,c(1,8,4,9,7)]
        
      }
      
      #optim approach
      if(Fitting.approach=='optim')
      {
        #theta=c(inflx=quantile(d$FL,probs=.6),slop=-5)
        theta=c(L50=quantile(d$FL,probs=.6),L95=quantile(d$FL,probs=.8))
        #fit=optim(theta,objfun, hessian=T, control=c(trace=1, maxit=10000))
        fit=optim(theta,objfun, hessian=T, method="L-BFGS-B",
                  lower=c(quantile(d$FL,probs=.2),quantile(d$FL,probs=.4)),
                  upper=c(quantile(d$FL,probs=.9),quantile(d$FL,probs=.9)),control=c(maxit=10000))
        # fit=optim(theta,objfun, hessian=T, method="L-BFGS-B",lower=c(quantile(d$FL,probs=.2),-50),
        #           upper=c(quantile(d$FL,probs=.9),-1),control=c(maxit=10000))
        # 
        NM=unique(d$COMMON_NAME)[1]
        newD=seq(min(d$FL),max(d$FL))
        #PRED=fn.logis(newD,1,fit$par[1],fit$par[2])
        PRED=fn.logis_L50(newD,1,fit$par[1],fit$par[2])
        
        plot(d$FL,d$N,pch=19,ylab="",xlab="",main="")
        mtext(paste(SPnme," (n=",nrow(d),")",sep=""),3,0,cex=.85)
        lines(newD,PRED,lwd=2,col="grey50")
        v_ob=solve(fit$hessian)	#variance covariance matrix
        std_ob=sqrt(diag(v_ob))
        TABL=data.frame(SPECIES=NM,L50=fit$par[1],L50_SE=std_ob[1],L95=fit$par[2],L95_SE=std_ob[2])
        Store[[s]]=TABL
      }
      
      
      if(show.classification=="YES")
      {   classify_data = classify_mature(d, varNames = c("FL", "CLASPLENTH"), 
                                          varSex = "SEX", selectSex = "M", method = "ld")
      par(mfrow = c(2,2))
      plot(classify_data)
      plot(classify_data, xlab = "FL", ylab = "CLASPLENTH")
      plot(classify_data, xlab = "FL", ylab = "CLASPLENTH",    col = c(2, 3), pch = c(5, 6))
      
      plot(classify_data, xlab = "FL", ylab = "CLASPLENTH", col = c(2, 3), pch = c(5, 6), lty_lines = c(1, 2), lwd_lines = c(1, 3), 
           cex = c(1, 3), main = "Classification")
      my_ogive_bayes = morph_mature(classify_data, method = "bayes", niter = 1000)
      print(my_ogive_bayes)
      par(mfrow = c(2,2))
      plot(my_ogive_bayes, xlab = "Carapace width (mm.)", ylab = "Proportion mature", col = c("blue", "red"))
      }
    }
    mtext("Fork length (cm)",1,1.5,outer=T,cex=1.5)
    mtext("Proportion mature",2,1.2,outer=T,las=3,cex=1.5)
    dev.off()
    
    Store=do.call(rbind,Store)
    write.csv(Store,"L50_pars_male.csv",row.names=F)
    
    #3. Running sperm for those with calcified claspers
    Select.sp=with(subset(Look.at.repro,!is.na(RUN_SPERM) & SEX=="M"),table(SPECIES))
    Select.sp=Select.sp[Select.sp>=n.min]
    Sperm_sp=rownames(Select.sp)
    
    fn.fig('Fig_male_running_sperm',2400, 2400)    
    par(mfcol=c(4,4),mar=c(1,1,2,1.2),oma=c(3,3,.1,.1),las=1,mgp=c(.9,.6,0))
    for(s in 1:length(Sperm_sp))
    {
      d=subset(Look.at.repro,SPECIES==Sperm_sp[s] & !is.na(RUN_SPERM) & SEX=="M" & CLASP_CALC=="Y")
      d$FL=with(d,ifelse(is.na(FL),.85*TL,FL))
      d=subset(d,!is.na(FL))
      SPnme=unique(d$COMMON_NAME)[1]
      d$COL=with(d,ifelse(RUN_SPERM=="Y","black",ifelse(RUN_SPERM=="N","grey80",NA)))
      plot(d$Month,d$FL,pch=21,col='transparent',cex=1.25,ylab="",xlab="",xlim=c(0.9,13))
      d1=subset(d,RUN_SPERM=="Y")
      points(d1$Month,d1$FL,pch=21,bg=d1$COL,cex=1.15)
      d2=subset(d,RUN_SPERM=="N")
      points(d2$Month+1.2,d2$FL,pch=21,bg=d2$COL,cex=1.15)
      NN=nrow(d)
      mtext(paste(SPnme," (n=",NN,")",sep=""),3,0,cex=.85)
    }
    plot(1,ann=F,xaxt='n',yaxt='n',col="transparent")
    box(col="white")
    legend("center",c("Y","N"),pch=21,pt.bg=c("black","grey80"),bty='n',title='Sperm presence',cex=1.5)
    
    mtext("Month",1,1.5,outer=T,cex=1.5)
    mtext("Fork length (cm)",2,1.2,outer=T,las=3,cex=1.5)
    dev.off()
    
    #FEMALE ANALYSES   
    Fem.dat=subset(DATA.bio,SPECIES%in%Reprod.sp & SEX=="F")
    Fem.dat=subset(Fem.dat,!SPECIES%in%c("BT","PJ","SD","SH","SW","WB","WD","WS","WW"))
    
    #Embryo number 
    #note: NO_UNDEVELOPED refers to embryos that did not develop. Only 5 observations, don't report
    
    #if emblen_1 is na, see if data for other emblens..
    Fem.dat$EMBLEN_1=with(Fem.dat,ifelse(is.na(EMBLEN_1) & !is.na(EMBLEN_2),EMBLEN_2,EMBLEN_1))
    
    
    #fit linear model of the form num. emb = b*FL + a
    Select.sp=with(subset(Fem.dat,!is.na(NO_EMBRYOS) & SEX=="F"),table(SPECIES))
    Select.sp=Select.sp[Select.sp>=n.min]
    Embryo_sp=rownames(Select.sp)
    
    Store=vector('list',length(Embryo_sp))  
    names(Store)=Embryo_sp
    
    fn.fig('Fig_female_embryo_number',2400, 2400)    
    par(mfcol=c(4,3),mar=c(1,1,2,1.2),oma=c(3,3,.1,.1),las=1,mgp=c(.9,.6,0))
    for(s in 1:length(Embryo_sp))
    {
      d=subset(Fem.dat,SPECIES==Embryo_sp[s] & !is.na(NO_EMBRYOS))
      d$FL=with(d,ifelse(is.na(FL),.85*TL,FL))
      d=subset(d,!is.na(FL))
      SPnme=unique(d$COMMON_NAME)[1]
      plot(d$FL,d$NO_EMBRYOS,pch=19,cex=1.25,ylab="",xlab="",ylim=c(0,max(d$NO_EMBRYOS)))
      NN=nrow(d)
      mtext(paste(SPnme," (n=",NN,")",sep=""),3,0,cex=.85)
      mod=lm(NO_EMBRYOS~FL,d)
      PRED.d=data.frame(FL=seq(min(d$FL),max(d$FL)))
      PRED=predict(mod,newdata=PRED.d)
      lines(PRED.d$FL,PRED,lwd=2,col="grey50")
      
      TABL=as.data.frame(summary(mod)$coefficients[,c(1:2,4)])
      TABL$SPECIES=SPnme
      TABL=cbind(TABL[1,c(4,1:3)],TABL[2,1:3])
      names(TABL)[2:7]=c("inflex","inflex.SE","P","slope","slope.SE","P")
      TABL$"a (SE)"=paste(round(TABL$inflex,2)," (",round(TABL$inflex.SE,2),")",sep="")
      TABL$"b (SE)"=paste(round(TABL$slope,2)," (",round(TABL$slope.SE,2),")",sep="")
      Store[[s]]=TABL[,c(1,8,4,9,7)]
    }
    mtext("Fork length (cm)",1,1.5,outer=T,cex=1.5)
    mtext("Number of embryos",2,1.2,outer=T,las=3,cex=1.5)
    dev.off()
    Store=do.call(rbind,Store)
    write.csv(Store,"Female_pups_regression.csv",row.names=F) 
    
    #Embryo length vs time
    Fem.dat$EMBLEN_1=with(Fem.dat,ifelse(UTERINESTG==4,0,EMBLEN_1))
    Fem.dat$EMBLEN_1=with(Fem.dat,ifelse(!UTERINESTG==4 & EMBLEN_1==0,NA,EMBLEN_1))
    
    Select.sp=with(subset(Fem.dat,!is.na(EMBLEN_1) & SEX=="F"),table(SPECIES))
    Select.sp=Select.sp[Select.sp>=n.min]
    Embryo_sp=rownames(Select.sp)
    fn.fig('Fig_female_embryo_length',2400, 2400)    
    par(mfcol=c(4,4),mar=c(1,1,2,1.2),oma=c(3,3,.1,.1),las=1,mgp=c(.9,.6,0))
    for(s in 1:length(Embryo_sp))
    {
      d=subset(Fem.dat,SPECIES==Embryo_sp[s]&!is.na(EMBLEN_1))
      SPnme=unique(d$COMMON_NAME)[1]
      d$COL=with(d,ifelse(UTERINESTG==4,"grey65",ifelse(UTERINESTG==5,"black",NA)))
      plot(d$Month,d$EMBLEN_1,pch=21,bg=d$COL,cex=1.25,ylab="",xlab="",xlim=c(0.9,13),
           ylim=c(0,max(d$EMBLEN_1)))
      mtext(paste(SPnme," (n=",nrow(d),")",sep=""),3,0,cex=.9)
    }
    mtext("Month",1,1.5,outer=T,cex=1.5)
    mtext("Embryo length (cm)",2,1.2,outer=T,las=3,cex=1.5)
    plot(1,ann=F,xaxt='n',yaxt='n',col="transparent")
    box(col="white")
    legend('bottomleft',c("4","5"),pch=21,pt.bg=c("grey65","black"),bty='n',cex=1.5,
           title="Uterus condition")
    dev.off()
    
    
    #Embryo length vs Max ova diameter
    Select.sp=with(subset(Fem.dat,SPECIES%in%Embryo_sp&!is.na(EMBLEN_1) & !is.na(MAXOVRYDIA) & SEX=="F"),table(SPECIES))
    Select.sp=Select.sp[Select.sp>=n.min]
    Embryo_sp_ova=rownames(Select.sp)
    
    fn.fig('Fig_female_embryo_length_ova_diam',1200, 2400)    
    par(mfcol=c(3,1),mar=c(1,1,2,1.2),oma=c(3,3,.1,.1),las=1,mgp=c(.9,.6,0))
    for(s in 1:length(Embryo_sp_ova))
    {
      d=subset(Fem.dat,SPECIES==Embryo_sp_ova[s]&!is.na(EMBLEN_1) &!is.na(MAXOVRYDIA))
      SPnme=unique(d$COMMON_NAME)[1]
      d$COL=with(d,ifelse(UTERINESTG==4,"grey65",ifelse(UTERINESTG==5,"black",NA)))
      plot(d$EMBLEN_1,d$MAXOVRYDIA,pch=21,bg=d$COL,cex=1.5,ylab="",xlab="",cex.axis=1.25,
           ylim=c(0,max(d$MAXOVRYDIA)))
      mtext(paste(SPnme," (n=",nrow(d),")",sep=""),3,0,cex=1.1)
    }
    mtext("Embryo length (cm)",1,1.5,outer=T,cex=1.5)
    mtext("Maximum ova diameter (cm)",2,1.2,outer=T,las=3,cex=1.5)
    #plot(1,ann=F,xaxt='n',yaxt='n',col="transparent")
    #box(col="white")
    legend('topright',c("4","5"),pch=21,pt.bg=c("grey65","black"),bty='n',cex=1.5,title="Uterus condition")
    dev.off()
    
    
    #Embryo length -Ova diameter vs time    
    fn.fig('Fig_female_emblen_ova_diam_time',1400, 2400)    
    par(mfcol=c(3,1),mar=c(1,1,2,4),oma=c(3,3,.1,1),las=1,mgp=c(.9,.6,0))
    for(s in 1:length(Ova_sp))
    {
      d=subset(Fem.dat,SPECIES==Ova_sp[s] & SEX=="F")
      d$MAXOVRYDIA=d$MAXOVRYDIA/10   #convert to cm as MOD measured in mm
      SPnme=unique(d$COMMON_NAME)[1]
      with(subset(d,!is.na(MAXOVRYDIA)),plot(Month,MAXOVRYDIA,pch=21,bg="grey30",cex=1.75,ylab="",xlab="",xlim=c(0.9,13),
                                             ylim=c(0,max(MAXOVRYDIA)),cex.axis=1.25))
      par(new = T)
      with(subset(d,!is.na(EMBLEN_1)), plot(Month, EMBLEN_1, pch=23, axes=F, xlab=NA, ylab=NA,bg="grey90",cex=1.4,
                                            ylim=c(0,max(EMBLEN_1)),cex.axis=1.25))
      axis(side = 4,cex.axis=1.25)
      mtext(paste(SPnme," (n=",nrow(subset(d,!is.na(MAXOVRYDIA))),")",sep=""),3,0,cex=1)
    }
    legend('right',c("ova diam","emb len"),pch=c(21,23),pt.bg=c("grey30","grey90"),bty='n',cex=1.75)
    
    mtext("Month",1,1.5,outer=T,cex=1.5)
    mtext("Maximum ova diameter (cm)",2,1.2,outer=T,las=3,cex=1.5)
    mtext("Embryo length (cm)",side = 4, -.5,outer=T,las=3,cex=1.5)
    dev.off()
    
    
    
    #Embryo Ova diameter vs time
    Select.sp=with(subset(Fem.dat,!is.na(MAXOVRYDIA) & SEX=="F"),table(SPECIES))
    Select.sp=Select.sp[Select.sp>=n.min]
    Ova_sp=rownames(Select.sp)
    fn.fig('Fig_female_ova_diam',2400, 2400)    
    par(mfcol=c(4,2),mar=c(1,1,2,1.2),oma=c(3,3,.1,.1),las=1,mgp=c(.9,.6,0))
    for(s in 1:length(Ova_sp))
    {
      d=subset(Fem.dat,SPECIES==Ova_sp[s] & !is.na(MAXOVRYDIA) & SEX=="F")
      SPnme=unique(d$COMMON_NAME)[1]
      d$COL=with(d,ifelse(UTERINESTG==4,"grey65",ifelse(UTERINESTG==5,"black",
                                                        ifelse(UTERINESTG==1,"white",ifelse(UTERINESTG==2,"blue",
                                                                                            ifelse(UTERINESTG==3,"green",ifelse(UTERINESTG==6,"red",NA)))))))                                   
      
      
      plot(d$Month,d$MAXOVRYDIA,pch=21,bg=d$COL,cex=1.25,ylab="",xlab="",xlim=c(0.9,13),
           ylim=c(0,max(d$MAXOVRYDIA)))
      mtext(paste(SPnme," (n=",nrow(d),")",sep=""),3,0,cex=.85)
    }
    mtext("Month",1,1.5,outer=T,cex=1.5)
    mtext("Maximum ova diameter (cm)",2,1.2,outer=T,las=3,cex=1.5)
    plot(1,ann=F,xaxt='n',yaxt='n',col="transparent")
    box(col="white")
    legend('top',as.character(1:6),pch=21,pt.bg=c("white","blue","green","grey65","black","red")
           ,bty='n',cex=1.25,title="Uterus condition")
    dev.off()
    
    
    #Maturity 
    Select.sp=with(subset(Fem.dat,!is.na(UTERINESTG) & SEX=="F"),table(SPECIES))
    Select.sp=Select.sp[Select.sp>=n.min]
    MAT_sp=rownames(Select.sp)
    MAT_sp=MAT_sp[-match(c("CW","GU","MS","WC","TA"),MAT_sp)]
    Select.sp=with(subset(Fem.dat,!is.na(GON_STAGE) & SEX=="F"),table(SPECIES))
    id=which(!rownames(Select.sp)%in%MAT_sp)    #no species had GI info but no UI info
    
    #2. get L50
    dat.size.brth=rbind(dat.size.brth,
                        data.frame(SPECIES=c("HZ","NS","TA"),
                                   Birth=c(50,35,22)))
    Fem.dat=merge(Fem.dat,dat.size.brth,by="SPECIES",all.x=T)
    
    Store=vector('list',length(MAT_sp))  
    names(Store)=MAT_sp
    fn.fig('Fig_female_maturity_fit',1800, 2400)    
    par(mfcol=c(5,3),mar=c(1,1,2,.9),oma=c(3,3,.1,.1),las=1,mgp=c(.9,.6,0))
    for(s in 1:length(MAT_sp))
    {
      d=subset(Fem.dat,SPECIES==MAT_sp[s] & SEX=="F")
      d$FL=with(d,ifelse(is.na(FL),.85*TL,FL))
      d=subset(d,!is.na(FL))
      d$N=with(d,ifelse(UTERINESTG%in%3:6 ,1,ifelse(UTERINESTG%in%1:2,0,NA)))
      d$N=with(d,ifelse(is.na(N) & GON_STAGE==3,1,N))
      d$N=with(d,ifelse(is.na(N) & GON_STAGE%in%c(1,2),0,N))
      
      #add dummies to anchor model
      BirtH=dat.size.brth[s,]$Birth
      Extra.n=5
      d.optim=d[1:Extra.n,]
      d.optim[,]=NA
      d.optim$FL=seq(BirtH,BirtH*1.1,length.out=Extra.n)
      d.optim$N=rep(0,Extra.n)
      
      d=rbind(d,d.optim)
      d=subset(d,!is.na(N))
      
      if(MAT_sp[s]=="TG") d=subset(d,!(N==1 & FL<200))
      NM=unique(d$COMMON_NAME)[1]
      #fit model
      #binominal glm
      # mod=glm(N~FL,d,family=binomial())
      # a=predict(mod,newdata=data.frame(FL=seq(min(d$FL),max(d$FL))),type='response')
      # lines(seq(min(d$FL),max(d$FL)),a)
      
      #nls approach
      if(Fitting.approach=='nls')
      {
        mod <- nls(N~1/(1+exp((FL-inflex)/slope)), data=d,
                   start=list(inflex=quantile(d$FL,probs=.6), slope=-2),
                   nls.control(maxiter = 100, tol = 1e-06, minFactor = 1/1024,
                               printEval = FALSE, warnOnly = T))
        #L50 and L95 approach
        #note: L50 == inflex point
        # mod1 <- nls(N~1/(1+exp(-log(19)*((FL-L50)/(L95-L50)))), data=d, 
        #                 start=list(L50=quantile(d$FL,probs=.6), L95=quantile(d$FL,probs=.8)),
        #                 nls.control(maxiter = 100, tol = 1e-06, minFactor = 1/1024,
        #                             printEval = FALSE, warnOnly = T))
        
        
        newD=data.frame(FL=seq(min(d$FL),max(d$FL)))
        PRED=predict(mod,newdata=newD)
        # d$COL=with(d,ifelse(UTERINESTG==4,"grey50",ifelse(UTERINESTG==5,"grey70",
        #             ifelse(UTERINESTG==3,"grey30",ifelse(UTERINESTG==6,"grey85","white")))))                                   
        d$COL=1
        plot(d$FL,d$N,pch=21,ylab="",xlab="",main="",bg=d$COL,ylim=c(0,1))
        mtext(paste(NM," (n=",nrow(d),")",sep=""),3,0,cex=.85)
        lines(newD$FL,PRED,lwd=2,col="grey50")
        
        TABL=as.data.frame(summary(mod)$coefficients[,c(1:2,4)])
        TABL$SPECIES=NM
        TABL=cbind(TABL[1,c(4,1:3)],TABL[2,1:3])
        names(TABL)[2:7]=c("inflex","inflex.SE","P","slope","slope.SE","P")
        TABL$"a (SE)"=paste(round(TABL$inflex,2)," (",round(TABL$inflex.SE,2),")",sep="")
        TABL$"b (SE)"=paste(round(TABL$slope,2)," (",round(TABL$slope.SE,2),")",sep="")
        Store[[s]]=TABL[,c(1,8,4,9,7)]
        
      }
      
      #optim approach
      if(Fitting.approach=='optim')
      {
        #theta=c(inflx=quantile(d$FL,probs=.6),slop=-5)
        theta=c(L50=quantile(d$FL,probs=.6),L95=quantile(d$FL,probs=.8))
        #fit=optim(theta,objfun, hessian=T, control=c(trace=1, maxit=10000))
        fit=optim(theta,objfun, hessian=T, method="L-BFGS-B",
                  lower=c(quantile(d$FL,probs=.2),quantile(d$FL,probs=.4)),
                  upper=c(quantile(d$FL,probs=.9),quantile(d$FL,probs=.9)),control=c(maxit=10000))
        # fit=optim(theta,objfun, hessian=T, method="L-BFGS-B",lower=c(quantile(d$FL,probs=.2),-50),
        #           upper=c(quantile(d$FL,probs=.9),-1),control=c(maxit=10000))
        # 
        NM=unique(d$COMMON_NAME)[1]
        newD=seq(min(d$FL),max(d$FL))
        #PRED=fn.logis(newD,1,fit$par[1],fit$par[2])
        PRED=fn.logis_L50(newD,1,fit$par[1],fit$par[2])
        
        d$COL=with(d,ifelse(UTERINESTG==4,"grey50",ifelse(UTERINESTG==5,"grey70",
                                                          ifelse(UTERINESTG==3,"grey30",ifelse(UTERINESTG==6,"grey85","white")))))                                   
        
        plot(d$FL,d$N,pch=21,ylab="",xlab="",main="",bg=d$COL,ylim=c(0,1))
        
        mtext(paste(SPnme," (n=",nrow(d),")",sep=""),3,0,cex=.85)
        lines(newD,PRED,lwd=2,col="grey50")
        v_ob=solve(fit$hessian)	#variance covariance matrix
        std_ob=sqrt(diag(v_ob))
        TABL=data.frame(SPECIES=NM,L50=fit$par[1],L50_SE=std_ob[1],L95=fit$par[2],L95_SE=std_ob[2])
        Store[[s]]=TABL
      }
      
      
      if(show.classification=="YES")
      {   classify_data = classify_mature(d, varNames = c("FL", "CLASPLENTH"), 
                                          varSex = "SEX", selectSex = "M", method = "ld")
      par(mfrow = c(2,2))
      plot(classify_data)
      plot(classify_data, xlab = "FL", ylab = "CLASPLENTH")
      plot(classify_data, xlab = "FL", ylab = "CLASPLENTH",    col = c(2, 3), pch = c(5, 6))
      
      plot(classify_data, xlab = "FL", ylab = "CLASPLENTH", col = c(2, 3), pch = c(5, 6), lty_lines = c(1, 2), lwd_lines = c(1, 3), 
           cex = c(1, 3), main = "Classification")
      my_ogive_bayes = morph_mature(classify_data, method = "bayes", niter = 1000)
      print(my_ogive_bayes)
      par(mfrow = c(2,2))
      plot(my_ogive_bayes, xlab = "Carapace width (mm.)", ylab = "Proportion mature", col = c("blue", "red"))
      }
    }
    
    # plot(1,ann=F,xaxt='n',yaxt='n',col="transparent")
    # box(col="white")
    # legend('top',c("Immature","UI 3","UI 4","UI 5","UI 6"),pch=21,
    #        pt.bg=c("white","grey30","grey50","grey70","grey85"),bty='n',cex=1.25)
    mtext("Fork length (cm)",1,1.5,outer=T,cex=1.5)
    mtext("Proportion mature",2,1.2,outer=T,las=3,cex=1.5)
    dev.off()
    
    Store=do.call(rbind,Store)
    write.csv(Store,"L50_pars_female.csv",row.names=F)
  }



#3. DIET  
#note: diet analyses based on frequency of occurrence only (weight/number prey not available/reliable
#  MISSING: if using tropics only (survey), not enough samples; separate tropical from subtropical;
#             report size distribution; map of where samples from
#             Separate some species into large, small? (e.g dusky) see Papastamatiou et al 2006
#        minimum sample size: 10
#            Also, among species comparisons should be meaningful (e.g. same location/shot)

if(Do.diet=="YES")
{
  Diet=subset(Diet,!COMMON_NAME=="Parrotfish (general)")
  Dummy.2=Dummy
  Dummy.2[Dummy.2>0]=1
  Look=sort(colSums(Dummy.2))
  DieT=Dummy.2[match("STMCH_CONT",row.names(Dummy.2)),]
  DieT=sort(DieT)
  DieT=subset(DieT,DieT>0)
  DieT=names(DieT)
  Look.at.diet=subset(DATA.bio,SPECIES%in%DieT & !is.na(STMCH_CONT),select=c(SPECIES,SHEET_NO,TL,FL,SEX,STMCH_CONT,STMCH_FULL,
                                   COMMON_NAME,SCIENTIFIC_NAME,year))
  Look.at.diet=subset(Look.at.diet,year==Current.Yr)
  Look.at.diet=Look.at.diet[,-match("year",names(Look.at.diet))]
  if(User=="Matias") setwd("C:/Matias/Analyses/Surveys/Naturaliste_longline/outputs/Diet")
  if(User=="Dani") setwd("whatever")
  write.csv(Look.at.diet,"Look.at.diet_all_sp.csv",row.names=F)
  
  if(User=="Matias") HNDL.diet="C:/Matias/Analyses/Surveys/Naturaliste_longline/outputs/Diet/"
  if(User=="Dani") HNDL.diet="C:/Matias/Analyses/Surveys/Naturaliste_longline/outputs/Diet/"
  
  Not.prey=c("SHEET_NO","COMMON_NAME","SCIENTIFIC_NAME","SPECIES","SEX","TL","FL","STMCH_FULL","STMCH_CONT")
  Prey=colnames(Diet)[-match(Not.prey,colnames(Diet))]
  id=colSums(Diet[,Prey],na.rm=T)
  id=subset(id,id==0)
  if(length(id)>0)
  {
    Diet=Diet[,-match(names(id),names(Diet))]
    Prey=colnames(Diet)[-match(Not.prey,colnames(Diet))]
  }
  
  #Calculate frequency of occurrence by prey for each shark species
  Spec.names=read.csv("Spec.names.csv",stringsAsFactors=F)
  Diet=Diet[,-match("COMMON_NAME",names(Diet))]
  Diet=merge(Diet,Spec.names,by="SPECIES",all.x=T)
  Diet$SPECIES=as.character(Diet$SPECIES)
  
  #tropical
  Tropics=-23.43697
  Rel.dat=subset(DATA.bio,SHEET_NO%in%unique(Diet$SHEET_NO),select=c(SHEET_NO,Mid.Lat, Mid.Long, Lat.round, Long.round))
  Rel.dat=Rel.dat[!duplicated(Rel.dat$SHEET_NO),]
  Tropical=subset(Rel.dat,Lat.round>=(-23.5))
  #Diet=subset(Diet,SHEET_NO%in%unique(Tropical$SHEET_NO))
  
  Diet.dat=data.table(Diet[,match(c("SPECIES",Prey),names(Diet))])
  Table.diet=Diet.dat[, lapply(.SD, sum,na.rm=T), by=SPECIES]     #aggregate all prey by species
  N.stom=table(Diet$SPECIES)
  N.stom=N.stom[match(Table.diet$SPECIES,names(N.stom))]
  N.STOM=matrix(N.stom,nrow=nrow(Table.diet[,2:ncol(Table.diet)]),ncol=ncol(Table.diet[,2:ncol(Table.diet)]))
  Table.diet[,2:ncol(Table.diet)]=100*Table.diet[,2:ncol(Table.diet)]/N.STOM
  
  Table.diet.prey.group=Diet[,-match(c("SHEET_NO","COMMON_NAME",
            "SCIENTIFIC_NAME","SEX","TL","FL","STMCH_FULL","STMCH_CONT"),names(Diet))]
  ID_Prey_TL2=match(colnames(Table.diet.prey.group)[-1],Prey_cnvrsn$Prey)
  colnames(Table.diet.prey.group)[2:ncol(Table.diet.prey.group)]=as.character(Prey_cnvrsn$Prey_TL[ID_Prey_TL2])
  Prey.groups=unique(as.character(Prey_cnvrsn$Prey_TL))
  Agg.diet=as.data.frame(matrix(nrow=nrow(Table.diet.prey.group),ncol=(1+length(Prey.groups))))
  colnames(Agg.diet)=c("SPECIES",Prey.groups)
  Agg.diet$SPECIES=Table.diet.prey.group$SPECIES
  for(x in 1:length(Prey.groups))
  {
    d=Table.diet.prey.group[,which(names(Table.diet.prey.group)==Prey.groups[x])]
    if(is.data.frame((d)))dd=rowSums(d,na.rm=T)
    if(is.integer(d))dd=d
    dd[dd>0]=1
    Agg.diet[,x+1]=dd
  }
  Agg.diet=data.table(Agg.diet)
  Agg.diet=Agg.diet[, lapply(.SD, sum,na.rm=T), by=SPECIES]
  N.STOM=matrix(N.stom,nrow=nrow(Agg.diet[,2:ncol(Agg.diet)]),ncol=ncol(Agg.diet[,2:ncol(Agg.diet)]))
  Agg.diet[,2:ncol(Agg.diet)]=100*Agg.diet[,2:ncol(Agg.diet)]/N.STOM
  
  Diet.Table1=cbind(Table.diet,N.STOM[,1])
  names(Diet.Table1)[ncol(Diet.Table1)]="Sample size"
  Diet.Table1[,2:(ncol(Diet.Table1)-1)]=round(Diet.Table1[,2:(ncol(Diet.Table1)-1)],3)
  Diet.Table1=Diet.Table1[order(-Diet.Table1$"Sample size"),]
  
  SP.nm=Diet[,match(c("COMMON_NAME","SCIENTIFIC_NAME","SPECIES"),names(Diet))]
  SP.nm=SP.nm[!duplicated(SP.nm$SPECIES),]
  Diet.Table1=merge(SP.nm,Diet.Table1,by="SPECIES")
  
  write.csv(Diet.Table1,paste(HNDL.diet,"Table1.csv",sep=""),row.names=F)
  
  
  #Pie charts   
  
    #aggregated prey
  Agg.diet=Agg.diet[order(Agg.diet$Teleost,Agg.diet$Cephalopods),]
  SpEcis=as.character(unique(Agg.diet$SPECIES))
  
  #colfunc <- colorRampPalette(c("black", "white"))
  #all.ColS.pie=colfunc(length(Prey.groups))
  all.ColS.pie=c("grey60","grey20","grey40","grey50","white",
                 "grey80","grey70","grey28","grey90","black")
  all.ColS.pie=data.frame(col=all.ColS.pie,Prey=Prey.groups)
  
  fn.pie=function(SP)
  {
    A=as.data.frame(subset(Agg.diet,SPECIES==SP))
    N=subset(Diet.Table1,SPECIES==SP)
    SPec=as.character(N$COMMON_NAME)
    N=N$'Sample size'
    A=A[,-match("SPECIES",names(A))]
    id=which(A>0)
    ColS.pie=all.ColS.pie[id,]
    A=A[,id]
    pie(unlist(A), labels = "", col = as.character(ColS.pie$col),border="transparent",main="")
    mtext(paste(SPec," (n= ",N,") ",sep=""),3,-1,cex=.7)
    
  }
  fn.fig(paste(HNDL.diet,"Pie_occurrence",sep=""),2400, 2400)  
  par(mfrow=c(6,4), mar=c(.1,.1,.1,0.5), mgp = c(1.5, 0.3, 0), tck = -0.01)
  for(i in 1:length(SpEcis))fn.pie(SP=SpEcis[i])
  plot(1:1,ann=F,col="transparent",axes=F)
  legend("center", as.character(all.ColS.pie$Prey), cex = 1, 
                     fill = as.character(all.ColS.pie$col),bty='n')
  dev.off()
  
      #high resolution prey species
  # fn.pie=function(SP)
  # {
  #   A=subset(Diet.Table1,SPECIES==SP)
  #   N=A$'Sample size'
  #   SPec=A$COMMON_NAME
  #   A1=A[,-match(c("SPECIES","COMMON_NAME","SCIENTIFIC_NAME",'Sample size'),names(A))]
  #   ColS.pie=rainbow(ncol(A1))
  #   ColS.pie=data.frame(col=ColS.pie,Prey=colnames(A1))
  #   id=which(A1>0)
  #   ColS.pie=ColS.pie[id,]
  #   A1=A1[,id]
  #   pie(unlist(A1), labels = colnames(A1), col = as.character(ColS.pie$col),border="white",
  #       main=paste(SPec," (n= ",N,") ",sep=""))
  # }
  # par(mfrow=c(6,4), mar=c(.1,.1,.1,0.5), mgp = c(1.5, 0.3, 0), tck = -0.01)
  # for(i in 1:length(SpEcis))fn.pie(SP=SpEcis[i])
  

  #Calculate trophic level
  Prey_TL=subset(Prey_TL, Prey%in%names(Agg.diet)[-1])
  Prey_TL=Prey_TL[order(Prey_TL$TL),]
  Agg.diet=as.data.frame(Agg.diet)
  Agg.diet=Agg.diet[,match(c("SPECIES",as.character(Prey_TL$Prey)),names(Agg.diet))]
  PREDATORS=as.character(Agg.diet$SPECIES)
  Troph.rel=Agg.diet[,-match(c("SPECIES"),names(Agg.diet))]
  
  PREY=Prey_TL$Prey
  fn.TL=function(d)
  {
    d=d[,-1]
    Nms=colnames(d)
    d=data.frame(t(d))
    colnames(d)="FO"
    d$Prey=Nms
    d=merge(d,Prey_TL,by="Prey",all.x=T)
    d$Prop=d$FO/sum(d$FO)    #get proportion
    TL=1+ sum(d$Prop*d$TL)
    return(TL)
  }
  Predator.TL=data.frame(SPECIES=PREDATORS,TL=NA)
  for(i in 1:length(PREDATORS)) Predator.TL$TL[i]=fn.TL(d=subset(Agg.diet,SPECIES==PREDATORS[i]))
  
  a=subset(Diet.Table1,select=c(SPECIES,COMMON_NAME))
  a=a[!duplicated(a$SPECIES),]
  Predator.TL=merge(Predator.TL,a,by="SPECIES",all.x=T)
  
    
  
  # #use these files Predator_TL Prey_TL
  # TL.pred=PREDATORS               # Predator_TL hast theoretical TL for these species 
  # names(TL.pred)=PREDATORS
  # TL.pred=runif(length(TL.pred),3.5,5.5)
  # TL.prey=PREY                    #DUMMY, change accordingly:  use this file Prey_TL
  # names(TL.prey)=PREY
  # TL.prey=runif(length(TL.prey),1,3)
  
  #Show trophic links
  function.plot.TL=function(DieT,TL.pred,TL.prey,scaler)
  {
    #Locations on X axis
    MaxX=max(dim(DieT))
    Location.pred=seq(1,MaxX,length.out=nrow(DieT))
    names(Location.pred)=PREDATORS
    Location.prey=seq(1,MaxX,length.out=ncol(DieT))
    names(Location.prey)=colnames(Agg.diet)[-1]
    
    #set things up
    plot(1,1,col="transparent",xlim=c(0,MaxX*1.1),
        ylim=c(min(TL.prey$TL)*.9,1.1*max(c(TL.pred$TL,TL.prey$TL))),
         ylab="",xlab="",xaxt="n")
    TL.pred=TL.pred[match(PREDATORS,TL.pred$SPECIES),]
    TL.prey=TL.prey[match(names(Location.prey),TL.prey$Prey),]
    
    #plot trophic links based on FO
    for(r in 1:nrow(DieT))
    {
      dd=DieT[r,]*scaler
       from.x=Location.pred[r]
      from.y=TL.pred$TL[r]
      to.x=Location.prey
      to.y=TL.prey$TL
      
      CL=dd
      CL[CL>0]="black"
      CL[CL=="0"]="transparent"
      CL=unlist(CL)
      arrows(rep(from.x,length(to.x)),rep(from.y,length(to.y)),to.x,to.y,lwd=dd,col=CL,length = 0.15)
    }
    
    #plot predator
    text(Location.pred,TL.pred$TL,PREDATORS,pos=3,offset = .2,srt=45)
    
    #plot prey
    text(Location.prey,TL.prey$TL,TL.prey$Prey,pos=1,offset = .1,srt=-15)
    
    
  }
  fn.fig(paste(HNDL.diet,"Trophic_links",sep=""),2400, 1200)    
  par(mfcol=c(1,1),mar=c(1.5,2.5,.1,.1),oma=c(.1,.2,.1,.2),las=1,mgp=c(1.9,.55,0))
  function.plot.TL(DieT=Troph.rel/100,TL.pred=Predator.TL,TL.prey=Prey_TL,scaler=3)
  mtext("Trophic level",2,1.5,las=3,cex=1.5)
  mtext("Prey",1,0.5,cex=1.5)
  dev.off()
  
  
  #ACA
  
  #Dietary overlap based on Simplified Morisita index with diet expressed as proportions
  
  
  #Generalist-Specialist spectrum / Prey diversity indices / Omnivory index  / Size of shark  /home range ( Roff et al 2016)

  
  #Map density distribution of shots to put study in context
  
  
  #size distribution of species
  
  
  #Multivariate analysis
  #note: could use FL or TL to see ontogenetic effects (e.g. good size range for duskies, size range other species??)
  
}








