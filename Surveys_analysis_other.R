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
#library(glmmADMB) #zero-inflated models
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


# DATA SECTION -------------------------------------------

if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

#Source Sharks data base 
if(User=="Matias") source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R"))
if(User=="Dani") source("C:/Users/ddw/Documents/Projects/Naturaliste_longline/Source codes/Source_Shark_bio.R")
rm(DATA)
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/fn.fig.R"))
Do.tiff="NO"    #select figure extension
Do.jpeg="YES"

#Comments data  
if(User=="Matias") 
{
  Comnt_Alive=read.csv(handl_OneDrive("Data/Naturaliste/Comments_Alive.csv"))
  Comnt_Dead=read.csv(handl_OneDrive("Data/Naturaliste/Comments_Dead.csv"))
}
if(User=="Dani") 
{
  Comnt_Alive=read.csv("C:/Users/ddw/Documents/Projects/Naturaliste_longline/Source codes/Comments_Alive.csv")
  Comnt_Dead=read.csv("C:/Users/ddw/Documents/Projects/Naturaliste_longline/Source codes/Comments_Dead.csv")  
}

Biol.vars=c("UMBIL_SCAR","CLASPLENTH","CLASP_CALC","GON_STAGE","RUN_SPERM","MAXOVRYDIA",
            "NO_YOLKOVA","UTERINESTG","NO_EMBRYOS","NO_UNDEVELOPED","EMBLEN_1",
            "STMCH_FULL","STMCH_CONT")

#Diet
#note: Diet is a csv file that must be manually updated (see '#manually add new years data' below). Data is presence/absence
if(User=="Matias")
{
  Diet=read.csv(handl_OneDrive("Data/Naturaliste/Look.at.diet_1.csv"))
  Predator_TL=read.csv(handl_OneDrive("Data/Naturaliste/Predator_TL.csv"))
  Prey_TL=read.csv(handl_OneDrive("Data/Naturaliste/Prey_TL.csv"))
  Prey_cnvrsn=read.csv(handl_OneDrive("Data/Naturaliste/Prey.csv"))
  
}
if(User=="Dani")   Diet=read.csv("C:/Users/ddw/Documents/Projects/Naturaliste_longline/Diet/Look.at.diet_1.csv")



# CONTROL SECTION -------------------------------------------------------------------------

#control if doing survival analysis
#Do.survival="YES"
Do.survival="NO"

#control if doing reproduction and diet analysis
Do.rep="NO"
Do.diet="YES"

#control if doing ecosystem analysis
Do.ecosystem="NO"



# DATA MANIPULATION -------------------------------------------------------------------------

#Create data sets for PCS and biological analysis 
#note: we assume that any shark that was dissected was not sacrificed so it was dissected 
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


if(Do.rep=="YES"| Do.diet=="YES")
{
  #Exclude all bony fish from data
  n.min=5
  Exclude <- c("CD", "BN", "FR", "GB", "GG", "SB","SE", "ST", "WR") #Ostheichthyes 
  
  DATA.rep <-   filter(DATA.bio,year>=1990)%>%
    # dplyr::select(SHEET_NO,date, Day, Month, year, Lat.round, Long.round, SPECIES, SEX,TL, FL, COMMON_NAME, SCIENTIFIC_NAME, UMBIL_SCAR, 
    #               CLASPLENTH, CLASP_CALC, RUN_SPERM,
    #               UTERINESTG, NO_EMBRYOS,EMBLEN_1, NO_UNDEVELOPED,STMCH_FULL,STMCH_CONT, GON_STAGE, NO_YOLKOVA, MAXOVRYDIA,Biol.dat)%>%
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
                                       "TL","FL","Disc.width","SEX",Reprod.vars),names(Look.at.repro))]
  
  Fem.dat=subset(DATA.bio,!SPECIES%in%c("GM","TK","BW","WH") & SEX=="F")
  
  #Look.at.repro[Look.at.repro=="NA"]<-NA
  these.sp.rep=unique(Look.at.repro$SPECIES)
  n.sp=length(these.sp.rep)
  n.vars=length(Reprod.vars)
  
  if(User=="Matias") setwd(handl_OneDrive("Analyses/Surveys/Naturaliste_longline/outputs/Reproduction"))
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

#1. IMMEDIATE POST CAPTURE MORTALITY -------------------------------------------------------------------------
if(Do.survival=="YES")  
{
  #Combine blacktips
  DATA.PCS$SPECIES=with(DATA.PCS,ifelse(SPECIES=="BL","BT",SPECIES))
  DATA.PCS$COMMON_NAME=with(DATA.PCS,ifelse(SPECIES=="Common blacktip shark","Blacktip sharks",COMMON_NAME))
  DATA.PCS$SCIENTIFIC_NAME=with(DATA.PCS,ifelse(SPECIES=="Carcharhinus limbatus","Carcharhinus limbatus/tilstoni",SCIENTIFIC_NAME))
  
  N.min=5  #minimum number of records per species for PCS analysis
  
  if(User=="Matias") setwd(handl_OneDrive("Analyses/Surveys/Naturaliste_longline/outputs/Survival"))
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
    mtext("Longitude (?E)",1,-.8,outer=T,cex=2)
    mtext("Latitude (?S)",2,line=-0.6,outer=T,cex=2)
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

#2. Export latest diet to add to  -------------------------------------------------------------------------
THESE.YEARS=2018:2024  #years for manual diet update (add to 'Look.at.diet_1.csv')
Dummy.2=Dummy
Dummy.2[Dummy.2>0]=1
Look=sort(colSums(Dummy.2))
DieT=Dummy.2[match("STMCH_CONT",row.names(Dummy.2)),]
DieT=sort(DieT)
DieT=subset(DieT,DieT>0)
DieT=names(DieT)

Look.at.diet=subset(DATA.bio,SPECIES%in%DieT & !is.na(STMCH_CONT),select=c(SPECIES,SHEET_NO,TL,FL,SEX,STMCH_CONT,
                                                                           STMCH_FULL,COMMON_NAME,SCIENTIFIC_NAME,year))
#manually add new years data
Look.at.diet=subset(Look.at.diet,year%in%THESE.YEARS)
Look.at.diet=Look.at.diet[,-match("year",names(Look.at.diet))]
if(User=="Matias") setwd(handl_OneDrive("Analyses/Surveys/Naturaliste_longline/outputs/Diet"))
if(User=="Georgina") setwd("whatever")
Look.at.diet=Look.at.diet%>%
  relocate(SHEET_NO,COMMON_NAME,SCIENTIFIC_NAME,SPECIES,SEX,TL,FL,STMCH_FULL,STMCH_CONT)
dumi=Diet%>%
  dplyr::select(-names(Look.at.diet))
dumi=dumi[1:nrow(Look.at.diet),]
dumi[]=''
write.csv(cbind(Look.at.diet,dumi),"Look.at.diet_all_sp.csv",row.names=F)

Rel.dat=subset(DATA.bio,SHEET_NO%in%unique(Diet$SHEET_NO),select=c(SHEET_NO,date.x,Mid.Lat, Mid.Long, Lat.round, Long.round))


#3. Export Reproduction and Diet data for Georgina -------------------------------------------------------------------------
write.csv(Look.at.repro%>%left_join(DATA.bio%>%
                                      filter(SHEET_NO%in%unique(Look.at.repro$SHEET_NO))%>%
                                      distinct(SHEET_NO,Mid.Lat,Mid.Long,date.x),by='SHEET_NO'),
          handl_OneDrive('Analyses/Surveys/Naturaliste_longline/Data for Georgina/Look.at.repro.csv'),row.names = F)
write.csv(Fem.dat, handl_OneDrive('Analyses/Surveys/Naturaliste_longline/Data for Georgina/Fem.dat.csv'),row.names = F)
write.csv(Diet%>%left_join(DATA.bio%>%
                             filter(SHEET_NO%in%unique(Diet$SHEET_NO))%>%
                             distinct(SHEET_NO,Mid.Lat,Mid.Long,date.x),by='SHEET_NO'),
          handl_OneDrive('Analyses/Surveys/Naturaliste_longline/Data for Georgina/Diet.csv'),row.names = F)
write.csv(Rel.dat,handl_OneDrive('Analyses/Surveys/Naturaliste_longline/Data for Georgina/Rel.dat.csv'),row.names = F)


# Analsyes are done in 'Surveys_analysis_other_Georgina.r'