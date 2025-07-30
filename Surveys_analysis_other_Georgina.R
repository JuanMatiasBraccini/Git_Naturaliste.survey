#1. DATA and Control sections ------------------------------------------------------------------------- 
library(tidyverse)
library(grDevices)
library(data.table)

#Define user
User="Matias"
#User="Georgina"
  
if(!exists('Look.at.repro'))
{
  if(User=="Matias") handl_OneDrive=function(x)paste('C:/Users','myb','OneDrive - Department of Primary Industries and Regional Development/Matias',x,sep='/')
  if(User=="Georgina") handl_OneDrive=function(x)paste('add path to where files will be stored or saved',x,sep='/')
  
  Look.at.repro=read.csv(handl_OneDrive('Analyses/Surveys/Naturaliste_longline/Data for Georgina/Look.at.repro.csv'))
  Fem.dat=read.csv(handl_OneDrive('Analyses/Surveys/Naturaliste_longline/Data for Georgina/Fem.dat.csv'))
  dat.size.brth=read.csv(handl_OneDrive("Analyses/Surveys/Naturaliste_longline/Data for Georgina/dat.size.brth.csv"))
  Diet=read.csv(handl_OneDrive('Analyses/Surveys/Naturaliste_longline/Data for Georgina/Diet.csv'))
  Rel.dat=read.csv(handl_OneDrive('Analyses/Surveys/Naturaliste_longline/Data for Georgina/Rel.dat.csv'))
  Spec.names=read.csv(handl_OneDrive('Analyses/Surveys/Naturaliste_longline/Data for Georgina/Spec.names.csv'),stringsAsFactors=F)
  Predator_TL=read.csv(handl_OneDrive("Data/Naturaliste/Predator_TL.csv"))
  Prey_TL=read.csv(handl_OneDrive("Data/Naturaliste/Prey_TL.csv"))
  Prey_cnvrsn=read.csv(handl_OneDrive("Data/Naturaliste/Prey.csv"))
}

n.min=5  #min sample size per species
Do.rep="YES"
Do.diet="YES"

if(!exists('fn.fig'))
{
  #handy function for exporting figures
  fn.fig=function(NAME,Width,Height)
  {
    if(Do.tiff=="YES") tiff(file=paste(NAME,".tiff",sep=""),width=Width,height=Height,units="px",res=300,compression="lzw")
    if(Do.jpeg=="YES") jpeg(file=paste(NAME,".jpeg",sep=""),width=Width,height=Height,units="px",res=300)
  }
  #choose if doing .jpeg or .tiff figures
  Do.jpeg="YES"
  Do.tiff="NO"
}

#2. Reproduction -------------------------------------------------------------------------    
#MISSING: 
#         EMBLEN_ there's several so look at variability and report! (use mean..?)
#         TL-FL conversion, add to definition of variable 'Length'
#         Check if removing species or not, see 'remove.sp'
remove.sp=FALSE

if(Do.rep=="YES")
{
  if(User=="Matias") setwd(handl_OneDrive('Analyses/Surveys/Naturaliste_longline/outputs/Reproduction'))
  if(User=="Georgina") setwd(handl_OneDrive('xxx'))
  
  #remove species with ID issues
  if(remove.sp)
  {
    sp.ID.issues=c("Guitarfish & shovelnose ray","Wobbegong (general)","Western Wobbegong","Spotted wobbegong",
                   "Angel Shark (general)","Banded wobbegong","Blacktip sharks","Sawsharks","Cobbler Wobbegong")
    Look.at.repro=subset(Look.at.repro,!COMMON_NAME%in%sp.ID.issues)
    Fem.dat=subset(Fem.dat,!COMMON_NAME%in%sp.ID.issues)
  }

  
  #Some final data manipulations
  #note:Length is TL or DW (rays)
  Ray.sp=c('ER','SM','SR')
  SP.TL.rather.than.FL=Look.at.repro%>%filter(!is.na(TL) & is.na(FL))
  TL.sp=unique(SP.TL.rather.than.FL%>%filter(!SPECIES%in%Ray.sp)%>%pull(SPECIES))
  Look.at.repro=Look.at.repro%>%
                  mutate(Length=ifelse(SPECIES%in%TL.sp,TL,
                                ifelse( SPECIES%in%Ray.sp,Disc.width,
                                        FL)))
  Fem.dat=Fem.dat%>%
    mutate(Length=ifelse(SPECIES%in%TL.sp,TL,
                         ifelse( SPECIES%in%Ray.sp,Disc.width,
                                 FL)))
  
  
  #1. Umbilical scars    #data only reliable for a few species
  bin=10  #size bins
  ID.v=match("UMBIL_SCAR",names(Look.at.repro))
  Select.sp=with(subset(Look.at.repro,!UMBIL_SCAR=="N"),as.matrix(table(SPECIES,UMBIL_SCAR)))
  Select.sp=Select.sp[rowSums(Select.sp)>1,]
  UMB_sca_sp=rownames(Select.sp)
  fn.fig('Fig_Umb_scar',2000, 2400)
  par(mfcol=n2mfrow(length(UMB_sca_sp)),mar=c(1,4,3,1.5),oma=c(3,.1,.1,.1),las=1,mgp=c(.9,.6,0))
  for(s in 1:length(UMB_sca_sp))
  {
    d=subset(Look.at.repro,SPECIES==UMB_sca_sp[s]&!UMBIL_SCAR=="N")
    SPnme=unique(d$COMMON_NAME)
    d$Length=floor(d$Length/bin)*bin
    TAB=as.matrix(table(d[,ID.v],d$Length))
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
  mtext("Length (cm)",1,1.5,outer=T,cex=1.5)
  mtext("Frequency",2,-1.8,outer=T,las=3,cex=1.5)
  dev.off()
  
  
  #2. MALE ANALYSES
  
  #2. 1. Clasper length, calcification and maturity
  #note: logistic ogive of the form Prop Mature=pmax/(1+exp((dat-a)/b)); a=inflex; b=slope
  
  Look.at.repro$CLASP_CALC=with(Look.at.repro,ifelse(SPECIES=="GN"&FL<140,"N",CLASP_CALC))
  Select.sp=with(subset(Look.at.repro,!is.na(CLASPLENTH) & SEX=="M"),table(SPECIES))
  Select.sp=Select.sp[Select.sp>=n.min]
  CLas_cal_sp=rownames(Select.sp)
  fn.fig('Fig_male_clas_length',2400, 2200)    
  par(mfcol=n2mfrow(length(CLas_cal_sp)),mar=c(1,1,2,.9),oma=c(3,3,.1,.1),las=1,mgp=c(.9,.6,0))
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
    
    
    plot(d$Length,d$CLASPLENTH,pch=d$CLASP_CALC_pt,bg=d$CLASP_CALC_col,cex=1.25,
         ylab="",xlab="")
    mtext(paste(SPnme," (n=",nrow(d),")",sep=""),3,0,cex=.72)
  }
  plot(1,ann=F,xaxt='n',yaxt='n',col="transparent")
  box(col="white")
  legend('top',c("Y","P","N","Unknown"),pch=c(21,21,21,3),pt.bg=c("black","grey60","white","black"),
         bty='n',cex=1.25,title="Calcification")
  mtext("Length (cm)",1,1.5,outer=T,cex=1.5)
  mtext("Clasper length (cm)",2,1.2,outer=T,las=3,cex=1.5)
  dev.off()
  
  
  #2.2 get L50
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
  par(mfcol=n2mfrow(length(CLas_cal_sp)),mar=c(1,1,2,.9),oma=c(3,3,.1,.1),las=1,mgp=c(.9,.6,0))
  for(s in 1:length(CLas_cal_sp))
  {
    d=subset(Look.at.repro,SPECIES==CLas_cal_sp[s] & CLASPLENTH>0)
    d=subset(d,!is.na(Length))
    d$N=d$CLASPLENTH/max(d$CLASPLENTH)
    
    #consider calcification
    d$N=with(d,ifelse(CLASP_CALC=="N",0,N))
    
    #add dummies to anchor model
    BirtH=dat.size.brth[s,]$Birth
    Extra.n=5
    d.optim=d[1:Extra.n,]
    d.optim[,]=NA
    d.optim$Length=seq(BirtH,BirtH*1.1,length.out=Extra.n)
    d.optim$N=rep(0,Extra.n)
    
    d=rbind(d,d.optim)
    d=subset(d,!is.na(N))
    
    #fit model
    if(Fitting.approach=='nls')
    {
      mod <- nls(N~1/(1+exp((Length-inflex)/slope)), data=d,
                 start=list(inflex=quantile(d$Length,probs=.6), slope=-2),
                 nls.control(maxiter = 100, tol = 1e-06, minFactor = 1/1024,
                             printEval = FALSE, warnOnly = T))
      #L50 and L95 approach
      #note: L50 == inflex point
      # mod1 <- nls(N~1/(1+exp(-log(19)*((Length-L50)/(L95-L50)))), data=d, 
      #                 start=list(L50=quantile(d$Length,probs=.6), L95=quantile(d$Length,probs=.8)),
      #                 nls.control(maxiter = 100, tol = 1e-06, minFactor = 1/1024,
      #                             printEval = FALSE, warnOnly = T))
      
      NM=unique(d$COMMON_NAME)[1]
      newD=data.frame(Length=seq(min(d$Length),max(d$Length)))
      PRED=predict(mod,newdata=newD)
      plot(d$Length,d$N,pch=19,ylab="",xlab="",main="")
      mtext(paste(NM," (n=",nrow(d),")",sep=""),3,0,cex=.85)
      lines(newD$Length,PRED,lwd=2,col="grey50")
      
      TABL=as.data.frame(summary(mod)$coefficients[,c(1:2,4)])
      TABL$SPECIES=NM
      TABL=cbind(TABL[1,c(4,1:3)],TABL[2,1:3])
      names(TABL)[2:7]=c("inflex","inflex.SE","P","slope","slope.SE","P")
      TABL$"a (SE)"=paste(round(TABL$inflex,2)," (",round(TABL$inflex.SE,2),")",sep="")
      TABL$"b (SE)"=paste(round(TABL$slope,2)," (",round(TABL$slope.SE,2),")",sep="")
      Store[[s]]=TABL[,c(1,8,4,9,7)]
      
    }
    if(Fitting.approach=='optim')
    {
      #theta=c(inflx=quantile(d$Length,probs=.6),slop=-5)
      theta=c(L50=quantile(d$Length,probs=.6),L95=quantile(d$Length,probs=.8))
      #fit=optim(theta,objfun, hessian=T, control=c(trace=1, maxit=10000))
      fit=optim(theta,objfun, hessian=T, method="L-BFGS-B",
                lower=c(quantile(d$Length,probs=.2),quantile(d$Length,probs=.4)),
                upper=c(quantile(d$Length,probs=.9),quantile(d$Length,probs=.9)),control=c(maxit=10000))
      # fit=optim(theta,objfun, hessian=T, method="L-BFGS-B",lower=c(quantile(d$Length,probs=.2),-50),
      #           upper=c(quantile(d$Length,probs=.9),-1),control=c(maxit=10000))
      # 
      NM=unique(d$COMMON_NAME)[1]
      newD=seq(min(d$Length),max(d$Length))
      #PRED=fn.logis(newD,1,fit$par[1],fit$par[2])
      PRED=fn.logis_L50(newD,1,fit$par[1],fit$par[2])
      
      plot(d$Length,d$N,pch=19,ylab="",xlab="",main="")
      mtext(paste(SPnme," (n=",nrow(d),")",sep=""),3,0,cex=.85)
      lines(newD,PRED,lwd=2,col="grey50")
      v_ob=solve(fit$hessian)	#variance covariance matrix
      std_ob=sqrt(diag(v_ob))
      TABL=data.frame(SPECIES=NM,L50=fit$par[1],L50_SE=std_ob[1],L95=fit$par[2],L95_SE=std_ob[2])
      Store[[s]]=TABL
    }
    
    
    if(show.classification=="YES")
    {   classify_data = classify_mature(d, varNames = c("Length", "CLASPLENTH"), 
                                        varSex = "SEX", selectSex = "M", method = "ld")
    par(mfrow = c(2,2))
    plot(classify_data)
    plot(classify_data, xlab = "Length", ylab = "CLASPLENTH")
    plot(classify_data, xlab = "Length", ylab = "CLASPLENTH",    col = c(2, 3), pch = c(5, 6))
    
    plot(classify_data, xlab = "Length", ylab = "CLASPLENTH", col = c(2, 3), pch = c(5, 6), lty_lines = c(1, 2), lwd_lines = c(1, 3), 
         cex = c(1, 3), main = "Classification")
    my_ogive_bayes = morph_mature(classify_data, method = "bayes", niter = 1000)
    print(my_ogive_bayes)
    par(mfrow = c(2,2))
    plot(my_ogive_bayes, xlab = "Length", ylab = "Proportion mature", col = c("blue", "red"))
    }
  }
  mtext("Length (cm)",1,1.5,outer=T,cex=1.5)
  mtext("Proportion mature",2,1.2,outer=T,las=3,cex=1.5)
  dev.off()
  
  Store=do.call(rbind,Store)
  write.csv(Store,"L50_pars_male.csv",row.names=F)
  
  #2.3. Running sperm for those with calcified claspers
  Select.sp=with(subset(Look.at.repro,!is.na(RUN_SPERM) & SEX=="M"),table(SPECIES))
  Select.sp=Select.sp[Select.sp>=n.min]
  Sperm_sp=rownames(Select.sp)
  
  fn.fig('Fig_male_running_sperm',2400, 2400)    
  par(mfcol=n2mfrow(length(Sperm_sp)),mar=c(1,1,2,1.2),oma=c(3,3,.1,.1),las=1,mgp=c(.9,.6,0))
  for(s in 1:length(Sperm_sp))
  {
    d=subset(Look.at.repro,SPECIES==Sperm_sp[s] & !is.na(RUN_SPERM) & SEX=="M" & CLASP_CALC=="Y")
    d=subset(d,!is.na(Length))
    SPnme=unique(d$COMMON_NAME)[1]
    d$COL=with(d,ifelse(RUN_SPERM=="Y","black",ifelse(RUN_SPERM=="N","grey80",NA)))
    plot(d$Month,d$Length,pch=21,col='transparent',cex=1.25,ylab="",xlab="",xlim=c(0.9,13))
    d1=subset(d,RUN_SPERM=="Y")
    points(d1$Month,d1$Length,pch=21,bg=d1$COL,cex=1.15)
    d2=subset(d,RUN_SPERM=="N")
    points(d2$Month+1.2,d2$Length,pch=21,bg=d2$COL,cex=1.15)
    NN=nrow(d)
    mtext(paste(SPnme," (n=",NN,")",sep=""),3,0,cex=.85)
  }
  plot(1,ann=F,xaxt='n',yaxt='n',col="transparent")
  box(col="white")
  legend("center",c("Y","N"),pch=21,pt.bg=c("black","grey80"),bty='n',title='Sperm presence',cex=1.5)
  
  mtext("Month",1,1.5,outer=T,cex=1.5)
  mtext("Length (cm)",2,1.2,outer=T,las=3,cex=1.5)
  dev.off()
  
  #3. FEMALE ANALYSES  
 # Fem.dat=subset(Fem.dat,!SPECIES%in%c("BT","PJ","SD","SH","SW","WB","WD","WS","WW"))
  
  #3.1. Embryo number 
  #note: NO_UNDEVELOPED refers to embryos that did not develop. Only 5 observations, don't report
  
  #if emblen_1 is na, see if data for other emblens..
  Fem.dat$EMBLEN_1=with(Fem.dat,ifelse(is.na(EMBLEN_1) & !is.na(EMBLEN_2),EMBLEN_2,EMBLEN_1))
  
  
  #fit linear model of the form num. emb = b*Length + a
  Select.sp=with(subset(Fem.dat,!is.na(NO_EMBRYOS)),table(SPECIES))
  Select.sp=Select.sp[Select.sp>=n.min]
  Embryo_sp=rownames(Select.sp)
  
  Store=vector('list',length(Embryo_sp))  
  names(Store)=Embryo_sp
  
  fn.fig('Fig_female_embryo_number',2400, 2400)    
  par(mfcol=n2mfrow(length(Embryo_sp)),mar=c(1,1,2,1.2),oma=c(3,3,.1,.1),las=1,mgp=c(.9,.6,0))
  for(s in 1:length(Embryo_sp))
  {
    d=subset(Fem.dat,SPECIES==Embryo_sp[s] & !is.na(NO_EMBRYOS))
    d=subset(d,!is.na(Length))
    SPnme=unique(d$COMMON_NAME)[1]
    plot(d$Length,d$NO_EMBRYOS,pch=19,cex=1.25,ylab="",xlab="",ylim=c(0,max(d$NO_EMBRYOS)))
    NN=nrow(d)
    mtext(paste(SPnme," (n=",NN,")",sep=""),3,0,cex=.85)
    mod=lm(NO_EMBRYOS~Length,d)
    PRED.d=data.frame(Length=seq(min(d$Length),max(d$Length)))
    PRED=predict(mod,newdata=PRED.d)
    lines(PRED.d$Length,PRED,lwd=2,col="grey50")
    
    TABL=as.data.frame(summary(mod)$coefficients[,c(1:2,4)])
    TABL$SPECIES=SPnme
    TABL=cbind(TABL[1,c(4,1:3)],TABL[2,1:3])
    names(TABL)[2:7]=c("inflex","inflex.SE","P","slope","slope.SE","P")
    TABL$"a (SE)"=paste(round(TABL$inflex,2)," (",round(TABL$inflex.SE,2),")",sep="")
    TABL$"b (SE)"=paste(round(TABL$slope,2)," (",round(TABL$slope.SE,2),")",sep="")
    Store[[s]]=TABL[,c(1,8,4,9,7)]
  }
  mtext("Length (cm)",1,1.5,outer=T,cex=1.5)
  mtext("Number of embryos",2,1.2,outer=T,las=3,cex=1.5)
  dev.off()
  Store=do.call(rbind,Store)
  write.csv(Store,"Female_pups_regression.csv",row.names=F) 
  
  #Embryo length vs time
  Fem.dat$EMBLEN_1=with(Fem.dat,ifelse(UTERINESTG==4,0,EMBLEN_1))
  Fem.dat$EMBLEN_1=with(Fem.dat,ifelse(!UTERINESTG==4 & EMBLEN_1==0,NA,EMBLEN_1))
  
  Select.sp=with(subset(Fem.dat,!is.na(EMBLEN_1)),table(SPECIES))
  Select.sp=Select.sp[Select.sp>=n.min]
  Embryo_sp=rownames(Select.sp)
  fn.fig('Fig_female_embryo_length',2400, 2400)    
  par(mfcol=n2mfrow(length(Embryo_sp)),mar=c(1,1,2,1.2),oma=c(3,3,.1,.1),las=1,mgp=c(.9,.6,0))
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
  legend('bottomleft',c("4","5"),pch=21,pt.bg=c("grey65","black"),bty='n',cex=.9,
         title="Uterus condition")
  dev.off()
  
  
  #Embryo length vs Max ova diameter
  Select.sp=with(subset(Fem.dat,SPECIES%in%Embryo_sp&!is.na(EMBLEN_1) & !is.na(MAXOVRYDIA)),table(SPECIES))
  Select.sp=Select.sp[Select.sp>=n.min]
  Embryo_sp_ova=rownames(Select.sp)
  
  fn.fig('Fig_female_embryo_length_ova_diam',1200, 2400)    
  par(mfcol=n2mfrow(length(Embryo_sp_ova)),mar=c(1,1,2,1.2),oma=c(3,3,.1,.1),las=1,mgp=c(.9,.6,0))
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
  legend('topright',c("4","5"),pch=21,pt.bg=c("grey65","black"),bty='n',cex=0.9,title="Uterus condition")
  dev.off()

    
  #Embryo Ova diameter vs time
  Select.sp=with(subset(Fem.dat,!is.na(MAXOVRYDIA)),table(SPECIES))
  Select.sp=Select.sp[Select.sp>=n.min]
  Ova_sp=rownames(Select.sp)
  fn.fig('Fig_female_ova_diam',2400, 2400)    
  par(mfcol=n2mfrow(length(Ova_sp)),mar=c(1,1,2,1.2),oma=c(3,3,.1,.1),las=1,mgp=c(.9,.6,0))
  for(s in 1:length(Ova_sp))
  {
    d=subset(Fem.dat,SPECIES==Ova_sp[s] & !is.na(MAXOVRYDIA))
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
  
  
  #Embryo length -Ova diameter vs time
  Select.sp=with(subset(Fem.dat,!is.na(MAXOVRYDIA) & !is.na(EMBLEN_1)),table(SPECIES))
  Select.sp=Select.sp[Select.sp>=n.min]
  Ova_sp_emb.len=rownames(Select.sp)

  fn.fig('Fig_female_emblen_ova_diam_time',1400, 2400)    
  par(mfcol=n2mfrow(length(Ova_sp_emb.len)),mar=c(1,1,2,4),oma=c(3,3,.1,1),las=1,mgp=c(.9,.6,0))
  for(s in 1:length(Ova_sp_emb.len))
  {
    d=subset(Fem.dat,SPECIES==Ova_sp_emb.len[s] &!is.na(MAXOVRYDIA))
    d$MAXOVRYDIA=d$MAXOVRYDIA/10   #convert to cm as MOD measured in mm
    SPnme=unique(d$COMMON_NAME)[1]
    with(subset(d,!is.na(MAXOVRYDIA)),plot(Month,MAXOVRYDIA,pch=21,bg="grey30",cex=1.75,ylab="",xlab="",xlim=c(0.9,13),
                                           ylim=c(0,max(MAXOVRYDIA,na.rm=T)),cex.axis=1.25))
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
  
  
  
  #Maturity 
  Select.sp=with(subset(Fem.dat,!is.na(UTERINESTG)),table(SPECIES))
  Select.sp=Select.sp[Select.sp>=n.min]
  MAT_sp=rownames(Select.sp)
  #MAT_sp=MAT_sp[-match(c("CW","GU","MS","WC","TA"),MAT_sp)]
  Select.sp=with(subset(Fem.dat,!is.na(GON_STAGE)),table(SPECIES))
  id=which(!rownames(Select.sp)%in%MAT_sp)    #no species had GI info but no UI info
  
  #2. get L50
  Fem.dat=merge(Fem.dat,dat.size.brth,by="SPECIES",all.x=T)
  
  Store=vector('list',length(MAT_sp))  
  names(Store)=MAT_sp
  fn.fig('Fig_female_maturity_fit',1800, 2400)    
  par(mfcol=n2mfrow(length(MAT_sp)),mar=c(1,1,2,.9),oma=c(3,3,.1,.1),las=1,mgp=c(.9,.6,0))
  for(s in 1:length(MAT_sp))
  {
    d=subset(Fem.dat,SPECIES==MAT_sp[s])
    d=subset(d,!is.na(Length))
    d$N=with(d,ifelse(UTERINESTG%in%3:6 ,1,ifelse(UTERINESTG%in%1:2,0,NA)))
    d$N=with(d,ifelse(is.na(N) & GON_STAGE==3,1,N))
    d$N=with(d,ifelse(is.na(N) & GON_STAGE%in%c(1,2),0,N))
    
    #add dummies to anchor model
    BirtH=dat.size.brth[s,]$Birth
    Extra.n=5
    d.optim=d[1:Extra.n,]
    d.optim[,]=NA
    d.optim$Length=seq(BirtH,BirtH*1.1,length.out=Extra.n)
    d.optim$N=rep(0,Extra.n)
    
    d=rbind(d,d.optim)
    d=subset(d,!is.na(N))
    
    if(MAT_sp[s]=="TG") d=subset(d,!(N==1 & Length<200))
    NM=unique(d$COMMON_NAME)[1]
    #fit model
    #binominal glm
    # mod=glm(N~Length,d,family=binomial())
    # a=predict(mod,newdata=data.frame(Length=seq(min(d$Length),max(d$Length))),type='response')
    # lines(seq(min(d$Length),max(d$Length)),a)
    
    #nls approach
    if(Fitting.approach=='nls')
    {
      mod <- nls(N~1/(1+exp((Length-inflex)/slope)), data=d,
                 start=list(inflex=quantile(d$Length,probs=.6), slope=-2),
                 nls.control(maxiter = 100, tol = 1e-06, minFactor = 1/1024,
                             printEval = FALSE, warnOnly = T))
      #L50 and L95 approach
      #note: L50 == inflex point
      # mod1 <- nls(N~1/(1+exp(-log(19)*((Length-L50)/(L95-L50)))), data=d, 
      #                 start=list(L50=quantile(d$Length,probs=.6), L95=quantile(d$Length,probs=.8)),
      #                 nls.control(maxiter = 100, tol = 1e-06, minFactor = 1/1024,
      #                             printEval = FALSE, warnOnly = T))
      
      
      newD=data.frame(Length=seq(min(d$Length),max(d$Length)))
      PRED=predict(mod,newdata=newD)
      # d$COL=with(d,ifelse(UTERINESTG==4,"grey50",ifelse(UTERINESTG==5,"grey70",
      #             ifelse(UTERINESTG==3,"grey30",ifelse(UTERINESTG==6,"grey85","white")))))                                   
      d$COL=1
      plot(d$Length,d$N,pch=21,ylab="",xlab="",main="",bg=d$COL,ylim=c(0,1))
      mtext(paste(NM," (n=",nrow(d),")",sep=""),3,0,cex=.85)
      lines(newD$Length,PRED,lwd=2,col="grey50")
      
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
      #theta=c(inflx=quantile(d$Length,probs=.6),slop=-5)
      theta=c(L50=quantile(d$Length,probs=.6),L95=quantile(d$Length,probs=.8))
      #fit=optim(theta,objfun, hessian=T, control=c(trace=1, maxit=10000))
      fit=optim(theta,objfun, hessian=T, method="L-BFGS-B",
                lower=c(quantile(d$Length,probs=.2),quantile(d$Length,probs=.4)),
                upper=c(quantile(d$Length,probs=.9),quantile(d$Length,probs=.9)),control=c(maxit=10000))
      # fit=optim(theta,objfun, hessian=T, method="L-BFGS-B",lower=c(quantile(d$Length,probs=.2),-50),
      #           upper=c(quantile(d$Length,probs=.9),-1),control=c(maxit=10000))
      # 
      NM=unique(d$COMMON_NAME)[1]
      newD=seq(min(d$Length),max(d$Length))
      #PRED=fn.logis(newD,1,fit$par[1],fit$par[2])
      PRED=fn.logis_L50(newD,1,fit$par[1],fit$par[2])
      
      d$COL=with(d,ifelse(UTERINESTG==4,"grey50",ifelse(UTERINESTG==5,"grey70",
                                                        ifelse(UTERINESTG==3,"grey30",ifelse(UTERINESTG==6,"grey85","white")))))                                   
      
      plot(d$Length,d$N,pch=21,ylab="",xlab="",main="",bg=d$COL,ylim=c(0,1))
      
      mtext(paste(SPnme," (n=",nrow(d),")",sep=""),3,0,cex=.85)
      lines(newD,PRED,lwd=2,col="grey50")
      v_ob=solve(fit$hessian)	#variance covariance matrix
      std_ob=sqrt(diag(v_ob))
      TABL=data.frame(SPECIES=NM,L50=fit$par[1],L50_SE=std_ob[1],L95=fit$par[2],L95_SE=std_ob[2])
      Store[[s]]=TABL
    }
    
    
    if(show.classification=="YES")
    {   classify_data = classify_mature(d, varNames = c("Length", "CLASPLENTH"), 
                                        varSex = "SEX", selectSex = "M", method = "ld")
    par(mfrow = c(2,2))
    plot(classify_data)
    plot(classify_data, xlab = "Length", ylab = "CLASPLENTH")
    plot(classify_data, xlab = "Length", ylab = "CLASPLENTH",    col = c(2, 3), pch = c(5, 6))
    
    plot(classify_data, xlab = "Length", ylab = "CLASPLENTH", col = c(2, 3), pch = c(5, 6), lty_lines = c(1, 2), lwd_lines = c(1, 3), 
         cex = c(1, 3), main = "Classification")
    my_ogive_bayes = morph_mature(classify_data, method = "bayes", niter = 1000)
    print(my_ogive_bayes)
    par(mfrow = c(2,2))
    plot(my_ogive_bayes, xlab = "Length", ylab = "Proportion mature", col = c("blue", "red"))
    }
  }
  
  # plot(1,ann=F,xaxt='n',yaxt='n',col="transparent")
  # box(col="white")
  # legend('top',c("Immature","UI 3","UI 4","UI 5","UI 6"),pch=21,
  #        pt.bg=c("white","grey30","grey50","grey70","grey85"),bty='n',cex=1.25)
  mtext("Length (cm)",1,1.5,outer=T,cex=1.5)
  mtext("Proportion mature",2,1.2,outer=T,las=3,cex=1.5)
  dev.off()
  
  Store=do.call(rbind,Store)
  write.csv(Store,"L50_pars_female.csv",row.names=F)
}



# 3. Diet -------------------------------------------------------------------------
#note: diet analyses based on frequency of occurrence only (weight/number prey not available/reliable
#  MISSING: not enough samples if using tropics only (survey); separate tropical from subtropical;
#             report size distribution; map of where samples from
#             Separate some species into large, small? (e.g dusky) see Papastamatiou et al 2006
#        minimum sample size: 10
#            Also, among species comparisons should be meaningful (e.g. same location/shot)
if(Do.diet=="YES")
{
  Diet=subset(Diet,!COMMON_NAME=="Parrotfish (general)")
  
  if(User=="Matias") HNDL.diet=handl_OneDrive("Analyses/Surveys/Naturaliste_longline/outputs/Diet/")
  if(User=="Georgina") HNDL.diet="whatever"
  
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
  Diet=Diet[,-match("COMMON_NAME",names(Diet))]
  Diet=merge(Diet,Spec.names,by="SPECIES",all.x=T)
  Diet$SPECIES=as.character(Diet$SPECIES)
  
  #tropical
  Tropics=-23.43697
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
  par(mfrow=n2mfrow(length(SpEcis)), mar=c(.1,.1,.1,0.5), mgp = c(1.5, 0.3, 0), tck = -0.01)
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
  
  
  #Missing:
  
  #Dietary overlap based on Simplified Morisita index with diet expressed as proportions
  
  
  #Generalist-Specialist spectrum / Prey diversity indices / Omnivory index  / Size of shark  /home range ( Roff et al 2016)
  
  
  #Map density distribution of shots to put study in context
  
  
  #size distribution of species
  
  
}








