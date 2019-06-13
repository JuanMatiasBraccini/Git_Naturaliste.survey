#1.12.5.1 Monte Carlo approach
if(do.MC=="YES")
{
  fun.MC.yr.pred=function(d)
  {
    #create new data
    Covar.pos=as.matrix(vcov(d))
    set.seed(999);Pos.pars.rand=rmvnorm(niter,mean=coef(d),sigma=Covar.pos)    
    MC.preds=matrix(nrow=niter,ncol=length(NEWDATA$year))
    
    for(n in 1:niter)
    {
      model=d
      model$coefficients$count=Pos.pars.rand[n,1:length(model$coefficients$count)]
      model$coefficients$zero=Pos.pars.rand[n,1:length(model$coefficients$zero)]
      MC.preds[n,]=predict(model,newdata=NEWDATA,type='response')
    }
    
    PRD=data.frame(Year=NEWDATA$year,
                   MEAN=colMeans(MC.preds,na.rm=T),
                   SD=apply(MC.preds,2,function(x) sd(x,na.rm=T)),
                   LOW=apply(MC.preds, 2, function(x) quantile(x, 0.025,na.rm=T)),
                   UP=apply(MC.preds, 2, function(x) quantile(x, 0.975,na.rm=T)))
    PRD$CV=PRD$SD/PRD$MEAN*100
    
    return(PRD)
    
  }
  Store.MC=Store.MC.all.stations=vector('list',niter)
  
  system.time(for(i in Species.cpue)   # 0.06 secs per iteration
  {
    Store.MC[[i]]=fun.MC.yr.pred(d=Store[[i]]$Fit)
    Store.MC.all.stations[[i]]=fun.MC.yr.pred(d=Store.all.stations[[i]]$Fit)
  })
}

if(do.boot=="YES")
{
  #bootstrap the data sets
  fn.yr.boot=function(dd)
  {
    n=nrow(dd)
    return(dd[sample(1:n,n,replace=T),])
  }
  fn.boot=function(DAT)
  {
    pos.dat=subset(DAT,Catch.Target>0)
    YR=sort(unique(pos.dat$year))
    POS=vector('list',length(YR))
    for(y in 1:length(YR)) POS[[y]]=fn.yr.boot(dd=subset(pos.dat,year==YR[y]))
    POS=do.call(rbind,POS)
    
    zero.dat=subset(DAT,Catch.Target==0)
    YR=sort(unique(zero.dat$year))
    ZERO=vector('list',length(YR))
    for(y in 1:length(YR)) ZERO[[y]]=fn.yr.boot(dd=subset(zero.dat,year==YR[y]))
    ZERO=do.call(rbind,ZERO)
    
    return(rbind(POS,ZERO))
  }
  # fn.boot=function(DAT)
  # {
  #   n=nrow(DAT)
  #   id=sample(1:n,n,replace=T)
  #   DAT.resampled=DAT[id,]
  #   return(DAT.resampled)
  # }
  
  STORE.BOOT=vector('list',N.species)
  names(STORE.BOOT)=names(DATA.list)
  STORE.BOOT.all.stations=STORE.BOOT
  
  #fit models to bootstrapped data
  system.time(for(i in Species.cpue)  #takes 2.5 secs per iteration 
  {
    Store.boot=Store.boot.all.stations=vector('list',n.boot)
    for(k in 1:n.boot)
    {
      #fixed stations
      DAT=fn.boot(subset(DATA.list[[i]],FixedStation=="YES",
                         select=c(Catch.Target,year,BOTDEPTH,Mid.Lat,N.hooks.Fixed,SOAK.TIME)))
      mod=try(fit.best(DAT,BEST.model[[i]],ERROR.st[[i]]),TRUE)  #deal with dodgy boot data
      if(isTRUE(class(mod)=="try-error")) { next } else { Store.boot[[k]]=mod } 
      rm(mod)
      #All stations
      # DAT=fn.boot(subset(DATA.list[[i]],
      #                    select=c(Catch.Target,year,BOTDEPTH,Mid.Lat,N.hooks.Fixed,SOAK.TIME)))
      # mod=try(fit.best(DAT,BEST.model[[i]],ERROR.st[[i]]),TRUE)  #deal with dodgy boot data
      # if(isTRUE(class(mod)=="try-error")) { next } else { Store.boot.all.stations[[k]]=mod } 
    }
    
    STORE.BOOT[[i]]=Store.boot
    #STORE.BOOT.all.stations[[i]]=Store.boot.all.stations
  })
  
  #predict each boot run   
  #Year effect
  PRED.CPUE=vector('list',N.species)
  names(PRED.CPUE)=names(DATA.list)
  PRED.CPUE.all.station=PRED.CPUE
  for(i in Species.cpue)
  {
    MOD=STORE.BOOT[[i]]
    #MOD.all.stations=STORE.BOOT.all.stations[[i]]
    An.Ktch=matrix(nrow=length(MOD),ncol=length(NEWDATA$year))
    colnames(An.Ktch)=NEWDATA$year
    An.Ktch.all.stations=An.Ktch
    for (x in 1:length(MOD))
    {
      a=Pred.Model(MOD[[x]]$Fit,NEWDATA)$PREDS
      if(!is.null(a))An.Ktch[x,]=a
      
      # a=Pred.Model(MOD.all.stations[[x]]$Fit,NEWDATA)$PREDS
      # if(!is.null(a))An.Ktch.all.stations[x,]=a
    }
    
    PRD=data.frame(Year=NEWDATA$year,
                   MEAN=apply(An.Ktch, 2, function(x) median(x, na.rm=T)),
                   SD=apply(An.Ktch,2,function(x) sd(x,na.rm=T)),
                   LOW=apply(An.Ktch, 2, function(x) quantile(x, 0.025,na.rm=T)),
                   UP=apply(An.Ktch, 2, function(x) quantile(x, 0.975,na.rm=T)))
    PRD$CV=PRD$SD/PRD$MEAN*100
    PRED.CPUE[[i]]=PRD
    
    # PRD=data.frame(Year=NEWDATA$year,
    #           MEAN=apply(An.Ktch.all.stations, 2, function(x) median(x, na.rm=T)),
    #           SD=apply(An.Ktch.all.stations,2,function(x) sd(x,na.rm=T)),
    #           LOW=apply(An.Ktch.all.stations, 2, function(x) quantile(x, 0.025,na.rm=T)),
    #           UP=apply(An.Ktch.all.stations, 2, function(x) quantile(x, 0.975,na.rm=T)))
    # PRD$CV=PRD$SD/PRD$MEAN*100
    # PRED.CPUE.all.station[[i]]=PRD
  }
  
  #Latitude efect                       #ACA: predict for when there's data
  PRED.lat=vector('list',N.species)
  names(PRED.lat)=names(DATA.list)
  #PRED.lat.all.station=PRED.lat
  for(i in Species.cpue)
  {
    dat=Store[[i]]$DAT
    LAT.seq=range(dat$Mid.Lat)
    LAT.seq=seq(round(LAT.seq[1]),round(LAT.seq[2]),.5)
    Yr.ref=names(sort(table(dat$year)))
    Yr.ref=Yr.ref[length(Yr.ref)]
    NEWDATA=data.frame(year=factor(Yr.ref,levels(dat$year)),
                       BOTDEPTH=mean(dat$BOTDEPTH,na.rm=T),
                       Mid.Lat=LAT.seq,
                       N.hooks.Fixed=50,
                       SOAK.TIME=mean(dat$SOAK.TIME,na.rm=T))
    
    MOD=STORE.BOOT[[i]]
    #MOD.all.stations=STORE.BOOT.all.stations[[i]]
    An.Ktch=matrix(nrow=length(MOD),ncol=length(NEWDATA$Mid.Lat))
    colnames(An.Ktch)=NEWDATA$Mid.Lat
    An.Ktch.all.stations=An.Ktch
    for (x in 1:length(MOD))
    {
      a=Pred.Model(MOD[[x]]$Fit,NEWDATA)$PREDS
      if(!is.null(a))An.Ktch[x,]=a
      
      #       a=Pred.Model(MOD.all.stations[[x]]$Fit,NEWDATA)$PREDS
      #      if(!is.null(a))An.Ktch.all.stations[x,]=a
    }
    
    PRD=data.frame(Mid.Lat=NEWDATA$Mid.Lat,
                   MEAN=apply(An.Ktch, 2, function(x) median(x, na.rm=T)),
                   SD=apply(An.Ktch,2,function(x) sd(x,na.rm=T)),
                   LOW=apply(An.Ktch, 2, function(x) quantile(x, 0.025,na.rm=T)),
                   UP=apply(An.Ktch, 2, function(x) quantile(x, 0.975,na.rm=T)))
    PRD$CV=PRD$SD/PRD$MEAN*100
    PRED.lat[[i]]=PRD
    
    # PRD=data.frame(Mid.Lat=NEWDATA$Mid.Lat,
    #                MEAN=apply(An.Ktch.all.stations, 2, function(x) median(x, na.rm=T)),
    #                SD=apply(An.Ktch.all.stations,2,function(x) sd(x,na.rm=T)),
    #                LOW=apply(An.Ktch.all.stations, 2, function(x) quantile(x, 0.025,na.rm=T)),
    #                UP=apply(An.Ktch.all.stations, 2, function(x) quantile(x, 0.975,na.rm=T)))
    # PRD$CV=PRD$SD/PRD$MEAN*100
    # PRED.lat.all.station[[i]]=PRD
  }
  
  #Depth effect
  PRED.z=vector('list',N.species)
  names(PRED.z)=names(DATA.list)
  # PRED.z.all.station=PRED.z
  for(i in Species.cpue)
  {
    dat=Store[[i]]$DAT
    Z.seq=range(dat$BOTDEPTH,na.rm=T)
    Z.seq=round(Z.seq/10)*10
    Z.seq=seq(round(Z.seq[1]),round(Z.seq[2]),10)
    Yr.ref=names(sort(table(dat$year)))
    Yr.ref=Yr.ref[length(Yr.ref)]
    NEWDATA=data.frame(year=factor(Yr.ref,levels(dat$year)),
                       BOTDEPTH=Z.seq,
                       Mid.Lat=mean(dat$Mid.Lat,na.rm=T),
                       N.hooks.Fixed=50,
                       SOAK.TIME=mean(dat$SOAK.TIME,na.rm=T))
    
    MOD=STORE.BOOT[[i]]
    # MOD.all.stations=STORE.BOOT.all.stations[[i]]
    An.Ktch=matrix(nrow=length(MOD),ncol=length(NEWDATA$Mid.Lat))
    colnames(An.Ktch)=NEWDATA$Mid.Lat
    An.Ktch.all.stations=An.Ktch
    for (x in 1:length(MOD))
    {
      a=Pred.Model(MOD[[x]]$Fit,NEWDATA)$PREDS
      if(!is.null(a))An.Ktch[x,]=a
      
      # a=Pred.Model(MOD.all.stations[[x]]$Fit,NEWDATA)$PREDS
      # if(!is.null(a))An.Ktch.all.stations[x,]=a
    }
    
    PRD=data.frame(BOTDEPTH=NEWDATA$BOTDEPTH,
                   MEAN=apply(An.Ktch, 2, function(x) median(x, na.rm=T)),
                   SD=apply(An.Ktch,2,function(x) sd(x,na.rm=T)),
                   LOW=apply(An.Ktch, 2, function(x) quantile(x, 0.025,na.rm=T)),
                   UP=apply(An.Ktch, 2, function(x) quantile(x, 0.975,na.rm=T)))
    PRD$CV=PRD$SD/PRD$MEAN*100
    PRED.z[[i]]=PRD
    
    # PRD=data.frame(BOTDEPTH=NEWDATA$BOTDEPTH,
    #                MEAN=apply(An.Ktch.all.stations, 2, function(x) median(x, na.rm=T)),
    #                SD=apply(An.Ktch.all.stations,2,function(x) sd(x,na.rm=T)),
    #                LOW=apply(An.Ktch.all.stations, 2, function(x) quantile(x, 0.025,na.rm=T)),
    #                UP=apply(An.Ktch.all.stations, 2, function(x) quantile(x, 0.975,na.rm=T)))
    # PRD$CV=PRD$SD/PRD$MEAN*100
    # PRED.z.all.station[[i]]=PRD
  }
}


# #Simulate multivariate data
# N=100
# dd=data.frame(
#   a=c(runif(N,1,5),runif(N,50,100)),
#   b=runif(N*2,0,10),
#   c=c(runif(N,2,10),runif(N,10,30)),
#   d=c(runif(N,10,30),runif(N,2,10)),
#   e=c(runif(N,50,100),runif(N,1,5)),
#   f=runif(N*2,0,100),
#   g=runif(N*2,0,100),
#   h=runif(N*2,0,1),
#   i=runif(N*2,0,1)
# )
# d.res.var.p=dd/rowSums(dd)
# predics=data.frame(
#   year=factor(c(sample(2000:2005,N,replace=T),sample(2006:2010,N,replace=T)),levels=2000:2010),
#   dummy=factor(sample(LETTERS[1:5],N*2,replace=T)),
#   dummy1=factor(sample(LETTERS[6:11],N*2,replace=T)),
#   dummy2=factor(sample(LETTERS[12:20],N*2,replace=T))
# )
# #Traditional Multivariate
# MDS1 <- metaMDS(d.res.var.p, trace = FALSE,parallel=detectCores()-1)
# plot(MDS1, type = "t")
# stressplot(MDS1)
# #see year effect
# nice.MDS.plot(MDS=MDS1,pt.col='transparent',pt.bg="transparent",txt.col="transparent",txt.cex=.9,
#               PRED=predics$year,col.pred=1,pred.cex=.5)
# mod <- adonis(d.res.var.p~., predics)
# mod$aov.tab
# mod <-anosim(d.res.var.p, predics$year,parallel=Cores)
# summary(mod)
# plot(mod)



#library(mvabund)

#Simulate multivariate data
# Stations=1:5
# n.stations=length(Stations)
# Years=1:4
# dd=expand.grid(Station=Stations,year=Years)
# N=nrow(dd)
# dd$Species1=with(dd,ifelse(year==1,runif(n.stations,100,200),
#                            ifelse(year==2,runif(n.stations,100,150),
#                                   ifelse(year==3,runif(n.stations,50,100),
#                                          ifelse(year==4,runif(n.stations,0,50),NA)))))
# 
# dd$Species2=with(dd,ifelse(year==1,runif(n.stations,10,20),
#                            ifelse(year==2,runif(n.stations,10,15),
#                                   ifelse(year==3,runif(n.stations,5,10),
#                                          ifelse(year==4,runif(n.stations,0,5),NA)))))
# dd$Species3=runif(N,0.8,1.2)*(max(dd$Species1)-dd$Species1)
# dd$Species4=runif(N,0.8,1.2)*(max(dd$Species2)-dd$Species2)
# dd=round(dd)
# 
# #add dummy preds
# dd$year=as.factor(dd$year)
# dd$dummy1=runif(N,1,2)
# dd$dummy2=factor(sample(LETTERS[1:5],N,replace=T))
# dd$dummy3=factor(sample(LETTERS[6:10],N,replace=T))
# 
# boxplot(Species1~year,dd)
# boxplot(Species1~dummy2,dd)
# 
# #d.res.var.p=dd[,3:6]/rowSums(dd[,3:6])  #proportional data
# 
# 
# 
# dd.y=mvabund(dd[,3:6])
# boxplot(dd.y)
# plot(dd.y~dd$year)
# #look at mean variance plot
# meanvar.plot(dd.y)
# 
# null <- manyglm(dd.y~1, family="poisson")
# full <- manyglm(dd.y~year+dummy1+dummy2+dummy3, family="poisson", data=dd)
# #if using Catch as res var, can use effort as offset....
# 
# anova(full, nBoot=199, test="wald")
# fun.dev.exp=function(null,modl) (null-modl)/null
# fun.dev.exp(null$deviance,full$deviance)
# plot(full)
# coefplot.manyglm(full, which.Xcoef=2:5)
# 
# mod.fn=function(X) factor(names(rev(sort(table(X)))[1]),levels=levels(X))
# 
# 
# fn.pred=function(what,dd)
# {
#   prdktr=dd[,match(what,names(dd))]
#   if(is.factor(prdktr)) prdktr=factor(unique(prdktr),levels(prdktr)) else
#     prdktr=seq(min(prdktr),max(prdktr),length.out=20)
#   others=
#     NEW=data.frame(
#       year=factor(levels(dd$year)),
#       dummy1=mean(dd$dummy1),
#       dummy2=mod.fn(dd$dummy2),
#       dummy3=mod.fn(dd$dummy3))
#   PRED=predict(full,newdata=NEW,type='response', se.fit =T)
#   
#   Stand=apply(PRED$fit,2,max)
#   PRED.dot=t(t(PRED$fit)/Stand)
#   PRED.dot.min2SE=t(t(PRED$fit-2*PRED$se.fit)/Stand)
#   PRED.dot.plus2SE=t(t(PRED$fit+2+PRED$se.fit)/Stand)
#   plot(0,type='n',axes=FALSE,ann=FALSE,
#        ylim=c(1,ncol(PRED.dot)),xlim=c(1,length(unique(NEW$year))))
#   for(n in 1:nrow(NEW))
#   {
#     x=NEW[n,]
#     y=PRED.dot[n,]
#     #Mean
#     points(rep(x$year,length(y)),1:length(y),cex=PRED.dot[n,]*3,col="black")
#     #CI
#     points(rep(x$year,length(y)),1:length(y),cex=PRED.dot.min2SE[n,]*3,col="grey60") 
#     points(rep(x$year,length(y)),1:length(y),cex=PRED.dot.plus2SE[n,]*3,col="grey60") 
#   }
#   axis(2,at=1:ncol(PRED.dot),colnames(PRED.dot),las=1,cex.axis=.85)
#   axis(1,at=1:length(unique(NEW$year)),unique(NEW$year),las=1,cex.axis=.85)
#   box()
# }
# fn.pred(what="year",dd=dd)



# if(Do.multivariate=="YES")   #ACA: missing: go thru with Corey, review unbalanced design and grouping by station-year, show lat effect?
# {
#   
#   #create data matrix
#   fn.reshp=function(d,Y,TimeVar,IdVAR)
#   {
#     
#     DATA.wide=reshape(d[,match(c(Y,IdVAR,TimeVar),names(d))],v.names=Y,
#                       idvar=IdVAR,timevar=TimeVar,direction="wide")
#     DATA.wide[is.na(DATA.wide)]=0
#     colnames(DATA.wide)=gsub(paste(Y,".",sep=""), "", names(DATA.wide))
#     
#     DATA.wide.cpue=reshape(d[,match(c("cpue",IdVAR,TimeVar),names(d))],v.names="cpue",
#                            idvar=IdVAR,timevar=TimeVar,direction="wide")
#     DATA.wide.cpue[is.na(DATA.wide.cpue)]=0
#     colnames(DATA.wide.cpue)=gsub(paste("cpue",".",sep=""), "", names(DATA.wide.cpue))
#     
#     props=DATA.wide
#     props[,-which(IdVAR%in%names(props))]=props[,-which(IdVAR%in%names(props))]/rowSums(props[,-which(IdVAR%in%names(props))])
#     
#     return(list(catch=DATA.wide,cpue=DATA.wide.cpue,proportion=props))
#   }
#   Dat=fn.reshp(d=Numbers.Station.year,Y=ResVar,TimeVar=MultiVar,IdVAR=Predictors) 
#   #check what type of transformation is needed  
#   #source:  Clarke & Warwick 2001. Change in marine communities: an approach to 
#   #         statistical analysis and interpretation. Chapter 9.
#   Trans.matrix=data.frame(Transf=c("none","root2","root4","log"),
#                           Min=c(0,0.4,0.65,.875),
#                           Max=c(0.4,0.65,0.875,1))  
#   
#   Store.transf=Dat
#   fn.what.trans=function(d)
#   {
#     Log.Mean=log(apply(d,2,mean,na.rm=T))
#     Log.SD=log(apply(d,2,sd,na.rm=T))  
#     
#     mod=lm(Log.SD~Log.Mean)
#     Pred=predict(mod,Log.Mean=seq(min(Log.Mean),max(Log.Mean)),type='response')
#     plot(Log.Mean,Log.SD)
#     lines(Log.Mean,Pred,col=2)
#     SLP=round(coef(mod)[2],3)
#     text(mean(Log.Mean),quantile(Log.SD,probs=0.075),paste("slope=",SLP),col=2,cex=1.25)
#     Trans=as.character(Trans.matrix$Transf[with(Trans.matrix, Min <= SLP & Max >= SLP)])
#     legend("topleft",paste("res var=",names(Store.transf)[s]),bty='n')
#     legend("bottomright",paste("Trans=",Trans),bty='n')
#     return(Trans)
#   }
#   smart.par(n.plots=length(Dat),MAR=c(2,2,1,1),OMA=c(1,1.5,.1,.1),MGP=c(2.5,.7,0))
#   for(s in 1:length(Dat)) Store.transf[[s]]=fn.what.trans(d=Dat[[s]][,-which(IDVAR%in%names(Dat[[s]]))]) 
#   
#   #transform data accordingly based on Store.tranfs 
#   transf.fn=function(d,transf)
#   {
#     if(transf=="none")  d[,-which(IDVAR%in%names(d))]=d[,-which(IDVAR%in%names(d))]
#     if(transf=="root2") d[,-which(IDVAR%in%names(d))]=d[,-which(IDVAR%in%names(d))]^0.5
#     if(transf=="root4") d[,-which(IDVAR%in%names(d))]=d[,-which(IDVAR%in%names(d))]^0.25
#     if(transf=="log")   d[,-which(IDVAR%in%names(d))]=log(d[,-which(IDVAR%in%names(d))]+1e-6)
#     return(d)
#     
#   }
#   for(s in 1:length(Dat)) Dat[[s]]=transf.fn(d=Dat[[s]],transf=Store.transf[[s]])
#   
#   #Explore patterns
#   explore.fn=function(d)    
#   {
#     for(p in 1:length(Predictors))
#     {
#       x=d[,c(Predictors[p],names(d)[-which(IDVAR%in%names(d))])]  
#       if(!is.factor(x[,1]))
#       {
#         if(length(unique(x[,1]))>30)
#         {
#           x[,1]=cut(x[,1],10) 
#         }else x[,1]=as.factor(x[,1])
#         
#       }
#       agg=aggregate(formula(paste(".",names(x)[which(names(x)%in%IDVAR)],sep="~")), x, mean)
#       agg.sd=aggregate(formula(paste(".",names(x)[which(names(x)%in%IDVAR)],sep="~")), x, sd)
#       smart.par(n.plots=length(2:ncol(agg)),MAR=c(2,2,1,1),OMA=c(1,1.5,.1,.1),MGP=c(2.5,.7,0))
#       for(n in 2:ncol(agg))
#       {
#         xx=1: length(unique(agg[,1]))
#         agg.sd[,n][is.na(agg.sd[,n])]=0
#         plot(xx,agg[,n],pch=19,ylim=c(0,max(agg[,n]+agg.sd[,n])),ylab="",xlab="",xaxt='n',main=names(agg)[n],cex.axis=.85)
#         segments(xx,agg[,n]-agg.sd[,n],xx,agg[,n]+agg.sd[,n])
#         
#       }
#       mtext(names(Dat)[s],2,0.25,outer=T,las=3,cex=1.25)
#       mtext(Predictors[p],1,-0.5,outer=T,cex=1.25)
#       
#     }
#   }
#   for(s in 1:length(Dat)) explore.fn(d=Dat[[s]])
#   
#   #apply multivariate analyses   
#   nice.MDS.plot=function(MDS,pt.col,pt.bg,txt.col,txt.cex,PRED,col.pred,pred.cex)
#   {
#     sites <- scores(MDS, display = "sites")
#     spps  <- scores(MDS, display = "species")
#     xlim <- range(sites[,1], spps[,1])
#     ylim <- range(sites[,2], spps[,2])
#     plot(MDS, type = "n", xlim = xlim, ylim = ylim,ylab='',xlab='')
#     points(sites, bg = pt.bg, col=pt.col,pch = 21, cex = 1.2)
#     text(spps, labels = rownames(spps), col =txt.col, cex = txt.cex)
#     if(is.numeric(PRED)) PRED=round(PRED,1)
#     if(!is.null(PRED)) text(sites, labels = PRED, col =col.pred, cex = pred.cex)
#     legend("bottomright",paste("Stress=",round(MDS$stress,2)),bty='n')
#   }
#   interpret.multivariate=function(permanova)
#   {
#     PERMANOVA=as.data.frame(permanova$aov.tab)
#     Variation.explained=1-PERMANOVA$R2[match("Residuals",rownames(PERMANOVA))]
#     return(list(PERMANOVA=PERMANOVA,Variation.explained=Variation.explained))
#   }
#   pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = Best.indx, p.adjust.m ='bonferroni')
#   {
#     co = combn(unique(as.character(factors)),2)
#     pairs = c()
#     F.Model =c()
#     R2 = c()
#     p.value = c()
#     
#     
#     for(elem in 1:ncol(co))
#     {
#       if(sim.function == 'daisy'){
#         library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
#       } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
#       
#       ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
#       pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
#       F.Model =c(F.Model,ad$aov.tab[1,4]);
#       R2 = c(R2,ad$aov.tab[1,5]);
#       p.value = c(p.value,ad$aov.tab[1,6])
#     }
#     p.adjusted = p.adjust(p.value,method=p.adjust.m)
#     sig = c(rep('',length(p.adjusted)))
#     sig[p.adjusted <= 0.05] <-'.'
#     sig[p.adjusted <= 0.01] <-'*'
#     sig[p.adjusted <= 0.001] <-'**'
#     sig[p.adjusted <= 0.0001] <-'***'
#     
#     pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
#     #print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
#     return(pairw.res)
#     
#   } 
#   
#   fn.traditional.multivar=function(d,FORMULA)
#   {
#     require(vegan)
#     require(parallel)
#     require(MASS)
#     d.res.var=d[,-which(IDVAR%in%names(d))]
#     d.preds=d[,Predictors]
#     Cores=detectCores()-1
#     
#     #1. nMDS ordination to identify groups
#     Dis.indx=rankindex(d.preds, d.res.var, c("euc","man","gow","bray","jac","kul"))
#     Best.indx=names(rev(sort(Dis.indx)))[1]
#     MDS <- metaMDS(d.res.var, distance = Best.indx, autotransform = FALSE,trace = FALSE,parallel=Cores)
#     par(mfcol=c(1,1))
#     stressplot(MDS) #inspect stress with Shepard's plot
#     
#     #show species groups
#     nice.MDS.plot(MDS=MDS,pt.col=1,pt.bg=rgb(.1,.1,1,alpha=.2),txt.col=2,txt.cex=1,PRED=NULL,col.pred=1,pred.cex=.5)
#     
#     #show ordination for each predictor
#     for(p in 1:length(Predictors)) 
#     {
#       nice.MDS.plot(MDS=MDS,pt.col='transparent',pt.bg="transparent",txt.col="transparent",txt.cex=.9,
#                     PRED=d.preds[,match(Predictors[p],names(d.preds))],col.pred=1,pred.cex=.5)
#       mtext(Predictors[p],3,cex=1.5)
#     }
#     
#     
#     #2. Permanova
#     #note: adonis studies the differences in the group means
#     #     adonis is more robust than anosim
#     permanova <- adonis(FORMULA, d.preds,method = Best.indx,parallel=Cores)
#     
#     #2.1. interpret permanova
#     Permanova.table=interpret.multivariate(permanova)
#     
#     #2.2 Pairwise comparisons
#     Fctrs=d.preds[,sapply(d.preds,is.factor)]
#     if(is.data.frame(Fctrs))permanova.pairwise=vector('list',ncol(Fctrs)) else
#       permanova.pairwise=vector('list',1)
#     names(permanova.pairwise)=names(Fctrs)
#     for(ss in 1:length(permanova.pairwise))
#     {
#       if(is.data.frame(Fctrs))permanova.pairwise[[ss]] <- pairwise.adonis(d.res.var, Fctrs[,ss])else
#         permanova.pairwise[[ss]] <- pairwise.adonis(d.res.var, Fctrs)
#     }
#     
#     
#     
#     #3. Simper analysis to identify species that discriminate among groups
#     #note: if there are strong differences, then a few species should discreminate among groups
#     Simper=permanova.pairwise
#     for(ss in 1:length(Simper))Simper[[ss]] <- simper(d.res.var, Fctrs[,ss],parallel=Cores)
#     
#     
#     return(list(Best.indx=Best.indx,MDS=MDS,permanova=permanova,
#                 Permanova.table=Permanova.table,permanova.pairwise=permanova.pairwise,Simper=Simper))
#     
#   }
#   STore.multi.var.trad=Dat
#   for(s in 1:length(STore.multi.var.trad)) STore.multi.var.trad[[s]]=fn.traditional.multivar(d=Dat[[s]],
#                                                                                              FORMULA=Formula)
#   
# }





#Latitude

#note: only show species found to have a significant effect in glm
id.lat.count=which(unlist(lapply(PRED.lat.count,function(x) !(is.null(x)))))
id.lat.zero=which(unlist(lapply(PRED.lat.zero,function(x) !(is.null(x)))))
pooled.lat=unique(c(names(id.lat.count),names(id.lat.zero)))
fn.fig("Paper/Figure 3",1800,2400)
par(mfrow=c(length(pooled.lat),2),mai=c(.3,.4,.1,.15),oma=c(2,1,.1,1),las=1,mgp=c(.04,.6,0))
for(i in 1:length(pooled.lat))
{
  #count
  ii=match(pooled.lat[i],names(PRED.lat.count))
  if(!is.na(match(pooled.lat[i],names(id.lat.count)))) a=fun.plot.yr.pred(PRED.lat.count[[ii]],X="Mid.Lat",normalised="YES",REV="NO",n.seq=1,YLIM=NULL)
  if(is.na(match(pooled.lat[i],names(id.lat.count))))  plot.new()
  
  #zero part
  ii=match(pooled.lat[i],names(PRED.lat.zero))
  if(!is.na(match(pooled.lat[i],names(id.lat.zero)))) a=fun.plot.yr.pred(PRED.lat.zero[[ii]],X="Mid.Lat",normalised="NO",REV="NO",n.seq=1,YLIM=c(0,1))
  if(is.na(match(pooled.lat[i],names(id.lat.zero))))  plot.new()
  mtext(Tar.names[ii],4,line=.55,cex=1.2,las=3)
}
mtext("Relative CPUE (positive records)",2,outer=T,line=-0.75,cex=1.25,las=3)
mtext("Probability of zero catch",2,outer=T,line=-22.25,cex=1.25,las=3)
mtext("Latitude (�S)",1,outer=T,line=0.25,cex=1.25)
dev.off()



#Longitude
#note: only show species found to have a significant effect in glm
id.z.count=which(unlist(lapply(PRED.z.count,function(x) !(is.null(x)))))
id.z.zero=which(unlist(lapply(PRED.z.zero,function(x) !(is.null(x)))))
pooled.z=unique(c(names(id.z.count),names(id.z.zero)))
fn.fig("Paper/Figure 4",1800,2400)
par(mfrow=c(length(pooled.z),2),mai=c(.3,.4,.1,.15),oma=c(2,1,.1,1),las=1,mgp=c(.04,.6,0))
for(i in 1:length(pooled.z))
{
  #count
  ii=match(pooled.z[i],names(PRED.z.count))
  if(!is.na(match(pooled.z[i],names(id.z.count)))) a=fun.plot.yr.pred(PRED.z.count[[ii]],X="BOTDEPTH",normalised="YES",REV="NO",n.seq=20,YLIM=NULL)
  if(is.na(match(pooled.z[i],names(id.z.count))))  plot.new()
  
  #zero part
  ii=match(pooled.z[i],names(PRED.z.zero))
  if(!is.na(match(pooled.z[i],names(id.z.zero)))) a=fun.plot.yr.pred(PRED.z.zero[[ii]],X="BOTDEPTH",normalised="NO",REV="NO",n.seq=20,YLIM=c(0,1))
  if(is.na(match(pooled.z[i],names(id.z.zero))))  plot.new()
  mtext(Tar.names[ii],4,line=.55,cex=.85,las=3)
}
mtext("Relative CPUE (positive records)",2,outer=T,line=-0.75,cex=1.25,las=3)
mtext("Probability of zero catch",2,outer=T,line=-22.25,cex=1.25,las=3)
mtext("Depth (m)",1,outer=T,line=0.25,cex=1.25)
dev.off()







