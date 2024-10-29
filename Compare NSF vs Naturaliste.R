library(ggpubr)

#Define user
User="Matias"

if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')



dat=DATA%>%
      filter(Mid.Lat>(-26) & Method=='LL')%>%
      filter(!BOAT%in%c('FLIN','RV GANNET','RV SNIPE 2','HAM','HOU','RV BREAKSEA'))%>%
      mutate(Vessel=ifelse(BOAT=='NAT','Naturaliste','Commercial'))

p1=dat%>%
  ggplot( aes(x=SOAK.TIME, fill=Vessel)) +
  geom_histogram(color="#e9ecef", alpha=0.6,  position = 'identity')+
  xlab('Soak time (hours)')

p2=dat%>%
  ggplot( aes(x=N.hooks, fill=Vessel)) +
  geom_histogram(color="#e9ecef", alpha=0.6,  position = 'identity')+
  xlab('Number of hooks')

p3=dat%>%
  ggplot( aes(x=BOTDEPTH, fill=Vessel)) +
  geom_histogram(color="#e9ecef", alpha=0.6,  position = 'identity')+
  xlab('Bottom depth (m)')+xlim(0,500)

p4=dat%>%
  ggplot(aes(Mid.Long,Mid.Lat,color=Vessel ))+
  geom_point()

ggarrange(plotlist=list(p1,p2,p3,p4), ncol = 2, nrow = 2,common.legend=TRUE)
ggsave('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/Surveys/Naturaliste_longline/outputs/Longline_NSF.Commercial.vs.Naturaliste.tiff',
       width = 8,height = 8,compression = "lzw")



