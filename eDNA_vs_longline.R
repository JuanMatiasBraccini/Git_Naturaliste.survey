# DATA SECTION ------------------------------------------------------------

#Define user
User="Matias"

fn.user=function(x1,x2)paste(x1,Sys.getenv("USERNAME"),x2,sep='/')
if(!exists('handl_OneDrive')) source(fn.user(x1='C:/Users',
                                             x2='OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R'))
library(tidyverse)
library(readxl)

#Sharks data base 
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R"))

#eDNA samples collected
eDNA.samples=read.csv(handl_OneDrive('eDNA/eDNA_Samples.csv'))

#eDNA results
eDNA.results_taxonomic_matrix <- read_excel(handl_OneDrive('eDNA/Matias_shark_eDNA_taxonomic_matrix.xlsx'), sheet = 'unsubsampled_matrix',skip = 0)
eDNA.results_Common.name=read.csv(handl_OneDrive('eDNA/Common.name_eDNA results.csv'))

# Manipulate data ------------------------------------------------------------
eDNA.samples=eDNA.samples%>%
              mutate(Date=as.Date(Date,"%d/%m/%Y"))

eDNA.dat_longline=DATA%>%
          filter(BOAT=='NAT' & Lat.round>(-27) & SHEET_NO%in%unique(eDNA.samples$Sheet.Number))%>%
          mutate(SPECIES=ifelse(SPECIES=='GH' & SHEET_NO=='I00842','HG',SPECIES))%>%
          dplyr::select(SHEET_NO,date,Mid.Lat,Mid.Long,SOAK.TIME,BOTDEPTH,SPECIES,CAAB_code,
                        CAES_Code,COMMON_NAME,SCIENTIFIC_NAME,SEX,FL,TL,BAG_NO,Number,RECORDER)

working='online'  #remove all this when back in Oz
if(working=='offline')
{
  #write.csv(eDNA.samples,'C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Desktop/eDNA data/eDNA.samples.csv',row.names = F)
  #write.csv(eDNA.dat_longline,'C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Desktop/eDNA data/eDNA.dat_longline.csv',row.names = F)
  #write.csv(eDNA.results_taxonomic_matrix,'C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Desktop/eDNA data/eDNA.results_taxonomic_matrix.csv',row.names = F)
  #write.csv(eDNA.results_Common.name,'C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Desktop/eDNA data/eDNA.results_Common.name.csv',row.names = F)

  eDNA.samples=read.csv('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Desktop/eDNA data/eDNA.samples.csv')
  eDNA.dat_longline=read.csv('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Desktop/eDNA data/eDNA.dat_longline.csv')
  eDNA.results_taxonomic_matrix=read.csv('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Desktop/eDNA data/eDNA.results_taxonomic_matrix.csv')
  eDNA.results_Common.name=read.csv('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Desktop/eDNA data/eDNA.results_Common.name.csv')
}

Meta.data=eDNA.samples%>%left_join(eDNA.dat_longline%>%
                                     distinct(SHEET_NO,Mid.Lat,Mid.Long),
                                   by=c('Sheet.Number'='SHEET_NO'))
setwd(handl_OneDrive('Analyses/Surveys/Naturaliste_longline/outputs/eDNA_vs_longline'))
write.csv(Meta.data,"Meta.data_eDNA_samples.csv",row.names = F)
write.csv(eDNA.dat_longline,"Longline_samples.csv",row.names = F)


eDNA.dat_longline=eDNA.dat_longline%>%
          left_join(eDNA.samples%>%
                      dplyr::select(Sample.Time,Filter.Time,Water.Depth,Sample.Depth,
                                                 Sample.ID,Comments,Sheet.Number)%>%
                      filter(Sheet.Number%in%unique(eDNA.dat_longline$SHEET_NO))%>%
                      filter(grepl("B",Sample.ID)),  #ACA: why selection "B" only??? loosing Top and 2023 altogether
                    by=c('SHEET_NO'='Sheet.Number'))

Matrix=eDNA.dat_longline%>%
  group_by(Sample.ID,COMMON_NAME)%>%
  tally()%>%
  spread(COMMON_NAME,n,fill=0)%>%
  data.frame
write.csv(Matrix,"Longline_matrix.csv",row.names = F)

eDNA.results_taxonomic_matrix=eDNA.results_taxonomic_matrix%>%
                                gather(Sample.ID,n,-c(domain, phylum, class, order, family, genus, species, taxa))%>%
  left_join(eDNA.results_Common.name,by=c('taxa'))


#combine eDNA results with longline
#note: combine replicates (_A & _B)
  #Common name
full_join(eDNA.results_taxonomic_matrix%>%
                     filter(grepl('BOTTOM',Sample.ID))%>%
                     mutate(Sample.ID=str_remove(Sample.ID,"(=?_).+?"),
                            Sample.ID=str_replace(Sample.ID, "St", "Stn"),
                            Sample.ID=str_replace(Sample.ID, "BOTTOM", "B"),
                            Sample.ID=str_replace(Sample.ID, "L", "_"),
                            dummy=sub("_.*", "", Sample.ID),
                            dummy2=str_extract(sub(".*_", "", Sample.ID), "[[:digit:]]+"),
                            Sample.ID=paste(dummy,dummy2,sep='_B'))%>%
                     dplyr::select(-c(dummy,dummy2))%>%
                     group_by(COMMON_NAME,Sample.ID)%>%
                     summarise(n=sum(n))%>%
                     mutate(n=ifelse(n>0,1,0))%>%
                     rename(eDNA=n)%>%
                     filter(eDNA>0),
                   eDNA.dat_longline%>%group_by(COMMON_NAME,Sample.ID)%>%
                     tally()%>%
                     mutate(n=ifelse(n>0,1,0))%>%
                     rename(longline=n),
                   by=c('COMMON_NAME','Sample.ID'))%>%
                mutate(eDNA=ifelse(is.na(eDNA),0,eDNA),
                       longline=ifelse(is.na(longline),0,longline))%>%
  gather(Method,n,-c(COMMON_NAME,Sample.ID))%>%
  mutate(n=ifelse(n==0,'Absence','Presence'))%>%
  ggplot(aes(x=Method, y=COMMON_NAME, fill = factor(n))) +
  geom_point(size = 3, shape = 21)+
  facet_wrap(~Sample.ID,scales = 'free_y')+
  theme(legend.position = 'top',
        legend.title = element_blank())+
  ylab('')+xlab('')
ggsave("eDNA_vs_longline.tiff",width = 20,height = 10,compression = "lzw")

  #Scientific name
full_join(eDNA.results_taxonomic_matrix%>%
            filter(grepl('BOTTOM',Sample.ID))%>%
            rename(SCIENTIFIC_NAME=taxa)%>%
            mutate(Sample.ID=str_remove(Sample.ID,"(=?_).+?"),
                   Sample.ID=str_replace(Sample.ID, "St", "Stn"),
                   Sample.ID=str_replace(Sample.ID, "BOTTOM", "B"),
                   Sample.ID=str_replace(Sample.ID, "L", "_"),
                   dummy=sub("_.*", "", Sample.ID),
                   dummy2=str_extract(sub(".*_", "", Sample.ID), "[[:digit:]]+"),
                   Sample.ID=paste(dummy,dummy2,sep='_B'))%>%
            dplyr::select(-c(dummy,dummy2))%>%
            group_by(SCIENTIFIC_NAME,Sample.ID)%>%
            summarise(n=sum(n))%>%
            mutate(n=ifelse(n>0,1,0))%>%
            rename(eDNA=n)%>%
            filter(eDNA>0),
          eDNA.dat_longline%>%group_by(SCIENTIFIC_NAME,Sample.ID)%>%
            tally()%>%
            mutate(n=ifelse(n>0,1,0))%>%
            rename(longline=n),
          by=c('SCIENTIFIC_NAME','Sample.ID'))%>%
  mutate(eDNA=ifelse(is.na(eDNA),0,eDNA),
         longline=ifelse(is.na(longline),0,longline))%>%
  gather(Method,n,-c(SCIENTIFIC_NAME,Sample.ID))%>%
  mutate(n=ifelse(n==0,'Absence','Presence'))%>%
  ggplot(aes(x=Method, y=SCIENTIFIC_NAME, fill = factor(n))) +
  geom_point(size = 3, shape = 21)+
  facet_wrap(~Sample.ID,scales = 'free_y')+
  theme(legend.position = 'top',
        legend.title = element_blank())+
  ylab('')+xlab('')
ggsave("eDNA_vs_longline_scien.name.tiff",width = 20,height = 10,compression = "lzw")



#combine eDNA bottom and surface
#note: combine replicates (_A & _B)

  #Common name
full_join(eDNA.results_taxonomic_matrix%>%
                     filter(grepl('BOTTOM',Sample.ID))%>%
                     mutate(Sample.ID=str_remove(Sample.ID,"(=?_).+?"),
                            Sample.ID=str_remove(Sample.ID, "BOTTOM"))%>%
                     group_by(COMMON_NAME,Sample.ID)%>%
                     summarise(n=sum(n))%>%
                     mutate(n=ifelse(n>0,1,0))%>%
                     rename(eDNA.bottom=n)%>%
                     filter(eDNA.bottom>0),
                   eDNA.results_taxonomic_matrix%>%
                     filter(grepl('TOP',Sample.ID))%>%
                     mutate(Sample.ID=str_remove(Sample.ID,"(=?_).+?"),
                            Sample.ID=str_remove(Sample.ID,"TOP"))%>%
                     group_by(COMMON_NAME,Sample.ID)%>%
                     summarise(n=sum(n))%>%
                     mutate(n=ifelse(n>0,1,0))%>%
                     rename(eDNA.surface=n)%>%
                     filter(eDNA.surface>0),
                   by=c('COMMON_NAME','Sample.ID'))%>%
  mutate(eDNA.bottom=ifelse(is.na(eDNA.bottom),0,eDNA.bottom),
         eDNA.surface=ifelse(is.na(eDNA.surface),0,eDNA.surface))%>%
  gather(Method,n,-c(COMMON_NAME,Sample.ID))%>%
  mutate(n=ifelse(n==0,'Absence','Presence'))%>%
  ggplot(aes(x=Method, y=COMMON_NAME, fill = factor(n))) +
  geom_point(size = 3, shape = 21)+
  facet_wrap(~Sample.ID,scales = 'free_y')+
  theme(legend.position = 'top',
        legend.title = element_blank())+
  ylab('')+xlab('')
ggsave("eDNA.bottom_vs_eDNA.surface.tiff",width = 20,height = 10,compression = "lzw")

  #Scientific name
full_join(eDNA.results_taxonomic_matrix%>%
            filter(grepl('BOTTOM',Sample.ID))%>%
            rename(SCIENTIFIC_NAME=taxa)%>%
            mutate(Sample.ID=str_remove(Sample.ID,"(=?_).+?"),
                   Sample.ID=str_remove(Sample.ID, "BOTTOM"))%>%
            group_by(SCIENTIFIC_NAME,Sample.ID)%>%
            summarise(n=sum(n))%>%
            mutate(n=ifelse(n>0,1,0))%>%
            rename(eDNA.bottom=n)%>%
            filter(eDNA.bottom>0),
          eDNA.results_taxonomic_matrix%>%
            filter(grepl('TOP',Sample.ID))%>%
            rename(SCIENTIFIC_NAME=taxa)%>%
            mutate(Sample.ID=str_remove(Sample.ID,"(=?_).+?"),
                   Sample.ID=str_remove(Sample.ID,"TOP"))%>%
            group_by(SCIENTIFIC_NAME,Sample.ID)%>%
            summarise(n=sum(n))%>%
            mutate(n=ifelse(n>0,1,0))%>%
            rename(eDNA.surface=n)%>%
            filter(eDNA.surface>0),
          by=c('SCIENTIFIC_NAME','Sample.ID'))%>%
  mutate(eDNA.bottom=ifelse(is.na(eDNA.bottom),0,eDNA.bottom),
         eDNA.surface=ifelse(is.na(eDNA.surface),0,eDNA.surface))%>%
  gather(Method,n,-c(SCIENTIFIC_NAME,Sample.ID))%>%
  mutate(n=ifelse(n==0,'Absence','Presence'))%>%
  ggplot(aes(x=Method, y=SCIENTIFIC_NAME, fill = factor(n))) +
  geom_point(size = 3, shape = 21)+
  facet_wrap(~Sample.ID,scales = 'free_y')+
  theme(legend.position = 'top',
        legend.title = element_blank())+
  ylab('')+xlab('')
ggsave("eDNA.bottom_vs_eDNA.surface_scien.name.tiff",width = 20,height = 10,compression = "lzw")
