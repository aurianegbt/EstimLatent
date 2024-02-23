# Required packages
library(ggplot2)
library(ggpubr)
library(scales)
library(dplyr)
# Custom Theme
source("~/Travail/00_Theme.R")

# Args
data = c("AB","G1","G2")
files ="Mult"
Nbr_ind = 15

# Load data
sim.list = lapply(data,FUN=function(d){read.csv(paste0("Model/IreneModel",files,"/Simulation/simulation_y",d,".txt"))})

Dataset <- data.frame()
for(d in 1:length(data)){
  Dataset <- Dataset %>% rbind(sim.list[[d]] %>% mutate(obsid = paste0("y",data[d]),.before = group) %>% rename(obs = paste0("y",data[d]))
  )
}

Dataset <- Dataset %>%
  filter(rep <= 5) %>%
  filter(if(Nbr_ind==50){group=="simulationGroup1"}else if(Nbr_ind==15){group=="simulationGroup2"}) %>%
  mutate(id=if(Nbr_ind==15){id-50}else{id}) %>%
  mutate(rep_id=paste0(rep,"_",id),.before = "id") %>%
  select(-group,-id,-rep)

if("AB" %in% data){
  DataSet = Dataset %>% filter(obsid=="yAB")
  moy = aggregate(DataSet[,"obs"],by = list(DataSet[,"time"]),FUN=function(x){mean(x,na.rm=TRUE)})

  med = aggregate(DataSet[,"obs"],by = list(DataSet[,"time"]),FUN=function(x){median(x,na.rm=TRUE)})

  sd = merge(aggregate(DataSet[,"obs",drop=FALSE],by = list(DataSet[,"time"]),FUN=function(x){sd(x,na.rm=TRUE)}), moy,by="Group.1")
  sd = as.data.frame(cbind(TIME=sd$Group.1,lcb = sd$x-1.96*sd$obs, ucb = sd$x +1.96*sd$obs,mean=sd$x))

  color = colorRampPalette(cbPalette)(length(unique(DataSet[,"rep_id"])))

  ggplot()+
    geom_point(data=DataSet,mapping=aes(x=time,y=obs,group=rep_id,color=rep_id))+
    geom_line(data=DataSet,mapping=aes(x=time,y=obs,group=rep_id,color=rep_id),alpha=0.4)+
    geom_line(data=sd,mapping=aes(x=TIME,y=mean),lwd=1)+
    geom_errorbar(data=sd,aes(x=TIME,ymin=lcb, ymax=ucb),lwd=0.8)+
    xlab("Time (in days)")+
    ylab("Antibody concentrations")+
    theme(axis.text = element_text(size=8))+
    theme(axis.title = element_text(size=14))+
    theme(plot.title = element_blank())+
    theme(plot.subtitle = element_text(hjust=0.5))+
    theme(strip.text = element_text(size = 16))+
    scale_color_manual(values=color)+
    theme(legend.position = "none")

  if(!dir.exists(paste0("~/Travail/Presentation/Plot/PhD_reflective_Irene/",files))){dir.create(paste0("~/Travail/Presentation/Plot/PhD_reflective_Irene/",files))}
  ggsave(filename = paste0("~/Travail/Presentation/Plot/PhD_reflective_Irene/",files,"/Ab",Nbr_ind,"_errorBar.png"))
}

for(d in data[stringr::str_detect(data,"G")]){
  DataSet <- Dataset %>% filter(obsid==paste0("y",d))

  moy = aggregate(DataSet[,"obs"],by = list(DataSet[,"time"]),FUN=function(x){mean(x,na.rm=TRUE)})

  med = aggregate(DataSet[,"obs"],by = list(DataSet[,"time"]),FUN=function(x){median(x,na.rm=TRUE)})

  sd = merge(aggregate(DataSet[,"obs",drop=FALSE],by = list(DataSet[,"time"]),FUN=function(x){sd(x,na.rm=TRUE)}), moy,by="Group.1")

  sd = as.data.frame(cbind(TIME=sd$Group.1,lcb = sd$x-1.96*sd$obs, ucb = sd$x +1.96*sd$obs,mean=sd$x))


  color = colorRampPalette(cbPalette)(length(unique(DataSet[,"rep_id"])))


  ggplot()+
    geom_point(data=DataSet,mapping=aes(x=time,y=obs,group=rep_id,color=rep_id),alpha=0.3)+
    geom_line(data=DataSet,mapping=aes(x=time,y=obs,group=rep_id,color=rep_id),alpha=0.4)+
    geom_line(data=sd,mapping=aes(x=TIME,y=mean),lwd=1)+
    geom_errorbar(data=sd,aes(x=TIME,ymin=lcb, ymax=ucb),lwd=0.8)+
    xlab("Time (in days)")+
    ylab(paste0("Genes ",d," expression"))+
    theme(axis.text = element_text(size=8))+
    theme(axis.title = element_text(size=14))+
    theme(plot.title = element_blank())+
    theme(plot.subtitle = element_text(hjust=0.5))+
    theme(strip.text = element_text(size = 16))+
    scale_color_manual(values=color)+
    theme(legend.position = "none")
  ggsave(filename = paste0("~/Travail/Presentation/Plot/PhD_reflective_Irene/",files,"/",d,Nbr_ind,"_errorBar.png"))
}



