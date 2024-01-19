# Required packages
library(ggplot2)
library(ggpubr)
library(scales)
library(dplyr)
# Custom Theme
source("~/Travail/00_Theme.R")
# Load data
DataSet = read.csv("IreneModel2/Simulation/simulatedData_Ab_G1.csv") %>%
  filter(obsid=="yAB")

moy = aggregate(DataSet[,"obs"],by = list(DataSet[,"TIME"]),FUN=function(x){mean(x,na.rm=TRUE)})

med = aggregate(DataSet[,"obs"],by = list(DataSet[,"TIME"]),FUN=function(x){median(x,na.rm=TRUE)})

sd = merge(aggregate(DataSet[,"obs",drop=FALSE],by = list(DataSet[,"TIME"]),FUN=function(x){sd(x,na.rm=TRUE)}), moy,by="Group.1")

sd = as.data.frame(cbind(TIME=sd$Group.1,lcb = sd$x-1.96*sd$obs, ucb = sd$x +1.96*sd$obs,mean=sd$x))

color = colorRampPalette(cbPalette)(length(unique(DataSet[DataSet$rep==1,"ID"])))

ggplot()+
  geom_point(data=DataSet[DataSet$rep==1,],mapping=aes(x=TIME,y=obs,group=UID,color=UID))+
  geom_line(data=DataSet[DataSet$rep==1,],mapping=aes(x=TIME,y=obs,group=UID,color=UID),alpha=0.4)+
  geom_line(data=sd,mapping=aes(x=TIME,y=mean),lwd=1)+
  # geom_errorbar(data=sd,aes(x=TIME,ymin=lcb, ymax=ucb),lwd=0.8)+
  xlab("Time (in days)")+
  ylab("Antibody concentrations")+
  theme(axis.text = element_text(size=8))+
  theme(axis.title = element_text(size=14))+
  theme(plot.title = element_blank())+
  theme(plot.subtitle = element_text(hjust=0.5))+
  theme(strip.text = element_text(size = 16))+
  scale_color_manual(values=color)+
  theme(legend.position = "none")
########################################################################
DataSet = read.csv("IreneModel2/Simulation/simulatedData_Ab_G1.csv") %>%
  filter(obsid=="yG1")

moy = aggregate(DataSet[,"obs"],by = list(DataSet[,"TIME"]),FUN=function(x){mean(x,na.rm=TRUE)})

med = aggregate(DataSet[,"obs"],by = list(DataSet[,"TIME"]),FUN=function(x){median(x,na.rm=TRUE)})

sd = merge(aggregate(DataSet[,"obs",drop=FALSE],by = list(DataSet[,"TIME"]),FUN=function(x){sd(x,na.rm=TRUE)}), moy,by="Group.1")

sd = as.data.frame(cbind(TIME=sd$Group.1,lcb = sd$x-1.96*sd$obs, ucb = sd$x +1.96*sd$obs,mean=sd$x))

color = colorRampPalette(cbPalette)(length(unique(DataSet[DataSet$rep==1,"ID"])))

ggplot()+
  geom_point(data=DataSet[DataSet$rep==1,],mapping=aes(x=TIME,y=obs,group=UID,color=UID),alpha=0.3)+
  geom_line(data=DataSet[DataSet$rep==1,],mapping=aes(x=TIME,y=obs,group=UID,color=UID),alpha=0.4)+
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

color = colorRampPalette(cbPalette)(5)

ggplot()+
  geom_point(data=DataSet[DataSet$rep==1 & DataSet$ID <=5,],mapping=aes(x=TIME,y=obs,group=UID,color=UID),alpha=0.3)+
  geom_line(data=DataSet[DataSet$rep==1 & DataSet$ID <=5,],mapping=aes(x=TIME,y=obs,group=UID,color=UID))+
  xlab("Time (in days)")+
  ylab("Antibody concentrations")+
  theme(axis.text = element_text(size=8))+
  theme(axis.title = element_text(size=14))+
  theme(plot.title = element_blank())+
  theme(plot.subtitle = element_text(hjust=0.5))+
  theme(strip.text = element_text(size = 16))+
  scale_color_manual(values=color)+
  theme(legend.position = "none")
