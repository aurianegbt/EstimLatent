library(dplyr)
library(flextable)
## DATA INFORMATION
data = "Ab_G1_G3"
type="default"
files="Mult"

load(paste0("Results/Results",files,"_",data,"/bootstrap",type,"/estim_1.RData"))
nBoot = length(resBoot)

# nrep = length(list.files(paste0("Results/Results",files,"_",data,"/bootstrap",type)))
nrep = 100
precision = abs(floor(log10(1/nrep)))

# Gather results

EstimatedStandardDeviation = data.frame() # tableau des standard errors estimées sur le réplicat par monolix
EstimatedPopulationParameters = data.frame() # tableau des paramètres estimés sur le réplicat par monolix

EstimatedStandardDeviationBOOT = data.frame() # tableau des écart-type estimées par bootstrap (sd des paramètres estimés sur chaque boot du réplicat)
EstimatedPopulationParametersBOOT = data.frame() # tableau du paramètre estimées par bootstrap (mean des paramètres estimées sur chaque boot du réplicat)

load(paste0("Results/Results",files,"_",data,"/",type,"/PopulationParameters_theo.RData"))

for(f in 1:nrep){
  if(exists("resBoot"))
    rm(resBoot)
  if(exists("res"))
    rm(res)
  suppressWarnings(tryCatch({
    pathToRes = paste0("Results/Results",files,"_",data,"/bootstrap",type,"/estim_",f,".RData")
    load(pathToRes)
    if(length(resBoot)!=500)
      cat("file n°",f," : ",length(resBoot),"\n")

      EstimatedStandardDeviation = rbind(EstimatedStandardDeviation, t(res[,"StandardError",drop=FALSE]))
      EstimatedPopulationParameters = rbind(EstimatedPopulationParameters,t(res[,"Parameters",drop=FALSE]))

    Est = Reduce(cbind,lapply(resBoot,function(dt){dt[,"Parameters",drop=FALSE]}))

    EstimatedStandardDeviationBOOT = rbind(EstimatedStandardDeviationBOOT, t(as.data.frame(apply(Est,1,sd))))
    EstimatedPopulationParametersBOOT = rbind(EstimatedPopulationParametersBOOT,t(as.data.frame(rowMeans(Est))))

    rownames(EstimatedPopulationParametersBOOT)[nrow(EstimatedPopulationParametersBOOT)] <-
      rownames(EstimatedStandardDeviationBOOT)[nrow(EstimatedStandardDeviationBOOT)] <-
      rownames(EstimatedPopulationParameters)[nrow(EstimatedPopulationParameters)] <-
      rownames(EstimatedStandardDeviation)[nrow(EstimatedStandardDeviation)] <-
      paste0("estim_",f)

  },error=function(e){cat(f,",")}))
}

EstimatedStandardDeviation <- t(EstimatedStandardDeviation)[names(PopulationParameters),]
EstimatedPopulationParameters <- t(EstimatedPopulationParameters)[names(PopulationParameters),]

EstimatedStandardDeviationBOOT <- t(EstimatedStandardDeviationBOOT)[names(PopulationParameters),]
EstimatedPopulationParametersBOOT <- t(EstimatedPopulationParametersBOOT)[names(PopulationParameters),]


# Simple stats

stats_estim = data.frame(target_values = setNames(as.numeric(PopulationParameters),names(PopulationParameters)))

stats_estim <- cbind(stats_estim,
                     mean=rowMeans(EstimatedPopulationParameters),
                     meanBoot = rowMeans(EstimatedPopulationParametersBOOT),
                     bias=as.numeric(((rowMeans(EstimatedPopulationParameters)-stats_estim[,'target_values']))/(PopulationParameters)),
                     biasBoot = as.numeric(((rowMeans(EstimatedPopulationParametersBOOT)-stats_estim[,'target_values']))/(PopulationParameters)),
                     sd_emp = apply(EstimatedPopulationParameters,1,sd),
                     sd_hess = rowMeans(apply(EstimatedStandardDeviation,MARGIN=2,FUN=as.numeric),na.rm = TRUE),
                     sd_boot = rowMeans(apply(EstimatedStandardDeviationBOOT,MARGIN=2,FUN=as.numeric),na.rm = TRUE))

countNA = sapply(apply(EstimatedStandardDeviation,MARGIN=1,FUN=function(x){sum(x=="NaN")}),FUN=function(x){if(x==0){"."}else{x}})
if(!all(countNA==".")){
  stats_estim <- cbind(stats_estim,countNA)[,c("target_values","mean","meanBoot","bias","biasBoot","sd_emp","sd_hess","countNA","sd_boot")]
}


# COVER
coverHess = rep(0,length(PopulationParameters)) # on calcule un taux de couverture pour chaque paramètres
n = ncol(EstimatedPopulationParameters) # nrep
nCount = rep(n,length(PopulationParameters)) # on compte le nombre de non NA pour chaque paramètres
for (j in 1:n){ # pour chaque réplicat
  sdhat = as.numeric(EstimatedStandardDeviation[,j]) # on regarde le vecteur contenant les se estimées pour chaque paramètre dans ce réplicat
  for(par in 1:length(sdhat)){ # pour chaque paramètre estimé dans ce réplicat
    if(sdhat[par]=="NaN"){ # si le sd n'est pas calculé (error IF), alors on réduit de 1 le nombre de valeur estimées au total
      nCount[par] <- nCount[par] -1
    }else{ # sinon, on regarde si le paramètre estimé pour ce paramètre, ce réplicat, est dans l'IC à 95%
      lw = EstimatedPopulationParameters[par,j] - 1.96*sdhat[par]
      up = EstimatedPopulationParameters[par,j] + 1.96*sdhat[par]
      bool = as.numeric(PopulationParameters[par] > lw & PopulationParameters[par] < up)
      coverHess[par] = coverHess[par] + bool
    }
  }
}
coverHess <- coverHess/nCount

coverBoot = rep(0,length(PopulationParameters))
n = ncol(EstimatedPopulationParameters)
nCount = rep(n,length(PopulationParameters))
for (j in 1:n){
  sdhat = as.numeric(EstimatedStandardDeviationBOOT[,j])
  for(par in 1:length(sdhat)){
    if(sdhat[par]=="NaN"){
      nCount[par] <- nCount[par] -1
    }else{
      lw = EstimatedPopulationParameters[par,j] - 1.96*sdhat[par]
      up = EstimatedPopulationParameters[par,j] + 1.96*sdhat[par]
      bool = as.numeric(PopulationParameters[par] > lw & PopulationParameters[par] < up)
      coverBoot[par] = coverBoot[par] + bool
    }
  }
}
coverBoot <- coverBoot/nCount

stats_estim <- cbind(stats_estim,cover_hess=coverHess,cover_boot = coverBoot)

if("countNA" %in% colnames(stats_estim)){
  colnames(stats_estim) <- c("Target Value","Mean of estimation","Mean of estimation by bootstrap","Relative Bias","Relative Bias by bootstrap","Empirical Standard Deviation","Estimated Standard Deviation","(NA)","Estimated Standard Deviation by bootstrap","Estimated Cover","Estimated Cover by bootstrap")
}else{
  colnames(stats_estim) <- c("Target Value","Mean of estimation","Mean of estimation by bootstrap","Relative Bias","Relative Bias by bootstrap","Empirical Standard Deviation","Estimated Standard Deviation","Estimated Standard Deviation by bootstrap","Estimated Cover","Estimated Cover by bootstrap")
}


if(files==""){
  if(data=="Ab"){
    tex2_names = c(
      " {f_M}_1 ",
      " \\delta_S ",
      " \\theta ",
      " \\delta_{Ab} ",
      " \\omega_{{f_M}_1} ",
      " \\omega_{\\delta_S} ",
      " \\omega_\\theta ",
      " \\sigma_{AB} "
    )
  }else if(data=="Ab_G1"){
    tex2_names = c(
      " {f_M}_1 ",
      " \\delta_S ",
      " \\theta ",
      " \\delta_{Ab} ",
      " \\alpha_1 ",
      " \\omega_{{f_M}_1} ",
      " \\omega_{\\delta_S} ",
      " \\omega_\\theta ",
      " \\sigma_{AB} ",
      " \\sigma_{G1} "
    )
  }else if(data=="Ab_G2"){
    tex2_names = c(
      " {f_M}_1 ",
      " \\delta_S ",
      " \\theta ",
      " \\delta_{Ab} ",
      " \\alpha_2 ",
      " \\omega_{{f_M}_1} ",
      " \\omega_{\\delta_S} ",
      " \\omega_\\theta ",
      " \\sigma_{AB} ",
      " \\sigma_{G2} "
    )
  }else{
    tex2_names = c(
      " {f_M}_1 ",
      " \\delta_S ",
      " \\theta ",
      " \\delta_{Ab} ",
      " \\alpha_1 ",
      " \\alpha_2 ",
      " \\omega_{{f_M}_1} ",
      " \\omega_{\\delta_S} ",
      " \\omega_\\theta ",
      " \\sigma_{AB} ",
      " \\sigma_{G1} ",
      " \\sigma_{G2} "
    )
  }
  if(type=="fixed_delta"){
    tex2_names <- setdiff(tex2_names,c(" \\delta_{Ab} "," \\delta_S "," \\omega_{\\delta_S} "))
  }else if(type=="default"){
    tex2_names <- setdiff(tex2_names,c(" \\omega_{\\delta_S} "))
  }else if(type=="only_S"){

    tex2_names <- setdiff(tex2_names,c(" \\delta_{Ab} "," \\omega_{\\delta_S} "))
  }
}else if(files=="2"){
  tex2_names = c(
    " {f_M}_2 ",
    " \\delta_S ",
    " \\theta_1 ",
    " \\theta_2 ",
    " \\alpha_1 ",
    " \\omega_{{f_M}_1} ",
    " \\omega_{{f_M}_2} ",
    " \\omega_{\\theta_1} ",
    " \\omega_{\\theta_2} ",
    " \\sigma_{AB} ",
    " \\sigma_{G1} "
  )
}else if(files=="Mult"){
  if(data=="Ab"){
    if(type=="default"){
      tex2_names = c(
        " {f_M}_1 ",
        " \\delta_S ",
        " \\theta ",
        " \\delta_{Ab} ",
        " \\omega_{{f_M}_1} ",
        " \\omega_\\theta ",
        " \\sigma_{AB} "
      )
    }else if(type=="poor"){
      tex2_names = c(
        " \\delta_S ",
        " \\theta ",
        " \\omega_{{f_M}_1} ",
        " \\omega_\\theta ",
        " \\sigma_{AB} "
      )
    }else if(type=="only_S"){
      tex2_names = c(
        " {f_M}_1 ",
        " \\delta_S ",
        " \\theta ",
        " \\omega_{{f_M}_1} ",
        " \\omega_\\theta ",
        " \\sigma_{AB} "
      )
    }
  }else if(data=="Ab_G1_G3"){
    if(type=="default"){
      tex2_names = c(
        " {f_M}_1 ",
        " \\delta_S ",
        " \\theta ",
        " \\delta_{Ab} ",
        " \\alpha_1 ",
        " \\alpha_3 ",
        " \\omega_{{f_M}_1} ",
        " \\omega_\\theta ",
        " \\sigma_{AB} ",
        " \\sigma_{G1} ",
        " \\sigma_{G3} "
      )
    }else if(type=="only_S"){
      tex2_names = c(
        " {f_M}_1 ",
        " \\delta_S ",
        " \\theta ",
        " \\alpha_1 ",
        " \\alpha_3 ",
        " \\omega_{{f_M}_1} ",
        " \\omega_\\theta ",
        " \\sigma_{AB} ",
        " \\sigma_{G1} ",
        " \\sigma_{G3} "
      )
    }
  }else if(data=="Ab_G1"){
    if(type=="default"){
      tex2_names = c(
        " {f_M}_1 ",
        " \\delta_S ",
        " \\theta ",
        " \\delta_{Ab} ",
        " \\alpha_1 ",
        " \\omega_{{f_M}_1} ",
        " \\omega_\\theta ",
        " \\sigma_{AB} ",
        " \\sigma_{G1} "
      )
    }else if(type=="only_S"){
      tex2_names = c(
        " {f_M}_1 ",
        " \\delta_S ",
        " \\theta ",
        " \\alpha_1 ",
        " \\omega_{{f_M}_1} ",
        " \\omega_\\theta ",
        " \\sigma_{AB} ",
        " \\sigma_{G1} "
      )
    }
  }else if(type=="default"){
    tex2_names = c(
      " {f_M}_1 ",
      " \\delta_S ",
      " \\theta ",
      " \\delta_{Ab} ",
      " \\alpha_1 ",
      " \\alpha_2 ",
      " \\alpha_3 ",
      " \\alpha_4 ",
      " \\alpha_5 ",
      " \\omega_{{f_M}_1} ",
      " \\omega_\\theta ",
      " \\sigma_{AB} ",
      " \\sigma_{G1} ",
      " \\sigma_{G2} ",
      " \\sigma_{G3} ",
      " \\sigma_{G4} ",
      " \\sigma_{G5} "
    )
  }else if(type=="only_S"){
    tex2_names = c(
      " {f_M}_1 ",
      " \\delta_S ",
      " \\theta ",
      " \\alpha_1 ",
      " \\alpha_2 ",
      " \\alpha_3 ",
      " \\alpha_4 ",
      " \\alpha_5 ",
      " \\omega_{{f_M}_1} ",
      " \\omega_\\theta ",
      " \\sigma_{AB} ",
      " \\sigma_{G1} ",
      " \\sigma_{G2} ",
      " \\sigma_{G3} ",
      " \\sigma_{G4} ",
      " \\sigma_{G5} "
    )
  }

}else if(files=="Mu"){
  if(type=="only_S"){
    tex2_names = c(
      " {f_M}_1 ",
      " \\delta_S ",
      " \\theta ",
      " \\mu_1 ",
      " \\mu_2 ",
      " \\mu_3 ",
      " \\mu_4 ",
      " \\mu_5 ",
      " \\alpha_1 ",
      " \\alpha_2 ",
      " \\alpha_3 ",
      " \\alpha_4 ",
      " \\alpha_5 ",
      " \\omega_{{f_M}_1} ",
      " \\omega_\\theta ",
      " \\sigma_{AB} ",
      " \\sigma_{G1} ",
      " \\sigma_{G2} ",
      " \\sigma_{G3} ",
      " \\sigma_{G4} ",
      " \\sigma_{G5} "
    )
  }
}else if(files=="Mu_i"){
  if(type=="only_S"){
    tex2_names = c(
      " {f_M}_1 ",
      " \\delta_S ",
      " \\theta ",
      " \\mu_1 ",
      " \\mu_2 ",
      " \\mu_3 ",
      " \\mu_4 ",
      " \\mu_5 ",
      " \\alpha_1 ",
      " \\alpha_2 ",
      " \\alpha_3 ",
      " \\alpha_4 ",
      " \\alpha_5 ",
      " \\omega_{{f_M}_1} ",
      " \\omega_\\theta ",
      " \\omega_{mu_1} ",
      " \\omega_{mu_2} ",
      " \\omega_{mu_3} ",
      " \\omega_{mu_4} ",
      " \\omega_{mu_5} ",
      " \\sigma_{AB} ",
      " \\sigma_{G1} ",
      " \\sigma_{G2} ",
      " \\sigma_{G3} ",
      " \\sigma_{G4} ",
      " \\sigma_{G5} "
    )
  }
}
options(scipen=999)

stats_estim[,!(colnames(stats_estim) %in% c("Target Value","(NA)"))] <- apply(stats_estim[,!(colnames(stats_estim) %in% c("Target Value","(NA)"))],2, FUN=function(x){as.character(round(x,digits=precision+1))})


cat("\n \n ----------- LATEX CODE FOR RESULTS -----------\n \n")
stats_estimLTX <- huxtable::huxtable(cbind(Parameter = paste0("$",tex2_names,"$"),stats_estim),add_rownames = F)[-1,]
stats_estimLTX <- huxtable::set_escape_contents(stats_estimLTX,FALSE)
print(xtable::xtable(stats_estimLTX,caption = paste0("Results of estimation, with bootstrap, for Irene model with",if(data=="Ab_G1-5"){" 5 genes"}else if(data=="Ab"){" anly antibodies measured"}else if(data=="Ab_G1"){" 1 gene"}else if(data=="Ab_G1_G3"){" 2 genes"}," using ",nBoot," bootstraps and ",nrep," replicates.\n")),sanitize.text.function=function(x){x},include.rownames = F)


stats_estimLTX <- cbind(Parameters=tex2_names,stats_estim) %>% flextable() %>%
  bold(par="header") %>%
  align(part="header",align="center") %>%
  align(part="body",align="right") %>%
  add_header_lines(paste0("Results of estimation, with bootstrap, for Irene model with",if(data=="Ab_G1-5"){" 5 genes"}else if(data=="Ab"){" anly antibodies measured"}else if(data=="Ab_G1"){" 1 gene"}else if(data=="Ab_G1_G3"){" 2 genes"}," using ",nBoot," bootstraps and ",nrep," replicates.\n")) %>%
  color(i=1,part="header",color="indianred") %>%
  fontsize(i=1,part="header",size=16)  %>%
  width(width=1,unit="cm") %>%
  compose(j="Parameters",value=as_paragraph(as_equation(Parameters,width=2,height=0.5))) %>%
  hline(part="body",border = fp_border_default(color = "grey", width = 0.6)) %>%
  vline(part="body",border = fp_border_default(color = "grey", width = 0.6)) %>%
  vline(part="header",border = fp_border_default(color = "grey", width = 0.6)) %>%
  vline(j=1,part="body", border = fp_border_default( width = 1))%>%
  vline(j=1,part="header", border = fp_border_default( width = 1))


  save_as_html(stats_estimLTX, path = paste0("Results/Results",files,"_",data,"/bootstrap",type,".html"),expand=20000)
  webshot2::webshot(url=paste0("Results/Results",files,"_",data,"/bootstrap",type,".html"),file=paste0("Results/Results",files,"_",data,"/bootstrap",type,".png"),quiet=TRUE)
  unlink(paste0("Results/Results",files,"_",data,"/bootstrap",type,".html"))

