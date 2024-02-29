library(dplyr)
library(flextable)
## DATA INFORMATION
data = c("AB")
type="default"
files="Mult"
Nind = "15"

load(paste0("Results/Results",files,if(Nind!=50){Nind},"_",paste0(data,collapse="_"),"/bootstrap",type,"/estim_1.RData"))
# nBoot = length(resBoot)
nBoot = 200

# nrep = length(list.files(paste0("Results/Results",files,"_",paste0(data,collapse="_"),"/bootstrap",type)))
nrep = 100
precision = abs(floor(log10(1/nrep)))

# Gather results

EstimatedStandardDeviation = data.frame() # tableau des standard errors estimées sur le réplicat par monolix
EstimatedPopulationParameters = data.frame() # tableau des paramètres estimés sur le réplicat par monolix

EstimatedStandardDeviationBOOT = data.frame() # tableau des écart-type estimées par bootstrap (sd des paramètres estimés sur chaque boot du réplicat)
EstimatedPopulationParametersBOOT = data.frame() # tableau du paramètre estimées par bootstrap (mean des paramètres estimées sur chaque boot du réplicat)

load(paste0("Results/PopulationParameters_theo_",files,".RData"))
if(!("G1" %in% data)){
  PopulationParameters <- PopulationParameters %>% select(-alpha_1_pop,-sigma_G1)
}
if(!("G2" %in% data)){
  PopulationParameters <- PopulationParameters %>% select(-alpha_2_pop,-sigma_G2)
}
if(type=="only_S"){
  PopulationParameters <- PopulationParameters %>% select(-delta_Ab_pop)
}


for(f in 1:nrep){
  if(exists("resBoot"))
    rm(resBoot)
  if(exists("res"))
    rm(res)
  suppressWarnings(tryCatch({
    pathToRes = paste0("Results/Results",files,if(Nind!=50){Nind},"_",paste0(data,collapse="_"),"/bootstrap",type,"/estim_",f,".RData")
    load(pathToRes)
    if(length(resBoot)<nBoot)
      cat("file n°",f," : ",length(resBoot),"\n")

      EstimatedStandardDeviation = rbind(EstimatedStandardDeviation, t(res[,"StandardError",drop=FALSE]))
      EstimatedPopulationParameters = rbind(EstimatedPopulationParameters,t(res[,"Parameters",drop=FALSE]))

    Est = Reduce(cbind,lapply(resBoot[1:nBoot],function(dt){dt[,"Parameters",drop=FALSE]}))

    EstimatedStandardDeviationBOOT = rbind(EstimatedStandardDeviationBOOT, t(as.data.frame(apply(Est,1,sd))))
    EstimatedPopulationParametersBOOT = rbind(EstimatedPopulationParametersBOOT,t(as.data.frame(rowMeans(Est))))

    rownames(EstimatedPopulationParametersBOOT)[nrow(EstimatedPopulationParametersBOOT)] <-
      rownames(EstimatedStandardDeviationBOOT)[nrow(EstimatedStandardDeviationBOOT)] <-
      rownames(EstimatedPopulationParameters)[nrow(EstimatedPopulationParameters)] <-
      rownames(EstimatedStandardDeviation)[nrow(EstimatedStandardDeviation)] <-
      paste0("estim_",f)

  },error=function(e){cat(f,",")}))
}
if("sigma_Ab" %in% colnames(EstimatedStandardDeviation)){
  EstimatedStandardDeviation <- EstimatedStandardDeviation %>%
    rename(sigma_AB=sigma_Ab)
  EstimatedStandardDeviationBOOT <- EstimatedStandardDeviationBOOT %>%
    rename(sigma_AB=sigma_Ab)
  EstimatedPopulationParameters <- EstimatedPopulationParameters%>%
    rename(sigma_AB=sigma_Ab)
  EstimatedPopulationParametersBOOT <- EstimatedPopulationParametersBOOT %>%
    rename(sigma_AB=sigma_Ab)
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
                     bias=as.numeric(((rowMeans(EstimatedPopulationParameters)-stats_estim[,'target_values']))/as.numeric((PopulationParameters))),
                     biasBoot = as.numeric(((rowMeans(EstimatedPopulationParametersBOOT)-stats_estim[,'target_values']))/as.numeric((PopulationParameters))),
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


tex2_names = unname(c(fM1_pop = "{f_M}_1",
                      delta_S_pop = "\\delta_S",
                      theta_pop = "\\theta",
                      delta_Ab_pop = "\\delta_{Ab}",
                      alpha_1_pop = "\\alpha_1",
                      alpha_2_pop = "\\alpha_2",
                      omega_fM1 = "\\omega_{{f_M}_1}",
                      omega_delta_S = "\\omega_{\\delta_S}",
                      omega_theta = "\\omega_\\theta",
                      sigma_AB = "\\sigma_{Ab}",
                      sigma_G1 = "\\sigma_{G1}",
                      sigma_G2 = "\\sigma_{G2}")[names(PopulationParameters)])
options(scipen=999)

stats_estim[,!(colnames(stats_estim) %in% c("Target Value","(NA)"))] <- apply(stats_estim[,!(colnames(stats_estim) %in% c("Target Value","(NA)"))],2, FUN=function(x){as.character(round(x,digits=precision+1))})


# cat("\n \n ----------- LATEX CODE FOR RESULTS -----------\n \n")
dir <- function(d){if(!dir.exists(d)){dir.create(d)}}
dir("Results/Table_Results")
dir(paste0("Results/Table_Results/",files,"/"))
dir(paste0("Results/Table_Results/",files,"/",paste0(data,collapse = "_")))
pathToResults = paste0("Results/Table_Results/",files,"/",paste0(data,collapse = "_"),"/bootstrap",type,"_",Nind,"ind")

stats_estimLTX <- huxtable::huxtable(cbind(Parameter = paste0("$",tex2_names,"$"),stats_estim),add_rownames = F)[-1,]
stats_estimLTX <- huxtable::set_escape_contents(stats_estimLTX,FALSE)
texTable = xtable::xtable(stats_estimLTX,caption = paste0("Results of estimation, with bootstrap, for Irene model with antibody Ab, ",if("G1" %in% data){" and a high noised gene,"}else if("G2" %in% data){" and a low noised gene,"}," using ",nBoot," bootstraps and ",nrep," replicates, with ",Nind," individuals.\n"))
# print(texTable,sanitize.text.function=function(x){x},include.rownames = F)
latex_code <- capture.output(print(texTable,sanitize.text.function=function(x){x},include.rownames = F))[-c(1,2)]
writeLines(latex_code, paste0(pathToResults,".txt"))




stats_estimLTX <- cbind(Parameters=tex2_names,stats_estim) %>% flextable() %>%
  bold(par="header") %>%
  align(part="header",align="center") %>%
  align(part="body",align="right") %>%
  add_header_lines(paste0("Results of estimation, with bootstrap, for Irene model with antibody Ab, ",if("G1" %in% data){" and a high noised gene,"}else if("G2" %in% data){" and a low noised gene,"}," using ",nBoot," bootstraps and ",nrep," replicates.\n")) %>%
  color(i=1,part="header",color="indianred") %>%
  fontsize(i=1,part="header",size=16)  %>%
  width(width=1,unit="cm") %>%
  compose(j="Parameters",value=as_paragraph(as_equation(Parameters,width=2,height=0.5))) %>%
  hline(part="body",border = fp_border_default(color = "grey", width = 0.6)) %>%
  vline(part="body",border = fp_border_default(color = "grey", width = 0.6)) %>%
  vline(part="header",border = fp_border_default(color = "grey", width = 0.6)) %>%
  vline(j=1,part="body", border = fp_border_default( width = 1))%>%
  vline(j=1,part="header", border = fp_border_default( width = 1))


  save_as_html(stats_estimLTX, path = paste0(pathToResults,".html"),expand=20000)
  webshot2::webshot(url=paste0(pathToResults,".html"),file=paste0(pathToResults,".png"),quiet=TRUE)
  unlink(paste0(pathToResults,".html"))

