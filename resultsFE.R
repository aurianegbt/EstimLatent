# default
# data
data = "Ab_G1-5"
type="only_S"
files="Mult"

EstimatedStandardDeviation = data.frame()
EstimatedPopulationParameters = data.frame()

load(paste0("Results/Results",files,"_",data,"/",type,"/PopulationParameters_theo.RData"))


for(f in 1:100){
  pathToRes = paste0("Results/Results",files,"_",data,"/",type,"/estim_",f,".RData")
  if(file.exists(pathToRes)){
    load(pathToRes)
    if(ncol(res)!=0 & nrow(res)!=0){
      EstimatedStandardDeviation = rbind(EstimatedStandardDeviation, t(res[,"StandardError",drop=FALSE]))
      EstimatedPopulationParameters = rbind(EstimatedPopulationParameters,t(res[,"Parameters",drop=FALSE]))

      rownames(EstimatedPopulationParameters)[nrow(EstimatedPopulationParameters)] <- rownames(EstimatedStandardDeviation)[nrow(EstimatedStandardDeviation)] <- paste0("estim_",f)
    }
  }
}
EstimatedStandardDeviation <- t(EstimatedStandardDeviation)[names(PopulationParameters),]
EstimatedPopulationParameters <- t(EstimatedPopulationParameters)[names(PopulationParameters),]

## A VOIR

stats_estim = data.frame(target_values = setNames(as.numeric(PopulationParameters),names(PopulationParameters)))

stats_estim <- cbind(stats_estim,mean=rowMeans(EstimatedPopulationParameters),
                     bias=as.numeric(((rowMeans(EstimatedPopulationParameters)-stats_estim[,'target_values']))/(PopulationParameters)),
                     sd_emp = apply(EstimatedPopulationParameters,1,sd),
                     sd_est = rowMeans(apply(EstimatedStandardDeviation,MARGIN=2,FUN=as.numeric),na.rm = TRUE))

countNA = sapply(apply(EstimatedStandardDeviation,MARGIN=1,FUN=function(x){sum(x=="NaN")}),FUN=function(x){if(x==0){"."}else{x}})
if(!all(countNA==".")){
  stats_estim <- cbind(stats_estim,countNA)
}

coverEst = rep(0,length(PopulationParameters))
n = ncol(EstimatedPopulationParameters)
nCount=rep(n,length(PopulationParameters))
for (j in 1:n){
  sdhat = as.numeric(EstimatedStandardDeviation[,j])
  for(par in 1:length(sdhat)){
    if(sdhat[par]=="NaN"){
      nCount[par] <- nCount[par] -1
    }else{
      lw = EstimatedPopulationParameters[par,j] - 1.96*sdhat[par]
      up = EstimatedPopulationParameters[par,j] + 1.96*sdhat[par]
      bool = as.numeric(PopulationParameters[par] > lw & PopulationParameters[par] < up)
      if(!is.nan(coverEst[par])){
        coverEst[par] = coverEst[par] + bool
      }
    }
  }
}
coverEst <- sapply(coverEst/nCount,FUN=function(x){round(x,digits=3)})

coverEmp <- rep(0,length(PopulationParameters))
for (j in 1:n){
  sdhat = stats_estim[,"sd_emp"]
  for(par in 1:length(sdhat)){
    if(sdhat[par]=="NaN"){
      coverEmp[par] <- NaN
    }else{
      lw = EstimatedPopulationParameters[par,j] - 1.96*sdhat[par]
      up = EstimatedPopulationParameters[par,j] + 1.96*sdhat[par]
      bool = as.numeric(PopulationParameters[par] > lw & PopulationParameters[par] < up)
      if(!is.nan(coverEmp[par])){
        coverEmp[par] = coverEmp[par] + bool
      }
    }
  }
}
coverEmp <- coverEmp/n

stats_estim <- cbind(stats_estim,cover_emp=coverEmp,cover_est = coverEst)

if("countNA" %in% colnames(stats_estim)){
  colnames(stats_estim) <- c("Target Value","Mean of estimation","Relative Bias","Empirical Standard Deviation","Estimated Standard Deviation","(NA)","Empirical Cover","Estimated Cover")
}else{
  colnames(stats_estim) <- c("Target Value","Mean of estimation","Relative Bias","Empirical Standard Deviation","Estimated Standard Deviation","Empirical Cover","Estimated Cover")
}



stats_estimLTX=stats_estim
if(files==""){
  if(data=="Ab"){
    tex2_names = c(
      "${f_M}_1$",
      "$\\delta_S$",
      "$\\theta$",
      "$\\delta_{Ab}$",
      "$\\omega_{{f_M}_1}$",
      "$\\omega_{\\delta_S}$",
      "$\\omega_\\theta$",
      "$\\sigma_{AB}$"
    )
  }else if(data=="Ab_G1"){
    tex2_names = c(
      "${f_M}_1$",
      "$\\delta_S$",
      "$\\theta$",
      "$\\delta_{Ab}$",
      "$\\alpha_1$",
      "$\\omega_{{f_M}_1}$",
      "$\\omega_{\\delta_S}$",
      "$\\omega_\\theta$",
      "$\\sigma_{AB}$",
      "$\\sigma_{G1}$"
    )
  }else if(data=="Ab_G2"){
    tex2_names = c(
      "${f_M}_1$",
      "$\\delta_S$",
      "$\\theta$",
      "$\\delta_{Ab}$",
      "$\\alpha_2$",
      "$\\omega_{{f_M}_1}$",
      "$\\omega_{\\delta_S}$",
      "$\\omega_\\theta$",
      "$\\sigma_{AB}$",
      "$\\sigma_{G2}$"
    )
  }else{
    tex2_names = c(
      "${f_M}_1$",
      "$\\delta_S$",
      "$\\theta$",
      "$\\delta_{Ab}$",
      "$\\alpha_1$",
      "$\\alpha_2$",
      "$\\omega_{{f_M}_1}$",
      "$\\omega_{\\delta_S}$",
      "$\\omega_\\theta$",
      "$\\sigma_{AB}$",
      "$\\sigma_{G1}$",
      "$\\sigma_{G2}$"
    )
  }
  if(type=="fixed_delta"){
    tex2_names <- setdiff(tex2_names,c("$\\delta_{Ab}$","$\\delta_S$","$\\omega_{\\delta_S}$"))
  }else if(type=="default"){
    tex2_names <- setdiff(tex2_names,c("$\\omega_{\\delta_S}$"))
  }else if(type=="only_S"){

    tex2_names <- setdiff(tex2_names,c("$\\delta_{Ab}$","$\\omega_{\\delta_S}$"))
  }
}else if(files=="2"){
  tex2_names = c(
    "${f_M}_2$",
    "$\\delta_S$",
    "$\\theta_1$",
    "$\\theta_2$",
    "$\\alpha_1$",
    "$\\omega_{{f_M}_1}$",
    "$\\omega_{{f_M}_2}$",
    "$\\omega_{\\theta_1}$",
    "$\\omega_{\\theta_2}$",
    "$\\sigma_{AB}$",
    "$\\sigma_{G1}$"
  )
}else if(files=="Mult"){
  if(type=="default"){
    tex2_names = c(
      "${f_M}_1$",
      "$\\delta_S$",
      "$\\theta$",
      "$\\delta_{Ab}$",
      "$\\alpha_1$",
      "$\\alpha_2$",
      "$\\alpha_3$",
      "$\\alpha_4$",
      "$\\alpha_5$",
      "$\\omega_{{f_M}_1}$",
      "$\\omega_\\theta$",
      "$\\sigma_{AB}$",
      "$\\sigma_{G1}$",
      "$\\sigma_{G2}$",
      "$\\sigma_{G3}$",
      "$\\sigma_{G4}$",
      "$\\sigma_{G5}$"
    )
  }else if(type=="only_S"){
    tex2_names = c(
      "${f_M}_1$",
      "$\\delta_S$",
      "$\\theta$",
      "$\\alpha_1$",
      "$\\alpha_2$",
      "$\\alpha_3$",
      "$\\alpha_4$",
      "$\\alpha_5$",
      "$\\omega_{{f_M}_1}$",
      "$\\omega_\\theta$",
      "$\\sigma_{AB}$",
      "$\\sigma_{G1}$",
      "$\\sigma_{G2}$",
      "$\\sigma_{G3}$",
      "$\\sigma_{G4}$",
      "$\\sigma_{G5}$"
    )
  }

}

cat("\n \n ----------- LATEX CODE FOR RESULTS -----------\n \n")
rownames(stats_estimLTX) <- tex2_names
print(xtable::xtable(stats_estimLTX,digits=3),sanitize.text.function=function(x){x})

