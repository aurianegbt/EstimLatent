estimateFS <- function(fileNumber,data,bootstrap=FALSE,type,files="",Nbr_ind=50){
  temporaryDirectory = paste0("/beegfs/agabaut/tmp",fileNumber,"_",paste0(data,collapse="_"))
  if(dir.exists(temporaryDirectory)){unlink(temporaryDirectory,recursive = T)}
  dir.create(temporaryDirectory)


  sim.list = lapply(data,FUN=function(d){read.csv(paste0("Model/IreneModel",files,"/Simulation/simulation_y",d,".txt"))})

  sim <- data.frame()
  for(d in 1:length(data)){
    sim <- sim %>% rbind(sim.list[[d]] %>% mutate(obsid = paste0("y",data[d]),.before = group) %>% rename(obs = paste0("y",data[d]))
      )
  }

  sim = sim %>%
    filter(rep==fileNumber) %>%
    filter(if(Nbr_ind==50){group=="simulationGroup1"}else if(Nbr_ind==15){group=="simulationGroup2"}) %>%
    mutate(id=if(Nbr_ind==15){id-50}else{id}) %>%
    select(-group,-rep)

  if(file.exists(paste0("Results/Results",files,if(Nbr_ind!=50){Nbr_ind},"_",paste0(data,collapse="_"),"/",if(bootstrap){"bootstrap"},type,"/estim_",fileNumber,".RData"))){
    load(paste0("Results/Results",files,if(Nbr_ind!=50){Nbr_ind},"_",paste0(data,collapse="_"),"/",if(bootstrap){"bootstrap"},type,"/estim_",fileNumber,".RData"))
    doRes = !exists("res")
    if(exists("resBoot")){
      bootMISS = setdiff(1:500,as.numeric(stringr::str_remove_all(names(resBoot),"boot_")))
    }
  }else{
    doRes=TRUE
    bootMISS=1:500
  }

  ## Normal part
  if(doRes){

    #project
    write.csv(sim,file=paste0(temporaryDirectory,"/sim.txt"),quote = F,row.names = F)

    newProject(data = list(dataFile = paste0(temporaryDirectory,"/sim.txt"),
                           headerTypes = c("id","time","observation","obsid")),
               modelFile = paste0("Files/model_generated_",paste0(data,collapse="_"),".txt"))

    # individual variability
    setIndividualParameterVariability(delta_V=FALSE,delta_Ab=FALSE) # never
    # alpha and mu
    for(d in data[stringr::str_detect(data,"G")]){
      d <- stringr::str_remove(d,"G")
      eval(parse(text=paste0("setIndividualParameterVariability(alpha_",d,"=FALSE)")))
      # if(files == "Mu"){
      #   eval(parse(text="setIndividualParameterVariability(mu",d,"=F)"))
      # }
    }
    # delta
    setPopulationParameterInformation(delta_V_pop=list(initialValue=2.7,method="FIXED"))
    setIndividualParameterVariability(delta_S=FALSE)
    if(type=="only_S"){
      setPopulationParameterInformation(delta_Ab_pop=list(initialValue=0.03,method="FIXED"))
    }
    # else if(type=="poor"){
    #   if(files %in% c("Mult","Mu","Mu_i")){
    #     setPopulationParameterInformation(
    #       fM1_pop=list(initialValue=4.5,method="FIXED"))
    #   }else if(files=="MP"){
    #     setPopulationParameterInformation(
    #       fM1_pop=list(initialValue=4.132851589,method="FIXED"))
    #   }
    # }

    # Error Model
    if("AB" %in% data){
      setErrorModel(yyAB="constant")
    }
    for(d in data[stringr::str_detect(data,"G")]){
      eval(parse(text=paste0("setErrorModel(yy",d,"='constant')")))
    }

    # Initialization (plafrim form)
    parameters = getPopulationParameterInformation()
    optimizedParameters = getFixedEffectsByAutoInit(parameters)
    setPopulationParameterInformation(optimizedParameters)

    # Estimated standard deviation
    scenario <- getScenario()
    scenario$tasks['standardErrorEstimation'] <- TRUE
    scenario$tasks['conditionalDistributionSampling'] <- FALSE
    scenario$tasks['conditionalModeEstimation'] <- FALSE
    scenario$tasks['plots'] <- FALSE
    setScenario(scenario)

    saveProject(projectFile = paste0(temporaryDirectory,"/MlxProj.mlxtran"))

    runScenario()

    # Results format
    if(type=="default"){
      keep = setdiff(names(getEstimatedPopulationParameters()),c("delta_V_pop"))
    }else if(type=="only_S"){
      keep = setdiff(names(getEstimatedPopulationParameters()),c("delta_V_pop","delta_Ab_pop"))
    }
    # else if(type=="poor"){
    #   keep = setdiff(names(getEstimatedPopulationParameters()),c("delta_V_pop","delta_Ab_pop","fM1_pop"))
    # }

    EstimatedPopulationParameters = getEstimatedPopulationParameters()[keep]
    EstimatedStandardErrors = getEstimatedStandardErrors()[[1]]


    res = data.frame(Parameters = EstimatedPopulationParameters,StandardError = EstimatedStandardErrors$se)
    if(identical("AB",data)){
      rownames(res)[which(rownames(res)=="a")] <- "sigma_Ab"
    }else{
      rownames(res) <-  stringr::str_replace_all(rownames(res),"ay","sigma_")
      if("sigma_AB" %in% rownames(res)){
        rownames(res)[which(rownames(res)=="sigma_AB")] <- "sigma_Ab"
      }
    }
  }

  doParallel::registerDoParallel(cl <- parallel::makeCluster(parallel::detectCores()))

  if(bootstrap){
    if(!exists("resBoot")){
      resBoot=list()
    }

    resBoot = append(resBoot,foreach::foreach(i = bootMISS) %dopar% {
      suppressMessages(library(lixoftConnectors))
      suppressMessages(library(dplyr))
      invisible(purrr::quietly(initializeLixoftConnectors)(software='monolix', force=T))

      # project
      simi <- data.frame()
      IDi = sample(unique(sim$id),replace=TRUE)
      for(k in 1:length(IDi)){
        simi <- rbind(simi,sim[sim$id==IDi[k],] %>% mutate(id=k) %>% mutate(original_id = IDi[k],.after="id"))
      }
      write.csv(simi,file=paste0(temporaryDirectory,"/sim_",i,".txt"),quote = FALSE,row.names = FALSE)

      newProject(data = list(dataFile = paste0(temporaryDirectory,"/sim_",i,".txt"),
                             headerTypes = c("id","ignore","time","observation","obsid")),
                 modelFile = paste0("Files/model_generated_",paste0(data,collapse="_"),".txt"))

      # individual variability
      setIndividualParameterVariability(delta_V=FALSE,delta_Ab=FALSE) # never
      # alpha and mu
      for(d in data[stringr::str_detect(data,"G")]){
        d <- stringr::str_remove(d,"G")
        eval(parse(text=paste0("setIndividualParameterVariability(alpha_",d,"=FALSE)")))
        # if(files == "Mu"){
        #   eval(parse(text="setIndividualParameterVariability(mu",d,"=F)"))
        # }
      }
      # delta
      setPopulationParameterInformation(delta_V_pop=list(initialValue=2.7,method="FIXED"))
      setIndividualParameterVariability(delta_S=FALSE)
      if(type=="only_S"){
        setPopulationParameterInformation(delta_Ab_pop=list(initialValue=0.03,method="FIXED"))
      }
      # else if(type=="poor"){
      #   if(files %in% c("Mult","Mu","Mu_i")){
      #     setPopulationParameterInformation(
      #       fM1_pop=list(initialValue=4.5,method="FIXED"))
      #   }else if(files=="MP"){
      #     setPopulationParameterInformation(
      #       fM1_pop=list(initialValue=4.132851589,method="FIXED"))
      #   }
      # }

      # Error Model
      if("AB" %in% data){
        setErrorModel(yyAB="constant")
      }
      for(d in data[stringr::str_detect(data,"G")]){
        eval(parse(text=paste0("setErrorModel(yy",d,"='constant')")))
      }

      # Initialization (plafrim form)
      parameters = getPopulationParameterInformation()
      optimizedParameters = getFixedEffectsByAutoInit(parameters)
      setPopulationParameterInformation(optimizedParameters)

      # Estimated standard deviation
      scenario <- getScenario()
      scenario$tasks['standardErrorEstimation'] <- TRUE
      scenario$tasks['conditionalDistributionSampling'] <- FALSE
      scenario$tasks['conditionalModeEstimation'] <- FALSE
      scenario$tasks['plots'] <- FALSE
      setScenario(scenario)

      saveProject(projectFile = paste0(temporaryDirectory,"/MlxProj",i,".mlxtran"))

      runScenario()

      # Results format
      if(type=="default"){
        keep = setdiff(names(getEstimatedPopulationParameters()),c("delta_V_pop"))
      }else if(type=="only_S"){
        keep = setdiff(names(getEstimatedPopulationParameters()),c("delta_V_pop","delta_Ab_pop"))
      }
      # else if(type=="poor"){
      #   keep = setdiff(names(getEstimatedPopulationParameters()),c("delta_V_pop","delta_Ab_pop","fM1_pop"))
      # }

      EstimatedPopulationParameters = getEstimatedPopulationParameters()[keep]
      EstimatedStandardErrors = getEstimatedStandardErrors()[[1]]

      resaux = data.frame(Parameters = EstimatedPopulationParameters,StandardError = EstimatedStandardErrors$se)
      if(identical("AB",data)){
        rownames(resaux)[which(rownames(resaux)=="a")] <- "sigma_Ab"
      }else{
        rownames(resaux) <-  stringr::str_replace_all(rownames(resaux),"ay","sigma_")
        if("sigma_AB" %in% rownames(resaux)){
          rownames(resaux)[which(rownames(resaux)=="sigma_AB")] <- "sigma_Ab"
        }
      }
      resaux
    })
    stopCluster(cl)

    resBoot <- setNames(resBoot,paste0("boot_",1:length(resBoot)))
  }

  unlink(temporaryDirectory,recursive=TRUE)
  if(bootstrap){
    save(res,resBoot,file=paste0("Results/Results",files,if(Nbr_ind!=50){Nbr_ind},"_",paste0(data,collapse=""),"/",if(bootstrap){"bootstrap"},type,"/estim_",fileNumber,".RData"))
  }else{
    save(res,file=paste0("Results/Results",files,if(Nbr_ind!=50){Nbr_ind},"_",paste0(data,collapse="_"),"/",if(bootstrap){"bootstrap"},type,"/estim_",fileNumber,".RData"))
  }
}
