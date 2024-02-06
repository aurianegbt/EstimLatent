estimateFS <- function(fileNumber,data,bootstrap=FALSE,type,files=""){
  temporaryDirectory = paste0("/beegfs/agabaut/tmp",fileNumber,"_",data)
  if(dir.exists(temporaryDirectory)){unlink(temporaryDirectory,recursive = T)}
  dir.create(temporaryDirectory)

  pathToSim = paste0("Model/IreneModel",files,"/Simulation/simulatedData_",data,".csv")

  sim = read.csv(pathToSim) %>%
    filter(rep==fileNumber) %>%
    select(ID,TIME,obs,obsid)

  if(file.exists(paste0("Results/Results",files,"_",data,"/",if(bootstrap){"bootstrap"},type,"/estim_",fileNumber,".RData"))){
    load(paste0("Results/Results",files,"_",data,"/",if(bootstrap){"bootstrap"},type,"/estim_",fileNumber,".RData"))
    doRes = !exists("res")
    if(exists("resBoot")){
      bootMISS = setdiff(1:200,as.numeric(stringr::str_remove_all(names(resBoot),"boot_")))
    }
  }else{
    doRes=TRUE
    bootMISS=1:200
  }

  ## Normal part
  if(doRes){

    #project
    write.csv(sim,file=paste0(temporaryDirectory,"/sim.txt"),quote = F,row.names = F)

    newProject(data = list(dataFile = paste0(temporaryDirectory,"/sim.txt"),
                           headerTypes = c("id","time","observation","obsid")),
               modelFile = paste0("Files/model_generated",files,"_",data,".txt"))

    # individual variability
    setIndividualParameterVariability(delta_V=FALSE,delta_Ab=FALSE) # never
    #mu
    if(files=="Mu"){
      setIndividualParameterVariability(mu_1=F,mu_2=F,mu_3=F,mu_4=F,mu_5=F) # never
    }
    # alpha
    if(data!="Ab"){
      if(files=="" | files=="2" | data %in% c("Ab_G1","Ab_G1_G3")){
        if(data=="Ab_G1"){
          setIndividualParameterVariability(alpha_1=FALSE)
        }else if(data=="Ab_G1_G3"){
          setIndividualParameterVariability(alpha_1=FALSE,alpha_3=FALSE)
        }else if(data=="Ab_G2"){
          setIndividualParameterVariability(alpha_2=FALSE)
        }
      }else if(files %in% c("Mult","Mu","Mu_i")){
        setIndividualParameterVariability(alpha_1=FALSE,alpha_2=FALSE,alpha_3=F,alpha_4=F,alpha_5=F)
        if(files=="Mu"){
          setIndividualParameterVariability(mu_1=FALSE,mu_2=FALSE,mu_3=F,mu_4=F,mu_5=F)
        }
      }
    }
    # delta
    if(type=="default"){
      setIndividualParameterVariability(delta_S=FALSE)
      setPopulationParameterInformation(
        delta_V_pop=list(initialValue=2.7,method="FIXED"))
    }else if(type=="default_S"){
      setPopulationParameterInformation(
        delta_V_pop=list(initialValue=2.7,method="FIXED"))
    }else if(type=="fixed_delta"){
      setIndividualParameterVariability(delta_S=FALSE)
      if(files %in% c("","Mult","Mu","Mu_i")){
        setPopulationParameterInformation(
          delta_V_pop=list(initialValue=2.7,method="FIXED"),
          delta_S_pop=list(initialValue=0.01,method="FIXED"),
          delta_Ab_pop=list(initialValue=0.03,method="FIXED"))
      }else if(files=="2"){
        setPopulationParameterInformation(
          delta_V_pop=list(initialValue=2.7,method="FIXED"),
          delta_S_pop=list(initialValue=0.0322,method="FIXED"),
          delta_Ab_pop=list(initialValue=0.08,method="FIXED"))
      }
    }else if(type=="only_S"){
      setIndividualParameterVariability(delta_S=FALSE)
      if(files %in% c("","Mult","Mu","Mu_i")){
        setPopulationParameterInformation(
          delta_V_pop=list(initialValue=2.7,method="FIXED"),
          delta_Ab_pop=list(initialValue=0.03,method="FIXED"))
      }else if(files=="2"){
        setPopulationParameterInformation(
          delta_V_pop=list(initialValue=2.7,method="FIXED"),
          delta_Ab_pop=list(initialValue=0.08,method="FIXED"))
      }
    }else if(type=="poor"){
      setIndividualParameterVariability(delta_S=FALSE)
      if(files %in% c("","Mult","Mu","Mu_i")){
        setPopulationParameterInformation(
          fM1_pop=list(initialValue=4.5,method="FIXED"),
          delta_V_pop=list(initialValue=2.7,method="FIXED"),
          delta_Ab_pop=list(initialValue=0.03,method="FIXED"))
      }else if(files=="2"){
        setPopulationParameterInformation(
          fM1_pop=list(initialValue=4.5,method="FIXED"),
          delta_V_pop=list(initialValue=2.7,method="FIXED"),
          delta_Ab_pop=list(initialValue=0.08,method="FIXED"))
      }
    }

    # Error Model
    if(data=="Ab"){
      setErrorModel(yyAB="constant")
    }else if(data=="Ab_G1"){
      setErrorModel(yyAB="constant",yyG1="constant")
    }else if(data=="Ab_G1_G3"){
      setErrorModel(yyAB="constant",yyG1="constant",yyG3="constant")
    }else if(data=="Ab_G2"){
      setErrorModel(yyAB="constant",yyG2="constant")
    }else if(data=="Ab_G1-5"){
      setErrorModel(yyAB="constant",yyG1="constant",yyG2="constant",yyG3="constant",yyG4="constant",yyG5="constant")
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
    if(type=="default" | type=="default_S"){
      keep = setdiff(names(getEstimatedPopulationParameters()),c("delta_V_pop"))
    }else if(type=="fixed_delta"){
      keep = setdiff(names(getEstimatedPopulationParameters()),c("delta_V_pop","delta_Ab_pop","delta_S_pop"))
    }else if(type=="only_S"){
      keep = setdiff(names(getEstimatedPopulationParameters()),c("delta_V_pop","delta_Ab_pop"))
    }else if(type=="poor"){
      keep = setdiff(names(getEstimatedPopulationParameters()),c("delta_V_pop","delta_Ab_pop","fM1_pop"))
    }
    if(files=="2"){
      keep <- setdiff(keep,"fM1_pop")
    }

    EstimatedPopulationParameters = getEstimatedPopulationParameters()[keep]
    EstimatedStandardErrors = getEstimatedStandardErrors()[[1]]


    res = data.frame(Parameters = EstimatedPopulationParameters,StandardError = EstimatedStandardErrors$se)
    if(data=="Ab"){
      rownames(res)[which(rownames(res)=="a")] <- "sigma_Ab"
    }else{
      rownames(res) <-  stringr::str_replace_all(rownames(res),"ay","sigma_")
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
      IDi = sample(unique(sim$ID),replace=TRUE)
      for(k in 1:length(IDi)){
        simi <- rbind(simi,sim[sim$ID==IDi[k],] %>% mutate(ID=k) %>% mutate(original_ID = IDi[k],.after="ID"))
      }
      write.csv(simi,file=paste0(temporaryDirectory,"/sim_",i,".txt"),quote = FALSE,row.names = FALSE)

      newProject(data = list(dataFile = paste0(temporaryDirectory,"/sim_",i,".txt"),
                             headerTypes = c("id","ignore","time","observation","obsid")),
                 modelFile = paste0("Files/model_generated",files,"_",data,".txt"))

      # individual variability
      setIndividualParameterVariability(delta_V=FALSE,delta_Ab=FALSE) # never
      # alpha
      if(data!="Ab"){
        if(files=="" | files=="2" | data %in% c("Ab_G1","Ab_G1_G3")){
          if(data=="Ab_G1"){
            setIndividualParameterVariability(alpha_1=FALSE)
          }else if(data=="Ab_G1_G3"){
            setIndividualParameterVariability(alpha_1=FALSE,alpha_3=FALSE)
          }else if(data=="Ab_G2"){
            setIndividualParameterVariability(alpha_2=FALSE)
          }
        }else if(files %in% c("Mult","Mu","Mu_i")){
          setIndividualParameterVariability(alpha_1=FALSE,alpha_2=FALSE,alpha_3=F,alpha_4=F,alpha_5=F)
          if(files=="Mu"){
            setIndividualParameterVariability(mu_1=FALSE,mu_2=FALSE,mu_3=F,mu_4=F,mu_5=F)
          }
        }
      }
      # fM1
      if(files=="2"){
        setPopulationParameterInformation()
      }
      # delta
      if(type=="default"){
        setIndividualParameterVariability(delta_S=FALSE)
        setPopulationParameterInformation(
          delta_V_pop=list(initialValue=2.7,method="FIXED"))
      }else if(type=="default_S"){
        setPopulationParameterInformation(
          delta_V_pop=list(initialValue=2.7,method="FIXED"))
      }else if(type=="fixed_delta"){
        setIndividualParameterVariability(delta_S=FALSE)
        if(files %in% c("","Mult","Mu","Mu_i")){
          setPopulationParameterInformation(
            delta_V_pop=list(initialValue=2.7,method="FIXED"),
            delta_S_pop=list(initialValue=0.01,method="FIXED"),
            delta_Ab_pop=list(initialValue=0.03,method="FIXED"))
        }else if(files=="2"){
          setPopulationParameterInformation(
            delta_V_pop=list(initialValue=2.7,method="FIXED"),
            delta_S_pop=list(initialValue=0.0322,method="FIXED"),
            delta_Ab_pop=list(initialValue=0.08,method="FIXED"))
        }
      }else if(type=="only_S"){
        setIndividualParameterVariability(delta_S=FALSE)
        if(files %in% c("","Mult","Mu","Mu_i")){
          setPopulationParameterInformation(
            delta_V_pop=list(initialValue=2.7,method="FIXED"),
            delta_Ab_pop=list(initialValue=0.03,method="FIXED"))
        }else if(files=="2"){
          setPopulationParameterInformation(
            delta_V_pop=list(initialValue=2.7,method="FIXED"),
            delta_Ab_pop=list(initialValue=0.08,method="FIXED"))
        }
      }else if(type=="poor"){
        setIndividualParameterVariability(delta_S=FALSE)
        if(files %in% c("","Mult","Mu","Mu_i")){
          setPopulationParameterInformation(
            fM1_pop=list(initialValue=4.5,method="FIXED"),
            delta_V_pop=list(initialValue=2.7,method="FIXED"),
            delta_Ab_pop=list(initialValue=0.03,method="FIXED"))
        }else if(files=="2"){
          setPopulationParameterInformation(
            fM1_pop=list(initialValue=4.5,method="FIXED"),
            delta_V_pop=list(initialValue=2.7,method="FIXED"),
            delta_Ab_pop=list(initialValue=0.08,method="FIXED"))
        }
      }

      # Error Model
      if(data=="Ab"){
        setErrorModel(yyAB="constant")
      }else if(data=="Ab_G1"){
        setErrorModel(yyAB="constant",yyG1="constant")
      }else if(data=="Ab_G1_G3"){
        setErrorModel(yyAB="constant",yyG1="constant",yyG3="constant")
      }else if(data=="Ab_G2"){
        setErrorModel(yyAB="constant",yyG2="constant")
      }else if(data=="Ab_G1-5"){
        setErrorModel(yyAB="constant",yyG1="constant",yyG2="constant",yyG3="constant",yyG4="constant",yyG5="constant")
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
      if(type=="default" | type=="default_S"){
        keep = setdiff(names(getEstimatedPopulationParameters()),c("delta_V_pop"))
      }else if(type=="fixed_delta"){
        keep = setdiff(names(getEstimatedPopulationParameters()),c("delta_V_pop","delta_Ab_pop","delta_S_pop"))
      }else if(type=="only_S"){
        keep = setdiff(names(getEstimatedPopulationParameters()),c("delta_V_pop","delta_Ab_pop"))
      }else if(type=="poor"){
        keep = setdiff(names(getEstimatedPopulationParameters()),c("delta_V_pop","delta_Ab_pop","fM1_pop"))
      }
      if(files=="2"){
        keep <- setdiff(keep,"fM1_pop")
      }

      EstimatedPopulationParameters = getEstimatedPopulationParameters()[keep]
      EstimatedStandardErrors = getEstimatedStandardErrors()[[1]]

      resaux = data.frame(Parameters = EstimatedPopulationParameters,StandardError = EstimatedStandardErrors$se)
      if(data=="Ab"){
        rownames(resaux)[which(rownames(resaux)=="a")] <- "sigma_Ab"
      }else{
        rownames(resaux) <-  stringr::str_replace_all(rownames(resaux),"ay","sigma_")
      }

      resaux
    })
    stopCluster(cl)

    resBoot <- setNames(resBoot,paste0("boot_",1:length(resBoot)))
  }

  unlink(temporaryDirectory,recursive=TRUE)
  if(bootstrap){
    save(res,resBoot,file=paste0("Results/Results",files,"_",data,"/",if(bootstrap){"bootstrap"},type,"/estim_",fileNumber,".RData"))
  }else{
    save(res,file=paste0("Results/Results",files,"_",data,"/",if(bootstrap){"bootstrap"},type,"/estim_",fileNumber,".RData"))
  }
}
