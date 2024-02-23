# empty the memory
rm(list=ls())

# Required Packages :
suppressWarnings(suppressMessages(library(lixoftConnectors)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(doParallel)))
suppressWarnings(suppressMessages(library(parallel)))
suppressWarnings(suppressMessages(library(foreach)))
suppressWarnings(suppressMessages(library(slurmR)))

# Load Functions and softwares
invisible(purrr::quietly(initializeLixoftConnectors)(software='monolix', force=T))
source("estimateFS.R")


# Arguments of batch
fileNumber <- as.numeric(Slurm_env(x='SLURM_ARRAY_TASK_ID'))

args = commandArgs(trailingOnly=TRUE)
args=c("c('AB','G1')","FALSE","only_S","Mult","15")
eval(parse(text=paste0("data <-",args[1])))
bootstrap <- as.logical(args[2])
type = args[3]
files = args[4]
Nbr_ind = as.numeric(args[5])

## Presentation of work :
cat("===========================================================================\n
    - - - - - - - - - - - Work Information - - - - - - - - - - -       ")
cat(paste0("\nWorking on file n°",fileNumber,", with ",if(bootstrap){" bootstrap, "},if(type=="only_S"){"only delta_S being estimated among delta parameters"}else if(type=="default"){" with only delta_V being fixed among delta parameters"},", from ",files," simulation",if(Nbr_ind!=50){paste0(" with ",Nbr_ind," individuals")},"."))
cat("\n===========================================================================\n")
dir <- function(x){
  if(!dir.exists(x)){dir.create(x)}
}
dir(paste0("Results/Results",files,if(Nbr_ind!=50){Nbr_ind},"_",paste0(data,collapse="_")))
dir(paste0("Results/Results",files,if(Nbr_ind!=50){Nbr_ind},"_",paste0(data,collapse="_"),"/",if(bootstrap){"bootstrap"},type))

estimateFS(fileNumber,data,bootstrap,type,files,Nbr_ind)

cat("===========================================================================\n")


