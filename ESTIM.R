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
args=c("Ab_G1","FALSE","only_S","")
data <- args[1]
bootstrap <- as.logical(args[2])
type = args[3]
files = args[4]

## Presentation of work :
cat("===========================================================================\n
    - - - - - - - - - - - Work Information - - - - - - - - - - -       ")
cat(paste0("\nWorking on file n°",fileNumber,", with",if(bootstrap){" bootstrap, "},if(type=="only_S"){"only delta_S being estimated among delta parameters"}else if(type=="default"){" with onlye delta_V being fixed among delta parameters"},", from ",files," simulation."))
cat("\n===========================================================================\n")
dir <- function(x){
  if(!dir.exists(x)){dir.create(x)}
}
dir(paste0("Results/Results",files,"_",data))
dir(paste0("Results/Results",files,"_",data,"/",if(bootstrap){"bootstrap"},type))

estimateFS(fileNumber,data,bootstrap,type,files)

cat("===========================================================================\n")


