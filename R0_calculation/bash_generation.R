#!/usr/bin/env Rscript
### Generation of bash script to run transmissions from index case
library(tidyverse, quietly = T)

# Script header
header <- "#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10"

# Script footer
footer <- "exit 0"

# Paths 
rlib = "/pasteur/appa/homes/maylayan/Pacri/Hospitals_modes_transmission/github_codes/cpp/"
rinput = "/pasteur/appa/homes/maylayan/Pacri/Hospitals_modes_transmission/github_codes/out/"

# Input parameters
threshold = 60
model = "linear"
schemes = data.frame(
  scheme =c("Scheme1", "Scheme2", "Scheme3", "Scheme4", "Scheme5"),
  beta_e = c('1-45', '1-60', '1-70', '1-100', '1-150'),
  beta_c = c(0.75, 1, 1.25, 1.5, 1.75)
)

n_sim = 5000
params = crossing(schemes, expand.grid(nSim = n_sim, network = c("icu1", "icu2")))

# Generate bash files 
for (curr in 1:nrow(params)) {
  # Redirect output and error messages
  reOut <- paste0("#SBATCH -o simu_", paste0(t(params[curr, ]), collapse = "_"), ".log")
  reErr <- paste0("#SBATCH -e simu_", paste0(t(params[curr, ]), collapse = "_"), ".err")
  
  # Command to run the script
  command <- paste('srun Rscript R0-calculation.R ', rlib, " ", rinput, " ", 
    paste0(t(params[curr,]), collapse = " "))
  
  # Create script to submit
  scriptToLaunch <- paste0("simu_", paste0(t(params[curr, ]), collapse = "_"), ".sh")
  
  script <- paste0(header, "\n",
                   reOut, "\n", 
                   reErr, "\n\n",
                   command, "\n\n",
                   paste0("rm ", reOut, "\n"),
                   paste0("rm ", reErr, "\n"),
                   paste0("rm ", scriptToLaunch, "\n"),
                   footer)
  
  write(script, scriptToLaunch)
  
  # Submit to cluster queue
  system(paste("sbatch", scriptToLaunch))
}


