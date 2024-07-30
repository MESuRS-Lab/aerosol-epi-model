#!/bin/bash

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=olivier.gaufres@pasteur.fr
#SBATCH --output=/pasteur/appa/homes/ogaufres/logs/intervention-analysis/simulation-%j.out
#SBATCH --error=/pasteur/appa/homes/ogaufres/logs/intervention-analysis/simulation-%j.err

module purge
module load R/4.4.0

# vars from launch-intervention-analysis.sh
beta=${beta:-"ERROR-NO-BETA"}
nu=${nu:-"ERROR-NO-NU"}
sim_id=${sim_id:-"ERROR-NO-ID"}

# pass vars to R script
Rscript $HOME/intervention-analysis/run_simulation_intervention.R $sim_id $beta $nu