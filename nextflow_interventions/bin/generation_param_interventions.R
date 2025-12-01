### GENERATION DES PARAMETRES ###
args <- commandArgs(trailingOnly = TRUE)

library(dplyr)
library(lubridate)

# Arguments 
pathways = 1:5
networks = c("icu1", "icu2")
threshold = 60
models = "linear" 
interventions = c("None", "Hand hygiene", "Symptomatic masking", "Universal masking", "Improved ventilation patients", "Improved ventilation hcws", "Mixed1", "Mixed2")
nSim = 1:500

# Generate random starting seeds and random index cases 
counterfactuals = data.frame()

for (n in networks) {
  for (i in nSim) {
    for (p in pathways) {
      s = sample(10000000:99999999, size=1) 
      counterfactual = data.frame(
        seed = s,
        sim = i,
        network = n,
        pathway = p
      )
      counterfactuals = bind_rows(counterfactuals, counterfactual)
    }
  }

}

# Dataframe of parameters
params = expand.grid(
  sim = nSim,
  pathway = pathways,
  threshold = threshold,
  network = networks,
  model = models,
  intervention = interventions
  ) %>%
  left_join(., counterfactuals, by = c("network", "sim", "pathway")) %>%
  mutate(
    beta_c = case_when(
      pathway == 1 ~ "0.75",
      pathway == 2 ~ "1",
      pathway == 3 ~ "1.25",
      pathway == 4 ~ "1.5",
      pathway == 5 ~ "1.75"
    ),
    beta_e = case_when(
      pathway == 1 ~ "1/45",
      pathway == 2 ~ "1/60",
      pathway == 3 ~ "1/70",
      pathway == 4 ~ "1/100",
      pathway == 5 ~ "1/150"
    ),
    .before = threshold
  ) %>%
  select(-pathway)
  

 write.table(params, file="param_grid_interventions.txt", row.names = F)
