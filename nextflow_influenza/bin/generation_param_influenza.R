### GENERATION DES PARAMETRES ###
args <- commandArgs(trailingOnly = TRUE)

library(dplyr)

nSim = 1:500
pathways = 1:5
networks = c("icu1", "icu2")
interventions = c("None", "Hand hygiene", "Symptomatic masking", "Universal masking", "Improved ventilation patients", "Improved ventilation hcws", "Mixed1", "Mixed2")

params = expand.grid(
  sim = nSim,
  pathway = pathways,
  network = networks,
  intervention = interventions
) %>%
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
    .before = network
  ) %>%
  select(-pathway)
  

 write.table(params, file="param_grid_influenza.txt", row.names = F)
