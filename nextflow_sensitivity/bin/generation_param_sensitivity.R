### GENERATION DES PARAMETRES ###
args <- commandArgs(trailingOnly = TRUE)

library(dplyr)

nSim = 1:500
pathways = 1:5
networks = c("icu1", "icu2")

params1 = expand.grid(
  sim = nSim,
  pathway = pathways,
  network = networks,
  sensitivity =  c("low-nu", "high-nu", "low-mu", "high-mu", "low-hh", "high-hh")
  ) %>%
  mutate(
      intervention = case_when(
        grepl("-nu", sensitivity) ~ "None",
        grepl("-mu", sensitivity) ~ "None",
        grepl("-hh", sensitivity) ~ "Hand hygiene"
      )
    )

params_masks = expand.grid(
  sim = nSim,
  pathway = pathways,
  network = networks,
  sensitivity = c("low-mw", "high-mw"),
  intervention = c("Symptomatic masking", "Universal masking", "Mixed1", "Mixed2")
  )

params_ventilation = expand.grid(
  sim = nSim,
  pathway = pathways,
  network = networks,
  sensitivity = c("low-ach", "high-ach"),
  intervention = c("Improved ventilation hcws", "Improved ventilation patients", "Mixed1", "Mixed2")
  )


params = bind_rows(params1, params_masks, params_ventilation) %>%
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
  

 write.table(params, file="param_grid_sensitivity.txt", row.names = F)
