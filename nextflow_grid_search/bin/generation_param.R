### GENERATION DES PARAMETRES ###
# args <- commandArgs(trailingOnly = TRUE)

library(dplyr)

#beta_c = c('0', '0.25', '0.5', '0.75', '1', '1.25', '1.5', '1.75', '2', '2.25', '2.5', '2.75', '3')
#beta_e = c('0', '1/120000', '1/110000', '1/100000', '1/90000', '1/80000', '1/70000', '1/60000', '1/50000', '1/40000', '1/30000', '1/28000', '1/26000', '1/24000', '1/22000', '1/20000', '1/18000','1/16000', '1/14000', '1/12000', '1/10000', '1/8000', '1/6000', '1/4000', '1/2000')
networks = c("herriot", "poincare")
# threshold = c(120,180) #60
# models = c("log-linear", "exponential")#"linear"


params = bind_rows( 
  expand.grid(
    scenario = 1:5,
    threshold = 60,
    network = networks,
    model = c("exponential", "log-linear")
    ),
  expand.grid(
    scenario = 1:5,
    threshold = c(120, 180),
    network = networks,
    model = "linear"
  ),
  )%>%
  mutate(
    beta_c = case_when(
      scenario == 1 ~ "0.5",
      scenario == 2 ~ "0.75",
      scenario == 3 ~ "1",
      scenario == 4 ~ "1.25",
      scenario == 5 ~ "1.5"
    ),
    beta_e = case_when(
      scenario == 1 ~ "1/20000",
      scenario == 2 ~ "1/26000",
      scenario == 3 ~ "1/30000",
      scenario == 4 ~ "1/40000",
      scenario == 5 ~ "1/60000"
    ),
    .before = threshold
  ) %>%
  select(-scenario)

write.table(params, file="param_grid.txt", row.names = F)