################################################################################
##
##          Calculation of R0 for the different transmission schemes
##
################################################################################

## Working environment----------------------------------------------------------
## R Packages
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(gt)
library(Rcpp)
library(foreach)
library(doParallel)

## Compile Rcpp code
sourceCpp("cpp/model-r0.cpp")

## Model parameterisation-------------------------------------------------------
n_sim = 100
threshold = 60
model = "linear"
schemes = data.frame(
  Scheme = paste("Scheme", 1:5),
  beta_e = c('1/45', '1/60', '1/70', '1/100', '1/150'),
  beta_c = c(0.75, 1, 1.25, 1.5, 1.75)
)

input_data = crossing(schemes, expand.grid(nSim = 1:n_sim, network = c("icu1", "icu2")))

## Load input data--------------------------------------------------------------
network_data = vector("list", 2)
names(network_data) = c("icu1", "icu2")

for (network in c("icu1", "icu2")) {
  
  load(paste0("out/parameters-synthetic-", network, "-", threshold, ".rda"))
  
  # Trim input data to 2880*7 + 1
  global_data = global_data[1:(2880*7 + 1)]
  global_environment = global_environment[1:(2880*7 + 1)]
  global_interaction = global_interaction[1:(2880*7 + 1)]
  
  firstDay = as_date(floor_date(begin_date, "day"))
  admission = admission %>% filter(firstDate <= firstDay+8)
  global_status = global_status %>% filter(id %in% admission$id)
  
  # Get individual IDs
  id_paramedical = admission$id[admission$cat == "Paramedical"]
  id_medical = admission$id[admission$cat == "Medical"]
  id_patient = admission$id[admission$cat == "Patient"] 
  
  # Save into final object
  network_data[[network]] = list(
    "admission" = admission, "global_data" = global_data, 
    "global_environment" = global_environment, "global_interaction" = global_interaction, 
    "global_status" = global_status, "B" = B, "firstDay" = firstDay, 
    "beta_c" = beta_c, "beta_e" = beta_e, "dt" = dt, "model" = model, 
    "mu" = mu, "nu" = nu, "id_medical" = id_medical, "id_paramedical" = id_paramedical, 
    "id_patient" = id_patient
    )
  
  # Erase unique objects
  rm(list = c("admission", "global_data", "global_environment", "global_interaction", 
              "global_status", "B", "begin_date", "firstDay", "beta_c", "beta_e",
              "dt", "end_date", "env_model", "mu", "n_days", "n_subdivisions", 
              "nu", "t_begin", "t_end", "id_medical", "id_paramedical", "id_patient"))
}

## Launch simulations-----------------------------------------------------------
set.seed(20251013)
registerDoParallel(2)

start = Sys.time()
r0 = foreach (r=1:nrow(input_data), .combine=bind_rows) %dopar% {
  
  # Arguments
  network = input_data$network[r]
  s = input_data$nSim[r]
  admission = network_data[[network]]$admission
  firstDay = network_data[[network]]$firstDay
  
  # Randomly select index case 
  id_index <- sample(x = admission$id[
    admission$firstDate == firstDay & admission$firstDate != admission$lastDate
    ], size = 1)
  
  # Update global status (index is infectious when the simulation starts)
  global_status_tmp <- network_data[[network]]$global_status %>%
    mutate(t_inf = ifelse(id == id_index,
                          as.integer(1),
                          t_inf),
           t_incub = ifelse(id == id_index,
                            as.integer(1),
                            t_recover),
           t_recover = ifelse(id == id_index,
                              as.integer((t_incub)  + runif(1, min = 2880*3, max = 2880*7)),
                              t_recover),
           inf_by = ifelse(id == id_index,
                           "INDEX",
                           inf_by))
  
  # Simulate transmission
  result <- simulation(
    global_interaction = network_data[[network]]$global_interaction,
    global_environment = network_data[[network]]$global_environment,
    global_data = network_data[[network]]$global_data,
    global_status = global_status_tmp,
    beta_c = input_data$beta_c[r],
    beta_e = eval(parse(text = input_data$beta_e[r])),
    B = network_data[[network]]$B,
    nu = network_data[[network]]$nu,
    mu = network_data[[network]]$mu,
    env_model = network_data[[network]]$model,
    dt = network_data[[network]]$dt,
    intervention = "None"
  )
    
  # Extract R0 and category of the index case 
  r0_tmp = result %>%
    mutate(
      network = network,
      Scheme = input_data$Scheme[r],
      nSim = s
    ) %>%
    filter(inf_by != "") %>%
    group_by(network, Scheme, nSim) %>%
    summarise(
      R0 = n(), 
      index = .$id[.$inf_by == "INDEX"],
      index = case_when(index %in% network_data[[network]]$id_patient ~ "Patient", 
                        index %in% network_data[[network]]$id_medical ~ "Medical staff",
                        index %in% network_data[[network]]$id_paramedical ~ "Paramedical staff"
      ),
      .groups = "drop"
    ) 
  r0_tmp
}
end = Sys.time()
end-start

## Reprensent correspondance----------------------------------------------------
# Summary statistics
r0_ss = r0 %>% 
  group_by(network, Scheme, index) %>%
  summarise(
    R0_mid = median(R0),
    R0_low = quantile(R0, 0.025),
    R0_up = quantile(R0, 0.975),
    .groups = "drop"
  )

# Plot
ggplot(r0, aes(x = Scheme, y = R0, col = index)) + 
  geom_violin() +
  geom_pointrange(data = r0_ss, aes(y = R0_mid, ymin = R0_low, ymax = R0_up), 
                  position = position_dodge(width = 0.9)) +
  facet_grid(cols = vars(network)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  expand_limits(y = 0) +
  labs(x = "", col = "Index case", y = "R0 (median, 2.5 and 97.5 percentiles)")

# Table
r0_ss %>%
  mutate(
    R0 = paste0(R0_mid, " (", R0_low, "-", R0_up, ")"),
    network = toupper(network)
    ) %>%
  select(-matches("R0_")) %>%
  pivot_wider(names_from = index, values_from = R0) %>%
  group_by(network) %>%
  gt(.) %>%
  tab_style(
    style = "vertical-align:middle; font-weight:bold",
    locations = list(cells_column_labels(), cells_row_groups())
  )


