### SIMULATION MODEL AVEC NEXTFLOW
library(dplyr)
library(lubridate)
library(tidyr)
library(purrr)
library(MASS)
library(Rcpp)

# On recupere la ligne avec le set de parametre
args <- commandArgs(trailingOnly = TRUE)
cat("\nParameters:\n")
cat(paste0(args, collapse = "\n"))
cat("\n")

# Modèle Rcpp
rlib=args[1] 
sourceCpp(rlib)
cat("\n")

# Fichier avec les codes pour l'analyse de la simulation
rlib2 = args[2]
source(rlib2)

# Other arguments
n_sim = as.numeric(args[4])
b_c_type <- args[5]
b_c <- as.numeric(args[5])
b_e_type <- gsub("/", "-", as.character(args[6]))
b_e <- eval(parse(text=args[6]))
network <- args[7]
threshold = 60
model = "linear"

# Load synthetic data
rinput = args[3] # dossier avec les fichiers input
data = paste0(rinput, "/parameters-synthetic-", network, "-", threshold, ".rda")
load(data)

# Set up intervention
intervention = args[9]
mu_int = mu
nu_int = nu 
rel_trans_risk = 1.0

if (intervention == "Hand hygiene") {
  rel_trans_risk = 1.0-0.17
}

if (intervention == "Symptomatic masking") {
  rel_trans_risk = 1.0-0.7
  nu_int = nu * (1-0.94) 
}

if (intervention == "Universal masking") {
  rel_trans_risk = 1.0-0.7
  nu_int = nu * (1-0.94) 
}

if (grepl("Improved ventilation", intervention)) {
  mu_int = (log(100)/(1/8)) * 24
}

if (intervention == "Mixed1") {
  mu_int = (log(100)/(1/8)) * 24
  nu_int = nu * (1-0.94)
  rel_trans_risk = 1.0-0.7
}

if (intervention == "Mixed2") {
  mu_int = (log(100)/(1/8)) * 24
  nu_int = nu * (1-0.94)
  rel_trans_risk = 1.0-0.7
}

# Sensitivity analysis
sensitivity = args[8]

if (sensitivity == "low-nu") nu = 1e3 * 2 * 24
if (sensitivity == "high-nu") nu = 1.8e7 * 2 * 24

if (sensitivity == "low-mu") mu = (log(2)/2.64) * 24
if (sensitivity == "high-mu") mu = (log(2)/0.64) * 24

if (sensitivity == "low-hh") rel_trans_risk = 1.0 - 0.05
if (sensitivity == "high-hh") rel_trans_risk = 1.0 - 0.24

if (sensitivity == "low-ach") mu_int = (log(100)/(1/2)) * 24
if (sensitivity == "high-ach") mu_int = (log(100)/(1/15)) * 24

if (sensitivity == "low-mw") {
  rel_trans_risk = 1.0 - 0.38
  nu_int = nu * (1.0-0.8)
}
if (sensitivity == "high-mw") {
  rel_trans_risk = 1.0 - 0.97
  nu_int = nu * (1.0-0.99) 
}

# Get ids by individual category 
id_paramedical = admission$id[admission$cat == "Paramedical"]
id_medical = admission$id[admission$cat == "Medical"]
id_patient = admission$id[admission$cat == "Patient"]

## RANDOM INDEX CASE
id_index <- sample(x = admission$id[admission$firstDate == as_date(floor_date(begin_date," day")) & 
                                      admission$firstDate != admission$lastDate], size = 1)
# id_index <- sample(x= admission %>% filter(id %in% admission$id) %>% distinct(id) %>% pull(), size = 1)

# Update global status (index is infectious when the simulation starts)
global_status_tmp <- global_status %>%
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

## SIMULATION USING RCPP
result <- simulation(global_interaction = global_interaction,
                     global_environment = global_environment,
                     global_data = global_data,
                     global_status = global_status_tmp,
                     beta_c = b_c,
                     beta_e = b_e,
                     B = B,
                     nu = nu,
                     mu = mu,
                     env_model = model,
                     dt = dt,
                     intervention = intervention,
                     mu_int = mu_int,
                     nu_int = nu_int,
                     rel_trans_risk = rel_trans_risk
                     )

# Save simulation summary statistics
stats = compute_SAR(result, id_paramedical, id_medical, id_patient) %>%
    mutate(
      Epidemic_duration = compute_epidemic_duration(result),
      Time_to_peak = get_time_to_peak_raw(result),
      network = network,
      model = model,
      threshold = threshold,
      beta_c = b_c_type,
      beta_e = b_e_type,
      intervention = intervention,
      sensitivity = sensitivity,
      nSim = n_sim
    )

write.csv2(stats, 
           paste0("summary_stat_", network, "_", b_c_type, "_", b_e_type, "_", intervention, "_", sensitivity, "_", n_sim, ".csv"), 
           row.names = F)


