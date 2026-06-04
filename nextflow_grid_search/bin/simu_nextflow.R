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
b_c_type <- gsub("/", "-", as.character(args[5]))
b_c <- eval(parse(text = args[5]))
b_e_type <- gsub("/", "-", as.character(args[6]))
b_e <- eval(parse(text=args[6]))
threshold = as.numeric(args[7])
network <- args[8]
model = args[9]

# Load synthetic data
rinput = args[3] # dossier avec les fichiers input
data = paste0(rinput, "/parameters-synthetic-", network, "-", threshold, ".rda")
load(data)

# Get ids by individual category 
id_paramedical = admission$id[admission$cat == "Paramedical"]
id_medical = admission$id[admission$cat == "Medical"]
id_patient = admission$id[admission$cat == "Patient"]

## RANDOM INDEX CASE
firstDay = as_date(floor_date(begin_date, "day"))
stats = data.frame()

for (i in 1:n_sim) {
  id_index <- sample(x = admission$id[admission$firstDate == firstDay & admission$firstDate != admission$lastDate], size = 1)
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
                       dt = dt)
  
  # assign(paste0("sim_", sim_id), result)
  # save_path <- file.path('grid-search', model, threshold, network, paste0('sim_', beta_c_type, '_', beta_e_type ))
  # if (!dir.exists(save_path)) dir.create(save_path, recursive = T, showWarnings = F)
  # 
  # # Save simulation results
  # save(list = paste0("sim_", sim_id), 
  #      file = file.path(save_path, paste0('sim_', sim_id, ".rda")))
  
  # Save simulation summary statistics
  stats = bind_rows(
    stats,
    compute_SAR(result, id_paramedical, id_medical, id_patient) %>%
      mutate(
        Epidemic_duration = compute_epidemic_duration(result),
        Time_to_peak = get_time_to_peak_raw(result),
        network = network,
        model = model,
        threshold = threshold,
        beta_c = b_c,
        beta_e = b_e,
        nSim = i
      )
  )
  
  # Print basic information
  cat(paste0("Simulation ", i, ": ", id_index, " - ", sum(!result$inf_by %in% c("", "INDEX")), "\n"))
}

write.csv2(stats, 
           paste0("summary_stat_", model, "_", threshold, "_", network, "_", b_c_type, "_", b_e_type, ".csv"), 
           row.names = F)


