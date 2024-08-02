######################
## READ RDATA FILES ##
######################

## READ RData FILES -> STRUCTURE OF FILE: sim_<beta>_<nu>_<sim_id>.RData 
## WARNING BETA = 1/3 WILL BE WRITTEN AS 1-3
## IN THE RDATA FILE, THERE WILL BE A DATAFRAME NAMED: sim_<beta>_<nu>_<sim_id>
load_rdata_to_list <- function(file, list_sim, dir) {
  load(file.path(intervention_path, 'out', dir, 'results', file))
  file_name <- basename(file)
  parts <- strsplit(file_name, "_")[[1]]
  sim_id <- parts[length(parts)] ## last number is the sim id
  sim_id <- as.integer(sub("\\.RData$", "", sim_id))  # Remove file extension
  # beta <- eval(parse(text = gsub(pattern = '-', replacement = '/', x = parts[2])))
  # nu <- eval(parse(text = gsub(pattern = '-', replacement = '/', x = parts[3])))
  sim_object <- get(sub("\\.RData$", "", file_name))
  list_sim[[dir]][[sim_id]] <<- sim_object
}


#####################
## SAR COMPUTATION ##
#####################

##################### SAR FOR ONE SIMULATION
compute_SAR <- function(global_status) {
  ## WE NEED INDEX TYPE FOR THE SAR BY TYPE COMPUTATION 
  index_id <- global_status %>% filter(inf_by == 'INDEX') %>% pull(id)
  index_type <- case_when(
    index_id %in% id_paramedical ~ 'PARAMEDICAL',
    index_id %in% id_medical ~ 'MEDICAL',
    index_id %in% id_patient ~ 'PATIENT'
  )
  ## SUSCEPTIBLE INDIVIDUALS DURING THE STUDY
  n_individual <- length(global_status %>% distinct(id) %>% pull()) -1
  n_medical <- ifelse(index_type == 'MEDICAL', length(id_medical) - 1, length(id_medical))
  n_paramedical <- ifelse(index_type == 'PARAMEDICAL', length(id_paramedical) - 1, length(id_paramedical))
  n_patient <- ifelse(index_type == 'PATIENT', length(id_patient) - 1, length(id_patient))
  n_hcw <- ifelse(index_type == 'PATIENT', length(id_hcw) - 1, length(id_hcw) - 1)
  
  ## SAR BY TYPE
  SAR_global <- (length(global_status %>% filter(t_inf != -1) %>% pull(id)) -1 )/(n_individual)
  SAR_patient <- length(global_status %>% filter(id %in% id_patient & t_inf != -1 & inf_by != 'INDEX') %>% pull(id)) / (n_patient)
  SAR_medical <- length(global_status %>% filter(id %in% id_medical & t_inf != -1 & inf_by != 'INDEX') %>% pull(id)) / n_medical
  SAR_paramedical <- length(global_status %>% filter(id %in% id_paramedical & t_inf != -1 & inf_by != 'INDEX') %>% pull(id)) / n_paramedical
  SAR_hcw <- length(global_status %>% filter(id %in% id_hcw & t_inf != -1 & inf_by != 'INDEX') %>% pull(id)) / n_hcw
  ## CLOSE CONTACT / ENVIRONMENT (SUM EQUALS TO GLOBAL)
  SAR_environment <- length(global_status %>% filter(grepl('ENVIRONMENT', inf_by)) %>% pull(id)) / n_individual
  SAR_contact <- length(global_status %>% filter(grepl('CONTACT', inf_by)) %>% pull(id)) / n_individual
  
  SAR_patient_e <- length(global_status %>% filter(id %in% id_patient & t_inf != -1 & grepl(pattern = 'ENVIRONMENT', inf_by)) %>% pull(id)) / (n_patient)
  SAR_patient_c <- length(global_status %>% filter(id %in% id_patient & t_inf != -1 & grepl(pattern = 'CONTACT', inf_by)) %>% pull(id)) / (n_patient)
  SAR_medical_e <- length(global_status %>% filter(id %in% id_medical & t_inf != -1 & grepl(pattern = 'ENVIRONMENT', inf_by)) %>% pull(id)) / n_medical
  SAR_medical_c <- length(global_status %>% filter(id %in% id_medical & t_inf != -1 & grepl(pattern = 'CONTACT', inf_by)) %>% pull(id)) / n_medical
  SAR_paramedical_e <- length(global_status %>% filter(id %in% id_paramedical & t_inf != -1 & grepl(pattern = 'ENVIRONMENT', inf_by)) %>% pull(id)) / n_paramedical
  SAR_paramedical_c <- length(global_status %>% filter(id %in% id_paramedical & t_inf != -1 & grepl(pattern = 'CONTACT', inf_by)) %>% pull(id)) / n_paramedical
  ##DF
  SAR_df <- data.frame(
    Global = SAR_global,
    HCW = SAR_hcw,
    Patient = SAR_patient,
    Paramedical = SAR_paramedical,
    Medical = SAR_medical,
    Environment = SAR_environment,
    Contact = SAR_contact,
    Patient_Contact = SAR_patient_c,
    Patient_Environment = SAR_patient_e,
    Paramedical_Contact = SAR_paramedical_c,
    Paramedical_Environment = SAR_paramedical_e,
    Medical_Contact = SAR_medical_c,
    Medical_Environment = SAR_medical_e)
  
  return(SAR_df)
}

##################### GET ALL SAR FOR ONE COUPLE (N SIMULATIONS)
get_SAR_couple <- function(couple, list_sim){
  lapply(list_sim[[couple]], function(global_status) {
    compute_SAR(global_status = global_status)
  })
}

##################### GET SAR FOR ALL BETA/NU COUPLES
get_all_SAR <- function(list_sim, id_patient, id_hcw, id_paramedical, id_medical) {
  list_SAR <- list()
  for(couple in names(list_sim)){
    list_SAR[[couple]] <- get_SAR_couple(couple, list_sim)
  }
  return(list_SAR)
}


#################
## SAR METRICS ##
#################

################# SAR METRICS FOR ONE COUPLE
get_SAR_g_metrics <- function(couple, list_SAR){
  couple_parts <- strsplit(couple, "_")[[1]]
  SAR_df <- rbindlist(list_SAR[[couple]])
  SAR_df_metrics <- data.frame(couple = couple,
                               beta = eval(parse(text = gsub(pattern = '-', replacement = '/', x = couple_parts[2]))),
                               nu = eval(parse(text = gsub(pattern = '-', replacement = '/', x = couple_parts[3]))),
                               SAR_min = min(SAR_df$Global),
                               SAR_max = max(SAR_df$Global),
                               SAR_median = median(SAR_df$Global),
                               SAR_mean = mean(SAR_df$Global),
                               SAR_sd = sd(SAR_df$Global))
  return(SAR_df_metrics)
}

################# SAR METRICS FOR EVERY COUPLES
get_all_SAR_g_metrics <- function(list_SAR){
  all_SAR_metrics <- data.frame()
  for(couple in names(list_SAR)){
    all_SAR_g_metrics <- rbind(all_SAR_metrics, get_SAR_g_metrics(couple = couple, list_SAR = list_SAR) )
  }
  return(all_SAR_g_metrics)
  
}

compute_all_SAR_metrics <- function(list_SAR) {
  all_metrics <- list()
  for (couple in names(list_SAR)) {
    SAR_dfs <- list_SAR[[couple]]
    
    SAR_couple <- rbindlist(SAR_dfs)
    
    SAR_metrics <- data.frame(
      couple = couple,
      Global_min = min(SAR_couple$Global),
      Global_max = max(SAR_couple$Global),
      Global_mean = mean(SAR_couple$Global),
      Global_median = median(SAR_couple$Global),
      Global_sd = sd(SAR_couple$Global),
      
      HCW_min = min(SAR_couple$HCW),
      HCW_max = max(SAR_couple$HCW),
      HCW_mean = mean(SAR_couple$HCW),
      HCW_median = median(SAR_couple$HCW),
      HCW_sd = sd(SAR_couple$HCW),
      
      Patient_min = min(SAR_couple$Patient),
      Patient_max = max(SAR_couple$Patient),
      Patient_mean = mean(SAR_couple$Patient),
      Patient_median = median(SAR_couple$Patient),
      Patient_sd = sd(SAR_couple$Patient),
      
      Paramedical_min = min(SAR_couple$Paramedical),
      Paramedical_max = max(SAR_couple$Paramedical),
      Paramedical_mean = mean(SAR_couple$Paramedical),
      Paramedical_median = median(SAR_couple$Paramedical),
      Paramedical_sd = sd(SAR_couple$Paramedical),
      
      Medical_min = min(SAR_couple$Medical),
      Medical_max = max(SAR_couple$Medical),
      Medical_mean = mean(SAR_couple$Medical),
      Medical_median = median(SAR_couple$Medical),
      Medical_sd = sd(SAR_couple$Medical),
      
      Environment_min = min(SAR_couple$Environment),
      Environment_max = max(SAR_couple$Environment),
      Environment_mean = mean(SAR_couple$Environment),
      Environment_median = median(SAR_couple$Environment),
      Environment_sd = sd(SAR_couple$Environment),
      
      Contact_min = min(SAR_couple$Contact),
      Contact_max = max(SAR_couple$Contact),
      Contact_mean = mean(SAR_couple$Contact),
      Contact_median = median(SAR_couple$Contact),
      Contact_sd = sd(SAR_couple$Contact)
    )
    
    # Append the metrics data frame to the list
    all_metrics[[couple]] <- SAR_metrics
  }
  
  # Combine all metrics into one data frame
  all_metrics_df <- do.call(rbind, all_metrics)
  
  return(all_metrics_df)
}

########################
## EPIDEMIC DURATIONS ##
########################

######################## FOR ONE SIMULATION
compute_epidemic_duration <- function(global_status) {
  ##LAST INFECTIOUS INDIVIDUAL TRANSITION TO RECOVERED
  n_time_step <- max(global_status$t_recover)
  ## DURATION IN DAYS
  epidemic_duration <- (n_time_step *30)/(3600*24)
  return(epidemic_duration)
}


######################## FOR ONE COUPLE (N SIMULATIONS)
get_epidemic_duration_couple <- function(couple, list_sim){
  epidemic_duration_df <- data.frame()
  results <- lapply(seq_along(list_sim[[couple]]), function(i) {
    duration <- compute_epidemic_duration(global_status = list_sim[[couple]][[i]])
    epidemic_duration_df <- rbind(epidemic_duration_df, data.frame(couple = couple, id_sim = i, duration = duration))
  })
}

get_all_epidemic_duration <- function(list_sim){
  list_epidemic_duration <- list()
  for(couple in names(list_sim)){
    list_epidemic_duration[[couple]] <-  get_epidemic_duration_couple(couple = couple, list_sim = list_sim)
  }
  return(list_epidemic_duration)
}


###############################
## EPIDEMIC DURATIONS METRICS##
###############################

################# EPIDEMIC DURATIONS METRICS FOR ONE COUPLE
get_epidemic_duration_metrics <- function(couple, list_epidemic_duration){
  couple_parts <- strsplit(couple, "_")[[1]]
  ED_df <- rbindlist(list_epidemic_duration[[couple]])
  ED_df_metrics <- data.frame(couple = couple,
                               beta = eval(parse(text = gsub(pattern = '-', replacement = '/', x = couple_parts[2]))),
                               nu = eval(parse(text = gsub(pattern = '-', replacement = '/', x = couple_parts[3]))),
                               ED_min = min(ED_df$duration),
                               ED_max = max(ED_df$duration),
                               ED_median = median(ED_df$duration),
                               ED_mean = mean(ED_df$duration),
                               ED_sd = sd(ED_df$duration))
  return(ED_df_metrics)
}

################# EPIDEMIC DURATIONS METRICS FOR EVERY COUPLES
get_all_epidemic_duration_metrics <- function(list_epidemic_duration){
  all_ED_metrics <- data.frame()
  for(couple in names(list_epidemic_duration)){
    all_ED_metrics <- rbind(all_ED_metrics, get_epidemic_duration_metrics(couple, list_epidemic_duration) )
  }
  return(all_ED_metrics)
  
}

#####################################
## MAX N INFECTED DURING EPIDEMIC  ##
#####################################


##########
## SEIR ##
##########

compute_SEIR <- function(global_status, n_subdivisions) {
  n_individual <- length(global_status %>% distinct(id) %>% pull())
  ## Time steps
  times <- seq(1, n_subdivisions, 100)
  ##DF
  counts_df <- data.frame(
    Susceptible = colSums(outer(global_status$t_inf, times, function(x, y) x == -1 | x > y)),
    Exposed = colSums(outer(global_status$t_inf, times, function(x, y) x <= y) & outer(global_status$t_incub, times, function(x, y) x > y)),
    Infectious = colSums(outer(global_status$t_incub, times, function(x, y) x <= y) & outer(global_status$t_recover, times, function(x, y) x > y)),
    Recovered = colSums(outer(global_status$t_recover, times, function(x, y) x <= y) & outer(global_status$t_inf, times, function(x, y) x != -1)),
    time = times
  )
  ## Pivot for easier plot
  counts_df <- counts_df %>%
    pivot_longer(cols = c("Susceptible", "Exposed", "Infectious", "Recovered"), names_to = "status", values_to = "count") %>%
    mutate(proportion = count * 100 / n_individual)
  ## ## S E I R as factor 
  counts_df$status <- factor(counts_df$status, levels = c("Susceptible", "Exposed", "Infectious", "Recovered"))
  
  return(counts_df)
}


##################### GET ALL SEIR FOR ONE COUPLE (N SIMULATIONS)
get_SEIR_couple <- function(couple, list_sim, n_subdivisions){
  lapply(list_sim[[couple]], function(global_status) {
    compute_SEIR(global_status = global_status, n_subdivisions = n_subdivisions)
  })
}

##################### GET SEIR FOR ALL BETA/NU COUPLES
get_all_SEIR <- function(list_sim, n_subdivisions) {
  list_SEIR <- list()
  for(couple in names(list_sim)){
    list_SEIR[[couple]] <- get_SEIR_couple(couple = couple, list_sim = list_sim, n_subdivisions = n_subdivisions)
  }
  return(list_SEIR)
}
