################################################################################
##                Dictionaries, color parlettes and variables 
##                      to analyze interaction data
################################################################################

# Dictionaries------------------------------------------------------------------
dict_cat = c("aux nurse" = "Paramedical", 
             "nurse" = "Paramedical",
             "student nurse" = "Paramedical",
             "reeducation staff" = "Paramedical",
             "ext physician" = "Medical",
             "physician" = "Medical",
             "patient" = "Patient")

dict_cat_initial = c("administration" = "Other",
                     "aux nurse" = "HCW", 
                     "ext physician" = "HCW",
                     "investigation" = "Other",
                     "logistic" = "Other",
                     "nurse" = "HCW",
                     "other" = "Other",
                     "physician" = "HCW",
                     "reeducation staff" = "HCW",
                     "student nurse" = "HCW",
                     "visitor" = "Visitor"
)

dict_cat_num = c(
  "1" = "nurse",
  "2" = "aux nurse",
  "3" = "reeducation staff",
  "4" = "physician",
  "5" = "ext physician",
  "6" = "administration",
  "7" = "logistic",
  "8" = "investigation",
  "9" = "patient",
  "10" = "visitor",
  "11" = "other",
  "12" = "student nurse"
)

dict_scenarios = c("sim_1-4_20" = "Scenario 1",
                   "sim_1-2_16" = "Scenario 2",
                   "sim_3-4_12" = "Scenario 3",
                   "sim_1_8" = "Scenario 4",
                   "sim_3-2_5" = "Scenario 5"
)

dict_interventions = c(
  "None" = "None", 
  "Hand hygiene" = "Hand hygiene", 
  "Symptomatic masking" = "Targeted masking", 
  "Universal masking" = "Universal masking", 
  "Improved ventilation patients" = "Ventilation patient rooms", 
  "Improved ventilation hcws" = "Ventilation HCW rooms", 
  "Mixed1" = "Targeted masking +\nVentilation patient rooms", 
  "Mixed2" = "Targeted masking +\nVentilation HCW rooms"
)

dict_sensitivity = c("low-nu", "intermediate-nu", "high-nu",
                     "low-sta", "intermediate-sta", "high-sta",
                     "low-mw", "intermediate-mw", "high-mw",
                     "low-ach", "intermediate-ach", "high-ach",
                     "low-hh", "intermediate-hh", "high-hh")

# Plots-------------------------------------------------------------------------
pal = c('Medical' = "#5CD6D6", 'Paramedical' = "#A9B9E8", 'Patient' = "#FFA766", 'HCW' = "grey50", "All" = "darkblue", 'Room' = "#666699")
network_pal = c('ICU1' = "darkcyan", "ICU2" = "chocolate4")
env_pal = c('Short-range' = 'darkorange', 'Long-range' =  'orchid')
algo_syn_pal = c("Observed" = "darkorchid", "Reconstructed"= "darkorange")
model_pal = c("exponential" = "cadetblue4", "linear" = "goldenrod1", "log-linear" = "cadetblue1")
threshold_pal = c("60" = "#FC766AFF", "120" = "#B0B8B4FF", "180" = "#184A44FF")
intervention_pal = c("None" = "grey90", 
                     "Hand hygiene" = "darkmagenta", 
                     "Universal masking" = "turquoise3", 
                     "Targeted masking" = "cadetblue1", 
                     "Ventilation patient rooms" = "darkorange", 
                     "Ventilation HCW rooms" = "gold", 
                     "Targeted masking +\nVentilation patient rooms" = "deeppink3",
                     "Targeted masking +\nVentilation HCW rooms" = "pink"
                     )
sensitivity_pal = c("High" = "#2c7fb8", "Intermediate (main analysis)" = "#7fcdbb", "Low" = "#edf8b1")
pairs_pal = c(
  "Patient to Patient" = "#95a757",
  "Patient to HCW" = "#4b662d",
  "HCW to Patient" = "#FEFAE0",
  "HCW to HCW" = "#DDA15E"
)

# Variables---------------------------------------------------------------------
noon_day1 = as_datetime("2020-05-06 12:00:00")
midnight_day1 = as_datetime("2020-05-06 00:00:00")
noon_day2 = as_datetime("2020-05-07 12:00:00")
midnight_day2 = as_datetime("2020-05-07 00:00:00")
noon_last_day = noon_day1 + 90*3600*24
midnight_last_day = midnight_day1 + 90*3600*24