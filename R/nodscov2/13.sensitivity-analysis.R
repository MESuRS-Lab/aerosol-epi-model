################################################################################
##
##                            Sensitvity analyses
##
################################################################################

## Working environment----------------------------------------------------------
## R packages
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(gtsummary)
library(DescTools)
library(rstatix)
library(rcompanion)
library(performance)

## Load data--------------------------------------------------------------------
# Number of individuals
ind_numbers = data.frame()
load("out/parameters-synthetic-icu1-60.rda")
ind_numbers = bind_rows(ind_numbers, 
                        global_status %>% 
                          summarise(n_all = n(), n_patient = sum(grepl("^PA", id)), n_hcw = sum(grepl("^PE", id))) %>%
                          mutate(network = "ICU1"))
rm(list=setdiff(ls(), "ind_numbers"))

load("out/parameters-synthetic-icu2-60.rda")
ind_numbers = bind_rows(ind_numbers, 
                        global_status %>% 
                          summarise(n_all = n(), n_patient = sum(grepl("^PA", id)), n_hcw = sum(grepl("^PE", id))) %>%
                          mutate(network = "ICU2"))
rm(list=setdiff(ls(), "ind_numbers"))

# Load functions and dictionaries
source('R/nodscov2/helper-functions-simulations.R')
source("R/nodscov2/helper-functions.R")
source('R/nodscov2/dictionaries.R')

# Load epidemic data from sensitivity analyses
stats_df = read.csv2("nextflow_sensitivity/results/resu_sensitivity_all.txt", header = T) %>%
    mutate(sensitivity = case_when(
      sensitivity == "low-mu" ~ "high-sta",
      sensitivity == "high-mu" ~ "low-sta",
      .default = sensitivity)) %>% 
  mutate(
    beta_e = factor(gsub("-", "/", beta_e), c('1/150', '1/100', '1/70', '1/60', '1/45')),
    Pathway = case_when(
      beta_e == "1/45" & beta_c == 0.75 ~ "Pathway 1",
      beta_e == "1/60" & beta_c == 1 ~ "Pathway 2",
      beta_e == "1/70" & beta_c == 1.25 ~ "Pathway 3",
      beta_e == "1/100" & beta_c == 1.5 ~ "Pathway 4",
      beta_e == "1/150" & beta_c == 1.75 ~ "Pathway 5",
      .default = NA
    ),
    network = ifelse(network == "icu1", "ICU1", "ICU2"),
    intervention = factor(recode(intervention, !!!dict_interventions), dict_interventions),
    sensitivity = factor(sensitivity, dict_sensitivity),
    parameter_value = str_to_title(gsub("-.*", "", sensitivity)),
    sensitivity2 = case_when(
      grepl("-nu", sensitivity) ~ "Exhalation",
      grepl("-sta", sensitivity) ~ "Aerostability",
      grepl("-mw", sensitivity) ~ paste0("Masks\n", intervention),
      grepl("-hh", sensitivity) ~ intervention,
      grepl("-ach", sensitivity) ~ paste0("ACH\n", intervention),
    )
  ) 
head(stats_df)

# Data from main analysis
stats_df_main = read.csv2("nextflow_interventions/results/resu_interventions_all.txt", header = T) %>%
  mutate(
    beta_e = factor(gsub("-", "/", beta_e), c('1/150', '1/100', '1/70', '1/60', '1/45')),
    Pathway = case_when(
      beta_e == "1/45" & beta_c == 0.75 ~ "Pathway 1",
      beta_e == "1/60" & beta_c == 1 ~ "Pathway 2",
      beta_e == "1/70" & beta_c == 1.25 ~ "Pathway 3",
      beta_e == "1/100" & beta_c == 1.5 ~ "Pathway 4",
      beta_e == "1/150" & beta_c == 1.75 ~ "Pathway 5",
      .default = NA
    ),
    intervention = factor(recode(intervention, !!!dict_interventions), dict_interventions),
    network = toupper(network)
  ) 

# Compute median cumulative incidence in baseline scenario
median_none_main = stats_df_main %>%
  filter(intervention == "None") %>%
  mutate(I_hcw = I_medical + I_paramedical) %>%
  pivot_longer(c(I_global, I_patient, I_medical, I_paramedical, I_hcw), names_to = "Category", values_to = "I") %>%
  group_by(network, Pathway, Category) %>%
  summarise(ref = median(I), .groups = "drop") 

## Plot parameters of epidemic dynamics for all Pathways-----------------------
# Global SAR
stats_df %>%
  ggplot(., aes(x = Pathway, y = Global, fill = parameter_value)) +
  geom_violin(position = position_dodge(width = 0.7), width = 0.5) +
  facet_grid(cols = vars(network), rows = vars(sensitivity2)) +
  scale_fill_manual(values = sensitivity_pal) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = element_blank(),
        panel.grid.major = element_line(linewidth = 11/22, color = "grey90")) +
  labs(fill = "Parameter value", y = "Global SAR")

# SAR stratified by individual category
for (s in unique(stats_df$sensitivity2)) {
  p = stats_df %>%
    filter(sensitivity2 == s) %>%
    pivot_longer(c(Patient, Paramedical, Medical), names_to = "Category", values_to = "SAR") %>%
    ggplot(., aes(x = Pathway, y = SAR, fill = parameter_value)) +
    geom_violin() +
    # geom_boxplot(position = position_dodge(width = 0.7), width = 0.5, outliers = F) +
    facet_grid(rows = vars(network), cols = vars(Category)) +
    scale_fill_manual(values = sensitivity_pal) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
          panel.border = element_rect(linewidth = 1),
          axis.line = element_blank()) +
    labs(title = s, fill = "Parameter value", y = "SAR")
  print(p)
} 

# SAR stratified by source
for (s in unique(stats_df$sensitivity2)) {
  p = stats_df %>%
    filter(sensitivity2 == s) %>%
    pivot_longer(c(Contact, Environment), names_to = "Source", values_to = "SAR") %>%
    ggplot(., aes(x = Pathway, y = SAR, fill = parameter_value)) +
    geom_violin() +
    # geom_boxplot(position = position_dodge(width = 0.7), width = 0.5, outliers = F) +
    facet_grid(rows = vars(network), cols = vars(Source)) +
    scale_fill_manual(values = sensitivity_pal) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
          panel.border = element_rect(linewidth = 1),
          axis.line = element_blank()) +
    labs(title = s, fill = "Parameter value", y = "SAR")
  print(p)
} 

## Combined plot on epidemic dynamics-------------------------------------------
# Epidemic duration
p1 = stats_df %>%
  filter(!is.na(Epidemic_duration)) %>%
  ggplot(., aes(x = Pathway, y = Epidemic_duration, col = parameter_value)) +
  geom_violin(position = position_dodge(width = 0.7), width = 0.7) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.05, outliers = F) +
  facet_grid(cols = vars(network), rows = vars(sensitivity2)) +
  scale_color_manual(values = sensitivity_pal) +
  expand_limits(ymin = 0) +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
        panel.border = element_rect(linewidth = 1),
        axis.line = element_blank()) +
  labs(col = "Sensitivity analysis", y = "Epidemic duration (days)") 
p1
ggsave("fig/sensitivity_analysis/epidemic_duration.png", p1, height = 8, width = 7)

# Time to peak
p2 = stats_df %>%
  filter(!is.na(Time_to_peak)) %>%
  ggplot(., aes(x = Pathway, y = Time_to_peak, col = parameter_value)) +
  geom_violin(position = position_dodge(width = 0.7), width = 0.7) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.05, outliers = F) +
  facet_grid(cols = vars(network), rows = vars(sensitivity2)) +
  scale_color_manual(values = sensitivity_pal) +
  expand_limits(ymin = 0) +
  theme_classic() +
  theme(axis.title.x = element_blank(), 
        panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
        panel.border = element_rect(linewidth = 1),
        axis.line = element_blank()) +
  labs(col = "Sensitivity analysis", y = "Time to the peak (days)") 
p2
ggsave("fig/sensitivity_analysis/time_to_peak.png", p2, height = 8, width = 7)

# Extinction probability
p3 = stats_df %>%
  group_by(model, threshold, network, Pathway, sensitivity2, parameter_value) %>%
  summarise(p = sum(Global == 0) / n(), .groups = "drop") %>%
  ggplot(., aes(x = Pathway, y = p, fill = parameter_value)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.5) +
  facet_grid(cols = vars(network), rows = vars(sensitivity2)) +
  scale_fill_manual(values = sensitivity_pal) +
  expand_limits(ymin = 0) +
  theme_classic() +
  theme(axis.title.x = element_blank(), 
        panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
        panel.border = element_rect(linewidth = 1),
        axis.line = element_blank()) +
  labs(fill = "Sensitivity analysis", y = "Extinction probability") 
p3
ggsave("fig/sensitivity_analysis/extinction_proba.png", p3, height = 8, width = 7)

## Plot final figures for the effect of interventions---------------------------
# Global effect on SAR
bind_rows(
  stats_df,
  stats_df_main %>% filter(intervention %in% c("Hand hygiene")) %>% mutate(parameter_value = "Intermediate (main analysis)", sensitivity2 = "Hand hygiene"),
  stats_df_main %>% filter(intervention %in% c("Targeted masking", "Universal masking", "Targeted masking +\nVentilation patient rooms", "Targeted masking +\nVentilation HCW rooms")) %>% mutate(parameter_value = "Intermediate (main analysis)", sensitivity2 = paste0("Masks\n", intervention)),
  stats_df_main %>% filter(intervention %in% c("Ventilation HCW rooms", "Ventilation patient rooms", "Targeted masking +\nVentilation patient rooms", "Targeted masking +\nVentilation HCW rooms")) %>% mutate(parameter_value = "Intermediate (main analysis)", sensitivity2 = paste0("ACH\n", intervention))
  ) %>%
  filter(!sensitivity2 %in% c("Aerostability", "Exhalation")) %>%
  ggplot(., aes(x = Pathway, y = Global, fill = parameter_value)) +
  geom_violin(position = position_dodge(width = 0.7), width = 0.5) +
  facet_grid(cols = vars(network), rows = vars(sensitivity2)) +
  scale_fill_manual(values = sensitivity_pal) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), 
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1),
        axis.line = element_blank(),
        panel.grid.major = element_line(linewidth = 11/22, color = "grey90")) +
  labs(fill = "Parameter value", y = "Global SAR") 
ggsave("fig/sensitivity_analysis/sar_interventions.png", height = 14, width = 7)

# Hand hygiene
bind_rows(
  stats_df %>% filter(intervention == "Hand hygiene"),
  stats_df_main %>% filter(intervention != "Hand hygiene") %>% mutate(parameter_value = "Low"),
  stats_df_main %>% mutate(parameter_value = "Intermediate (main analysis)"),
  stats_df_main %>% filter(intervention != "Hand hygiene") %>% mutate(parameter_value = "High")
  ) %>%
  group_by(network, Pathway, parameter_value) %>%
  wilcox_test(I_global ~ intervention, ref.group = "None", paired = F, alternative = "two.sided", detailed = T) %>% 
  left_join(., median_none_main %>% filter(Category == "I_global"), by = c("network", "Pathway")) %>%
  mutate(
    estimate = -estimate/ref,
    conf.high = -conf.high/ref,
    conf.low = -conf.low/ref,
    group2 = factor(group2, dict_interventions)
  ) %>%
  nest(.by = c(network, Pathway, parameter_value)) %>%
  mutate(out = map(data, rank_interventions)) %>%
  select(-data) %>%
  unnest(out) %>%
  ggplot(., aes(x = Pathway, y = group2, fill = estimate)) +
  geom_tile() +
  geom_text(aes(label = rank), fontface = "bold", color = "white") +
  facet_grid(rows = vars(parameter_value), cols = vars(network)) +
  scale_fill_viridis(labels = scales::percent) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(), axis.title.y = element_blank(), 
    legend.position = "bottom", legend.title.position = "top", 
    legend.title = element_text(hjust = 0.5),
    panel.border = element_rect(linewidth = 1),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  labs(fill = "Relative difference of the cumulative incidence in all individuals")
ggsave("fig/sensitivity_analysis/rank_sensitivity_hand-hygiene.png", height = 8, width = 6)

# Mask wearing
bind_rows(
  stats_df %>% filter(grepl("-mw$", sensitivity)),
  stats_df_main %>% filter(intervention %in% c("None", "Ventilation patient rooms", "Ventilation HCW rooms", "Hand hygiene")) %>% mutate(parameter_value = "Low"),
  stats_df_main %>% mutate(parameter_value = "Intermediate (main analysis)"),
  stats_df_main %>% filter(intervention  %in% c("None", "Ventilation patient rooms", "Ventilation HCW rooms", "Hand hygiene")) %>% mutate(parameter_value = "High")
) %>%
  group_by(network, Pathway, parameter_value) %>%
  wilcox_test(I_global ~ intervention, ref.group = "None", paired = F, alternative = "two.sided", detailed = T) %>% 
  left_join(., median_none_main %>% filter(Category == "I_global"), by = c("network", "Pathway")) %>%
  mutate(
    estimate = -estimate/ref,
    conf.high = -conf.high/ref,
    conf.low = -conf.low/ref,
    group2 = factor(group2, dict_interventions)
  ) %>%
  nest(.by = c(network, Pathway, parameter_value)) %>%
  mutate(out = map(data, rank_interventions)) %>%
  select(-data) %>%
  unnest(out) %>%
  ggplot(., aes(x = Pathway, y = group2, fill = estimate)) +
  geom_tile() +
  geom_text(aes(label = rank), fontface = "bold", color = "white") +
  facet_grid(rows = vars(parameter_value), cols = vars(network)) +
  scale_fill_viridis(labels = scales::percent) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(), axis.title.y = element_blank(), 
    legend.position = "bottom", legend.title.position = "top", 
    legend.title = element_text(hjust = 0.5),
    panel.border = element_rect(linewidth = 1),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  labs(fill = "Relative difference of the cumulative incidence in all individuals")
ggsave("fig/sensitivity_analysis/rank_sensitivity_mask-wearing.png", height = 8, width = 6)

# Ventilation
bind_rows(
  stats_df %>% filter(grepl("ach$", sensitivity)),
  stats_df_main %>% filter(intervention %in% c("None", "Hand hygiene", "Universal masking", "Targeted masking")) %>% mutate(parameter_value = "Low"),
  stats_df_main %>% mutate(parameter_value = "Intermediate (main analysis)"),
  stats_df_main %>% filter(intervention %in% c("None", "Hand hygiene", "Universal masking", "Targeted masking")) %>% mutate(parameter_value = "High")
) %>%
  group_by(network, Pathway, parameter_value) %>%
  wilcox_test(I_global ~ intervention, ref.group = "None", paired = F, alternative = "two.sided", detailed = T) %>% 
  left_join(., median_none_main %>% filter(Category == "I_global"), by = c("network", "Pathway")) %>%
  mutate(
    estimate = -estimate/ref,
    conf.high = -conf.high/ref,
    conf.low = -conf.low/ref,
    group2 = factor(group2, dict_interventions)
  ) %>%
  nest(.by = c(network, Pathway, parameter_value)) %>%
  mutate(out = map(data, rank_interventions)) %>%
  select(-data) %>%
  unnest(out) %>%
  ggplot(., aes(x = Pathway, y = group2, fill = estimate)) +
  geom_tile() +
  geom_text(aes(label = rank), fontface = "bold", color = "white") +
  facet_grid(rows = vars(parameter_value), cols = vars(network)) +
  scale_fill_viridis(labels = scales::percent) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(), axis.title.y = element_blank(), 
    legend.position = "bottom", legend.title.position = "top", 
    legend.title = element_text(hjust = 0.5),
    panel.border = element_rect(linewidth = 1),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  labs(fill = "Relative difference of the cumulative incidence in all individuals")
ggsave("fig/sensitivity_analysis/rank_sensitivity_ventilation.png", height = 8, width = 6)

## Plot final figures for the effect of exhalation and aerostability------------
# Baseline Pathways
bind_rows(
  stats_df %>%filter(sensitivity2 %in% c("Aerostability", "Exhalation")), 
  stats_df_main %>% filter(intervention == "None") %>% mutate(parameter_value = "Intermediate (main analysis)", sensitivity2 = "Aerostability"),
  stats_df_main %>% filter(intervention == "None") %>% mutate(parameter_value = "Intermediate (main analysis)", sensitivity2 = "Exhalation")
  ) %>%
  ggplot(., aes(x = Pathway, y = Global, fill = parameter_value)) +
  geom_violin(position = position_dodge(width = 0.7), width = 0.5) +
  facet_grid(cols = vars(network), rows = vars(sensitivity2)) +
  scale_fill_manual(values = sensitivity_pal) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), 
        axis.title.x = element_blank(),
        panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
        panel.border = element_rect(linewidth = 1),
        axis.line = element_blank()) +
  labs(fill = "Parameter value", y = "Global SAR") 
ggsave("fig/paper/figure6.png", height = 4, width = 7)

