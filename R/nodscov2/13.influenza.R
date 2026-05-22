################################################################################
##
##                          Intervention evaluation
##
################################################################################

## Working environment----------------------------------------------------------
## R packages
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(ggpattern)
library(gtsummary)
library(DescTools)
library(rstatix)
library(rcompanion)
library(performance)
library(viridis)

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
source('R/nodscov2/helper-functions.R')
source('R/nodscov2/dictionaries.R')

# All conditions
stats_df = read.csv2("nextflow_influenza/results/resu_influenza_all.txt", header = T) %>%
  mutate(
    beta_e = factor(gsub("-", "/", beta_e), c('1/150', '1/100', '1/70', '1/60', '1/45')),
    Scheme = case_when(
      beta_e == "1/45" & beta_c == 0.75 ~ "Scheme 1",
      beta_e == "1/60" & beta_c == 1 ~ "Scheme 2",
      beta_e == "1/70" & beta_c == 1.25 ~ "Scheme 3",
      beta_e == "1/100" & beta_c == 1.5 ~ "Scheme 4",
      beta_e == "1/150" & beta_c == 1.75 ~ "Scheme 5",
      .default = NA
    ),
    intervention = factor(recode(intervention, !!!dict_interventions), dict_interventions),
    network = toupper(network)
  ) 
head(stats_df)
stats_df %>% count(intervention)

## Plot parameters of epidemic dynamics for all scenarios-----------------------
# Global SAR
stats_df %>%
  ggplot(., aes(x = Scheme, y = Global, fill = intervention)) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.5) +
  facet_grid(rows = vars(network)) +
  scale_fill_manual(values = intervention_pal) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  labs(fill = "Intervention", y = "SAR")

# SAR stratified by individual category
stats_df %>%
  mutate(intervention = factor(intervention, dict_interventions)) %>%
  pivot_longer(c(Patient, Paramedical, Medical), names_to = "Category", values_to = "SAR") %>%
  ggplot(., aes(x = Scheme, y = SAR, fill = intervention)) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.5, outliers = F) +
  facet_grid(rows = vars(network), cols = vars(Category)) +
  scale_fill_manual(values = intervention_pal) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  labs(fill = "Intervention", y = "SAR")

# SAR stratified by source
stats_df %>%
  mutate(intervention = factor(intervention, dict_interventions)) %>%
  pivot_longer(c(Contact, Environment), names_to = "Source", values_to = "SAR") %>%
  ggplot(., aes(x = Scheme, y = SAR, fill = intervention)) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.5, outliers = F) +
  facet_grid(rows = vars(network), cols = vars(Source)) +
  theme_bw() +
  scale_fill_manual(values = intervention_pal) +
  theme(axis.title.x = element_blank()) +
  labs(fill = "Intervention", y = "SAR")

## Combined plot on epidemic dynamics-------------------------------------------
# Epidemic duration
p1 = stats_df %>%
  mutate(intervention = factor(intervention, dict_interventions)) %>%
  mutate(Epidemic_duration = ifelse(is.na(Epidemic_duration), 0, Epidemic_duration)) %>%
  ggplot(., aes(x = Scheme, y = Epidemic_duration, fill = intervention)) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.5, outliers = F) +
  facet_grid(rows = vars(network)) +
  scale_fill_manual(values = intervention_pal) +
  expand_limits(ymin = 0) +
  theme_classic() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 30, hjust = 1), 
        legend.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
        panel.border = element_rect(linewidth = 1),
        axis.line = element_blank()) +
  guides(fill=guide_legend(title.position="top", nrow=2)) + 
  labs(fill = "Intervention", y = "Epidemic duration (days)") 

# Time to peak
p2 = stats_df %>%
  mutate(intervention = factor(intervention, dict_interventions)) %>%
  filter(!is.na(Time_to_peak)) %>%
  ggplot(., aes(x = Scheme, y = Time_to_peak, col = intervention)) +
  geom_violin(position = position_dodge(width = 0.7), width = 0.7) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.05, outliers = F) +
  facet_grid(rows = vars(network)) +
  scale_color_manual(values = intervention_pal) +
  expand_limits(ymin = 0) +
  theme_classic() +
  theme(axis.title.x = element_blank(), panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
        panel.border = element_rect(linewidth = 1),
        axis.line = element_blank()) +
  labs(col = "Intervention", y = "Time to the peak (days)") 

# Extinction probability
p3 = stats_df %>%
  mutate(intervention = factor(intervention, dict_interventions)) %>%
  group_by(model, threshold, network, Scheme, intervention) %>%
  summarise(p = sum(Global == 0) / n(), .groups = "drop") %>%
  ggplot(., aes(x = Scheme, y = p, fill = intervention)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.5) +
  facet_grid(rows = vars(network)) +
  scale_fill_manual(values = intervention_pal) +
  expand_limits(ymin = 0) +
  theme_classic() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom",
        panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
        panel.border = element_rect(linewidth = 1),
        axis.line = element_blank()
  ) +
  labs(fill = "Intervention", y = "Extinction probability") 

ggarrange(p1, p2, p3, nrow = 3, align = "hv", common.legend = T, legend = "bottom")
ggsave("fig/influenza/epidemic_dynamics.png", height = 10, width = 6)

## Plot cumulative incidence reduction compared to baseline scenario------------
# Generate databases with estimates and confidence intervals of relative 
# cumulative incidence 
median_none = stats_df %>%
  filter(intervention == "None") %>%
  mutate(I_hcw = I_medical + I_paramedical) %>%
  pivot_longer(c(I_global, I_patient, I_medical, I_paramedical, I_hcw), names_to = "Category", values_to = "I") %>%
  group_by(network, Scheme, Category) %>%
  summarise(ref = median(I), .groups = "drop") 

sim_none = stats_df %>%
  filter(intervention == "None") %>%
  select(network, Scheme, nSim, matches("I_[a-z]+$"), -c(I_environment, I_contact)) %>%
  mutate(I_hcw = I_medical + I_paramedical) %>%
  pivot_longer(c(I_global, I_patient, I_medical, I_paramedical, I_hcw), names_to = "Category", values_to = "I_none")

## !!! THE RELATIVE ESTIMATES ARE NOT VALID WHEN THE MEDIANE IN THE BASELINE SCENARIO IS NULL !!!!
median_interventions = stats_df %>%
  select(network, Scheme, intervention, nSim, matches("I_[a-z]+$"), -c(I_environment, I_contact)) %>%
  arrange(network, Scheme, intervention, nSim) %>%
  mutate(I_hcw = I_medical + I_paramedical) %>%
  pivot_longer(c(I_global, I_patient, I_medical, I_paramedical, I_hcw), names_to = "Category", values_to = "I") %>%
  group_by(network, Scheme, Category) %>%
  wilcox_test(I ~ intervention, ref.group = "None", paired = F,
              p.adjust.method = "BH", conf.level = 0.95, alternative = "two.sided",
              detailed = T) %>%
  left_join(., median_none, by = c("network", "Scheme", "Category"))

# Plot of the reduction of cumulative incidence for all individuals, patients, and HCWs
median_interventions %>%
  mutate(
    estimate = ifelse(ref == 0, NA, -estimate/ref),
    conf.high = ifelse(ref == 0, NA, -conf.high/ref),
    conf.low = ifelse(ref == 0, NA, -conf.low/ref),
    group2 = factor(group2, dict_interventions),
    Category = case_when(
      Category == "I_global" ~ "All",
      Category == "I_patient" ~ "Patient",
      Category == "I_medical" ~ "Medical",
      Category == "I_paramedical" ~ "Paramedical",
      Category == "I_hcw" ~ "HCW"
    )
  ) %>%
  ggplot(., aes(x = Scheme, y = estimate, ymin = conf.high, ymax = conf.low, fill = Category)) +
  geom_hline(yintercept = -1, linetype = "dashed", col = "grey70") +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.5) +
  geom_errorbar(position = position_dodge(width = 0.7), width = 0.2) +
  ggh4x::facet_grid2(cols = vars(network), rows = vars(group2)) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(labels = scales::percent, breaks = c(0, -0.25, -0.5, -0.75, -1)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
        panel.border = element_rect(linewidth = 1),
        axis.line = element_blank()
  ) +
  guides(fill=guide_legend(title.position="top", position = "bottom", nrow=1)) +
  labs(fill = "Individual category", y = "Relative difference of the cumulative incidence over 90 days")

# Plot of the reduction of cumulative incidence for patients only
p4_patients = median_interventions %>%
  mutate(
    estimate = ifelse(ref == 0, NA, -estimate/ref),
    conf.high = ifelse(ref == 0, NA, -conf.high/ref),
    conf.low = ifelse(ref == 0, NA, -conf.low/ref),
    group2 = factor(group2, dict_interventions),
    Category = case_when(
      Category == "I_global" ~ "All",
      Category == "I_patient" ~ "Patient",
      Category == "I_medical" ~ "Medical",
      Category == "I_paramedical" ~ "Paramedical",
      Category == "I_hcw" ~ "HCW"
    )
  ) %>%
  filter(Category %in% c("Patient", "HCW")) %>%
  ggplot(., aes(x = Scheme, y = estimate, ymin = conf.high, ymax = conf.low)) +
  geom_hline(yintercept = -1, linetype = "dashed", col = "grey70") +
  # geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.4) +
  geom_col_pattern(
    aes(fill = network, pattern=Category), 
    pattern_fill = "black",
    pattern_density = 0.3,
    pattern_spacing = 0.05,
    # pattern_key_scale_factor = 0.6,
    position = position_dodge(width = 0.7), 
    width = 0.6
  ) +
  geom_errorbar(aes(group = interaction(Category, network)), 
                position = position_dodge(width = 0.7), width = 0.2) +
  facet_grid(rows = vars(group2)) +
  scale_fill_manual(values = network_pal) +
  # scale_pattern_fill_manual(values = network_pal) +
  scale_pattern_manual(values = c(HCW = "stripe", Patient = "none")) +
  scale_y_continuous(breaks = c(0, -.25, -.5, -.75, -1), labels = scales::percent) +
  theme_classic() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
        panel.border = element_rect(linewidth = 1),
        axis.line = element_blank()) +
  guides(
    fill = guide_legend(title.position="top", position = "bottom", nrow=1, 
                        override.aes = list(pattern = c("none", "none"))),
    pattern = guide_legend(title.position="top", position = "bottom", nrow=1, 
                           override.aes = list(fill = c("#00000000", "#00000000"),
                                               colour=c("black", "black")
                           ))
  ) +
  labs(pattern_fill = "Network", fill = "Network", 
       y = "Relative difference of the cumulative incidence among patients compared to the no-intervention scenario")

# 
median_interventions %>% filter(Category %in% c("I_patient", "I_hcw"), ref == 0)

# Final figure for paper
ggarrange(
  p4_patients,
  ggarrange(
    p1  + theme(legend.position = "none"), 
    p3 + theme(legend.position = "none"),
    get_legend(p1),
    nrow = 3,
    labels = c("B", "C", ""), 
    heights = c(1, 1, 0.2)
  ),
  labels = c("A", ""),
  ncol = 2,
  widths = c(0.6, 1)
)
ggsave("fig/influenza/figure5_influenza.png", height = 12, width = 11)

## Plot of intervention ranking-------------------------------------------------
# Ranking that displays cumulative incidence reduction when considering 
# all individuals
median_interventions %>%
  filter(Category == "All") %>%
  nest(.by = c(network, Scheme, Category)) %>%
  mutate(out = map(data, rank_interventions)) %>%
  unnest(out) %>%
  ggplot(., aes(x = Scheme, y = group2)) +
  geom_tile(aes(fill = estimate)) +  
  geom_text(aes(label = rank), color = "white", fontface='bold') +
  scale_fill_viridis(labels = scales::percent) +
  facet_grid(cols = vars(network)) +
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
ggsave("fig/influenza/rank_interventions.png", height = 5, width = 8)
