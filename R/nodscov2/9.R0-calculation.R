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

## Represent correspondance-----------------------------------------------------
# Load input data 
r0 = read.csv2("R0_calculation/r0_all.txt", header = T)

# Global R0
r0_global = r0 %>%
  group_by(network, Scheme, R0) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(index = 'Global')

# Plot
r0 %>%
  group_by(network, Scheme, index, R0) %>%
  summarise(n = n(), .groups = "drop") %>%
  bind_rows(r0_global) %>%
  ggplot(., aes(x = R0, y = n, fill = index)) + 
  geom_bar(stat ="identity", position = position_dodge(width = 0.8)) +
  facet_grid(cols = vars(network), rows = vars(Scheme)) +
  theme_bw() + 
  theme(strip.background = element_rect(fill = "#00000000")) +
  expand_limits(y = 0) +
  labs(y = "Frequency", col = "Index case", x = "R0 (median, 2.5 and 97.5 percentiles)")

# Global R0 - median, 2.5 and 97.5 percentiles
r0_global_sumstat = r0 %>%
  group_by(network, Scheme) %>%
  summarise(
    R0_mid = round(mean(R0), 2),
    R0_low = quantile(R0, 0.025),
    R0_up = quantile(R0, 0.975),
    .groups = "drop"
  ) %>%
  mutate(index = 'Global')

# Median, 2.5 and 97.7 percentiles
r0_sumstat = r0 %>%
  group_by(network, Scheme, index) %>%
  summarise(
    R0_mid = round(mean(R0), 2),
    R0_low = quantile(R0, 0.025),
    R0_up = quantile(R0, 0.975),
    .groups = "drop"
  ) %>%
  bind_rows(., r0_global_sumstat)

# Table
tab = r0_sumstat %>%
  mutate(
    R0_low = round(R0_low, 2),
    R0_up = round(R0_up, 2),
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
tab
gtsave(tab, "fig/baseline_scenario/R0_correspondance_tab.png")
gtsave(tab, "fig/baseline_scenario/R0_correspondance_tab.docx")

