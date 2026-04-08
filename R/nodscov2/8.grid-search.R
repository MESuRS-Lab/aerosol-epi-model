################################################################################
##
##                  Grid-search procedure
##
################################################################################

## Working environment----------------------------------------------------------
## DATA MANAGEMENT
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(gtsummary)
library(DescTools)
library(rstatix)
library(rcompanion)
library(performance)
library(igraph)
library(ggnetwork)

## LOAD FUNCTIONS
source('R/nodscov2/helper-functions-simulations.R')
source('R/nodscov2/dictionaries.R')

# Concatenation of the results when multiple nextflow rounds-------------------
# f_final = lapply(list.files("grid_search/unique_files", pattern = "summary.*.csv", full.names = T),
# function(f_name) { read.csv2(f_name, header = T)})
# f_final = do.call("bind_rows", f_final)
# f_final %>% count(beta_e, beta_c, network) %>% nrow()
# write.csv2(f_final, "grid_search/resu_simu_all.txt", row.names = F)

## Load data--------------------------------------------------------------------
# All conditions
stats = read.csv2("grid_search/resu_simu_all.txt", header = T) %>%
  mutate(
    network = ifelse(network == "icu1", "ICU1", "ICU2"),
    beta_e = factor(
      ifelse(beta_e == 0, "0", paste0("1/", gsub(" ", "", format(round(1/beta_e), scientific = F)))), 
      c('0', '1/200', '1/150', '1/100', '1/90', '1/80', '1/70', '1/60', '1/50', '1/45', '1/40', '1/30', '1/20', '1/10')))
head(stats)
nrow(stats)

# Dimensions
stats %>% count(beta_e, beta_c) %>% filter(n != 200)

# Sub-selection corresponding to studied scenarios
stats_sub = stats %>%
  mutate(
    Pathway = case_when(
      beta_e == "1/45" & beta_c == 0.75 ~ "Pathway 1",
      beta_e == "1/60" & beta_c == 1 ~ "Pathway 2",
      beta_e == "1/70" & beta_c == 1.25 ~ "Pathway 3",
      beta_e == "1/100" & beta_c == 1.5 ~ "Pathway 4",
      beta_e == "1/150" & beta_c == 1.75 ~ "Pathway 5",
      .default = NA
    )
  ) %>%
  filter(!is.na(Pathway))

stats_sub_principal = stats_sub %>%
  filter(model == "linear", threshold == 60)

## Plots------------------------------------------------------------------------
# Summary statistics
stats_sum = stats %>%
  filter(model == "linear", threshold == 60) %>%
  group_by(beta_e, beta_c, network, model, threshold) %>%
  summarise(Mean = mean(Global), Median = median(Global), .groups = "drop") %>%
  pivot_longer(c(Mean, Median), values_to = "val", names_to = "metric")

stats_sum_sel = stats_sum %>%
  filter(
    beta_e == "1/45" & beta_c == 0.75 | 
    beta_e == "1/60" & beta_c == 1 | 
    beta_e == "1/70" & beta_c == 1.25  | 
    beta_e=="1/100" & beta_c==1.5 | 
    beta_e == "1/150" & beta_c == 1.75
    )

ggplot(stats_sum, aes(x = beta_c, col = beta_e, y = val, group = beta_e)) +
  geom_point() +
  geom_line(linewidth = 0.3) +
  geom_point(data = stats_sum_sel, col = "black", shape = 4, size = 4) +
  facet_grid(rows = vars(metric), cols = vars(network)) +
  viridis::scale_colour_viridis(option = "magma", discrete = T, direction = -1) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
    panel.border = element_rect(linewidth = 1),
    axis.line = element_blank()
    ) +
  labs(y = "Global SAR", x = expression(beta[c]), col = expression(beta[e])) +
  expand_limits(y = c(0,1))
ggsave("fig/grid_search/grid-search.png", width = 8, height = 5)

# Stratified SAR 
  # by individual category
pa = stats_sub_principal %>%
  pivot_longer(c(Patient, Paramedical, Medical), names_to = "Category", values_to = "SAR") %>%
  ggplot(., aes(x = Pathway, y = SAR, group = interaction(Pathway, Category))) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.4, outliers = F) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.1, jitter.height = 0), aes(col = Category), alpha = 0.5) +
  scale_color_manual(values = pal) +
  facet_grid(cols = vars(network)) +
  labs() +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))

  # by transmission route
pb = stats_sub_principal %>%
  pivot_longer(c(Contact, Environment), names_to = "Route", values_to = "SAR") %>%
  mutate(Route = ifelse(Route == "Contact", "Short-range", "Long-range")) %>%
  ggplot(., aes(x = Pathway, y = SAR, group = interaction(Pathway, Route))) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.4, outliers = ) +
  geom_jitter(aes(col = Route), position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.1, jitter.height = 0), alpha = 0.5) +
  facet_grid(cols = vars(network)) +
  scale_color_manual(values = env_pal) +
  labs() +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))

p_all = ggarrange(pa, pb, ncol = 1, nrow = 2, align = "hv")
p_all
ggsave("fig/grid_search/sar_stratified.png", p_all, height = 6, width = 8)

# SAR by transmission route and individual category
p_jonc = stats_sub_principal %>%
  select(network, Pathway, matches("[A-Z][a-z]+_Contact|[A-Z][a-z]+_Environment")) %>%
  pivot_longer(matches(".*_Contact|.*_Environment"), names_to = c("Category", "Route"), names_pattern = "(.*)_(.*)", values_to = "SAR") %>%
  arrange(network, Category, Route, Pathway) %>%
  nest(.by = c("network", "Category", "Route")) %>%
  mutate(data = map(data, jonckheere_test)) %>%
  unnest(data) %>%
  filter(p <= 0.05) %>%
  mutate(p = paste0("trend p-value=", base::signif(p, digits = 1)),
         Route = ifelse(Route == "Contact", "Close proximity", "Long-distance"), 
         SAR = 1, 
         Pathway = NA)

stats_sub_principal %>%
  pivot_longer(matches("Patient_|Paramedical_|Medical_"), names_to = c("Category", "Route"), names_pattern = "(.*)_(.*)", values_to = "SAR") %>%
  mutate(Route = ifelse(Route == "Contact", "Close proximity", "Long-distance")) %>%
  ggplot(., aes(x = Route, y = SAR, group = interaction(Category, Route, Pathway))) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.4, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05, jitter.height = 0), 
              alpha = 0.5, aes(col = Category)) +
  geom_text(data = p_jonc, aes(label = p)) +
  facet_grid(cols = vars(network), rows = vars(Category)) +
  scale_color_manual(values = pal) +
  expand_limits(y = c(0,1)) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank())
ggsave("fig/grid_search/sar_stratified_trends.png", height = 8, width = 8)

## Show scenario selection results----------------------------------------------
# Network ICU1
admission_nodscov2 = read.csv2("data/data-nodscov2/clean/admission_cleaned_icu1.csv")
graph_data = read.csv2("data/data-synthetic-graphs/input/interactions_icu1.csv") %>%
  mutate(date_posix = floor_date(as_date(date_posix), "day")) %>%
  group_by(from, to, date_posix) %>%
  summarise(length = sum(length), .groups = "drop") %>%
  arrange(date_posix) %>%
  filter(date_posix == min(date_posix))
graph_example = simplify(graph_from_data_frame(graph_data, directed = F))

vertex_atts = data.frame(id = vertex_attr(graph_example, "name")) %>%
  left_join(admission_nodscov2, "id")

graph_example = graph_example %>%
  set_vertex_attr("cat", value = vertex_atts$cat)

p1 = ggplot(ggnetwork(graph_example, layout = igraph::layout_with_kk(graph_example)),
            aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(alpha = 0.5, linewidth = 0.5) +
  geom_nodes(aes(colour = cat), size = 2) +
  theme_blank() +
  scale_colour_discrete(type = pal) +
  labs(colour = "Category")

# Network ICU2
admission_nodscov2 = read.csv2("data/data-nodscov2/clean/admission_cleaned_icu2.csv")
graph_data = read.csv2("data/data-synthetic-graphs/input/interactions_icu2.csv") %>%
  mutate(date_posix = floor_date(as_date(date_posix), "day")) %>%
  group_by(from, to, date_posix) %>%
  summarise(length = sum(length), .groups = "drop") %>%
  arrange(date_posix) %>%
  filter(date_posix == min(date_posix))
graph_example = simplify(graph_from_data_frame(graph_data, directed = F))

vertex_atts = data.frame(id = vertex_attr(graph_example, "name")) %>%
  left_join(admission_nodscov2, "id")

graph_example = graph_example %>%
  set_vertex_attr("cat", value = vertex_atts$cat)

p2 = ggplot(ggnetwork(graph_example, layout = igraph::layout_with_kk(graph_example)),
            aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(alpha = 0.5, linewidth = 0.5) +
  geom_nodes(aes(colour = cat), size = 2) +
  theme_blank() +
  scale_colour_discrete(type = pal) +
  labs(colour = "Category")

# Global SAR
p_sar = stats_sub_principal %>%
  group_by(network) %>%
  wilcox_test(Global ~ Pathway, p.adjust.method = "BH") %>%
  add_xy_position(x = "Pathway") %>%
  filter(p.adj <= 0.05)
p3 = stats_sub_principal %>%
  arrange(network, Pathway) %>%
  ggboxplot(., x = "Pathway", y = "Global", width = 0.4, outlier.shape = NA) +
  stat_pvalue_manual(p_sar, label = "p.adj.signif") +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, aes(col = network)) +
  facet_grid(cols = vars(network)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(y = "SAR (in %)", col = "Network") +
  scale_color_manual(values = network_pal) +
  expand_limits(y = c(0,1))

# SAR by source
p4 = stats_sub_principal %>%
  pivot_longer(c(Contact, Environment), names_to = "Route", values_to = "SAR") %>%
  mutate(Route = ifelse(Route == "Contact", "Short-range", "Long-range")) %>%
  ggplot(., aes(x = Pathway, y = SAR, group = interaction(Pathway, Route))) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.4, outliers = ) +
  geom_jitter(aes(col = Route), position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.1, jitter.height = 0), 
              alpha = 0.5) +
  facet_grid(cols = vars(network)) +
  scale_color_manual(values = env_pal) +
  scale_y_continuous(labels = scales::percent) +
  labs() +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))

# Save figure
p_all = ggarrange(
  ggarrange(p1, p2, ncol = 2, labels = c("A", ""), common.legend = T, legend = "right"), 
  ggarrange(p3, p4, ncol = 1, nrow = 2, legend = "right", align = "hv", labels = c("B", "C")),
  ncol = 1, nrow = 2, heights = c(1/3, 2/3))
p_all
ggsave("fig/grid_search/figure3_on_grid_search_simulations.png", p_all, height = 7, width = 7)

## Global epidemic metrics------------------------------------------------------
# Extinction probability
p_proba = stats_sub_principal %>%
  mutate(Extinction = ifelse(Global == 0, "Yes", "No")) %>%
  count(Pathway, network, Extinction) %>%
  nest(.by = network) %>%
  mutate(data = map(data, chisq_test_df)) %>%
  unnest(data) %>%
  mutate(y = y+0.05) %>%
  filter(p.adj <= 0.05)
stats_ext = stats_sub_principal %>% 
  group_by(network, Pathway) %>% 
  summarise(prop = sum(Global == 0) / n(), .groups = "drop") %>%
  mutate(Category = rep(names(pal)[1:3], n())[1:n()])
p1 = ggplot(stats_ext) +
  geom_bar(stat = "identity", aes(x = Pathway, y = prop, fill = network), width = 0.4) +
  geom_point(x = 0, y = 0, aes(col = Category)) +
  stat_pvalue_manual(p_proba, label = "p.adj.signif", y.position = "y") +
  facet_grid(cols = vars(network)) +
  expand_limits(y = c(0,1)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 20, hjust = 1), 
        legend.position = "bottom") +
  scale_fill_manual(values = network_pal) +
  scale_color_manual(values = pal[1:3]) +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5),
         col = guide_legend(title.position="top", title.hjust = 0.5)) +
  labs(y = "Extinction probability", fill = "Network", col = "Category")

# Epidemic duration
p_duration = stats_sub_principal %>%
  group_by(network) %>%
  wilcox_test(Epidemic_duration ~ Pathway, p.adjust.method = "BH") %>%
  add_xy_position(x = "Pathway") %>%
  filter(p.adj <= 0.05)
p2 = stats_sub_principal %>%
  filter(!is.na(Epidemic_duration)) %>%
  arrange(network, Pathway) %>%
  ggboxplot(., x = "Pathway", y = "Epidemic_duration", width = 0.4, outlier.shape = NA) +
  stat_pvalue_manual(p_duration, label = "p.adj.signif") +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, aes(col = network)) +  
  facet_grid(cols = vars(network)) +
  labs(y = "Epidemic duration (in days)") +
  theme_bw() +
  scale_color_manual(values = network_pal) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1)) +
  expand_limits(ymin = 0)

# Time to the peak with statistical tests
p_peak = stats_sub_principal %>%
  group_by(network) %>%
  wilcox_test(Time_to_peak ~ Pathway, p.adjust.method = "BH") %>%
  add_xy_position(x = "Pathway") %>%
  filter(p.adj <= 0.05)
p3 = stats_sub_principal %>%
  arrange(network, Pathway) %>%
  filter(!is.na(Time_to_peak)) %>%
  ggboxplot(., x = "Pathway", y = "Time_to_peak", width = 0.4, outlier.shape = NA) +
  stat_pvalue_manual(p_peak, label = "p.adj.signif") +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, aes(col = network)) +
  facet_grid(cols = vars(network)) +
  labs(y = "Time to the peak (in days)") +
  theme_bw() +
  scale_color_manual(values = c(network_pal, env_pal), drop = F) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1)) +
  expand_limits(ymin = 0)

# SAR difference between the two routes for the different individual categories
p4 = stats_sub_principal %>%
  pivot_longer(matches("[A-Z][a-z]+_Environment|[A-Z][a-z]+_Contact"), names_to = c("Category", "Route"), names_pattern = "(.*)_(.*)", values_to = "SAR") %>%
  group_by(Pathway, network, Category) %>%
  wilcox_test(SAR ~ Route, paired = T, ref.group = "Contact", alternative = "two.sided", detailed = T) %>%
  ggplot(., aes(x = Pathway, y = estimate, ymin = conf.low, ymax = conf.high, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.4) +
  geom_errorbar(position = "dodge", width = 0.4) +
  facet_grid(cols = vars(network)) +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(y = "Contact - Environment SAR\n(Wilcoxon test estimate)")

# Save figure 
p_all = ggarrange(p1, p2, p3, p4, 
                  ncol = 2, nrow = 2, align = "hv", 
                  common.legend = T, legend = "bottom",
                  labels = c("A", "B", "C", "D"))
p_all
ggsave("fig/grid_search/figure4_on_grid_search_data.png", p_all, height = 6, width = 10)

## Statistical comparisons between scenarios------------------------------------
# Comparison of basic epidemic metrics
tab_global = stats_sub_principal %>%
  mutate(Extinction = ifelse(Global == 0, 1, 0)) %>%
  distinct(network, nSim, Pathway, Extinction, Global, Epidemic_duration, Time_to_peak) %>%
  pivot_wider(names_from = network, values_from = c(Extinction, Global, Epidemic_duration, Time_to_peak)) %>%
  select(Pathway, matches("_ICU1"), matches("_ICU2")) %>%
  tbl_summary(
    by = Pathway, 
    missing = "no",
    list(
      Global_ICU1 ~ "Global SAR (in %)",
      Extinction_ICU1 ~ "Extinction probability (in %)",
      Time_to_peak_ICU1 ~ "Time to the peak (in days)",
      Epidemic_duration_ICU1 ~ "Epidemic duration (in days)",
      Global_ICU2 ~ "Global SAR (in %)",
      Extinction_ICU2 ~ "Extinction proportion",
      Time_to_peak_ICU2 ~ "Time to the peak (in days)",
      Epidemic_duration_ICU2 ~ "Epidemic duration (in days)"
      )
    ) %>%
  add_variable_group_header(header = "ICU1", variables = matches("_ICU1")) %>%
  add_variable_group_header(header = "ICU2", variables = matches("_ICU2")) %>%
  add_p() %>%
  bold_levels() %>%
  bold_p() %>%
  modify_abbreviation("ICU = Intensive Care Unit; SAR = Secondary Attack Rate")

as_gt(tab_global)
gt::gtsave(as_gt(tab_global), "fig/grid_search/tab_global_metrics.png")

# Participant categories are not affected by the same route 
for (n in c("ICU1", "ICU2")) {
  out = stats_sub_principal %>%
    filter(network == n) %>%
    select(matches("[A-Z].+_Contact|[A-Z].+_Environment"), Pathway, nSim) %>%
    distinct() %>%
    pivot_longer(dplyr::matches("_Contact|_Environment"), names_to = c("Category", "Route"), names_pattern = "(.*)_(.*)", values_to = "SAR") %>%
    arrange(Pathway, nSim, Category) %>%
    mutate(Pathway_Category = paste0(gsub(" ", "", Pathway), "_", Category)) %>%
    dplyr::select(-c(Pathway, Category)) %>%
    pivot_wider(names_from = Pathway_Category, values_from = SAR) %>%
    select(nSim, Route, ends_with("_Medical"), ends_with("_Paramedical"), ends_with("_Patient")) %>%
    tbl_summary(
      by = Route, 
      include = -nSim,
      label = list(Pathway1_Patient = "Pathway 1",
                   Pathway1_Medical = "Pathway 1",
                   Pathway1_Paramedical = "Pathway 1",
                   Pathway2_Patient = "Pathway 2",
                   Pathway2_Medical = "Pathway 2",
                   Pathway2_Paramedical = "Pathway 2",
                   Pathway3_Patient = "Pathway 3",
                   Pathway3_Medical = "Pathway 3",
                   Pathway3_Paramedical = "Pathway 3",
                   Pathway4_Patient = "Pathway 4",
                   Pathway4_Medical = "Pathway 4",
                   Pathway4_Paramedical = "Pathway 4",
                   Pathway5_Patient = "Pathway 5",
                   Pathway5_Medical = "Pathway 5",
                   Pathway5_Paramedical = "Pathway 5")
    ) %>%
    add_variable_group_header(header = "Medical staff", variables = ends_with("Medical")) %>%
    add_variable_group_header(header = "Paramedical staff", variables = ends_with("Paramedical")) %>%
    add_variable_group_header(header = "Patients", variables = ends_with("Patient")) %>%
    modify_header(label = paste0("**", n, "**")) %>%
    bold_levels() %>%
    add_p() %>%
    bold_p()
  gt::gtsave(as_gt(out), paste0("fig/grid_search/tab_stratified_", tolower(n), ".png"))
  print(as_gt(out))
}

# # Sensitivity analyses on the dose-response function----------------------------
# # Dose-response model (linear, log-linear, exponential) and time spent out of the 
# # ward (60, 120, 180 minutes)
# stats_sub %>%
#   filter(model == "linear") %>%
#   mutate(threshold = factor(threshold, c("60", "120", "180"))) %>%
#   ggplot(., aes(x=Pathway, y=Global, group = interaction(Pathway, model, threshold))) +
#   geom_boxplot(position = position_dodge(width = 0.7), width = 0.4, outliers = F) +
#   geom_jitter(aes(col = threshold), alpha = 0.2, position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0.7)) +
#   facet_grid(rows = vars(network), cols = vars(model)) +
#   scale_color_manual(values = threshold_pal) + 
#   theme_bw() +
#   expand_limits(y = c(0,1)) +
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom") +
#   labs(y = "SAR", col = "Threshold before leaving the ward (in min)")
# ggsave("fig/grid_search/sensitivity.png", height = 5, width = 6)
# 
# # Does the threshold changes the results when stratified by individual 
# # category
# stats_sub %>%
#   filter(model == "linear") %>%
#   pivot_longer(matches("Patient_|Paramedical_|Medical_"), names_to = c("Category", "Route"), names_pattern = "(.*)_(.*)", values_to = "SAR") %>%
#   mutate(Route = ifelse(Route == "Contact", "Close proximity", "Long-distance"),
#          threshold = factor(threshold, c("60", "120", "180"))) %>%
#   ggplot(., aes(x = Route, y = SAR, group = interaction(Category, threshold, Route, Pathway))) +
#   geom_boxplot(position = position_dodge(width = 0.7), width = 0.4, outlier.shape = NA) +
#   geom_jitter(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.05, jitter.height = 0), 
#               alpha = 0.5, aes(col = threshold)) +
#   facet_grid(cols = vars(network), rows = vars(Category)) +
#   scale_color_manual(values = threshold_pal) +
#   expand_limits(y = c(0,1)) +
#   theme_bw() +
#   theme(legend.position = "bottom", axis.title.x = element_blank()) +
#   labs(col = "Threshold before leaving the ward (in min)")
# 
