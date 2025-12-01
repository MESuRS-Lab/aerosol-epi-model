################################################################################
##
##            Figures of epidemic dynamics in baseline scenario
##
################################################################################

## Working environment----------------------------------------------------------
## R Packages
rm(list=ls())
library(tidyverse)
library(hms)
library(ggpubr)
library(gtsummary)
library(DescTools)
library(Rcpp)
library(rstatix)
library(rcompanion)
library(performance)
library(igraph)
library(ggnetwork)

## Load data--------------------------------------------------------------------
# Number of individuals
ind_numbers = data.frame()
load("out/parameters-synthetic-icu2-60.rda")
ind_numbers = bind_rows(ind_numbers, 
                        global_status %>% 
                          summarise(n_all = n(), n_patient = sum(grepl("^PA", id)), n_hcw = sum(grepl("^PE", id))) %>%
                          mutate(network = "ICU2"))
rm(list=setdiff(ls(), "ind_numbers"))

load("out/parameters-synthetic-icu1-60.rda")
ind_numbers = bind_rows(ind_numbers, 
                        global_status %>% 
                          summarise(n_all = n(), n_patient = sum(grepl("^PA", id)), n_hcw = sum(grepl("^PE", id))) %>%
                          mutate(network = "ICU1"))
dict_rooms_simplified = global_environment[[1]] %>% distinct(room, id_room) %>% mutate(room = ifelse(room == as.character(id_room), paste0("Patient Room ", id_room), room))
dict_rooms_simplified = SetNames(dict_rooms_simplified$room, dict_rooms_simplified$id_room)
rm(list=setdiff(ls(), c("ind_numbers", "admission", "begin_date", "end_date", "dict_rooms_simplified")))
   
# Movement trajectories 
load("data/data-synthetic-graphs/loc/icu1-simulated-reconstructed-locations-60.rda")
rm(list=setdiff(ls(), c("ind_numbers", "paths", "admission", "begin_date", "end_date", "dict_rooms_simplified")))

# Load functions and dictionaries
source('R/nodscov2/helper-functions-simulations.R')
source("R/nodscov2/helper-functions.R")
source('R/nodscov2/dictionaries.R')

# Epidemics for all conditions
stats_df = read.csv2("nextflow_interventions/resu_interventions_all.txt", header = T) %>%
  filter(model == "linear", threshold == 60, intervention == "None") %>%
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
    network = toupper(network)
  ) 
head(stats_df)

## Plot individual movement trajectories and aerosol generation-----------------
# Movement trajectories 
set.seed(20251013)

agenda = read.csv2("data/data-synthetic-graphs/input/agenda_icu1.csv")
individuals_to_show = c(
  sample(admission$id[admission$cat == "Paramedical"], 1),
  sample(admission$id[admission$cat == "Medical"], 1),
  sample(admission$id[admission$cat == "Patient"], 1)
)

individual_paths = data.frame()
for (y in individuals_to_show) {
  if (grepl("^PA-", y)) {
    days_of_presence = seq(admission %>% filter(id == y) %>% pull(firstDate),
                           admission %>% filter(id == y) %>% pull(lastDate),
                           1)
  } else {
    days_of_presence = agenda %>%
      filter(id == y) %>%
      mutate(k = 1:n(), firstDate = as_date(floor_date(as_datetime(firstDate), "day")), lastDate = as_date(floor_date(as_datetime(lastDate), "day"))) %>%
      nest(.by = k) %>%
      mutate(data = map(data, unroll_days)) %>%
      unnest(data) %>%
      pull(data)
  }
  
  ind_cat = as.character(admission$cat[admission$id == y])
  
  random_day = sample(unique(days_of_presence), 1)
  individual_path = data.frame(
    ti = seq(begin_date, end_date-30, 30),
    loc = as.vector(paths[, y]),
    cat = ind_cat,
    id = paste0(ind_cat, ifelse(ind_cat != "Patient", " staff", ""), " n°", gsub("PE-|PA-", "", y))
  ) %>%
    filter(floor_date(ti, "day") == random_day) %>%
    mutate(
      loc = recode(loc, !!!c(dict_rooms_simplified, "-1" = "Absent")),
      ti = as_hms(ti)
      ) %>%
    filter(!is.na(ti)) # Remove midnight
  individual_paths = bind_rows(individual_paths, individual_path)
}

p1 = individual_paths %>%
  mutate(loc = factor(loc, c("Medical Staff Room", "Paramedical Staff Room", 
                             "Office", "Nursing station", "Corridor", 
                             paste0("Patient Room ", c(2,3,4,5,6,8,10,13,15,16,17)), 
                             "Absent"))) %>%
    mutate(la = lag(loc), le = lead(loc)) %>%
    ungroup() %>%
    mutate(xstart = case_when(is.na(la) | loc != la | (!is.na(la) & cat != lag(cat)) ~ ti, .default = NA), 
           xend = case_when(is.na(le) | loc != le | (!is.na(le) & cat != lead(cat)) ~ ti, .default = NA)) %>%
    filter(!is.na(xstart) | !is.na(xend)) %>%
    mutate(
      xend_le = lead(xend),
      xend = case_when(
        cat == lead(cat) & loc == le ~ xend_le,
        .default = xend
    )) %>%
    filter(!is.na(xstart)) %>%
  ggplot(., aes(x = xstart, xend = xend, y = loc, col = cat)) +
    geom_point(aes(x=xend, y=loc), size = 0.5, shape = 20) +
    geom_point(aes(x=xstart, y=loc), size = 0.5, shape = 20) +
    geom_segment(linewidth = 2) +
  facet_wrap(facets = vars(id), ncol = 1, scales = "free_x") + 
  scale_color_manual(values = pal) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = 11/22, color = "grey90"),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(linewidth = 1),
    axis.line = element_blank()
    ) +
  labs(x = "Time of the day")

# Generate one simulation
sourceCpp("cpp/dev-model-nodscov2.cpp")
load("out/parameters-synthetic-icu1-60.rda")
beta_c <- 1.25
beta_e <- 1/70
intervention = "None"

id_index <- sample(x = admission$id[admission$firstDate == as_date(floor_date(begin_date," day")) & 
                                      admission$firstDate != admission$lastDate], size = 1)
global_status <- global_status %>%
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

result_noint <- simulation(global_interaction = global_interaction,
                     global_environment = global_environment,
                     global_data = global_data,
                     global_status = global_status,
                     beta_c = beta_c,
                     beta_e = beta_e,
                     B = B,
                     nu = nu,
                     mu = mu,
                     env_model = env_model,
                     dt = dt,
                     intervention = intervention
)

# Within-room aerosol dynamics 
inf_status_noint = result_noint$global_status %>% 
  filter(t_inf>=0) %>% 
  select(id, t_infectious_start, t_recover) %>% 
  mutate(t_infectious_start = as.POSIXct("2020-01-01 00:00:30") + 30 * t_infectious_start,
         t_recover = as.POSIXct("2020-01-01 00:00:30") + 30 * t_recover) %>%
  nest(.by = id) %>% 
  mutate(date_hour = map(data, function(df) seq(df$t_infectious_start, df$t_recover, by = 30))) %>% 
  unnest(date_hour) %>% 
  select(id, date_hour)

one_day_loc_noint = do.call("rbind", mapply(function(x, y) {x$ti = y; return(x)}, x=result_noint$global_data[1:(24*60*2)], y=1:(24*60*2), SIMPLIFY = F)) %>%
  mutate(date_hour = as.POSIXct("2020-01-01 00:00:30") + 30 * ti) %>%
  filter(location_ti >= 50 | location_ti == 5) %>%
  mutate(room = case_when(
    location_ti == 5 ~ "Patient Room 5",
    location_ti == 54 ~ "Corridor",
    location_ti == 50 ~ "Medical Staff Room",
    location_ti == 52 ~ "Nursing station",
    location_ti == 53 ~ "Office",
    location_ti == 51 ~ "Paramedical Staff Room"
  )) %>% 
  inner_join(., inf_status_noint, by = c("date_hour", "id")) %>%
  group_by(room, date_hour) %>%
  summarise(n_inf = n(), .groups = "drop")

one_day_all_noint = do.call("rbind", mapply(function(x, y) {x$ti = y; return(x)}, x=result_noint$global_environment[1:(24*60*2)], y=1:(24*60*2), SIMPLIFY = F)) %>%
  mutate(date_hour = as.POSIXct("2020-01-01 00:00:30") + 30 * ti) %>%
  filter(id_room >= 50 | id_room==5) %>%
  mutate(room = ifelse(room==5, paste0("Patient Room ", room), room)) %>%
  left_join(., one_day_loc_noint, by = c("room", "date_hour")) %>%
  mutate(n_inf = ifelse(is.na(n_inf), 0, n_inf))

p2 = ggplot(one_day_all_noint, aes(x = date_hour, y = env/volume)) +
  geom_line(aes(y = n_inf*1000), col = "red") +
  geom_line() +
  scale_y_continuous(
    "Aerosol concentration (per m3)",
    sec.axis = sec_axis(~ . / 1000)
  ) +
  facet_wrap(facets = vars(room), ncol = 1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank(),
    axis.line.y.right = element_line(color = "red"), 
    axis.ticks.y.right = element_line(color = "red"),
    axis.text.y.right = element_text(color = "red"),
    panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
    panel.border = element_rect(linewidth = 1),
    axis.line = element_blank()
    )

# With intervention
intervention = "Improved ventilation patients"
result_int <- simulation(global_interaction = global_interaction,
                         global_environment = global_environment,
                         global_data = global_data,
                         global_status = global_status,
                         beta_c = beta_c,
                         beta_e = beta_e,
                         B = B,
                         nu = nu,
                         mu = mu,
                         env_model = env_model,
                         dt = dt,
                         intervention = intervention,
                         mu_int = (log(100)/(1/8)) * 24
)

inf_status = result_int$global_status %>% 
  filter(t_inf>=0) %>% 
  select(id, t_infectious_start, t_recover) %>% 
  mutate(t_infectious_start = as.POSIXct("2020-01-01 00:00:30") + 30 * t_infectious_start,
         t_recover = as.POSIXct("2020-01-01 00:00:30") + 30 * t_recover) %>%
  nest(.by = id) %>% 
  mutate(date_hour = map(data, function(df) seq(df$t_infectious_start, df$t_recover, by = 30))) %>% 
  unnest(date_hour) %>% 
  select(id, date_hour)

one_day_loc = do.call("rbind", mapply(function(x, y) {x$ti = y; return(x)}, x=result_int$global_data[1:(24*60*2)], y=1:(24*60*2), SIMPLIFY = F)) %>%
  mutate(date_hour = as.POSIXct("2020-01-01 00:00:30") + 30 * ti) %>%
  filter(location_ti >= 50 | location_ti == 5) %>%
  mutate(room = case_when(
    location_ti == 5 ~ "Patient Room 5",
    location_ti == 54 ~ "Corridor",
    location_ti == 50 ~ "Medical Staff Room",
    location_ti == 52 ~ "Nursing station",
    location_ti == 53 ~ "Office",
    location_ti == 51 ~ "Paramedical Staff Room"
  )) %>% 
  inner_join(., inf_status, by = c("date_hour", "id")) %>%
  group_by(room, date_hour) %>%
  summarise(n_inf = n(), .groups = "drop")

one_day_all = do.call("rbind", mapply(function(x, y) {x$ti = y; return(x)}, x=result_int$global_environment[1:(24*60*2)], y=1:(24*60*2), SIMPLIFY = F)) %>%
  mutate(date_hour = as.POSIXct("2020-01-01 00:00:30") + 30 * ti) %>%
  filter(id_room >= 50 | id_room==5) %>%
  mutate(room = ifelse(room==5, paste0("Patient Room ", room), room)) %>%
  left_join(., one_day_loc, by = c("room", "date_hour")) %>%
  mutate(n_inf = ifelse(is.na(n_inf), 0, n_inf))

p3 = ggplot(one_day_all, aes(x = date_hour, y = env/volume)) +
  geom_line(aes(y = n_inf*1000), col = "red") +
  geom_line() +
  scale_y_continuous(
    "",
    sec.axis = sec_axis(~ . / 1000, name = "Number of infectious individuals")
  ) +
  expand_limits(y = c(0, 3800)) +
  facet_wrap(facets = vars(room), ncol = 1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank(),
    axis.line.y.right = element_line(color = "red"), 
    axis.ticks.y.right = element_line(color = "red"),
    axis.text.y.right = element_text(color = "red"),
    panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
    panel.border = element_rect(linewidth = 1),
    axis.line = element_blank()
    )

# Combined plot
p_final = ggarrange(p1, p2, p3, ncol = 3, labels = c("A", "B", "C"), 
                    widths = c(1, 0.5, 0.5))
p_final
ggsave("fig/paper/figure2.png", p_final, height = 8, width = 12)

## Plot SAR for baseline scenario without interventions-------------------------
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
  geom_edges(alpha = 0.2, linewidth = 0.2) +
  geom_nodes(aes(fill = cat), color = "black", size = 2, shape = 21) +
  theme_blank() +
  scale_fill_discrete(type = pal) +
  labs(fill = "Category")

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
  geom_edges(alpha = 0.2, linewidth = 0.2) +
  geom_nodes(aes(fill = cat), color = "black", size = 2, shape = 21) +
  theme_blank() +
  scale_fill_discrete(type = pal) +
  labs(fill = "Category")

# Global SAR
# p_sar = stats_df %>%
#   group_by(network) %>%
#   wilcox_test(Global ~ Pathway, p.adjust.method = "BH") %>%
#   add_xy_position(x = "Pathway") %>%
#   filter(p.adj <= 0.05)

p3 = stats_df %>%
  arrange(network, Pathway) %>%
  # ggboxplot(., x = "Pathway", y = "Global", width = 0.4, outlier.shape = NA, col = "network") +
  # stat_pvalue_manual(p_sar, label = "p.adj.signif") +
  ggplot(., aes(x = Pathway, y = Global, col = network)) +
  geom_boxplot(width = 0.4, outliers = F) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.1, size = 0.1, aes(col = network)) +
  facet_grid(cols = vars(network)) +
  scale_color_manual(values = network_pal) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Global SAR", col = "Network") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
    panel.border = element_rect(linewidth = 1),
    axis.line = element_blank()
    ) +
  expand_limits(y = c(0,1))

stats_df %>%
  group_by(network) %>%
  summarise(
    m = round(mean(Global)*100),
    me = round(median(Global)*100),
    iqr = round(IQR(Global)*100),
    .groups = "drop"
  )

stats_df %>%
  nest(.by = network) %>%
  mutate(kruskall_wallis_p = map_dbl(data, ~kruskal.test(.$Global ~ .$Pathway)$p.value))

# SAR by source
p4 = stats_df %>%
  pivot_longer(c(Contact, Environment), names_to = "Route", values_to = "SAR") %>%
  mutate(Route = ifelse(Route == "Contact", "Short-range", "Long-range")) %>%
  ggplot(., aes(x = Pathway, y = SAR, col = Route, group = interaction(Pathway, Route))) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.4, outliers = F) +
  geom_jitter(aes(col = Route), position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.3, jitter.height = 0), alpha = 0.1, size = 0.1) +
  facet_grid(cols = vars(network)) +
  scale_color_manual(values = env_pal) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "SAR by transmission route") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
    panel.border = element_rect(linewidth = 1),
    axis.line = element_blank()
    )

# Save figure
p_all = ggarrange(
  ggarrange(p1, p2, ncol = 2, labels = c("A", ""), common.legend = T, legend = "right", hjust = -1.5), 
  ggarrange(p3, p4, ncol = 1, nrow = 2, align = "hv", labels = c("B", "C"), hjust = -1.5),
  ncol = 1, nrow = 2, heights = c(1/3, 2/3))
p_all
ggsave("fig/paper/figure3.png", p_all, height = 7, width = 7)

## Plot global epidemic metrics for baseline scenario without intervention------
# Extinction probability
p_proba = stats_df %>%
  mutate(Extinction = ifelse(Global == 0, "Yes", "No")) %>%
  count(Pathway, network, Extinction) %>%
  nest(.by = network) %>%
  mutate(data = map(data, chisq_test_df)) %>%
  unnest(data) %>%
  mutate(y = y+0.05) %>%
  filter(p.adj <= 0.05)

stats_ext = stats_df %>% 
  group_by(network, Pathway) %>% 
  summarise(prop = sum(Global == 0) / n(), .groups = "drop") %>%
  mutate(Category = rep(names(pal)[1:3], n())[1:n()])

p1 = ggplot(stats_ext) +
  geom_bar(stat = "identity", aes(x = Pathway, y = prop, fill = network), width = 0.4) +
  geom_point(x = 0, y = 0, aes(col = Category)) +
  facet_grid(cols = vars(network)) +
  expand_limits(y = c(0,1)) +
  theme_classic() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 20, hjust = 1), 
        legend.position = "bottom", 
        panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
        panel.border = element_rect(linewidth = 1),
        axis.line = element_blank()) +
  scale_fill_manual(values = network_pal) +
  scale_color_manual(values = pal[1:3]) +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5),
         col = guide_legend(title.position="top", title.hjust = 0.5)) +
  labs(y = "Extinction probability", fill = "Network", col = "Category")

# p-value of extinction probability comparison across transmission scenarios
# within a network
stats_df %>%
  mutate(Extinction = ifelse(Global == 0, "Yes", "No")) %>%
  count(network, Extinction, Pathway) %>%
  pivot_wider(names_from = Extinction, values_from = n) %>%
  select(-Pathway) %>%
  nest(.by = network) %>% 
  mutate(chisq_p = map_dbl(data, ~chisq.test(.)$p.value))

# Epidemic duration
p_duration = stats_df %>%
  group_by(network) %>%
  wilcox_test(Epidemic_duration ~ Pathway, p.adjust.method = "BH") %>%
  add_xy_position(x = "Pathway") %>%
  filter(p.adj <= 0.05)

p2 = stats_df %>%
  filter(!is.na(Epidemic_duration)) %>%
  arrange(network, Pathway) %>%
  ggboxplot(., x = "Pathway", y = "Epidemic_duration", width = 0.4, outlier.shape = NA, col = "network") +
  geom_jitter(width = 0.2, height = 0, alpha = 0.1, aes(col = network), size = 0.1) +  
  facet_grid(cols = vars(network)) +
  labs(y = "Epidemic duration (in days)", col = "Network") +
  theme_classic() +
  scale_color_manual(values = network_pal) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 20, hjust = 1),
        panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
        panel.border = element_rect(linewidth = 1),
        axis.line = element_blank()) +
  expand_limits(ymin = 0)

# Kruskall-Wallis test to compare time to the peak across transmission pathways for 
# each transmission network
stats_df %>%
  nest(.by = network) %>%
  mutate(kruskall_wallis_p = map_dbl(data, ~kruskal.test(.$Epidemic_duration ~ .$Pathway)$p.value))

# Time to the peak with statistical tests
p_peak = stats_df %>%
  group_by(network) %>%
  wilcox_test(Time_to_peak ~ Pathway, p.adjust.method = "BH") %>%
  add_xy_position(x = "Pathway") %>%
  filter(p.adj <= 0.05)

p3 = stats_df %>%
  arrange(network, Pathway) %>%
  filter(!is.na(Time_to_peak)) %>%
  ggboxplot(., x = "Pathway", y = "Time_to_peak", width = 0.4, outlier.shape = NA, col = "network") +
  stat_pvalue_manual(p_peak, label = "p.adj.signif") +
  geom_jitter(width = 0.2, height = 0, alpha = 0.1, size = 0.1, aes(col = network)) +
  facet_grid(cols = vars(network)) +
  labs(y = "Time to the peak (in days)", col ="Network") +
  theme_classic() +
  scale_color_manual(values = c(network_pal, env_pal), drop = F) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 20, hjust = 1),
        panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
        panel.border = element_rect(linewidth = 1),
        axis.line = element_blank()) +
  expand_limits(ymin = 0)

# Kruskall-Wallis test to compare epidemic duration across transmission scenarios for 
# each transmission network
stats_df %>%
  nest(.by = network) %>%
  mutate(kruskall_wallis_p = map_dbl(data, ~kruskal.test(.$Time_to_peak ~ .$Pathway)$p.value))

# SAR difference between the two routes for the different individual categories
p4 = stats_df %>%
  pivot_longer(matches("[A-Z][a-z]+_Environment|[A-Z][a-z]+_Contact"), names_to = c("Category", "Route"), names_pattern = "(.*)_(.*)", values_to = "SAR") %>%
  group_by(Pathway, network, Category) %>%
  wilcox_test(SAR ~ Route, paired = T, ref.group = "Contact", alternative = "two.sided", detailed = T) %>%
  ggplot(., aes(x = Pathway, y = estimate, ymin = conf.low, ymax = conf.high, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.4) +
  geom_errorbar(position = "dodge", width = 0.4) +
  facet_grid(cols = vars(network)) +
  scale_fill_manual(values = pal) +
  theme_classic() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 20, hjust = 1),
        panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
        panel.border = element_rect(linewidth = 1),
        axis.line = element_blank()) +
  labs(y = "Absolute difference between\nshort- and long-range SAR")

stats_df %>%
  pivot_longer(c(Patient, HCW), values_to = "val", names_to = "cat") %>%
  group_by(network, cat) %>%
  summarise(
    m = round(mean(val*100)),
    me = round(median(val)*100), 
    iqr = round(IQR(val)*100),
    .groups = "drop"
  )

# Save figure 
p_all = ggarrange(p1, p2, p3, p4, 
                  ncol = 2, nrow = 2, align = "hv", 
                  common.legend = T, legend = "bottom",
                  labels = c("A", "B", "C", "D"))
p_all
ggsave("fig/paper/figure4.png", p_all, height = 6, width = 10)

## Plot cumulative incidence by infector-infectee pairs-------------------------
stats_df %>%
  filter(intervention == "None") %>%
  pivot_longer(c(I_patient_to_patient, I_patient_to_hcw, I_hcw_to_patient, I_hcw_to_hcw), names_to = "Pair", values_to = "I") %>%
  mutate(
    Pair = case_when(
      Pair == "I_patient_to_patient" ~ "Patient to Patient",
      Pair == "I_patient_to_hcw" ~ "Patient to HCW",
      Pair == "I_hcw_to_patient" ~ "HCW to Patient",
      Pair == "I_hcw_to_hcw" ~ "HCW to HCW"
    )
  ) %>% 
  left_join(., ind_numbers, by = c("network")) %>%
  mutate(SAR = case_when(
    grepl("to Patient", Pair) ~ I/n_patient,
    grepl("to HCW", Pair) ~ I/n_hcw
    )
  ) %>%
  ggplot(., aes(x = Pathway, y = SAR, fill = Pair)) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.6, 
               outlier.size = 0.5) +
  facet_grid(cols = vars(network)) +
  scale_fill_manual(values = pairs_pal) +
  scale_y_continuous(labels = scales::percent) +
  expand_limits(y = c(0,1)) +
  theme_classic() +
  theme(axis.title.x = element_blank(), legend.position = "bottom", 
        panel.grid.major = element_line(linewidth = 11/22, color = "grey90"),
        panel.border = element_rect(linewidth = 1),
        axis.line = element_blank()
        ) +
  labs(y = "Secondary attack rates", fill = "Transmission pairs")
ggsave("fig/baseline_scenario/pairs_transmission.png", height = 3, width = 8)
