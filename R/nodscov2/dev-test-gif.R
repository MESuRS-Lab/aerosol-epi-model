library(ggplot2)
library(gganimate)
library(gifski)
library(dplyr)
library(ggrepel)
library(tidyr)

# Function to distribute points in a circle
distribute_points <- function(n) {
  theta <- seq(0, 2 * pi, length.out = n + 1)[-1]
  radius <- 0.3  # Adjust as needed
  x <- radius * cos(theta)
  y <- radius * sin(theta)
  return(data.frame(offset_x = x, offset_y = y))
}

# Prepare data
chambres <- rooms %>%
  mutate(
    x = case_when(
      grepl("001", id) ~ id_room,
      grepl("Corridor", room) ~ 1,
      room == "Office" ~ 1,
      room == "Nursing station" ~ 1,
      grepl("Restroom", room) ~ (id_room %% 2) + 1,
      TRUE ~ NA_real_
    ),
    y = case_when(
      grepl("001", id) ~ 1,
      grepl("Corridor", room) ~ 2,
      room == "Office" ~ 3,
      room == "Nursing station" ~ 3,
      grepl("Restroom", room) ~ 4,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(localization = id_room) %>%
  distinct(id_room, .keep_all = TRUE) %>%
  select(localization, x, y)

all_data <- do.call(rbind, global_localization) %>% filter(time < 150)
all_data <- all_data %>%
  filter(localization != -1) %>%
  left_join(chambres, by = c("localization")) %>%
  left_join(admission, by = c("id")) %>%
  mutate(cat = ifelse(is.na(cat), "patient", cat)) %>%
  filter(cat != "patient")



# Determine the maximum number of individuals in the same localization at the same time
max_n <- all_data %>%
  group_by(localization, time) %>%
  summarise(n = n(), .groups = 'drop') %>%
  summarise(max_n = max(n)) %>%
  pull(max_n)

test <- all_data %>%
  group_by(localization, time) %>%
  mutate(count = n()) %>%
  mutate(offset = list(distribute_points(count))) %>%
  unnest(cols = offset) %>%
  mutate(x_adj = x + offset_x, 
         y_adj = y + offset_y) %>%
  ungroup()

# Plot
p <- ggplot(all_data, aes(x = x_adj, y = y_adj, color = cat, group = id)) +
  geom_point(size = 5) +
  geom_path(aes(group = id), alpha = 1) +
  geom_text_repel(aes(label = id), size = 3, box.padding = 0.5, point.padding = 0.5, max.overlaps = Inf) +
  scale_color_discrete(name = "Category") +
  scale_y_continuous(breaks = 1:4, labels = c("Chambre", "Corridor", "Office / Nursing station", "Restroom")) +
  scale_x_continuous(breaks = 1:sum(chambres$y == 1), labels = c(paste0("Chambre ", 1:sum(chambres$y == 1)))) +
  labs(title = 'Déplacement des individus au temps: {frame_time}', x = 'X', y = 'Type de pièce') +
  theme_minimal()

animation <- p +
  transition_time(time) +
  ease_aes('linear')

animate(animation, renderer = gifski_renderer("deplacement.gif"), width = 2000, height = 1333, duration = 20)
