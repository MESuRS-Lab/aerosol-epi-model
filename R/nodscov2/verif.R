library(tidyverse)
admission_p = read.csv2("data/data-synthetic-graphs/input/admission_poincare.csv")
admission_h = read.csv2("data/data-synthetic-graphs/input/admission_herriot.csv")
admission_old = read.csv2("data/data-synthetic-graphs/input/toy_admission.csv")

interaction_p = read.csv2("data/data-synthetic-graphs/input/interactions_poincare.csv")
interaction_h = read.csv2("data/data-synthetic-graphs/input/interactions_herriot.csv")
interaction_old = read.csv2("data/data-synthetic-graphs/input/toy_mat_ctc.csv")

agenda_p = read.csv2("data/data-synthetic-graphs/input/agenda_poincare.csv")
agenda_h = read.csv2("data/data-synthetic-graphs/input/agenda_herriot.csv")
agenda_old = read.csv2("data/data-synthetic-graphs/input/toy_agenda.csv")

max(admission$lastDate)

c(min(admission_p$lastDate), max(admission_p$lastDate))

c(min(admission_h$lastDate), max(admission_h$lastDate))

all(admission_h$id %in% c(interaction_h$from, interaction_h$to))
all(c(interaction_h$from, interaction_h$to) %in% admission_h$id)
all(agenda_h$id %in% admission_h$id)

all(c(interaction_p$from, interaction_p$to) %in% admission_p$id)
all(agenda_p$id %in% admission_p$id)

identical(colnames(admission_h), colnames(admission_old))
identical(colnames(admission_p), colnames(admission_old))

identical(colnames(agenda_h), colnames(agenda_old))
identical(colnames(agenda_p), colnames(agenda_old))

identical(colnames(interaction_h), colnames(interaction_old))
identical(colnames(interaction_p), colnames(interaction_old))

identical(agenda_p, agenda_old)
identical(interaction_p, interaction_old)
identical(admission_p, admission_old)

head(admission_p)
head(admission_old)

admission_p %>%
  anti_join(.,admission_old) %>%
  count(status, hospitalization, cat, sex)

admission_p %>%
  inner_join(.,admission_old) %>%
  filter(status == "PA")

i = "PA-001-0117-M-R"
interaction_p %>% filter(from == i | to == i) 
interaction_old %>% summarise(min(length), max(length))
interaction_p %>% summarise(min(length), max(length))

interaction_p %>% 
  anti_join(., interaction_old)

interaction_old %>% 
  anti_join(., interaction_p)



