library(pacman)
p_load(sf, terra, dplyr, purrr, tidytable, data.table, vtable, skimr, rmapshaper, tigris)
sf_use_s2(FALSE)
Sys.setenv("TIGRIS_CACHE_DIR" = sprintf("%s/tigris_cache/", Sys.getenv("HOME")))

## data path determination
## compute_mode 1 (wine mount), compute_mode 2 (hpc compute node), 3 (container internal)
COMPUTE_MODE <- 1
path_pegs <-
  ifelse(COMPUTE_MODE == 1,
  "/Volumes/PEGS/",
  ifelse(COMPUTE_MODE == 2,
  "/ddn/gs1/project/controlled/PEGS/",
  ifelse(COMPUTE_MODE == 3,
  "/opt/", stop("COMPUTE_MODE should be one of 1, 2, or 3.\n"))))


list.files(paste0(path_pegs, "Data_Freezes/latest"), 
        pattern = "*.(csv|RData|rdata|rds|RDS|CSV|txt)$",
        recursive = TRUE, full.names = TRUE)


rdatas <- list.files(
  path = paste0(path_pegs, "Data_Freezes/latest"),
  pattern = "*.(RData)$",
  recursive = TRUE,
  full.names = TRUE)

for (rdata in rdatas) {
  load(rdata)
}

obj.list <- sapply(rdatas, load)
obj.names <- obj.list %>%
  t %>%
  data.frame(fpath = rownames(.), .) %>%
  tidyr::pivot_longer(cols = 2:3)
obj.sizes <- obj.list %>%
  lapply(function(x) Reduce(rbind, sapply(x, function(y) dim(get(y)), simplify = FALSE))) %>%
  Reduce(rbind, .) %>%
  as.data.frame


epr.gis <- epr.gis %>%
  mutate(gis_event_date = as.Date(gis_event_date, format = "%M/%d/%Y"))
is.date <- function(x) inherits(x, 'Date')


## outcomes
### ATC
epr_atc_psych <-
  epr.atc %>%
  as_tibble %>%
  # filter(grepl("DEPRESS", atc_level3_name)) %>%
  mutate(flag_antidep_1p = ifelse(grepl("DEPRESS", atc_level3_name), 1, 0),
         epr_number = as.character(epr_number))
length(unique(epr_atc_psych$epr_number))
unique(epr_atc_psych$atc_level3_name)
## 563 participants reported that they took antidepressants

## How many participants took 2+ antidepressants?
epr_atc_psych2p <-
  epr.atc %>%
  as_tibble %>%
  filter(grepl("DEPRESS", atc_level3_name)) %>%
  group_by(epr_number) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  mutate(flag_antidep_2p = 1) %>%
  mutate(epr_number = as.character(epr_number))
length(unique(epr_atc_psych2p$epr_number))
unique(epr_atc_psych2p$atc_level3_name)
## 86 participants took 2+ types of antidepressants
## todo: is the list cumulative or the full list at the survey period?

epr_atc_stim <-
  epr.atc %>%
  as_tibble %>%
  filter(grepl("PSYCHOSTIMUL", atc_level3_name))
length(unique(epr_atc_stim$epr_number))
unique(epr_atc_stim$atc_level3_name)
## 93 participants reported that they took psychostimulants (for ADHD)

epr_atc_sl <-
  epr.atc %>%
  as_tibble %>%
  filter(grepl("HYPNOTIC", atc_level3_name))
length(unique(epr_atc_sl$epr_number))
unique(epr_atc_sl$atc_level3_name)
## 79 participants took hypnotics (benzodiazepine; sleeping pill)

epr.atc %>%
  as_tibble %>%
  filter(grepl("venlafaxine", atc_level5_name))

## mental health status: u210a-k, u211-u214.
epr_he_mental <-
  epr.he %>%
  as_tibble %>%
  select(1, starts_with("he_u21")) %>%
  mutate(epr_number = as.character(epr_number)) %>%
  mutate(across(-1, ~as.integer(.))) %>%
  mutate(across(-1, ~ifelse(is.na(.), -999, .))) %>%
#   left_join(epr_atc_psych %>% select(1, flag_antidep_1p), by = "epr_number") %>%
#   left_join(epr_atc_psych2p %>% select(1, flag_antidep_2p) %>% unique, by = "epr_number") %>%
  mutate(he_flsum = rowSums(across(he_u210a_unusual_hyper_PARQ:he_u214_bipolar)))# %>%
#   select(1, he_flsum, flag_antidep_1p, flag_antidep_2p, 2:(ncol(.)-4)) %>%
#   mutate(across(starts_with("flag_antidep"), ~ifelse(is.na(.), 0, .)))

epr_eb_mental <-
  epr.eb %>%
  as_tibble %>%
  select(1, starts_with("eb_e")) %>%
  mutate(epr_number = as.character(epr_number)) %>%
  mutate(across(-1, ~as.integer(.))) %>%
  mutate(across(-1, ~ifelse(is.na(.), -999, .)))



### GIS / EJI
# temporal inconsistency check
epr_gis <- epr.gis %>%
  mutate(event_code = case_when(
    gis_study_event == "longest_lived_child" ~ "A",
    gis_study_event == "longest_lived_adult" ~ "B",
    gis_study_event == "enrollment" ~ "C",
    gis_study_event == "health_and_exposure" ~ "D",
    gis_study_event == "current_address_exposome_a" ~ "E",
    TRUE ~ NA
  )) %>%
  arrange(epr_number, event_code) %>%
  mutate(epr_number = as.character(epr_number))

epr_gis_last <- epr_gis %>%
  group_by(epr_number) %>%
  filter(event_code == event_code[n()]) %>%
  ungroup() %>%
  filter(!is.na(epr_number))

epr_eji <- epr.eji %>%
  mutate(event_code = case_when(
    gis_study_event == "longest_lived_child" ~ "A",
    gis_study_event == "longest_lived_adult" ~ "B",
    gis_study_event == "enrollment" ~ "C",
    gis_study_event == "health_and_exposure" ~ "D",
    gis_study_event == "current_address_exposome_a" ~ "E",
    TRUE ~ NA
  )) %>%
  arrange(epr_number, event_code) %>%
  group_by(epr_number) %>%
  filter(event_code == event_code[n()]) %>%
  ungroup()
epr_eji$epr_number <- as.character(epr_eji$epr_number)


epr_earth <- epr.earthdata %>%
  mutate(event_code = case_when(
    gis_study_event == "longest_lived_child" ~ "A",
    gis_study_event == "longest_lived_adult" ~ "B",
    gis_study_event == "enrollment" ~ "C",
    gis_study_event == "health_and_exposure" ~ "D",
    gis_study_event == "current_address_exposome_a" ~ "E",
    TRUE ~ NA
  )) %>%
  arrange(epr_number, event_code) %>%
  mutate(epr_number = as.character(epr_number)) %>%
  group_by(epr_number) %>%
  filter(event_code == event_code[n()]) %>%
  ungroup()


epr_svi <- epr.svi %>%
  mutate(event_code = case_when(
    gis_study_event == "longest_lived_child" ~ "A",
    gis_study_event == "longest_lived_adult" ~ "B",
    gis_study_event == "enrollment" ~ "C",
    gis_study_event == "health_and_exposure" ~ "D",
    gis_study_event == "current_address_exposome_a" ~ "E",
    TRUE ~ NA
  ),
    epr_number = as.character(epr_number)) %>%
  arrange(epr_number, event_code) %>%
  group_by(epr_number) %>%
  filter(event_code == event_code[n()]) %>%
  ungroup()


epr_hazards <- epr.hazards %>%
  mutate(event_code = case_when(
    gis_study_event == "longest_lived_child" ~ "A",
    gis_study_event == "longest_lived_adult" ~ "B",
    gis_study_event == "enrollment" ~ "C",
    gis_study_event == "health_and_exposure" ~ "D",
    gis_study_event == "current_address_exposome_a" ~ "E",
    TRUE ~ NA
  )) %>%
  mutate(epr_number = as.character(epr_number)) %>%
  arrange(epr_number, event_code) %>%
  group_by(epr_number) %>%
  filter(event_code == event_code[n()]) %>%
  ungroup()


## join
epr_atc_psychs <- epr_atc_psych %>%
  select(epr_number) %>%
  unique() %>%
  as.data.table()
epr_console <- 
  epr_atc_psychs %>%
  left_join(epr_gis_last) %>%
  left_join(epr_atc_psych %>% select(1, flag_antidep_1p) %>% unique()) %>%
  left_join(epr_he_mental) %>%
  left_join(epr_eji) %>%
  left_join(epr_svi) %>%
  left_join(epr_hazards) %>%
  left_join(epr_earth)


