
pkgs <- c("sf", "terra", "dplyr", "purrr", "tidytable",
          "data.table", "vtable", "skimr", "rmapshaper", "tigris")
invisible(sapply(pkgs, library, quietly = TRUE, character.only = TRUE))
sf_use_s2(FALSE)
Sys.setenv("TIGRIS_CACHE_DIR" = sprintf("%s/tigris_cache/", Sys.getenv("HOME")))

## data path determination
## compute_mode 1 (wine mount), compute_mode 2 (hpc compute node), 3 (container internal)
COMPUTE_MODE <- 2
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

## GIS record history
epr.gis.hist <- epr.gis |>
  mutate(across(ends_with("_date"), ~as.POSIXct(as.character(.), format = "%Y-%m-%d"))) |>
  mutate(gis_duration = difftime(gis_end_date, gis_start_date, units = "days"))

epr.gis.hist |>
  group_by(gis_study_event) |>
  summarize(mean_duration = mean(gis_duration, na.rm = T))
epr.gis.hist |>
  group_by(gis_study_event) |>
  summarize(N = n())
# How many participants reported child/adult residence?
epr.gis.hist |>
  filter(endsWith(gis_study_event, "exposome_a")) |>
  group_by(epr_number) |>
  summarize(N = (n() == 3)) |>
  ungroup() |>
  _$N |>
  sum()

## distance calc (+ excluding one-per-participant records)
## + indiv. char
epr.bcbb.ind <- epr.bcbb.map |>
  select(1:42, 57:65, 68:92) |>
  mutate(epr_number = as.character(epr_number))

epr_reloc <-
  epr.gis.hist |>
  mutate(epr_number = as.character(epr_number)) |>
  group_by(epr_number) |>
  filter(n() != 1L) |>
  arrange(epr_number, gis_start_date) |>
  mutate(gis_lat_lag = lag(gis_latitude),
         gis_lon_lag = lag(gis_longitude)) |>
  rowwise() |>
  mutate(
    reloc_dist =
    geosphere::distGeo(
      c(gis_longitude, gis_latitude),
      c(gis_lon_lag, gis_lat_lag)
    )
  ) |>
  ungroup() |>
  left_join(epr.bcbb.ind, by = c("epr_number"))

# 660 withdrawn
table(epr_reloc$withdrawal_date) |> _[-1] |> sum()
table(epr_reloc$sex_derived)

#
epr.gis.he <- epr.gis %>%
  mutate(epr_number = as.character(epr_number)) %>%
  filter(grepl("health_and_exposure", gis_study_event)) %>%
  inner_join(epr.he, by = "epr_number") %>%
  filter(gis_quality == "residential")



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

# eb survey:
# 
# B032a-B041a: past one year and medication history (anxiety, depression)
# B052-053: current status, continuing, and suspended are distinguishable
# E090-E099: stress (E): past month
# G group: sleep quality + behavior / past month
# K group: genetic history (proxy of potential concerns)
epr_eb_mental <-
  epr.eb %>%
  as_tibble %>%
  select(1, starts_with("eb_e"), starts_with("eb_b03"), starts_with("eb_b04"),
    starts_with("eb_b05")) %>%
  mutate(epr_number = (epr_number)) %>%
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
epr_gis_ea <- epr_gis %>%
  group_by(epr_number) %>%
  filter(event_code == "E") %>%
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
  # 12062023: valid values only
  filter(!is.na(hazards_pm_25_ugm3)) %>%
  mutate(epr_number = as.character(epr_number)) %>%
  arrange(epr_number, event_code) #%>%
  #group_by(epr_number) %>%
  #filter(event_code == "E") %>%
  #ungroup()

# mental health related outcomes
epr_he_mental <-
  epr.he %>%
  as_tibble %>%
  mutate(he_date_completed = as.Date(he_date_completed)) %>%
  group_by(epr_number) %>%
  filter(he_date_completed == max(he_date_completed)) %>%
  ungroup() %>%
  select(1, starts_with("he_u20"), starts_with("he_u21")) %>%
  mutate(epr_number = (epr_number)) %>%
  mutate(across(-1, ~as.integer(.))) %>%
  mutate(across(-1, ~ifelse(is.na(.), -999, .))) %>%
  # left_join(epr_atc_psych %>% select(1, flag_antidep_1p), by = "epr_number") %>%
  # left_join(epr_atc_psych2p %>% select(1, flag_antidep_2p) %>% unique, by = "epr_number") %>%
  mutate(he_flsum = rowSums(across(he_u210a_unusual_hyper_PARQ:he_u214_bipolar))) %>%
  # select(1, he_flsum, flag_antidep_1p, flag_antidep_2p, 2:(ncol(.)-3)) %>%
  mutate(across(starts_with("flag_antidep"), ~ifelse(is.na(.), 0, .)))

epr_he_behav <-
  epr.he %>%
  as_tibble %>%
  select(1, starts_with("he_s17"), starts_with("he_s19")) %>%
  mutate(epr_number = (epr_number)) %>%
  mutate(across(-1, ~as.integer(.))) %>%
  mutate(across(-1, ~ifelse(is.na(.), -999, .)))

epr_he_family <-
  epr.he %>%
  as_tibble %>%
  select(1, starts_with("he_t19")) %>%
  mutate(epr_number = (epr_number)) %>%
  mutate(across(-1, ~as.integer(.))) %>%
  mutate(across(-1, ~ifelse(is.na(.), -999, .)))



## Basic personal information
epr_bcbb <- epr.bcbb.map |>
  mutate(across(-1:-2, ~as.integer(.))) |>
  transmute(
    epr_number = epr_number,
    gender_legacy = factor(gender_legacy, levels = 0:1, labels = c("Male", "Female")),
    race = factor(race, levels = 1:6,
        labels = c("AIAN", "Asian", "Black", "NHPI", "White", "Multiple")),
    birth_date = birth_date#as.Date(birth_date, format = "%m/%d/%y")
  )
epr_bcbb


## join
eprs <- 
  cbind(obj.names, obj.sizes) %>%
  filter(V1 > 1000) %>%
  .$value %>%
  grep("epr", ., value = TRUE) %>%
  .[!grepl("meta", .)]
eprs_vals <- eprs %>%
  lapply(get) %>%
  lapply(as_tibble) %>%
  lapply(function(x) x %>% select(epr_number)) %>%
  lapply(unlist) 

epr_number_allpresence <- Reduce(intersect, eprs_vals)

epr_allp_sf <- epr.gis %>%
  dplyr::filter(epr_number %in% as.integer(epr_number_allpresence)) %>%
  st_as_sf(coords = c("gis_longitude", "gis_latitude"), crs = 4326)


epr_allp_sf_ea <-
  epr_allp_sf %>%
  dplyr::mutate(gis_event_date = as.Date(gis_event_date, format = "%m/%d/%Y")) %>%
  dplyr::filter(gis_study_event == "current_address_exposome_a") %>%
  dplyr::filter(!gis_state %in% c("HI", "AK", "GU", "VI", "PR")) %>%
  dplyr::filter(!is.na(gis_event_date)) %>%
  sf::st_transform("EPSG:5070")

epr_allp_hms <- cbind(epr_allp_sf_ea, hms_smoke_exp) |>
  mutate(epr_number = as.character(epr_number))
#     mutate(geoid10 = as.character(gis_geoid10)) |>
#     left_join(firerisk, by = c("geoid10" = "geoid10"))




epr_atc_psychs <- epr_atc_psych %>%
  select(epr_number) %>%
  unique() %>%
  as.data.table()
epr_console <-
  epr_allp_hms %>%
  # epr_atc_psychs %>%
  # left_join(epr_gis_ea) %>%
  # left_join(epr_atc_psych %>% select(1, flag_antidep_1p) %>% unique()) %>%
  inner_join(epr_eb_mental) %>%
  inner_join(epr_he_mental) %>%
  inner_join(epr_eji) %>%
  inner_join(epr_svi) %>%
  inner_join(epr_hazards) %>%
  inner_join(epr_earth) %>%
  inner_join(epr_bcbb)


