# library(pacman)
# p_load(sf, terra, dplyr, purrr, data.table, vtable, skimr, rmapshaper, tigris)
source("./R/load_packages.r")
# Uncomment the line below when running the code on triton
# .libPaths("/ddn/gs1/home/songi2/r-libs")
# sf_use_s2(FALSE)
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
# cbind(obj.names, obj.sizes) %>%
#   write.csv("./output/object_list.csv", row.names = FALSE)

## Key info & glossary
## epr_number: EPR (previous name of PEGS) unique id
## Data value flags
## Basically all data values are supposed to character/factor
##    .M is missing
##    .N is nullified by QC
##    .S is skipped
##    -88888 is don't know
## ATC: anatomical therapeutic chemical
## gis: coordinates (latest)
## [dia]betes
## [ecz]ema
## [e]xposome survey part [a, b]
##      [ea]: external exposome
##      [eb]: internal exposome --> deriving atc
##      Supplements/medicine intake (current/history)
##      Life history (housing, water source, urbanity, pet, surrounding, hobby, sunscreen)
##          Recreational/Occupational/Residential exposure as confounding/mediating factors
## [h]ealth & [e]xposure study: multiton, reported/known health status (*mental*), 
##      employment, regular exposure to stressors
##      mental health status: u210a-k, u211-u214.
## eji
## svi
## hazards (environmental)
## bcbb.map: changes over multiple surveys, gender identification, withdrawal, death
##    has_address_d: last known address indicator

# grep("meta", obj.names$value, value = TRUE) %>%
#   lapply(get)
epr.gis <- epr.gis %>%
  mutate(gis_event_date = as.Date(gis_event_date, format = "%M/%d/%Y"))
# Timeflag incompleteness
# Spatial coord incompleteness


is.date <- function(x) inherits(x, 'Date')

## summary tables
# sink(file = "./output/summary_stats_objects.txt")
# options(width = 240)
# lapply(obj.names$value[!grepl("meta", obj.names$value)],
#   function(x) {
#     print(sprintf("Summary statistics of %s", x))
#     get(x) %>%
#     mutate_if(is.character, ~as.factor(.)) %>%
#     mutate_if(is.date, ~as.numeric(.)) %>%
#     skimr::skim()
#   })
# sink(NULL)
# options(width = 120)
# summary(epr.he)


## Location maps
# p_load(ggplot2, tigris)
us_cnty <- tigris::counties(year = 2020) %>%
    dplyr::filter(!STATEFP %in% c("02", "72", "78", "15", "60", "69", "66")) %>%
    st_transform(4326)
epr_sf <- epr.gis %>%
    st_as_sf(coords = c("gis_longitude", "gis_latitude"), crs = 4326)
epr.gis$gis_state %>% table

epr_sfs <- epr_sf[us_cnty,] %>%
    st_jitter(0.3)
ggindloc <- ggplot() +
    geom_sf(data = us_cnty, linewidth = 0.8) +
    geom_sf(data = epr_sfs, pch = 19, cex = 0.6, col = 'red') +
    labs(caption = "Coordinates were jittered.")
# ggsave("./output/distr.png", ggindloc, width = 15, height = 9, units = "in", dpi = 300)






## Detailed metadata tables -- failed
## Metadata tables are seemingly generated from SAS
# codebooks <- list.files(paste0(path_pegs, "Data_Freezes/latest"),
#             pattern = "*_codebook_([0-9]+{1,2})*.*.pdf$",
#             recursive = TRUE,
#             full.names = TRUE)
# codebooks <- codebooks[!grepl("appendix", codebooks)]

# codebook_tables <- lapply(codebooks[2], 
#     function(x) tabulizer::extract_tables(x, method = "lattice", output = "data.frame") %>% map(as_tibble))

# pdf_tables <- pdf_data |> 
#   extract_tables(output = "data.frame") |>
#   map(as_tibble)


## check the number of participants who were "all available" in surveys
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
# 1,897 participants are available across attribute tables with 1000+ rows

## Key reference dataset!!! ####
## true number of participants, not total recorded locations
## 18462
length(unique(epr.gis$epr_number))

epr_allp_sf = epr.gis %>%
  dplyr::filter(epr_number %in% as.integer(epr_number_allpresence)) %>%
  st_as_sf(coords = c("gis_longitude", "gis_latitude"), crs = 4326)
# In NC
epr_allp_sf[us_cnty %>% dplyr::filter(STATEFP == 37),]
# types of enrolment
unique(epr_allp_sf$gis_study_event)
# 1,648 participants have full coordinate information
# across five key events (exposome, enrollment, he, ...)
# 1,849 participants have four key events
# Exposome a and b surveys were done (almost) at the same time
epr_allp_sf %>%
  st_drop_geometry() %>%
  dplyr::group_by(epr_number) %>%
  dplyr::summarize(Nevents = length(unique(gis_study_event))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(Nevents >= 4)

epr_allp_sf %>%
  st_drop_geometry() %>%
  dplyr::group_by(epr_number) %>%
  dplyr::summarize(Nevents = length(unique(gis_study_event))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(Nevents >= 4)

us_cnty_s <- rmapshaper::ms_simplify(us_cnty, 0.1, keep_shapes = TRUE)
epr_allp_sf_sub <- epr_allp_sf[us_cnty_s,]
ggindallploc = ggplot() +
  geom_sf(data = us_cnty_s, linewidth = 0.8) +
  geom_sf(data = epr_allp_sf_sub, pch = 19, cex = 0.6, col = 'red')
ggindallploc


### ATC
epr_atc_psych <-
  epr.atc %>%
  as_tibble %>%
  # filter(grepl("DEPRESS", atc_level3_name)) %>%
  mutate(flag_antidep_1p = ifelse(grepl("DEPRESS", atc_level3_name), 1, 0)) %>%
  mutate(epr_number = as.integer(epr_number))
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
  mutate(epr_number = as.integer(epr_number))
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

as.data.frame(epr_atc_psych[501:520,])


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

names(epr.he)



## Look into the temporal consistency and order
## Strong assumption: no potential reverse causality
## Mental health status should follow our exposure of interest

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
  mutate(epr_number = as.integer(epr_number)) %>%
  mutate(across(-1, ~as.integer(.))) %>%
  mutate(across(-1, ~ifelse(is.na(.), -999, .)))

epr_eb_q_psych <- inner_join(epr_eb_mental, epr_atc_psych)
antidep_cope <- table(epr_eb_q_psych$flag_antidep_1p, epr_eb_q_psych$eb_e095_unable_to_cope)
antidep_cope/rowSums(antidep_cope)

antidep_cope <- table(epr_eb_q_psych$flag_antidep_1p, epr_eb_q_psych$eb_e099_cannot_overcome)
antidep_cope/rowSums(antidep_cope)

antidep_cope <- table(epr_eb_q_psych$flag_antidep_1p, epr_eb_q_psych$eb_e094_going_your_way)
antidep_cope/rowSums(antidep_cope)


# he survey:
# s177-s193: health behaviors (drinking, sleeping, ...)
# t192-199: family history (well-off, etc.) in youth
# t200-207: current family history (their own)
# u208*: self perception and social behaviors
# u209: compounded occurrences of u208* entries
# u210: attribution to difficulties
# u211-212: bipolar history in blood siblings
## mental health status: u210a-k, u211-u214.
epr_he_mental <-
  epr.he %>%
  as_tibble %>%
  select(1, starts_with("he_u20"), starts_with("he_u21")) %>%
  mutate(epr_number = as.integer(epr_number)) %>%
  mutate(across(-1, ~as.integer(.))) %>%
  mutate(across(-1, ~ifelse(is.na(.), -999, .))) %>%
  left_join(epr_atc_psych %>% select(1, flag_antidep_1p), by = "epr_number") %>%
  left_join(epr_atc_psych2p %>% select(1, flag_antidep_2p) %>% unique, by = "epr_number") %>%
  mutate(he_flsum = rowSums(across(he_u210a_unusual_hyper_PARQ:he_u214_bipolar))) %>%
  select(1, he_flsum, flag_antidep_1p, flag_antidep_2p, 2:(ncol(.)-3)) %>%
  mutate(across(starts_with("flag_antidep"), ~ifelse(is.na(.), 0, .)))

epr_he_behav <-
  epr.he %>%
  as_tibble %>%
  select(1, starts_with("he_s17"), starts_with("he_s19")) %>%
  mutate(epr_number = as.integer(epr_number)) %>%
  mutate(across(-1, ~as.integer(.))) %>%
  mutate(across(-1, ~ifelse(is.na(.), -999, .)))

epr_he_family <-
  epr.he %>%
  as_tibble %>%
  select(1, starts_with("he_t19")) %>%
  mutate(epr_number = as.integer(epr_number)) %>%
  mutate(across(-1, ~as.integer(.))) %>%
  mutate(across(-1, ~ifelse(is.na(.), -999, .)))

skim(epr_he_mental)
skim(epr_he_behav)
skim(epr_he_family)

names(epr_he_mental)

table(epr_he_mental$flag_antidep_1p, epr_he_mental$he_flsum)
table(epr_he_mental$flag_antidep_2p, epr_he_mental$he_flsum)
# questionnaire gives high frequency of self-reported negative mental health
# in participants who do not take antidepressants
# indicative of self-recognition of negative mental health status or
# the willingness to get treatment


# Checking temporal inconsistency across survey/event sequences ####
epr_gis <- epr.gis %>%
  mutate(event_code = case_when(
    gis_study_event == "longest_lived_child" ~ "A",
    gis_study_event == "longest_lived_adult" ~ "B",
    gis_study_event == "enrollment" ~ "C",
    gis_study_event == "health_and_exposure" ~ "D",
    gis_study_event == "current_address_exposome_a" ~ "E",
    TRUE ~ NA
  )) %>%
  arrange(epr_number, event_code)

epr_gis %>%
  group_by(epr_number) %>%
  #mutate(lag_evdate = lag(gis_event_date)) %>%
  #filter(!is.na(lag_evdate)) %>%
  mutate(evdate_diff = c(NA, diff(gis_event_date))) %>%
  filter(!is.na(evdate_diff)) %>%
  filter(evdate_diff < 0) %>%
  ungroup() %>%
  # .$gis_study_event %>%
  # unique()
  filter(gis_study_event == "health_and_exposure") %>%
  select(epr_number, evdate_diff)
## in some cases, health_and_exposure completion date precedes the enrolment
## the number is negligible (18/18K+, <0.1 %)
## heard that it is not unusual


## EB
## eb_b and eb_e: medication and mental health status
