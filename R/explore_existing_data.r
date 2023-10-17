library(pacman)
p_load(sf, terra, dplyr, purrr, data.table, vtable, skimr)
sf_use_s2(FALSE)


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
cbind(obj.names, obj.sizes) %>%
  write.csv("./output/object_list.csv", row.names = FALSE)

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

grep("meta", obj.names$value, value = TRUE) %>%
  lapply(get)
epr.gis <- epr.gis %>%
  mutate(gis_event_date = as.Date(gis_event_date, format = "%M/%d/%Y"))

is.date <- function(x) inherits(x, 'Date')

## summary tables
lapply(obj.names$value[!grepl("meta", obj.names$value)],
  function(x) {
    get(x) %>%
    mutate_if(is.character, ~as.factor(.)) %>%
    mutate_if(is.date, ~as.numeric(.)) %>%
    skimr::skim()
  })
# epr.svi %>%
#     vtable::sumtable(vars = names(.), 
#                      add.median = TRUE,
#                      file = "./output/summary_svi.html")

summary(epr.he)


## Location maps
p_load(ggplot2, tigris)
us_cnty <- tigris::counties(year = 2020) %>%
    filter(!STATEFP %in% c("02", "72", "78", "15", "60", "69", "66")) %>%
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
ggsave("./output/distr.png", ggindloc, width = 15, height = 9, units = "in", dpi = 300)






## Detailed metadata tables -- failed
## Metadata tables are seemingly generated from SAS
codebooks <- list.files(paste0(path_pegs, "Data_Freezes/latest"),
            pattern = "*_codebook_([0-9]+{1,2})*.*.pdf$",
            recursive = TRUE,
            full.names = TRUE)
codebooks <- codebooks[!grepl("appendix", codebooks)]

codebook_tables <- lapply(codebooks[2], 
    function(x) tabulizer::extract_tables(x, method = "lattice", output = "data.frame") %>% map(as_tibble))

pdf_tables <- pdf_data |> 
  extract_tables(output = "data.frame") |>
  map(as_tibble)


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

## true number of participants, not total recorded locations
## 18462
length(unique(epr.gis$epr_number))

epr_allp_sf = epr.gis %>%
  filter(epr_number %in% as.integer(epr_number_allpresence)) %>%
  st_as_sf(coords = c("gis_longitude", "gis_latitude"), crs = 4326)
epr_allp_sf[us_cnty %>% filter(STATEFP == 37),]

ggindallploc = ggplot() +
  geom_sf(data = us_cnty, linewidth = 0.8) +
  geom_sf(data = epr_allp_sf, pch = 19, cex = 0.6, col = 'red')



### ATC
epr_atc_psych <-
  epr.atc %>%
  as_tibble %>%
  filter(grepl("DEPRESS", atc_level3_name)) %>%
  mutate(flag_antidep_1p = 1) %>%
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

## mental health status: u210a-k, u211-u214.
epr_he_mental <-
  epr.he %>%
  as_tibble %>%
  select(1, starts_with("he_u21")) %>%
  mutate(epr_number = as.integer(epr_number)) %>%
  mutate(across(-1, ~as.integer(.))) %>%
  mutate(across(-1, ~ifelse(is.na(.), -999, .))) %>%
  left_join(epr_atc_psych %>% select(1, flag_antidep_1p), by = "epr_number") %>%
  left_join(epr_atc_psych2p %>% select(1, flag_antidep_2p) %>% unique, by = "epr_number") %>%
  mutate(he_flsum = rowSums(across(he_u210a_unusual_hyper_PARQ:he_u214_bipolar))) %>%
  select(1, he_flsum, flag_antidep_1p, flag_antidep_2p, 2:(ncol(.)-4)) %>%
  mutate(across(starts_with("flag_antidep"), ~ifelse(is.na(.), 0, .)))

epr_he_mental

table(epr_he_mental$flag_antidep_1p, epr_he_mental$he_flsum)
table(epr_he_mental$flag_antidep_2p, epr_he_mental$he_flsum)
# questionnaire gives high frequency of self-reported negative mental health
# in participants who do not take antidepressants
# result of positiveness?
