if (dir.exists("~/r-libs")) {
  .libPaths("~/r-libs")
}

source("./R/load_packages.r")
source("./R/prep_data.r")

firerisk <- "~/data_input/fsf_fire_tract_summary.csv"
firerisk <- read.csv(firerisk) |>
    mutate(geoid10 = ifelse(nchar(fips) == 10, paste0("0", fips), as.character(fips)))

library(NLinteraction)

# join smoke exposure
hms_30day <- readRDS("./output/smoke_days_30days_ea_nona_addr.rds")
hms_60day <- readRDS("./output/smoke_days_60days_ea_nona_addr.rds")
hms_120day <- readRDS("./output/smoke_days_120days_ea_nona_addr.rds")
hms_365day <- readRDS("./output/smoke_days_365days_ea_nona_addr.rds")

hms_smoke_exp <-
  data.frame(
    smoke_d60 = hms_60day,
    smoke_d120 = hms_120day,
    smoke_d365 = hms_365day
  )

epr.gis <- epr.gis %>%
  mutate(gis_event_date = as.Date(gis_event_date, format = "%M/%d/%Y"))
epr_sf <- epr.gis %>%
    st_as_sf(coords = c("gis_longitude", "gis_latitude"), crs = 4326)

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

epr_allp_sf = epr.gis %>%
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



epr_console_sub <- epr_console %>%
  filter(gis_study_event == "current_address_exposome_a") %>%
  select(epr_number, flag_antidep_1p, 18:34, 36:593) %>%
  mutate(across(-1, ~ifelse(is.na(.), 0, as.numeric(.))))
epr_console_mv <- sapply(epr_console_sub, function(x) length(unique(x)) > 1)
epr_console_sub <- epr_console_sub[, ..epr_console_mv] #%>%
  # mutate(across(seq_len(ncol(.))[length(unique(.)) > 2], ~as.vector(scale(.))))

epr_console_sub_noe <-
  epr_console_sub %>%
  mutate(epr_number = as.character(epr_number)) %>%
  select(-contains("_E_"), -contains("EPL")) %>%
  left_join(epr_he_behav, by = "epr_number") %>%
  left_join(epr_he_family, by = "epr_number") 
  # mutate(geoid10 = as.character(gis_geoid10)) %>%
  # left_join(firerisk, by = c("geoid10" = "geoid10"))

pegs_all_join <- right_join(epr_console_sub_noe, epr_allp_hms)


# too many NAs (except for one record, actually)
# why were there too many NAs in GIS covariates?
pegs_X <- pegs_all_join |>
  st_drop_geometry() |>
  select(svi_2018_EP_GROUPQ, smoke_d365) |>
  as.matrix()
pegs_X[,1] <- rgamma(nrow(pegs_X), 1, 0.3)
summary(pegs_X)
pegs_C <- pegs_all_join |>
  st_drop_geometry() |>
  select(all_of(seq(53, 72))) |>
  as.matrix()
pegs_Y <- pegs_all_join |>
  st_drop_geometry() |>
  select(flag_antidep_1p) |>
  as.matrix()

# nltest
pegs_nlint <- NLint(pegs_Y, pegs_X, pegs_C, nChains = 4, nIter = 2000, nBurn = 500, thin = 5, ns = 2)
