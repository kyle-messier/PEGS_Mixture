if (dir.exists("~/r-libs")) {
  .libPaths("~/r-libs")
}

source("./R/load_packages.r")
source("./R/explore_existing_data.r")

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

epr_allp_hms <- cbind(epr_allp_sf_ea, hms_smoke_exp) |>
    mutate(geoid10 = as.character(gis_geoid10)) |>
    left_join(firerisk, by = c("geoid10" = "geoid10"))


