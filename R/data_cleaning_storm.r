## Load packages ####
pkgs <- c("data.table", "sf", "stars", "sftime", "terra")
invisible(sapply(pkgs, library, quietly = TRUE, character.only = TRUE))
options(sf_use_s2 = FALSE)

#
library(rvest)
turl <- "https://www.ncei.noaa.gov/pub/data/swdi/stormevents/csvfiles/"
thtml <- read_html(turl)
tlinks <- thtml |>
  html_table()
tlinks <- tlinks[[1]]$Name
tlinks <- grep("StormEvents", tlinks, value = TRUE)
tlinkspre2000 <- grep("(d20*)", tlinks, value = TRUE, invert = TRUE)
tlinkspre2000

tlinkspre2000e <- sprintf("%s%s", turl, tlinkspre2000)

# for (i in tlinkspre2000) {
#   download.file(
#     sprintf("%s%s", turl, i),
#     sprintf("%s%s", "./input/noaastorm/", i),
#     mode = "wb"
#   )
#   Sys.sleep(0.5)
# }

# types
# details & fatalities -> 1950-2023
# locations -> 1996-2023 (197206 is the start, only with two records)

quicklist <- \(d, q) list.files(d, q, FALSE, TRUE)
dir_input <- "input/noaastorm/csv"
csv_details <- quicklist(dir_input, "details") |>
  lapply(read.csv) |>
  data.table::rbindlist()
csv_fatal <- quicklist(dir_input, "fatalities") |>
  lapply(read.csv) |>
  data.table::rbindlist()
csv_locs <- quicklist(dir_input, "locations") |>
  lapply(read.csv) |>
  data.table::rbindlist()

# fwrite(csv_details, "input/noaastorm/storm_details_all.csv.gz")
# fwrite(csv_fatal, "input/noaastorm/storm_fatal_all.csv.gz")
# fwrite(csv_locs, "input/noaastorm/storm_locs_all.csv.gz")

# load file
csv_details <- readRDS("output/noaa_stormdb.rds")

unique(csv_details$EVENT_TYPE)
csv_details[EVENT_TYPE == "Heat", ]
csv_details[EVENT_TYPE == "Flood", ]
unique(csv_details$DAMAGE_PROPERTY)

# property damage cleaning
# TODO: check if as-then or standardized dollar values
csv_details2 <- csv_details[,
list(damage_property = substr(DAMAGE_PROPERTY, 1, nchar(DAMAGE_PROPERTY) - 1),
     damage_property_unit = substr(DAMAGE_PROPERTY, nchar(DAMAGE_PROPERTY), nchar(DAMAGE_PROPERTY)),
     damage_property = as.numeric(damage_property),
     damage_property_unit = ifelse())]

# unusual cases: property unit ends with H or h (houses?)
csv_details |>
  dplyr::filter(endsWith(DAMAGE_PROPERTY, "H"))
csv_details |>
  dplyr::filter(!grepl("[A-Za-z]$", DAMAGE_PROPERTY) & !DAMAGE_PROPERTY %in% c("", "0"))

# temporal trends
csv_details[EVENT_TYPE == "Drought", "YEAR"] |>
  table()
csv_details[grepl("(S|s)torm", EVENT_TYPE), ][STATE == "NORTH CAROLINA", c("YEAR", "MAGNITUDE")] %>%
  table()
csv_details[grepl("(T|t)orna", EVENT_TYPE), ] %>%
  .[STATE == "NORTH CAROLINA", c("YEAR", "TOR_F_SCALE")] %>%
  table()
csv_details[grepl("(F|f)ire", EVENT_TYPE), "YEAR"] |>
  table()

csv_details[!is.na(TOR_WIDTH) & TOR_WIDTH > 0, ] %>%
  terra::vect(geom = c("END_LON", "END_LAT"), crs = "EPSG:4326", keepgeom = TRUE) %>%
  terra::buffer(width = .$TOR_WIDTH) %>%
  plot(xlim = c(-84, -80), ylim = c(33, 37))

# Tornados
# "blocked" pattern; still a few in NC
csv_details[!is.na(TOR_WIDTH) & TOR_WIDTH > 0,
  c("BEGIN_LON", "BEGIN_LAT", "END_LON", "END_LAT")] %>%
  as.matrix() %>%
  terra::vect(crs = "EPSG:4326", type = "lines") %>%
  plot()
  # plot(xlim = c(-90, -80), ylim = c(32, 38))


## Intersection of noaa and fema
years <- seq(2000, 2023)
map_noaa <- csv_details[, c("BEGIN_LON", "BEGIN_LAT", "YEAR", "EVENT_TYPE")] %>%
  terra::vect(geom = c("BEGIN_LON", "BEGIN_LAT"), crs = "EPSG:4326") %>%
  na.omit(geom = TRUE)
par(mfcol = c(5, 6))
for (i in years) {
  mapsf::mf_map(
    sf::st_as_sf(map_noaa[map_noaa$YEAR == i & grepl("(Storm|storm|STORM)", map_noaa$EVENT_TYPE), ]),
    type = "base",
    cex = 0.01,
    xlim = c(-128, -64), ylim = c(22, 52))
}


map_fema <- data.table::fread("input/fema/IndividualAssistanceHousingRegistrantsLargeDisasters.csv")
map_fema <- data.table::fread("input/fema/HousingAssistanceRenters.csv")
map_fema <- data.table::fread("input/fema/DisasterDeclarationsSummaries.csv")
map_fema[grepl("(S|s)torm", incidentType), c("fyDeclared")] %>%
  table
map_fema[, c("fyDeclared")] %>%
  table
