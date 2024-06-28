## Load packages ####
pkgs <- c("data.table", "sf", "stars", "sftime", "terra", "dplyr", "tidyr", "collapse", "qs")
invisible(sapply(pkgs, library, quietly = TRUE, character.only = TRUE))
options(sf_use_s2 = FALSE)

#
# library(rvest)
# turl <- "https://www.ncei.noaa.gov/pub/data/swdi/stormevents/csvfiles/"
# thtml <- read_html(turl)
# tlinks <- thtml |>
#   html_table()
# tlinks <- tlinks[[1]]$Name
# tlinks <- grep("StormEvents", tlinks, value = TRUE)
# tlinkspre2000 <- grep("(d20*)", tlinks, value = TRUE, invert = TRUE)
# tlinkspre2000

# tlinkspre2000e <- sprintf("%s%s", turl, tlinkspre2000)

# for (i in tlinkspre2000) {
#   download.file(
#     sprintf("%s%s", turl, i),
#     sprintf("%s%s", "./input/noaastorm/", i),
#     mode = "wb"
#   )
#   Sys.sleep(0.5)
# }

# load file
csv_details <- qs::qread("output/noaa_stormdb.qs")

# unique(csv_details$EVENT_TYPE)
# csv_details[EVENT_TYPE == "Heat", ]
# csv_details[EVENT_TYPE == "Flood", ]
# unique(csv_details$DAMAGE_PROPERTY)

# property damage cleaning
# TODO: check if as-then or standardized dollar values
csv_details2 <- csv_details[,
`:=`(damage_property =
       substr(DAMAGE_PROPERTY, 1, nchar(DAMAGE_PROPERTY) - 1),
     damage_property_unit =
       substr(DAMAGE_PROPERTY, nchar(DAMAGE_PROPERTY), nchar(DAMAGE_PROPERTY)))][,
`:=`(damage_property =
       as.numeric(damage_property), 
     damage_property_unit =
       plyr::mapvalues(damage_property_unit, c("h", "H", "K", "M", "B"), c(1e3, 1e3, 1e3, 1e6, 1e9)))]



### NOAA event types reclassification ####
trans_hazard <-
  dplyr::tribble(
    ~type,  ~group,
    "Tornado",  "Tornado",
    "TORNADO/WATERSPOUT", "Tornado",
    "TORNADOES, TSTM WIND, HAIL", "Tornado",
    "Tropical Storm", "Storm",
    "Hurricane (Typhoon)",  "Storm",
    "Tropical Depression",  "Storm",
    "Marine Tropical Storm",  "Storm",
    "Marine Hurricane/Typhoon", "Storm",
    "Marine Tropical Depression", "Storm",
    "Hurricane",  "Storm",
    "Thunderstorm Wind",  "Thunderstorm",
    "THUNDERSTORM WINDS/FLOODING", "Flood",
    "THUNDERSTORM WINDS/FLASH FLOOD", "Flood",
    "THUNDERSTORM WIND/ TREE", "Thunderstorm",
    "THUNDERSTORM WIND/ TREES", "Thunderstorm",
    "THUNDERSTORM WINDS FUNNEL CLOU", "Thunderstorm",
    "THUNDERSTORM WINDS/HEAVY RAIN", "Thunderstorm",
    "THUNDERSTORM WINDS/ FLOOD", "Flood",
    "THUNDERSTORM WINDS HEAVY RAIN", "Thunderstorm",
    "Flash Flood",  "Flood",
    "Winter Storm", "Storm",
    "Ice Storm",  "Storm",
    "Flood",  "Flood",
    "Coastal Flood",  "Flood",
    "Wildfire",  "Wildfire",
    "Drought",  "Drought",
    "Dust Storm", "Storm",
    "Storm Surge/Tide", "Storm",
    "Extreme Cold/Wind Chill",  "Cold",
    "Excessive Heat", "Heat",
    "Marine Thunderstorm Wind", "Thunderstorm",
    "Dense Smoke",  "Wildfire",
    "Lakeshore Flood",  "Flood",
    "Frost/Freeze", "Cold",
    "Cold/Wind Chill",  "Cold",
    "HAIL FLOODING",  "Flood",
    "Blizzard", "Cold",
    "Heavy Snow", "Cold"
  )

csv_details_sub <-
  csv_details[EVENT_TYPE %in% trans_hazard$type, ]
csv_details_sub <- csv_details_sub[,
    `:=`(event_re = plyr::mapvalues(EVENT_TYPE, trans_hazard$type, trans_hazard$group))]
csv_details_sub2 <- csv_details_sub[!is.na(BEGIN_LON), ][
  (BEGIN_LON <= -64 & BEGIN_LON >= -128) & (BEGIN_LAT >= 20 & BEGIN_LAT <= 52),
][YEAR >= 2000,]


map_fema <- data.table::fread("input/fema/IndividualAssistanceHousingRegistrantsLargeDisasters.csv")
map_fema <- data.table::fread("input/fema/HousingAssistanceRenters.csv")
map_fema <- data.table::fread("input/fema/DisasterDeclarationsSummaries.csv")
map_fema[grepl("(S|s)torm", incidentType), c("fyDeclared")] %>%
  table
map_fema[, c("fyDeclared")] %>%
  table
