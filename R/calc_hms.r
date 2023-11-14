pkgs <- c("sf", "terra", "dplyr")
invisible(sapply(pkgs, library, character.only = TRUE, quietly = TRUE))
options(sf_use_s2 = FALSE)

args <- commandArgs(trailingOnly = TRUE)

# Mitchell's hms smoke data processing code
# patched
source("https://github.com/Spatiotemporal-Exposures-and-Toxicology/NRTAPmodel/raw/mm_noaa_smoke_patch_1103/input/Rinput/download_functions/download_noaa_hms_smoke_data.R")

# source("https://github.com/Spatiotemporal-Exposures-and-Toxicology/NRTAPmodel/raw/main/input/Rinput/download_functions/download_noaa_hms_smoke_data.R")
# source("https://github.com/Spatiotemporal-Exposures-and-Toxicology/NRTAPmodel/raw/main/input/Rinput/process_noaa_hms_smoke_data.R")


# process_noaa_hms_smoke_data
# Data available from 2005-08-06
# Smoke density available from 2008-04-11

# as.Date("2020-12-31") - as.Date("2008-04-11")
startdate <- args[1]
enddate <- args[2]

print(startdate)
print(enddate)


targdir <- "./hmsdata/"
if (!dir.exists(targdir)) {
    dir.create(targdir)
}

# "https://satepsanone.nesdis.noaa.gov/pub/FIRE/web/HMS/Smoke_Polygons/Shapefile/2008/04/hms_smoke20080411.zip"
startdate <- "2022-07-01"
enddate <- "2022-12-31"


download_noaa_hms_smoke_data(
    date_start = startdate,
    date_end = enddate,
    directory_to_download = "/Users/songi2/Documents/input/noaa_hms/raw/",
    directory_to_save = "/Users/songi2/Documents/input/noaa_hms/shapefile/",
    data_download_acknowledgement = TRUE,
    remove_download = FALSE,
    time_wait_download = 1L
)


# Downloading requested files...
# Requested files downloaded.
# Unzipping shapefiles to /Users/songi2/Documents/input/noaa_hms/shapefile/...
# Files unzipped and saved in/Users/songi2/Documents/input/noaa_hms/shapefile/.
# There were 50 or more warnings (use warnings() to see the first 50)
# > warnings()
# Warning messages:
# 1: In download.file(file_urls, download_names, method = "libcurl",  ... :
#   cannot open URL 'https://satepsanone.nesdis.noaa.gov/pub/FIRE/web/HMS/Smoke_Polygons/Shapefile/2008/12/hms_smoke20080821.zip': HTTP status was '404 Not Found'
# 2: In download.file(file_urls, download_names, method = "libcurl",  ... :
#   cannot open URL 'https://satepsanone.nesdis.noaa.gov/pub/FIRE/web/HMS/Smoke_Polygons/Shapefile/2008/12/hms_smoke20080902.zip': HTTP status was '404 Not Found'
# 3: In download.file(file_urls, download_names, method = "libcurl",  ... :
#   cannot open URL 'https://satepsanone.nesdis.noaa.gov/pub/FIRE/web/HMS/Smoke_Polygons/Shapefile/2008/12/hms_smoke20080501.zip': HTTP status was '404 Not Found'
# 4: In download.file(file_urls, download_names, method = "libcurl",  ... :
#   cannot open URL 'https://satepsanone.nesdis.noaa.gov/pub/FIRE/web/HMS/Smoke_Polygons/Shapefile/2008/12/hms_smoke20080502.zip': HTTP status was '404 Not Found'
# 5: In download.file(file_urls, download_names, method = "libcurl",  ... :
#   cannot open URL 'https://satepsanone.nesdis.noaa.gov/pub/FIRE/web/HMS/Smoke_Polygons/Shapefile/2008/12/hms_smoke20080505.zip': HTTP status was '404 Not Found'
# 6: In download.file(file_urls, download_names, method = "libcurl",  ... :
#   cannot open URL 'https://satepsanone.nesdis.noaa.gov/pub/FIRE/web/HMS/Smoke_Polygons/Shapefile/2008/12/hms_smoke20080506.zip': HTTP status was '404 Not Found'
# 7: In download.file(file_urls, download_names, method = "libcurl",  ... :
#   cannot open URL 'https://satepsanone.nesdis.noaa.gov/pub/FIRE/web/HMS/Smoke_Polygons/Shapefile/2008/12/hms_smoke20080507.zip': HTTP status was '404 Not Found'
# [truncated to save space ...]