source("./R/load_packages.r")

dpath <- "~/Documents/input/noaa_hms"
zips <- list.files(path = paste0(dpath, "/raw"), pattern = "*.zip$", full.names = TRUE)
targets <- paste0(dpath, "/shapefile")
# invisible(sapply(zips, unzip, exdir = targets))

filedates <- gsub("*.*hms_smoke_Shapefile_", "", zips)
filedates <- gsub("\\.(zip)", "", filedates)
filedates <- as.Date(filedates, format = "%Y%m%d")

shps <- list.files(path = targets, pattern = "*.shp$", full.names = TRUE)
# version: by density
# version: by total extent (union)
# single worker to process single file
worker_density <- function(path, time = NULL,
  sub_dens = "Heavy", template = baserast) {
  # parse time
  smoke <- terra::vect(path)
  smoke <- terra::project(smoke, "EPSG:5070")
  # smoke$Density <- factor(smoke$Density, levels = c("Heavy", "Medium", "Light"))
  # smoke <- smoke[order(smoke$Density, smoke$End),]
  formattext <- "%Y%j %H%M"
  smoke$time_start <- as.POSIXct(smoke$Start, format = formattext)
  smoke$time_end <- as.POSIXct(smoke$End, format = formattext)
  smoke$timediff_h <- (as.numeric(smoke$time_end) - as.numeric(smoke$time_start)) / 3600
  smoke <- smoke[smoke$Density == sub_dens, ]
  smoke_r <- terra::rasterize(smoke, template, field = "timediff_h")
  return(smoke_r)
}


worker_all <- function(path, time, template = baserast) {
  # parse time
  smoke <- terra::vect(path)
  smoke <- terra::project(smoke, "EPSG:5070")
  smoke <- terra::aggregate(smoke)
  smoke_r <- terra::rasterize(smoke, template)
  return(smoke_r)
}


# blank raster
target_ext <-
  terra::vect(terra::ext(c(-126, -70, 25, 50)),
    crs = "EPSG:4326") |>
  terra::project("EPSG:5070")


baserast <- terra::rast(
  extent = target_ext,
  resolution = 4e3L,
  crs = "EPSG:5070"
)

targrasts <- gsub("raw/hms_smoke_Shapefile", "processed_heavy/smoke_heavy", zips)
targrasts <- gsub("zip", "nc", targrasts)

plan(multicore, workers = 8)
doFuture::registerDoFuture()
foreach(
  x = seq_along(shps),
  .export = c("baserast", "zips", "shps", "targrasts", "worker_density"),
  .packages = c("dplyr", "terra")
) %dopar% {
  pcs <- worker_density(shps[x], template = baserast)
  terra::writeRaster(pcs, filename = targrasts[x])
}

targrasts <- gsub("raw/hms_smoke_Shapefile", "processed_medium/smoke_medium", zips)
targrasts <- gsub("zip", "nc", targrasts)

foreach(
  x = seq_along(shps),
  .export = c("baserast", "zips", "shps", "targrasts", "worker_density"),
  .packages = c("dplyr", "terra")
) %dopar% {
  pcs <- worker_density(shps[x], template = baserast, sub_dens = "Medium")
  terra::writeRaster(pcs, filename = targrasts[x])
}

targrasts <- gsub("raw/hms_smoke_Shapefile", "processed_light/smoke_light", zips)
targrasts <- gsub("zip", "nc", targrasts)

foreach(
  x = seq_along(shps),
  .export = c("baserast", "zips", "shps", "targrasts", "worker_density"),
  .packages = c("dplyr", "terra")
) %dopar% {
  pcs <- worker_density(shps[x], template = baserast, sub_dens = "Light")
  terra::writeRaster(pcs, filename = targrasts[x])
}


# testcode
test <- worker_density(shps[1])
plot(test)

test1 <- lapply(shps[1:100], worker_density, sub_dens = "Light")
test1m <- do.call(c, test1)
plot(test1m[,,80:95])


ff <- terra::vect(shps[124])
tf <- "%Y%j %H%M"
indx <- 7
dd <- as.POSIXct(ff$Start[indx], format = tf)
de <- as.POSIXct(ff$End[indx], format = tf)
unclass(de-dd)
diff(de, dd)
