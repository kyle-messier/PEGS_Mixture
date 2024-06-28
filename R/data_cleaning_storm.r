
source("R/load_storm_data.r")
source("R/prep_add_snps.r")
# this file will return
# csv_details_sub2: recoded event types
# csv_details2: minimally cleaned details
# map_fema: FEMA declared disaster data

collapse::set_collapse(mask = "manip")

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

# qs::qsave(csv_details, "input/noaastorm/storm_details_all.qs")
# qs::qsave(csv_fatal, "input/noaastorm/storm_fatal_all.qs")
# qs::qsave(csv_locs, "input/noaastorm/storm_locs_all.qs")

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

# tornado width around NC
csv_details2[!is.na(TOR_WIDTH) & TOR_WIDTH > 0, ] %>%
  terra::vect(geom = c("END_LON", "END_LAT"), crs = "EPSG:4326", keepgeom = TRUE) %>%
  terra::buffer(width = .$TOR_WIDTH) %>%
  plot(xlim = c(-84, -80), ylim = c(33, 37))

# Tornados
# "blocked" pattern; still a few in NC
csv_details2[!is.na(TOR_WIDTH) & TOR_WIDTH > 0,
  c("BEGIN_LON", "BEGIN_LAT", "END_LON", "END_LAT")] %>%
  as.matrix() %>%
  .[!is.na(.[,3]), ] %>%
  terra::vect(crs = "EPSG:4326", type = "lines") %>%
  plot()
  #plot(xlim = c(-90, -80), ylim = c(32, 38))


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





## nc simulation
ncpath <- system.file("gpkg/nc.gpkg", package = "sf")
nc <- sf::st_read(ncpath)
nc <- sf::st_transform(nc, "EPSG:5070")
ncp <- sf::st_sample(nc, 5000, type = "Thomas", kappa = 8e-10, mu = 50, scale = 20000)
ncp$pid <- sprintf("ID-%04d", seq_len(nrow(ncp)))
nrow(ncp)
plot(ncp)

unique(csv_details$EVENT_TYPE)



## Exploratory plots ####
years <- seq(2006, 2023)
map_noaa_l <- copy(csv_details_sub2) %>%
  tidytable::filter((!is.na(BEGIN_LON) & !is.na(END_LON)) & !(BEGIN_LON == END_LON & BEGIN_LAT == END_LAT)) %>%
  tidytable::select(tidytable::all_of(c("EVENT_ID", "BEGIN_LON", "BEGIN_LAT", "END_LON", "END_LAT", "BEGIN_RANGE", "END_RANGE", "YEAR", "event_re"))) %>%
  tidytable::filter(!is.na(BEGIN_RANGE) & !is.na(END_RANGE)) %>%
  tidytable::filter((BEGIN_RANGE > 0) & (END_RANGE > 0)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(geometry = list(sf::st_linestring(matrix(c(BEGIN_LON, BEGIN_LAT, END_LON, END_LAT), ncol = 2, byrow = TRUE)))) %>%
  dplyr::ungroup() %>%
  sf::st_as_sf(sf_column_name = "geometry") %>%
  sf::st_set_crs("EPSG:4326") %>%
  dplyr::mutate(geom_type = "line")

noaa_lm <-
  Map(function(x1, y1, x2, y2, br, er) {
  conv_point_trajectory(
    c(x1, y1),
    c(x2, y2),
    wd_start = br * 1609,
    wd_end = er * 1609,
  )
  }, map_noaa_l$BEGIN_LON, map_noaa_l$BEGIN_LAT, map_noaa_l$END_LON, map_noaa_l$END_LAT, map_noaa_l$BEGIN_RANGE, map_noaa_l$END_RANGE)

  dplyr::mutate(
    geometry =
    list(sf::st_as_sfc(
      terra::geom(
        conv_point_trajectory(
          c(BEGIN_LON, BEGIN_LAT),
          c(END_LON, END_LAT),
          wd_start = BEGIN_RANGE * 1609,
          wd_end = END_RANGE * 1609,
        ),
        wkt = TRUE
      )
    )
    )
  ) %>%
  dplyr::ungroup() %>%
  sf::st_as_sf(sf_column_name = "geometry") %>%
  sf::st_set_crs("EPSG:4326") %>%
  dplyr::mutate(geom_type = "line")
map_noaa_p <- copy(csv_details_sub2) %>%
  tidytable::filter((!is.na(BEGIN_LON) & is.na(END_LON)) | (BEGIN_LON == END_LON & BEGIN_LAT == END_LAT)) %>%
  tidytable::select(tidytable::all_of(c("EVENT_ID", "BEGIN_LON", "BEGIN_LAT", "BEGIN_RANGE", "YEAR", "event_re"))) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(geometry = list(sf::st_point(matrix(c(BEGIN_LON, BEGIN_LAT), ncol = 2, byrow = TRUE)))) %>%
  dplyr::ungroup() %>%
  sf::st_as_sf(sf_column_name = "geometry") %>%
  sf::st_set_crs("EPSG:4326") %>%
  dplyr::mutate(geom_type = "point")


map_noaa <-
  dplyr::bind_rows(map_noaa_l, map_noaa_p)


## point + line
line_sub <- copy(csv_details_sub)[(!is.na(BEGIN_LON) & !is.na(END_LON)) & !(BEGIN_LON == END_LON & BEGIN_LAT == END_LAT), ]
line_sub_v <- line_sub[, c("BEGIN_LON", "BEGIN_LAT", "END_LON", "END_LAT")] %>%
  as.matrix() %>%
  terra::vect(crs = "EPSG:4326", type = "lines")
line_sub_v <- cbind(line_sub_v, line_sub)

point_sub2 <- copy(csv_details_sub)[(BEGIN_LON == END_LON & BEGIN_LAT == END_LAT), ]
point_sub <- copy(csv_details_sub)[!is.na(BEGIN_LON) & is.na(END_LON), ]
point_sub_v <- point_sub[, c("BEGIN_LON", "BEGIN_LAT")] %>%
  as.matrix() %>%
  terra::vect(crs = "EPSG:4326", type = "points")
point_sub_v <- cbind(point_sub_v, point_sub)



# Assuming you have the starting point `start_point` and the end point `end_point` as terra spatvector objects,
# and you have set the sampling interval `interval`:
# Assuming you have the coordinates of the starting point and the end point


#' Convert start end end points to funnel-shaped
#' trajectory polygon
#' 
#' Extreme weather observation points are often
#' containing the widths at the start and end points.
#' This function creates a funnel-shaped trajectory
#' polygon based on the information with interval.
#' The direct line between the start and end points
#' is densified and buffered with increasing radius
#' at each interval. Please note that radius is linearly
#' increased between the start and end points.
#' 
#' @param start_point A numeric vector of length 2 representing the
#'  coordinates of the starting point.
#' @param end_point A numeric vector of length 2 representing the
#' coordinates of the ending point.
#' @param interval A numeric value representing the sampling
#' interval between the start and end points. It follows
#' the linear unit of the CRS.
#' @param wd_start A numeric value representing the width at the
#' starting point.
#' @param wd_end A numeric value representing the width at the
#' ending point.
#' @param fun A function to be applied to the radius.
#' This function will be applied to the serial index of the
#' radius. The default is NULL.
#' @returns A SpatVector object representing the funnel-shaped
#' trajectory polygon.
#' @export
#' @examples
#' start_point <- c(-80, 35)
#' end_point <- c(-78, 37)
#' interval <- 1e2
#'
#' conv_point_trajectory(start_point, end_point, interval, 1, 10)
#' conv_point_trajectory(start_point, end_point, interval, 1, 10, function(x) -0.5 * (x-4)^2))
#'
conv_point_trajectory <-
  function(
    start_point = numeric(2),
    end_point = numeric(2),
    interval = 1e2,
    wd_start = numeric(1),
    wd_end = numeric(1),
    fun = NULL
  ) {

    start_point <-
      terra::vect(
        matrix(start_point, nrow = 1),
        crs = "EPSG:4326",
        type = "points"
      )
    end_point <-
      terra::vect(
        matrix(end_point, nrow = 1),
        crs = "EPSG:4326",
        type = "points"
      )

    # Create a spatvector object with the two points
    spatvector <- terra::vect(c(start_point, end_point))
    line <- terra::as.lines(spatvector)
    line$id <- 1

    linedens <- terra::densify(line, interval)
    linedensp <- terra::as.points(linedens)
    # return(linedensp)
    radius <- seq(wd_start, wd_end, length.out = nrow(linedensp))
    if (!is.null(fun)) {
      radius <- fun(seq_along(radius))
    }
    linedenspb <-
      terra::buffer(
        linedensp,
        radius,
        quadsegs = 90L
      )
    linedenspbm <- terra::aggregate(linedenspb)
    return(linedenspbm)
  }

start_point <- c(-80, 35)
end_point <- c(-78, 37)
interval <- 1e2

t1 <- conv_point_trajectory(start_point, end_point, 100, 100, 50000)
t2 <-
  conv_point_trajectory(
    start_point, end_point,
    interval,
    100, 10000, function(x) 10000* log10(x+1))
plot(t1)

t2 <-
  conv_point_trajectory(
    start_point, end_point,
    interval,
    1000, 10000)
plot(t2)

library(data.table)
library(collapse)


map_noaa_l



# par(mfcol = c(1, 1))
# map_noaa_l %>%
#   dplyr::filter(event_re == "Tornado") %>%
#   mapsf::mf_map(lwd = 0.8, xlim = c(-128, -64), ylim = c(22, 50))
# map_noaa_l %>%
#   dplyr::filter(event_re == "Thunderstorm") %>%
#   mapsf::mf_map(lwd = 0.8, xlim = c(-128, -64), ylim = c(22, 50))


# ##
# # png("output/noaa_tornado.png", width = 12, height = 10, units= "in", res = 508)
# par(mfcol = c(5, 4))
# for (i in years) {
#   mapsf::mf_map(
#     map_noaa_p %>% dplyr::filter(YEAR == i & event_re == "Tornado"),
#     type = "base",
#     cex = 0.2,
#     col = "red",
#     pch = 19,
#     xlim = c(-128, -64), ylim = c(22, 52)
#   )
#   mapsf::mf_map(
#     map_noaa_l %>% dplyr::filter(YEAR == i & event_re == "Tornado"),
#     add = TRUE,
#     lwd = 0.8)
# }
# # dev.off()

# # png("output/noaa_thunderstorm.png", width = 12, height = 10, units= "in", res = 508)
# par(mfcol = c(5, 4))
# for (i in years) {
#   mapsf::mf_map(
#     map_noaa_p %>% dplyr::filter(YEAR == i & event_re == "Thunderstorm"),
#     type = "base",
#     cex = 0.2,
#     col = "red",
#     pch = 19,
#     xlim = c(-128, -64), ylim = c(22, 52)
#   )
#   mapsf::mf_map(
#     map_noaa_l %>% dplyr::filter(YEAR == i & event_re == "Thunderstorm"),
#     add = TRUE,
#     lwd = 0.8)
# }
# # dev.off()

# sts <- tigris::states(year = 2020, cb = TRUE) %>%
#   sf::st_transform("EPSG:4326")

# plot(
#   map_noaa_l %>% dplyr::filter(event_re == "Storm") %>% .$geometry,
#     # st_intersection(., st_as_sfc(st_bbox(c(xmin = -84, xmax = -80, ymin = 33, ymax = 37), crs = 4326))),
#   #type = "base",
#   xlim = c(-84, -76),
#   ylim = c(34, 37.5),
#   lwd = 0.8
# )
# plot(sts$geometry, add = TRUE, border = "purple", lwd = 1.5)

# plot(
#   map_noaa_l %>% dplyr::filter(event_re == "Tornado") %>% .$geometry,
#     # st_intersection(., st_as_sfc(st_bbox(c(xmin = -84, xmax = -80, ymin = 33, ymax = 37), crs = 4326))),
#   #type = "base",
#   xlim = c(-84, -76),
#   ylim = c(34, 37.5),
#   lwd = 0.8
# )
# plot(sts$geometry, add = TRUE, border = "purple", lwd = 1.5)

# plot(
#   map_noaa %>% dplyr::filter(event_re == "Storm") %>% .$geometry,
#   xlim = c(-84, -76),
#   ylim = c(34, 37.5),
#   cex = 0.6,
#   pch = 19,
#   lwd = 0.8
# )
# plot(sts$geometry, add = TRUE, border = "purple", lwd = 1.5)


## count
# ncp <- sf::st_set_crs(ncp, "EPSG:5070")
source("R/prep_add_snps.r")
pegs_main_psssf <-
  pegs_main_pss %>%
  as.data.frame %>%
  sf::st_as_sf(sf_column_name = "geometry", crs = 5070)

map_noaa_lt <- sf::st_transform(map_noaa_l, "EPSG:5070") %>%
  dplyr::mutate(
    endrange = ifelse(is.na(END_RANGE) | END_RANGE == 0,
    min(END_RANGE[which(END_RANGE != 0)], na.rm = TRUE),
    END_RANGE)
  ) %>%
  dplyr::ungroup()
map_noaa_pt <- sf::st_transform(map_noaa_p, "EPSG:5070") %>%
  dplyr::mutate(
    endrange = ifelse(is.na(BEGIN_RANGE) | BEGIN_RANGE == 0,
    min(BEGIN_RANGE[which(BEGIN_RANGE != 0)], na.rm = TRUE),
    BEGIN_RANGE)
  ) %>%
  dplyr::ungroup()


map_noaa_ltb <- map_noaa_lt %>%
  sf::st_buffer(units::set_units(.$endrange, "mi"))
map_noaa_ptb <- map_noaa_pt %>%
  sf::st_buffer(units::set_units(.$endrange, "mi"))

map_noaa_lptb <- dplyr::bind_rows(map_noaa_ltb, map_noaa_ptb)

nc <- sf::st_read(system.file("gpkg/nc.gpkg", package = "sf"))
nc <- sf::st_transform(nc, "EPSG:5070")
map_noaa_lptbnc <- map_noaa_lptb[nc, ]
plot(map_noaa_lptbnc$geometry)

map_noaa_lptbnc <- qs::qread("output/noaa_nc_polygons.qs")
vct_noaa_nc <- terra::vect(map_noaa_lptbnc)
vct_noaa_nc$occurrence <- 1
rst_nc_template <- terra::rast(terra::ext(vct_noaa_nc), resolution = 100)
rst_noaa_nc <- terra::rasterize(vct_noaa_nc, rst_nc_template, field = "occurrence", fun = "sum")


pegs_nc_occ <- sf::st_buffer(pegs_main_psssf, units::set_units(10, "km"))
pegs_nc_occ_count <- exactextractr::exact_extract(
  rst_noaa_nc, pegs_nc_occ, fun = "mean"
)

pegs_main_psssfx <- pegs_main_psssf %>%
  dplyr::mutate(
    n_ewe = pegs_nc_occ_count
  ) %>%
  dplyr::filter(!is.nan(n_ewe)) %>%
  dplyr::mutate(age_2020 = 2021 - as.POSIXlt(birth_date)$year)
dim(pegs_main_psssfx)

explanatory <- c(687, 688, 691:703, 706, 707)
# 494:533, 
mainform <- reformulate(
  termlabels = names(pegs_main_psssfx)[explanatory],
  response = "pss_class"
)

pegs_main_pss <-
  pegs_main_pss |>
  select(-geometry)
pegs_main_pss_com <- pegs_main_pss[complete.cases(pegs_main_pss), ]

mainmult <- MASS::polr(mainform, data = pegs_main_psssfx, Hess = TRUE)


## previous way: set circular buffer from residences and
##   count intersecting features
ncpb <- sf::st_buffer(pegs_main_psssf, units::set_units(10, "km")) %>%
  dplyr::mutate(n_l_ewe = lengths(sf::st_intersects(geometry, map_noaa_lt))) %>%
  dplyr::mutate(n_p_ewe = lengths(sf::st_intersects(geometry, map_noaa_pt)))

years <- seq(2001, 2020)
lyears <- split(years, years) %>%
  lapply(function(tg) {
    na0 <- function(x) ifelse(is.na(x), 0, x)
    nck <-
      sf::st_buffer(pegs_main_psssf, units::set_units(20, "km")) %>%
      dplyr::mutate(
        n_l_ewe = lengths(sf::st_intersects(geometry, map_noaa_lt %>% dplyr::filter(YEAR == tg)))
      ) %>%
      dplyr::mutate(
        n_p_ewe = lengths(sf::st_intersects(geometry, map_noaa_pt %>% dplyr::filter(YEAR == tg)))
      ) %>%
      dplyr::mutate(n_ewe = na0(n_l_ewe) + na0(n_p_ewe))
    names(nck)[ncol(nck)] <- sprintf("n_ewe_%d", tg)
    return(nck)
  })


ncbox <- st_bbox(c(xmin = -84.5, xmax = -76, ymin = 34.5, ymax = 37.5), crs = 4326)
ncbox <- st_as_sfc(ncbox)
ncbox <- st_transform(ncbox, 5070)

png("output/yearly_exposure.png", width = 12, height = 10, units= "in", res = 508)
par(mfcol = c(4, 5))
for (y in years) {
  mapsf::mf_init(x = ncbox)
  mapsf::mf_choro(sf::st_centroid(lyears[[y - 2000]]), sprintf("n_ewe_%d", y),
  breaks = c(0, 25, 50, 75, 100, 150, 200), cex = 0.5, pch = 19,
  add = TRUE)
}
dev.off()

par(mfrow = c(1, 1))
mapsf::mf_init(x = ncbox)
mapsf::mf_choro(lyears[[1]], "n_ewe_2001")

# TODO: by intensity