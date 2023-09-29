library(pacman)
p_load(sf, terra, dplyr, purrr, data.table, vtable)
sf_use_s2(F)


## data path determination
## compute_mode 1 (wine mount), compute_mode 2 (hpc compute node), 3 (container internal)
COMPUTE_MODE = 1
path_pegs = 
    ifelse(COMPUTE_MODE == 1,
    "/Volumes/PEGS/",
    ifelse(COMPUTE_MODE == 2,
    "/ddn/gs1/project/controlled/PEGS/",
    ifelse(COMPUTE_MODE == 3,
    "/opt/", stop("COMPUTE_MODE should be one of 1, 2, or 3.\n"))))


list.files(paste0(path_pegs, "Data_Freezes/latest"), 
        pattern = "*.(csv|RData|rdata|rds|RDS|CSV|txt)$",
        recursive = TRUE, full.names = TRUE)


rdatas = list.files(paste0(path_pegs, "Data_Freezes/latest"), 
        pattern = "*.(RData)$",
        recursive = TRUE, full.names = TRUE)

for (rdata in rdatas) {
    load(rdata)
}

obj.list = sapply(rdatas, load)
obj.names = obj.list %>%
    t %>%
    data.frame(fpath = rownames(.), .) %>%
    tidyr::pivot_longer(cols = 2:3)
obj.sizes = obj.list %>%
    lapply(function(x) Reduce(rbind, sapply(x, function(y) dim(get(y)), simplify = F))) %>%
    Reduce(rbind, .) %>%
    as.data.frame
cbind(obj.names, obj.sizes) %>%
    write.csv("./output/object_list.csv", row.names = F)

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
epr.gis = epr.gis %>%
    mutate(gis_event_date = as.Date(gis_event_date, format = "%M/%d/%Y"))


is.date <- function(x) inherits(x, 'Date')

## summary tables
lapply(obj.names$value[!grepl("meta", obj.names$value)],
    function(x) {
        get(x) %>%
            mutate_if(is.character, ~as.factor(.)) %>%
            mutate_if(is.date, ~as.numeric(.)) %>%
            vtable::sumtable(vars = names(.), 
                            add.median = TRUE,
                            #out = "csv",
                            file = sprintf("./output/summary_%s.html", x))
    })
epr.svi %>%
    vtable::sumtable(vars = names(.), 
                     add.median = TRUE,
                     file = "./output/summary_svi.html")

summary(epr.he)


## Location maps
p_load(ggplot2, tigris)
us_cnty = tigris::counties(year = 2020) %>%
    filter(!STATEFP %in% c("02", "72", "78", "15", "60", "69", "66")) %>%
    st_transform(4326)
epr_sf = epr.gis %>%
    st_as_sf(coords = c("gis_longitude", "gis_latitude"), crs = 4326)
epr.gis$gis_state %>% table

epr_sfs = epr_sf[us_cnty,] %>%
    st_jitter(0.3)
ggindloc = ggplot() +
    geom_sf(data = us_cnty, linewidth = 0.8) +
    geom_sf(data = epr_sfs, pch = 19, cex = 0.6, col = 'red') +
    labs(caption = "Coordinates were jittered.")
ggsave("./output/distr.png", ggindloc, width = 15, height = 9, units = "in", dpi = 300)






## Detailed metadata tables -- failed
## Metadata tables are seemingly generated from SAS
codebooks = list.files(paste0(path_pegs, "Data_Freezes/latest"),
            pattern = "*_codebook_([0-9]+{1,2})*.*.pdf$",
            recursive = TRUE,
            full.names = TRUE)
codebooks = codebooks[!grepl("appendix", codebooks)]

codebook_tables = lapply(codebooks[2], 
    function(x) tabulizer::extract_tables(x, method = "lattice", output = "data.frame") %>% map(as_tibble))

pdf_tables <- pdf_data |> 
  extract_tables(output = "data.frame") |>
  map(as_tibble)