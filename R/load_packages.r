pkgs <- c("qs", "dplyr", "data.table", "tidytable", "terra", "here",
    "lubridate", "future", "future.apply", "foreach",
    "skimr", "tigris", "vtable", "sf", "ggplot2")
suppressMessages(invisible(sapply(pkgs, library, quietly = TRUE, character.only = TRUE)))
options(sf_use_s2 = FALSE)
