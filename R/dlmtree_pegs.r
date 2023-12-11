.libPaths("~/r-libs")
source("./R/load_packages.r")

if (!require(pak)) {
  install.packages("pak")
  library(pak)
}
if (!require(dlmtree)) {
  pak::pak("danielmork/dlmtree")
  library(dlmtree)
}

source("./R/prep_data.r")



