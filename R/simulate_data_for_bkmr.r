### Read packages ###
pkgs = c('dplyr', 'data.table', 'sf', 'terra', 'skimr', 'rgenoud', 'MatchIt', 'simstudy', 'fields')
invisible(lapply(pkgs, library, character.only = TRUE, quietly = TRUE, logical.return = TRUE))

Sys.setenv(sf_use_s2 = FALSE)
set.seed(202310)

## Part 1 : generate raster datasets with fields
gr_list <- list(x = seq(1, 100, 1), y = seq(1, 100, 1))
gr_exgr <- gr_list

# range and nu values
smth_rng <- expand.grid(
  range = seq(5, 30, 5),
  nu = c(seq(0.25, 2, 0.25), 2.5)
)
n_smth_rng <- nrow(smth_rng)
seq_smth_rng <- seq(1, n_smth_rng)


gr_mat_sim <- matern.image.cov(grid = gr_list, setup = TRUE, aRange = 20, smoothness = 0.6)
base_grid <- expand.grid(
  x = gr_mat_sim$grid$x,
  y = gr_mat_sim$grid$y
  ) |>
  as.matrix()

look <- sim.rf(gr_mat_sim)
# circulantEmbeddingSetup(matrix(nrow = 100, ncol = 100), cov.function = "Matern")
# set.panel(1,1)
image.plot(gr_list$x, gr_list$y, look)
lk_vgram <- fields::vgram(base_grid, as.vector(t(look)), d = seq(0, 0.5, 0.05), type = 'variogram')
plot(lk_vgram)


grid <- list(x = seq(0, 5, , 100), y = seq(0, 5, , 100))
obj <- fields::Exp.image.cov(grid = grid, theta = .5, setup = TRUE)
look <- fields::sim.rf(obj)
# Now simulate another ...
look2 <- fields::sim.rf(obj)
# take a look
set.panel(2, 1)
image.plot(grid$x, grid$y, look)
title("simulated gaussian field")
image.plot(grid$x, grid$y, look2)
title("another (independent) realization ...")

##
lookv <- as.vector(look)
def <- defData(varname = "C", formula = "..lookv", dist = "nonrandom")
# def <- defData(varname = "C", formula = 0.4, dist = "binary")
def <- defData(def, "X1", formula = "0.3 + 0.4 * C", dist = "nonrandom")
def <- defData(def, "X2", dist = "gamma", formula = 1.25, variance = 2)
def <- defData(def, "X3", dist = "gamma", formula = 3, variance = 5)
def <- defData(def, "X4", dist = "negBinomial", formula = 2, variance = 4)
def <- defData(def, "X5", dist = "negBinomial", formula = 4, variance = 4.2)
def <- defData(def, "X6", dist = "normal", formula = 1.5, variance = 3.6)
# def <- defData(def, "X", formula = "0.3 + 0.4 * C", dist = "binary")
def <- defData(def, "e", formula = 0, variance = 2, dist = "normal")
def <- defData(def, "Y0", formula = "2 * C + e", dist="nonrandom")
def <- defData(def, "Y1", formula = "1 + 2 * C + e", dist="nonrandom")
def <- defData(def, "Y_obs", formula = "Y0 + (Y1 - Y0) * X1", dist = "nonrandom")

# set.seed(202211)
gr <- expand.grid(long = seq(1, 100), lat = seq(1, 100))

easyfit_bkmr <- function(
  indata,
  yindex,
  zindex,
  xindex,
  family = c("gaussian", "binomial")
) {
  match.arg(family)
  indata <- as.matrix(indata)
  inY <- indata[, yindex]
  inZ <- indata[, zindex]
  if (is.null(xindex)) {inX <- NULL} else {inX <- indata[, xindex]}
  invisible(bkmr::kmbayes(X = inX, y = inY, Z = inZ, family = family, iter = 3000))
}

dtg <- genData(10000, def)

# not fit in 20 minutes
sim_bkmr <- easyfit_bkmr(
  dtg, yindex = 12, zindex = seq(3, 8), xindex = NULL, family = "gaussian"
)