source("./R/prep_bkmr_data.r")

pkgs <- c("bkmr", "data.table", "dplyr", "tidytable", "sf", "broom")
invisible(sapply(pkgs, library, character.only = TRUE))

sf_use_s2(FALSE)
Sys.setenv("TIGRIS_CACHE_DIR" = sprintf("%s/tigris_cache/", Sys.getenv("HOME")))

## base data is epr_console
names(epr_console)

epr_console_sub <- epr_console %>%
  select(flag_antidep_1p, 18:34, 36:593) %>%
  mutate(across(everything(), ~ifelse(is.na(.), 0, as.numeric(.))))
epr_console_mv <- sapply(epr_console_sub, function(x) length(unique(x)) > 1)
epr_console_sub <- epr_console_sub[, ..epr_console_mv] %>%
  mutate(across(seq_len(ncol(.))[length(unique(.)) > 2], ~as.vector(scale(.))))

epr_form <- reformulate(
  termlabels = names(epr_console_sub)[-1],
  response = names(epr_console_sub)[1]
)

## no-brainer logistic reg
epr_glm <- glm(epr_form, binomial, data.frame(epr_console_sub))
broom::glance(epr_glm)
epr_glm_res <- broom::tidy(epr_glm)
tail(epr_glm_res)


## TODO: tidy codes
epr_console_sub_noe <-
  epr_console_sub %>%
  select(-contains("_E_"), -contains("EPL"))
names(epr_console_sub_noe)

y_in <- epr_console_sub_noe[,1] %>%
  as.matrix()
targ_interaction <- seq(81, 84)
Z_in <- epr_console_sub_noe[, ..targ_interaction] %>%
  as.matrix()
X_in <- epr_console_sub_noe[, -..targ_interaction] %>%
  .[,-1] %>%
  as.matrix() %>%
  .[,-apply(., 2, \(x) any(is.na(x)))] %>%
  .[,-17] %>%
  .[,-apply(., 2, \(x) all(x == 0))]

# Error: the leading minor of order 29 is not positive
# check the data (i.e., X_in)
testbkmr <- bkmr::kmbayes(y = y_in, Z = Z_in, X = X_in, family = "binomial")
