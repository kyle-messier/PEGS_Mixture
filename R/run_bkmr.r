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
  .[,!apply(., 2, \(x) any(is.na(x)))] %>%
  .[,-17] %>%
  .[,!apply(., 2, \(x) all(x == 0))]

X_in <- epr_console_sub_noe[, -..targ_interaction] %>%
  .[,-1] %>%
  as.matrix() %>%
  .[,!apply(., 2, \(x) any(is.na(x)))] %>%
  # remove covariates with "excessive" zeros (i.e., 2000+ / 2751)
  # reduces the # covariates dramatically to five
  .[,apply(., 2, \(x) sum(x == 0) < 2700)]
dim(X_in)

# check the data (i.e., X_in)
testbkmr <- bkmr::kmbayes(y = y_in, Z = Z_in, X = X_in, family = "binomial", varsel = TRUE)
# 10%: 1.12396 hours...

# Fitted object of class 'bkmrfit'
# Iterations: 1000 
# Outcome family: binomial (probit link) 
# Model fit on: 2023-10-24 17:10:49.467435 
# Running time:  5.86523 hours 

# Acceptance rates for Metropolis-Hastings algorithm:
#               param      rate
# 1            lambda 0.2432432
# 2 r/delta (overall) 0.4474474
# 3 r/delta  (move 1) 0.6641509
# 4 r/delta  (move 2) 0.2025586

# Parameter estimates (based on iterations 501-1000):
#        param     mean      sd    q_2.5   q_97.5
# 1      beta1  0.00003 0.00036 -0.00067  0.00070
# 2      beta2 -0.00016 0.00028 -0.00071  0.00041
# 3      beta3  0.00027 0.00034 -0.00038  0.00103
# 4      beta4  0.00023 0.00006  0.00012  0.00035
# 5      beta5  0.00007 0.00009 -0.00012  0.00025
# 6  sigsq.eps  1.00000 0.00000  1.00000  1.00000
# 7         r1  0.01168 0.01865  0.00000  0.06973
# 8         r2  0.03308 0.15350  0.00000  0.08467
# 9         r3  0.00854 0.01471  0.00000  0.05890
# 10        r4  0.06164 0.38266  0.00000  0.25699
# 11    lambda  3.92528 5.23059  0.75361 19.66435
# 12    ystar1 -1.08635 0.73720 -2.74390 -0.05360
# 13    ystar2  0.61001 0.48697  0.02041  1.76516
# 14 ystar2751  0.56680 0.46340  0.02678  1.78262

# single parent or group quarter inclusion p > 0.5
# Posterior inclusion probabilities:
#              variable   PIP
# 1    svi_2010_P_AGE17 0.442
# 2  svi_2010_P_SNGPRNT 0.610
# 3 svi_2010_P_MINORITY 0.394
# 4   svi_2010_P_GROUPQ 0.562
