source("R/prep_participant_data.r")

dirlist <- read.csv("~/.songi2init")
snps <- read.csv(file.path(dirlist$val[8], "snp_presence_absence_binary.csv"))
snps$FID <- as.character(snps$FID)

# use epr_main
pegs_main <- epr_main |>
  inner_join(snps, by = c("epr_number" = "FID"))
# 1,575


# perceived stress score
# reverse scale 4,5,7,8
# 1469 excluding no response in any of eb_e09* questions
pegs_main_pss <- pegs_main %>%
  # 1-5 to 0-4
  mutate(
    across(starts_with("eb_e09"), ~. - 1)
  ) %>%
  mutate(
    eb_e093_confidence = 4 - eb_e093_confidence,
    eb_e094_going_your_way = 4 - eb_e094_going_your_way,
    eb_e096_irritated = 4 - eb_e096_irritated,
    eb_e097_on_top_of_things = 4 - eb_e097_on_top_of_things
  ) %>%
  mutate(pss = rowSums(select(., starts_with("eb_e09")))) %>%
  mutate(pss_class = cut(pss, c(0, 16, 26, 40), labels = c("low", "medium", "high"))) %>%
  filter(!is.na(pss_class))
  # filter(across(everything(), ~!is.na(.)))

pegs_main_psso <- pegs_main_pss
pegs_main_pss <- pegs_main_pss[, sapply(pegs_main_pss, function(x) sum(x == 0) < 1e3)]

dim(pegs_main_pss)
summary(pegs_main_pss$pss_class)


##
library(qgcompint)
library(nnet)

explanatory <- c(688, 689, 494:533, 691:703)

mainform <- reformulate(
  termlabels = names(pegs_main_pss)[explanatory],
  response = "pss_class"
)

pegs_main_pss <-
  pegs_main_pss |>
  select(-geometry)
pegs_main_pss_com <- pegs_main_pss[complete.cases(pegs_main_pss), ]

mainmult <- MASS::polr(mainform, data = pegs_main_pss_com, Hess = TRUE)
# summary(mainmult)

# gtsummary::tbl_regression(mainmult)
# mainmult
# broom::tidy(mainmult)
