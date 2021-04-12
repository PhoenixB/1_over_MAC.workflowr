## ---------------------------
##
## Script name: create_models.R
##
## Purpose of script: Create models for further use
##
## Author: Pascal Jerney
##
## Date Created: 2021-03-25
##
## Copyright (c) Pascal Jerney, 2021
## Email: pascal@jerney.ch
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)

if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(
  "conflicted",
  "tidyverse",
  "readr",
  "here",
  "glue",
  "fs",
  "nlme",
  "broom.mixed",
  "data.table",
  "tictoc"
)

## ---------------------------

setwd(here())                                        # Project working directory

## ---------------------------

set.seed(12345)                        # Set pseudo-random number generator seed
options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation
suppressWarnings(memory.limit(30000000)) # this is needed on some PCs to increase memory allowance, but has no impact on macs.
conflict_prefer("select", "dplyr", quiet = TRUE)
conflict_prefer("filter", "dplyr", quiet = TRUE)
conflict_prefer("collapse", "dplyr", quiet = TRUE)

## ---------------------------

## load up our functions into memory

`%!in%` <- Negate(`%in%`)

## ---------------------------

## load data

one_over_mac_filtered <-
  read_rds(file = path_wd("output", "one_over_mac_filtered", ext = "rds"))

setkey(one_over_mac_filtered, pid, time)

## ---------------------------

# Create resamples

tictoc::tic(msg = "Resampling data", quiet = FALSE)

n_resamples <- 100
prop <- 0.01

resamples <-
  data.table(seed = 1:n_resamples)

setkey(resamples, seed)

resamples[, resample := map(seed,
                      function(seed, data, prop) {
                        set.seed(seed)
                        data[, slice_sample(.SD, prop = prop), by = pid]
                        },
                      data = one_over_mac_filtered,
                      prop = prop)][]

resamples[, resample := map(resample, ~ setkey(.x, pid, time))]

tictoc::toc(log = TRUE)

tictoc::tic(msg = "Saving resamples", quiet = FALSE)

write_rds(resamples, file = path_wd("output", "resamples", ext = "rds"))

tictoc::toc(log = TRUE)

# Calculate average subsample size
mean_subsample_size <-
  resamples[, .(subsample_size = map_dbl(resample, ~ nrow(.x)))
            ][, mean(subsample_size)]

# Calculate effective sample fraction
eff_sample_frac <- mean_subsample_size / nrow(one_over_mac_filtered)

# Save effective sample fraction
write_rds(eff_sample_frac,
          file = path_wd("output", "eff_sample_frac", ext = "rds"))

# Null model

tictoc::tic(msg = "Modelling model 0", quiet = FALSE)

model0 <-
  resamples[,
            .(seed, resample, model = map(resample,
                         ~ gls(
                           log(peak_alpha) ~ 1,
                           data = .x,
                           method = "REML",
                           control = glsControl())))
  ][,
    .(seed, resample, model,
      tidy = map(model,
                 ~ setDT(broom.mixed::tidy(.x))))
  ][,
    .(seed, resample, model, tidy,
      performance = map(model,
                        ~ data.table(
                          rmse = performance::performance_rmse(.x),
                          mse = performance::performance_mse(.x))))
  ][,
    .(seed, resample, model, tidy, performance,
      glance = map(model,
                   ~ setDT(broom.mixed::glance(.x))))
  ]

tictoc::toc(log = TRUE)

tictoc::tic(msg = "Saving model 0", quiet = FALSE)

write_rds(model0, file = path_wd("output", "model0", ext = "rds"))

tictoc::toc(log = TRUE)

# Linear model with only centered ce_mac as predictor

tictoc::tic(msg = "Modelling model 1", quiet = FALSE)

model1 <-
  resamples[,
            .(seed, resample,
              model = map(resample,
                          ~ gls(
                            log(peak_alpha) ~ ce_mac_pid_centered,
                            data = .x,
                            method = "REML",
                            control = glsControl())))
  ][,
    .(seed, resample, model,
      tidy = map(model,
                 ~ setDT(broom.mixed::tidy(.x))))
  ][,
    .(seed, resample, model, tidy,
      performance = map(model,
                        ~ data.table(
                          rmse = performance::performance_rmse(.x),
                          mse = performance::performance_mse(.x))))
  ][,
    .(seed, resample, model, tidy, performance,
      glance = map(model,
                   ~ setDT(broom.mixed::glance(.x))))
  ]

tictoc::toc(log = TRUE)

tictoc::tic(msg = "Saving model 1", quiet = FALSE)

write_rds(model1, file = path_wd("output", "model1", ext = "rds"))

tictoc::toc(log = TRUE)

# Mixed model with centered ce_mac as fixed effect and pid as random intercept

tictoc::tic(msg = "Modelling model 2", quiet = FALSE)

model2 <-
  resamples[,
            .(seed, resample,
              model = map(resample,
                          ~ lme(
                            log(peak_alpha) ~ ce_mac_pid_centered,
                            data = .x,
                            random = ~ 1 | pid,
                            method = "REML",
                            control = lmeControl(maxIter = 100)
                            )))
  ][,
    .(seed, resample, model,
      tidy = map(model,
                 ~ setDT(broom.mixed::tidy(.x))))
  ][,
    .(seed, resample, model, tidy,
      performance = map(model,
                        ~ data.table(
                          rmse = performance::performance_rmse(.x),
                          mse = performance::performance_mse(.x))))
  ][,
    .(seed, resample, model, tidy, performance,
      glance = map(model,
                   ~ setDT(broom.mixed::glance(.x))))
  ]

tictoc::toc(log = TRUE)

tictoc::tic(msg = "Saving model 2", quiet = FALSE)

write_rds(model2, file = path_wd("output", "model2", ext = "rds"))

tictoc::toc(log = TRUE)

# Mixed model with centered ce_mac as fixed effect, pid as random intercept and
# centered ce_mac as random slope

tictoc::tic(msg = "Modelling model 3", quiet = FALSE)

model3 <-
  resamples[,
            .(seed, resample,
              model = map(resample,
                          ~ lme(
                            log(peak_alpha) ~ ce_mac_pid_centered,
                            data = .x,
                            random = ~ ce_mac_pid_centered | pid,
                            method = "REML",
                            control = lmeControl(maxIter = 100)
                            )))
  ][,
    .(seed, resample, model,
      tidy = map(model,
                 ~ setDT(broom.mixed::tidy(.x))))
  ][,
    .(seed, resample, model, tidy,
      performance = map(model,
                        ~ data.table(
                          rmse = performance::performance_rmse(.x),
                          mse = performance::performance_mse(.x))))
  ][,
    .(seed, resample, model, tidy, performance,
      glance = map(model,
                   ~ setDT(broom.mixed::glance(.x))))
  ]

tictoc::toc(log = TRUE)

tictoc::tic(msg = "Saving model 3", quiet = FALSE)

write_rds(model3, file = path_wd("output", "model3", ext = "rds"))

tictoc::toc(log = TRUE)

# Mixed model with centered ce_mac as fixed effect, pid as random intercept,
# centered ce_mac as random slope and AR1 process within pid

tictoc::tic(msg = "Modelling model 4", quiet = FALSE)

model4 <-
  resamples[,
            .(seed, resample,
              model = map(resample,
                          ~ lme(
                            log(peak_alpha) ~ ce_mac_pid_centered,
                            data = .x,
                            random = ~ ce_mac_pid_centered | pid,
                            correlation = corAR1(form = ~ time),
                            method = "REML",
                            control = lmeControl(maxIter = 100)
                            )))
  ][,
    .(seed, resample, model,
      tidy = map(model,
                 ~ setDT(broom.mixed::tidy(.x))))
  ][,
    .(seed, resample, model, tidy,
      performance = map(model,
                        ~ data.table(
                          rmse = performance::performance_rmse(.x),
                          mse = performance::performance_mse(.x))))
  ][,
    .(seed, resample, model, tidy, performance,
      glance = map(model,
                   ~ setDT(broom.mixed::glance(.x))))
  ]

tictoc::toc(log = TRUE)

tictoc::tic(msg = "Saving model 4", quiet = FALSE)

write_rds(model4, file = path_wd("output", "model4", ext = "rds"))

tictoc::toc(log = TRUE)

# Mixed model with centered ce_mac as fixed effect, pid as random intercept,
# centered ce_mac as random slope and continuous AR1 process within pid

tictoc::tic(msg = "Modelling model 5", quiet = FALSE)

model5 <-
  resamples[,
            .(seed, resample,
              model = map(resample,
                          ~ lme(
                            log(peak_alpha) ~ ce_mac_pid_centered,
                            data = .x,
                            random = ~ ce_mac_pid_centered | pid,
                            correlation = corCAR1(form = ~ time),
                            method = "REML",
                            control = lmeControl(maxIter = 100)
                            )))
  ][,
    .(seed, resample, model,
      tidy = map(model,
                 ~ setDT(broom.mixed::tidy(.x))))
  ][,
    .(seed, resample, model, tidy,
      performance = map(model,
                        ~ data.table(
                          rmse = performance::performance_rmse(.x),
                          mse = performance::performance_mse(.x))))
  ][,
    .(seed, resample, model, tidy, performance,
      glance = map(model,
                   ~ setDT(broom.mixed::glance(.x))))
  ]

tictoc::toc(log = TRUE)

tictoc::tic(msg = "Saving model 5", quiet = FALSE)

write_rds(model5, file = path_wd("output", "model5", ext = "rds"))

tictoc::toc(log = TRUE)

# Mixed model with centered ce_mac as fixed effect, pid as random intercept,
# centered ce_mac as random slope and MA2 process within pid

tictoc::tic(msg = "Modelling model 6", quiet = FALSE)

model6 <-
  resamples[,
            .(seed, resample,
              model = map(resample,
                          ~ lme(
                            log(peak_alpha) ~ ce_mac_pid_centered,
                            data = .x,
                            random = ~ ce_mac_pid_centered | pid,
                            correlation = corARMA(form = ~ time, q = 2),
                            method = "REML",
                            control = lmeControl(maxIter = 100)
                          )))
  ][,
    .(seed, resample, model,
      tidy = map(model,
                 ~ setDT(broom.mixed::tidy(.x))))
  ][,
    .(seed, resample, model, tidy,
      performance = map(model,
                        ~ data.table(
                          rmse = performance::performance_rmse(.x),
                          mse = performance::performance_mse(.x))))
  ][,
    .(seed, resample, model, tidy, performance,
      glance = map(model,
                   ~ setDT(broom.mixed::glance(.x))))
  ]

tictoc::toc(log = TRUE)

tictoc::tic(msg = "Saving model 6", quiet = FALSE)

write_rds(model6, file = path_wd("output", "model6", ext = "rds"))

tictoc::toc(log = TRUE)

# Mixed model with centered ce_mac as fixed effect, pid as random intercept,
# centered ce_mac as random slope and ARMA1,1 process within pid

tictoc::tic(msg = "Modelling model 7", quiet = FALSE)

model7 <-
  resamples[,
            .(seed, resample,
              model = map(resample,
                          ~ lme(
                            log(peak_alpha) ~ ce_mac_pid_centered,
                            data = .x,
                            random = ~ ce_mac_pid_centered | pid,
                            correlation = corARMA(form = ~ time, p = 1,
                                                  q = 1, fixed = TRUE),
                            method = "REML",
                            control = lmeControl(maxIter = 100)
                          )))
  ][,
    .(seed, resample, model,
      tidy = map(model,
                 ~ setDT(broom.mixed::tidy(.x))))
  ][,
    .(seed, resample, model, tidy,
      performance = map(model,
                        ~ data.table(
                          rmse = performance::performance_rmse(.x),
                          mse = performance::performance_mse(.x))))
  ][,
    .(seed, resample, model, tidy, performance,
      glance = map(model,
                   ~ setDT(broom.mixed::glance(.x))))
  ]

tictoc::toc(log = TRUE)

tictoc::tic(msg = "Saving model 7", quiet = FALSE)

write_rds(model7, file = path_wd("output", "model7", ext = "rds"))

tictoc::toc(log = TRUE)

tictoc::tic.log()
