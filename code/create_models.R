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

# Project working directory
setwd(here())

## ---------------------------

# Set pseudo-random number generator seed
set.seed(12345)

# I prefer to view outputs in non-scientific notation
options(scipen = 6, digits = 4)

# This is needed on some PCs to increase memory allowance, but has no impact on
# macs
suppressWarnings(memory.limit(30000000))

# Resolve function conflicts
conflict_prefer("select", "dplyr", quiet = TRUE)
conflict_prefer("filter", "dplyr", quiet = TRUE)
conflict_prefer("collapse", "dplyr", quiet = TRUE)

## ---------------------------

## load up our functions into memory

`%!in%` <- Negate(`%in%`)

## ---------------------------

## Load resamples

resamples <-
  read_rds(file = path_wd("output", "resamples", ext = "rds"))

resamples[, resample := map(resample, ~ setkey(.x, pid, time))]

## ---------------------------

# Null model

tictoc::tic(msg = "Modelling model 0", quiet = FALSE)

model0 <-
  resamples[,
            .(seed, model = map(resample,
                                function(data, ...) {
                                  args <- list(...)
                                  model <- gls(data = data, ...)
                                  model$call$model <- args$model
                                  model$call$control <- args$control
                                  model$data <- data
                                  model
                                  },
                                model = log(peak_alpha) ~ 1,
                                method = "ML",
                                control = glsControl()))
  ]

model0_specs <-
  model0[,
         .(seed,
           tidy = map(model,
                      ~ setDT(broom.mixed::tidy(.x))),
           performance = map(model,
                             ~ setDT(performance::model_performance(.x))),
           residuals = map(model,
                           ~ data.table(
                             raw = residuals(.x, type = "response"),
                             pearson = residuals(.x, type = "pearson"),
                             normalized = residuals(.x, type = "normalized"))),
           fitted = map(model,
                        ~ data.table(
                          fitted = fitted(.x))))
  ]

tictoc::toc(log = TRUE)

tictoc::tic(msg = "Saving model 0", quiet = FALSE)

# write_rds(model0, file = path_wd("output", "model0", ext = "rds"))
write_rds(model0_specs, file = path_wd("output", "model0_specs", ext = "rds"))

tictoc::toc(log = TRUE)

# Linear model with only centered ce_mac as predictor

tictoc::tic(msg = "Modelling model 1", quiet = FALSE)

model1 <-
  resamples[,
            .(seed, model = map(resample,
                                function(data, ...) {
                                  args <- list(...)
                                  model <- gls(data = data, ...)
                                  model$call$model = args$model
                                  model$call$control = args$control
                                  model
                                },
                                model = log(peak_alpha) ~ ce_mac_pid_centered,
                                method = "ML",
                                control = glsControl()))
  ]

model1_specs <-
  model1[,
         .(seed,
           tidy = map(model,
                      ~ setDT(broom.mixed::tidy(.x))),
           performance = map(model,
                             ~ setDT(performance::model_performance(.x))),
           residuals = map(model,
                           ~ data.table(
                             raw = residuals(.x, type = "response"),
                             pearson = residuals(.x, type = "pearson"),
                             normalized = residuals(.x, type = "normalized"))),
           fitted = map(model,
                        ~ data.table(
                          fitted = fitted(.x))))
  ]

tictoc::toc(log = TRUE)

tictoc::tic(msg = "Saving model 1", quiet = FALSE)

# write_rds(model1, file = path_wd("output", "model1", ext = "rds"))
write_rds(model1_specs, file = path_wd("output", "model1_specs", ext = "rds"))

tictoc::toc(log = TRUE)

# Mixed model with centered ce_mac as fixed effect and pid as random intercept

tictoc::tic(msg = "Modelling model 2", quiet = FALSE)

model2 <-
  resamples[,
            .(seed,
              model = map(resample,
                          function(data, ...) {
                            args = list(...)
                            model <- lme(data = data, ...)
                            model$call$fixed <- args$fixed
                            model$call$random <- args$random
                            model$call$control <- args$control
                            model
                          },
                          fixed = log(peak_alpha) ~ ce_mac_pid_centered,
                          random = ~ 1 | pid,
                          method = "ML",
                          control = lmeControl(maxIter = 100,
                                               returnObject = TRUE)))
  ]

model2_specs <-
  model2[,
         .(seed,
           tidy = map(model,
                      ~ setDT(broom.mixed::tidy(.x))),
           performance = map(model,
                             ~ setDT(performance::model_performance(.x))),
           residuals = map(model,
                           ~ data.table(
                             raw = residuals(.x, type = "response"),
                             pearson = residuals(.x, type = "pearson"),
                             normalized = residuals(.x, type = "normalized"))),
           fitted = map(model,
                        ~ data.table(
                          fitted = fitted(.x, level = 0:1))),
           ranef = map(model,
                        ~ data.table(
                          ranef = ranef(.x))))
  ]

tictoc::toc(log = TRUE)

tictoc::tic(msg = "Saving model 2", quiet = FALSE)

# write_rds(model2, file = path_wd("output", "model2", ext = "rds"))
write_rds(model2_specs, file = path_wd("output", "model2_specs", ext = "rds"))

tictoc::toc(log = TRUE)

# Mixed model with centered ce_mac as fixed effect, pid as random intercept and
# centered ce_mac as random slope

tictoc::tic(msg = "Modelling model 3", quiet = FALSE)

model3 <-
  resamples[,
            .(seed,
              model = map(resample,
                          function(data, ...) {
                            args = list(...)
                            model <- lme(data = data, ...)
                            model$call$fixed <- args$fixed
                            model$call$random <- args$random
                            model$call$control <- args$control
                            model
                          },
                          fixed = log(peak_alpha) ~ ce_mac_pid_centered,
                          random = ~ ce_mac_pid_centered | pid,
                          method = "ML",
                          control = lmeControl(maxIter = 100,
                                               returnObject = TRUE)))
  ]

model3_specs <-
  model3[,
         .(seed,
           tidy = map(model,
                      ~ setDT(broom.mixed::tidy(.x))),
           performance = map(model,
                             ~ setDT(performance::model_performance(.x))),
           residuals = map(model,
                           ~ data.table(
                             raw = residuals(.x, type = "response"),
                             pearson = residuals(.x, type = "pearson"),
                             normalized = residuals(.x, type = "normalized"))),
           fitted = map(model,
                        ~ data.table(
                          fitted = fitted(.x, level = 0:1))),
           ranef = map(model,
                       ~ data.table(
                         ranef = ranef(.x))))
  ]

tictoc::toc(log = TRUE)

tictoc::tic(msg = "Saving model 3", quiet = FALSE)

# write_rds(model3, file = path_wd("output", "model3", ext = "rds"))
write_rds(model3_specs, file = path_wd("output", "model3_specs", ext = "rds"))

tictoc::toc(log = TRUE)

# Mixed model with centered ce_mac as fixed effect, pid as random intercept,
# centered ce_mac as random slope and AR1 process within pid

tictoc::tic(msg = "Modelling model 4", quiet = FALSE)

model4 <-
  resamples[,
            .(seed,
              model = map(resample,
                          function(data, ...) {
                            args = list(...)
                            model <- lme(data = data, ...)
                            model$call$fixed <- args$fixed
                            model$call$random <- args$random
                            model$call$control <- args$control
                            model$call$correlation <- args$correlation
                            model
                          },
                          fixed = log(peak_alpha) ~ ce_mac_pid_centered,
                          random = ~ ce_mac_pid_centered | pid,
                          correlation = corAR1(form = ~ time),
                          method = "ML",
                          control = lmeControl(maxIter = 100,
                                               returnObject = TRUE)))
  ]

model4_specs <-
  model4[,
         .(seed,
           tidy = map(model,
                      ~ setDT(broom.mixed::tidy(.x))),
           performance = map(model,
                             ~ setDT(performance::model_performance(.x))),
           residuals = map(model,
                           ~ data.table(
                             raw = residuals(.x, type = "response"),
                             pearson = residuals(.x, type = "pearson"),
                             normalized = residuals(.x, type = "normalized"))),
           fitted = map(model,
                        ~ data.table(
                          fitted = fitted(.x, level = 0:1))),
           ranef = map(model,
                       ~ data.table(
                         ranef = ranef(.x))))
  ]

tictoc::toc(log = TRUE)

tictoc::tic(msg = "Saving model 4", quiet = FALSE)

# write_rds(model4, file = path_wd("output", "model4", ext = "rds"))
write_rds(model4_specs, file = path_wd("output", "model4_specs", ext = "rds"))

tictoc::toc(log = TRUE)

# Mixed model with centered ce_mac as fixed effect, pid as random intercept,
# centered ce_mac as random slope and continuous AR1 process within pid

tictoc::tic(msg = "Modelling model 5", quiet = FALSE)

model5 <-
  resamples[,
            .(seed,
              model = map(resample,
                          function(data, ...) {
                            args = list(...)
                            model <- lme(data = data, ...)
                            model$call$fixed <- args$fixed
                            model$call$random <- args$random
                            model$call$control <- args$control
                            model$call$correlation <- args$correlation
                            model
                          },
                          fixed = log(peak_alpha) ~ ce_mac_pid_centered,
                          random = ~ ce_mac_pid_centered | pid,
                          correlation = corCAR1(form = ~ time),
                          method = "ML",
                          control = lmeControl(maxIter = 100,
                                               returnObject = TRUE)))
  ]

model5_specs <-
  model5[,
         .(seed,
           tidy = map(model,
                      ~ setDT(broom.mixed::tidy(.x))),
           performance = map(model,
                             ~ setDT(performance::model_performance(.x))),
           residuals = map(model,
                           ~ data.table(
                             raw = residuals(.x, type = "response"),
                             pearson = residuals(.x, type = "pearson"),
                             normalized = residuals(.x, type = "normalized"))),
           fitted = map(model,
                        ~ data.table(
                          fitted = fitted(.x, level = 0:1))),
           ranef = map(model,
                       ~ data.table(
                         ranef = ranef(.x))))
  ]

tictoc::toc(log = TRUE)

tictoc::tic(msg = "Saving model 5", quiet = FALSE)

write_rds(model5, file = path_wd("output", "model5", ext = "rds"))
write_rds(model5_specs, file = path_wd("output", "model5_specs", ext = "rds"))

tictoc::toc(log = TRUE)

# Mixed model with centered ce_mac as fixed effect, pid as random intercept,
# centered ce_mac as random slope and MA2 process within pid

tictoc::tic(msg = "Modelling model 6", quiet = FALSE)

model6 <-
  resamples[,
            .(seed,
              model = map(resample,
                          function(data, ...) {
                            args = list(...)
                            model <- lme(data = data, ...)
                            model$call$fixed <- args$fixed
                            model$call$random <- args$random
                            model$call$control <- args$control
                            model$call$correlation <- args$correlation
                            model
                          },
                          fixed = log(peak_alpha) ~ ce_mac_pid_centered,
                          random = ~ ce_mac_pid_centered | pid,
                          correlation = corARMA(form = ~ time, q = 2),
                          method = "ML",
                          control = lmeControl(maxIter = 100,
                                               returnObject = TRUE)))
  ]

model6_specs <-
  model6[,
         .(seed,
           tidy = map(model,
                      ~ setDT(broom.mixed::tidy(.x))),
           performance = map(model,
                             ~ setDT(performance::model_performance(.x))),
           residuals = map(model,
                           ~ data.table(
                             raw = residuals(.x, type = "response"),
                             pearson = residuals(.x, type = "pearson"),
                             normalized = residuals(.x, type = "normalized"))),
           fitted = map(model,
                        ~ data.table(
                          fitted = fitted(.x, level = 0:1))),
           ranef = map(model,
                       ~ data.table(
                         ranef = ranef(.x))))
  ]

tictoc::toc(log = TRUE)

tictoc::tic(msg = "Saving model 6", quiet = FALSE)

# write_rds(model6, file = path_wd("output", "model6", ext = "rds"))
write_rds(model6_specs, file = path_wd("output", "model6_specs", ext = "rds"))

tictoc::toc(log = TRUE)

# Mixed model with centered ce_mac as fixed effect, pid as random intercept,
# centered ce_mac as random slope and ARMA1,1 process within pid

tictoc::tic(msg = "Modelling model 7", quiet = FALSE)

model7 <-
  resamples[,
            .(seed,
              model = map(resample,
                          function(data, ...) {
                            args = list(...)
                            model <- lme(data = data, ...)
                            model$call$fixed <- args$fixed
                            model$call$random <- args$random
                            model$call$control <- args$control
                            model$call$correlation <- args$correlation
                            model
                          },
                          fixed = log(peak_alpha) ~ ce_mac_pid_centered,
                          random = ~ ce_mac_pid_centered | pid,
                          correlation = corARMA(form = ~ time, p = 1,
                                                q = 1, fixed = TRUE),
                          method = "ML",
                          control = lmeControl(maxIter = 100,
                                               returnObject = TRUE)))
  ]

model7_specs <-
  model7[,
         .(seed,
           tidy = map(model,
                      ~ setDT(broom.mixed::tidy(.x))),
           performance = map(model,
                             ~ setDT(performance::model_performance(.x))),
           residuals = map(model,
                           ~ data.table(
                             raw = residuals(.x, type = "response"),
                             pearson = residuals(.x, type = "pearson"),
                             normalized = residuals(.x, type = "normalized"))),
           fitted = map(model,
                        ~ data.table(
                          fitted = fitted(.x, level = 0:1))),
           ranef = map(model,
                       ~ data.table(
                         ranef = ranef(.x))))
  ]

tictoc::toc(log = TRUE)

tictoc::tic(msg = "Saving model 7", quiet = FALSE)

# write_rds(model7, file = path_wd("output", "model7", ext = "rds"))
write_rds(model7_specs, file = path_wd("output", "model7_specs", ext = "rds"))

tictoc::toc(log = TRUE)

tictoc::tic.log()
