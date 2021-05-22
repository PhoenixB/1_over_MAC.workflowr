## ---------------------------
##
## Script name: create_resamples.R
##
## Purpose of script: Create resamples for modelling
##
## Author: Pascal Jerney
##
## Date Created: 2021-05-18
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

## load data

one_over_mac_filtered <-
  read_rds(file = path_wd("output", "one_over_mac_filtered", ext = "rds"))

setkey(one_over_mac_filtered, pid, time)

## ---------------------------

# Create resamples

tictoc::tic(msg = "Resampling data", quiet = FALSE)

n_resamples <- 100
n_obs <- 100

resamples <-
  data.table(seed = 1:n_resamples)

setkey(resamples, seed)

resamples[, resample := map(seed,
                      function(seed, data, prop) {
                        set.seed(seed)
                        data[, slice_sample(.SD, n = n_obs), by = pid]
                        },
                      data = one_over_mac_filtered,
                      prop = prop)][]

resamples[, resample := map(resample, ~ setkey(.x, pid, time))]

tictoc::toc(log = TRUE)

tictoc::tic(msg = "Saving resamples", quiet = FALSE)

write_rds(resamples, file = path_wd("output", "resamples", ext = "rds"))

tictoc::toc(log = TRUE)
