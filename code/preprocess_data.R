## ---------------------------
##
## Script name: preprocess_data.R
##
## Purpose of script: Preprocess EEG data for further use
##
## Author: Pascal Jerney
##
## Date Created: 2021-02-26
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
  "tidymodels",
  "readxl",
  "fs",
  "readr",
  "here",
  "lubridate",
  "janitor"
)

## ---------------------------

## set working directory

setwd(here()) # Project working directory

## ---------------------------

set.seed(12345) # Set pseudo-random number generator seed
options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation
suppressWarnings(memory.limit(30000000)) # this is needed on some PCs to increase memory allowance, but has no impact on macs.
conflict_prefer("filter", "dplyr") # Set preferred functions of conflicting packages
conflict_prefer("discard", "purrr")

## ---------------------------

## load up our functions into memory

source("code/group_sampling.R")

rangediff <- function(x) {
  range(x)[2] - range(x)[1]
}

`%!in%` <- Negate(`%in%`)

## ---------------------------

## Create helpers ----
helper.colnames <- c(
  "time",
  "ce_mac",
  "peak_alpha",
  "osc_alpha_power",
  "osc_alpha"
)

helper.colspec <- cols(
  time = col_integer(),
  ce_mac = col_double(),
  peak_alpha = col_double(),
  osc_alpha_power = col_double(),
  osc_alpha = col_logical()
)

helper.data_path = path_wd(
  "data"
)

helper.rejection_files <- dir_ls(path = helper.data_path, regexp = "to_reject_from_analysis_.*?\\.xlsx$")

tibble.rejection_list <-
  helper.rejection_files %>%
  map(function(file_path) {
    read_excel(
      file_path,
      col_names = c("pid"),
      # col_names = helper.colnames,
      # col_types = helper.colspec,
      na = c("", "NA", "NaN"),
      skip = 1
    )
  }) %>%
  bind_rows(.id = "filename") %>%
  mutate(file_date = dmy(str_replace(filename, ".*?_(((((0?[1-9])|(1\\d)|(2[0-8]))\\.((0?[1-9])|(1[0-2])))|((31\\.((0[13578])|(1[02])))|((29|30)\\.((0?[1,3-9])|(1[0-2])))))\\.((20[0-9][0-9]))|(29\\.0?2\\.20(([02468][048])|([13579][26]))))\\.xlsx$", "\\1"))) %>%
  # modify_at(vars(file_date), as_date) %>%
  select(file_date, pid) %>%
  modify_at(vars(pid), as_factor) %>%
  filter(file_date == max(file_date))

helper.data_files <- dir_ls(path = helper.data_path, regexp = "\\.csv$")

## Load patient data from csv files into tibble and save to temporary RDS file ----
tibble.1overmac_list <-
  helper.data_files %>%
  map(function(file_path) {
    read_csv(
      file_path,
      col_names = helper.colnames,
      col_types = helper.colspec,
      na = c("", "NA", "NaN")
    )
  }) %>%
  discard(~ nrow(.x) == 0)

## Create recipe for data tidying ----
recipe.peakalpha_base <-
  recipe(
    tibble.1overmac_list[[1]]
  )

recipe.peakalpha_movingmean <-
  recipe.peakalpha_base %>%
  step_window(peak_alpha, size = 3, statistic = "mean", names = "peak_alpha_mean_3pt") %>%
  step_window(peak_alpha, size = 5, statistic = "mean", names = "peak_alpha_mean_5pt") %>%
  step_window(peak_alpha, size = 7, statistic = "mean", names = "peak_alpha_mean_7pt") %>%
  step_window(peak_alpha, size = 15, statistic = "mean", names = "peak_alpha_mean_15pt") %>%
  step_window(ce_mac, size = 5, statistic = "mean", names = "ce_mac_mean_5pt") %>%
  step_window(ce_mac, size = 15, statistic = "mean", names = "ce_mac_mean_15pt")

recipe.peakalpha_prepped <-
  recipe.peakalpha_movingmean %>%
  prep(training = tibble.1overmac_list[[1]])

tibble.1overmac_raw <-
  tibble.1overmac_list %>%
  map(~ bake(object = recipe.peakalpha_prepped, new_data = .x)) %>%
  imap(~ mutate(.x, filename = as.character(.y))) %>%
  bind_rows()

tibble.1overmac_tidy <-
  tibble.1overmac_raw %>%
  mutate(pid = str_replace(filename, ".*?_([0-9]+)\\.csv$", "\\1")) %>%
  modify_at(vars(pid), as_factor) %>%
  modify_at(vars(pid), fct_inseq) %>%
  select(pid, time, ce_mac, ce_mac_mean_5pt, ce_mac_mean_15pt, peak_alpha, peak_alpha_mean_3pt, peak_alpha_mean_5pt, peak_alpha_mean_7pt, peak_alpha_mean_15pt) %>%
  group_by(pid) %>%
  mutate(
    peak_alpha_mz = (0.6745 * (peak_alpha - median(peak_alpha, na.rm = TRUE))) / (mad(peak_alpha, na.rm = TRUE)),
    peak_alpha_outlier = as.logical(abs(peak_alpha_mz) > 3.5)
  ) %>%
  filter(ce_mac != 0) %>%
  arrange(pid, time)

# write_rds(tibble.1overmac_tidy, file = path_wd("output", "1overmac_tidy", ext = "rds"))

## Exclude patients from exclusion list ----
tibble.1overmac_excluded <-
  tibble.1overmac_tidy %>%
  anti_join(tibble.rejection_list, by = "pid")

# write_rds(tibble.1overmac_excluded, file = path_wd("output", "1overmac_excluded", ext = "rds"))

## Create helper for the median of ce_mac_range ----
helper.cemac_range_median <-
  tibble.1overmac_excluded %>%
  group_by(pid) %>%
  summarise(across(ce_mac, list(range = rangediff)), .groups = "drop_last") %>%
  summarise(ce_mac_median = median(ce_mac_range)) %>%
  pull(ce_mac_median)

## Create helper for splitting data at the median of ce_mac_range ----
helper.cemac_range_split <-
  tibble.1overmac_excluded %>%
  group_by(pid) %>%
  summarise(across(ce_mac, list(range = rangediff)), .groups = "drop_last") %>%
  filter(ce_mac_range >= helper.cemac_range_median) %>%
  pull(pid)

tibble.1overmac_filtered <-
  tibble.1overmac_excluded %>%
  filter(pid %in% helper.cemac_range_split) %>%
  # filter(!peak_alpha_outlier) %>%
  select(-peak_alpha_outlier, -peak_alpha_mz) %>%
  drop_na(peak_alpha) %>%
  modify_at(vars(pid), fct_drop)

write_rds(tibble.1overmac_filtered, file = path_wd("output", "1overmac_filtered", ext = "rds"))

## Create helper for training and test set ----
helper.training_test <-
  tibble.1overmac_filtered %>%
  distinct(pid, .keep_all = FALSE) %>%
  initial_split(prop = 1/2)

# helper.training_test_pid <-
#   tibble.1overmac_thinned_15 %>%
#   distinct(pid, .keep_all = TRUE) %>%
#   rownames_to_column() %>%
#   filter(rowname %in% helper.training_test) %>%
#   pull(pid)

## Create training and test sets ----
tibble.1overmac_train <-
  tibble.1overmac_filtered %>%
  filter(as.character(pid) %in% (training(helper.training_test) %>% pluck("pid"))) %>%
  modify_at(vars(pid), fct_drop)

write_rds(tibble.1overmac_train, file = path_wd("output", "1overmac_train", ext = "rds"))

tibble.1overmac_test <-
  tibble.1overmac_filtered %>%
  filter(as.character(pid) %in% (testing(helper.training_test) %>% pluck("pid"))) %>%
  modify_at(vars(pid), fct_drop)

write_rds(tibble.1overmac_test, file = path_wd("output", "1overmac_test", ext = "rds"))

## Thin training set ----
tibble.1overmac_train_thinned_5 <-
  tibble.1overmac_train %>%
  filter(row_number() %% 5 == 1) %>%
  modify_at(vars(pid), fct_drop)

write_rds(tibble.1overmac_train_thinned_5, file = path_wd("output", "1overmac_train_thinned_5", ext = "rds"))

tibble.1overmac_train_thinned_15 <-
  tibble.1overmac_train %>%
  filter(row_number() %% 15 == 1) %>%
  modify_at(vars(pid), fct_drop)

write_rds(tibble.1overmac_train_thinned_15, file = path_wd("output", "1overmac_train_thinned_15", ext = "rds"))
