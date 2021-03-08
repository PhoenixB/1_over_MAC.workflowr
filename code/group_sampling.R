## ---------------------------
##
## Script name: group_sampling.R
##
## Purpose of script: Create functions for sampling within groups
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
  "tidyverse",
  "conflicted",
  "here"
)

## ---------------------------

setwd(here()) # Project working directory

## ---------------------------

sample_n_groups <- function(df, grp_var, n, replace = FALSE, weight = NULL) {
  grp_var <- enquo(grp_var)
  random_grp <- df %>%
    group_by(!!grp_var) %>%
    summarise(.groups = "drop_last") %>%
    slice_sample(n = n, weight_by = weight, replace = replace) %>%
    mutate(unique_id = 1:NROW(.))
  df %>%
    right_join(random_grp, by=quo_name(grp_var)) %>%
    group_by(!!grp_var)
}

sample_prop_groups <- function(df, grp_var, prop, replace = FALSE, weight = NULL) {
  grp_var <- enquo(grp_var)
  random_grp <- df %>%
    group_by(!!grp_var) %>%
    summarise(.groups = "drop_last") %>%
    slice_sample(prop = prop, weight_by = weight, replace = replace) %>%
    mutate(unique_id = 1:NROW(.))
  df %>%
    right_join(random_grp, by=quo_name(grp_var)) %>%
    group_by(!!grp_var)
}
