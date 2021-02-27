## ---------------------------
##
## Script name: groupwise_resampling.R
##
## Purpose of script: Create function for groupwise resampling
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

## load up the packages we will need:

if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(
  "dplyr",
  "purrr",
  "forcats",
  "rlang"
)

## ---------------------------

setwd(here()) # Project working directory

## ---------------------------

resample_groupwise <- function(.data, grp_var, n, prop, times = 1, weight_by = NULL, replace = FALSE) {

  n_int_digits = function(x) {
    result = floor(log10(abs(x))) + 1
    result[!is.finite(result)] = 0
    result
  }

  if (missing(grp_var)) rlang::abort("`grp_var` must be a valid grouping variable")
  grp_var <- enquo(grp_var)
  size <- dplyr:::check_slice_size(n, prop, "slice_sample")
  idx <- switch(
    size$type,
    n = function(x, n) dplyr:::sample_int(n, size$n, replace = replace, wt = x),
    prop = function(x, n) dplyr:::sample_int(n, size$prop * n, replace = replace, wt = x),
  )
  .data <- group_by(.data, !!grp_var)
  .data <- replicate(
    times,
    slice(
      .data,
      idx({
        {
          weight_by
        }
      }, dplyr::n()),
      .preserve = TRUE
    ),
    simplify = FALSE
  )
  .data <- map_dfr(
    .data,
    modify_at,
    vars(!!grp_var),
    fct_drop,
    .id = "Resample"
  )
  .data <- group_by(.data, Resample)
  .data <- nest(.data)
  .data <- mutate(.data, data = map(data, ~ group_by(.x, !!grp_var)))
  ungroup(.data)
}
