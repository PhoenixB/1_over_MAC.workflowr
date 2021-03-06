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
  "janitor",
  "data.table",
  "glue"
)

## ---------------------------

setwd(here())                                        # Project working directory

## ---------------------------

set.seed(12345)                        # Set pseudo-random number generator seed
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

col_names <- c(                                 # Create helper for column names
  "time",
  "ce_mac",
  "peak_alpha",
  "osc_alpha_power",
  "osc_alpha"
)

data_path = path_wd(      # Set path of data files relative to working directory
  "data"
)

rejection_files <-                        # Directory listing of rejection files
  dir_ls(path = data_path, regexp = "to_reject_from_analysis_.*?\\.xlsx$")

patients_to_reject <-
  rejection_files %>%              # Read rejection files from directory listing
  map(~ read_excel(
    path = .x,
    col_names = c("pid"),
    na = c("", "NA", "NaN"),
    skip = 1)
  ) %>%
  bind_rows(.id = "file_name") %>%                         # Bind lists together
  mutate(file_date = dmy(                   # Set date column of rejection files
    str_replace(
      file_name,
      glue(
        ".*?_(((((0?[1-9])|(1\\d)|(2[0-8]))\\.((0?[1-9])|(1[0-2])))|((31\\.((0",
        "[13578])|(1[02])))|((29|30)\\.((0?[1,3-9])|(1[0-2])))))\\.((20[0-9]",
        "[0-9]))|(29\\.0?2\\.20(([02468][048])|([13579][26]))))\\.xlsx$"),
      "\\1"
    )
  )) %>%
  select(file_date, pid) %>%
  modify_at(vars(pid), as_factor) %>%
  filter(file_date == max(file_date))            # Filter only newest rejections

data_files <-                                  # Directory listing of data files
  dir_ls(path = data_path, regexp = "\\.csv$")

one_over_mac <-
  data_files %>%
  map(~ fread(                            # Read data files into data.table list
    file = .x,
    col.names = col_names,
    na.strings = c("", "NA", "NaN")
  )) %>%
  discard( ~ nrow(.x) == 0) %>%         # Reject data.tables with 0 observations
  rbindlist(idcol = "file_name")                            # Bind list together

one_over_mac[,                                        # Set osc_alpha to logical
                osc_alpha := lapply(.SD, as.logical), .SDcols = "osc_alpha"]

one_over_mac[,                                     # Create pid column as factor
  pid := fct_inseq(factor(str_replace(file_name, ".*?_([0-9]+)\\.csv$", "\\1")))
  ]

one_over_mac[,                                         # Remove file name column
  file_name := NULL
]

setkey(one_over_mac, pid, time)              # Set keys and sort by pid and time

# one_over_mac[,                  # Create rolling means for ce_mac and peak_alpha
#                 paste0(c("ce_mac", "peak_alpha"), "_mean_5pt") :=
#                   lapply(
#                     .SD,
#                     frollmean,
#                     n = 5,
#                     align = "center",
#                     na.rm = FALSE,
#                     hasNA = TRUE
#                   ),
#                 by = pid,
#                 .SDcols = c("ce_mac", "peak_alpha")]

one_over_mac[,                               # Mark patients from rejection list
  rejected := pid %in% pluck(patients_to_reject, "pid")
]

one_over_mac[,               # Calculate modified z-score from peak_alpha by pid
  peak_alpha_mz := (0.6745 *
    (peak_alpha - median(peak_alpha, na.rm = TRUE))) /
    (mad(peak_alpha, na.rm = TRUE)),
  by = pid
  ]

one_over_mac[,                      # Mark absolute z-score above 3.5 as outlier
  peak_alpha_outlier := as.logical(abs(peak_alpha_mz) > 3.5)
  ]


helper.upper_boundary <- 1.5
helper.lower_boundary <- 0.5

# Mark ce_mac values above 1.5 and below 0.5, which is clinically sensible
one_over_mac[,
  ce_mac_outlier := (ce_mac < helper.lower_boundary |
                       ce_mac > helper.upper_boundary)
  ]

# Create helper for the median of the range of ce_mac
helper.ce_mac_range_median <-
  one_over_mac[,
    .(range = rangediff(ce_mac)),
    by = pid
    ][,
      median(range)
      ]

one_over_mac[,                       # Create range column for ce_mac by patient
  ce_mac_range := rangediff(ce_mac),
  by = pid
  ]

one_over_mac[,             # filter patients with lower-than-median ce_mac range
  range_outlier := ce_mac_range <= helper.ce_mac_range_median
  ]

one_over_mac_filtered <-                               # Create filtered dataset
  na.omit(
    one_over_mac[
      ce_mac_outlier != TRUE & rejected != TRUE &
        peak_alpha_outlier != TRUE & range_outlier != TRUE
    ]
  )

# Set class to tibble
setattr(one_over_mac, "class", c("tbl", "tbl_df", "data.frame"))
setattr(one_over_mac_filtered, "class", c("tbl", "tbl_df", "data.frame"))

# Save processed data to disk
write_rds(one_over_mac, file = path_wd("output", "one_over_mac", ext = "rds"))
write_rds(one_over_mac_filtered,
          file = path_wd("output", "one_over_mac_filtered", ext = "rds"))
