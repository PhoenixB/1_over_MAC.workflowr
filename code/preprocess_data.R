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

pacman::p_load("profvis", "htmlwidgets")

preprocess_profile <- profvis({

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
conflict_prefer("filter", "dplyr", quiet = TRUE) # Set preferred functions of conflicting packages
conflict_prefer("discard", "purrr", quiet = TRUE)
conflict_prefer("col_factor", "readr", quiet = TRUE)
conflict_prefer(":=", "data.table", quiet = TRUE)

## ---------------------------

## load up our functions into memory

source("code/group_sampling.R")

rangediff <- function(x) {
  range(x)[2] - range(x)[1]
}

`%!in%` <- Negate(`%in%`)

scale_vec <- function(x, scale = TRUE, attributes = TRUE) {
  scaled_x <- scale(x, scale = scale)
  scaled_x_vec <- as.vector(scaled_x)
  if (attributes) attributes(scaled_x_vec) <- attributes(scaled_x)
  scaled_x_vec
  }

## ---------------------------

# Create helper for column names
col_names <- c(
  "time",
  "ce_mac",
  "peak_alpha",
  "osc_alpha_power",
  "osc_alpha"
)

# Set path of data files relative to working directory
data_path = path_wd(
  "data"
)

# Directory listing of rejection files
rejection_files <-
  dir_ls(path = data_path, regexp = "to_reject_from_analysis_.*?\\.xlsx$")

# Read rejection files from directory listing
patients_to_reject <-
  rejection_files %>%
  map(~ read_excel(
    path = .x,
    col_names = c("pid"),
    na = c("", "NA", "NaN"),
    skip = 1)
  ) %>%
  map(~ setDT(.x)) %>%
  rbindlist(idcol = "file_name")                           # Bind lists together

# Convert `pid` to factor
patients_to_reject[, pid := lapply(.SD, as_factor), .SDcols = "pid"]

# Create file date column
patients_to_reject[,
  file_date := dmy(str_replace(
    file_name,
    glue(
      ".*?_(((((0?[1-9])|(1\\d)|(2[0-8]))\\.((0?[1-9])|(1[0-2])))|((31\\.((0",
      "[13578])|(1[02])))|((29|30)\\.((0?[1,3-9])|(1[0-2])))))\\.((20[0-9]",
      "[0-9]))|(29\\.0?2\\.20(([02468][048])|([13579][26]))))\\.xlsx$"
      ),
    "\\1"
    ))]

# Remove file name column
patients_to_reject[, file_name := NULL]

# Set key on `pid`
setkey(patients_to_reject, pid)

# Directory listing of data files
data_files <-
  dir_ls(
    path = data_path,
    regexp = glue(
      "alpha_freq_(?:Jan(?:uary)?|Feb(?:ruary)?|Mar(?:ch)?|Apr(?:il)?|May",
      "|Jun(?:e)?|Jul(?:y)?|Aug(?:ust)?|Sep(?:t)?(?:ember)?|Oct(?:ober)?|",
      "Nov(?:ember)?|Dec(?:ember)?)(?:19[7-9]\\d|2\\d{{3}})_\\d*?\\.csv$"
    )
  )

# Read data files into data.table list
one_over_mac <-
  data_files %>%
  map(~ fread(
    file = .x,
    col.names = col_names,
    na.strings = c("", "NA", "NaN")
  )) %>%
  discard( ~ nrow(.x) == 0) %>%         # Reject data.tables with 0 observations
  rbindlist(idcol = "file_name")                            # Bind list together

# Set `osc_alpha` to logical
one_over_mac[,
  osc_alpha := lapply(.SD, as.logical), .SDcols = "osc_alpha"]

# Create `pid` column as factor
one_over_mac[,
  pid := fct_inseq(factor(str_replace(file_name, ".*?_([0-9]+)\\.csv$", "\\1")))
  ]

# Remove file name column
one_over_mac[,
  file_name := NULL
]

# Set keys and sort by `pid` and `time`
setkey(one_over_mac, pid, time)

# Directory listing of age files
age_files <-
  dir_ls(path = data_path, regexp = "alpha_freq_ages_.*?\\.csv$")

# Read age files into data.table list
patients_ages <-
  age_files %>%
  map(~ fread(
    file = .x,
    col.names = c("pid", "age"),
    colClasses = c("factor", "integer"),
    na.strings = c("", "NA", "NaN")
  )) %>%
  rbindlist(idcol = "file_name")                            # Bind list together

# Reorder `pid` levels by numeric value
patients_ages[, pid := lapply(.SD, fct_inseq), .SDcols = "pid"]

# Create file date column
patients_ages[,
  file_date := my(str_replace(
    file_name,
    glue(
      ".*?_(Jan(?:uary)?|Feb(?:ruary)?|Mar(?:ch)?|Apr(?:il)?|May|Jun(?:e)?|J",
      "ul(?:y)?|Aug(?:ust)?|Sep(?:tember)?|Oct(?:ober)?|Nov(?:ember)?|Dec(?:",
      "ember)?)(19[7-9]\\d|2\\d{{3}})\\.csv$"),
    "\\1 \\2"))
  ]

# Remove file name column
patients_ages[, file_name := NULL]

# Set key on `pid`
setkey(patients_ages, pid)

# Full outer join age and data tables
one_over_mac <-
  merge(
    patients_ages[
      file_date == max(file_date),                    # Use only newest age data
      .(pid, age)
      ],
    one_over_mac, by = c("pid"), all.y = TRUE)

# Reset keys and sort by `pid` and `time` after join
setkey(one_over_mac, pid, time)

# Calculate age-adjusted `ce_mac`
one_over_mac[, aa_ce_mac := ce_mac / 10 ^ (-0.00269 * (age - 40))]

# # Create rolling means for `ce_mac`, `aa_ce_mac` and `peak_alpha`
# one_over_mac[,
#              paste0(rep(c("ce_mac", "aa_ce_mac", "peak_alpha"), each = 2),
#                     rep(c("_mean_5pt", "_mean_15pt"), times = 2)) :=
#                frollmean(.SD, n = c(5, 15), align = "center", na.rm = FALSE,
#                          hasNA = TRUE),
#              by = pid,
#              .SDcols = c("ce_mac", "aa_ce_mac", "peak_alpha")]

# Mark patients from rejection list
one_over_mac[,
  rejected := pid %in% patients_to_reject[, pid]
]

# Calculate modified z-score from `peak_alpha` by `pid`
one_over_mac[,
  peak_alpha_mz := (0.6745 *
    (peak_alpha - median(peak_alpha, na.rm = TRUE))) /
    (mad(peak_alpha, na.rm = TRUE)),
  by = pid
  ]

# Mark absolute z-score above 3.5 as outlier
one_over_mac[,
  peak_alpha_outlier := as.logical(abs(peak_alpha_mz) > 3.5)
  ]

ce_mac_upper_boundary <- 1.5
ce_mac_lower_boundary <- 0.5

# Mark `ce_mac` values above 1.5 and below 0.5, which is clinically sensible
one_over_mac[,
  ce_mac_outlier := (ce_mac < ce_mac_lower_boundary |
                       ce_mac > ce_mac_upper_boundary)
  ]

# Create helper for the median of the range of `ce_mac`
helper.ce_mac_range_median <-
  one_over_mac[,
    .(range = rangediff(ce_mac)),
    by = pid
    ][,
      median(range)
      ]

# Create range column for `ce_mac` by patient
one_over_mac[,
  ce_mac_range := rangediff(ce_mac),
  by = pid
  ]

# Filter patients with lower-than-median `ce_mac` range
one_over_mac[,
  range_outlier := ce_mac_range <= helper.ce_mac_range_median
  ]

# Create filtered dataset
one_over_mac_filtered <-
  na.omit(
    one_over_mac[
      ce_mac_outlier != TRUE & rejected != TRUE &
        peak_alpha_outlier != TRUE & range_outlier != TRUE
    ]
  )

# Set keys and sort by `pid` and `time`
setkey(one_over_mac_filtered, pid, time)

# Drop unused factor levels
one_over_mac_filtered[,
                      pid := lapply(.SD, fct_drop), .SDcols = c("pid")
]

# Create helper for centering columns
center_cols <-
  c("ce_mac", "aa_ce_mac")

# Create centered column for `ce_mac` and `aa_ce_mac`
one_over_mac[,
  paste0(center_cols, "_centered") :=
    lapply(.SD, scale_vec, scale = FALSE, attributes = FALSE),
  .SDcols = center_cols
]
one_over_mac_filtered[,
             paste0(center_cols, "_centered") :=
               lapply(.SD, scale_vec, scale = FALSE, attributes = FALSE),
             .SDcols = center_cols
]

# Create within-group centered column for `ce_mac` and `aa_ce_mac`
one_over_mac[,
  paste0(center_cols, "_pid_centered") :=
    lapply(.SD, scale_vec, scale = FALSE, attributes = FALSE),
  by = pid,
  .SDcols = center_cols
]
one_over_mac_filtered[,
  paste0(center_cols, "_pid_centered") :=
    lapply(.SD, scale_vec, scale = FALSE, attributes = FALSE),
  by = pid,
  .SDcols = center_cols
]

# # Create scaled column for `ce_mac`
# one_over_mac[,
#   ce_mac_scaled := scale(ce_mac)
# ]

# # Create within-group scaled column for `ce_mac`
# one_over_mac[,
#   ce_mac_pid_scaled := scale(ce_mac),
#   by = pid
# ]

# Create mean column for `ce_mac` and `aa_ce_mac`
one_over_mac[,
  paste0(center_cols, "_mean") := lapply(.SD, mean),
  .SDcols = center_cols
]
one_over_mac_filtered[,
  paste0(center_cols, "_mean") := lapply(.SD, mean),
  .SDcols = center_cols
]

# Create within-group mean column for `ce_mac` and `aa_ce_mac`
one_over_mac[,
  paste0(center_cols, "_pid_mean") := lapply(.SD, mean),
  by = pid,
  .SDcols = center_cols
]
one_over_mac_filtered[,
  paste0(center_cols, "_pid_mean") := lapply(.SD, mean),
  by = pid,
  .SDcols = center_cols
]

# Save processed data to disk
write_rds(one_over_mac, file = path_wd("output", "one_over_mac", ext = "rds"))
write_rds(one_over_mac_filtered,
          file = path_wd("output", "one_over_mac_filtered", ext = "rds"))

})

saveWidget(preprocess_profile, file = path_wd("output", "preprocess_profile", ext = "html"))
