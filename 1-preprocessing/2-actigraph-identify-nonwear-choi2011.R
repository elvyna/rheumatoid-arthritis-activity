################################################################################################################################################
## 00. LOAD LIBRARIES AND SET WORKING DIRECTORY
################################################################################################################################################

library(zoo)

DIRECTORY <- '/PUT-THE-WORKING-DIRECTORY/'
setwd(DIRECTORY)

################################################################################################################################################
## 01. DEFINE FUNCTIONS
## 01-1. IDENTIFY NONWEAR PERIOD
## Reference:
## Choi, L., Liu, Z., Matthews, C. E., & Buchowski, M. S. (2011). 
## Validation of accelerometer wear and nonwear time classification algorithm. 
## Medicine and science in sports and exercise, 43(2), 357–364. https://doi.org/10.1249/MSS.0b013e3181ed61a3
################################################################################################################################################

## input: 1 dataframe, 2 numeric vectors: start and end indices of consecutive zeroes
## for each pair, slice the dataframe indices, set is_nonwear to TRUE
set_nonwear_period_full_window <- function(series, start_idx, end_idx) {
  if (length(start_idx) != length(end_idx)) {
    stop("Start and end indices do not have similar length!")
  } else if ((length(start_idx) == 0) || (is.na(start_idx))) {
    stop("Index is empty")
  } else {
    idx_pair_sequence <- 1:length(start_idx)
    is_nonwear <- rep(FALSE, length.out = length(series))
    for (i in idx_pair_sequence) {
      mask_nonwear_idx <- (start_idx[i] : end_idx[i])
      is_nonwear[mask_nonwear_idx] <- TRUE
    }
  }
  is_nonwear
}

## input: 
## 1. series of is_nonwear data
## 2. gap: number of consecutive zeroes
## 3. window length of consecutive zeroes, default = 90
set_nonwear_period_with_break_allowance <- function(series, 
                                                    nonzero_idx,
                                                    gap, 
                                                    consecutive_zeroes_window = 90,
                                                    break_allowance_duration = 2) {
  ## Use zoo::rollsum
  ### ASSUMPTION: total consecutive zeroes must be >= 88 (90 min - 2 artifactual movement interval)
  ### then, check the downstream & upstream, i.e., it must have >= 30 consecutive zeroes
  ### if true, then the period will be considered as nonwear
  MINIMUM_CONSECUTIVE_ZEROES <- consecutive_zeroes_window
  ARTIFICIAL_MOVEMENT_ALLOWANCE <- break_allowance_duration
  ## index of the consecutive zeroes with break (but is it 2 minute break?) TO DO: CHECK
  rollsum_geq_90min <- which(rollsum(gap, 3) >= MINIMUM_CONSECUTIVE_ZEROES - ARTIFICIAL_MOVEMENT_ALLOWANCE)
  ## based on the 3-window rolling sum, 
  ## we need to find the 3-window which have 0 consecutive zero in the middle (that's why we add +1) 
  ## then we add 1 to convert the index back to original value (due to diff)
  break_allowance_2min <- which(gap[(rollsum_geq_90min + 1)] == 0)  
  break_allowance_2min_idx <- rollsum_geq_90min[break_allowance_2min] + 1
  ## get the index of the break in the actual series (axis1) (note: the break is up to 2 minutes)
  allowance_2min_idx <- nonzero_idx[break_allowance_2min_idx]
  
  ## check the downstream 30 minutes before the break
  DOWNSTREAM_PERIOD <- 30
  UPSTREAM_PERIOD <- 30
  # print(rollsum_geq_90min)
  # print(break_allowance_2min)
  # print(break_allowance_2min_idx)
  # print(allowance_2min_idx)
  # 
  ## should add 1 minute to search upstream (since we allow 2 minute breaks, and current index is the first break)
  for (i in seq_along(allowance_2min_idx)) {
    downstream_start_idx <- (allowance_2min_idx[i] - DOWNSTREAM_PERIOD)
    upstream_end_idx <- (allowance_2min_idx[i] + UPSTREAM_PERIOD + 1)
    check_downstream_consecutive_zeroes <- sum(series[downstream_start_idx : allowance_2min_idx[i]] == 0)
    check_upstream_consecutive_zeroes <- sum(series[allowance_2min_idx[i] : upstream_end_idx] == 0)
    
    ## If both of them equals to DOWNSTREAM_PERIOD and UPSTREAM_PERIOD, 
    ## then set is_nonwear to TRUE during the consecutive zeroes + 2 minute allowance
    if ((check_downstream_consecutive_zeroes >= DOWNSTREAM_PERIOD) & (check_upstream_consecutive_zeroes >= UPSTREAM_PERIOD)) {
      ## Get the start of downstream period & end of upstream period, set those period as nonwear
      downstream_zeroes_idx <- rollsum_geq_90min[break_allowance_2min[i]]
      ## Start of downstream period. add 1 to convert it back to the original axis
      downstream_zeroes_actual_idx <- nonzero_idx[downstream_zeroes_idx] + 1 
      
      upstream_zeroes_idx <- rollsum_geq_90min[break_allowance_2min[i]] + 2
      ## start of upstream period (may include 1 minute with nonzero CPM)
      upstream_zeroes_actual_idx <- nonzero_idx[upstream_zeroes_idx] 
      ## Get the end of upstream period (based on length)
      upstream_length <- gap[upstream_zeroes_idx]
      if (length(upstream_zeroes_actual_idx) == length(upstream_length)) {
        upstream_end_period_idx <- upstream_zeroes_actual_idx + upstream_length
      }
      
      ## Determine is_nonwear = TRUE
      series[downstream_zeroes_actual_idx:upstream_end_period_idx] <- TRUE
    }
  }
  series
}

## Nonwear period identification based on the proposed algorithm from Choi et al. (2011).
## Choi, L., Liu, Z., Matthews, C. E., & Buchowski, M. S. (2011). 
## Validation of accelerometer wear and nonwear time classification algorithm. 
## Medicine and science in sports and exercise, 43(2), 357–364. 
## https://doi.org/10.1249/MSS.0b013e3181ed61a3
detect_nonwear <- function(series,
                           consecutive_zeroes_window = 90,
                           break_allowance_duration = 2) {
  ## Get index of nonzero cpm
  nonzero_idx <- which(series != 0)
  if (is.na(nonzero_idx) || length(nonzero_idx) == 0) {
    return(0)
  }
  
  ## Set minimum consecutive zeroes count to be considered as nonwear (in minute)
  MINIMUM_CONSECUTIVE_ZEROES <- consecutive_zeroes_window 
  ARTIFICIAL_MOVEMENT_ALLOWANCE <- break_allowance_duration
  ## Calculate the length of consecutive zeroes. if 0, it means no gap (no zero cpm)
  gap <- (diff(nonzero_idx) - 1)
  ## Get indices where we find the END of long consecutive zeroes (90 min window):
  ## Add 1 to get the index of gap, 
  ## Then minus 1 on nonzero_idx to get the last consecutive zero 0
  gap_idx_geq_90 <- which(gap >= MINIMUM_CONSECUTIVE_ZEROES) + 1
  geq90_consecutive_zeroes_idx_end <- nonzero_idx[gap_idx_geq_90] - 1
  ## Get indices where we find the START of long consecutive zeroes (the first zero)
  ## Add 1 to get the start of zero (not the index of nonzero value before consecutive zeroes)
  geq90_consecutive_zeroes_idx_start <- geq90_consecutive_zeroes_idx_end - gap[gap >= MINIMUM_CONSECUTIVE_ZEROES] + 1
  ## Determine wear / nonwear period based on 90 minute consecutive windows (no break)
  if (is.na(gap_idx_geq_90) || length(gap_idx_geq_90) == 0) {
    return(0)
  } else {
    is_nonwear <- set_nonwear_period_full_window(series, geq90_consecutive_zeroes_idx_start, geq90_consecutive_zeroes_idx_end)
  }
  ## Determine wear / nonwear period based on 90 minute consecutive windows with 
  is_nonwear <- set_nonwear_period_with_break_allowance(is_nonwear, 
                                                        nonzero_idx, 
                                                        gap, 
                                                        MINIMUM_CONSECUTIVE_ZEROES, 
                                                        ARTIFICIAL_MOVEMENT_ALLOWANCE)
  is_nonwear
}

# ## test in one patient
# df <- read.csv("INDIVIDUAL_PATIENT_DATA.csv")
# df$is_nonwear <- detect_nonwear(df$axis1)
# nrow(df[df$is_nonwear,])

################################################################################################################################################
## 01. DEFINE FUNCTIONS
## 01-2. CLASSIFY ACTIVITY TYPE
################################################################################################################################################

## Classify activity type.
## Input: series of vertical axis movement (axis1)
## Reference:
## Meiring, R. M., Frimpong, E., Mokete, L., Pietrzak, J., Van Der Jagt, D., Tikly, M., & McVeigh, J. A. (2016). 
## Rationale, design and protocol of a longitudinal study assessing the effect of total knee arthroplasty 
## on habitual physical activity and sedentary behavior in adults with osteoarthritis. 
## BMC musculoskeletal disorders, 17(1), 281.
determine_activity_type <- function(series) {
  SEDENTARY_THRESHOLD <- 100
  LIGHT_THRESHOLD <- 1952
  MODERATE_THRESHOLD <- 5725
  ifelse(
    series < SEDENTARY_THRESHOLD, 
    'sedentary_behaviour',
    ifelse(
      (series >= SEDENTARY_THRESHOLD) & (series < LIGHT_THRESHOLD),
      'light_intensity',
      ifelse(
        (series >= LIGHT_THRESHOLD) & (series < MODERATE_THRESHOLD),
        'moderate',
        ifelse(
          (series >= MODERATE_THRESHOLD),
          'vigorous',
          NA
        )
      ) 
    )
  )
}

################################################################################################################################################
## 02. DATA IMPORT AND PREPROCESSING
################################################################################################################################################

## Identify non-wear in all patients data
ACTIGRAPH_FILE <- "../dataset/actigraph_combined_with_meta.csv"
df_actigraph <- read.csv(ACTIGRAPH_FILE)
head(df_actigraph)

## Split data frame per patient-period group, then run the function
pid_period_groups <- split(df_actigraph, list(df_actigraph$patient_id, df_actigraph$directory_name))

is_nonwear_output <- sapply(
  pid_period_groups, 
  function(x) { 
    detect_nonwear(x$axis1)
    }
  )

idx <- names(is_nonwear_output)
for (i in idx) {
  pid_period_groups[[i]][,'is_nonwear'] <- is_nonwear_output[i]
}

## Combine the lists back to a data frame
df_actigraph <- do.call(rbind, pid_period_groups)

## Classify activity type based on vertical axis movement
df_actigraph$activity_type <- determine_activity_type(df_actigraph$axis1)

## Determine final activity classification
## If is_nonwear is true, change activity_type to non_wear
df_actigraph$behavior_final <- ifelse(
  df_actigraph$is_nonwear == 1,
  'nonwear',
  df_actigraph$activity_type
)

################################################################################################################################################
## 03. ASSESS RESULTS
################################################################################################################################################

## check number of nonwear occurrence per patient
with(
  df_actigraph, 
  tapply(is_nonwear, list(patient_id, directory_name), sum)
)

## check percentage of nonwear per patient
with(
  df_actigraph, 
  tapply(is_nonwear, list(patient_id, directory_name), sum)
) / 
  with(
    df_actigraph, 
    tapply(is_nonwear, list(patient_id, directory_name), length)
  )

## check percentage of nonwear by day of week & hour
df_actigraph$day_of_week_string <- factor(
  df_actigraph$day_of_week_string, 
  levels = c('Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday')
  )

actigraph_nonwear_time_of_day <- with(
  df_actigraph, 
  tapply(is_nonwear, list(day_of_week_string, hour), sum)
) / 
  with(
    df_actigraph, 
    tapply(is_nonwear, list(day_of_week_string, hour), length)
  )

## Show as percentage
t(actigraph_nonwear_time_of_day) * 100

## assess results in one patient & obs period
PID_INPUT <- 2
sub <- df_actigraph[
  (df_actigraph$patient_id == PID_INPUT) &
    (df_actigraph$directory_name == 'Baseline'), 
  ]
with(sub, tapply(is_nonwear, list(day_of_week_string, hour), sum))

## Check percentage of nonwear based on classified activity
## As expected, most nonwear are classified as SB or LIPA
with(
  df_actigraph, 
  tapply(is_nonwear, list(patient_id, activity_type), sum)
) / 
  with(
    df_actigraph, 
    tapply(is_nonwear, list(patient_id, activity_type), length)
  )

## Check final behavior classification after is_nonwear overrides classified activity type
with(df_actigraph, table(day_of_week, behavior_final))
with(df_actigraph, table(hour, behavior_final))

# write.csv(
#   df_actigraph,
#   paste0(DIRECTORY, "../dataset/actigraph_preprocessed_choi.csv"),
#   row.names = FALSE
#   )