################################################################################################################################################
## 00. LOAD LIBRARIES
## ALSO SET WORKING DIRECTORY
################################################################################################################################################

library(dplyr)
library(ggthemes) # viz
library(ggplot2) # viz
theme_set(theme_economist())

## set working directory  
DIRECTORY <- '/PUT-THE-WORKING-DIRECTORY/'
setwd(DIRECTORY)

################################################################################################################################################
## 01. DATA IMPORT AND PREPARATION
## 01-1. CREATE HOURLY DATA
################################################################################################################################################

## Import combined data after using SAS code from Winkler et al. (2016)
ACTIVPAL_FILE <- "../dataset/activpal_combined_processed_in_sas_with_meta.csv"
df_activpal <- read.csv(ACTIVPAL_FILE)

## Identify invalid observations
## Based on preprocessing results
df_activpal <- df_activpal %>% 
  mutate(
    is_invalid = case_when(
      (
        (include_activpal_data == 1) &
          (validday == 1) &
          ((AW_steps_n + SL_steps_n) == steps_n) &
          (sleepboutall == 0) 
      ) ~ 0,
      TRUE ~ 1
    )
  )

df_activpal$day_of_week_string <- weekdays(as.Date(df_activpal$date))
df_activpal$day_of_week_string <- factor(
  df_activpal$day_of_week_string, 
  levels = c('Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday')
)

## Observe % of invalid data based on patients who complete both the baseline and follow-up period
## Only take patients with valid data in both observation periods
activpal_patient_obs_check <- df_activpal %>% 
  filter(
    (visit_info %in% c('Baseline','Follow up'))
  ) %>% 
  group_by(patient_id) %>% 
  summarise(
    observation_count = length(unique(visit_info))
  ) %>% 
  filter(observation_count == 2)
ACTIVPAL_PATIENT_ID_BOTH_PERIOD <- activpal_patient_obs_check$patient_id

## Get observations of patients who participated in both period ## with valid data in both period
sub_activpal <- df_activpal %>% 
  filter(patient_id %in% ACTIVPAL_PATIENT_ID_BOTH_PERIOD)

## Based on day and hour
activpal_invalid_time_of_day <- with(
  sub_activpal, 
  tapply(is_invalid, list(day_of_week_string, hour), sum)
) / 
  with(
    sub_activpal, 
    tapply(is_invalid, list(day_of_week_string, hour), length)
  )

df_activpal_invalid <- data.frame(
  t(activpal_invalid_time_of_day), 
  "hour" = row.names(t(activpal_invalid_time_of_day)), 
  "device" = "activpal"
)

## Based on day and period
activpal_invalid_day_period <- with(
  sub_activpal, 
  tapply(is_invalid, list(day_of_week_string, visit_info), sum)
) / 
  with(
    sub_activpal, 
    tapply(is_invalid, list(day_of_week_string, visit_info), length)
  )

df_activpal_invalid <- data.frame(
  t(activpal_invalid_day_period), 
  "period" = row.names(t(activpal_invalid_day_period)), 
  "device" = "activpal"
)

################################################################################################################################################
## 01. DATA IMPORT AND PREPARATION
## 01-2. VISUALISE PROPORTION OF INVALID DATA -- COMPARE WITH ACTIGRAPH PREPROCESSING RESULT
## RUN 3-actigraph-identify-nonwear-choi2011-sanity-check.R FIRST
################################################################################################################################################

## Only use this to compare % of invalid observation vs. actigraph
## Run 2020-09-05-actigraph-identify-nonwear-choi2011-sanity-check.R first
df_invalid <- rbind(df_actigraph_invalid, df_activpal_invalid)
rm(df_actigraph_invalid, df_activpal_invalid)

## Create a dataframe with long format
df_invalid_long <- df_invalid %>% 
  reshape2::melt(
    measure.vars = c("Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday"),
    id.vars = c("device","period"),
    variable.name = "day_of_week",
    value.name = "proportion_invalid_obs"
  )

## Visualise
ggplot(data = df_invalid_long) +
  geom_bar(
    aes(
      x = day_of_week,
      y = proportion_invalid_obs,
      fill = device
    ),
    position = "dodge",
    stat = "identity"
  ) +
  facet_wrap(. ~ period) +
  scale_fill_colorblind(
    name = ""
  ) +
  scale_x_discrete(
    name = "day of week"
  ) +
  scale_y_continuous(
    name = 'proportion of invalid records',
    breaks = seq(0,1,.2),
    labels = scales::percent
  ) +
  coord_flip() +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)
  )

################################################################################################################################################
## 01. DATA IMPORT AND PREPARATION
## 01-3. REMOVE INVALID DATA
################################################################################################################################################

## Only use "valid" observations
## Remove sleep bouts
df_activpal <- df_activpal %>% 
  filter(is_invalid == 0)

agg <- df_activpal %>% 
  group_by(
    patient_id, 
    visit_info, 
    datetime_hour
    ) %>% 
  summarise(
    sleep_total = sum(SL_all_t),
    awake_total = sum(AW_all_t),
    awake_sitting_total = sum(AW_sitting_t),
    awake_standing_total = sum(AW_standing_t),
    awake_stepping_total = sum(AW_stepping_t),
    awake_steps_total = sum(AW_steps_n),
    duration_total = sum(SL_all_t + AW_all_t),
    sleep_proportion = sum(SL_all_t) / sum(SL_all_t + AW_all_t)
    )

################################################################################################################################################
## 01. DATA IMPORT AND PREPARATION
## 01-4. ATTRIBUTE DURATION BASED ON EACH HOUR, NOT ON THE START OF THE EVENT
################################################################################################################################################

## Generate full record of patient_id per hour, based on their observation timeframe
patient_timeframe <- agg %>% 
  group_by(patient_id, visit_info) %>% 
  summarise(
    min_datetime = min(as.Date(datetime_hour)), 
    max_datetime = max(as.Date(datetime_hour))
    )

## Create empty dataframe to store the populated time sequence per patient ID
df_patient_obs <- data.frame()

for (pid in unique(patient_timeframe$patient_id)){
  patient_sub <- patient_timeframe %>% filter(patient_id == pid)
  for (period in patient_sub$visit_info){
    patient_sub_period <- patient_sub %>% filter(visit_info == period)
    tmp_seq <- seq(
      as.POSIXct(patient_sub_period$min_datetime), 
      as.POSIXct(patient_sub_period$max_datetime),
      'hour')
    tmp_df <- data.frame(
      patient_id = pid,
      visit_info = period,
      datetime_hour = tmp_seq
    )
    df_patient_obs <- rbind(
      df_patient_obs,
      tmp_df
    )
  }
}

agg$datetime_hour <- as.POSIXct(agg$datetime_hour)
df_patient_obs <- df_patient_obs %>% 
  left_join(
    agg %>% 
      select(
        patient_id, 
        visit_info, 
        datetime_hour,
        sleep_total, 
        awake_total, 
        awake_steps_total,
        awake_sitting_total, 
        awake_standing_total, 
        awake_stepping_total
        ),
    by = c(
      'patient_id' = 'patient_id',
      'visit_info' = 'visit_info',
      'datetime_hour' = 'datetime_hour'
      )
    )

## Remove rows with full NAs
df_patient_obs <- df_patient_obs[!is.na(df_patient_obs$patient_id), ]

## generate hourly duration of the awake and sleep time, attribute to the corresponding hour
df_patient_obs <- df_patient_obs %>% 
  mutate(
    sleep_total_hour = sleep_total / 3600,
    sleep_modulo_hour_second = sleep_total %% 3600,
    awake_total_hour = awake_total / 3600,
    awake_modulo_hour_second = awake_total %% 3600,
    awake_sitting_total_hour = awake_sitting_total / 3600,
    awake_sitting_modulo_hour_second = awake_sitting_total %% 3600,
    awake_standing_total_hour = awake_standing_total / 3600,
    awake_standing_modulo_hour_second = awake_standing_total %% 3600,
    awake_stepping_total_hour = awake_stepping_total / 3600,
    awake_stepping_modulo_hour_second = awake_stepping_total %% 3600
    )

## Fill out awake time per hour
for (idx in 1:nrow(df_patient_obs)){
  if ((df_patient_obs[idx, 'awake_total_hour'] > 1) && (!is.na(df_patient_obs[idx, 'awake_total_hour']))){
    awake_dur_hour <- floor(df_patient_obs[idx, 'awake_total_hour'])
    additional_awake_obs <- awake_dur_hour - 1
    
    if (additional_awake_obs >= 1){
      for (next_record in 1:additional_awake_obs) {
        df_patient_obs[idx + next_record, 'awake_additional_hour'] <- 1
      }
    }
    ## add the modulo
    if (
      (additional_awake_obs < 1) &&
      (df_patient_obs[idx, 'awake_modulo_hour_second'] > 0) && 
      (!is.na(df_patient_obs[idx, 'awake_modulo_hour_second']))
    ) {
      df_patient_obs[idx + 1, 'awake_additional_second'] <- df_patient_obs[idx, 'awake_modulo_hour_second']
    }
    else {
      df_patient_obs[idx + next_record + 1, 'awake_additional_second'] <- df_patient_obs[idx, 'awake_modulo_hour_second']
    }
  }
}

## Aggregate the actual awake duration per hour
for (idx in 1:nrow(df_patient_obs)){
  ## replace NA as 0, so that we can calculate the aggregated value
  df_patient_obs[idx, ] <- replace(df_patient_obs[idx, ], is.na(df_patient_obs[idx, ]), 0)
  if ((df_patient_obs[idx, 'awake_total_hour'] >= 1) && (!is.na(df_patient_obs[idx, 'awake_total_hour']))) {
    df_patient_obs[idx, 'awake_duration_second_clean'] <- (1 + df_patient_obs[idx, 'awake_additional_hour']) * 3600 + df_patient_obs[idx, 'awake_additional_second']
  }
  else {
    df_patient_obs[idx, 'awake_duration_second_clean'] <- df_patient_obs[idx, 'awake_modulo_hour_second'] + (df_patient_obs[idx, 'awake_additional_hour'] * 3600) + df_patient_obs[idx, 'awake_additional_second']
  }
}

## Fill out awake_sitting time per hour
for (idx in 1:nrow(df_patient_obs)){
  if ((df_patient_obs[idx, 'awake_sitting_total_hour'] > 1) && (!is.na(df_patient_obs[idx, 'awake_sitting_total_hour']))){
    awake_sitting_dur_hour <- floor(df_patient_obs[idx, 'awake_sitting_total_hour'])
    additional_awake_sitting_obs <- awake_sitting_dur_hour - 1
    if (additional_awake_sitting_obs >= 1){
      for (next_record in 1:additional_awake_sitting_obs){
        df_patient_obs[idx + next_record, 'awake_sitting_additional_hour'] <- 1
      }
    }
    ## Add the modulo
    if (
      (additional_awake_sitting_obs < 1) &&
      (df_patient_obs[idx, 'awake_sitting_modulo_hour_second'] > 0) && 
      (!is.na(df_patient_obs[idx, 'awake_sitting_modulo_hour_second']))
    ) {
      df_patient_obs[idx + 1, 'awake_sitting_additional_second'] <- df_patient_obs[idx, 'awake_sitting_modulo_hour_second']
    }
    else {
      df_patient_obs[idx + next_record + 1, 'awake_sitting_additional_second'] <- df_patient_obs[idx, 'awake_sitting_modulo_hour_second']
    }
  }
}

## Aggregate the actual awake_sitting duration per hour
for (idx in 1:nrow(df_patient_obs)) {
  ## replace NA as 0, so that we can calculate the aggregated value
  df_patient_obs[idx, ] <- replace(df_patient_obs[idx, ], is.na(df_patient_obs[idx, ]), 0)
  if ((df_patient_obs[idx, 'awake_sitting_total_hour'] >= 1) && (!is.na(df_patient_obs[idx, 'awake_sitting_total_hour']))) {
    df_patient_obs[idx, 'awake_sitting_duration_second_clean'] <- (1 + df_patient_obs[idx, 'awake_sitting_additional_hour']) * 3600 + df_patient_obs[idx, 'awake_sitting_additional_second']
  }
  else {
    df_patient_obs[idx, 'awake_sitting_duration_second_clean'] <- df_patient_obs[idx, 'awake_sitting_modulo_hour_second'] + (df_patient_obs[idx, 'awake_sitting_additional_hour'] * 3600) + df_patient_obs[idx, 'awake_sitting_additional_second']
  }
}

## Standing and stepping time = awake.total - awake.sitting time
## Fill out awake.standing time per hour
df_patient_obs$awake_standing_additional_hour <- 0
for (idx in 1:nrow(df_patient_obs)) {
  if ((df_patient_obs[idx, 'awake_standing_total_hour'] > 1) && (!is.na(df_patient_obs[idx, 'awake_standing_total_hour']))){
    awake_standing_dur_hour <- floor(df_patient_obs[idx, 'awake_standing_total_hour'])
    additional_awake_standing_obs <- awake_standing_dur_hour - 1
    
    if (additional_awake_standing_obs >= 1){
      for (next_record in 1:additional_awake_standing_obs) {
        df_patient_obs[idx + next_record, 'awake_standing_additional_hour'] <- 1
      }
    }
    ## Add the modulo
    if (
      (additional_awake_standing_obs < 1) &&
      (df_patient_obs[idx, 'awake_standing_modulo_hour_second'] > 0) && 
      (!is.na(df_patient_obs[idx, 'awake_standing_modulo_hour_second']))
    ) {
      df_patient_obs[idx + 1, 'awake_standing_additional_second'] <- df_patient_obs[idx, 'awake_standing_modulo_hour_second']
    }
    else {
      df_patient_obs[idx + next_record + 1, 'awake_standing_additional_second'] <- df_patient_obs[idx, 'awake_standing_modulo_hour_second']
    }
  }
}

## Aggregate the actual awake.standing duration per hour
for (idx in 1:nrow(df_patient_obs)) {
  ## replace NA as 0, so that we can calculate the aggregated value
  df_patient_obs[idx, ] <- replace(df_patient_obs[idx, ], is.na(df_patient_obs[idx, ]), 0)
  if ((df_patient_obs[idx, 'awake_standing_total_hour'] >= 1) && (!is.na(df_patient_obs[idx, 'awake_standing_total_hour']))) {
    df_patient_obs[idx, 'awake_standing_duration_second_clean'] <- (1 + df_patient_obs[idx, 'awake_standing_additional_hour']) * 3600 + df_patient_obs[idx, 'awake_standing_additional_second']
  }
  else {
    df_patient_obs[idx, 'awake_standing_duration_second_clean'] <- df_patient_obs[idx, 'awake_standing_modulo_hour_second'] + (df_patient_obs[idx, 'awake_standing_additional_hour'] * 3600) + df_patient_obs[idx, 'awake_standing_additional_second']
  }
}

## Fill out awake_stepping time per hour
df_patient_obs$awake_stepping_additional_hour <- 0
df_patient_obs$awake_stepping_additional_second <- 0
for (idx in 1:nrow(df_patient_obs)) {
  if ((df_patient_obs[idx, 'awake_stepping_total_hour'] > 1) && (!is.na(df_patient_obs[idx, 'awake_stepping_total_hour']))) {
    awake_stepping_dur_hour <- floor(df_patient_obs[idx, 'awake_stepping_total_hour'])
    additional_awake_stepping_obs <- awake_stepping_dur_hour - 1
    
    if (additional_awake_stepping_obs >= 1) {
      for (next_record in 1:additional_awake_stepping_obs){
        df_patient_obs[idx + next_record, 'awake_stepping_additional_hour'] <- 1
      }
    }
    ## Add the modulo
    if (
      (additional_awake_stepping_obs < 1) &&
      (df_patient_obs[idx, 'awake_stepping_modulo_hour_second'] > 0) && 
      (!is.na(df_patient_obs[idx, 'awake_stepping_modulo_hour_second']))
    ) {
      df_patient_obs[idx + 1, 'awake_stepping_additional_second'] <- df_patient_obs[idx, 'awake_stepping_modulo_hour_second']
    }
    else {
      df_patient_obs[idx + next_record + 1, 'awake_stepping_additional_second'] <- df_patient_obs[idx, 'awake_stepping_modulo_hour_second']
    }
  }
}

## Aggregate the actual awake.stepping duration per hour
for (idx in 1:nrow(df_patient_obs)) {
  ## replace NA as 0, so that we can calculate the aggregated value
  df_patient_obs[idx, ] <- replace(df_patient_obs[idx, ], is.na(df_patient_obs[idx, ]), 0)
  if ((df_patient_obs[idx, 'awake_stepping_total_hour'] >= 1) && (!is.na(df_patient_obs[idx, 'awake_stepping_total_hour']))) {
    df_patient_obs[idx, 'awake_stepping_duration_second_clean'] <- (1 + df_patient_obs[idx, 'awake_stepping_additional_hour']) * 3600 + df_patient_obs[idx, 'awake_stepping_additional_second']
  }
  else {
    df_patient_obs[idx, 'awake_stepping_duration_second_clean'] <- df_patient_obs[idx, 'awake_stepping_modulo_hour_second'] + (df_patient_obs[idx, 'awake_stepping_additional_hour'] * 3600) + df_patient_obs[idx, 'awake_stepping_additional_second']
  }
}

## To reset calculation:
# df_patient_obs <- df_patient_obs[,!names(df_patient_obs) %in% c('sleep_additional_hour',
#                                                                 'sleep_additional_second',
#                                                                 'awake_additional_hour',
#                                                                 'awake_additional_second')]

## Leave the step counts as it is; they're attributed to the start time (we don't know how many steps are conducted per hour)
## Assuming that they conduct it within an hour period

################################################################################################################################################
## 02. EXPORT PREPROCESSED DATA
################################################################################################################################################

## Merge back with patient metadata from df_activpal
METADATA_FILE <- '../dataset/patient_metadata_trimmed.csv'
df_patient_metadata <- read.csv(METADATA_FILE)
df_patient_metadata$visit_info <- ifelse(
  df_patient_metadata$visit == 1, 
  'Baseline', 
  ifelse(
    df_patient_metadata$visit == 2,
    'Follow up',
    NA
    )
  )
df_patient_metadata$gender <- ifelse(
  df_patient_metadata$gender == 1, 
  'female', 
  ifelse(
    df_patient_metadata$gender == 0,
    'male',
    NA
  )
)

df_patient_metadata <- df_patient_metadata %>% 
  rename(
    sf_physical_functioning = sf_pf,
    sf_role_physical = sf_role,
    sf_bodily_pain = sf_pain,
    sf_social_functioning = sf_social,
    sf_mental_health = sf_mental,
    sf_role_emotional = sf_emot,
    sf_vitality = sf_vitality,
    sf_general_health = sf_gen_health
  ) 

df_patient_obs_with_meta <- df_patient_obs %>% 
  left_join(
    df_patient_metadata %>% select(-c('visit')),
    by = c(
      'patient_id' = 'patient_id',
      'visit_info' = 'visit_info'
    )
  )

# write.csv(
#   df_patient_obs_with_meta,
#   file = '../dataset/activpal_hourly_cleaned_duration_without_slnw_with_meta.csv', # '../dataset/activpal_hourly_cleaned_duration.csv',
#   row.names = FALSE
#   )

## Note: there are patients without follow up period observations
## should be removed later before further analysis