################################################################################################################################################
## 00. LOAD LIBRARIES
## ALSO SET WORKING DIRECTORY
################################################################################################################################################

library(dplyr)
library(reshape2)
library(ggthemes) # viz
library(ggplot2) # viz
library(lubridate)
theme_set(theme_economist())

## set working directory  
DIRECTORY <- '/PUT-THE-WORKING-DIRECTORY/'
setwd(DIRECTORY)

## define colorblind palette
cbbPalette <- c("#000000", "#E69F00", "#0072B2", "#66FF00", "#FF33FF", "#56B4E9", "#D55E00", "#CC79A7")

################################################################################################################################################
## 01. DATA IMPORT AND PREPARATION
## 01-1. REMOVE PATIENTS WHO DO NOT HAVE BOTH BASELINE AND FOLLOW-UP PERIOD
################################################################################################################################################

df_actigraph <- read.csv('../dataset/actigraph_preprocessed_choi.csv')

## only use patients with pre and post surgery data in this model
actigraph_patient_obs_check <- df_actigraph %>% 
  filter(
    directory_name %in% c('Baseline','Follow up')
    ) %>% 
  group_by(patient_id) %>% 
  summarise(
    observation_count = length(unique(directory_name))
    )

print(
  paste(
    'Number of patients - Actigraph:', 
    length(unique(df_actigraph$patient_id)),
    sep = ' '
  )
)

actigraph_patient_obs_check <- actigraph_patient_obs_check %>% filter(observation_count == 2)
print(
  paste(
    'Number of patients - Actigraph (have both baseline and follow-up data):', 
    length(actigraph_patient_obs_check$patient_id),
    sep = ' '
  )
)

df_actigraph_full_visit <- df_actigraph %>% 
  filter(
    df_actigraph$patient_id %in% actigraph_patient_obs_check$patient_id
    ) %>% 
  mutate(
    day_hour = paste(day_of_week, hour, sep='')
    )

df_actigraph_full_visit$day_hour <- factor(
  df_actigraph_full_visit$day_hour, 
  levels = c(
    '00','01','02','03','04','05','06','07','08','09',
    '010','011','012','013','014','015','016','017','018','019',
    '020','021','022','023',
    '10','11','12','13','14','15','16','17','18','19',
    '110','111','112','113','114','115','116','117','118','119',
    '120','121','122','123',
    '20','21','22','23','24','25','26','27','28','29',
    '210','211','212','213','214','215','216','217','218','219',
    '220','221','222','223',
    '30','31','32','33','34','35','36','37','38','39',
    '310','311','312','313','314','315','316','317','318','319',
    '320','321','322','323',
    '40','41','42','43','44','45','46','47','48','49',
    '410','411','412','413','414','415','416','417','418','419',
    '420','421','422','423',
    '50','51','52','53','54','55','56','57','58','59',
    '510','511','512','513','514','515','516','517','518','519',
    '520','521','522','523',
    '60','61','62','63','64','65','66','67','68','69',
    '610','611','612','613','614','615','616','617','618','619',
    '620','621','622','623'
  )
)

df_actigraph_full_visit$datetime <- as.POSIXct(df_actigraph_full_visit$datetime)
df_actigraph_full_visit$minute <- minute(df_actigraph_full_visit$datetime)

rm(df_actigraph)

################################################################################################################################################
## 01. DATA IMPORT AND PREPARATION
## 01-2. REMOVE PATIENTS WHO DO NOT HAVE SUFFICIENT VALID DAYS
## VALID DAYS: >= 10 WEAR HOURS FOR AT LEAST 4 DAYS
## THERE ARE MISSING DATA IN TUE 10AM-1PM (BASELINE) AND THU 10AM-1PM (FOLLOW-UP)
## THOSE ARE CLOSE TO THE END OF OBSERVATION PERIOD -- THESE PERIODS ARE EXCLUDED FROM VALID HOUR CALCULATION
## NOTE: based on visualisation in the latter section
################################################################################################################################################

### ACTIGRAPH - check number of waking wear hour
### Only take patients that have >= 10 wear hours in >= 4 days
VALID_DAY_WEAR_HOUR_THRESHOLD <- 10
VALID_DAY_THRESHOLD <- 4
check_valid_day <- df_actigraph_full_visit %>% 
  filter(
    (behavior_final != 'nonwear')
  ) %>% 
  group_by(patient_id, directory_name, date, day_hour) %>% 
  count() %>% 
  group_by(patient_id, directory_name, date) %>% 
  summarise(
    wear_hour = sum(n) / 60
  ) %>% 
  filter(wear_hour >= VALID_DAY_WEAR_HOUR_THRESHOLD) %>%
  group_by(patient_id, directory_name) %>%
  count() %>% 
  filter(n >= VALID_DAY_THRESHOLD) %>% 
  group_by(patient_id) %>% 
  summarise(
    is_valid = n_distinct(directory_name)
  )

PATIENT_ID_VALID_DAYS <- check_valid_day %>%
  filter(is_valid == 2) %>%
  select(patient_id)
PATIENT_ID_VALID_DAYS <- PATIENT_ID_VALID_DAYS$patient_id

print(
  paste(
    'Number of patients - Actigraph (have both baseline and follow-up data) & valid wear days:', 
    length(PATIENT_ID_VALID_DAYS),
    sep = ' '
  )
)

## Observe % of invalid data based on patients who complete both the baseline and follow-up period
df_actigraph_full_visit <- df_actigraph_full_visit %>% 
  mutate(
    is_invalid_days_or_nonwear = case_when(
      (
        (is_nonwear == 1) |
        !(patient_id %in% PATIENT_ID_VALID_DAYS)
      ) ~ 1,
      TRUE ~ 0
    )
  )

## Based on day and hour
actigraph_invalid_time_of_day <- with(
  df_actigraph_full_visit, 
  tapply(is_invalid_days_or_nonwear, list(day_of_week_string, hour), sum)
) / 
  with(
    df_actigraph_full_visit, 
    tapply(is_invalid_days_or_nonwear, list(day_of_week_string, hour), length)
  )

df_actigraph_invalid <- data.frame(
  t(actigraph_invalid_time_of_day), 
  "hour" = row.names(t(actigraph_invalid_time_of_day)), 
  "device" = "actigraph"
)

## Based on day and period
actigraph_invalid_day_period <- with(
  df_actigraph_full_visit, 
  tapply(is_invalid_days_or_nonwear, list(day_of_week_string, directory_name), sum)
) / 
  with(
    df_actigraph_full_visit, 
    tapply(is_invalid_days_or_nonwear, list(day_of_week_string, directory_name), length)
  )

df_actigraph_invalid <- data.frame(
  t(actigraph_invalid_day_period), 
  "period" = row.names(t(actigraph_invalid_day_period)), 
  "device" = "actigraph"
)

## Remove patients with insufficient number of valid days
df_actigraph_full_visit <- df_actigraph_full_visit %>% 
  filter(patient_id %in% PATIENT_ID_VALID_DAYS) 

################################################################################################################################################
## 01. DATA IMPORT AND PREPARATION
## 01-3. CREATE TIME SEQUENCE PER DAY-HOUR
################################################################################################################################################

## currently the time_seq is reset after each day_hour
df_actigraph_full_visit <- df_actigraph_full_visit %>% 
  group_by(day_hour, directory_name) %>% 
  mutate(
    time_seq = ceiling(cur_group_id() / 2)
  )

df_actigraph_full_visit$behavior_final <- factor(
  df_actigraph_full_visit$behavior_final,
  levels = c('nonwear','sedentary_behaviour','light_intensity','moderate','vigorous')
)

## Export data for further use
# write.csv(
#   df_actigraph_full_visit,
#   paste0(DIRECTORY, "../dataset/actigraph_preprocessed_choi_valid_patients53.csv"),
#   row.names = FALSE
#   )

################################################################################################################################################
## 02. VISUALISATION
## 02-1. CHECK IDENTIFIED BEHAVIOUR
################################################################################################################################################

## Sanity check of the identified behavior (per minute)
## Generate full record of patient_id per minute, based on their observation timeframe
patient_timeframe <- df_actigraph_full_visit %>% 
  group_by(patient_id, visit_info) %>% summarise(
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
      'min')
    tmp_df <- data.frame(
      patient_id = pid,
      visit_info = period,
      datetime_minute = tmp_seq
    )
    df_patient_obs <- rbind(
      df_patient_obs,
      tmp_df
    )
  }
}

df_patient_obs$datetime_minute <- as.POSIXct(df_patient_obs$datetime_minute)
df_patient_obs$hour <- hour(df_patient_obs$datetime_minute)
df_patient_obs$minute <- minute(df_patient_obs$datetime_minute)
df_patient_obs$day_of_week_string <- weekdays(df_patient_obs$datetime_minute)
## day_of_week 0 = Sunday, 6 = Saturday
df_patient_obs$day_of_week <- as.numeric(format(df_patient_obs$datetime_minute,'%w'))
df_patient_obs$day_of_week <- ifelse(
  df_patient_obs$day_of_week_string == 'Sunday', 
  6,
  ifelse(
    df_patient_obs$day_of_week_string == 'Monday',
    0,
    df_patient_obs$day_of_week - 1
  )
)

## Create time sequence per day-hour-minute
df_patient_obs <- df_patient_obs %>% 
  mutate(
    day_hour_minute = paste(day_of_week, hour, minute, sep='')
    )

df_patient_obs <- df_patient_obs %>% 
  group_by(day_of_week, hour, minute, visit_info) %>% 
  mutate(
    time_seq = ceiling(cur_group_id() / 2)
  )

## Merge the full period and the observed data
df_actigraph_minute_behavior <- df_patient_obs %>% 
  left_join(
    df_actigraph_full_visit %>% select(patient_id, visit_info, day_of_week_string, hour, minute, behavior_final, axis1),
    by = c(
      'patient_id' = 'patient_id',
      'visit_info' = 'visit_info',
      'day_of_week_string' = 'day_of_week_string',
      'hour' = 'hour',
      'minute' = 'minute'
      )
    )

## Ensure factor levels
df_actigraph_minute_behavior$behavior_final <- factor(
  df_actigraph_minute_behavior$behavior_final, 
  levels = c('nonwear','sedentary_behaviour','light_intensity','moderate','vigorous')
)

## Observe activity sequence per patient-period
ggplot(data = df_actigraph_minute_behavior) +
  geom_point(
    aes(
      x = time_seq,
      y = patient_id,
      col = behavior_final
    )
  ) + 
  scale_color_manual(
    name = 'Behavior Type',
    labels = c('Non-wear','Sedentary', 'Light', 'Moderate', 'Vigorous'),
    values = cbbPalette
  ) +
  facet_grid(
    visit_info ~ .
  ) + 
  scale_x_continuous(
    name = '',
    labels = c('Mon 12AM','Tue 12AM','Wed 12AM','Thu 12AM','Fri 12AM','Sat 12AM','Sun 12AM','Mon 12AM'),
    breaks = seq(0, 10080, 1440)
  ) + 
  scale_y_continuous(
    name = 'patient ID'
  )

## Sanity check, plot one patient's data
PID <- 54
PERIOD <- "Baseline"
pid_sub <- df_actigraph_minute_behavior[
  (df_actigraph_minute_behavior$patient_id == PID) &
    (df_actigraph_minute_behavior$directory_name == PERIOD),
  ]

ggplot(data = pid_sub,
       aes(
         x = datetime_minute,
         y = axis1,
         col = behavior_final
       )
       ) +
  geom_point()

################################################################################################################################################
## 02. VISUALISATION
## 02-2. CHECK PROPORTION OF ACTIVITY
################################################################################################################################################

### Overall activity proportion
ggplot(
  data = df_actigraph_full_visit %>% filter(
    (behavior_final != 'nonwear')
  ),
  aes(
    x = time_seq,
    fill = behavior_final
  )
) + 
  stat_bin(
    geom = 'bar',
    position = 'fill',
    binwidth = 4
    ) + 
  scale_fill_brewer(
    name = 'Behavior Type',
    labels = c('Sedentary', 'Light', 'Moderate', 'Vigorous'),
    palette = 'Spectral'
    ) +
  facet_grid(
    directory_name ~ .
    ) +
  scale_x_continuous(
    name = '',
    labels = c('Mon 12AM','Tue 12AM','Wed 12AM','Thu 12AM','Fri 12AM','Sat 12AM','Sun 12AM','Mon 12AM'),
    breaks = seq(0, 168, 24)
    ) +
  scale_y_continuous(
    name = 'proportion of activity',
    breaks = seq(0,1,.2),
    labels = scales::percent
    ) + 
  ggtitle(
    label = 'Fewer physical activities are observed during midnight in the follow-up period',
    subtitle = 'Almost none vigorous activity is observed'
    )

## Observe activity proportion after removing non-wear time and period closer to end of observation
df_actigraph_minute_behavior_clean <- df_actigraph_minute_behavior %>% 
  filter(
    (behavior_final != 'nonwear')
    )

## Check activity proportion per patient-period
ggplot(data = df_actigraph_minute_behavior_clean) +
  geom_bar(
    aes(
      x = patient_id,
      fill = behavior_final
    ),
    position = 'fill'
  ) +
  scale_fill_brewer(
    name = '',
    label = c(
      'Sedentary','Light','Moderate','Vigorous'
    ),
    palette = 'Spectral'
  )  + 
  coord_flip() +
  facet_grid(. ~ visit_info) +
  scale_x_continuous(
    name = 'patient ID'
  ) +
  scale_y_continuous(
    name = 'proportion of activity',
    breaks = seq(0,1,.2),
    labels = scales::percent
  ) +
  ggtitle(
    'Proportion of sedentary behavior looks similar',
    subtitle = 'but some patients perform more light physical activity'
  )  +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  )

## Export data for further use, e.g., compositional analysis
df_actigraph_minute_behavior_out <- df_patient_obs %>% 
  left_join(
    df_actigraph_full_visit %>%
      select(
        patient_id, visit_info, day_of_week_string, hour, minute, behavior_final, sf36_total, haq, gender, age, bmi, steps, date
      ),
    by = c(
      'patient_id' = 'patient_id',
      'visit_info' = 'visit_info',
      'day_of_week_string' = 'day_of_week_string',
      'hour' = 'hour',
      'minute' = 'minute'
    )
  )

df_actigraph_minute_behavior_out <- df_actigraph_minute_behavior_out %>% 
  filter(
    (behavior_final != 'nonwear')
  )

# write.csv(
#   df_actigraph_minute_behavior_out %>% select(
#     date, patient_id, visit_info, gender, age, bmi, sf36_total, haq, steps,
#     behavior_final
#     ),
#   '../dataset/actigraph_preprocessed_choi_valid_patients53_epoch_minute_clean.csv',
#   row.names = FALSE
#   )