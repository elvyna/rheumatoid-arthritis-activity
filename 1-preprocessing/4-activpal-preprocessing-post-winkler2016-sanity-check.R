################################################################################################################################################
## 00. LOAD LIBRARIES
## ALSO SET WORKING DIRECTORY
################################################################################################################################################

library(dplyr)
library(ggthemes) # viz
library(ggplot2) # viz
library(lubridate)
theme_set(theme_economist())

## set working directory  
DIRECTORY <- '/PUT-THE-WORKING-DIRECTORY/'
setwd(DIRECTORY)

SLEEP_TIME <- c(21,22,23,0,1,2,3,4,5)

################################################################################################################################################
## 01. DATA IMPORT AND PREPARATION
## 01-1. REMOVE PATIENTS WHO DON'T HAVE BOTH BASELINE AND FOLLOW UP PERIOD
################################################################################################################################################

## Import combined data after using SAS code from Winkler et al. (2016)
## Previously: 'activpal_combined_processed_in_sas_fixed_with_meta.csv"
ACTIVPAL_FILE <- "../dataset/activpal_hourly_cleaned_duration_without_slnw_with_meta.csv"

### ACTIVPAL - LOAD DATASET -- AFTER CLEANUP 2020-04-21 activpal-preparation.R
df_activpal_hourly <- read.csv(ACTIVPAL_FILE)

df_activpal_hourly$datetime_hour <- as.POSIXct(df_activpal_hourly$datetime_hour)
df_activpal_hourly$hour <- as.factor(hour(df_activpal_hourly$datetime_hour))
df_activpal_hourly$day_of_week_string <- weekdays(df_activpal_hourly$datetime_hour)
## day_of_week 0 = Sunday, 6 = Saturday
df_activpal_hourly$day_of_week <- as.numeric(format(df_activpal_hourly$datetime_hour,'%w'))

## create time sequence
df_activpal_hourly <- df_activpal_hourly %>%
  mutate(
    day_hour = paste(day_of_week, hour, sep='')
    )

## set the level to start on Monday, so that it's similar to Actigraph
## note: day_of_week 0 = Sunday, 6 = Saturday
df_activpal_hourly$day_hour <- factor(
  df_activpal_hourly$day_hour, 
  levels = c(
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
    '620','621','622','623',
    '00','01','02','03','04','05','06','07','08','09',
    '010','011','012','013','014','015','016','017','018','019',
    '020','021','022','023'
  )
)

df_activpal_hourly <- df_activpal_hourly %>%
  group_by(day_hour, visit_info) %>% 
  mutate(
    time_seq = ceiling(cur_group_id() / 2)
  )

## Create dummy variable to indicate baseline or follow up period
df_activpal_hourly$is_followup_period <- as.numeric(df_activpal_hourly$visit_info == 'Follow up')
## Create dummy variable to indicate gender of the patient
df_activpal_hourly$is_female <- as.numeric(df_activpal_hourly$gender == 'female')
## Add is_sleep_time as dummy attribute
## Assumption: most patients sleep between 9 PM to 5 AM (inclusive)
df_activpal_hourly$is_sleep_time <- as.numeric(df_activpal_hourly$hour %in% SLEEP_TIME)
df_activpal_hourly$is_awake_time <- (df_activpal_hourly$is_sleep_time != 1)

## There are patients without follow up period observations
## Should be removed later before further analysis
activpal_patient_obs_check <- df_activpal_hourly %>% 
  group_by(patient_id) %>% 
  summarise(
    observation_count = length(unique(visit_info))
  )

print(
  paste(
    "Number of patients with valid Activpal data:",
    nrow(activpal_patient_obs_check)
  )
)

print(
  paste(
    "Number of patients with valid Activpal data AND have both baseline & follow-up data:",
    sum(activpal_patient_obs_check$observation_count == 2)
  )
)

activpal_patient_obs_check <- activpal_patient_obs_check %>% filter(observation_count == 2)

df_activpal_hourly_full_obs <- df_activpal_hourly %>% filter(patient_id %in% activpal_patient_obs_check$patient_id)
rm(df_activpal_hourly)

################################################################################################################################################
## 02. EXPORT PREPROCESSED DATA
################################################################################################################################################

## Hourly data, only patients with both baseline and follow up period
# write.csv(
#   df_activpal_hourly_full_obs,
#   file = '../dataset/activpal_hourly_cleaned_duration_without_slnw_with_meta_patients33.csv',
#   row.names = FALSE
# )

################################################################################################################################################
## 03. EXPLORE DATA
## 03-1. PROPORTION OF ACTIVITY - OVERALL
################################################################################################################################################

library(reshape2)
df_activpal_hourly_full_obs_melt <- melt(
  df_activpal_hourly_full_obs[, c(
    'time_seq','gender','visit_info',
    'awake_sitting_duration_second_clean','awake_standing_duration_second_clean','awake_stepping_duration_second_clean')],
  id = c('time_seq','visit_info','gender'),
  measure.vars = c('awake_sitting_duration_second_clean','awake_standing_duration_second_clean',
                   'awake_stepping_duration_second_clean')
)

df_activpal_hourly_full_obs_melt$variable <- factor(
  df_activpal_hourly_full_obs_melt$variable,
  levels = c('awake_sitting_duration_second_clean','awake_standing_duration_second_clean','awake_stepping_duration_second_clean')
)

ggplot(
  data = df_activpal_hourly_full_obs_melt,
  aes(
    x = time_seq,
    y = value,
    fill = variable
  )
) + geom_bar(
  position = 'fill',
  stat = 'identity'
) +
  facet_grid(visit_info ~ .) +
  scale_x_continuous(
    name = '',
    labels = c('Mon\n 12AM','Tue\n 12AM','Wed\n 12AM','Thu\n 12AM','Fri\n 12AM','Sat\n 12AM','Sun\n 12AM','Mon\n 12AM'),
    breaks = seq(0, 168, 24)
  ) +
  scale_fill_brewer(
    name = '',
    label = c('Sitting','Standing','Stepping'), # 'Sleeping'
    palette = 'Spectral'
  ) +
  scale_y_continuous(
    name = 'proportion of activity',
    breaks = seq(0,1,.2),
    labels = scales::percent
  )

################################################################################################################################################
## 03. EXPLORE DATA
## 03-2. PROPORTION OF ACTIVITY - PER PATIENT-PERIOD
################################################################################################################################################

### ACTIVPAL - observe changes in behavior proportion per patient (daily proportion)
df_activpal_daily_full_obs_melt <- melt(
  df_activpal_hourly_full_obs[, c(
    'day_of_week','gender','visit_info','patient_id',
    'awake_sitting_duration_second_clean','awake_standing_duration_second_clean','awake_stepping_duration_second_clean')],
  id = c('day_of_week','visit_info','patient_id'),
  measure.vars = c('awake_sitting_duration_second_clean','awake_standing_duration_second_clean',
                   'awake_stepping_duration_second_clean')
)

df_activpal_daily_full_obs_melt$variable <- factor(
  df_activpal_daily_full_obs_melt$variable,
  levels = c('awake_sitting_duration_second_clean','awake_standing_duration_second_clean','awake_stepping_duration_second_clean')
)

df_activpal_daily_full_obs_melt$device <- 'activpal'

ggplot(
  data = df_activpal_daily_full_obs_melt,
  aes(
    x = patient_id,
    y = value,
    fill = variable
  )
) + geom_bar(
  position = 'fill',
  stat = 'identity'
) + coord_flip() +
  facet_grid(. ~ visit_info) +
  scale_x_continuous(
    name = 'patient ID'
  ) +
  scale_fill_brewer(
    name = '',
    label = c('Sitting','Standing','Stepping'), # 'Sleeping'
    palette = 'Spectral'
  ) +
  scale_y_continuous(
    name = 'proportion of activity',
    breaks = seq(0,1,.2),
    labels = scales::percent
  )

################################################################################################################################################
## 03. EXPLORE DATA
## 03-3. PROPORTION OF ACTIVITY - PER PATIENT-PERIOD
## USE THIS CODE ALONGSIDE 2020-09-06-actigraph-preprocessing-choi2011-sanity-check.R
## TO PLOT THE PROPORTIONS AS A SINGLE GRAPH OBJECT
################################################################################################################################################

sub_actigraph_behavior_patient <- df_actigraph_minute_behavior_clean %>% select(day_of_week, visit_info, patient_id, behavior_final)
sub_actigraph_behavior_patient <- sub_actigraph_behavior_patient[, !names(sub_actigraph_behavior_patient) %in% c('hour','minute')]
sub_actigraph_behavior_patient$value <- 1
sub_actigraph_behavior_patient$device <- 'actigraph'
colnames(sub_actigraph_behavior_patient) <- c('day_of_week','visit_info','patient_id','variable','value','device')

df_activpal_daily_full_obs_melt$device <- 'activpal'
sub_all_behavior_patient_daily <- dplyr::bind_rows(sub_actigraph_behavior_patient, df_activpal_daily_full_obs_melt)

ggplot(
  data = sub_all_behavior_patient_daily,
  aes(
    x = patient_id,
    y = value,
    fill = variable
  )
) + 
  geom_bar(
    position = 'fill',
    stat = 'identity'
  ) +
  facet_grid(device ~ visit_info) +
  scale_x_continuous(
    name = 'patient ID'
  ) +
  scale_fill_brewer(
    name = '',
    label = c(
      'Sedentary','Light','Moderate','Vigorous',
      'Sitting','Standing','Stepping'
    ),
    palette = 'Dark2'
  ) +
  scale_y_continuous(
    name = 'proportion of activity',
    breaks = seq(0,1,.2),
    labels = scales::percent
  ) +
  coord_flip() +
  theme(
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 8)
  )
