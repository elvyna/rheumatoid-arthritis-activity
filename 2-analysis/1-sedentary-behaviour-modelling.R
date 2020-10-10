################################################################################################################################################
## 00. LOAD LIBRARIES
## ALSO SET WORKING DIRECTORY
################################################################################################################################################

library(dplyr)
library(tidyr)
library(ggthemes)
library(ggplot2)
library(lubridate)
library(compositions)
theme_set(theme_economist())

DIRECTORY <- '/PUT-THE-WORKING-DIRECTORY/'
setwd(DIRECTORY)

## assumption: most patients sleep between 9 PM to 5 AM (inclusive)
SLEEP_TIME <- c(21,22,23,0,1,2,3,4,5)

################################################################################################################################################
## 01. DATA IMPORT AND PREPARATION
## 01-1. ACTIVPAL
################################################################################################################################################

ACTIVPAL_FILE <- '../dataset/activpal_hourly_cleaned_duration_without_slnw_with_meta_patients33.csv'
df_activpal_hourly <- read.csv(ACTIVPAL_FILE)

################################################################################################################################################
## 01. DATA IMPORT AND PREPARATION
## 01-2. ACTIGRAPH
################################################################################################################################################

ACTIGRAPH_FILE <- '../dataset/actigraph_preprocessed_choi_valid_patients53.csv'
df_actigraph <- read.csv(ACTIGRAPH_FILE)

## create dummy variable to indicate baseline or follow up period
df_actigraph$is_followup_period <- as.numeric(df_actigraph$directory_name == 'Follow up')
## create dummy variable to indicate gender of the patient
df_actigraph$is_female <- as.numeric(df_actigraph$gender == 'female')
## add is_sleep_time as dummy attribute
## assumption: most patients sleep between 9 PM to 5 AM (inclusive)
df_actigraph$is_sleep_time <- as.numeric(df_actigraph$hour %in% SLEEP_TIME)
df_actigraph$is_awake_time <- (df_actigraph$is_sleep_time != 1)

################################################################################################################################################
## 01. DATA IMPORT AND PREPARATION
## 01-2. ACTIGRAPH
## 01-2A. STORE HOURLY RECORD IN A DATAFRAME, TO BE USED ON GLM
################################################################################################################################################

## create new time variable, so that the aggregation is based on day-hour of week
## not the sequence from the start of each patient's observation period
## filter out identified "non-wear" time
## put "behavior_final" in the group by, since it is the indicator of the sedentary or non-sedentary
df_actigraph_hourly <- df_actigraph %>% 
  filter(behavior_final != 'nonwear') %>% 
  group_by(
    date,
    day_of_week, 
    day_of_week_string, 
    hour, 
    day_hour,
    time_seq,
    patient_id, 
    is_female,
    age,
    height,
    weight,
    bmi,
    sf_physical_functioning,
    sf_role_physical,
    sf_bodily_pain,
    sf_social_functioning,
    sf_mental_health,
    sf_role_emotional,
    sf_vitality,
    sf_general_health,
    sf36_total,
    haq,
    ntx,
    oc,
    directory_name,
    is_followup_period,
    is_sleep_time,
    is_awake_time,
    behavior_final
  ) %>% 
  summarise(
    total_step = sum(steps),
    duration_minute = n()
  )

## convert "behavior_final" into separate variables
df_actigraph_hourly <- df_actigraph_hourly %>% 
  spread(
    key = behavior_final,
    value = duration_minute
  ) %>% 
  group_by(
    date,
    day_hour,
    time_seq,
    patient_id, 
    is_female,
    age,
    height,
    weight,
    bmi,
    sf_physical_functioning,
    sf_role_physical,
    sf_bodily_pain,
    sf_social_functioning,
    sf_mental_health,
    sf_role_emotional,
    sf_vitality,
    sf_general_health,
    sf36_total,
    haq,
    ntx,
    oc,
    directory_name,
    is_followup_period, 
    is_sleep_time,
    is_awake_time
  ) %>% 
  summarise(
    sedentary_behaviour = sum(sedentary_behaviour, na.rm = TRUE),
    light_intensity = sum(light_intensity, na.rm = TRUE),
    moderate = sum(moderate, na.rm = TRUE),
    vigorous = sum(vigorous, na.rm = TRUE),
    total_step = sum(total_step)
  )

################################################################################################################################################
## 01. DATA IMPORT AND PREPARATION
## 01-2. ACTIGRAPH
## 01-2B. STORE DAILY RECORD IN A DATAFRAME, 
##        TO BE USED ON COMPOSITIONAL ANALYSIS
################################################################################################################################################

df_actigraph_daily <- df_actigraph %>% 
  filter(behavior_final != 'nonwear') %>% 
  group_by(
    date,
    day_of_week, 
    day_of_week_string, 
    day_hour,
    patient_id, 
    is_female,
    age,
    height,
    weight,
    bmi,
    sf_physical_functioning,
    sf_role_physical,
    sf_bodily_pain,
    sf_social_functioning,
    sf_mental_health,
    sf_role_emotional,
    sf_vitality,
    sf_general_health,
    sf36_total,
    haq,
    ntx,
    oc,
    is_followup_period,
    directory_name,
    behavior_final
  ) %>% 
  summarise(
    total_step = sum(steps),
    duration.minute = n()
  )

df_actigraph_daily <- df_actigraph_daily %>% 
  group_by(
    date,
    day_of_week, 
    day_of_week_string, 
    patient_id, 
    is_female,
    age,
    height,
    weight,
    bmi,
    sf_physical_functioning,
    sf_role_physical,
    sf_bodily_pain,
    sf_social_functioning,
    sf_mental_health,
    sf_role_emotional,
    sf_vitality,
    sf_general_health,
    sf36_total,
    haq,
    ntx,
    oc,
    is_followup_period,
    directory_name,
    behavior_final
  )  %>% 
  summarise(
    total_step = sum(total_step),
    duration.minute = sum(duration.minute)
  )

## convert "behavior_final" into separate variables
df_actigraph_daily <- df_actigraph_daily %>% 
  spread(
    key = behavior_final,
    value = duration.minute
  ) %>% 
  group_by(
    date,
    patient_id, 
    is_female,
    age,
    height,
    weight,
    bmi,
    sf_physical_functioning,
    sf_role_physical,
    sf_bodily_pain,
    sf_social_functioning,
    sf_mental_health,
    sf_role_emotional,
    sf_vitality,
    sf_general_health,
    sf36_total,
    haq,
    ntx,
    oc,
    directory_name,
    is_followup_period
  ) %>% 
  summarise(
    sedentary_behaviour = sum(sedentary_behaviour, na.rm = TRUE),
    light_intensity = sum(light_intensity, na.rm = TRUE),
    moderate = sum(moderate, na.rm = TRUE),
    vigorous = sum(vigorous, na.rm = TRUE),
    total_step = sum(total_step)
  )

################################################################################################################################################
## 02. INITIAL PROPORTION DATA MODELLING USING LOGISTIC REGRESSION
## 02-1. ACTIVPAL
################################################################################################################################################

### ACTIVPAL - calculate hourly composition
ACTIVPAL_ACTIVITY_DURATION_COLUMNS <- c('awake_sitting_duration_second_clean','awake_standing_duration_second_clean',
                                        'awake_stepping_duration_second_clean')

activpal_hourly_composition <- acomp(
  df_activpal_hourly[, ACTIVPAL_ACTIVITY_DURATION_COLUMNS] 
)
activpal_hourly_composition <- as.data.frame(activpal_hourly_composition)
colnames(activpal_hourly_composition) <- c('sitting_proportion','standing_proportion','stepping_proportion')

df_activpal_hourly_with_prop <- dplyr::bind_cols(df_activpal_hourly, activpal_hourly_composition)
df_activpal_hourly_with_prop$sedentary_proportion <- df_activpal_hourly_with_prop$sitting_proportion
## remove record with NA
df_activpal_hourly_with_prop <- df_activpal_hourly_with_prop %>% 
  filter(!is.na(sedentary_proportion) & !is.na(haq) & !is.na(sf36_total))

df_activpal_hourly_with_prop$sedentary_minute <- df_activpal_hourly_with_prop$awake_sitting_duration_second_clean / 60
df_activpal_hourly_with_prop$non_sedentary_minute <- (df_activpal_hourly_with_prop$awake_standing_duration_second_clean + df_activpal_hourly_with_prop$awake_stepping_duration_second_clean) / 60

## modeling, to find effect of haq and sf36_total on sb
## input of GLM with binomial family (logistic regression) is binary data - can be represented as m x 2 matrix
## where the first column is the count of success (sedentary) and the second one is the count of fail (non-sedentary)
activpal_sb_hourly_matrix <- matrix(data = c(df_activpal_hourly_with_prop$sedentary_minute, df_activpal_hourly_with_prop$non_sedentary_minute), 
                                    ncol = 2,
                                    nrow = length(df_activpal_hourly_with_prop$sedentary_minute)
)
activpal_sb_hourly_matrix <- round(activpal_sb_hourly_matrix)
activpal_sb_glm <- glm(
  activpal_sb_hourly_matrix ~ is_awake_time + is_female + awake_steps_total + bmi +
    is_followup_period * sf36_total + is_followup_period*haq, 
  data = df_activpal_hourly_with_prop,
  family = binomial(link = 'logit')
)

summary(activpal_sb_glm)
## check odd ratios
exp(cbind(odd.ratio = coef(activpal_sb_glm), confint(activpal_sb_glm)))
## calculate R-squared from glm, i.e, proportion of the explained variance  
# with(summary(activpal_sb_glm), 1 - deviance/null.deviance) 
with(summary(activpal_sb_glm), (null.deviance - deviance) / null.deviance)
## use inverse logit function to transform the coefficient back to probability
## we only use it for the intercept
## for other covariates, we interpret their effect in terms of odd ratios

## show model summary as latex table
# library(xtable)
# xtable(activpal_sb_glm, label = "table:ch3-glm-activpal-summary")

df_activpal_hourly_with_prop$sb_prop_pred <- predict(activpal_sb_glm, type='response')
## use inverse logit function 
df_activpal_hourly_with_prop$sb_prop_pred <- 1 / (1 + exp(-df_activpal_hourly_with_prop$sb_prop_pred))

ggplot(data = df_activpal_hourly_with_prop) +
  geom_point(aes(x = time_seq, y = sedentary_proportion), col = 'orange', alpha = 0.5) +
  geom_point(aes(x = time_seq, y = sb_prop_pred), col = 'green', alpha = 0.5) +
  facet_grid(visit_info ~ .) +
  scale_x_continuous(
    name = '',
    labels = c('Mon\n 12AM','Tue\n 12AM','Wed\n 12AM','Thu\n 12AM','Fri\n 12AM','Sat\n 12AM','Sun\n 12AM','Mon\n 12AM'),
    breaks = seq(0, 168, 24)
  ) +
  scale_y_continuous(
    name = 'proportion of sedentary behavior',
    breaks = seq(0,1,.2),
    labels = scales::percent
  ) +
  ggtitle('Activpal')

################################################################################################################################################
## 02. INITIAL PROPORTION DATA MODELLING USING LOGISTIC REGRESSION
## 02-2. ACTIGRAPH
################################################################################################################################################

ACTIGRAPH_ACTIVITY_DURATION_COLUMNS <- c('sedentary_behaviour','light_intensity','moderate','vigorous')
ACTIGRAPH_PA_DURATION_COLUMNS <- c('light_intensity','moderate','vigorous')
df_actigraph_hourly$non_sedentary <- apply(df_actigraph_hourly[, ACTIGRAPH_PA_DURATION_COLUMNS], 1, sum)

actigraph_hourly_composition <- acomp(
  df_actigraph_hourly[, ACTIGRAPH_ACTIVITY_DURATION_COLUMNS]
)
actigraph_hourly_composition <- as.data.frame(actigraph_hourly_composition)
colnames(actigraph_hourly_composition) <- c('sedentary_proportion','lipa_proportion','moderate_proportion','vigorous_proportion')

df_actigraph_hourly_with_prop <- dplyr::bind_cols(df_actigraph_hourly, actigraph_hourly_composition)

## remove record with NA
df_actigraph_hourly_with_prop <- df_actigraph_hourly_with_prop %>% 
  filter(
    !is.na(haq) & !is.na(sf36_total)
  )

## modeling, to find effect of haq and sf36_total on sb
actigraph_sb_hourly_matrix <- matrix(
  data = c(df_actigraph_hourly_with_prop$sedentary_behaviour, df_actigraph_hourly_with_prop$non_sedentary), 
  ncol = 2,
  nrow = length(df_actigraph_hourly_with_prop$sedentary_behaviour)
)
actigraph_sb_glm <- glm(
  actigraph_sb_hourly_matrix ~ is_awake_time + is_female + total_step + bmi +
    is_followup_period * sf36_total + is_followup_period*haq, 
  data = df_actigraph_hourly_with_prop,
  family = binomial(link = 'logit')
)

summary(actigraph_sb_glm)
## check odd ratios
exp(cbind(odd.ratio = coef(actigraph_sb_glm), confint(actigraph_sb_glm)))
## calculate R-squared from glm
with(summary(actigraph_sb_glm), 1 - deviance/null.deviance) 

## show model summary as latex table
# library(xtable)
# xtable(actigraph_sb_glm, label = "table:ch3-glm-actigraph-summary")

df_actigraph_hourly_with_prop$sb_prop_pred <- predict(actigraph_sb_glm, type='response')
## use inverse logit function 
df_actigraph_hourly_with_prop$sb_prop_pred <- 1 / (1 + exp(-df_actigraph_hourly_with_prop$sb_prop_pred))

ggplot(data = df_actigraph_hourly_with_prop) +
  geom_point(aes(x = time_seq, y = sedentary_proportion), col = 'orange', alpha = 0.5) +
  geom_point(aes(x = time_seq, y = sb_prop_pred), col = 'green', alpha = 0.5) +
  facet_grid(directory_name ~ .) +
  scale_x_continuous(
    name = '',
    labels = c('Mon\n 12AM','Tue\n 12AM','Wed\n 12AM','Thu\n 12AM','Fri\n 12AM','Sat\n 12AM','Sun\n 12AM','Mon\n 12AM'),
    breaks = seq(0, 168, 24)
  ) +
  scale_y_continuous(
    name = 'proportion of sedentary behavior',
    breaks = seq(0,1,.2),
    labels = scales::percent
  ) +
  ggtitle('Actigraph')

################################################################################################################################################
## 02. INITIAL PROPORTION DATA MODELLING USING LOGISTIC REGRESSION
## 02-3. VISUALISE BOTH ACTIVPAL AND ACTIGRAPH RESULTS
################################################################################################################################################

## Combine hourly steps plot between Actigraph and Activpal
sub_activpal <- df_activpal_hourly_with_prop[, c('time_seq', 'visit_info', 'sedentary_proportion', 'sb_prop_pred')]
sub_activpal$device <- 'activpal'
colnames(sub_activpal) <- c('time_seq','visit_info','sb_proportion','predicted_sb_proportion','device')

## Ensure we have comparable patients between actigraph and activpal
sub_actigraph <- df_actigraph_hourly_with_prop[
  (df_actigraph_hourly_with_prop$patient_id %in% df_activpal_hourly$patient_id), 
  c('time_seq', 'directory_name', 'sedentary_proportion', 'sb_prop_pred')
  ]
sub_actigraph$device <- 'actigraph'
colnames(sub_actigraph) <- c('time_seq','visit_info','sb_proportion','predicted_sb_proportion','device')
sub_all <- dplyr::bind_rows(sub_activpal, sub_actigraph)
rm(sub_actigraph, sub_activpal)

sub_all_avg <- sub_all %>% 
  group_by(
    time_seq,
    device,
    visit_info
  ) %>% 
  summarise(
    avg_sb_proportion = mean(sb_proportion, na.rm = TRUE),
    avg_predicted_sb_proportion = mean(predicted_sb_proportion, na.rm = TRUE)
  )

sub_all <- sub_all %>% 
  left_join(
    sub_all_avg,
    by = c(
      'time_seq' = 'time_seq',
      'device' = 'device',
      'visit_info' = 'visit_info'
    )
  )

ggplot(data = sub_all) +
  geom_point(
    mapping = aes(
      x = time_seq,
      y = sb_proportion
    ),
    color = 'orange',
    alpha = 0.25
  ) + geom_point(
    mapping = aes(
      x = time_seq,
      y = predicted_sb_proportion
    ),
    color = 'green',
    alpha = 0.25
  ) +
  geom_line(
    mapping = aes(
      x = time_seq,
      y = avg_predicted_sb_proportion
    ),
    color = 'grey8',
    alpha = 0.5
  ) +
  scale_x_continuous(
    name = '',
    labels = c('Mon\n 12AM','Tue\n 12AM','Wed\n 12AM','Thu\n 12AM','Fri\n 12AM','Sat\n 12AM','Sun\n 12AM','Mon\n 12AM'),
    breaks = seq(0, 168, 24)
  ) +
  scale_y_continuous(
    name = 'proportion of sedentary behavior',
    breaks = seq(0,1,.2),
    labels = scales::percent
  ) + 
  facet_grid(
    visit_info ~ device
  ) +
  theme(
    axis.text = element_text(size = 8)
  )

################################################################################################################################################
## 03. EXPORT DATASET FOR FURTHER USE
## 03-1. ACTIVPAL
################################################################################################################################################

# write.csv(
#   df_activpal_hourly_with_prop,
#   "../dataset/activpal_hourly_cleaned_duration_without_slnw_with_meta_patients33_with_activity_proportion.csv",
#   row.names = FALSE
# )

################################################################################################################################################
## 03. EXPORT DATASET FOR FURTHER USE
## 03-2. ACTIGRAPH
################################################################################################################################################

# write.csv(
#   df_actigraph_hourly_with_prop,
#   '../dataset/actigraph_preprocessed_choi_valid_patients53_hourly_with_activity_proportion.csv',
#   row.names = FALSE
# )
