################################################################################################################################################
## 00. LOAD LIBRARIES
## ALSO SET WORKING DIRECTORY
################################################################################################################################################

library(dplyr)
library(mgcv) # for gam
library(mgcv.helper) ## helper library to show confint  
library(ggthemes) # viz
library(ggplot2) # viz
library(lubridate)
theme_set(theme_economist())

## set working directory  
DIRECTORY <- '/PUT-THE-WORKING-DIRECTORY/'
setwd(DIRECTORY)

## Define colorblind palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Define sleep time
SLEEP_TIME <- c(21,22,23,0,1,2,3,4,5)

################################################################################################################################################
## 01. DATA IMPORT AND PREPARATION
## 01-1. ACTIVPAL
################################################################################################################################################

ACTIVPAL_FILE <- "../dataset/activpal_hourly_cleaned_duration_without_slnw_with_meta_patients33.csv"
df_activpal_hourly <- read.csv(ACTIVPAL_FILE)

## create sine and cosine component of the assumed activity cycle
xc1 <- cos(2 * pi * df_activpal_hourly$time_seq / 24) ## cycle per 24 hours
xs1 <- sin(2 * pi * df_activpal_hourly$time_seq / 24) ## cycle per 24 hours
xc2 <- cos(4 * pi * df_activpal_hourly$time_seq / 24) ## cycle per 12 hours
xs2 <- sin(4 * pi * df_activpal_hourly$time_seq / 24) ## cycle per 12 hours
sine_cosine <- data.frame(cbind(xc1, xs1, xc2, xs2))
df_activpal_hourly <- data.frame(df_activpal_hourly, sine_cosine)

## Illustrate sine-cosine components
ggplot(data = df_activpal_hourly) +
  geom_line(
    mapping = aes(
      x = time_seq,
      y = xs1,
      color = cbbPalette[1]
    ),
    alpha = 0.75
  ) +
  geom_line(
    mapping = aes(
      x = time_seq,
      y = xs2,
      color = cbbPalette[2]
    ),
    alpha = 0.75
  ) +
  geom_line(
    mapping = aes(
      x = time_seq,
      y = xc1,
      color = cbbPalette[3]
    ),
    alpha = 0.75
  ) +
  geom_line(
    mapping = aes(
      x = time_seq,
      y = xc2,
      color = cbbPalette[4]
    ),
    alpha = 0.75
  ) +
  scale_x_continuous(
    name = '',
    labels = c('Mon\n 12AM','Tue\n 12AM','Wed\n 12AM','Thu\n 12AM','Fri\n 12AM','Sat\n 12AM','Sun\n 12AM','Mon\n 12AM'),
    breaks = seq(0, 168, 24)
  ) +
  scale_y_continuous(
    name = 'value'
  ) +
  scale_colour_manual(
    values = cbbPalette[1:4],
    name = '',
    labels = c('sin 24 hours', 'sin 12 hours', 'cos 24 hours', 'cos 12 hours')
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

################################################################################################################################################
## 01. DATA IMPORT AND PREPARATION
## 01-2. ACTIGRAPH
################################################################################################################################################

ACTIGRAPH_FILE <- "../dataset/actigraph_preprocessed_choi_valid_patients53.csv"
df_actigraph <- read.csv(ACTIGRAPH_FILE)

## create new time variable, so that the aggregation is based on day-hour of week
## not the sequence from the start of each patient's observation period
## filter out identified "non-wear" time
df_actigraph_hourly <- df_actigraph %>% 
  filter(behavior_final != 'nonwear') %>% 
  group_by(
    day_of_week, 
    day_of_week_string, 
    hour, 
    day_hour,
    time_seq,
    patient_id, 
    gender,
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
    directory_name
  ) %>% 
  summarise(
    total_step = sum(steps),
    mean_step = mean(steps),
    total_axis1 = sum(axis1),
    total_axis2 = sum(axis2),
    total_axis3 = sum(axis3),
    total_vector_magnitude = sum(vector_magnitude),
    mean_vector_magnitude = mean(vector_magnitude)
  )

## create dummy variable to indicate baseline or follow up period
df_actigraph_hourly$is_followup_period <- as.numeric(df_actigraph_hourly$directory_name == 'Follow up')
## create dummy variable to indicate gender of the patient
df_actigraph_hourly$is_female <- as.numeric(df_actigraph_hourly$gender == 'female')
## add is_sleep_time as dummy attribute
## assumption: most patients sleep between 9 PM to 5 AM (inclusive)
df_actigraph_hourly$is_sleep_time <- as.numeric(df_actigraph_hourly$hour %in% SLEEP_TIME)
df_actigraph_hourly$is_awake_time <- (df_actigraph_hourly$is_sleep_time != 1)

################################################################################################################################################
## 01. DATA IMPORT AND PREPARATION
## 01-2b. ACTIGRAPH - DECOMPOSE TIME SERIES COMPONENTS
################################################################################################################################################

## create sine-cosine covariates to indicate 12-hourly and 24-hourly cycle
xc1 <- cos(2 * pi * df_actigraph_hourly$time_seq / 24) ## cycle per 24 hours
xs1 <- sin(2 * pi * df_actigraph_hourly$time_seq / 24) ## cycle per 24 hours
xc2 <- cos(4 * pi * df_actigraph_hourly$time_seq / 24) ## cycle per 12 hours
xs2 <- sin(4 * pi * df_actigraph_hourly$time_seq / 24) ## cycle per 12 hours
sine_cosine <- data.frame(cbind(xc1, xs1, xc2, xs2))
df_actigraph_hourly <- data.frame(df_actigraph_hourly, sine_cosine)

################################################################################################################################################
## 02. EXPLORATORY ANALYSIS
## 02-1. CHECK AVG HOURLY STEPS
## 02-1A. ACTIVPAL
################################################################################################################################################

### Plot the mean of hourly step counts
df_activpal_hourly %>%
  group_by(time_seq, visit_info) %>% 
  summarise(avg_step = mean(awake_steps_total)) %>%
  ggplot() + 
  geom_line(
    aes(
      x = time_seq, 
      y = avg_step
    )
  ) +
  scale_x_continuous(
    name = '',
    labels = c('Mon\n 12AM','Tue\n 12AM','Wed\n 12AM','Thu\n 12AM','Fri\n 12AM','Sat\n 12AM','Sun\n 12AM','Mon\n 12AM'),
    breaks = seq(0, 168, 24)
  ) +
  facet_grid(
    visit_info ~ .
  )

################################################################################################################################################
## 02. EXPLORATORY ANALYSIS
## 02-1. CHECK AVG HOURLY STEPS
## 02-1B. ACTIGRAPH
################################################################################################################################################

df_actigraph_hourly_extended %>% 
  group_by(time_seq, directory_name) %>% 
  summarise(avg_step = mean(total_step)) %>% 
  ggplot() +
  geom_line(
    aes(
      x = time_seq,
      y = avg_step
    )
  ) +
  scale_x_continuous(
    name = '',
    labels = c('Mon\n 12AM','Tue\n 12AM','Wed\n 12AM','Thu\n 12AM','Fri\n 12AM','Sat\n 12AM','Sun\n 12AM','Mon\n 12AM'),
    breaks = seq(0, 168, 24)
  ) +
  facet_grid(
    directory_name ~ .
  )

ggplot(
  df_actigraph_hourly_extended
) +
  geom_point(
    aes(
      x = time_seq,
      y = total_step
    ),
    alpha = 0.25
  ) +
  scale_x_continuous(
    name = '',
    labels = c('Mon\n 12AM','Tue\n 12AM','Wed\n 12AM','Thu\n 12AM','Fri\n 12AM','Sat\n 12AM','Sun\n 12AM','Mon\n 12AM'),
    breaks = seq(0, 168, 24)
  ) +
  facet_grid(
    directory_name ~ .
  )

################################################################################################################################################
## 03. STEP COUNT MODELLING
## 03-1. USE SINE-COSINE AS COVARIATES
## 03-1A. ACTIVPAL
################################################################################################################################################

### ACTIVPAL - INITIAL MODEL
## note: the step counts are attributed to the start time (we don't know how many steps are conducted per hour)
## we assume the patients do not walk continuously for more than an hour

activpal_model_gam <- gam(
  awake_steps_total + 1 ~ s(time_seq, bs = 'cc')  + is_followup_period + is_female + xs1 + xc1 + xs2 + xc2 +
    bmi + sf_role_emotional + age + is_followup_period*is_female, 
  family = Gamma(link='log'),
  data = df_activpal_hourly
)

plot(activpal_model_gam)
summary(activpal_model_gam)
# plot(activpal_model_gam$residuals)
confint(activpal_model_gam)
cbind(exp(confint(activpal_model_gam)[, c("Estimate","2.5%","97.5%")]))

## show model summary as latex table
# library(itsadug)
# gamtabs(activpal_model_gam, label = "table:ch3-activpal-step-count-model-initial")

## store prediction
df_activpal_hourly$predicted <- exp(
  predict(
    activpal_model_gam,
    newdata = df_activpal_hourly
  )
)

## plot the prediction
ggplot(data = df_activpal_hourly) +
  geom_point(
    mapping = aes(
      x = time_seq,
      y = awake_steps_total
    ),
    color = 'orange',
    alpha = 0.25
  ) +
  geom_point(
    mapping = aes(
      x = time_seq,
      y = predicted
    ),
    color = 'green',
    alpha = 0.25
  ) +
  scale_x_continuous(
    name = '',
    labels = c('Mon\n 12AM','Tue\n 12AM','Wed\n 12AM','Thu\n 12AM','Fri\n 12AM','Sat\n 12AM','Sun\n 12AM','Mon\n 12AM'),
    breaks = seq(0, 168, 24)
  ) +
  scale_y_continuous(
    name = 'step count'
  ) +
  facet_grid(visit_info ~ .) +
  theme(
    axis.text = element_text(size = 8)
  )

################################################################################################################################################
## 03. STEP COUNT MODELLING
## 03-1. USE SINE-COSINE AS COVARIATES
## 03-1B. ACTIGRAPH
################################################################################################################################################

actigraph_model_gam <- gam(
  total_step + 1 ~ s(time_seq, bs = 'cc') + is_followup_period * is_female + bmi + xs1 + xc1 + xs2 + xc2 +
    is_followup_period*sf_role_emotional,
  family = Gamma(link='log'),
  data = df_actigraph_hourly_extended
)

plot(actigraph_model_gam)
summary(actigraph_model_gam)
confint(actigraph_model_gam)
cbind(exp(confint(actigraph_model_gam)[, c("Estimate","2.5%","97.5%")]))

## show model summary as latex table
# library(itsadug)
# gamtabs(actigraph_model_gam, label = "table:ch3-actigraph-step-count-model-initial")

## store the prediction result
df_actigraph_hourly_extended$predicted <- exp(predict(actigraph_model_gam, newdata = df_actigraph_hourly_extended))

ggplot(data = df_actigraph_hourly_extended) +
  # ggplot(data = df_actigraph_hourly_extended %>% filter(patient_id %in% activpal_patient_obs_check$patient_id)) + 
  geom_point(
    mapping = aes(
      x = time_seq,
      y = total_step
    ), 
    color = 'orange',
    alpha = 0.25
  ) + 
  geom_point(
    mapping = aes(
      x = time_seq,
      y = predicted
    ),
    color = 'green',
    alpha = 0.25
  ) +
  ggtitle(
    ''
  ) +
  scale_x_continuous(
    name = '',
    labels = c('Mon\n 12AM','Tue\n 12AM','Wed\n 12AM','Thu\n 12AM','Fri\n 12AM','Sat\n 12AM','Sun\n 12AM','Mon\n 12AM'),
    breaks = seq(0, 168, 24)
  ) +
  scale_y_continuous(
    name = 'total steps',
    breaks = seq(0, 7000, 1000),
    limits = c(0, 7000)
  ) + 
  facet_grid(
    directory_name ~ .
  ) +
  theme(
    axis.text = element_text(size = 8)
  )

################################################################################################################################################
## 03. STEP COUNT MODELLING
## 03-1. USE SINE-COSINE AS COVARIATES
## 03-1C. VISUALISE PREDICTION IN BOTH DEVICES
## PLOT OBSERVED DATA, OVERLAYED BY PREDICTED VALUE AND AVG PREDICTED VALUE
################################################################################################################################################

## Combine hourly steps plot between Actigraph and Activpal; to be visualized in the same plot object
sub_activpal <- df_activpal_hourly %>% select(time_seq, hour, awake_steps_total, predicted, visit_info, gender)
sub_activpal$device <- 'activpal'
colnames(sub_activpal) <- c('time_seq','hour','step_total','predicted_step_total','visit_info','gender','device')

sub_actigraph <- df_actigraph_hourly_extended %>%
  select(time_seq, hour, total_step, predicted, directory_name, gender)
sub_actigraph$device <- 'actigraph'
colnames(sub_actigraph) <- c('time_seq','hour','step_total','predicted_step_total','visit_info','gender','device')

sub_all <- rbind(sub_activpal, sub_actigraph)
rm(sub_actigraph, sub_activpal)

sub_all_avg <- sub_all %>% 
  group_by(
    time_seq,
    device,
    visit_info
  ) %>% 
  summarise(
    avg_step_total = mean(step_total, na.rm = TRUE),
    avg_predicted_step_total = mean(predicted_step_total, na.rm = TRUE)
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
      y = step_total
    ), 
    color = 'orange',
    alpha = 0.25
  ) + 
  geom_point(
    mapping = aes(
      x = time_seq,
      y = predicted_step_total
    ),
    color = 'green',
    alpha = 0.25
  ) +
  geom_line(
    mapping = aes(
      x = time_seq,
      y = avg_predicted_step_total
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
    name = 'total steps',
    breaks = seq(0, 7000, 1000),
    limits = c(0, 7000)
  ) + 
  facet_grid(
    visit_info ~ device
  ) +
  theme(
    axis.text = element_text(size = 8)
  )

################################################################################################################################################
## 03. STEP COUNT MODELLING
## 03-2. USE IS_AWAKE_TIME DUMMY VARIABLE TO REPLACE SINCE-COSINE COVARIATES
## 03-2A. ACTIVPAL
################################################################################################################################################

activpal_model_gam <- gam(
  awake_steps_total + 1 ~ s(time_seq, bs = 'cc') + is_followup_period + is_female + age + bmi +
    sf_role_emotional + is_awake_time*is_followup_period, 
  family = Gamma(link='log'),
  data = df_activpal_hourly
)

plot(activpal_model_gam)
summary(activpal_model_gam)
# plot(activpal_model_gam$residuals)
confint(activpal_model_gam)
cbind(exp(confint(activpal_model_gam)[, c("Estimate","2.5%","97.5%")]))

## show model summary as latex table
# library(itsadug)
# gamtabs(activpal_model_gam, label = "table:ch3-activpal-step-count-model-current")

## store prediction
df_activpal_hourly$predicted <- exp(
  predict(
    activpal_model_gam,
    newdata = df_activpal_hourly
  )
)

## plot the prediction
ggplot(data = df_activpal_hourly) +
  geom_point(
    mapping = aes(
      x = time_seq,
      y = awake_steps_total
    ),
    color = 'orange',
    alpha = 0.25
  ) +
  geom_point(
    mapping = aes(
      x = time_seq,
      y = predicted
    ),
    color = 'green',
    alpha = 0.25
  ) +
  scale_x_continuous(
    name = '',
    labels = c('Mon\n 12AM','Tue\n 12AM','Wed\n 12AM','Thu\n 12AM','Fri\n 12AM','Sat\n 12AM','Sun\n 12AM','Mon\n 12AM'),
    breaks = seq(0, 168, 24)
  ) +
  scale_y_continuous(
    name = 'step count'
  ) +
  facet_grid(visit_info ~ .) +
  theme(
    axis.text = element_text(size = 8)
  )

################################################################################################################################################
## 03. STEP COUNT MODELLING
## 03-2. USE IS_AWAKE_TIME DUMMY VARIABLE TO REPLACE SINCE-COSINE COVARIATES
## 03-2B. ACTIGRAPH
################################################################################################################################################

actigraph_model_gam <- gam(
  total_step + 1 ~ s(time_seq, bs = 'cc') + is_followup_period * is_female + bmi + # upper_bound_step_seasonality +
    is_awake_time  + is_awake_time*is_followup_period,
  family = Gamma(link='log'),
  data = df_actigraph_hourly_extended
)

plot(actigraph_model_gam)
summary(actigraph_model_gam)
confint(actigraph_model_gam)
cbind(exp(confint(actigraph_model_gam)[, c("Estimate","2.5%","97.5%")]))

## show model summary as latex table
# library(itsadug)
# gamtabs(actigraph_model_gam, label = "table:ch3-actigraph-step-count-model-current")

## store the prediction result
df_actigraph_hourly_extended$predicted <- exp(predict(actigraph_model_gam, newdata = df_actigraph_hourly_extended))

ggplot(data = df_actigraph_hourly_extended) +
  # ggplot(data = df_actigraph_hourly_extended %>% filter(patient_id %in% activpal_patient_obs_check$patient_id)) + 
  geom_point(
    mapping = aes(
      x = time_seq,
      y = total_step
    ), 
    color = 'orange',
    alpha = 0.25
  ) + 
  geom_point(
    mapping = aes(
      x = time_seq,
      y = predicted
    ),
    color = 'green',
    alpha = 0.25
  ) +
  ggtitle(
    ''
  ) +
  scale_x_continuous(
    name = '',
    labels = c('Mon\n 12AM','Tue\n 12AM','Wed\n 12AM','Thu\n 12AM','Fri\n 12AM','Sat\n 12AM','Sun\n 12AM','Mon\n 12AM'),
    breaks = seq(0, 168, 24)
  ) +
  scale_y_continuous(
    name = 'total steps',
    breaks = seq(0, 7000, 1000),
    limits = c(0, 7000)
  ) + 
  facet_grid(
    directory_name ~ .
  ) +
  theme(
    axis.text = element_text(size = 8)
  )

################################################################################################################################################
## 03. STEP COUNT MODELLING
## 03-2. USE IS_AWAKE_TIME DUMMY VARIABLE TO REPLACE SINCE-COSINE COVARIATES
## 03-2C. VISUALISE PREDICTION IN BOTH DEVICES
## PLOT OBSERVED DATA, OVERLAYED BY PREDICTED VALUE AND AVG PREDICTED VALUE
################################################################################################################################################

## Combine hourly steps plot between Actigraph and Activpal; to be visualized in the same plot object
sub_activpal <- df_activpal_hourly %>% select(time_seq, hour, awake_steps_total, predicted, visit_info, gender)
sub_activpal$device <- 'activpal'
colnames(sub_activpal) <- c('time_seq','hour','step_total','predicted_step_total','visit_info','gender','device')

sub_actigraph <- df_actigraph_hourly_extended %>%
  # filter(patient_id %in% df_activpal_hourly$patient_id) %>% ## use this only for EDA, to get an apple to apple comparison with activpal
  select(time_seq, hour, total_step, predicted, directory_name, gender)
sub_actigraph$device <- 'actigraph'
colnames(sub_actigraph) <- c('time_seq','hour','step_total','predicted_step_total','visit_info','gender','device')

sub_all <- rbind(sub_activpal, sub_actigraph)
rm(sub_actigraph, sub_activpal)

sub_all_avg <- sub_all %>% 
  group_by(
    time_seq,
    device,
    visit_info
  ) %>% 
  summarise(
    avg_step_total = mean(step_total, na.rm = TRUE),
    avg_predicted_step_total = mean(predicted_step_total, na.rm = TRUE)
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
      y = step_total
    ), 
    color = 'orange',
    alpha = 0.25
  ) + 
  geom_point(
    mapping = aes(
      x = time_seq,
      y = predicted_step_total
    ),
    color = 'green',
    alpha = 0.25
  ) +
  geom_line(
    mapping = aes(
      x = time_seq,
      y = avg_predicted_step_total
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
    name = 'total steps',
    breaks = seq(0, 7000, 1000),
    limits = c(0, 7000)
  ) + 
  facet_grid(
    visit_info ~ device
  ) +
  theme(
    axis.text = element_text(size = 8)
  )

################################################################################################################################################
## 04. ADDITIONAL EXPLORATION
## 04-1. PLOT OBSERVED DATA, OVERLAYED BY HOURLY AVG
################################################################################################################################################

### Combine hourly steps plot between Actigraph and Activpal; to be visualized in the same plot object
sub_activpal <- df_activpal_hourly %>% select(time_seq, hour, awake_steps_total, predicted, visit_info, gender)
sub_activpal$device <- 'activpal'
colnames(sub_activpal) <- c('time_seq','hour','step_total','predicted_step_total','visit_info','gender','device')

## Only select patients with both actigraph and activpal data, to compare the records
sub_actigraph <- df_actigraph_hourly_extended %>%
  filter(patient_id %in% df_activpal_hourly$patient_id) %>% 
  select(time_seq, hour, total_step, predicted, directory_name, gender)
sub_actigraph$device <- 'actigraph'
colnames(sub_actigraph) <- c('time_seq','hour','step_total','predicted_step_total','visit_info','gender','device')

sub_all <- rbind(sub_activpal, sub_actigraph)
rm(sub_actigraph, sub_activpal)

sub_all_avg <- sub_all %>% 
  group_by(
    time_seq,
    device,
    visit_info
  ) %>% 
  summarise(
    avg_step_total = mean(step_total, na.rm = TRUE),
    avg_predicted_step_total = mean(predicted_step_total, na.rm = TRUE)
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
      y = step_total
    ), 
    color = 'orange',
    alpha = 0.25
  ) + 
  geom_line(
    mapping = aes(
      x = time_seq,
      y = avg_step_total
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
    name = 'total steps',
    breaks = seq(0, 7000, 1000),
    limits = c(0, 7000)
  ) + 
  facet_grid(
    visit_info ~ device
  ) +
  theme(
    axis.text = element_text(size = 8)
  )

################################################################################################################################################
## 04. STEP COUNT MODELLING
## 04-2. PLOT BOXPLOT OF OBSERVED DATA, FOR EDA SECTION
################################################################################################################################################

## Boxplot, to estimate sleep time
ggplot(data = sub_all) + 
  geom_boxplot(
    mapping = aes(
      x = hour,
      y = step_total,
      group = hour,
      fill = 'orange',
      alpha = 0.25
    )
  ) +
  scale_x_continuous(
    name = 'hour of day',
    breaks = seq(0, 23, 1)
  ) +
  scale_y_continuous(
    name = 'total steps',
    breaks = seq(0, 7000, 1000),
    limits = c(0, 7000)
  ) + 
  facet_grid(
    device ~ . 
  ) +
  theme(
    legend.position = 'none',
    axis.ticks = element_blank()
  )