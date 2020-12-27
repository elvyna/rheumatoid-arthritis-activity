################################################################################################################################################
## 00. LOAD LIBRARIES
## LOAD deltacomp-adjustment.R, which contains adjusted functions to accommodate linear mixed effects model
## ALSO SET WORKING DIRECTORY
################################################################################################################################################

library(dplyr)
library(tidyr)
library(reshape2)
library(ggthemes)
library(ggplot2)
library(lubridate)
library(deltacomp)
library(compositions)
library(lme4)
theme_set(theme_minimal())

DIRECTORY <- '/media/elvyna/DATA/uoa/compsci791-dissertation/ap-rheumatoid-arthritis/script/'
setwd(DIRECTORY)
source('deltacomp-function-adjustment.R')

## assumption: most patients sleep between 9 PM to 5 AM (inclusive)
SLEEP_TIME <- c(21,22,23,0,1,2,3,4,5)

################################################################################################################################################
## 01. DATA IMPORT AND PREPARATION
## 01-1. ACTIVPAL
################################################################################################################################################

ACTIVPAL_FILE <- "../dataset/activpal_hourly_cleaned_duration_without_slnw_with_meta_patients33_with_activity_proportion.csv"
df_activpal_hourly_with_prop <- read.csv(ACTIVPAL_FILE)

df_activpal_hourly_with_prop$date <- as.Date(df_activpal_hourly_with_prop$datetime_hour, '%Y-%m-%d')
df_activpal_hourly_with_prop$is_followup_period <- as.numeric(df_activpal_hourly_with_prop$visit_info == 'Follow up')

################################################################################################################################################
## 01. DATA IMPORT AND PREPARATION
## 01-2. ACTIGRAPH
################################################################################################################################################

ACTIGRAPH_HOURLY_FILE <- '../dataset/actigraph_preprocessed_choi_valid_patients53_hourly_with_activity_proportion.csv'
df_actigraph_hourly_with_prop <- read.csv(ACTIGRAPH_HOURLY_FILE)
df_actigraph_hourly_with_prop$is_followup_period <- as.numeric(df_actigraph_hourly_with_prop$directory_name == 'Follow up')

################################################################################################################################################
## 02. LINEAR MIXED EFFECTS MODEL ON THE COMPOSITIONAL DATA
## 02-1. ACTIVPAL
## DATA PREPARATION
################################################################################################################################################

df_activpal_daily_comp_awake <- df_activpal_hourly_with_prop %>% 
  filter(is_awake_time == 1) %>% 
  group_by(patient_id,
           visit_info,
           is_followup_period,
           date,
           age,
           bmi,
           sf36_total,
           sf_physical_functioning,
           sf_role_physical,
           sf_bodily_pain,
           sf_social_functioning,
           sf_mental_health,
           sf_role_emotional,
           sf_vitality,
           sf_general_health,
           haq) %>%
  summarise(
    sitting_minute = sum(awake_sitting_duration_second_clean) / 60,
    standing_minute = sum(awake_standing_duration_second_clean) / 60,
    stepping_minute = sum(awake_stepping_duration_second_clean) / 60
  )
activpal_comp_temp <- as.data.frame(
  acomp(df_activpal_daily_comp_awake[, c("sitting_minute", "standing_minute","stepping_minute")])
)
colnames(activpal_comp_temp) <- c('sitting_proportion','standing_proportion','stepping_proportion')
df_activpal_daily_comp_awake <- dplyr::bind_cols(df_activpal_daily_comp_awake, activpal_comp_temp)

## compute daily average per observation period
df_activpal_daily_comp_awake <- df_activpal_daily_comp_awake %>% 
  group_by(
    patient_id, 
    visit_info, 
    is_followup_period,
    age, 
    bmi, 
    sf36_total, 
    sf_physical_functioning,
    sf_role_physical,
    sf_bodily_pain,
    sf_social_functioning,
    sf_mental_health,
    sf_role_emotional,
    sf_vitality,
    sf_general_health,
    haq
  ) %>% 
  summarise(
    sitting_proportion = mean(sitting_proportion),
    standing_proportion = mean(standing_proportion),
    stepping_proportion = mean(stepping_proportion)
  )
df_activpal_daily_comp_awake <- as.data.frame(df_activpal_daily_comp_awake)


################################################################################################################################################
## 02. LINEAR MIXED EFFECTS MODEL ON THE COMPOSITIONAL DATA
## 02-1. ACTIVPAL
## 02-1A. RESPONSE VARIABLE: SF36_TOTAL
## random effects on the subjects (patient_id)
################################################################################################################################################

################################################################################################################################################
## 02-1A-1. ACTIVPAL - ILR1: SITTING
################################################################################################################################################

activpal_sf36_out1 <- predict_delta_comps(
  dataf = df_activpal_daily_comp_awake[
    ,
    c('sf36_total', 'sitting_proportion', 'standing_proportion', 'stepping_proportion', 'is_followup_period', 'bmi','patient_id')
    ],
  y = 'sf36_total',
  comps = c('sitting_proportion','standing_proportion','stepping_proportion'), 
  ## note: we can switch the order of covariate to get information of relative amount of time spent in 1st behavior
  covars = c('bmi','is_followup_period','patient_id'),
  random_effect = '(1|patient_id)',
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

plot_delta_comp(
  activpal_sf36_out1,
  comp_total = 24 * 60,
  units_lab = "min"
)

################################################################################################################################################
## 02-1A-2. ACTIVPAL - ILR1: STANDING
################################################################################################################################################

activpal_sf36_out2 <- predict_delta_comps(
  dataf = df_activpal_daily_comp_awake[
    ,
    c('sf36_total', 'sitting_proportion', 'standing_proportion', 'stepping_proportion', 'is_followup_period', 'bmi','patient_id')
    ],
  y = 'sf36_total',
  comps = c('standing_proportion','stepping_proportion','sitting_proportion'), 
  ## note: we can switch the order of covariate to get information of relative amount of time spent in 1st behavior
  covars = c('bmi','is_followup_period','patient_id'),
  random_effect = '(1|patient_id)',
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

################################################################################################################################################
## 02-1A-3. ACTIVPAL - ILR1: STEPPING
################################################################################################################################################

activpal_sf36_out3 <- predict_delta_comps(
  dataf = df_activpal_daily_comp_awake[
    ,
    c('sf36_total', 'sitting_proportion', 'standing_proportion', 'stepping_proportion', 'is_followup_period', 'bmi','patient_id')
    ],
  y = 'sf36_total',
  comps = c('stepping_proportion','sitting_proportion','standing_proportion'), 
  ## note: we can switch the order of covariate to get information of relative amount of time spent in 1st behavior
  covars = c('bmi','is_followup_period','patient_id'),
  random_effect = '(1|patient_id)',
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

################################################################################################################################################
## 02. LINEAR MIXED EFFECTS MODEL ON THE COMPOSITIONAL DATA
## 02-1. ACTIVPAL
## 02-1B. RESPONSE VARIABLE: HAQ
## random effects on the subjects (patient_id)
################################################################################################################################################

################################################################################################################################################
## 02-1B-1. ACTIVPAL - ILR1: SITTING
################################################################################################################################################

activpal_haq_out1 <- predict_delta_comps(
  dataf = df_activpal_daily_comp_awake[
    ,
    c('haq', 'sitting_proportion', 'standing_proportion', 'stepping_proportion', 'is_followup_period', 'bmi','patient_id')
    ],
  y = 'haq',
  comps = c('sitting_proportion','standing_proportion','stepping_proportion'), 
  ## note: we can switch the order of covariate to get information of relative amount of time spent in 1st behavior
  covars = c('bmi','is_followup_period','patient_id'),
  random_effect = '(1|patient_id)',
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

plot_delta_comp(
  activpal_haq_out1,
  comp_total = 24 * 60,
  units_lab = "min"
)

################################################################################################################################################
## 02-1B-2. ACTIVPAL - ILR1: STANDING
################################################################################################################################################

activpal_haq_out2 <- predict_delta_comps(
  dataf = df_activpal_daily_comp_awake[
    ,
    c('haq', 'sitting_proportion', 'standing_proportion', 'stepping_proportion', 'is_followup_period', 'bmi','patient_id')
    ],
  y = 'haq',
  comps = c('standing_proportion','stepping_proportion','sitting_proportion'), 
  ## note: we can switch the order of covariate to get information of relative amount of time spent in 1st behavior
  covars = c('bmi','is_followup_period','patient_id'),
  random_effect = '(1|patient_id)',
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

################################################################################################################################################
## 02-1B-3. ACTIVPAL - ILR1: STEPPING
################################################################################################################################################

activpal_haq_out3 <- predict_delta_comps(
  dataf = df_activpal_daily_comp_awake[
    ,
    c('haq', 'sitting_proportion', 'standing_proportion', 'stepping_proportion', 'is_followup_period', 'bmi','patient_id')
    ],
  y = 'haq',
  comps = c('stepping_proportion','sitting_proportion','standing_proportion'), 
  ## note: we can switch the order of covariate to get information of relative amount of time spent in 1st behavior
  covars = c('bmi','is_followup_period','patient_id'),
  random_effect = '(1|patient_id)',
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

################################################################################################################################################
## 02. LINEAR MIXED EFFECTS MODEL ON THE COMPOSITIONAL DATA
## 02-2. ACTIGRAPH
## DATA PREPARATION
################################################################################################################################################

df_actigraph_daily_comp_awake <- df_actigraph_hourly_with_prop %>% 
  filter(is_awake_time == 1) %>%
  group_by(patient_id, 
           directory_name, 
           is_followup_period,
           date, 
           age, 
           bmi, 
           sf36_total, 
           sf_physical_functioning,
           sf_role_physical,
           sf_bodily_pain,
           sf_social_functioning,
           sf_mental_health,
           sf_role_emotional,
           sf_vitality,
           sf_general_health,
           haq
  ) %>% 
  summarise(
    sedentary_behaviour = sum(sedentary_behaviour),
    light_intensity = sum(light_intensity),
    mvpa = sum(moderate + vigorous)
  )
actigraph_comp_temp <- as.data.frame(
  acomp(df_actigraph_daily_comp_awake[, c("sedentary_behaviour", "light_intensity","mvpa")])
)

colnames(actigraph_comp_temp) <- c('sedentary_proportion','lipa_proportion','mvpa_proportion')
df_actigraph_daily_comp_awake <- dplyr::bind_cols(df_actigraph_daily_comp_awake, actigraph_comp_temp)

## compute daily average per observation period
df_actigraph_daily_comp_awake <- df_actigraph_daily_comp_awake %>% 
  group_by(
    patient_id, 
    directory_name, 
    is_followup_period,
    age, 
    bmi, 
    sf36_total, 
    sf_physical_functioning,
    sf_role_physical,
    sf_bodily_pain,
    sf_social_functioning,
    sf_mental_health,
    sf_role_emotional,
    sf_vitality,
    sf_general_health,
    haq
  ) %>% 
  summarise(
    sedentary_proportion = mean(sedentary_proportion),
    lipa_proportion = mean(lipa_proportion),
    mvpa_proportion = mean(mvpa_proportion),
    sedentary_behaviour_raw = mean(sedentary_behaviour),
    light_intensity_raw = mean(light_intensity),
    mvpa_raw = mean(mvpa)
  )
df_actigraph_daily_comp_awake <- as.data.frame(df_actigraph_daily_comp_awake)

################################################################################################################################################
## 02. LINEAR MIXED EFFECTS MODEL ON THE COMPOSITIONAL DATA
## 02-2. ACTIGRAPH
## 02-2A. RESPONSE VARIABLE: SF36_TOTAL
## random effects on the subjects (patient_id)
################################################################################################################################################

################################################################################################################################################
## 02-2A-1. ACTIGRAPH - ILR1: SEDENTARY BEHAVIOUR
################################################################################################################################################

actigraph_sf36_out1 <- predict_delta_comps(
  dataf = df_actigraph_daily_comp_awake[
    ,
    c('sf36_total', 'sedentary_proportion','lipa_proportion','mvpa_proportion','is_followup_period','bmi','patient_id')
    ],
  y = 'sf36_total',
  comps = c('sedentary_proportion','lipa_proportion','mvpa_proportion'), 
  ## note: we can switch the order of covariate to get information of relative amount of time spent in 1st behavior
  covars = c('bmi','is_followup_period','patient_id'),
  random_effect = '(1|patient_id)',
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

plot_delta_comp(
  actigraph_sf36_out1,
  comp_total = 24 * 60,
  units_lab = "min"
)

################################################################################################################################################
## 02-2A-2. ACTIGRAPH - ILR1: LIPA
################################################################################################################################################

actigraph_sf36_out2 <- predict_delta_comps(
  dataf = df_actigraph_daily_comp_awake[
    ,
    c('sf36_total', 'sedentary_proportion','lipa_proportion','mvpa_proportion','is_followup_period','bmi','patient_id')
    ],
  y = 'sf36_total',
  comps = c('lipa_proportion','mvpa_proportion','sedentary_proportion'), 
  ## note: we can switch the order of covariate to get information of relative amount of time spent in 1st behavior
  covars = c('bmi','is_followup_period','patient_id'),
  random_effect = '(1|patient_id)',
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

################################################################################################################################################
## 02-2A-3. ACTIGRAPH - ILR1: MVPA
################################################################################################################################################

actigraph_sf36_out3 <- predict_delta_comps(
  dataf = df_actigraph_daily_comp_awake[
    ,
    c('sf36_total', 'sedentary_proportion','lipa_proportion','mvpa_proportion','is_followup_period','bmi','patient_id')
    ],
  y = 'sf36_total',
  comps = c('mvpa_proportion','sedentary_proportion','lipa_proportion'), 
  ## note: we can switch the order of covariate to get information of relative amount of time spent in 1st behavior
  covars = c('bmi','is_followup_period','patient_id'),
  random_effect = '(1|patient_id)',
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

################################################################################################################################################
## 02. LINEAR MIXED EFFECTS MODEL ON THE COMPOSITIONAL DATA
## 02-2. ACTIGRAPH
## 02-2B. RESPONSE VARIABLE: HAQ
## random effects on the subjects (patient_id)
################################################################################################################################################

################################################################################################################################################
## 02-2B-1. ACTIGRAPH - ILR1: SEDENTARY BEHAVIOUR
################################################################################################################################################

actigraph_haq_out1 <- predict_delta_comps(
  dataf = df_actigraph_daily_comp_awake[
    ,
    c('haq', 'sedentary_proportion','lipa_proportion','mvpa_proportion','is_followup_period','bmi','patient_id')
    ],
  y = 'haq',
  comps = c('sedentary_proportion','lipa_proportion','mvpa_proportion'), 
  ## note: we can switch the order of covariate to get information of relative amount of time spent in 1st behavior
  covars = c('bmi','is_followup_period','patient_id'),
  random_effect = '(1|patient_id)',
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

plot_delta_comp(
  actigraph_haq_out1,
  comp_total = 24 * 60,
  units_lab = "min"
)

################################################################################################################################################
## 02-2B-2. ACTIGRAPH - ILR1: LIPA
################################################################################################################################################

actigraph_haq_out2 <- predict_delta_comps(
  dataf = df_actigraph_daily_comp_awake[
    ,
    c('haq', 'sedentary_proportion','lipa_proportion','mvpa_proportion','is_followup_period','bmi','patient_id')
    ],
  y = 'haq',
  comps = c('lipa_proportion','mvpa_proportion','sedentary_proportion'), 
  ## note: we can switch the order of covariate to get information of relative amount of time spent in 1st behavior
  covars = c('bmi','is_followup_period','patient_id'),
  random_effect = '(1|patient_id)',
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

################################################################################################################################################
## 02-2B-3. ACTIGRAPH - ILR1: MVPA
################################################################################################################################################

actigraph_haq_out3 <- predict_delta_comps(
  dataf = df_actigraph_daily_comp_awake[
    ,
    c('haq', 'sedentary_proportion','lipa_proportion','mvpa_proportion','is_followup_period','bmi','patient_id')
    ],
  y = 'haq',
  comps = c('mvpa_proportion','sedentary_proportion','lipa_proportion'), 
  ## note: we can switch the order of covariate to get information of relative amount of time spent in 1st behavior
  covars = c('bmi','is_followup_period','patient_id'),
  random_effect = '(1|patient_id)',
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

################################################################################################################################################
## 03. OBSERVE MORE COMPLETE MODEL SUMMARY
## 03-1. ACTIVPAL
################################################################################################################################################

library(lmerTest) ## to show p-values in summary()
## to plot random effects
library(sjPlot)
library(sjmisc)
library(glmmTMB)

## sf36
activpal_sf36_out1_mdl <- predict_delta_comps(
  dataf = df_activpal_daily_comp_awake[
    ,
    c('sf36_total', 'sitting_proportion', 'standing_proportion', 'stepping_proportion', 'is_followup_period', 'bmi','patient_id')
    ],
  y = 'sf36_total',
  comps = c('sitting_proportion','standing_proportion','stepping_proportion'),
  # comps = c('standing_proportion','stepping_proportion','sitting_proportion'),
  # comps = c('stepping_proportion','sitting_proportion','standing_proportion'),
  ## note: we can switch the order of covariate to get information of relative amount of time spent in 1st behavior
  covars = c('bmi','is_followup_period','patient_id'),
  random_effect = '(1|patient_id)',
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05,
  return_model = TRUE
)
plot_model(activpal_sf36_out1_mdl, type = 're', title = 'random effects - activpal - sf36 total')

## haq
activpal_haq_out1_mdl <- predict_delta_comps(
  dataf = df_activpal_daily_comp_awake[
    ,
    c('haq', 'sitting_proportion', 'standing_proportion', 'stepping_proportion', 'is_followup_period', 'bmi','patient_id')
    ],
  y = 'haq',
  # comps = c('sitting_proportion','standing_proportion','stepping_proportion'),
  # comps = c('standing_proportion','stepping_proportion','sitting_proportion'),
  comps = c('stepping_proportion','sitting_proportion','standing_proportion'),
  ## note: we can switch the order of covariate to get information of relative amount of time spent in 1st behavior
  covars = c('bmi','is_followup_period','patient_id'),
  random_effect = '(1|patient_id)',
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05,
  return_model = TRUE
)
plot_model(activpal_haq_out1_mdl, type = 're', title = 'random effects - activpal - haq')

################################################################################################################################################
## 03. OBSERVE MORE COMPLETE MODEL SUMMARY
## 03-2. ACTIGRAPH
################################################################################################################################################

## sf36
actigraph_sf36_out1_mdl <- predict_delta_comps(
  dataf = df_actigraph_daily_comp_awake[
    ,
    c('sf36_total', 'sedentary_proportion','lipa_proportion','mvpa_proportion','is_followup_period','bmi','patient_id')
    ],
  y = 'sf36_total',
  comps = c('sedentary_proportion','lipa_proportion','mvpa_proportion'),
  # comps = c('lipa_proportion','sedentary_proportion','mvpa_proportion'), 
  # comps = c('mvpa_proportion','lipa_proportion','sedentary_proportion'), 
  ## note: we can switch the order of covariate to get information of relative amount of time spent in 1st behavior
  covars = c('bmi','is_followup_period','patient_id'),
  random_effect = '(1|patient_id)',
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05,
  return_model = TRUE
)
plot_model(actigraph_sf36_out1_mdl, type = 're', title = 'random effects - actigraph - sf36 total')

## haq
actigraph_haq_out1_mdl <- predict_delta_comps(
  dataf = df_actigraph_daily_comp_awake[
    ,
    c('haq', 'sedentary_proportion','lipa_proportion','mvpa_proportion','is_followup_period','bmi','patient_id')
    ],
  y = 'haq',
  comps = c('sedentary_proportion','lipa_proportion','mvpa_proportion'),
  # comps = c('lipa_proportion','sedentary_proportion','mvpa_proportion'), 
  # comps = c('mvpa_proportion','lipa_proportion','sedentary_proportion'), 
  ## note: we can switch the order of covariate to get information of relative amount of time spent in 1st behavior
  covars = c('bmi','is_followup_period','patient_id'),
  random_effect = '(1|patient_id)',
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05,
  return_model = TRUE
)
plot_model(actigraph_haq_out1_mdl, type = 're', title = 'random effects - actigraph - haq')

set_theme(
  base = theme_classic(),
  axis.tickslen = 0, # hides tick marks
  axis.textsize = .8,
  geom.label.size = 3.5
)

################################################################################################################################################
## 03. OBSERVE MORE COMPLETE MODEL SUMMARY
## 03-3. PLOT THE RANDOM EFFECTS PER PATIENT
################################################################################################################################################

p_ap_sf36 <- plot_model(activpal_sf36_out1_mdl, type = 're', title = " ")
p_ap_haq <- plot_model(activpal_haq_out1_mdl, type = 're', title = " ")
p_ag_sf36 <- plot_model(actigraph_sf36_out1_mdl, type = 're', title = " ")
p_ag_haq <- plot_model(actigraph_haq_out1_mdl, type = 're', title = " ")

pl <- list(p_ag_sf36, p_ag_haq, p_ap_sf36, p_ap_haq)

## reference: https://stackoverflow.com/questions/45473843/put-row-and-column-titles-using-grid-arrange-in-r
library(grid)
library(gridExtra)
N <- length(pl)
nr <- 2
nc <- 2

combine <- rbind(tableGrob(t(c('sf36_total','haq')), theme = ttheme_minimal(), rows = ""), 
                 cbind(tableGrob(c('actigraph','activpal'), theme = ttheme_minimal()), 
                       arrangeGrob(grobs = pl),  size = "last"), size = "last")
grid.newpage()
grid.draw(combine)