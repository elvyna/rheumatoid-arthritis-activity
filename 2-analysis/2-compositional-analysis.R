################################################################################################################################################
## 00. LOAD LIBRARIES
## ALSO SET WORKING DIRECTORY
################################################################################################################################################

library(dplyr)
library(tidyr)
library(reshape2)
library(ggthemes)
library(ggplot2)
library(lubridate)
library(MASS)
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

ACTIVPAL_FILE <- "../dataset/activpal_hourly_cleaned_duration_without_slnw_with_meta_patients33_with_activity_proportion.csv"
df_activpal_hourly_with_prop <- read.csv(ACTIVPAL_FILE)

df_activpal_hourly_with_prop$date <- as.Date(df_activpal_hourly_with_prop$datetime_hour, '%Y-%m-%d')

################################################################################################################################################
## 01. DATA IMPORT AND PREPARATION
## 01-2. ACTIGRAPH
################################################################################################################################################

ACTIGRAPH_HOURLY_FILE <- '../dataset/actigraph_preprocessed_choi_valid_patients53_hourly_with_activity_proportion.csv'
df_actigraph_hourly_with_prop <- read.csv(ACTIGRAPH_HOURLY_FILE)

################################################################################################################################################
## 02. VISUALISE COMPOSITIONAL DATA IN 3-DIMENSIONS USING TERNARY PLOT
## 02-1. ACTIVPAL
################################################################################################################################################

library(ggtern)

### USE AVG OF DAILY DATA
df_activpal_daily_comp_all <- df_activpal_hourly_with_prop %>% 
  group_by(patient_id, visit_info, date, age, bmi, sf36_total, haq, is_awake_time) %>% 
  summarise(
    sitting_minute = sum(awake_sitting_duration_second_clean) / 60,
    standing_minute = sum(awake_standing_duration_second_clean) / 60,
    stepping_minute = sum(awake_stepping_duration_second_clean) / 60
  )
activpal_comp_temp <- as.data.frame(
  acomp(df_activpal_daily_comp_all[, c("sitting_minute", "standing_minute","stepping_minute")])
)
colnames(activpal_comp_temp) <- c('sitting_proportion','standing_proportion','stepping_proportion')
df_activpal_daily_comp_all <- dplyr::bind_cols(df_activpal_daily_comp_all, activpal_comp_temp)

df_activpal_daily_comp_all <- df_activpal_daily_comp_all %>% 
  group_by(patient_id, visit_info, age, bmi, sf36_total, haq, is_awake_time) %>% 
  summarise(
    sitting_proportion = mean(sitting_proportion),
    standing_proportion = mean(standing_proportion),
    stepping_proportion = mean(stepping_proportion)
  )
df_activpal_daily_comp_all <- as.data.frame(df_activpal_daily_comp_all)
df_activpal_daily_comp_all$time <- ifelse(df_activpal_daily_comp_all$is_awake_time == 0, 'Night', 'Day')

ggtern(
  data = df_activpal_daily_comp_all, 
  aes(
    x = sitting_proportion, 
    y = standing_proportion, 
    z = stepping_proportion, 
    col = visit_info
  )
) + 
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("#FF69B4", "#7FFF00"),
                     name = "Period",
                     breaks = c("Baseline", "Follow up"),
                     labels = c("Baseline", "Follow up")
  ) +
  facet_grid(. ~ time) +
  ggtitle('Activpal') +
  labs(x = 'sit', y = 'stand', z = 'step') +
  theme_bluelight()

################################################################################################################################################
## 02. VISUALISE COMPOSITIONAL DATA IN 3-DIMENSIONS USING TERNARY PLOT
## 02-2. ACTIGRAPH
################################################################################################################################################

### USE AVG OF DAILY DATA
df_actigraph_daily_comp_all <- df_actigraph_hourly_with_prop %>% 
  group_by(patient_id, directory_name, date, age, bmi, sf36_total, haq, is_awake_time) %>% 
  summarise(
    sedentary_behaviour = sum(sedentary_behaviour),
    light_intensity = sum(light_intensity),
    mvpa = sum(moderate + vigorous)
  )
actigraph_comp_temp <- as.data.frame(
  acomp(df_actigraph_daily_comp_all[, c("sedentary_behaviour", "light_intensity","mvpa")])
)
colnames(actigraph_comp_temp) <- c('sedentary_proportion','lipa_proportion','mvpa_proportion')
df_actigraph_daily_comp_all <- dplyr::bind_cols(df_actigraph_daily_comp_all, actigraph_comp_temp)

df_actigraph_daily_comp_all <- df_actigraph_daily_comp_all %>% 
  group_by(patient_id, directory_name, age, bmi, sf36_total, haq, is_awake_time) %>% 
  summarise(
    sedentary_proportion = mean(sedentary_proportion),
    lipa_proportion = mean(lipa_proportion),
    mvpa_proportion = mean(mvpa_proportion)
  )
df_actigraph_daily_comp_all <- as.data.frame(df_actigraph_daily_comp_all)
df_actigraph_daily_comp_all$time <- ifelse(df_actigraph_daily_comp_all$is_awake_time == 0, 'Night', 'Day')

library(ggtern)
ggtern(
  data = df_actigraph_daily_comp_all, 
  aes(
    x = sedentary_proportion, 
    y = lipa_proportion, 
    z = mvpa_proportion, 
    col = directory_name
  )
) + 
  geom_point(alpha = 0.8) +
  scale_color_manual(values=c("#FF69B4", "#7FFF00"),
                     name="Period",
                     breaks=c("Baseline", "Follow up"),
                     labels=c("Baseline", "Follow up")
  ) +
  facet_grid(.~time) +
  ggtitle('Actigraph') +
  labs(x = 'SB', y = 'LIPA', z = 'MVPA') +
  theme_bluelight()

################################################################################################################################################
## 02. VISUALISE COMPOSITIONAL DATA IN 3-DIMENSIONS USING TERNARY PLOT
## 02-3. SHOW DAILY AVERAGE COMPOSITION OF ACTIGRAPH AND ACTIVPAL IN ONE PLOT OBJECT
################################################################################################################################################

library(ggtern)

p1 <- ggtern(
  data = df_actigraph_daily_comp_all, 
  aes(
    x = sedentary_proportion, 
    y = lipa_proportion, 
    z = mvpa_proportion, 
    col = directory_name
  )
) + 
  geom_point(alpha = 0.8) +
  scale_color_manual(values=c("#FF69B4", "#7FFF00"),
                     name="Period",
                     breaks=c("Baseline", "Follow up"),
                     labels=c("Baseline", "Follow up")
  ) +
  facet_grid(.~time) +
  ggtitle('actigraph') +
  labs(x = 'SB', y = 'LIPA', z = 'MVPA') +
  theme_bluelight() +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  )

p2 <- ggtern(
  data = df_activpal_daily_comp_all, 
  aes(
    x = sitting_proportion, 
    y = standing_proportion, 
    z = stepping_proportion, 
    col = visit_info
  )
) + 
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("#FF69B4", "#7FFF00"),
                     name = "Period",
                     breaks = c("Baseline", "Follow up"),
                     labels = c("Baseline", "Follow up")
  ) +
  facet_grid(. ~ time) +
  ggtitle('activpal') +
  labs(x = 'sit', y = 'stand', z = 'step') +
  theme_bluelight() +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  )

grid.arrange(p1, p2)

unloadNamespace('ggtern')

################################################################################################################################################
## 03. HYPOTHESIS TEST TO COMPARE DURATION (IN COMPOSITION VALUE)
## 03-1. ACTIVPAL
################################################################################################################################################

df_activpal_daily_comp_awake <- df_activpal_hourly_with_prop %>% 
  filter(is_awake_time == 1) %>% 
  group_by(patient_id, 
           visit_info, 
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

df_activpal_daily_comp_awake <- df_activpal_daily_comp_awake %>% 
  group_by(
    patient_id, 
    visit_info, 
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

### STATISTICAL TEST: Multivariate (input: avg composition)
## cannot perform manova on the daily avg proportion. use shrinkage estimate
library(mvabund)
ap_values <- as.matrix(df_activpal_daily_comp_awake[, c('sitting_proportion','standing_proportion','stepping_proportion')])
# ap_manova_daily <- manova(ap_values ~ df_activpal_daily_comp_awake$visit_info)
ap_manova_daily <- manylm(ap_values ~ df_activpal_daily_comp_awake$visit_info, cor.type = 'shrink')
summary(ap_manova_daily)

################################################################################################################################################
## 03. HYPOTHESIS TEST TO COMPARE DURATION (IN COMPOSITION VALUE)
## 03-2. ACTIGRAPH
################################################################################################################################################

### USE AVG OF DAILY DATA
df_actigraph_daily_comp_awake <- df_actigraph_hourly_with_prop %>% 
  filter(is_awake_time == 1) %>%
  group_by(patient_id, 
           directory_name, 
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

df_actigraph_daily_comp_awake <- df_actigraph_daily_comp_awake %>% 
  group_by(
    patient_id, 
    directory_name, 
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

### STATISTICAL TEST: Multivariate (input: avg composition)
## cannot perform manova on the daily avg proportion. use shrinkage estimate
library(mvabund)
ag_values <- as.matrix(df_actigraph_daily_comp_awake[, c('sedentary_proportion', 'lipa_proportion', 'mvpa_proportion')])
# ag_manova_daily <- manova(ag_values ~ df_actigraph_daily_comp_awake$directory_name)
ag_manova_daily <- manylm(ag_values ~ df_actigraph_daily_comp_awake$directory_name, cor.type = 'shrink')
summary(ag_manova_daily)

################################################################################################################################################
## 04. COMPOSITIONAL MULTIPLE REGRESSION
## 04-1. ACTIVPAL
## 04-1A. SF36_TOTAL
################################################################################################################################################

library(deltacomp)

## Predict sf36_total
## Baseline data
ap_out_sf36_pre <- predict_delta_comps(
  dataf = df_activpal_daily_comp_awake[
    (df_activpal_daily_comp_awake$visit_info == 'Baseline'),
    c('sf36_total', 'sitting_proportion', 'standing_proportion', 'stepping_proportion', 'age', 'bmi')
    ],
  y = 'sf36_total',
  comps = c('sitting_proportion','standing_proportion','stepping_proportion'), 
  ## note: we can switch the order of covariate to get information of relative amount of time spent in 1st behavior
  # covars = c('bmi'),
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

plot_delta_comp(
  ap_out_sf36_pre,
  comp_total = 24 * 60,
  units_lab = "min"
)

## Follow-up data
ap_out_sf36_post <- predict_delta_comps(
  dataf = df_activpal_daily_comp_awake[
    (df_activpal_daily_comp_awake$visit_info == 'Follow up'),
    c('sf36_total', 'sitting_proportion', 'standing_proportion', 'stepping_proportion', 'age', 'bmi')
    ],
  y = 'sf36_total',
  comps = c('sitting_proportion','standing_proportion','stepping_proportion'), 
  # covars = c('bmi'),
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

plot_delta_comp(
  ap_out_sf36_post, 
  comp_total = 24 * 60,
  units_lab = "min" 
)

################################################################################################################################################
## 04. COMPOSITIONAL MULTIPLE REGRESSION
## 04-1. ACTIVPAL
## 04-1B. HAQ
################################################################################################################################################

## Predict haq
## Baseline data
ap_out_haq_pre <- predict_delta_comps(
  dataf = df_activpal_daily_comp_awake[
    (df_activpal_daily_comp_awake$visit_info == 'Baseline'),
    c('haq', 'sitting_proportion', 'standing_proportion', 'stepping_proportion', 'age', 'bmi')
    ],
  y = 'haq',
  comps = c('sitting_proportion','standing_proportion','stepping_proportion'),
  # covars = c('bmi'),
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

plot_delta_comp(
  ap_out_haq_pre,
  comp_total = 24 * 60, 
  units_lab = "min"
)

## Follow-up data
ap_out_haq_post <- predict_delta_comps(
  dataf = df_activpal_daily_comp_awake[
    (df_activpal_daily_comp_awake$visit_info == 'Follow up'),
    c('haq', 'sitting_proportion', 'standing_proportion', 'stepping_proportion', 'age', 'bmi')
    ],
  y = 'haq',
  comps = c('sitting_proportion','standing_proportion','stepping_proportion'),
  # covars = c('bmi'),
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

plot_delta_comp(
  ap_out_haq_post,
  comp_total = 24 * 60,
  units_lab = "min"
)

################################################################################################################################################
## 04. COMPOSITIONAL MULTIPLE REGRESSION
## 04-2. ACTIGRAPH
## 04-2A. SF36_TOTAL
################################################################################################################################################

library(deltacomp)

## baseline
ag_out_sf36_pre <- predict_delta_comps(
  dataf = df_actigraph_daily_comp_awake[
    (df_actigraph_daily_comp_awake$directory_name == 'Baseline'),
    c('sf36_total','sedentary_proportion','lipa_proportion','mvpa_proportion','age','bmi')
    ],
  y = 'sf36_total',
  comps = c("sedentary_proportion", "lipa_proportion", "mvpa_proportion"),
  # covars = c('bmi'),
  comparisons = 'prop-realloc',
  # deltas = c(0, 10, 20)/(24 * 60),
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

plot_delta_comp(
  ag_out_sf36_pre,
  comp_total = 24 * 60,
  units_lab = "min"
)

## follow-up
ag_out_sf36_post <- predict_delta_comps(
  dataf = df_actigraph_daily_comp_awake[
    (df_actigraph_daily_comp_awake$directory_name == 'Follow up'),
    c('sf36_total','sedentary_proportion','lipa_proportion','mvpa_proportion','age','bmi')
    ],
  y = 'sf36_total',
  comps = c("sedentary_proportion", "lipa_proportion", "mvpa_proportion"),
  # covars = c('bmi'),
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

plot_delta_comp(
  ag_out_sf36_post,
  comp_total = 24 * 60,
  units_lab = "min"
)

################################################################################################################################################
## 04. COMPOSITIONAL MULTIPLE REGRESSION
## 04-2. ACTIGRAPH
## 04-2B. HAQ
################################################################################################################################################

## predict haq
## baseline
ag_out_haq_pre <- predict_delta_comps(
  dataf = df_actigraph_daily_comp_awake[
    (df_actigraph_daily_comp_awake$directory_name == 'Baseline'),
    c('haq','sedentary_proportion','lipa_proportion','mvpa_proportion','age','bmi')
    ],
  y = 'haq',
  comps = c("sedentary_proportion", "lipa_proportion", "mvpa_proportion"),
  # covars = c('bmi'),
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

plot_delta_comp(
  ag_out_haq_pre,
  comp_total = 24 * 60,
  units_lab = "min"
)

## follow-up
ag_out_haq_post <- predict_delta_comps(
  dataf = df_actigraph_daily_comp_awake[
    (df_actigraph_daily_comp_awake$directory_name == 'Follow up'),
    c('haq','sedentary_proportion','lipa_proportion','mvpa_proportion','age','bmi')
    ],
  y = 'haq',
  comps = c("sedentary_proportion", "lipa_proportion", "mvpa_proportion"),
  # covars = c('bmi'),
  comparisons = 'prop-realloc',
  deltas =  seq(-60, 60, by = 10) / (24 * 60),
  alpha = 0.05
)

plot_delta_comp(
  ag_out_haq_post,
  comp_total = 24 * 60,
  units_lab = "min"
)

unloadNamespace('deltacomp')