################################################################################################################################################
## 00. LOAD LIBRARIES
## ALSO SET WORKING DIRECTORY
################################################################################################################################################

library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggfortify)

DIRECTORY <- 'PUT-THE-WORKING-DIRECTORY'
setwd(DIRECTORY)

## define colorblind palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

################################################################################################################################################
## 01. DATA IMPORT AND PREPARATION
################################################################################################################################################

METADATA_FILE <- '../dataset/patient_metadata_trimmed.csv'
df_meta <- read.csv(METADATA_FILE) %>% 
  as_tibble()

df_meta$visit <- ifelse(
  df_meta$visit == 1, 
  'Baseline', 
  ifelse(
    df_meta$visit == 2,
    'Follow up',
    NA
  )
)

ACTIGRAPH_FILE <- "../dataset/actigraph_preprocessed_choi_valid_patients53.csv"
df_actigraph <- read.csv(ACTIGRAPH_FILE)
VALID_PATIENT_ID <- unique(df_actigraph$patient_id)
rm(df_actigraph)

df_meta_health <- df_meta %>% 
  filter(!is.na(haq) & !is.na(sf36_total) & patient_id %in% VALID_PATIENT_ID) %>% 
  dplyr::select(c(starts_with("sf"), 'haq'))

colnames(df_meta_health) <- c(
  'sf_physical_functioning',
  'sf_role_physical',
  'sf_bodily_pain',
  'sf_social_functioning',
  'sf_mental_health',
  'sf_role_emotional',
  'sf_vitality',
  'sf_general_health',
  'sf36_total',
  'haq'
)

################################################################################################################################################
## 02. PAIRED T-TEST
## 02-1. SF36_TOTAL
################################################################################################################################################

PATIENT_ID_COMPLETE_META <- df_meta %>% 
  filter(!is.na(haq) & !is.na(sf36_total) & patient_id %in% VALID_PATIENT_ID) %>% 
  group_by(patient_id) %>% 
  count() %>% 
  filter(n == 2) %>% 
  select(patient_id)
PATIENT_ID_COMPLETE_META <- PATIENT_ID_COMPLETE_META$patient_id

df_meta_valid <- df_meta %>% 
  filter(!is.na(haq) & 
           !is.na(sf36_total) & 
           patient_id %in% VALID_PATIENT_ID &
           patient_id %in% PATIENT_ID_COMPLETE_META
  ) %>% 
  dplyr::select('patient_id', 'visit', 'sf36_total', 'haq')

t.test(
  x = subset(df_meta_valid, visit == 'Baseline')$sf36_total, 
  y = subset(df_meta_valid, visit == 'Follow up')$sf36_total,
  paired = TRUE, 
  alternative = "two.sided",
  var.equal = FALSE
)

plot_sf36 <- ggplot(df_meta_valid) +
  geom_density(
    aes(
      sf36_total,
      group = visit,
      fill = visit
    ),
    alpha = 0.5
  ) +
  scale_x_continuous(
    name = 'SF36 total',
    limits = c(0,100),
    breaks = seq(0,100,10)
  ) +
  scale_fill_manual(
    values = cbbPalette
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, colour = "grey80")
  ) +
  labs(fill = 'Period')

################################################################################################################################################
## 02. PAIRED T-TEST
## 02-2. HAQ
################################################################################################################################################

t.test(
  x = subset(df_meta_valid, visit == 'Baseline')$haq,
  y = subset(df_meta_valid, visit == 'Follow up')$haq,
  paired = TRUE, 
  alternative = "two.sided",
  var.equal = FALSE
)

plot_haq <- ggplot(df_meta_valid) +
  geom_density(
    aes(
      haq,
      group = visit,
      fill = visit
    ),
    alpha = 0.5
  ) +
  scale_x_continuous(
    name = 'HAQ',
    limits = c(0,3),
    breaks = seq(0,3,.25)
  ) +
  scale_fill_manual(
    values = cbbPalette
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, colour = "grey80")
  ) +
  labs(fill = 'Period')

library(gridExtra)
grid.arrange(plot_sf36, plot_haq)

################################################################################################################################################
## 03. CHECK CORRELATION
################################################################################################################################################

df_meta_health %>% 
  cor() %>% 
  as_tibble() %>% 
  pivot_longer(cols = colnames(df_meta_health))

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)] <- NA
  return (cormat)
}

health_cor <- round(cor(df_meta_health), 2)
health_cor_melt <- melt(get_upper_tri(health_cor))
colnames(health_cor_melt) <- c('var1', 'var2', 'correlation')
ggplot(health_cor_melt,
       aes(
         x = var1,
         y = var2,
         fill = correlation
       )
) +
  geom_tile(
    colour = 'white'
  ) +
  geom_text(
    aes(
      x = var1, 
      y = var2, 
      label = correlation
    ), 
    color = "black", 
    size = 4
  ) +
  scale_x_discrete(
    name = ''
  ) +
  scale_y_discrete(
    name = ''
  ) +
  scale_fill_gradient2(low = "white", high = "blue", mid = "cyan",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Pearson\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1,
                                   size = 10, hjust = 1)) +
  coord_fixed()

################################################################################################################################################
## 04. CONSTRUCT PRINCIPAL COMPONENTS
################################################################################################################################################

pc <- prcomp(df_meta_health, scale = TRUE)
total_variation <- sum(pc$sdev^2)
(pc$sdev^2) / total_variation

df_pc <- data.frame(
  component = 1:ncol(pc$rotation),
  variance = pc$sdev^2,
  variance_pct = (pc$sdev^2) / total_variation
)

## prettier
ggplot(data = df_pc) +
  geom_line(
    aes(
      x = component,
      y = variance_pct
    )
  ) +
  scale_x_continuous(
    name = 'n-th component',
    limits = c(1,10),
    breaks = seq(1,10,1)
  ) +
  scale_y_continuous(
    name = '% of explained variance',
    labels = scales::percent
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, colour = "grey80")
  )

print(paste(
  "Variation represented using 2 PCs: ",
  round(100 * sum((pc$sdev^2)[1:2]) / total_variation, 2),
  "%",
  sep = ""
))    

print(paste(
  "Variation represented using 3 PCs: ",
  round(100 * sum((pc$sdev^2)[1:3]) / total_variation, 2),
  "%",
  sep = ""
))

pc_rotated <- varimax(pc$rotation[, 1:2])
pc_rotated

## prettier; no overlapping variable names
library(ggbiplot)
ggbiplot(
  pc,
  alpha = 0.5,
  varname.adjust = 1.25,
  varname.size = 3
) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, colour = "grey80")
  ) +
  scale_x_continuous(
    name = 'PC1 (51.9% explained variance)'
  ) +
  scale_y_continuous(
    name = 'PC2 (11.1% explained variance)'
  )