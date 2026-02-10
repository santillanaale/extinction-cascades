rm(list=ls())
library(bipartite)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(lme4)
library(broom.mixed)
library(lmerTest)
library(ggeffects)
library(car)
set.seed(123)
setwd("~/")
source("lab_paths.R")
local.path

# ---- Load Networks ----
setwd(file.path(local.path, "skyIslands"))
load("data/networks/YearSR_PlantPollinator_Bees.Rdata") #Site x Year x Survey Round

## --- Make net key ----
net_meta <- tibble(net_id = names(nets)) %>%
  separate(net_id, into = c("Site", "Year", "SampleRound"), sep = "\\.") %>%
  mutate(
    Year = as.numeric(Year),
    SampleRound = as.numeric(SampleRound)
  )

bad_net <- "PL.2017.1"

keep <- names(nets) != bad_net

nets <- nets[keep]
net_meta <- net_meta[keep, ]

# ---- TCM: Combining all scenarios into one loop ----
## ---- Calculating Robustness ----
### ---- Load modified extinction functions ----
setwd("~/")
setwd(file.path(local.path, "extinction-cascades"))

source("modified_extinction.R")
source("modified_second_extinct.R")

### ---- Unified helper function ----
run_robustness <- function(web, scenario, nrep = 100) {
  
  stopifnot(scenario$method %in% c("abundance", "degree"))
  
  ext <- modified.second.extinct(
    web = web,
    participant = "lower",      # plants
    method = scenario$method,
    reverse = scenario$reverse,
    nrep = nrep,
    details = FALSE
  )
  
  robustness(ext)
}

### ---- Define scenarios ----
scenarios <- tibble(
  scenario = c(
    "abundance_low", # remove rare first
    "abundance_high", # remove dominant first
    "degree_low", # remove specialists first
    "degree_high" # remove generalists first
  ),
  method = c(
    "abundance",
    "abundance",
    "degree",
    "degree"
  ),
  reverse = c(
    FALSE,
    TRUE,
    FALSE,
    TRUE
  )
)

### ---- Robustness Loop (Site × Year × SampleRound × scenario × Robustness----
robustness_results <- map_dfr(seq_len(nrow(scenarios)), function(s) {
  
  scenario <- scenarios[s, ]
  
  message("Running scenario: ", scenario$scenario)
  
  map_dfr(seq_along(nets), function(i) {
    
    web <- nets[[i]]
    
    rob <- run_robustness(
      web = web,
      scenario = scenario,
      nrep = 100
    )
    
    tibble(
      Site = net_meta$Site[i],
      Year = net_meta$Year[i],
      SampleRound = net_meta$SampleRound[i],
      scenario = scenario$scenario,
      Robustness = rob
    )
  })
})

## ---- Merge Climate ----
setwd("~/")
setwd(file.path(local.path, "skyIslands_saved"))

climate <- read.csv("data/relational/original/climate.csv")

tcm_climate <- robustness_results %>%
  left_join(
    climate,
    by = c("Site", "Year", "SampleRound")
  )

setwd("~/")
setwd(file.path(local.path, "extinction-cascades"))
save(
  robustness_results,
  tcm_climate,
  file = "analysis/network/robustness_metrics/YearSR_robustness_climate.Rdata"
)

## ---- Modeling ----
### ---- Standardize predictors ----
# tcm_climate <- tcm_climate %>%
#   mutate(across(
#     c(SpringPrecip, CumulativePrecip, RoundPrecip,
#       SpringTmean, CumulativeTmean, RoundTmean, 
#       SpringTmeanAnom, RoundTmeanAnom, CumulativeTmeanAnom),
#     ~ scale(.)[,1],
#     .names = "z_{.col}"
#   ))

### ---- Core models ----
#### ---- Precipitation only ----
# mod_precip <- lmer(
#   Robustness ~ 
#     #scale(SpringPrecip) +
#     scale(log(CumulativePrecip)) +
#     #scale(RoundPrecip) +
#     #scenario +
#     (1 | Site) +
#     (1 | Year)
#   ,data = tcm_climate[tcm_climate$scenario=="abundance_high",]
# )
# 
# summary(mod_precip)

#### ---- Temperature only (absolute, TMean) ----
# mod_tmean <- lmer(
#   Robustness ~ 
#     # scale(SpringTmean) +
#     scale(log(CumulativeTmean)) +
#     # scale(RoundTmean) +
#     # scenario +
#     (1 | Site) +
#     (1 | Year),
#   data = tcm_climate[tcm_climate$scenario=="abundance_high",]
# )
# 
# summary(mod_tmean)

#### ---- Temperature anomalies ----
# mod_tmean_anom <- lmer(
#   Robustness ~ 
#     # scale(SpringTmeanAnom) +
#     scale(CumulativeTmeanAnom) +
#     # scale(RoundTmeanAnom) +
#     # scenario +
#     (1 | Site) +
#     (1 | Year),
#   data = tcm_climate[tcm_climate$scenario=="abundance_high",]
# )
# 
# summary(mod_tmean_anom)

#### ---- Full Climate Model (Cumulative Precip and TmeanAnom) ----
mod_tcm_all <- lmer(
  Robustness ~ 
    scale(CumulativePrecip) +
    scale(CumulativeTmeanAnom) +
    scenario +
    (1 | Site) +
    (1 | Year),
  data = tcm_climate
)

summary(mod_tcm_all)

vif(mod_tcm_all)

##### ---- Model-based predictions ----
tcm_preds_precip <- ggpredict(
  mod_tcm_all,
  terms = c("CumulativePrecip [all]", "scenario")
)

tcm_preds_temp <- ggpredict(
  mod_tcm_all,
  terms = c("CumulativeTmeanAnom [all]", "scenario")
)

##### ---- Plot TCM Effects ----
p_tcm_precip <- ggplot(tcm_preds_precip,
                       aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(
    x = "Cumulative precipitation",
    y = "TCM robustness",
    color = "Extinction scenario",
    fill = "Extinction scenario"
  )

p_tcm_temp <- ggplot(tcm_preds_temp,
                     aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(
    x = "Cumulative temperature anomaly",
    y = NULL,
    color = "Extinction scenario",
    fill = "Extinction scenario"
  )

(p_tcm_precip | p_tcm_temp) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(
  "analysis/network/figures/TCM_Prediction_Robustness_by_CumulativeClimate_RoundScale.pdf",
  width = 12, height = 4
)

# Interpretation framework (this is paper-ready)
# 
# Spring terms → legacy / phenological setup
# 
# Cumulative terms → seasonal resource tracking
# 
# Round terms → short-term sensitivity
# 
# Scenario effects → robustness mechanism dependence
# 
# Anomalies vs absolute → physiological vs climatological drivers

# ---- SCM ----
## ---- Calculating Robustness ----
source("SCM/netcascade_mod.R")
source("SCM/scm_robustness_helpers.R")

### ---- Define SCM scenarios ----
scm_scenarios <- tibble(
  scenario = c(
    "scm_random_plant",
    "scm_dominant_plant"
  ),
  TargetGuild = c("rows", "rows"),
  TargetSpecies = c("random_binary", "random_abundance")
)

### ---- Run SCM robustness in the same Site × Year x Sample Round loop ----
scm_results <- map_dfr(seq_len(nrow(scm_scenarios)), function(s) {
  
  scenario <- scm_scenarios[s, ]
  
  message("Running SCM scenario: ", scenario$scenario)
  
  map_dfr(seq_along(nets), function(i) {
    
    web <- as.matrix(nets[[i]])
    print(names(nets)[i])
    
    R_rows <- rep(0.5, nrow(web))  # or sample like VA2015
    R_cols <- rep(0.5, ncol(web))
    
    out <- run_scm_robustness(
      web = web,
      R_val_rows = R_rows,
      R_val_cols = R_cols,
      TargetGuild = scenario$TargetGuild,
      TargetSpecies = scenario$TargetSpecies,
      nrep = 100
    )
    
    tibble(
      Site = net_meta$Site[i],
      Year = net_meta$Year[i],
      SampleRound = net_meta$SampleRound[i],
      scenario = scenario$scenario,
      prop_cascade = out$prop_cascade,
      SCM_Robustness = out$SCM_robustness
    )
  })
})

## ---- Merge with Climate ----
scm_climate <- scm_results %>%
  left_join(climate, by = c("Site", "Year", "SampleRound"))

## ---- Modeling ----
### ---- Standardize predictors ----
# scm_climate <- scm_climate %>%
#   mutate(across(
#     c(SpringPrecip, CumulativePrecip, RoundPrecip,
#       SpringTmean, CumulativeTmean, RoundTmean, 
#       SpringTmeanAnom, RoundTmeanAnom, CumulativeTmeanAnom),
#     ~ scale(.)[, 1],
#     .names = "z_{.col}"
#   ))

### ---- Core Models ---- 
#### ---- Precipitation only ----
# mod_scm_precip <- lmer(
#   SCM_Robustness ~ 
#     # z_SpringPrecip +
#     scale(CumulativePrecip) *
#     # z_RoundPrecip +
#     scenario +
#     (1 | Site) +
#     (1 | Year),
#   data = scm_climate
# )
# 
# summary(mod_scm_precip)

#### ---- Temperature only (absolute, TMean) ----
# mod_scm_tmean <- lmer(
#   SCM_Robustness ~ 
#     # z_SpringTmean +
#     scale(CumulativeTmean) *
#     # z_RoundTmean +
#     scenario +
#     (1 | Site) +
#     (1 | Year),
#   data = scm_climate
# )
# 
# summary(mod_scm_tmean)

#### ---- Temperature anomalies ----
# mod_scm_tmean_anom <- lmer(
#   SCM_Robustness ~
#     # z_SpringTmeanAnom +
#     scale(CumulativeTmeanAnom) *
#     # z_RoundTmeanAnom +
#     scenario +
#     (1 | Site) +
#     (1 | Year),
#   data = scm_climate
# )
# 
# summary(mod_scm_tmean_anom)

#### ---- Full Climate Model (Cumulative Precip and TmeanAnom) ----
mod_scm_all <- lmer(
  SCM_Robustness ~ 
    scale(CumulativePrecip) *
    scenario +
    scale(CumulativeTmeanAnom) *
    scenario +
    (1 | Site) +
    (1 | Year),
  data = scm_climate
)

summary(mod_scm_all)

##### ---- Model-based predictions ----
scm_preds_precip <- ggpredict(
  mod_scm_all,
  terms = c("CumulativePrecip [all]", "scenario")
)

scm_preds_temp <- ggpredict(
  mod_scm_all,
  terms = c("CumulativeTmeanAnom [all]", "scenario")
)

##### ---- Plot SCM Effects ----
p_scm_precip <- ggplot(scm_preds_precip,
                       aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(
    x = "Cumulative precipitation",
    y = "SCM robustness",
    color = "SCM scenario",
    fill = "SCM scenario"
  )

p_scm_temp <- ggplot(scm_preds_temp,
                     aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(
    x = "Cumulative temperature anomaly",
    y = NULL
  )

(p_scm_precip | p_scm_temp) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(
  "analysis/network/figures/SCM_Prediction_Robustness_by_CumulativeClimate_RoundScale.pdf",
  width = 12, height = 4
)

# ---- Plotting ----
## ---- TCM ----
### ---- Precipitation ----
# p1 <- ggplot(tcm_climate,
#              aes(x = SpringPrecip, y = Robustness)) +
#   geom_point(aes(color = Site), alpha = 0.6) +
#   geom_smooth(method = "lm", se = FALSE) +
#   theme_minimal() +
#   labs(title = "Spring Precipitation",
#        x = "Spring precipitation",
#        y = "Robustness")
# 
p1 <- ggplot(tcm_climate,
             aes(x = CumulativePrecip, y = Robustness)) +
  geom_point(aes(color = Site), alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Cumulative Precipitation",
       x = "Cumulative precipitation")
# 
# p3 <- ggplot(tcm_climate,
#              aes(x = RoundPrecip, y = Robustness)) +
#   geom_point(aes(color = Site), alpha = 0.6) +
#   geom_smooth(method = "lm", se = FALSE) +
#   theme_minimal() +
#   labs(title = "Round Precipitation",
#        x = "Round precipitation")
# 
# (p1 | p2 | p3)
# 
# ggsave(
#   "analysis/network/figures/TCM_Robustness_by_Precip_RoundScale.pdf",
#   width = 12, height = 4
# )

ggsave(
  "analysis/network/figures/TCM_Robustness_by_CumulativePrecip_RoundScale.pdf",
  width = 6, height = 4
)

### ---- Mean Temperature (absolute) ----
# p1 <- ggplot(tcm_climate,
#              aes(x = SpringTmean, y = Robustness)) +
#   geom_point(aes(color = Site), alpha = 0.6) +
#   geom_smooth(method = "lm", se = FALSE) +
#   theme_minimal() +
#   labs(title = "Spring Mean Temperature")

# p2 <- ggplot(tcm_climate,
#              aes(x = CumulativeTmean, y = Robustness)) +
#   geom_point(aes(color = Site), alpha = 0.6) +
#   geom_smooth(method = "lm", se = FALSE) +
#   theme_minimal() +
#   labs(title = "Cumulative Mean Temperature")

# p3 <- ggplot(tcm_climate,
#              aes(x = RoundTmean, y = Robustness)) +
#   geom_point(aes(color = Site), alpha = 0.6) +
#   geom_smooth(method = "lm", se = FALSE) +
#   theme_minimal() +
#   labs(title = "Round Mean Temperature")
# 
# (p1 | p2 | p3)
# 
# ggsave(
#   "analysis/network/figures/TCM_Robustness_by_Temperature_RoundScale.pdf",
#   width = 12, height = 4
# )

# ggsave(
#   "analysis/network/figures/TCM_Robustness_by_CumulativeTmean_RoundScale.pdf",
#   width = 6, height = 4
# )

### ---- Mean Temperature (Anomaly) ----
# p1 <- ggplot(tcm_climate,
#              aes(x = SpringTmeanAnom, y = Robustness)) +
#   geom_point(aes(color = Site), alpha = 0.6) +
#   geom_smooth(method = "lm", se = FALSE) +
#   theme_minimal() +
#   labs(title = "Spring Mean Temperature")

p2 <- ggplot(tcm_climate,
             aes(x = CumulativeTmeanAnom, y = Robustness)) +
  geom_point(aes(color = Site), alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Cumulative Mean Temperature Anomaly")

# p3 <- ggplot(tcm_climate,
#              aes(x = RoundTmeanAnom, y = Robustness)) +
#   geom_point(aes(color = Site), alpha = 0.6) +
#   geom_smooth(method = "lm", se = FALSE) +
#   theme_minimal() +
#   labs(title = "Round Mean Temperature")
# 
# (p1 | p2 | p3)
# 
# ggsave(
#   "analysis/network/figures/TCM_Robustness_by_TemperatureAnomaly_RoundScale.pdf",
#   width = 12, height = 4
# )

ggsave(
  "analysis/network/figures/TCM_Robustness_by_CumulativeTmeanAnom_RoundScale.pdf",
  width = 6, height = 4
)

### ---- Combining Robustness by Climate ----

(p1 | p2)

ggsave(
  "analysis/network/figures/TCM_Robustness_by_CumulativeClimate_RoundScale.pdf",
  width = 12, height = 4
)

### ---- By extinction scenario to show whether extinction order changes climate sensitivity ----
#### ---- Precip ----
# ggplot(tcm_climate,
#        aes(x = SpringPrecip, y = Robustness)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm", se = FALSE) +
#   facet_wrap(~ scenario) +
#   theme_minimal()

p1 <- ggplot(tcm_climate,
       aes(x = CumulativePrecip, y = Robustness)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ scenario) +
  theme_minimal()

ggsave(
  "analysis/network/figures/TCM_ExtinctionScenario_CumulativePrecip_RoundScale.pdf",
  width = 6, height = 4
)

# ggplot(tcm_climate,
#        aes(x = RoundPrecip, y = Robustness)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm", se = FALSE) +
#   facet_wrap(~ scenario) +
#   theme_minimal()

#### ---- Temp ----
# ggplot(tcm_climate,
#        aes(x = SpringTmean, y = Robustness)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm", se = FALSE) +
#   facet_wrap(~ scenario) +
#   theme_minimal()

# p2 <- ggplot(tcm_climate,
#        aes(x = CumulativeTmean, y = Robustness)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm", se = FALSE) +
#   facet_wrap(~ scenario) +
#   theme_minimal()
# 
# ggsave(
#   "analysis/network/figures/TCM_ExtinctionScenario_CumulativeTmean_RoundScale.pdf",
#   width = 6, height = 4
# )

# ggplot(tcm_climate,
#        aes(x = RoundTmean, y = Robustness)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm", se = FALSE) +
#   facet_wrap(~ scenario) +
#   theme_minimal()

#### ---- TempAnom ----
# ggplot(tcm_climate,
#        aes(x = SpringTmeanAnom, y = Robustness)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm", se = FALSE) +
#   facet_wrap(~ scenario) +
#   theme_minimal()

p2 <- ggplot(tcm_climate,
       aes(x = CumulativeTmeanAnom, y = Robustness)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ scenario) +
  theme_minimal()

ggsave(
  "analysis/network/figures/TCM_ExtinctionScenario_CumulativeTmeanAnom_RoundScale.pdf",
  width = 6, height = 4
)

# ggplot(tcm_climate,
#        aes(x = RoundTmeanAnom, y = Robustness)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm", se = FALSE) +
#   facet_wrap(~ scenario) +
#   theme_minimal()

#### ---- Combining Extinction Scenarios by Climate ----
(p1 | p2)

ggsave(
  "analysis/network/figures/TCM_ExtinctionScenario_Climate_RoundScale.pdf",
  width = 12, height = 4
)


## ---- SCM ----
### ---- Precipitation ----
# p1 <- ggplot(scm_climate,
#              aes(x = SpringPrecip, y = SCM_Robustness)) +
#   geom_point(aes(color = Site), alpha = 0.6) +
#   geom_smooth(method = "lm", se = FALSE) +
#   theme_minimal() +
#   labs(title = "Spring Precipitation",
#        x = "Spring precipitation",
#        y = "SCM Robustness")

p1 <- ggplot(scm_climate,
             aes(x = CumulativePrecip, y = SCM_Robustness)) +
  geom_point(aes(color = Site), alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Cumulative Precipitation")

# p3 <- ggplot(scm_climate,
#              aes(x = RoundPrecip, y = SCM_Robustness)) +
#   geom_point(aes(color = Site), alpha = 0.6) +
#   geom_smooth(method = "lm", se = FALSE) +
#   theme_minimal() +
#   labs(title = "Round Precipitation")
# 
# (p1 | p2 | p3)
# 
# ggsave(
#   "analysis/network/figures/SCM_Robustness_by_Precip_RoundScale.pdf",
#   width = 12, height = 4
# )

ggsave(
  "analysis/network/figures/SCM_Robustness_by_CumulativePrecip_RoundScale.pdf",
  width = 6, height = 4
)

### ---- Mean Temperature (absolute) ----
# p1 <- ggplot(scm_climate,
#              aes(x = SpringTmean, y = SCM_Robustness)) +
#   geom_point(aes(color = Site), alpha = 0.6) +
#   geom_smooth(method = "lm", se = FALSE) +
#   theme_minimal() +
#   labs(title = "Spring Mean Temperature")

# p2 <- ggplot(scm_climate,
#              aes(x = CumulativeTmean, y = SCM_Robustness)) +
#   geom_point(aes(color = Site), alpha = 0.6) +
#   geom_smooth(method = "lm", se = FALSE) +
#   theme_minimal() +
#   labs(title = "Cumulative Mean Temperature")

# p3 <- ggplot(scm_climate,
#              aes(x = RoundTmean, y = SCM_Robustness)) +
#   geom_point(aes(color = Site), alpha = 0.6) +
#   geom_smooth(method = "lm", se = FALSE) +
#   theme_minimal() +
#   labs(title = "Round Mean Temperature")
# 
# (p1 | p2 | p3)
# 
# ggsave(
#   "analysis/network/figures/SCM_Robustness_by_Temperature_RoundScale.pdf",
#   width = 12, height = 4
# )

# ggsave(
#   "analysis/network/figures/SCM_Robustness_by_CumulativeTmean_RoundScale.pdf",
#   width = 6, height = 4
# )

### ---- Mean Temperature (Anomaly) ----
# p1 <- ggplot(scm_climate,
#              aes(x = SpringTmeanAnom, y = SCM_Robustness)) +
#   geom_point(aes(color = Site), alpha = 0.6) +
#   geom_smooth(method = "lm", se = FALSE) +
#   theme_minimal()

p2 <- ggplot(scm_climate,
             aes(x = CumulativeTmeanAnom, y = SCM_Robustness)) +
  geom_point(aes(color = Site), alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Cumulative Mean Temperature Anomaly")

# p3 <- ggplot(scm_climate,
#              aes(x = RoundTmeanAnom, y = SCM_Robustness)) +
#   geom_point(aes(color = Site), alpha = 0.6) +
#   geom_smooth(method = "lm", se = FALSE) +
#   theme_minimal()
# 
# (p1 | p2 | p3)
# 
# ggsave(
#   "analysis/network/figures/SCM_Robustness_by_TemperatureAnomaly_RoundScale.pdf",
#   width = 12, height = 4
# )

ggsave(
  "analysis/network/figures/SCM_Robustness_by_CumulativeTmeanAnom_RoundScale.pdf",
  width = 6, height = 4
)

### ---- Combining Robustness by Climate ----
(p1 | p2)

ggsave(
  "analysis/network/figures/SCM_Robustness_by_CumulativeClimate_RoundScale.pdf",
  width = 12, height = 4
)

### ---- By extinction scenario to show whether extinction order changes climate sensitivity ----
#### ---- Precip ----
# ggplot(scm_climate,
#        aes(x = SpringPrecip, y = SCM_Robustness)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm", se = FALSE) +
#   facet_wrap(~ scenario) +
#   theme_minimal()

p1 <- ggplot(scm_climate,
       aes(x = CumulativePrecip, y = SCM_Robustness)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ scenario) +
  theme_minimal()

# ggplot(scm_climate,
#        aes(x = RoundPrecip, y = SCM_Robustness)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm", se = FALSE) +
#   facet_wrap(~ scenario) +
#   theme_minimal()

#### ---- Temp ----
# ggplot(scm_climate,
#        aes(x = SpringTmean, y = SCM_Robustness)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm", se = FALSE) +
#   facet_wrap(~ scenario) +
#   theme_minimal()

# p2 <- ggplot(scm_climate,
#        aes(x = CumulativeTmean, y = SCM_Robustness)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm", se = FALSE) +
#   facet_wrap(~ scenario) +
#   theme_minimal()

# ggplot(scm_climate,
#        aes(x = RoundTmean, y = SCM_Robustness)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm", se = FALSE) +
#   facet_wrap(~ scenario) +
#   theme_minimal()

#### ---- TempAnom ----
# ggplot(scm_climate,
#        aes(x = SpringTmeanAnom, y = SCM_Robustness)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm", se = FALSE) +
#   facet_wrap(~ scenario) +
#   theme_minimal()

p2 <- ggplot(scm_climate,
       aes(x = CumulativeTmeanAnom, y = SCM_Robustness)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ scenario) +
  theme_minimal()

# ggplot(scm_climate,
#        aes(x = RoundTmeanAnom, y = SCM_Robustness)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm", se = FALSE) +
#   facet_wrap(~ scenario) +
#   theme_minimal()

#### ---- Combining Extinction Scenarios by Climate ----
(p1 | p2)

ggsave(
  "analysis/network/figures/SCM_ExtinctionScenario_Climate_RoundScale.pdf",
  width = 12, height = 4
)

# ---- Harmonize TCM and SCM ----
## ---- Combine into one dataframe ----
comparison_df <- bind_rows(
  tcm_climate %>%
    transmute(
      Model = "TCM",
      Site,
      Year,
      SampleRound,
      scenario,
      # SpringPrecip,
      CumulativePrecip,
      # RoundPrecip,
      # SpringTmean,
      # CumulativeTmean,
      # RoundTmean,
      # SpringTmeanAnom,
      CumulativeTmeanAnom,
      # RoundTmeanAnom,
      Robustness = Robustness
    ),
  
  scm_climate %>%
    transmute(
      Model = "SCM",
      Site,
      Year,
      SampleRound,
      scenario,
      # SpringPrecip,
      CumulativePrecip,
      # RoundPrecip,
      # SpringTmean,
      # CumulativeTmean,
      # RoundTmean,
      # SpringTmeanAnom,
      CumulativeTmeanAnom,
      # RoundTmeanAnom,
      Robustness = SCM_Robustness
    )
)

## ---- Core comparison plots ----
### ---- Cumulative precipitation ----
p_compare_precip <- ggplot(
  comparison_df,
  aes(x = CumulativePrecip, y = Robustness)
) +
  geom_point(aes(color = Site), alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(~ Model, nrow = 1) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Climate sensitivity of robustness depends on extinction mechanism",
    subtitle = "Comparison of topological (TCM) vs stochastic (SCM) plant removals",
    x = "Cumulative precipitation",
    y = "Robustness"
  )

ggsave(
  "analysis/network/figures/TCM_vs_SCM_CumulativePrecip.pdf",
  p_compare_precip,
  width = 10,
  height = 4
)

### ---- Cumulative TmeanAnom ----
p_compare_tmeananom <- ggplot(
  comparison_df,
  aes(x = CumulativeTmeanAnom, y = Robustness)
) +
  geom_point(aes(color = Site), alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(~ Model, nrow = 1) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Climate sensitivity of robustness depends on extinction mechanism",
    subtitle = "Comparison of topological (TCM) vs stochastic (SCM) plant removals",
    x = "Cumulative temperature anomaly",
    y = "Robustness"
  )

ggsave(
  "analysis/network/figures/TCM_vs_SCM_CumulativeTmeanAnom.pdf",
  p_compare_tmeananom,
  width = 10,
  height = 4
)


### ---- Multi-climate panel ----
#### ---- Precip ----
# comparison_long <- comparison_df %>%
#   pivot_longer(
#     cols = c(SpringPrecip, CumulativePrecip, RoundPrecip),
#     names_to = "ClimateWindow",
#     values_to = "ClimateValue"
#   )
# 
# p_compare_multi <- ggplot(
#   comparison_long,
#   aes(x = ClimateValue, y = Robustness)
# ) +
#   geom_point(alpha = 0.35) +
#   geom_smooth(method = "lm", se = FALSE, color = "black") +
#   facet_grid(Model ~ ClimateWindow, scales = "free_x") +
#   theme_minimal(base_size = 12) +
#   labs(
#     x = "Climate variable",
#     y = "Robustness"
#   )
# 
# ggsave(
#   "analysis/network/figures/TCM_vs_SCM_Precip_AllWindows.pdf",
#   p_compare_multi,
#   width = 12,
#   height = 6
# )
# 
# #### ---- Temp ----
# comparison_long <- comparison_df %>%
#   pivot_longer(
#     cols = c(SpringTmean, CumulativeTmean, RoundTmean),
#     names_to = "ClimateWindow",
#     values_to = "ClimateValue"
#   )
# 
# p_compare_multi <- ggplot(
#   comparison_long,
#   aes(x = ClimateValue, y = Robustness)
# ) +
#   geom_point(alpha = 0.35) +
#   geom_smooth(method = "lm", se = FALSE, color = "black") +
#   facet_grid(Model ~ ClimateWindow, scales = "free_x") +
#   theme_minimal(base_size = 12) +
#   labs(
#     x = "Climate variable",
#     y = "Robustness"
#   )
# 
# ggsave(
#   "analysis/network/figures/TCM_vs_SCM_Temp_AllWindows.pdf",
#   p_compare_multi,
#   width = 12,
#   height = 6
# )
# 
# #### ---- TempAnom ----
# comparison_long <- comparison_df %>%
#   pivot_longer(
#     cols = c(SpringTmeanAnom, CumulativeTmeanAnom, RoundTmeanAnom),
#     names_to = "ClimateWindow",
#     values_to = "ClimateValue"
#   )
# 
# p_compare_multi <- ggplot(
#   comparison_long,
#   aes(x = ClimateValue, y = Robustness)
# ) +
#   geom_point(alpha = 0.35) +
#   geom_smooth(method = "lm", se = FALSE, color = "black") +
#   facet_grid(Model ~ ClimateWindow, scales = "free_x") +
#   theme_minimal(base_size = 12) +
#   labs(
#     x = "Climate variable",
#     y = "Robustness"
#   )
# 
# ggsave(
#   "analysis/network/figures/TCM_vs_SCM_TempAnom_AllWindows.pdf",
#   p_compare_multi,
#   width = 12,
#   height = 6
# )

## ---- Model Comparison ----
### ---- Precip ----
mod_compare_precip <- lmer(
  Robustness ~ CumulativePrecip * Model +
    (1 | Site) +
    (1 | Year),
  data = comparison_df
)
summary(mod_compare_precip)


#### ---- Effect Sizes ----
# newdat <- expand.grid(
#   SpringPrecip = seq(
#     min(comparison_df$SpringPrecip, na.rm = TRUE),
#     max(comparison_df$SpringPrecip, na.rm = TRUE),
#     length.out = 100
#   ),
#   Model = c("TCM", "SCM")
# )
# 
# newdat$pred <- predict(
#   mod_compare_precip,
#   newdata = newdat,
#   re.form = NA
# )
# 
# p_effect <- ggplot(
#   newdat,
#   aes(x = CumulativePrecip, y = pred, color = Model)
# ) +
#   geom_line(linewidth = 1.2) +
#   theme_minimal(base_size = 14) +
#   scale_color_manual(
#     values = c("TCM" = "#1b9e77", "SCM" = "#d95f02")
#   ) +
#   labs(
#     title = "Climate sensitivity of network robustness depends on extinction mechanism",
#     subtitle = "Predicted effects of cumulative precipitation from mixed-effects model",
#     x = "Cumulative precipitation",
#     y = "Predicted robustness",
#     color = "Model"
#   )
# 
# p_effect


# ## ---- Model x Climate interaction
# mod_compare_precip <- lmer(
#   Robustness ~ SpringPrecip * Model +
#     # scenario +
#     (1 | Site) +
#     (1 | Year),
#   data = comparison_df
# )
# 
# summary(mod_compare_precip)
# 
# ## ---- Effect-size
# slopes_precip <- tidy(mod_compare_precip, effects = "fixed") %>%
#   filter(str_detect(term, "SpringPrecip"))
# 
# slopes_precip

# ---- Community Resistance ----

## ---- Join climate to community metrics ----
load('analysis/network/saved/corMets_PlantPollinator_YearSR.Rdata')

# cor.dats: Site x Year x SampleRound network metrics

cor.dats <- cor.dats %>%
  separate(
    Site,
    into = c("Site", "Year", "SampleRound"),
    sep = "\\.",
    convert = TRUE,   # converts Year and SampleRound to numeric
    remove = FALSE    # keeps original Site column if you want it
  )


cor.dats <- cor.dats %>%
  left_join(
    climate %>%
      select(
        Site, Year, SampleRound,
        CumulativePrecip,
        CumulativeTmean,
        CumulativeTmeanAnom
      ),
    by = c("Site", "Year", "SampleRound")
  )


cor.dats <- cor.dats %>%
  mutate(
    Site = factor(Site),
    Year = factor(Year),
    SampleRound = factor(SampleRound)
  )

## ---- Define metrics and labels ----
# metrics <- c(
#   "FunRedundancy.Pols",
#   "FunRedundancy.Plants",
#   "functional.complementarity.HL",
#   "functional.complementarity.LL",
#   "mean.number.of.links.HL",
#   "mean.number.of.links.LL"
# )

metrics <- c(
  # Pollinators
  "functional.complementarity.HL",
  "FunRedundancy.Pols",
  "mean.number.of.links.HL",
  "number.of.species.HL",
  
  # Plants
  "functional.complementarity.LL",
  "FunRedundancy.Plants",
  "mean.number.of.links.LL",
  "number.of.species.LL"
)


metric_labels <- c(
  "functional.complementarity.HL" = "Pollinator Complementarity",
  "FunRedundancy.Pols"             = "Pollinator Redundancy",
  "mean.number.of.links.HL"        = "Pollinator Generalization",
  "number.of.species.HL"           = "Pollinator Richness",
  
  "functional.complementarity.LL" = "Plant Complementarity",
  "FunRedundancy.Plants"           = "Plant Redundancy",
  "mean.number.of.links.LL"        = "Plant Generalization",
  "number.of.species.LL"           = "Plant Richness"
)


## ---- Model ----
fit_climate_models <- function(response_vars, climate_var, data) {
  
  mods <- lapply(response_vars, function(y) {
    formula <- as.formula(
      paste0(y, " ~ scale(", climate_var, ") + (1 | Site) + (1 | Year)")
    )
    lmer(formula, data = data, REML = FALSE)
  })
  
  names(mods) <- response_vars
  return(mods)
}

mods_precip <- fit_climate_models(metrics, "CumulativePrecip", cor.dats)
# mods_tmean  <- fit_climate_models(metrics, "CumulativeTmean", cor.dats)
mods_anom   <- fit_climate_models(metrics, "CumulativeTmeanAnom", cor.dats)

### ---- Extract Results ----
extract_climate_results <- function(mods, climate_var) {
  lapply(names(mods), function(name) {
    broom.mixed::tidy(mods[[name]], effects = "fixed") %>%
      filter(term == paste0("scale(", climate_var, ")")) %>%
      mutate(
        Metric = name,
        ClimateVar = climate_var
      )
  }) %>%
    bind_rows()
}

results_all <- bind_rows(
  extract_climate_results(mods_precip, "CumulativePrecip"),
  # extract_climate_results(mods_tmean,  "CumulativeTmean"),
  extract_climate_results(mods_anom,   "CumulativeTmeanAnom")
)

results_all

sig_table <- results_all %>%
  mutate(
    signif = if_else(p.value < 0.05, "significant", "ns")
  ) %>%
  select(Metric, ClimateVar, signif)

sig_table

## ---- Predictions ----
predict_metric <- function(mod, climate_var, values) {
  new_data <- data.frame(x = values)
  colnames(new_data) <- climate_var
  
  mm <- model.matrix(
    as.formula(paste0("~ scale(", climate_var, ")")),
    new_data
  )
  
  fit <- mm %*% fixef(mod)
  se  <- sqrt(diag(mm %*% vcov(mod) %*% t(mm)))
  
  new_data %>%
    mutate(
      fit = as.vector(fit),
      lower = fit - 1.96 * se,
      upper = fit + 1.96 * se
    )
}

## ---- Plot ----
plot_climate_effects <- function(mods, climate_var, data) {
  
  climate_range <- seq(
    min(data[[climate_var]], na.rm = TRUE),
    max(data[[climate_var]], na.rm = TRUE),
    length.out = 100
  )
  
  preds <- lapply(names(mods), function(m) {
    predict_metric(mods[[m]], climate_var, climate_range) %>%
      mutate(metric = m)
  }) %>%
    bind_rows() %>%
    left_join(
      sig_table %>%
        filter(ClimateVar == climate_var),
      by = c("metric" = "Metric")
    )
  
  raw <- data %>%
    pivot_longer(
      cols = all_of(metrics),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(
        metric,
        levels = c(
          # Pollinators
          "functional.complementarity.HL",
          "FunRedundancy.Pols",
          "mean.number.of.links.HL",
          "number.of.species.HL",
          
          # Plants
          "functional.complementarity.LL",
          "FunRedundancy.Plants",
          "mean.number.of.links.LL",
          "number.of.species.LL"
        )
      )
    )
  
  preds <- preds %>%
    mutate(metric = factor(metric, levels = levels(raw$metric)))
  
  ggplot(preds, aes(x = .data[[climate_var]], y = fit)) +
    geom_line(aes(linetype = signif), color = "black") +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = "gray80", alpha = 0.4
    ) +
    geom_point(
      data = raw,
      aes(
        x = .data[[climate_var]],
        y = value,
        color = Site,
        shape = SampleRound
      ),
      inherit.aes = FALSE,
      alpha = 0.6
    ) +
    scale_linetype_manual(
      values = c(significant = "solid", ns = "dashed"),
      guide = "none"
    ) +
    facet_wrap(
      ~ metric,
      nrow = 2,
      scales = "free_y",
      labeller = labeller(metric = metric_labels)
    ) +
    theme_minimal(base_size = 14) +
    labs(
      x = climate_var,
      y = "Predicted community resistance metric",
      color = "Site"
    ) +
    guides(shape = "none")
}

p_precip <- plot_climate_effects(mods_precip, "CumulativePrecip", cor.dats)

ggsave("analysis/network/figures/CommunityResistance_CumulativePrecip.pdf",
       p_precip, width = 11, height = 8)

# p_tmean <- plot_climate_effects(mods_tmean, "CumulativeTmean", cor.dats)
# 
# ggsave("analysis/network/figures/CommunityResistance_CumulativeTmean.pdf",
#        p_tmean, width = 11, height = 8)

p_anom <- plot_climate_effects(mods_anom, "CumulativeTmeanAnom", cor.dats)

ggsave("analysis/network/figures/CommunityResistance_CumulativeTmeanAnom.pdf",
       p_anom, width = 11, height = 8)
