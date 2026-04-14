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

# ---- Load Data ----
setwd(file.path(local.path, "extinction-cascades"))
load("analysis/network/robustness_metrics/YearSR_robustness_climate.Rdata")

setwd(file.path(local.path, "skyIslands_saved"))
climate <- read.csv("data/relational/original/climate.csv")

setwd(file.path(local.path, "extinction-cascades"))
load("analysis/network/saved/corMets_PlantPollinator_YearSR.Rdata")

## *******************************************************************
# ---- TCM Robustness Models ----
## *******************************************************************

## ---- Merge climate into robustness results ----
tcm_climate <- robustness_results %>%
  left_join(climate, by = c("Site", "Year", "SampleRound"))

## ---- Full Climate Model ----
mod_tcm <- lmer(
  Robustness ~
    scale(APi) +
    scale(WindowTmeanAnom) +
    scenario +
    (1 | Site) +
    (1 | Year),
  data = tcm_climate
)

summary(mod_tcm)
vif(mod_tcm)

## ---- Model-based predictions ----
tcm_preds_api <- ggpredict(
  mod_tcm,
  terms = c("APi [all]", "scenario")
)

tcm_preds_anom <- ggpredict(
  mod_tcm,
  terms = c("WindowTmeanAnom [all]", "scenario")
)

## ---- Plots ----
p_tcm_api <- ggplot(tcm_preds_api,
                    aes(x = x, y = predicted,
                        color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(
    x = "Antecedent Precipitation Index",
    y = "TCM Robustness",
    color = "Extinction scenario",
    fill = "Extinction scenario"
  )

p_tcm_anom <- ggplot(tcm_preds_anom,
                     aes(x = x, y = predicted,
                         color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(
    x = "Seasonal Temperature Anomaly",
    y = NULL,
    color = "Extinction scenario",
    fill = "Extinction scenario"
  )

(p_tcm_api | p_tcm_anom) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(
  "analysis/network/figures/TCM_Robustness_APi_WindowTmeanAnom.pdf",
  width = 12, height = 4
)

## *******************************************************************
# ---- SCM Robustness Models ----
## *******************************************************************

scm_climate <- scm_results %>%
  left_join(climate, by = c("Site", "Year", "SampleRound"))

mod_scm <- lmer(
  SCM_Robustness ~
    scale(APi) * scenario +
    scale(WindowTmeanAnom) * scenario +
    (1 | Site) +
    (1 | Year),
  data = scm_climate
)

summary(mod_scm)

## ---- Predictions ----
scm_preds_api <- ggpredict(
  mod_scm,
  terms = c("APi [all]", "scenario")
)

scm_preds_anom <- ggpredict(
  mod_scm,
  terms = c("WindowTmeanAnom [all]", "scenario")
)

## ---- Plots ----
p_scm_api <- ggplot(scm_preds_api,
                    aes(x = x, y = predicted,
                        color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(
    x = "Antecedent Precipitation Index",
    y = "SCM Robustness",
    color = "SCM scenario",
    fill = "SCM scenario"
  )

p_scm_anom <- ggplot(scm_preds_anom,
                     aes(x = x, y = predicted,
                         color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(
    x = "Seasonal Temperature Anomaly",
    y = NULL
  )

(p_scm_api | p_scm_anom) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(
  "analysis/network/figures/SCM_Robustness_APi_WindowTmeanAnom.pdf",
  width = 12, height = 4
)

## *******************************************************************
# ---- Community Resistance Models ----
## *******************************************************************

cor.dats <- cor.dats %>%
  separate(
    Site,
    into = c("Site", "Year", "SampleRound"),
    sep = "\\.",
    convert = TRUE,
    remove = FALSE
  ) %>%
  left_join(
    climate %>%
      select(Site, Year, SampleRound, APi, WindowTmeanAnom),
    by = c("Site", "Year", "SampleRound")
  ) %>%
  mutate(
    Site = factor(Site),
    Year = factor(Year),
    SampleRound = factor(SampleRound)
  )

## ---- Define metrics ----
metrics <- c(
  "functional.complementarity.HL",
  "FunRedundancy.Pols",
  "mean.number.of.links.HL",
  "number.of.species.HL",
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
  "functional.complementarity.LL"  = "Plant Complementarity",
  "FunRedundancy.Plants"           = "Plant Redundancy",
  "mean.number.of.links.LL"        = "Plant Generalization",
  "number.of.species.LL"           = "Plant Richness"
)

## ---- Fit models ----
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

mods_api  <- fit_climate_models(metrics, "APi", cor.dats)
mods_anom <- fit_climate_models(metrics, "WindowTmeanAnom", cor.dats)

## ---- Extract results ----
extract_climate_results <- function(mods, climate_var) {
  lapply(names(mods), function(name) {
    broom.mixed::tidy(mods[[name]], effects = "fixed") %>%
      filter(term == paste0("scale(", climate_var, ")")) %>%
      mutate(Metric = name, ClimateVar = climate_var)
  }) %>%
    bind_rows()
}

results_all <- bind_rows(
  extract_climate_results(mods_api,  "APi"),
  extract_climate_results(mods_anom, "WindowTmeanAnom")
)

sig_table <- results_all %>%
  mutate(signif = if_else(p.value < 0.05, "significant", "ns")) %>%
  select(Metric, ClimateVar, signif)

## ---- Plot resistance effects ----
# (plot_climate_effects function stays identical to your original,
#  it uses climate_var as a string so it works with any variable name)
source("plot_climate_effects.R")

### ---- Predictions ----
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

p_api  <- plot_climate_effects(mods_api,  "APi",            cor.dats)
p_anom <- plot_climate_effects(mods_anom, "WindowTmeanAnom", cor.dats)

ggsave("analysis/network/figures/CommunityResistance_APi.pdf",
       p_api,  width = 11, height = 8)

ggsave("analysis/network/figures/CommunityResistance_WindowTmeanAnom.pdf",
       p_anom, width = 11, height = 8)