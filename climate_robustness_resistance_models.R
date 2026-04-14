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
library(performance)
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
tcm_climate <- robustness_results %>%
  left_join(climate, by = c("Site", "Year", "SampleRound"))

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
r2(mod_tcm)

## ---- Predictions ----
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
  geom_line(linewidth = 1, linetype = "solid") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(
    x = "Antecedent Precipitation Index",
    y = "TCM Robustness",
    color = "Extinction scenario",
    fill = "Extinction scenario"
  )

p_tcm_api <- p_tcm_api +
  annotate("text",
           x = Inf, y = Inf,
           label = "Marginal R² = 0.83\nConditional R² = 0.89",
           hjust = 1.1, vjust = 1.5,
           size = 3.5)

p_tcm_anom <- ggplot(tcm_preds_anom,
                     aes(x = x, y = predicted,
                         color = group, fill = group)) +
  geom_line(linewidth = 1, linetype = "dashed") +
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
vif(mod_scm)
r2(mod_scm)

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
  plot_layout(guides = "collect") +
  plot_annotation(
    caption = "Marginal R² = 0.26; Conditional R² = 0.78"
  ) &
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

### ---- Extract r2 ----
extract_r2 <- function(mods, climate_var) {
  lapply(names(mods), function(name) {
    r2_vals <- performance::r2(mods[[name]])
    
    data.frame(
      Metric = name,
      ClimateVar = climate_var,
      Marginal_R2 = r2_vals$R2_marginal,
      Conditional_R2 = r2_vals$R2_conditional
    )
  }) %>%
    bind_rows()
}

r2_api  <- extract_r2(mods_api, "APi")
r2_anom <- extract_r2(mods_anom, "WindowTmeanAnom")

r2_all <- bind_rows(r2_api, r2_anom)

results_all <- results_all %>%
  left_join(r2_all, by = c("Metric", "ClimateVar"))

sig_table <- results_all %>%
  mutate(signif = if_else(p.value < 0.05, "significant", "ns")) %>%
  select(Metric, ClimateVar, signif)

## ---- Plot resistance effects ----
r2_labels <- r2_all %>%
  mutate(
    r2_text = paste0("R²m = ", round(Marginal_R2, 2),
                     "\nR²c = ", round(Conditional_R2, 2))
  ) %>%
  select(Metric, ClimateVar, r2_text)


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
    mutate(metric = factor(metric, levels = levels(raw$metric))) %>%
    left_join(
      r2_labels %>% filter(ClimateVar == climate_var),
      by = c("metric" = "Metric")
    )
  
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
    geom_text(
      data = preds,
      aes(x = -Inf, y = Inf, label = r2_text),
      hjust = -0.05, vjust = 1.1,
      size = 3,
      inherit.aes = FALSE
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
