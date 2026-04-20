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
# ---- TCM Robustness Models (one model per scenario) ----
## *******************************************************************

tcm_climate <- robustness_results %>%
  left_join(climate, by = c("Site", "Year", "SampleRound"))

# TCM scenarios
tcm_scenarios <- c("abundance_low", "abundance_high", "degree_low", "degree_high")

# Fit one model per TCM scenario
tcm_models <- lapply(tcm_scenarios, function(s) {
  dat <- filter(tcm_climate, scenario == s)
  lmer(
    Robustness ~
      scale(APi) +
      scale(WindowTmeanAnom) +
      # (1 | Year) +
      (1 | Site),
    data = dat
  )
})

names(tcm_models) <- tcm_scenarios

# Summaries and diagnostics
lapply(tcm_scenarios, function(s) {
  cat("\n\n============================\n")
  cat("TCM scenario:", s, "\n")
  cat("============================\n")
  print(summary(tcm_models[[s]]))
  cat("\nVIF:\n")
  print(vif(tcm_models[[s]]))
  cat("\nR2:\n")
  print(r2(tcm_models[[s]]))
})

## ---- Extract R2 and p-values per scenario x climate variable ----
tcm_stats <- lapply(tcm_scenarios, function(s) {
  mod     <- tcm_models[[s]]
  r2_vals <- r2(mod)
  coefs   <- broom.mixed::tidy(mod, effects = "fixed")

  p_api  <- coefs %>% filter(term == "scale(APi)")             %>% pull(p.value)
  p_anom <- coefs %>% filter(term == "scale(WindowTmeanAnom)") %>% pull(p.value)

  data.frame(
    scenario       = s,
    Marginal_R2    = round(r2_vals$R2_marginal,    2),
    Conditional_R2 = round(r2_vals$R2_conditional, 2),
    p_APi          = p_api,
    p_Anom         = p_anom
  )
}) %>%
  bind_rows() %>%
  mutate(
    # Legend label: scenario name + both R² values
    legend_label = paste0(
      scenario,
      "\n  R²m = ", Marginal_R2,
      ", R²c = ", Conditional_R2
    )
  )

## ---- Predictions per scenario ----
predict_tcm <- function(s, climate_var) {
  ggpredict(tcm_models[[s]], terms = paste0(climate_var, " [all]")) %>%
    as.data.frame() %>%
    mutate(scenario = s, climate_var = climate_var)
}

tcm_preds <- bind_rows(
  lapply(tcm_scenarios, predict_tcm, climate_var = "APi"),
  lapply(tcm_scenarios, predict_tcm, climate_var = "WindowTmeanAnom")
) %>%
  left_join(
    tcm_stats %>% select(scenario, legend_label, p_APi, p_Anom),
    by = "scenario"
  ) %>%
  mutate(
    # Significance flag per row based on which climate variable this prediction is for
    sig = case_when(
      climate_var == "APi"              & p_APi  < 0.05 ~ "significant",
      climate_var == "WindowTmeanAnom"  & p_Anom < 0.05 ~ "significant",
      TRUE ~ "ns"
    ),
    climate_var = factor(
      climate_var,
      levels = c("APi", "WindowTmeanAnom"),
      labels = c("Antecedent Precipitation Index", "Seasonal Temperature Anomaly")
    ),
    # Order legend by scenario; label carries R² info
    legend_label = factor(legend_label, levels = tcm_stats$legend_label)
  )

## ---- Plot: facet_wrap(~ climate_var), lines colored by scenario ----
p_tcm <- ggplot(tcm_preds, aes(x = x, y = predicted,
                                color = legend_label,
                                fill  = legend_label)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.15, color = NA) +
  geom_line(aes(linetype = sig), linewidth = 1) +
  scale_linetype_manual(
    values = c(significant = "solid", ns = "dashed"),
    guide  = "none"
  ) +
  facet_wrap(~ climate_var, scales = "free_x") +
  theme_minimal(base_size = 13) +
  labs(
    x      = NULL,
    y      = "TCM Robustness",
    color  = "Extinction ordering",
    fill   = "Extinction ordering"
  ) +
  theme(
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.text      = element_text(size = 9),
    legend.key.height = unit(1.2, "lines")
  )

ggsave(
  "analysis/network/figures/TCM_Robustness_by_climate.pdf",
  p_tcm,
  width = 11, height = 5
)

## *******************************************************************
# ---- SCM Robustness Models (one model per scenario) ----
## *******************************************************************

scm_climate <- scm_results %>%
  left_join(climate, by = c("Site", "Year", "SampleRound"))

# SCM scenarios
scm_scenarios <- c("scm_random_plant", "scm_dominant_plant")

# Fit one model per SCM scenario
scm_models <- lapply(scm_scenarios, function(s) {
  dat <- filter(scm_climate, scenario == s)
  lmer(
    SCM_Robustness ~
      scale(APi) +
      scale(WindowTmeanAnom) +
      (1 | Year) +
      (1 | Site),
    data = dat
  )
})

names(scm_models) <- scm_scenarios

# Summaries and diagnostics
lapply(scm_scenarios, function(s) {
  cat("\n\n============================\n")
  cat("SCM scenario:", s, "\n")
  cat("============================\n")
  print(summary(scm_models[[s]]))
  cat("\nVIF:\n")
  print(vif(scm_models[[s]]))
  cat("\nR2:\n")
  print(r2(scm_models[[s]]))
})

## ---- Extract R2 and p-values per scenario x climate variable ----
scm_stats <- lapply(scm_scenarios, function(s) {
  mod     <- scm_models[[s]]
  r2_vals <- r2(mod)
  coefs   <- broom.mixed::tidy(mod, effects = "fixed")

  p_api  <- coefs %>% filter(term == "scale(APi)")             %>% pull(p.value)
  p_anom <- coefs %>% filter(term == "scale(WindowTmeanAnom)") %>% pull(p.value)

  data.frame(
    scenario       = s,
    Marginal_R2    = round(r2_vals$R2_marginal,    2),
    Conditional_R2 = round(r2_vals$R2_conditional, 2),
    p_APi          = p_api,
    p_Anom         = p_anom
  )
}) %>%
  bind_rows() %>%
  mutate(
    legend_label = paste0(
      scenario,
      "\n  R²m = ", Marginal_R2,
      ", R²c = ", Conditional_R2
    )
  )

## ---- Predictions per scenario ----
predict_scm <- function(s, climate_var) {
  ggpredict(scm_models[[s]], terms = paste0(climate_var, " [all]")) %>%
    as.data.frame() %>%
    mutate(scenario = s, climate_var = climate_var)
}

scm_preds <- bind_rows(
  lapply(scm_scenarios, predict_scm, climate_var = "APi"),
  lapply(scm_scenarios, predict_scm, climate_var = "WindowTmeanAnom")
) %>%
  left_join(
    scm_stats %>% select(scenario, legend_label, p_APi, p_Anom),
    by = "scenario"
  ) %>%
  mutate(
    sig = case_when(
      climate_var == "APi"             & p_APi  < 0.05 ~ "significant",
      climate_var == "WindowTmeanAnom" & p_Anom < 0.05 ~ "significant",
      TRUE ~ "ns"
    ),
    climate_var = factor(
      climate_var,
      levels = c("APi", "WindowTmeanAnom"),
      labels = c("Antecedent Precipitation Index", "Seasonal Temperature Anomaly")
    ),
    legend_label = factor(legend_label, levels = scm_stats$legend_label)
  )

## ---- Plot: facet_wrap(~ climate_var), lines colored by scenario ----
p_scm <- ggplot(scm_preds, aes(x = x, y = predicted,
                                color = legend_label,
                                fill  = legend_label)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.15, color = NA) +
  geom_line(aes(linetype = sig), linewidth = 1) +
  scale_linetype_manual(
    values = c(significant = "solid", ns = "dashed"),
    guide  = "none"
  ) +
  facet_wrap(~ climate_var, scales = "free_x") +
  theme_minimal(base_size = 13) +
  labs(
    x     = NULL,
    y     = "SCM Robustness",
    color = "Extinction ordering",
    fill  = "Extinction ordering"
  ) +
  theme(
    strip.text        = element_text(face = "bold"),
    panel.grid.minor  = element_blank(),
    legend.text       = element_text(size = 9),
    legend.key.height = unit(1.2, "lines")
  )

ggsave(
  "analysis/network/figures/SCM_Robustness_by_climate.pdf",
  p_scm,
  width = 11, height = 5
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
