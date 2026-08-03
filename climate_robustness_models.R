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
geography <- read.csv("data/relational/original/geography.csv")

setwd(file.path(local.path, "extinction-cascades"))
exposure <- read_csv("exposure.csv")
load("analysis/network/saved/corMets_PlantPollinator_YearSR.Rdata")

## *******************************************************************
# ---- TCM Robustness Models (one model per scenario) ----
## *******************************************************************
tcm_climate <- robustness_results %>%
  left_join(exposure, by = c("Site", "Year", "SampleRound"))

# TCM scenarios
tcm_scenarios <- c("abundance_low", "abundance_high", "degree_low", "degree_high")

## *******************************************************************
## ---- Model Set 1 Mechanistic climate variables APi + GDD ----
## *******************************************************************
# Fit one model per TCM scenario
tcm_api_models <- lapply(tcm_scenarios, function(s) {
  dat <- filter(tcm_climate, scenario == s)
  lmer(
    Robustness ~
      scale(APi) +
      scale(GDD) +
      (1 | Year) +
      (1 | Site),
    data = dat
  )
})

names(tcm_api_models) <- tcm_scenarios

tcm_api_results <-
  bind_rows(
    lapply(names(tcm_api_models), function(s){
      
      broom.mixed::tidy(tcm_api_models[[s]]) %>%
        mutate(
          scenario = s,
          model = "APi_GDD"
        )
    })
  )

# Summaries and diagnostics
lapply(tcm_scenarios, function(s) {
  cat("\n\n============================\n")
  cat("TCM scenario:", s, "\n")
  cat("============================\n")
  print(summary(tcm_api_models[[s]]))
  cat("\nVIF:\n")
  print(vif(tcm_api_models[[s]]))
  cat("\nR2:\n")
  print(r2(tcm_api_models[[s]]))
})

### ---- Extract R2 and p-values per scenario x climate variable ----
tcm_api_stats <- lapply(tcm_scenarios, function(s) {
  mod     <- tcm_api_models[[s]]
  r2_vals <- r2(mod)
  coefs   <- broom.mixed::tidy(mod, effects = "fixed")
  
  p_api  <- coefs %>% filter(term == "scale(APi)")             %>% pull(p.value)
  p_gdd <- coefs %>% filter(term == "scale(GDD)") %>% pull(p.value)
  
  data.frame(
    scenario       = s,
    Marginal_R2    = round(r2_vals$R2_marginal,    2),
    Conditional_R2 = round(r2_vals$R2_conditional, 2),
    p_APi          = p_api,
    p_GDD          = p_gdd
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

### ---- Predictions per scenario ----
predict_api_tcm <- function(s, climate_var) {
  ggpredict(tcm_api_models[[s]], terms = paste0(climate_var, " [all]")) %>%
    as.data.frame() %>%
    mutate(scenario = s, climate_var = climate_var)
}

tcm_api_preds <- bind_rows(
  lapply(tcm_scenarios, predict_api_tcm, climate_var = "APi"),
  lapply(tcm_scenarios, predict_api_tcm, climate_var = "GDD")
) %>%
  left_join(
    tcm_api_stats %>% select(scenario, legend_label, p_APi, p_GDD),
    by = "scenario"
  ) %>%
  mutate(
    # Significance flag per row based on which climate variable this prediction is for
    sig = case_when(
      climate_var == "APi"              & p_APi  < 0.05 ~ "significant",
      climate_var == "GDD"              & p_GDD  < 0.05 ~ "significant",
      TRUE ~ "ns"
    ),
    climate_var = factor(
      climate_var,
      levels = c("APi", "GDD"),
      labels = c("Antecedent Precipitation Index", "Growing Degree Days")
    ),
    # Order legend by scenario; label carries R² info
    legend_label = factor(legend_label, levels = tcm_api_stats$legend_label)
  )

### ---- Plot: facet_wrap(~ climate_var), lines colored by scenario ----
p_api_tcm <- ggplot(tcm_api_preds, aes(x = x, y = predicted,
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
  p_api_tcm,
  width = 11, height = 5
)

## *******************************************************************
## ---- Model Set 2 Multivariate climate exposure PC1 + PC2 ----
## *******************************************************************
# Fit models to PCA
tcm_pca_models <-
  lapply(tcm_scenarios, function(s){
    dat <- filter(tcm_climate,
                  scenario==s)
    lmer(
      Robustness ~
        scale(PC1)+
        scale(PC2)+
        (1|Year)+
        (1|Site),
      data=dat
    )
  })

names(tcm_pca_models) <- tcm_scenarios

tcm_pca_results <-
  bind_rows(
    lapply(names(tcm_pca_models), function(s){
      
      broom.mixed::tidy(tcm_pca_models[[s]]) %>%
        mutate(
          scenario = s,
          model = "PCA"
        )
      
    })
  )

# PCA Summaries and diagnostics
lapply(tcm_scenarios, function(s) {
  cat("\n\n============================\n")
  cat("TCM scenario:", s, "\n")
  cat("============================\n")
  print(summary(tcm_pca_models[[s]]))
  cat("\nVIF:\n")
  print(vif(tcm_pca_models[[s]]))
  cat("\nR2:\n")
  print(r2(tcm_pca_models[[s]]))
})

### ---- Extract R2 and p-values per scenario x climate variable ----
# PCA Extract R2 and p-values per scenario x climate variable
tcm_pca_stats <- lapply(tcm_scenarios, function(s) {
  mod     <- tcm_pca_models[[s]]
  r2_vals <- r2(mod)
  coefs <- tidy(mod,effects="fixed")
  
  p_PC1  <- coefs %>% filter(term=="scale(PC1)") %>% pull(p.value)
  p_PC2  <- coefs %>% filter(term=="scale(PC2)") %>% pull(p.value)
  
  data.frame(
    scenario       = s,
    Marginal_R2    = round(r2_vals$R2_marginal,    2),
    Conditional_R2 = round(r2_vals$R2_conditional, 2),
    p_PC1          = p_PC1,
    p_PC2          = p_PC2
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

### ---- Predictions per scenario ----
predict_pca_tcm <- function(s, climate_var) {
  ggpredict(tcm_pca_models[[s]], terms = paste0(climate_var, " [all]")) %>%
    as.data.frame() %>%
    mutate(scenario = s, climate_var = climate_var)
}

tcm_pca_preds <- bind_rows(
  lapply(tcm_scenarios, predict_pca_tcm, climate_var = "PC1"),
  lapply(tcm_scenarios, predict_pca_tcm, climate_var = "PC2")
) %>%
  left_join(
    tcm_pca_stats %>% select(scenario, legend_label, p_PC1, p_PC2),
    by = "scenario"
  ) %>%
  mutate(
    # Significance flag per row based on which climate variable this prediction is for
    sig = case_when(
      climate_var == "PC1"              & p_PC1  < 0.05 ~ "significant",
      climate_var == "PC2"              & p_PC2 < 0.05 ~ "significant",
      TRUE ~ "ns"
    ),
    climate_var = factor(
      climate_var,
      levels = c("PC1", "PC2"),
      labels = c("Principal Component 1", "Principal Component 2")
    ),
    # Order legend by scenario; label carries R² info
    legend_label = factor(legend_label, levels = tcm_pca_stats$legend_label)
  )

### ---- Plot: facet_wrap(~ climate_var), lines colored by scenario ----
p_tcm_pca <- ggplot(tcm_pca_preds, aes(x = x, y = predicted,
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
  "analysis/network/figures/TCM_Robustness_by_PCAclimate.pdf",
  p_tcm_pca,
  width = 11, height = 5
)

## ---- Compare models ----
tcm_compare_results <- bind_rows(
  tcm_api_results,
  tcm_pca_results
)

## *******************************************************************
# ---- SCM Robustness Models (one model per scenario) ----
## *******************************************************************
scm_climate <- scm_results %>%
  left_join(exposure, by = c("Site", "Year", "SampleRound"))

# SCM scenarios
scm_scenarios <- c("scm_random_plant", "scm_dominant_plant")

## *******************************************************************
## ---- Model Set 1 Mechanistic climate variables APi + GDD ----
## *******************************************************************
# Fit one model per SCM scenario
scm_api_models <- lapply(scm_scenarios, function(s) {
  dat <- filter(scm_climate, scenario == s)
  lmer(
    SCM_Robustness ~
      scale(APi) +
      scale(GDD) +
      (1 | Year) +
      (1 | Site),
    data = dat
  )
})

names(scm_api_models) <- scm_scenarios

scm_api_results <-
  bind_rows(
    lapply(names(scm_api_models), function(s){
      
      broom.mixed::tidy(scm_api_models[[s]]) %>%
        mutate(
          scenario = s,
          model = "APi_GDD"
        )
    })
  )

# Summaries and diagnostics
lapply(scm_scenarios, function(s) {
  cat("\n\n============================\n")
  cat("SCM scenario:", s, "\n")
  cat("============================\n")
  print(summary(scm_api_models[[s]]))
  cat("\nVIF:\n")
  print(vif(scm_api_models[[s]]))
  cat("\nR2:\n")
  print(r2(scm_api_models[[s]]))
})

### ---- Extract R2 and p-values per scenario x climate variable ----
scm_api_stats <- lapply(scm_scenarios, function(s) {
  mod     <- scm_api_models[[s]]
  r2_vals <- r2(mod)
  coefs   <- broom.mixed::tidy(mod, effects = "fixed")

  p_api  <- coefs %>% filter(term == "scale(APi)")             %>% pull(p.value)
  p_gdd  <- coefs %>% filter(term == "scale(GDD)")             %>% pull(p.value)

  data.frame(
    scenario       = s,
    Marginal_R2    = round(r2_vals$R2_marginal,    2),
    Conditional_R2 = round(r2_vals$R2_conditional, 2),
    p_APi          = p_api,
    p_GDD          = p_gdd
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

### ---- Predictions per scenario ----
predict_api_scm <- function(s, climate_var) {
  ggpredict(scm_api_models[[s]], terms = paste0(climate_var, " [all]")) %>%
    as.data.frame() %>%
    mutate(scenario = s, climate_var = climate_var)
}

scm_api_preds <- bind_rows(
  lapply(scm_scenarios, predict_api_scm, climate_var = "APi"),
  lapply(scm_scenarios, predict_api_scm, climate_var = "GDD")
) %>%
  left_join(
    scm_api_stats %>% select(scenario, legend_label, p_APi, p_GDD),
    by = "scenario"
  ) %>%
  mutate(
    sig = case_when(
      climate_var == "APi"             & p_APi  < 0.05 ~ "significant",
      climate_var == "GDD"             & p_GDD < 0.05 ~ "significant",
      TRUE ~ "ns"
    ),
    climate_var = factor(
      climate_var,
      levels = c("APi", "GDD"),
      labels = c("Antecedent Precipitation Index", "Growing Degree Days")
    ),
    legend_label = factor(legend_label, levels = scm_api_stats$legend_label)
  )

### ---- Plot: facet_wrap(~ climate_var), lines colored by scenario ----
p_api_scm <- ggplot(scm_api_preds, aes(x = x, y = predicted,
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
  p_api_scm,
  width = 11, height = 5
)

## *******************************************************************
## ---- Model Set 2 Multivariate climate exposure PC1 + PC2 ----
## *******************************************************************
# Fit one model per SCM scenario
scm_pca_models <- lapply(scm_scenarios, function(s) {
  dat <- filter(scm_climate, scenario == s)
  lmer(
    SCM_Robustness ~
      scale(PC1) +
      scale(PC2) +
      (1 | Year) +
      (1 | Site),
    data = dat
  )
})

names(scm_pca_models) <- scm_scenarios

scm_pca_results <-
  bind_rows(
    lapply(names(scm_pca_models), function(s){
      
      broom.mixed::tidy(scm_pca_models[[s]]) %>%
        mutate(
          scenario = s,
          model = "PCA"
        )
    })
  )

# Summaries and diagnostics
lapply(scm_scenarios, function(s) {
  cat("\n\n============================\n")
  cat("SCM scenario:", s, "\n")
  cat("============================\n")
  print(summary(scm_pca_models[[s]]))
  cat("\nVIF:\n")
  print(vif(scm_pca_models[[s]]))
  cat("\nR2:\n")
  print(r2(scm_pca_models[[s]]))
})

### ---- Extract R2 and p-values per scenario x climate variable ----
scm_pca_stats <- lapply(scm_scenarios, function(s) {
  mod     <- scm_pca_models[[s]]
  r2_vals <- r2(mod)
  coefs   <- broom.mixed::tidy(mod, effects = "fixed")
  
  p_PC1  <- coefs %>% filter(term=="scale(PC1)") %>% pull(p.value)
  p_PC2  <- coefs %>% filter(term=="scale(PC2)") %>% pull(p.value)
  
  data.frame(
    scenario       = s,
    Marginal_R2    = round(r2_vals$R2_marginal,    2),
    Conditional_R2 = round(r2_vals$R2_conditional, 2),
    p_PC1          = p_PC1,
    p_PC2          = p_PC2
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

### ---- Predictions per scenario ----
predict_pca_scm <- function(s, climate_var) {
  ggpredict(scm_pca_models[[s]], terms = paste0(climate_var, " [all]")) %>%
    as.data.frame() %>%
    mutate(scenario = s, climate_var = climate_var)
}

scm_pca_preds <- bind_rows(
  lapply(scm_scenarios, predict_pca_scm, climate_var = "PC1"),
  lapply(scm_scenarios, predict_pca_scm, climate_var = "PC2")
) %>%
  left_join(
    scm_pca_stats %>% select(scenario, legend_label, p_PC1, p_PC2),
    by = "scenario"
  ) %>%
  mutate(
    sig = case_when(
      climate_var == "PC1"             & p_PC1  < 0.05 ~ "significant",
      climate_var == "PC2"             & p_PC2  < 0.05 ~ "significant",
      TRUE ~ "ns"
    ),
    climate_var = factor(
      climate_var,
      levels = c("PC1", "PC2"),
      labels = c("Principal Component 1", "Principal Component 2")
    ),
    legend_label = factor(legend_label, levels = scm_pca_stats$legend_label)
  )

### ---- Plot: facet_wrap(~ climate_var), lines colored by scenario ----
p_pca_scm <- ggplot(scm_pca_preds, aes(x = x, y = predicted,
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
  "analysis/network/figures/SCM_Robustness_by_PCAclimate.pdf",
  p_pca_scm,
  width = 11, height = 5
)

## *******************************************************************
# ---- Network Structural Properties Models ----
# These models characterize how eight network structural metrics
# co-vary with climate predictors. Following Ponisio (2020), these
# properties (redundancy, complementarity, generalization, richness)
# are theorized to determine a community's capacity to withstand
# perturbation, so they reveal the structural landscape that robustness
# results emerge from. They are NOT a direct measure of resistance.
## *******************************************************************
## ---- Order sites by latitude ----
site_order <- geography %>%
  group_by(Site) %>%
  summarize(
    Lat = mean(Lat, na.rm = TRUE)
  ) %>%
  arrange(desc(Lat)) %>%   # north -> south
  pull(Site)

cor.dats <- cor.dats %>%
  separate(
    Site,
    into = c("Site", "Year", "SampleRound"),
    sep = "\\.",
    convert = TRUE,
    remove = FALSE
  ) %>%
  left_join(
    exposure %>%
      select(
        Site,
        Year,
        SampleRound,
        APi,
        GDD,
        PC1,
        PC2
      ),
    by = c("Site", "Year", "SampleRound")
  ) %>%
  mutate(
    Site = factor(Site, levels = site_order),
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

## ---- Fit one model per structural metric x climate variable ----
fit_structural_models <- function(response_vars, climate_var, data) {
  mods <- lapply(response_vars, function(y) {
    formula <- as.formula(
      paste0(y, " ~ scale(", climate_var, ") + (1 | Site) + (1 | Year)")
    )
    lmer(formula, data = data, REML = FALSE)
  })
  names(mods) <- response_vars
  return(mods)
}

mods_api  <- fit_structural_models(metrics, "APi", cor.dats)
mods_gdd <- fit_structural_models(metrics, "GDD", cor.dats)

## ---- Extract results ----
extract_structural_results <- function(mods, climate_var) {
  lapply(names(mods), function(name) {
    broom.mixed::tidy(mods[[name]], effects = "fixed") %>%
      filter(term == paste0("scale(", climate_var, ")")) %>%
      mutate(Metric = name, ClimateVar = climate_var)
  }) %>%
    bind_rows()
}

results_all <- bind_rows(
  extract_structural_results(mods_api,  "APi"),
  extract_structural_results(mods_gdd, "GDD")
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
r2_gdd <- extract_r2(mods_gdd, "GDD")

r2_all <- bind_rows(r2_api, r2_gdd)

results_all <- results_all %>%
  left_join(r2_all, by = c("Metric", "ClimateVar"))

sig_table <- results_all %>%
  mutate(signif = if_else(p.value < 0.05, "significant", "ns")) %>%
  select(Metric, ClimateVar, signif)

## ---- Predictions and Plot structural property effects ----
r2_labels <- r2_all %>%
  mutate(
    r2_text = paste0("R²m = ", round(Marginal_R2, 2),
                     "\nR²c = ", round(Conditional_R2, 2))
  ) %>%
  select(Metric, ClimateVar, r2_text)


plot_structural_effects <- function(mods, climate_var, data) {
  
  # Predictions via ggpredict (same method used for the PCA plots below),
  # holding other model terms at their mean
  preds <- lapply(names(mods), function(m) {
    ggpredict(mods[[m]], terms = paste0(climate_var, " [all]")) %>%
      as.data.frame() %>%
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
  
  r2_text_by_panel <- preds %>% distinct(metric, r2_text)
  
  ggplot(preds, aes(x = x, y = predicted)) +
    geom_ribbon(
      aes(ymin = conf.low, ymax = conf.high),
      fill = "gray80", alpha = 0.4
    ) +
    geom_line(aes(linetype = signif), color = "black", linewidth = 1) +
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
      data = r2_text_by_panel,
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
      y = "Network structural property",
      color = "Site (north to south)"
    ) +
    guides(shape = "none")
}

p_api  <- plot_structural_effects(mods_api,  "APi",            cor.dats)
p_gdd  <- plot_structural_effects(mods_gdd, "GDD", cor.dats)

ggsave("analysis/network/figures/NetworkStructure_APi.pdf",
       p_api,  width = 11, height = 8)

ggsave("analysis/network/figures/NetworkStructure_GDD.pdf",
       p_gdd, width = 11, height = 8)

## *******************************************************************
# ---- Network Structural Properties: PCA Exposure Models ----
## *******************************************************************
## ---- Fit one model that includes both PCs simultaneously x climate variable ----
fit_structural_models_pca <- function(response_vars, data){
  mods <- lapply(response_vars, function(y){
    formula <- as.formula(
      paste0(
        y,
        " ~ scale(PC1) + scale(PC2) + (1 | Site) + (1 | Year)"
      )
    )
    lmer(formula, data=data, REML=FALSE)
  })
  names(mods) <- response_vars
  mods
}

mods_pca <- fit_structural_models_pca(metrics, cor.dats)

## ---- Extract results ----
extract_structural_results_pca <- function(mods){
  lapply(names(mods), function(name){
    tidy(mods[[name]], effects="fixed") %>%
      filter(term %in% c("scale(PC1)", "scale(PC2)")) %>%
      mutate(Metric=name)
  }) %>%
    bind_rows()
}

results_pca <- extract_structural_results_pca(mods_pca)

### ---- Extract r2 ----
r2_pca <- extract_r2(mods_pca, "PCA")

## ---- Predictions ----
predict_structural_pca <- function(metric, climate_var){
  ggpredict(
    mods_pca[[metric]],
    terms = paste0(climate_var, " [all]")
  ) %>%
    as.data.frame() %>%
    mutate(
      metric = metric,
      climate_var = climate_var
    )
}

preds_pca <- bind_rows(
  lapply(metrics,
         predict_structural_pca,
         climate_var = "PC1"),
  lapply(metrics,
         predict_structural_pca,
         climate_var = "PC2")
)

### ---- Add significance ----
sig_pca <- results_pca %>%
  mutate(
    climate_var = case_when(
      term == "scale(PC1)" ~ "PC1",
      term == "scale(PC2)" ~ "PC2"
    ),
    signif = if_else(p.value < 0.05,
                     "significant",
                     "ns")
  ) %>%
  select(Metric, climate_var, signif)

preds_pca <- preds_pca %>%
  left_join(
    sig_pca,
    by = c("metric" = "Metric",
           "climate_var")
  )

## ---- Plot structural property effects ----
# R2 labels (same joint model underlies both PC1 and PC2 panels)
r2_pca_labels <- r2_pca %>%
  mutate(
    r2_text = paste0(
      "R²m = ", round(Marginal_R2, 2),
      "\nR²c = ", round(Conditional_R2, 2)
    )
  ) %>%
  select(Metric, r2_text)

# Raw data in long format, ordered/labeled the same way as the APi/GDD plots
raw_pca <- cor.dats %>%
  pivot_longer(
    cols = all_of(metrics),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = names(metric_labels)))

plot_structural_effects_pca <- function(climate_var){
  
  preds <- preds_pca %>%
    filter(climate_var == !!climate_var) %>%
    mutate(metric = factor(metric, levels = names(metric_labels))) %>%
    left_join(r2_pca_labels, by = c("metric" = "Metric"))
  
  r2_text_by_panel <- preds %>% distinct(metric, r2_text)
  
  ggplot(preds, aes(x = x, y = predicted)) +
    
    geom_ribbon(
      aes(ymin = conf.low, ymax = conf.high),
      fill = "gray80", alpha = 0.4
    ) +
    
    geom_line(
      aes(linetype = signif),
      color = "black", linewidth = 1
    ) +
    
    geom_point(
      data = raw_pca,
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
      data = r2_text_by_panel,
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
      y = "Network structural property",
      color = "Site (north to south)"
    ) +
    
    guides(shape = "none")
}

p_pc1 <- plot_structural_effects_pca("PC1")
p_pc2 <- plot_structural_effects_pca("PC2")

ggsave("analysis/network/figures/NetworkStructure_PC1.pdf",
       p_pc1, width = 11, height = 8)

ggsave("analysis/network/figures/NetworkStructure_PC2.pdf",
       p_pc2, width = 11, height = 8)

## *******************************************************************
# ---- Network Structural Properties: PC1 + Seasonal GDD Model ----
## *******************************************************************
## Does this season's climate anomaly (PC1) predict network structure,
## once we've accounted for how far along the season is
## (GDD, allowed to bend via a quadratic term)?
## ---- Fit one model per structural metric: PC1 + poly(GDD, 2) ----
fit_structural_models_pc1gdd <- function(response_vars, data){
  mods <- lapply(response_vars, function(y){
    formula <- as.formula(
      paste0(
        y,
        " ~ scale(PC1) + poly(GDD, 2) + (1 | Site) + (1 | Year)"
      )
    )
    lmer(formula, data = data, REML = FALSE)
  })
  names(mods) <- response_vars
  mods
}

mods_pc1gdd <- fit_structural_models_pc1gdd(metrics, cor.dats)

## ---- Extract results ----
extract_structural_results_pc1gdd <- function(mods){
  lapply(names(mods), function(name){
    tidy(mods[[name]], effects = "fixed") %>%
      filter(term %in% c("scale(PC1)", "poly(GDD, 2)1", "poly(GDD, 2)2")) %>%
      mutate(Metric = name)
  }) %>%
    bind_rows()
}

results_pc1gdd <- extract_structural_results_pc1gdd(mods_pc1gdd)

### ---- Extract r2 ----
r2_pc1gdd <- extract_r2(mods_pc1gdd, "PC1_GDD")

## ---- Predictions ----
# For the PC1 panel: PC1 varied across its range, GDD held at its mean
# For the GDD panel: GDD varied across its range, PC1 held at its mean
predict_structural_pc1gdd <- function(metric, climate_var){
  ggpredict(
    mods_pc1gdd[[metric]],
    terms = paste0(climate_var, " [all]")
  ) %>%
    as.data.frame() %>%
    mutate(
      metric = metric,
      climate_var = climate_var
    )
}

preds_pc1gdd <- bind_rows(
  lapply(metrics,
         predict_structural_pc1gdd,
         climate_var = "PC1"),
  lapply(metrics,
         predict_structural_pc1gdd,
         climate_var = "GDD")
)

### ---- Add significance ----
# PC1 significance comes from the single scale(PC1) term.
# GDD significance is flagged if either the linear or quadratic
# poly(GDD, 2) term is significant, since "GDD's effect" now spans two terms.
sig_pc1gdd <- results_pc1gdd %>%
  mutate(
    climate_var = case_when(
      term == "scale(PC1)" ~ "PC1",
      term %in% c("poly(GDD, 2)1", "poly(GDD, 2)2") ~ "GDD"
    )
  ) %>%
  group_by(Metric, climate_var) %>%
  summarise(
    signif = if_else(any(p.value < 0.05), "significant", "ns"),
    .groups = "drop"
  )

preds_pc1gdd <- preds_pc1gdd %>%
  left_join(
    sig_pc1gdd,
    by = c("metric" = "Metric",
           "climate_var")
  )

## ---- Plot structural property effects ----
r2_pc1gdd_labels <- r2_pc1gdd %>%
  mutate(
    r2_text = paste0(
      "R²m = ", round(Marginal_R2, 2),
      "\nR²c = ", round(Conditional_R2, 2)
    )
  ) %>%
  select(Metric, r2_text)

raw_pc1gdd <- cor.dats %>%
  pivot_longer(
    cols = all_of(metrics),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = names(metric_labels)))

plot_structural_effects_pc1gdd <- function(climate_var){
  
  preds <- preds_pc1gdd %>%
    filter(climate_var == !!climate_var) %>%
    mutate(metric = factor(metric, levels = names(metric_labels))) %>%
    left_join(r2_pc1gdd_labels, by = c("metric" = "Metric"))
  
  r2_text_by_panel <- preds %>% distinct(metric, r2_text)
  
  ggplot(preds, aes(x = x, y = predicted)) +
    
    geom_ribbon(
      aes(ymin = conf.low, ymax = conf.high),
      fill = "gray80", alpha = 0.4
    ) +
    
    geom_line(
      aes(linetype = signif),
      color = "black", linewidth = 1
    ) +
    
    geom_point(
      data = raw_pc1gdd,
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
      data = r2_text_by_panel,
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
      y = "Network structural property",
      color = "Site (north to south)"
    ) +
    
    guides(shape = "none")
}

p_pc1_adj <- plot_structural_effects_pc1gdd("PC1")
p_gdd_adj <- plot_structural_effects_pc1gdd("GDD")

ggsave("analysis/network/figures/NetworkStructure_PC1_adjGDD.pdf",
       p_pc1_adj, width = 11, height = 8)

ggsave("analysis/network/figures/NetworkStructure_GDD_adjPC1.pdf",
       p_gdd_adj, width = 11, height = 8)
