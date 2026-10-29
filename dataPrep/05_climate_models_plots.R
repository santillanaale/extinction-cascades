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

# ---- Load data ----
setwd(file.path(local.path, "extinction-cascades"))
load("analysis/network/robustness_metrics/YearSR_robustness_climate.Rdata")

setwd(file.path(local.path, "skyIslands_saved"))
geography <- read.csv("data/relational/original/geography.csv")

setwd(file.path(local.path, "extinction-cascades"))
exposure <- read_csv("data/exposure.csv")
# Load audited Q2 metrics (61 networks, 52 columns)
cor.dats <- readRDS(
  "analysis/network/saved/q2_metrics_table_expanded_61_N100.rds"
)

stopifnot(
  nrow(cor.dats) == 61,
  ncol(cor.dats) == 52,
  !anyDuplicated(cor.dats$Network)
)

## *******************************************************************
# ---- Helper functions ----
## *******************************************************************

## Fit one lmer per response variable, same climate formula each time.
## `response` can be a vector (structural metrics) or a single string
## (TCM/SCM, used via fit_scenario_models below).
fit_models <- function(response, data, climate_formula){
  mods <- lapply(response, function(y){
    form <- as.formula(
      paste0(y, " ~ ", climate_formula, " + (1 | Year) + (1 | Site)")
    )
    lmer(form, data = data, REML = FALSE,
         control = lmerControl(optimizer = "bobyqa"))
  })
  names(mods) <- response
  mods
}

## Fit one model per scenario (filters data to scenario == s, fits a
## single response variable). Used for TCM/SCM robustness models.
fit_scenario_models <- function(response, data, scenarios,
                                climate_formula, scenario_col = "scenario"){
  mods <- lapply(scenarios, function(s){
    dat <- filter(data, .data[[scenario_col]] == s)
    fit_models(response, dat, climate_formula)[[1]]
  })
  names(mods) <- scenarios
  mods
}

## Full fixed-effects table, for CSV export.
extract_results <- function(mods, model_name, id_col = "id"){
  bind_rows(lapply(names(mods), function(nm){
    broom.mixed::tidy(mods[[nm]], effects = "fixed") %>%
      mutate(!!id_col := nm, model = model_name)
  }))
}

## R2 table, for CSV export.
extract_r2 <- function(mods, model_name, id_col = "id"){
  bind_rows(lapply(names(mods), function(nm){
    r <- performance::r2(mods[[nm]])
    tibble(
      !!id_col := nm,
      model = model_name,
      Marginal_R2 = r$R2_marginal,
      Conditional_R2 = r$R2_conditional
    )
  }))
}

## R2 + p-value (min across any term matching that climate variable,
## so poly(GDD, 2)'s linear + quadratic terms both count) for one model.
get_stats <- function(mod, climate_vars){
  r2_vals <- performance::r2(mod)
  coefs   <- broom.mixed::tidy(mod, effects = "fixed")
  bind_rows(lapply(climate_vars, function(cv){
    terms_cv <- coefs %>% filter(str_detect(term, fixed(cv)))
    tibble(
      climate_var    = cv,
      p_min          = if (nrow(terms_cv) > 0) min(terms_cv$p.value, na.rm = TRUE) else NA_real_,
      Marginal_R2    = r2_vals$R2_marginal,
      Conditional_R2 = r2_vals$R2_conditional
    )
  }))
}

## get_stats() across every model in a named list (scenarios or metrics).
get_stats_all <- function(mods, climate_vars){
  bind_rows(lapply(names(mods), function(nm){
    get_stats(mods[[nm]], climate_vars) %>% mutate(id = nm)
  })) %>%
    mutate(signif = if_else(p_min < 0.05, "significant", "ns"))
}

## ggpredict() across every model x every climate variable.
predict_all <- function(mods, climate_vars){
  bind_rows(lapply(names(mods), function(nm){
    bind_rows(lapply(climate_vars, function(cv){
      ggpredict(mods[[nm]], terms = paste0(cv, " [all]")) %>%
        as.data.frame() %>%
        mutate(id = nm, climate_var = cv)
    }))
  }))
}

## TCM/SCM style plot: facet by climate_var, color/legend by scenario
## (legend label carries R2), no raw data overlay.
plot_scenario_climate <- function(mods, climate_vars, var_labels, y_lab, out_file,
                                  width = 11, height = 5){
  
  stats <- get_stats_all(mods, climate_vars) %>%
    mutate(
      legend_label = paste0(
        id, "\n  R²m = ", round(Marginal_R2, 2),
        ", R²c = ", round(Conditional_R2, 2)
      )
    )
  
  preds <- predict_all(mods, climate_vars) %>%
    left_join(
      stats %>% select(id, climate_var, legend_label, signif),
      by = c("id", "climate_var")
    ) %>%
    mutate(
      legend_label = factor(legend_label, levels = unique(stats$legend_label)),
      climate_var  = factor(climate_var, levels = names(var_labels), labels = var_labels)
    )
  
  p <- ggplot(preds, aes(x = x, y = predicted, color = legend_label, fill = legend_label)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
    geom_line(aes(linetype = signif), linewidth = 1) +
    scale_linetype_manual(values = c(significant = "solid", ns = "dashed"), guide = "none") +
    facet_wrap(~ climate_var, scales = "free_x") +
    theme_minimal(base_size = 13) +
    labs(x = NULL, y = y_lab, color = "Extinction ordering", fill = "Extinction ordering") +
    theme(
      strip.text        = element_text(face = "bold"),
      panel.grid.minor  = element_blank(),
      legend.text       = element_text(size = 9),
      legend.key.height = unit(1.2, "lines")
    )
  
  ggsave(out_file, p, width = width, height = height)
  p
}

## Structural-property style plot: facet by metric, one climate_var at a
## time, raw data points colored by Site, R2 text annotation per panel.
plot_metric_climate <- function(mods, climate_var, data, metrics, metric_labels,
                                climate_lab = climate_var, out_file,
                                width = 11, height = 8, nrow = 2){
  
  stats <- get_stats_all(mods, climate_var) %>%
    mutate(r2_text = paste0("R²m = ", round(Marginal_R2, 2),
                            "\nR²c = ", round(Conditional_R2, 2)))
  
  preds_all <- predict_all(mods, climate_var) %>%
    left_join(stats %>% select(id, signif, r2_text), by = "id")
  
  raw_all <- data %>%
    pivot_longer(cols = all_of(metrics), names_to = "id", values_to = "value")
  
  ## Build one standalone ggplot per structural property, each with its own
  ## real x-axis title (climate variable, same text on every panel) and
  ## y-axis title (that property's label). No facet_wrap, so there is no
  ## strip/facet title to remove and no risk of shared/collapsed titles.
  panel_list <- lapply(names(metric_labels), function(m){
    
    preds  <- filter(preds_all, id == m)
    raw    <- filter(raw_all, id == m)
    r2_lab <- preds$r2_text[1]
    
    ggplot(preds, aes(x = x, y = predicted)) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "gray80", alpha = 0.4) +
      geom_line(aes(linetype = signif), color = "black", linewidth = 1) +
      geom_point(
        data = raw,
        aes(x = .data[[climate_var]], y = value, color = Site, shape = SampleRound),
        inherit.aes = FALSE, alpha = 0.6
      ) +
      scale_linetype_manual(values = c(significant = "solid", ns = "dashed"), guide = "none") +
      # annotate("text", x = -Inf, y = Inf, label = r2_lab,
      # hjust = -0.05, vjust = 1.1, size = 3) +
      theme_minimal(base_size = 14) +
      labs(x = climate_lab, y = metric_labels[[m]], color = "Site (N to S)") +
      guides(shape = "none")
  })
  
  p <- wrap_plots(panel_list, nrow = nrow, guides = "collect")
  
  ggsave(out_file, p, width = width, height = height)
  p
}

## Optional diagnostics (summary/VIF/R2 per model) - uncomment to run.
print_diagnostics <- function(mods, label){
  for (nm in names(mods)){
    cat("\n\n============================\n")
    cat(label, ":", nm, "\n")
    cat("============================\n")
    print(summary(mods[[nm]]))
    cat("\nVIF:\n"); print(vif(mods[[nm]]))
    cat("\nR2:\n");  print(r2(mods[[nm]]))
  }
}
# lapply(names(tcm_pc1gdd_models), function(s) print_diagnostics(tcm_pc1gdd_models[s], "TCM"))


## *******************************************************************
# ---- TCM robustness ----
## *******************************************************************
tcm_climate <- tcm_results %>%
  left_join(exposure, by = c("Site", "Year", "SampleRound"))

tcm_scenarios <- c("abundance_low", "abundance_high", "degree_low", "degree_high")

tcm_api_models    <- fit_scenario_models("Robustness", tcm_climate, tcm_scenarios, "scale(APi)+scale(GDD)")
tcm_pca_models    <- fit_scenario_models("Robustness", tcm_climate, tcm_scenarios, "scale(PC1)+scale(PC2)")
tcm_pc1gdd_models <- fit_scenario_models("Robustness", tcm_climate, tcm_scenarios, "scale(PC1)+poly(GDD,2)")

p_tcm_api <- plot_scenario_climate(
  tcm_api_models, c("APi", "GDD"),
  var_labels = c(APi = "Antecedent Precipitation Index", GDD = "Growing Degree Days"),
  y_lab = "TCM Robustness",
  out_file = "analysis/network/figures/TCM_Robustness_by_climate.pdf"
)

p_tcm_pca <- plot_scenario_climate(
  tcm_pca_models, c("PC1", "PC2"),
  var_labels = c(PC1 = "Principal Component 1", PC2 = "Principal Component 2"),
  y_lab = "TCM Robustness",
  out_file = "analysis/network/figures/TCM_Robustness_by_PCAclimate.pdf"
)

p_tcm_pc1gdd <- plot_scenario_climate(
  tcm_pc1gdd_models, c("PC1", "GDD"),
  var_labels = c(PC1 = "Principal Component 1", GDD = "Growing Degree Days (seasonal progression)"),
  y_lab = "TCM Robustness",
  out_file = "analysis/network/figures/TCM_Robustness_by_PC1GDD.pdf"
)


## *******************************************************************
# ---- SCM robustness ----
## *******************************************************************
scm_climate <- scm_results %>%
  left_join(exposure, by = c("Site", "Year", "SampleRound"))

scm_scenarios <- c("scm_random_plant", "scm_dominant_plant")

scm_api_models    <- fit_scenario_models("SCM_Robustness", scm_climate, scm_scenarios, "scale(APi)+scale(GDD)")
scm_pca_models    <- fit_scenario_models("SCM_Robustness", scm_climate, scm_scenarios, "scale(PC1)+scale(PC2)")
scm_pc1gdd_models <- fit_scenario_models("SCM_Robustness", scm_climate, scm_scenarios, "scale(PC1)+poly(GDD,2)")

p_scm_api <- plot_scenario_climate(
  scm_api_models, c("APi", "GDD"),
  var_labels = c(APi = "Antecedent Precipitation Index", GDD = "Growing Degree Days"),
  y_lab = "SCM Robustness",
  out_file = "analysis/network/figures/SCM_Robustness_by_climate.pdf"
)

p_scm_pca <- plot_scenario_climate(
  scm_pca_models, c("PC1", "PC2"),
  var_labels = c(PC1 = "Principal Component 1", PC2 = "Principal Component 2"),
  y_lab = "SCM Robustness",
  out_file = "analysis/network/figures/SCM_Robustness_by_PCAclimate.pdf"
)

p_scm_pc1gdd <- plot_scenario_climate(
  scm_pc1gdd_models, c("PC1", "GDD"),
  var_labels = c(PC1 = "Principal Component 1", GDD = "Growing Degree Days (seasonal progression)"),
  y_lab = "SCM Robustness",
  out_file = "analysis/network/figures/SCM_Robustness_by_PC1GDD.pdf"
)


## *******************************************************************
# ---- Network structural properties ----
## *******************************************************************
site_order <- geography %>%
  group_by(Site) %>%
  summarize(Lat = mean(Lat, na.rm = TRUE)) %>%
  arrange(desc(Lat)) %>%
  pull(Site)

cor.dats <- cor.dats %>%
  left_join(
    exposure %>% select(Site, Year, SampleRound, APi, GDD, PC1, PC2),
    by = c("Site", "Year", "SampleRound")
  ) %>%
  mutate(
    Site = factor(Site, levels = site_order),
    Year = factor(Year),
    SampleRound = factor(SampleRound)
  )

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

mods_api    <- fit_models(metrics, cor.dats, "scale(APi)")
mods_gdd    <- fit_models(metrics, cor.dats, "scale(GDD)")
mods_pca    <- fit_models(metrics, cor.dats, "scale(PC1)+scale(PC2)")
mods_pc1gdd <- fit_models(metrics, cor.dats, "scale(PC1)+poly(GDD,2)")

p_api <- plot_metric_climate(
  mods_api, "APi", cor.dats, metrics, metric_labels,
  out_file = "analysis/network/figures/NetworkStructure_APi.pdf"
)

p_gdd <- plot_metric_climate(
  mods_gdd, "GDD", cor.dats, metrics, metric_labels,
  out_file = "analysis/network/figures/NetworkStructure_GDD.pdf"
)

p_pc1 <- plot_metric_climate(
  mods_pca, "PC1", cor.dats, metrics, metric_labels,
  out_file = "analysis/network/figures/NetworkStructure_PC1.pdf"
)

p_pc2 <- plot_metric_climate(
  mods_pca, "PC2", cor.dats, metrics, metric_labels,
  out_file = "analysis/network/figures/NetworkStructure_PC2.pdf"
)

p_pc1_adj <- plot_metric_climate(
  mods_pc1gdd, "PC1", cor.dats, metrics, metric_labels,
  out_file = "analysis/network/figures/NetworkStructure_PC1_adjGDD.pdf"
)

p_gdd_adj <- plot_metric_climate(
  mods_pc1gdd, "GDD", cor.dats, metrics, metric_labels,
  out_file = "analysis/network/figures/NetworkStructure_GDD_adjPC1.pdf"
)


## *******************************************************************
# ---- Functional diversity (bee community, SiteYearSR level) ----
## *******************************************************************
## Bee functional divergense (BeeFDiv), dispersion (BeeFDis) and evenness (BeeFEve) come from
## src/calcFuncUniqOrig.R / 5_funcTraits.R and are stored on spec.net,
## duplicated across every specimen row collected in that site-year-
## round. Collapse to one row per SiteYearSR before joining to climate,
## or the models below would treat every bee as an independent
## observation of its community's functional diversity.

setwd(file.path(local.path, "skyIslands"))
load("data/spec_traits.Rdata")  # spec.net w/ BeeFDisSiteYearSR, BeeFEveSiteYearSR, BeeFDivSiteYearSR
setwd(file.path(local.path, "extinction-cascades"))

funcdiv <- spec.net %>%
  distinct(SiteYearSR, Site, Year, SampleRound,
           BeeFDisSiteYearSR, BeeFEveSiteYearSR, BeeFDivSiteYearSR) %>%
  mutate(
    Year        = as.numeric(Year),
    SampleRound = as.numeric(SampleRound)
  )

funcdiv_climate <- funcdiv %>%
  left_join(
    exposure %>% select(Site, Year, SampleRound, APi, GDD, PC1, PC2),
    by = c("Site", "Year", "SampleRound")
  ) %>%
  mutate(
    Site        = factor(Site, levels = site_order),
    Year        = factor(Year),
    SampleRound = factor(SampleRound)
  )

fd_metrics <- c("BeeFDisSiteYearSR", "BeeFEveSiteYearSR", "BeeFDivSiteYearSR")

fd_metric_labels <- c(
  "BeeFDisSiteYearSR" = "Bee Functional Dispersion",
  "BeeFEveSiteYearSR" = "Bee Functional Evenness",
  "BeeFDivSiteYearSR" = "Bee Functional Divergence"
)

mods_fd_api    <- fit_models(fd_metrics, funcdiv_climate, "scale(APi)")
mods_fd_gdd    <- fit_models(fd_metrics, funcdiv_climate, "scale(GDD)")
mods_fd_pca    <- fit_models(fd_metrics, funcdiv_climate, "scale(PC1)+scale(PC2)")
mods_fd_pc1gdd <- fit_models(fd_metrics, funcdiv_climate, "scale(PC1)+poly(GDD,2)")

p_fd_api <- plot_metric_climate(
  mods_fd_api, "APi", funcdiv_climate, fd_metrics, fd_metric_labels,
  out_file = "analysis/network/figures/FuncDiv_APi.pdf",
  width = 11, height = 4, nrow = 1
)

p_fd_gdd <- plot_metric_climate(
  mods_fd_gdd, "GDD", funcdiv_climate, fd_metrics, fd_metric_labels,
  out_file = "analysis/network/figures/FuncDiv_GDD.pdf",
  width = 11, height = 4, nrow = 1
)

p_fd_pc1 <- plot_metric_climate(
  mods_fd_pca, "PC1", funcdiv_climate, fd_metrics, fd_metric_labels,
  climate_lab = "Principal Component 1",
  out_file = "analysis/network/figures/FuncDiv_PC1.pdf",
  width = 11, height = 4, nrow = 1
)

p_fd_pc2 <- plot_metric_climate(
  mods_fd_pca, "PC2", funcdiv_climate, fd_metrics, fd_metric_labels,
  climate_lab = "Principal Component 2",
  out_file = "analysis/network/figures/FuncDiv_PC2.pdf",
  width = 11, height = 4, nrow = 1
)

p_fd_pc1_adj <- plot_metric_climate(
  mods_fd_pc1gdd, "PC1", funcdiv_climate, fd_metrics, fd_metric_labels,
  climate_lab = "Principal Component 1 (adj. for GDD)",
  out_file = "analysis/network/figures/FuncDiv_PC1_adjGDD.pdf",
  width = 11, height = 4, nrow = 1
)

p_fd_gdd_adj <- plot_metric_climate(
  mods_fd_pc1gdd, "GDD", funcdiv_climate, fd_metrics, fd_metric_labels,
  climate_lab = "Growing Degree Days (adj. for PC1)",
  out_file = "analysis/network/figures/FuncDiv_GDD_adjPC1.pdf",
  width = 11, height = 4, nrow = 1
)


## *******************************************************************
# ---- Export model summaries ----
## *******************************************************************
model_results <- bind_rows(
  extract_results(tcm_api_models,    "TCM_APi_GDD",   "scenario"),
  extract_results(tcm_pca_models,    "TCM_PCA",       "scenario"),
  extract_results(tcm_pc1gdd_models, "TCM_PC1_GDD",   "scenario"),
  
  extract_results(scm_api_models,    "SCM_APi_GDD",   "scenario"),
  extract_results(scm_pca_models,    "SCM_PCA",       "scenario"),
  extract_results(scm_pc1gdd_models, "SCM_PC1_GDD",   "scenario"),
  
  extract_results(mods_api,    "Network_APi",     "Metric"),
  extract_results(mods_gdd,    "Network_GDD",     "Metric"),
  extract_results(mods_pca,    "Network_PCA",     "Metric"),
  extract_results(mods_pc1gdd, "Network_PC1_GDD", "Metric"),
  
  extract_results(mods_fd_api,    "FuncDiv_APi",     "Metric"),
  extract_results(mods_fd_gdd,    "FuncDiv_GDD",     "Metric"),
  extract_results(mods_fd_pca,    "FuncDiv_PCA",     "Metric"),
  extract_results(mods_fd_pc1gdd, "FuncDiv_PC1_GDD", "Metric")
)

model_r2 <- bind_rows(
  extract_r2(tcm_api_models,    "TCM_APi_GDD",   "scenario"),
  extract_r2(tcm_pca_models,    "TCM_PCA",       "scenario"),
  extract_r2(tcm_pc1gdd_models, "TCM_PC1_GDD",   "scenario"),
  
  extract_r2(scm_api_models,    "SCM_APi_GDD",   "scenario"),
  extract_r2(scm_pca_models,    "SCM_PCA",       "scenario"),
  extract_r2(scm_pc1gdd_models, "SCM_PC1_GDD",   "scenario"),
  
  extract_r2(mods_api,    "Network_APi",     "Metric"),
  extract_r2(mods_gdd,    "Network_GDD",     "Metric"),
  extract_r2(mods_pca,    "Network_PCA",     "Metric"),
  extract_r2(mods_pc1gdd, "Network_PC1_GDD", "Metric"),
  
  extract_r2(mods_fd_api,    "FuncDiv_APi",     "Metric"),
  extract_r2(mods_fd_gdd,    "FuncDiv_GDD",     "Metric"),
  extract_r2(mods_fd_pca,    "FuncDiv_PCA",     "Metric"),
  extract_r2(mods_fd_pc1gdd, "FuncDiv_PC1_GDD", "Metric")
)

write_csv(model_results, "analysis/network/model_fixed_effects.csv")
write_csv(model_r2, "analysis/network/model_R2.csv")

## *******************************************************************
# ---- Check for Convergence ----
## *******************************************************************
check_convergence <- function(mods, label){
  msgs <- lapply(mods, function(m) m@optinfo$conv$lme4$messages)
  tibble(
    model_set = label,
    id = names(mods),
    message = sapply(msgs, function(x) if (is.null(x)) "OK" else paste(x, collapse = " | ")),
    failed_to_converge = sapply(msgs, function(x) any(grepl("failed to converge", x))),
    singular = sapply(msgs, function(x) any(grepl("singular", x)))
  )
}

convergence_report <- bind_rows(
  check_convergence(tcm_api_models,    "TCM_APi_GDD"),
  check_convergence(tcm_pca_models,    "TCM_PCA"),
  check_convergence(tcm_pc1gdd_models, "TCM_PC1_GDD"),
  
  check_convergence(scm_api_models,    "SCM_APi_GDD"),
  check_convergence(scm_pca_models,    "SCM_PCA"),
  check_convergence(scm_pc1gdd_models, "SCM_PC1_GDD"),
  
  check_convergence(mods_api,    "Network_APi"),
  check_convergence(mods_gdd,    "Network_GDD"),
  check_convergence(mods_pca,    "Network_PCA"),
  check_convergence(mods_pc1gdd, "Network_PC1_GDD"),
  
  check_convergence(mods_fd_api,    "FuncDiv_APi"),
  check_convergence(mods_fd_gdd,    "FuncDiv_GDD"),
  check_convergence(mods_fd_pca,    "FuncDiv_PCA"),
  check_convergence(mods_fd_pc1gdd, "FuncDiv_PC1_GDD")
)

# Just the real problems (non-convergence, not just singularity)
convergence_report %>% filter(failed_to_converge)