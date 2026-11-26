# 05_climate_models.R
# Climate (PC1, PC2) effects on:
#   Q1 - community abundance, richness, and Shannon diversity
#   Q2 - plant-pollinator network metrics
#   Q3 - extinction-cascade robustness
# Functions live in 05_climate_models_helpers.R.

rm(list = ls())
setwd("~/")
source("lab_paths.R")

library(dplyr)
library(lme4)
library(performance)
library(see)
library(ggplot2)
library(qqplotr)

dir.SI <- file.path(local.path, "skyIslands")
dir.EC <- file.path(local.path, "extinction-cascades")
dir.SI_saved <- file.path(local.path, "skyIslands_saved")

source(file.path(dir.EC, "dataPrep", "src", "05_climate_models_helpers.R"))

keys <- c("Site", "Year", "SampleRound")


# 0. Configuration ---------------------------------------------------------
# Everything you might edit (responses, labels, panel order) lives here.

## Q1 responses
q1_abundance_ys <- c(
  "BeeAbundance", "SyrphidAbundance", "HBAbundance", "BombusAbundance",
  "NonBombusHBAbundance", "ButterflyHawkAbundance", "BeeSyrphLepAbundance",
  "MeanFloralAbundance", "MeanFloweringPlantAbundance"
)

q1_richness_ys <- c(
  "BeeRichness", "SyrphidRichness", "BombusRichness",
  "ButterflyHawkRichness", "BeeSyrphLepRichness", "MeanFloralRichness"
)

# Each Shannon response is paired with the richness of the same group.
shannon_richness <- c(
  BeeDiversity = "BeeRichness",
  SyrphidDiversity = "SyrphidRichness",
  BombusDiversity = "BombusRichness",
  ButterflyHawkDiversity = "ButterflyHawkRichness",
  BeeSyrphLepDiversity = "BeeSyrphLepRichness",
  MeanFloralDiversity = "MeanFloralRichness"
)

q1_shannon_ys <- c(names(shannon_richness), "MeanFloweringPlantDiversity")

q1_fd_ys <- c("BeeFDis", "BeeFEve", "BeeFDiv")

q1_ys <- c(q1_abundance_ys, q1_richness_ys, q1_shannon_ys, q1_fd_ys)

# Q1 plot labels: "Bee" in "BeeAbundance" is the response key.
response_labels <- c(
  Bee = "Bees",
  Syrphid = "Hoverflies",
  HB = "Honey bees",
  Bombus = "Bumble bees",
  NonBombusHB = "Other bees",
  ButterflyHawk = "Butterflies and hawkmoths",
  BeeSyrphLep = "Focal pollinators",
  VisitedFloral = "Visited plants",
  MeanFloral = "Quadrat flowers",
  MeanFloweringPlant = "Quadrat flowering plants",
  BeeFDis = "Bee functional dispersion (FDis)",
  BeeFEve = "Bee functional evenness (FEve)",
  BeeFDiv = "Bee functional divergence (FDiv)"
)

q1_panel_order <- c(
  "BeeSyrphLep", "Bee", "Bombus", "HB", "NonBombusHB", "Syrphid",
  "ButterflyHawk", "MeanFloral", "MeanFloweringPlant", "BeeFDis",
  "BeeFEve", "BeeFDiv"
)

q1_suffix <- c(
  abundance = "abundance", richness = "richness", shannon = "Shannon diversity", fd = ""
)

# The influence plot fails for honey bee abundance, so skip it there.
q1_checks_without_influence <- c(
  "vif", "qq", "normality", "linearity", "homogeneity",
  "reqq", "pp_check", "binned_residuals", "overdispersion"
)

## Q2 responses
q2_labels <- c(
  zfunctional.complementarity.HL = "Pollinator complementarity",
  zfunctional.complementarity.LL = "Plant complementarity",
  zFunRedundancy.Pols = "Pollinator redundancy",
  zFunRedundancy.Plants = "Plant redundancy",
  zmean.number.of.links.HL = "Pollinator mean number of links",
  zmean.number.of.links.LL = "Plant mean number of links"
)

q2_ys <- c(
  "zfunctional.complementarity.HL", "zfunctional.complementarity.LL",
  "zFunRedundancy.Plants", "zFunRedundancy.Pols",
  "zmean.number.of.links.HL", "zmean.number.of.links.LL"
)

q2_panel_order <- c(
  "zfunctional.complementarity.HL", "zFunRedundancy.Pols",
  "zmean.number.of.links.HL", "zfunctional.complementarity.LL",
  "zFunRedundancy.Plants", "zmean.number.of.links.LL"
)

## Q3 scenarios
q3_labels <- c(
  abundance_low = "Rare plants first",
  abundance_high = "Dominant plants first",
  degree_low = "Specialist plants first",
  degree_high = "Generalist plants first",
  scm_random_plant = "Random plant removal",
  scm_dominant_plant = "Abundance-weighted plant removal"
)

# "BeeAbundance" -> "Bee"; used to look up labels and panel order.
response_key <- function(y) sub("(Abundance|Richness|Diversity)$", "", y)

q1_y_label <- function(y, suffix) {
  key <- response_key(y)
  if (!key %in% names(response_labels)) {
    stop("Missing plot label for ", y, " (key: ", key, ")")
  }
  trimws(paste(response_labels[[key]], suffix))
}


# 1. Read inputs -----------------------------------------------------------

community <- read_keyed_csv(file.path(dir.SI, "data", "sitestats.csv"))

load(file.path(
  dir.SI, "analysis", "network", "saved",
  "corMets_PlantPollinator_YearSR_BeesSyrphidsButterflyHawk_N100.Rdata"
))
network_metrics <- cor.dats

robustness_scores <- read_keyed_csv(
  file.path(dir.EC, "data", "YearSR_robustness_scores.csv")
) %>%
  mutate(Year = factor(Year))

climate <- read_keyed_csv(file.path(dir.EC, "data", "exposure.csv")) %>%
  select(all_of(keys), PC1, PC2)


# 2. Assemble the Q1/Q2 analysis table -------------------------------------

# Keep all community rows.
analysis <- community %>%
  left_join(network_metrics, by = keys, relationship = "one-to-one") %>%
  left_join(climate, by = keys, relationship = "one-to-one") %>%
  mutate(Site = factor(Site), Year = factor(Year))

cat("\n--- Analysis dimensions ---\n")
print(dim(analysis))

cat("\n--- RP 2018 round 1 ---\n")
print(
  analysis %>%
    filter(Site == "RP", Year == "2018", SampleRound == "1") %>%
    select(all_of(keys), PollAbundance, PC1, PC2)
)

stopifnot(
  nrow(analysis) == nrow(community),
  all(is.finite(analysis$PC1)),
  all(is.finite(analysis$PC2)),
  all(c(q2_ys, "PC1", "PC2", "Site", "Year") %in% names(analysis)),
  all(q2_ys %in% names(q2_labels))
)

analysis <- analysis %>%
  select(all_of(unique(c(
    keys, "PC1", "PC2", q1_ys, unname(shannon_richness), q2_ys
  ))))

cat("\n--- Analysis dimensions ---\n")
print(dim(analysis))

write.csv(analysis, file.path(dir.EC, "data", "analysis.csv"), row.names = FALSE)


# 3. Assemble the Q3 analysis table ----------------------------------------

robustness_analysis <- robustness_scores %>%
  left_join(analysis, by = keys, relationship = "many-to-one")

q3_needed <- c(
  "Robustness", "PC1", "PC2", "Site", "Year", "SampleRound",
  "cascade_model", "scenario"
)

stopifnot(
  all(q3_needed %in% names(robustness_analysis)),
  all(robustness_analysis$scenario %in% names(q3_labels))
)

q3_data <- robustness_analysis %>%
  select(all_of(q3_needed)) %>%
  filter(
    is.finite(Robustness), is.finite(PC1), is.finite(PC2),
    !is.na(Site), !is.na(Year), !is.na(SampleRound)
  ) %>%
  mutate(
    Site = factor(Site),
    Year = factor(Year),
    scenario = factor(scenario, levels = names(q3_labels)),
    CommunityID = interaction(Site, Year, SampleRound, drop = TRUE)
  )

# Check how many observations enter the model for each scenario.
print(count(q3_data, cascade_model, scenario))


# 4. Shared plot setup -----------------------------------------------------

geography <- read.csv(
  file.path(dir.SI_saved, "data", "relational", "original", "geography.csv")
)

# Order the site legend from north to south.
site_order <- geography %>%
  filter(!Site %in% c("JM", "CC", "SH")) %>%
  group_by(Site) %>%
  summarize(Lat = mean(Lat, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Lat)) %>%
  pull(Site)

year_levels <- sort(unique(as.character(analysis$Year[!is.na(analysis$Year)])))

lv <- list(site = site_order, year = year_levels)


# 5. Q1: abundance, richness, Shannon diversity ----------------------------

q1_dirs <- make_dirs(dir.EC, "Q1")

## 5a. Fit  ----------------------------
model_groups <- list(
  abundance = fit_models(q1_abundance_ys, analysis, "nb",       q1_dirs$diag, "Q1 abundance"),
  richness  = fit_models(q1_richness_ys,  analysis, "gaussian", q1_dirs$diag, "Q1 richness"),
  shannon   = fit_models(q1_shannon_ys,   analysis, "gaussian", q1_dirs$diag, "Q1 shannon"),
  fd = fit_models(q1_fd_ys, analysis, "gaussian", q1_dirs$diag, "Q1 functional diversity")
)

# Keep these names for code that uses the individual model lists.
abundance_models <- model_groups[["abundance"]]
richness_models  <- model_groups[["richness"]]
shannon_models   <- model_groups[["shannon"]]
fd_models <- model_groups[["fd"]]

# All Q1 models in one named list.
models <- do.call(c, unname(model_groups))

## 5b. Checks and diagnostics  ----------------------------
model_checks <- q1_model_checks(model_groups)
print(model_checks)
write_csv_rows(model_checks, file.path(q1_dirs$diag, "model_checks.csv"))

with_pdf(file.path(q1_dirs$diag, "Q1_diagnostics.pdf"), 11, 8.5, {
  for (y in names(models)) {
    checks <- if (y == "HBAbundance") q1_checks_without_influence else "all"
    diagnostics_page(models[[y]], y, checks)
  }
})

print_singular_fits(model_checks, models)

## 5c. Plot climate effects  ----------------------------
q1_pvalues <- pvalue_matrix(models)
print(q1_pvalues)

ncol_q1 <- 4

for (group in names(model_groups)) {
  group_models <- model_groups[[group]]
  ordered_ys <- names(group_models)[
    order(match(response_key(names(group_models)), q1_panel_order))
  ]
  
  test_line <- if (group == "abundance") {
    "Solid lines indicate slopes with two-sided Wald z-test p < 0.05."
  } else {
    "Solid lines indicate slopes with two-sided Satterthwaite p < 0.05."
  }
  
  write_scatter_pdf(
    file = file.path(q1_dirs$plot, paste0("Q1_", group, "_scatterplots.pdf")),
    models = group_models[ordered_ys],
    pvalues = q1_pvalues,
    lv = lv,
    ncol = ncol_q1,
    width = 14,
    height = max(8, 3.5 * ceiling(length(group_models) / ncol_q1)),
    caption = paste(
      ci_note,
      "The other climate PC is held at its mean among observations used in that model.",
      test_line,
      sep = "\n"
    ),
    y_label = function(y) q1_y_label(y, q1_suffix[[group]])
  )
}


# 6. Q2: network metrics ---------------------------------------------------

q2_dirs <- make_dirs(dir.EC, "Q2")

## 6a. Fit  ----------------------------
q2_models <- fit_models(
  q2_ys, analysis, "gaussian", q2_dirs$diag, "Q2", file_prefix = "Q2_"
)

## 6b. Checks and diagnostics  ----------------------------
q2_checks <- q2_model_checks(q2_models)
print(q2_checks)
write_csv_rows(q2_checks, file.path(q2_dirs$diag, "Q2_model_checks.csv"))

for (y in names(q2_models)) {
  write_diagnostics_pdf(
    q2_models[[y]],
    file.path(q2_dirs$diag, paste0("Q2_", y, "_diagnostics.pdf")),
    label = paste("Q2", y)
  )
}

## 6c. Plot climate effects  ----------------------------
q2_pvalues <- pvalue_matrix(q2_models)
print(q2_pvalues)  # which slopes meet the cutoff

write_scatter_pdf(
  file = file.path(q2_dirs$plot, "Q2_scatterplots.pdf"),
  models = q2_models[q2_panel_order],
  pvalues = q2_pvalues,
  lv = lv,
  ncol = 3,
  width = 14,
  height = 8,
  caption = paste(
    ci_note,
    "The other climate PC is held at its mean among observations used in each model.",
    "Solid fitted lines indicate slopes with two-sided Satterthwaite p < 0.05.",
    "Dashed fitted lines indicate slopes with p >= 0.05. The dotted horizontal line marks the null expectation (z = 0).",
    sep = "\n"
  ),
  y_label = function(y) "Null-standardized z-score",
  title = function(y) q2_labels[[y]],
  null_line = TRUE
)


# 7. Q3: extinction-cascade robustness -------------------------------------

q3_dirs <- make_dirs(dir.EC, "Q3")

## 7a. Fit  ----------------------------
q3_models <- fit_q3_models(q3_data, q3_dirs$diag)

## 7b. Checks and diagnostics  ----------------------------
q3_checks <- q3_model_checks(q3_models, q3_data)
print(q3_checks)
write_csv_rows(q3_checks, file.path(q3_dirs$diag, "Q3_model_checks.csv"))

for (cascade in names(q3_models)) {
  for (scenario_name in names(q3_models[[cascade]])) {
    write_diagnostics_pdf(
      q3_models[[cascade]][[scenario_name]],
      file.path(
        q3_dirs$diag,
        paste0("Q3_", cascade, "_", scenario_name, "_diagnostics.pdf")
      ),
      label = paste("Q3", cascade, scenario_name)
    )
  }
}

## 7c. Plot climate effects  ----------------------------
for (cascade in names(q3_models)) {
  q3_pvalues <- pvalue_matrix(q3_models[[cascade]])
  message("Q3 slope p-values | ", cascade)
  print(q3_pvalues)
  
  write_q3_scatter_pdf(
    file = file.path(q3_dirs$plot, paste0("Q3_", cascade, "_scatterplots.pdf")),
    cascade = cascade,
    models = q3_models[[cascade]],
    pvalues = q3_pvalues,
    lv = lv,
    labels = q3_labels,
    caption = paste(
      ci_note,
      "The other climate PC is held at its mean among observations used for that scenario.",
      "Solid lines indicate unadjusted Satterthwaite p < 0.05. Dashed lines indicate p >= 0.05.",
      sep = "\n"
    )
  )
}


# 8. Master climate slopes table -------------------------------------------

q1_results <- lapply(names(models), function(y) {
  slope_rows(models[[y]], "Q1", y)
})

q2_results <- lapply(names(q2_models), function(y) {
  slope_rows(q2_models[[y]], "Q2", q2_labels[[y]])
})

q3_results <- lapply(names(q3_models), function(cascade) {
  bind_rows(lapply(names(q3_models[[cascade]]), function(scenario_name) {
    slope_rows(
      q3_models[[cascade]][[scenario_name]], "Q3",
      paste(cascade, "robustness"), q3_labels[[scenario_name]]
    )
  }))
})

# Match the unadjusted p < 0.05 cutoff used in the climate plots.
master_results <- bind_rows(q1_results, q2_results, q3_results) %>%
  mutate(
    significance = case_when(
      is.na(p_value) ~ "Unavailable",
      p_value < 0.05 ~ "Significant",
      TRUE ~ "Not significant"
    )
  )

write.csv(
  master_results,
  file.path(dir.EC, "results", "climate_models", "master_climate_results.csv"),
  row.names = FALSE, na = ""
)
