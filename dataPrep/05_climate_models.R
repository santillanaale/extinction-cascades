rm(list=ls())
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

keys <- c("Site", "Year", "SampleRound")

# 1. Read community summaries ----------------------------------------

community <- read.csv(
  file.path(dir.SI, "data", "sitestats.csv"),
  colClasses = c(SampleRound = "character")
)

# 2. Read Q2, Q3, and climate outputs --------------------------------

load(file.path(
  dir.SI, "analysis", "network", "saved",
  "corMets_PlantPollinator_YearSR_BeesSyrphidsButterflyHawk_N100.Rdata" 
))

network_metrics <- cor.dats

# robustness_scores <- read.csv(
#   "PATH_TO_ROBUSTNESS_SCORES.csv",
#   colClasses = c(SampleRound = "character")
# )

climate <- read.csv(
  file.path(dir.EC, "data", "exposure.csv"),
  colClasses = c(SampleRound = "character")
) %>%
  select(all_of(keys), PC1, PC2)

# 3. Assemble Q1/Q2 analysis table -----------------------------------

# Keep all community rows.
analysis <- community %>%
  left_join(
    network_metrics,
    by = keys,
    relationship = "one-to-one"
  ) %>%
  left_join(
    climate,
    by = keys,
    relationship = "one-to-one"
  ) %>%
  mutate(
    Site = factor(Site),
    Year = factor(Year)
  )

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
  all(is.finite(analysis$PC2))
)

# Select the model inputs
q1_abundance_ys <- c(
  "BeeAbundance",
  "SyrphidAbundance",
  "HBAbundance",
  "BombusAbundance",
  "NonBombusHBAbundance",
  "ButterflyHawkAbundance",
  "BeeSyrphLepAbundance",
  "MeanFloralAbundance",
  "MeanFloweringPlantAbundance"
)

q1_richness_ys <- c(
  "BeeRichness",
  "SyrphidRichness",
  "BombusRichness",
  "ButterflyHawkRichness",
  "BeeSyrphLepRichness",
  "MeanFloralRichness"
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

q1_ys <- c(q1_abundance_ys, q1_richness_ys, q1_shannon_ys)

q2_ys <- c(
  "zfunctional.complementarity.HL",
  "zfunctional.complementarity.LL",
  "zFunRedundancy.Plants",
  "zFunRedundancy.Pols",
  "zmean.number.of.links.HL",
  "zmean.number.of.links.LL"
)

analysis_cols <- unique(c(
  keys, "PC1", "PC2",
  q1_ys, unname(shannon_richness),
  q2_ys
))

analysis <- analysis %>%
  select(all_of(analysis_cols))

cat("\n--- Analysis dimensions ---\n")
print(dim(analysis))

analysis_file <- file.path(dir.EC, "data", "analysis.csv")

write.csv(analysis, analysis_file, row.names = FALSE)

# 4. Assemble Q3 analysis table --------------------------------------

# robustness_analysis <- robustness_scores %>%
#   left_join(
#     analysis,
#     by = keys,
#     relationship = "many-to-one"
#   )

# 5. Fit and Plot Q1 Models  --------------------------------------------------

# Store fitted models and diagnostics separately from the figures.
model.dir <- file.path(dir.EC, "results", "climate_models", "Q1")
plot.dir <- file.path(model.dir, "plotting")
diag.dir <- file.path(model.dir, "diagnostics")

dir.create(model.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diag.dir, recursive = TRUE, showWarnings = FALSE)

## 5a. Fit the Q1 models --------------------------------------------------
responses <- stack(list(
  abundance = q1_abundance_ys,
  richness = q1_richness_ys,
  shannon = q1_shannon_ys
))

model_groups <- list(
  abundance = list(),
  richness = list(),
  shannon = list()
)

for (i in seq_len(nrow(responses))) {
  y <- responses$values[i]
  group <- as.character(responses$ind[i])
  needed <- c(y, "PC1", "PC2", "Site", "Year")
  dat <- analysis[complete.cases(analysis[, needed]), ]
  
  message(group, " | ", y, " | n = ", nrow(dat))
  equation <- as.formula(
    paste(y, "~ PC1 + PC2 + (1 | Site) + (1 | Year)")
  )
  
  if (group == "abundance") {
    model <- lme4::glmer.nb(equation, data = dat)
  } else {
    model <- lme4::lmer(equation, data = dat)
  }
  
  print(coef(summary(model)))
  model_groups[[group]][[y]] <- model
  saveRDS(model, file.path(diag.dir, paste0(y, "_model.rds")))
}

# Keep these names for code that uses the individual model lists.
abundance_models <- model_groups[["abundance"]]
richness_models <- model_groups[["richness"]]
shannon_models <- model_groups[["shannon"]]

## 5b. Checks and diagnostic PDFs -----------------------------------------
# Combine the groups into one list, keeping the response names.
models <- do.call(c, unname(model_groups))

# Summarize sample size, singularity, and the Pearson dispersion ratio for counts.
model_checks <- data.frame(
  response = names(models),
  group = rep(names(model_groups), lengths(model_groups)),
  n = vapply(models, nobs, numeric(1)),
  singular = vapply(models, performance::check_singularity, logical(1)),
  pearson_dispersion = NA_real_
)

for (i in which(model_checks$group == "abundance")) {
  model <- models[[i]]
  model_checks$pearson_dispersion[i] <-
    sum(residuals(model, type = "pearson")^2) / df.residual(model)
}

print(model_checks)
write.csv(
  model_checks,
  file.path(diag.dir, "model_checks.csv"),
  row.names = FALSE
)

# Write one diagnostic PDF per response. Close it even if plotting fails.
pdf(file.path(diag.dir, "Q1_diagnostics.pdf"), width = 11, height = 8.5)

tryCatch({
  for (y in names(models)) {
    message("Writing diagnostics for ", y)
    
    checks <- "all"
    if (y == "HBAbundance") {
      # Skip the influence plot that failed for honey bee abundance.
      checks <- c(
        "vif", "qq", "normality", "linearity", "homogeneity",
        "reqq", "pp_check", "binned_residuals", "overdispersion"
      )
    }
    
    plot(performance::check_model(models[[y]], panel = TRUE, check = checks))
  }
}, finally = dev.off())

# Assess singularities
singular_ys <- model_checks$response[model_checks$singular]

for (y in singular_ys) {
  cat("\n", y, "\n", sep = "")
  print(lme4::VarCorr(models[[y]]))
}

## 5c. Plot Q1 climate effects ------------------------------------------

geography <- read.csv(
  file.path(dir.SI_saved, "data", "relational", "original", "geography.csv"))

# Order the site legend from north to south.
site_order <- geography %>%
  filter(!Site %in% c("JM", "CC", "SH")) %>%
  group_by(Site) %>%
  summarize(Lat = mean(Lat, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Lat)) %>%
  pull(Site)

pollinator_labels <- c(
  Bee = "Bees",
  Syrphid = "Hoverflies",
  HB = "Honey bees",
  Bombus = "Bumble bees",
  NonBombusHB = "Other bees",
  ButterflyHawk = "Butterflies and hawkmoths",
  BeeSyrphLep = "Focal pollinators"
)

response_labels <- c(
  pollinator_labels,
  VisitedFloral = "Visited plants",
  MeanFloral = "Quadrat flowers",
  MeanFloweringPlant = "Quadrat flowering plants"
)

# Plot the observations used in each model and the fixed-effect fitted curve.
# The shaded band is an approximate pointwise 95% CI for the fitted mean.

panel_order <- c(
  "BeeSyrphLep", "Bee", "Bombus", "HB", "NonBombusHB", "Syrphid", "ButterflyHawk",
  "MeanFloral", "MeanFloweringPlant"
)

year_levels <- sort(unique(as.character(analysis$Year[!is.na(analysis$Year)])))
ncol_panels <- 4

for (group in names(model_groups)) {
  group_models <- model_groups[[group]]
  if (length(group_models) == 0L) next
  
  response_suffix <- switch(
    group,
    abundance = "abundance",
    richness = "richness",
    shannon = "Shannon diversity"
  )
  
  nrow_panels <- ceiling(length(group_models) / ncol_panels)
  
  pdf(
    file.path(plot.dir, paste0("Q1_", group, "_scatterplots.pdf")),
    width = 14, height = max(8, 3.5 * nrow_panels)
  )
  
  tryCatch({
    for (pc in c("PC1", "PC2")) {
      response_plots <- list()
      
      for (y in names(group_models)) {
        model <- group_models[[y]]
        dat <- model.frame(model)
        
        stopifnot(all(as.character(dat$Site) %in% site_order))
        dat$Site <- factor(dat$Site, levels = site_order)
        dat$Year <- factor(dat$Year, levels = year_levels)
        
        group_key <- sub("(Abundance|Richness|Diversity)$", "", y)
        if (!group_key %in% names(response_labels)) {
          stop("Missing plot label for ", y, " (key: ", group_key, ")")
        }
        y_label <- paste(response_labels[[group_key]], response_suffix)
        
        # Vary the plotted PC and hold the other PC at its model-data mean.
        scores <- seq(min(dat[[pc]]), max(dat[[pc]]), length.out = 100)
        new_data <- data.frame(
          PC1 = rep(mean(dat$PC1), length(scores)),
          PC2 = rep(mean(dat$PC2), length(scores))
        )
        new_data[[pc]] <- scores
        
        # Compute the fixed-effect fit and interval on the link scale.
        X <- model.matrix(~ PC1 + PC2, data = new_data)
        coef_names <- colnames(X)
        beta <- lme4::fixef(model)[coef_names]
        fixed_vcov <- as.matrix(vcov(model))[coef_names, coef_names, drop = FALSE]
        
        fit_link <- drop(X %*% beta)
        fit_se <- sqrt(rowSums((X %*% fixed_vcov) * X))
        inverse_link <- family(model)$linkinv
        
        pc_se <- sqrt(fixed_vcov[pc, pc])
        line_type <- if (abs(beta[[pc]]) > 1.96 * pc_se) "solid" else "dashed"
        
        fit_data <- data.frame(
          score = scores,
          fitted = inverse_link(fit_link),
          lower = inverse_link(fit_link - 1.96 * fit_se),
          upper = inverse_link(fit_link + 1.96 * fit_se)
        )
        
        response_plots[[y]] <- ggplot(dat) +
          geom_ribbon(
            data = fit_data,
            aes(x = score, ymin = lower, ymax = upper),
            fill = "gray80", alpha = 0.4
          ) +
          geom_point(
            aes(x = .data[[pc]], y = .data[[y]], color = Site, shape = Year),
            alpha = 0.7, size = 2,
            show.legend = TRUE
          ) +
          geom_line(
            data = fit_data,
            aes(x = score, y = fitted),
            color = "black", linewidth = 1, linetype = line_type
          ) +
          scale_color_discrete(limits = site_order, drop = FALSE) +
          scale_shape_discrete(limits = year_levels, drop = FALSE) +
          labs(x = pc, y = y_label, color = "Site (N to S)", shape = "Year") +
          theme_minimal(base_size = 12) +
          theme(
            panel.grid.minor = element_blank(),
            plot.margin = margin(8, 12, 8, 8)
          )
      }
      
      panel_groups <- sub(
        "(Abundance|Richness|Diversity)$", "", names(response_plots)
      )
      response_plots <- response_plots[order(match(panel_groups, panel_order))]
      
      combined <- patchwork::wrap_plots(
        response_plots, ncol = ncol_panels, guides = "collect"
      ) +
        patchwork::plot_annotation(
          caption = paste(
            "Shading shows approximate pointwise 95% confidence intervals for the fitted means.",
            "The other climate PC is held at its mean among observations used in that model.",
            "Solid lines indicate a PC coefficient whose 95% Wald interval excludes zero.",
            sep = "\n"
          )
        )
      
      print(combined & theme(legend.position = "right"))
    }
  }, finally = dev.off())
}

# 6. Fit and Plot Q2 Models --------------------------------------------------

model.dir <- file.path(dir.EC, "results", "climate_models", "Q2")
plot.dir <- file.path(model.dir, "plotting")
diag.dir <- file.path(model.dir, "diagnostics")

dir.create(model.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diag.dir, recursive = TRUE, showWarnings = FALSE)

## 6a. Fit the Q2 models ----------------------------------------------------

q2_labels <- c(
  zfunctional.complementarity.HL = "Pollinator complementarity",
  zfunctional.complementarity.LL = "Plant complementarity",
  zFunRedundancy.Pols = "Pollinator redundancy",
  zFunRedundancy.Plants = "Plant redundancy",
  zmean.number.of.links.HL = "Pollinator mean number of links",
  zmean.number.of.links.LL = "Plant mean number of links"
)

stopifnot(
  all(c(q2_ys, "PC1", "PC2", "Site", "Year") %in% names(analysis)),
  all(q2_ys %in% names(q2_labels))
)

q2_models <- list()

for (y in q2_ys) {
  needed <- c(y, "PC1", "PC2", "Site", "Year")
  keep <- complete.cases(analysis[, needed]) & is.finite(analysis[[y]])
  dat <- analysis[keep, needed, drop = FALSE]
  
  message("Q2 | ", y, " | n = ", nrow(dat))
  
  equation <- as.formula(
    paste(y, "~ PC1 + PC2 + (1 | Site) + (1 | Year)")
  )
  
  model <- lme4::lmer(equation, data = dat)
  
  print(coef(summary(model)))
  q2_models[[y]] <- model
  saveRDS(model, file.path(diag.dir, paste0("Q2_", y, "_model.rds")))
}

## 6b. Checks and diagnostic PDFs -------------------------------------------

q2_checks <- data.frame(
  response = names(q2_models),
  n = vapply(q2_models, nobs, numeric(1)),
  singular = vapply(q2_models, performance::check_singularity, logical(1)),
  site_sd = vapply(q2_models, function(model) {
    as.numeric(attr(lme4::VarCorr(model)[["Site"]], "stddev"))
  }, numeric(1)),
  year_sd = vapply(q2_models, function(model) {
    as.numeric(attr(lme4::VarCorr(model)[["Year"]], "stddev"))
  }, numeric(1))
)

print(q2_checks)
write.csv(
  q2_checks,
  file.path(diag.dir, "Q2_model_checks.csv"),
  row.names = FALSE
)

for (y in names(q2_models)) {
  message("Writing Q2 diagnostics for ", y)
  
  pdf(
    file.path(diag.dir, paste0("Q2_", y, "_diagnostics.pdf")),
    width = 11, height = 8.5
  )
  
  tryCatch({
    diagnostics <- performance::check_model(q2_models[[y]], panel = TRUE)
    print(plot(diagnostics))
  }, error = function(err) {
    warning("Q2 diagnostics failed for ", y, ": ", conditionMessage(err))
  }, finally = dev.off())
}

## 6c. Plot Q2 climate effects ----------------------------------------------

q2_panel_order <- c(
  "zfunctional.complementarity.HL",
  "zFunRedundancy.Pols",
  "zmean.number.of.links.HL",
  "zfunctional.complementarity.LL",
  "zFunRedundancy.Plants",
  "zmean.number.of.links.LL"
)

pdf(
  file.path(plot.dir, "Q2_scatterplots.pdf"),
  width = 14, height = 8
)

tryCatch({
  for (pc in c("PC1", "PC2")) {
    response_plots <- list()
    
    for (y in q2_panel_order)  {
      model <- q2_models[[y]]
      dat <- model.frame(model)
      
      stopifnot(all(as.character(dat$Site) %in% site_order))
      dat$Site <- factor(dat$Site, levels = site_order)
      dat$Year <- factor(dat$Year, levels = year_levels)
      
      # Vary the plotted PC and hold the other PC at its model-data mean.
      scores <- seq(min(dat[[pc]]), max(dat[[pc]]), length.out = 100)
      new_data <- data.frame(
        PC1 = rep(mean(dat$PC1), length(scores)),
        PC2 = rep(mean(dat$PC2), length(scores))
      )
      new_data[[pc]] <- scores
      
      X <- model.matrix(~ PC1 + PC2, data = new_data)
      coef_names <- colnames(X)
      beta <- lme4::fixef(model)[coef_names]
      fixed_vcov <- as.matrix(vcov(model))[coef_names, coef_names, drop = FALSE]
      
      fitted <- drop(X %*% beta)
      fit_se <- sqrt(rowSums((X %*% fixed_vcov) * X))
      pc_se <- sqrt(fixed_vcov[pc, pc])
      line_type <- if (abs(beta[[pc]]) > 1.96 * pc_se) "solid" else "dashed"
      
      fit_data <- data.frame(
        score = scores,
        fitted = fitted,
        lower = fitted - 1.96 * fit_se,
        upper = fitted + 1.96 * fit_se
      )
      
      response_plots[[y]] <- ggplot(dat) +
        geom_hline(yintercept = 0, linetype = "dotted", color = "gray55") +
        geom_ribbon(
          data = fit_data,
          aes(x = score, ymin = lower, ymax = upper),
          fill = "gray80", alpha = 0.4
        ) +
        geom_point(
          aes(x = .data[[pc]], y = .data[[y]], color = Site, shape = Year),
          alpha = 0.7, size = 2, show.legend = TRUE
        ) +
        geom_line(
          data = fit_data,
          aes(x = score, y = fitted),
          color = "black", linewidth = 1, linetype = line_type
        ) +
        scale_color_discrete(limits = site_order, drop = FALSE) +
        scale_shape_discrete(limits = year_levels, drop = FALSE) +
        labs(
          x = pc,
          y = "Null-standardized z-score",
          title = q2_labels[[y]],
          color = "Site (N to S)",
          shape = "Year"
        ) +
        theme_minimal(base_size = 12) +
        theme(panel.grid.minor = element_blank())
    }
    
    combined <- patchwork::wrap_plots(
      response_plots, ncol = 3, guides = "collect"
    ) +
      patchwork::plot_annotation(
        caption = paste(
          "Shading shows approximate pointwise 95% confidence intervals for the fitted means.",
          "The other climate PC is held at its mean among observations used in that model.",
          "Solid fitted lines indicate a PC coefficient whose 95% Wald interval excludes zero.",
          sep = "\n"
        )
      )
    
    print(combined & theme(legend.position = "right"))
  }
}, finally = dev.off())