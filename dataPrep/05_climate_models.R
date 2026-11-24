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
  "VisitedFloralRichness",
  "MeanFloralRichness"
)

# Each Shannon response is paired with the richness of the same group.
shannon_richness <- c(
  BeeDiversity = "BeeRichness",
  SyrphidDiversity = "SyrphidRichness",
  BombusDiversity = "BombusRichness",
  ButterflyHawkDiversity = "ButterflyHawkRichness",
  BeeSyrphLepDiversity = "BeeSyrphLepRichness",
  VisitedFloralDiversity = "VisitedFloralRichness",
  MeanFloralDiversity = "MeanFloralRichness"
)

q1_shannon_ys <- c(names(shannon_richness), "MeanFloweringPlantDiversity")

q1_ys <- c(q1_abundance_ys, q1_richness_ys, q1_shannon_ys)

q2_cols <- c(
  "number.of.species.HL",
  "number.of.species.LL",
  "functional.complementarity.HL",
  "functional.complementarity.LL",
  "FunRedundancy.Plants",
  "FunRedundancy.Pols",
  "mean.number.of.links.HL",
  "mean.number.of.links.LL",
  "znumber.of.species.HL",
  "znumber.of.species.LL",
  "zfunctional.complementarity.HL",
  "zfunctional.complementarity.LL",
  "zmean.number.of.links.HL",
  "zmean.number.of.links.LL",
  "zFunRedundancy.Plants",
  "zFunRedundancy.Pols",
  "pnumber.of.species.HL",
  "pnumber.of.species.LL",
  "pfunctional.complementarity.HL",
  "pfunctional.complementarity.LL",
  "pmean.number.of.links.HL",
  "pmean.number.of.links.LL",
  "pFunRedundancy.Plants",
  "pFunRedundancy.Pols"
)

analysis_cols <- unique(c(
  keys, "PC1", "PC2",
  q1_ys, unname(shannon_richness),
  q2_cols
))

analysis <- analysis %>%
  select(all_of(analysis_cols))

cat("\n--- Analysis dimensions ---\n")
print(dim(analysis))

analysis_file <- file.path(model.dir, "Q1_analysis.csv")

write.csv(analysis, analysis_file, row.names = FALSE)

# 4. Assemble Q3 analysis table --------------------------------------

# robustness_analysis <- robustness_scores %>%
#   left_join(
#     analysis,
#     by = keys,
#     relationship = "many-to-one"
#   )

# 5. Fit Q1 response models -----------------------------------------
# PC1 and PC2 are fixed effects; Site and Year each get a random intercept.

# Store fitted models and diagnostics separately from the figures.
model.dir <- file.path(dir.EC, "results", "climate_models", "Q1")
plot.dir <- file.path(model.dir, "plotting")
diag.dir <- file.path(model.dir, "diagnostics")

dir.create(model.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diag.dir, recursive = TRUE, showWarnings = FALSE)

## 5a. fit the abundance models -----------------------------------------
# Fit a negative binomial model for each abundance response.
abundance_models <- list()

for (y in q1_abundance_ys) {
  needed <- c(y, "PC1", "PC2", "Site", "Year")
  dat <- analysis[complete.cases(analysis[, needed]), ]
  message(y, " | n = ", nrow(dat))
  
  equation <- as.formula(
    paste(y, "~ PC1 + PC2 + (1 | Site) + (1 | Year)")
  )
  
  model <- glmer.nb(equation, data = dat)
  print(coef(summary(model)))
  
  abundance_models[[y]] <- model
  saveRDS(model, file.path(diag.dir, paste0(y, "_model.rds")))
}

## 5b. fit the richness models -----------------------------------------
# Fit a Poisson model for each richness response. 
richness_models <- list()

for (y in q1_richness_ys) {
  needed <- c(y, "PC1", "PC2", "Site", "Year")
  dat <- analysis[complete.cases(analysis[, needed]), ]
  message(y, " | n = ", nrow(dat))
  
  equation <- as.formula(
    paste(y, "~ PC1 + PC2 + (1 | Site) + (1 | Year)")
  )
  
  if (y == "BombusRichness") {
    # Year SD was estimated at zero. With Site only, COM-Poisson
    # had lower AIC than Poisson (214.71 vs. 218.36).
    model <- glmmTMB::glmmTMB(
      BombusRichness ~ PC1 + PC2 + (1 | Site),
      data = dat,
      family = glmmTMB::compois()
    )
    print(summary(model)$coefficients$cond)
  } else {
    model <- glmer(equation, data = dat, family = poisson)
    print(coef(summary(model)))
  }
  
  richness_models[[y]] <- model
  saveRDS(model, file.path(diag.dir, paste0(y, "_model.rds")))
}

## 5c. fit the Shannon models -----------------------------------------
# Fit Gaussian mixed models for Shannon diversity where the corresponding richness
# is positive. 
shannon_models <- list()

for (y in q1_shannon_ys) {
  if (y == "MeanFloweringPlantDiversity") {
    presence_col <- "MeanFloweringPlantAbundance"
  } else {
    presence_col <- shannon_richness[[y]]
  }
  
  needed <- c(y, presence_col, "PC1", "PC2", "Site", "Year")
  dat <- analysis[complete.cases(analysis[, needed]), ]
  dat <- dat[dat[[presence_col]] > 0, ]
  
  message(y, " | n = ", nrow(dat))
  
  equation <- as.formula(
    paste(y, "~ PC1 + PC2 + (1 | Site) + (1 | Year)")
  )
  
  model <- lmer(equation, data = dat)
  print(coef(summary(model)))
  
  shannon_models[[y]] <- model
  saveRDS(model, file.path(diag.dir, paste0(y, "_model.rds")))
}

## 5d. Q1 checks and diagnostic PDFs -----------------------------------------

model_groups <- list(
  abundance = abundance_models,
  richness = richness_models,
  shannon = shannon_models
)

models <- c(abundance_models, richness_models, shannon_models)

# Summarize sample size, singularity, and the Pearson dispersion ratio for counts.
model_checks <- data.frame(
  response = names(models),
  group = rep(names(model_groups), lengths(model_groups)),
  n = vapply(models, nobs, numeric(1)),
  singular = vapply(models, performance::check_singularity, logical(1)),
  pearson_dispersion = NA_real_
)

for (y in c(q1_abundance_ys, q1_richness_ys)) {
  model <- models[[y]]
  model_checks$pearson_dispersion[model_checks$response == y] <-
    sum(residuals(model, type = "pearson")^2) / df.residual(model)
}

print(model_checks)
write.csv(
  model_checks,
  file.path(diag.dir, "model_checks.csv"),
  row.names = FALSE
)

# Write one diagnostic PDF per response. Close the PDF even if plotting fails.
for (y in names(models)) {
  pdf(
    file.path(diag.dir, paste0(y, "_diagnostics.pdf")),
    width = 11, height = 8.5
  )
  
  tryCatch({
    plot(performance::check_model(models[[y]], panel = TRUE))
  }, finally = dev.off())
}

# 5e. Assess singularities
singular_ys <- model_checks$response[model_checks$singular]

for (y in singular_ys) {
  cat("\n", y, "\n", sep = "")
  print(lme4::VarCorr(models[[y]]))
}

# 6. Plot Q1 climate effects ------------------------------------------

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

## 6a. coefficient plots ----
# Exponentiate the count-model slopes and intervals to get expected count ratios.
# Shannon estimates are changes on the original scale.
response_labels <- pollinator_labels
response_labels[c("VisitedFloral", "MeanFloral", "MeanFloweringPlant")] <- c(
  "Visited plants",
  "Quadrat flowers",
  "Quadrat flowering plants"
)

axis_labels <- c(
  abundance = "Change in abundance per 1-unit increase in PC score",
  richness  = "Change in richness per 1-unit increase in PC score",
  shannon   = "Change in Shannon diversity per 1-unit increase in PC score"
)

fixed_coefs <- function(model) {
  if (inherits(model, "glmmTMB")) {
    return(summary(model)$coefficients$cond)
  }
  coef(summary(model))
}

for (group in names(model_groups)) {
  group_models <- model_groups[[group]]
  effects <- list()
  
  for (y in names(group_models)) {
    model <- group_models[[y]]
    coefs <- fixed_coefs(model)[c("PC1", "PC2"), , drop = FALSE]
    
    group_key <- sub("(Abundance|Richness|Diversity)$", "", y)
    if (!group_key %in% names(response_labels)) {
      stop("Missing plot label for ", y, " (key: ", group_key, ")")
    }
    
    link <- family(model)$link
    if (!link %in% c("identity", "log")) {
      stop("Unsupported link for ", y, ": ", link)
    }
    
    effects[[y]] <- data.frame(
      taxon = response_labels[[group_key]],
      component = rownames(coefs),
      estimate = coefs[, "Estimate"],
      lower = coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"],
      upper = coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"],
      link = link,
      row.names = NULL
    )
  }
  
  group_effects <- bind_rows(effects)
  
  for (link in unique(group_effects$link)) {
    plot_data <- group_effects[group_effects$link == link, ]
    plot_data$taxon <- factor(
      plot_data$taxon,
      levels = rev(unique(plot_data$taxon))
    )
    
    reference <- 0
    axis_label <- axis_labels[[group]]
    file_group <- group
    
    if (link == "log") {
      plot_data$estimate <- exp(plot_data$estimate)
      plot_data$lower <- exp(plot_data$lower)
      plot_data$upper <- exp(plot_data$upper)
      reference <- 1
      axis_label <- paste0(
        "Expected ", group, " ratio per 1-unit increase in PC score"
      )
    }
    
    if (group == "abundance" && link == "identity") {
      file_group <- "abundance_original_scale"
    }
    
    dodge <- position_dodge(width = 0.5)
    
    p <- ggplot(plot_data, aes(taxon, estimate, color = component)) +
      geom_hline(
        yintercept = reference, linetype = "dashed", color = "grey50"
      ) +
      geom_errorbar(
        aes(ymin = lower, ymax = upper),
        position = dodge, width = 0.15, linewidth = 0.6
      ) +
      geom_point(aes(shape = component), position = dodge, size = 2.8) +
      coord_flip() +
      scale_color_manual(values = c(PC1 = "#0072B2", PC2 = "#B64F00")) +
      scale_shape_manual(values = c(PC1 = 16, PC2 = 17)) +
      labs(
        x = NULL,
        y = axis_label,
        color = NULL,
        shape = NULL,
        caption = paste(
          "Points show fixed-effect estimates.",
          "Bars show approximate 95% Wald intervals."
        )
      ) +
      theme_classic(base_size = 12) +
      theme(
        legend.position = "top",
        plot.caption = element_text(hjust = 0)
      )
    
    print(p)
    ggsave(
      file.path(plot.dir, paste0("Q1_", file_group, "_effects.pdf")),
      plot = p, width = 8, height = 4.8
    )
  }
}

## 6b. scatterplots ----
# Plot observed points and a fixed-effect fitted line for each climate axis.
# The shaded band is an approximate pointwise 95% CI for the fitted mean.

panel_order <- c(
  "BeeSyrphLep", "Bee", "Syrphid", "ButterflyHawk",
  "HB", "Bombus", "NonBombusHB", "VisitedFloral",
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
        
        if (!"Year" %in% names(dat)) {
          # Recover the survey year for the rows used to fit the model.
          survey_rows <- match(rownames(dat), rownames(analysis))
          stopifnot(!anyNA(survey_rows))
          dat$Year <- analysis$Year[survey_rows]
        }
        
        stopifnot(all(as.character(dat$Site) %in% site_order))
        dat$Site <- factor(dat$Site, levels = site_order)
        dat$Year <- factor(dat$Year, levels = year_levels)
        
        group_key <- sub("(Abundance|Richness|Diversity)$", "", y)
        if (!group_key %in% names(response_labels)) {
          stop("Missing plot label for ", y, " (key: ", group_key, ")")
        }
        y_label <- paste(response_labels[[group_key]], response_suffix)
        
        # Vary the plotted PC across the range used in this model.
        scores <- seq(min(dat[[pc]]), max(dat[[pc]]), length.out = 100)
        new_data <- data.frame(
          PC1 = rep(mean(dat$PC1), length(scores)),
          PC2 = rep(mean(dat$PC2), length(scores))
        )
        new_data[[pc]] <- scores
        
        coefs <- fixed_coefs(model)
        fixed_vcov <- if (inherits(model, "glmmTMB")) {
          as.matrix(vcov(model)$cond)
        } else {
          as.matrix(vcov(model))
        }
        
        # Compute the fixed-effect fit and interval on the link scale.
        X <- model.matrix(~ PC1 + PC2, data = new_data)
        terms <- colnames(X)
        beta <- coefs[terms, "Estimate"]
        fixed_vcov <- fixed_vcov[terms, terms, drop = FALSE]
        
        fit_link <- drop(X %*% beta)
        fit_se <- sqrt(rowSums((X %*% fixed_vcov) * X))
        inverse_link <- family(model)$linkinv
        
        estimate <- coefs[pc, "Estimate"]
        cutoff <- 1.96 * coefs[pc, "Std. Error"]
        line_type <- if (abs(estimate) > cutoff) "solid" else "dashed"
        
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
            "Shading shows approximate pointwise 95% confidence intervals",
            "for the fitted means.",
            "Solid lines indicate a PC coefficient whose 95% Wald interval",
            "excludes zero."
          )
        )
      
      print(combined & theme(legend.position = "right"))
    }
  }, finally = dev.off())
}


# # 6. Fit Q2 response models -----------------------------------------
# 
# q2_ys <- c(
#   "functional.complementarity.HL", "functional.complementarity.LL",
#   "FunRedundancy.Pols", "FunRedundancy.Plants",
#   "mean.number.of.links.HL", "mean.number.of.links.LL"
# )

# ys <- q2_ys
#
# model.dir <- file.path(
#   dir.EC, "results", "climate_models", guild, "Q2"
# )
# dir.create(model.dir, recursive = TRUE, showWarnings = FALSE)

# models <- list()
# 
# for (y in ys) {
#   
#   print(y)
#   
#   # Keep complete observations for this particular model.
#   needed <- c(y, "PC1", "PC2", "Site", "Year")
#   dat <- analysis[complete.cases(analysis[, needed]), ]
#   
#   print(nrow(dat))
#   
#   equation <- as.formula(
#     paste(y, "~ PC1 + PC2 + (1 | Site) + (1 | Year)")
#   )
#   
#   model <- lmer(equation, data = dat)
#   
#   print(coef(summary(model)))       # Fixed-effect estimates and standard errors
#   print(as.data.frame(VarCorr(model)))  # Random-effect and residual variances
#   
#   models[[y]] <- model
#   
#   saveRDS(
#     model,
#     file.path(model.dir, paste0(y, "_model.rds"))
#   )
# }

# 7. Check model diagnostics -----------------------------------------

# pdf(file.path(model.dir, "model_diagnostics.pdf"), width = 11, height = 8.5)
# tryCatch({
#   for (y in names(models)) {
#     plot(performance::check_model(models[[y]], panel = TRUE))
#   }
# }, finally = dev.off())
# 
# model_checks <- data.frame(
#   response = names(models),
#   n = vapply(models, nobs, numeric(1)),
#   singular = vapply(models, lme4::isSingular, logical(1))
# )
# 
# print(model_checks)
# write.csv(
#   model_checks,
#   file.path(model.dir, "model_checks.csv"),
#   row.names = FALSE
# )