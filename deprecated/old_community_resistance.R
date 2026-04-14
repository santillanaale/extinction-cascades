# ---- Modeling Community Resistance ----
# Metrics: Network redundancy, complementarity, and generalization
rm(list=ls())
setwd("C:/Users/ale_s/University of Oregon Dropbox/Alejandro Santillana Fernandez/extinction-cascades/analysis/network")
source('src/initialize.R')
source("src/misc.R")
library(dplyr)
library(ggplot2)
library(tidyr)
library(broom.mixed)  # for tidy() and augment() with lme4 models
library(lme4)
load('saved/mods/metrics.Rdata')


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

## ---- By Antecedent Monsoon Precipitation ----

## Join Monsoon Precip into cor.dats
# Shift monsoon year forward so 2011 precip matches 2012 network data
monsoon_precip_data_corrected <- monsoon_precip_data %>%
  mutate(Year = Year + 1)  # Shift year forward by one

# Make sure Site and Year are the same type in both data frames
monsoon_precip_data_corrected$Year <- as.numeric(monsoon_precip_data_corrected$Year)
cor.dats$Year <- as.numeric(cor.dats$Year)

# Merge monsoon precipitation into main dataset
cor.dats <- left_join(cor.dats, monsoon_precip_data_corrected, by = c("Site", "Year"))

# Define the metrics to analyze
ys <- c("FunRedundancy.Pols",
        "FunRedundancy.Plants",
        "functional.complementarity.HL",
        "functional.complementarity.LL",
        "mean.number.of.links.HL",
        "mean.number.of.links.LL")

## Fit models: metric ~ Mean_Monsoon_Precip
mods.monsoon <- lapply(ys, function(y) {
  formula <- as.formula(paste(y, "~ Mean_Monsoon_Precip + (1 | Site)"))
  lmer(formula, data = cor.dats, REML = FALSE)  
})
summary(mods.monsoon[[1]])
names(mods.monsoon) <- ys

library(broom.mixed)

# Extract tidy results for all models
monsoon_results <- lapply(names(mods.monsoon), function(name) {
  broom.mixed::tidy(mods.monsoon[[name]], effects = "fixed") %>%
    mutate(Metric = name)
}) %>%
  bind_rows()

# Keep only the climate predictor rows
monsoon_results_clean <- monsoon_results %>%
  filter(term == "Mean_Monsoon_Precip") %>%
  select(Metric, estimate, std.error, statistic, p.value)

# Print a neat summary table
print(monsoon_results_clean)

## Prediction function using monsoon precipitation
predict_metric_precip <- function(mod, yname, precip_values) {
  new_data <- data.frame(Mean_Monsoon_Precip = precip_values)
  
  mm <- model.matrix(~ Mean_Monsoon_Precip, new_data)
  fit <- mm %*% fixef(mod)
  se <- sqrt(diag(mm %*% vcov(mod) %*% t(mm)))
  
  new_data %>%
    mutate(
      fit = fit,
      se = se,
      lower = fit - 1.96 * se,
      upper = fit + 1.96 * se,
      metric = yname
    )
}

##  Generate predictions across a reasonable precip range 
precip_range <- seq(min(cor.dats$Mean_Monsoon_Precip, na.rm = TRUE),
                    max(cor.dats$Mean_Monsoon_Precip, na.rm = TRUE),
                    length.out = 100)

predictions_precip <- lapply(ys, function(y) {
  predict_metric_precip(mods.monsoon[[y]], y, precip_range)
})
pred_precip_df <- bind_rows(predictions_precip)

## Plot: Metric ~ Monsoon Precip
metric_labels <- c(
  "FunRedundancy.Pols" = "Pollinator Redundancy",
  "FunRedundancy.Plants" = "Plant Redundancy",
  "functional.complementarity.HL" = "Pollinator Complementarity",
  "functional.complementarity.LL" = "Plant Complementarity",
  "mean.number.of.links.HL" = "Pollinator Generalization",
  "mean.number.of.links.LL" = "Plant Generalization"
)

p <- ggplot(pred_precip_df, aes(x = Mean_Monsoon_Precip, y = fit)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80", alpha = 0.4) +
  geom_point(
    data = cor.dats %>%
      pivot_longer(cols = all_of(ys), names_to = "metric", values_to = "value"),
    aes(x = Mean_Monsoon_Precip, y = value, color = Site, shape = as.factor(Year)),
    inherit.aes = FALSE,
    size = 2,
    alpha = 0.6
  ) +
  facet_wrap(~ metric, scales = "free_y", labeller = labeller(metric = metric_labels)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Monsoon Precipitation (mm)",
    y = "Predicted value",
    color = "Site",
    shape = "Year"
  ) +
  theme(strip.text = element_text(size = 12))

ggsave("figures/NetworkMetricsByMonsoonPrecip.pdf", plot = p, width = 10, height = 7)
print(p)

### Filter predictions and data for Plant Generalization
pred_monsoon_LL <- pred_precip_df %>%
  filter(metric == "mean.number.of.links.LL")

data_LL <- cor.dats %>%
  select(Site, Year, Mean_Monsoon_Precip, mean.number.of.links.LL) %>%
  rename(value = mean.number.of.links.LL)

# Extract from the combined results table
coef_LL <- monsoon_results_clean %>%
  filter(Metric == "mean.number.of.links.LL")

slope_LL <- coef_LL$estimate
p_LL <- coef_LL$p.value

label_text <- paste0(
  "Slope = ", round(slope_LL, 3),
  "\nP = ", signif(p_LL, 2)
)

### Plot only Plant Generalization
p_monsoon_LL <- ggplot(pred_monsoon_LL, aes(x = Mean_Monsoon_Precip, y = fit)) +
  geom_line(color = "black", linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "gray80", alpha = 0.4) +
  geom_point(data = data_LL,
             aes(x = Mean_Monsoon_Precip, y = value, color = Site, shape = as.factor(Year)),
             size = 2, alpha = 0.6) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Effect of Monsoon Precipitation on Plant Generalization",
    x = "Monsoon Precipitation (mm)",
    y = "Mean number of links (Plant Generalization)",
    color = "Site",
    shape = "Year"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 12)
  ) +
  annotate("label",
           x = max(pred_monsoon_LL$Mean_Monsoon_Precip) - 
             0.05 * diff(range(pred_monsoon_LL$Mean_Monsoon_Precip)),
           y = max(pred_monsoon_LL$fit) - 
             0.05 * diff(range(pred_monsoon_LL$fit)),
           label = label_text,
           hjust = 1, vjust = 1,
           size = 5,
           fill = "white",
           color = "black",
           alpha = 0.8)

### Print and save plot
print(p_monsoon_LL)

ggsave("figures/PlantGeneralizationByMonsoonPrecip.pdf", plot = p_monsoon_LL, width = 10, height = 7)

## ---- By Winter Precipitation ----

winter_precip_data <- winter_precip_data %>%
  rename(Year = Winter_Year)

# Ensure types match for join
winter_precip_data$Year <- as.numeric(winter_precip_data$Year)
cor.dats$Year <- as.numeric(cor.dats$Year)

# Join winter precip data to main dataset
cor.dats <- left_join(cor.dats, winter_precip_data, by = c("Site", "Year"))

# Define metrics to model
ys <- c("FunRedundancy.Pols",
        "FunRedundancy.Plants",
        "functional.complementarity.HL",
        "functional.complementarity.LL",
        "mean.number.of.links.HL",
        "mean.number.of.links.LL")

### Fit models: metric ~ Mean_Winter_Precip + (1 | Site)
mods.winter <- lapply(ys, function(y) {
  formula <- as.formula(paste(y, "~ Mean_Winter_Precip + (1 | Site)"))
  lmer(formula, data = cor.dats, REML = FALSE)
})
names(mods.winter) <- ys

# Extract tidy results for all models
winter_results <- lapply(names(mods.winter), function(name) {
  broom.mixed::tidy(mods.winter[[name]], effects = "fixed") %>%
    mutate(Metric = name)
}) %>%
  bind_rows()

# Keep only the climate predictor rows
winter_results_clean <- winter_results %>%
  filter(term == "Mean_Winter_Precip") %>%
  select(Metric, estimate, std.error, statistic, p.value)

# Print a neat summary table
print(winter_results_clean)

### Prediction function for winter precip
predict_metric_winter <- function(mod, yname, precip_values) {
  new_data <- data.frame(Mean_Winter_Precip = precip_values)
  
  mm <- model.matrix(~ Mean_Winter_Precip, new_data)
  fit <- mm %*% fixef(mod)
  se <- sqrt(diag(mm %*% vcov(mod) %*% t(mm)))
  
  new_data %>%
    mutate(
      fit = fit,
      se = se,
      lower = fit - 1.96 * se,
      upper = fit + 1.96 * se,
      metric = yname
    )
}

### Generate predictions
precip_range_winter <- seq(min(cor.dats$Mean_Winter_Precip, na.rm = TRUE),
                           max(cor.dats$Mean_Winter_Precip, na.rm = TRUE),
                           length.out = 100)

predictions_winter <- lapply(ys, function(y) {
  predict_metric_winter(mods.winter[[y]], y, precip_range_winter)
})
pred_winter_df <- bind_rows(predictions_winter)

### Plot
metric_labels <- c(
  "FunRedundancy.Pols" = "Pollinator Redundancy",
  "FunRedundancy.Plants" = "Plant Redundancy",
  "functional.complementarity.HL" = "Pollinator Complementarity",
  "functional.complementarity.LL" = "Plant Complementarity",
  "mean.number.of.links.HL" = "Pollinator Generalization",
  "mean.number.of.links.LL" = "Plant Generalization"
)

p_winter <- ggplot(pred_winter_df, aes(x = Mean_Winter_Precip, y = fit)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80", alpha = 0.4) +
  geom_point(
    data = cor.dats %>%
      pivot_longer(cols = all_of(ys), names_to = "metric", values_to = "value"),
    aes(x = Mean_Winter_Precip, y = value, color = Site, shape = as.factor(Year)),
    inherit.aes = FALSE,
    size = 2,
    alpha = 0.6
  ) +
  facet_wrap(~ metric, scales = "free_y", labeller = labeller(metric = metric_labels)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Winter Precipitation (mm)",
    y = "Predicted value",
    color = "Site",
    shape = "Year"
  ) +
  theme(strip.text = element_text(size = 12))

ggsave("figures/NetworkMetricsByWinterPrecip.pdf", plot = p_winter, width = 10, height = 7)
print(p_winter)



### Filter predictions and data for Plant Generalization
pred_winter_LL <- pred_winter_df %>%
  filter(metric == "mean.number.of.links.LL")

data_LL <- cor.dats %>%
  select(Site, Year, Mean_Winter_Precip, mean.number.of.links.LL) %>%
  rename(value = mean.number.of.links.LL)

# Extract from the combined results table
coef_LL <- winter_results_clean %>%
  filter(Metric == "mean.number.of.links.LL")

slope_LL <- coef_LL$estimate
p_LL <- coef_LL$p.value

label_text <- paste0(
  "Slope = ", round(slope_LL, 3),
  "\nP = ", signif(p_LL, 2)
)

### Plot only Plant Generalization
p_winter_LL <- ggplot(pred_winter_LL, aes(x = Mean_Winter_Precip, y = fit)) +
  geom_line(color = "black", linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "gray80", alpha = 0.4) +
  geom_point(data = data_LL,
             aes(x = Mean_Winter_Precip, y = value, color = Site, shape = as.factor(Year)),
             size = 2, alpha = 0.6) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Effect of Winter Precipitation on Plant Generalization",
    x = "Winter Precipitation (mm)",
    y = "Mean number of links (Plant Generalization)",
    color = "Site",
    shape = "Year"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 12)
  )

p_winter_LL <- p_winter_LL +
  annotate("label",
           x = min(pred_winter_LL$Mean_Winter_Precip) + 
             0.05 * diff(range(pred_winter_LL$Mean_Winter_Precip)),
           y = max(pred_winter_LL$fit) - 
             0.05 * diff(range(pred_winter_LL$fit)),
           label = label_text,
           hjust = 0, vjust = 1,
           size = 5,
           fill = "white",
           color = "black",
           alpha = 0.8)


### Print plot
p_winter_LL
ggsave("figures/PlantGeneralizationByWinterPrecip.pdf", plot = p_winter_LL, width = 10, height = 7)

## ---- By Degree Days ----
annual_gdd <- annual_gdd %>% rename(Site = site_name)
annual_gdd <- annual_gdd %>% rename(Year = year)

# Merge with your cor.dats
cor.dats <- left_join(cor.dats, annual_gdd, by = c("Site", "Year"))

ys <- c("FunRedundancy.Pols",
        "FunRedundancy.Plants",
        "functional.complementarity.HL",
        "functional.complementarity.LL",
        "mean.number.of.links.HL",
        "mean.number.of.links.LL")

# Model
mods.annual_gdd <- lapply(ys, function(y) {
  formula <- as.formula(paste(y, "~ total_gdd + (1 | Site)"))
  lmer(formula, data = cor.dats, REML = FALSE)
})
names(mods.annual_gdd) <- ys

# Extract tidy results for all models
gdd_results <- lapply(names(mods.annual_gdd), function(name) {
  broom.mixed::tidy(mods.annual_gdd[[name]], effects = "fixed") %>%
    mutate(Metric = name)
}) %>%
  bind_rows()

# Keep only the climate predictor rows
gdd_results_clean <- gdd_results %>%
  filter(term == "mean_annual_gdd") %>%
  select(Metric, estimate, std.error, statistic, p.value)

# Print a neat summary table
print(gdd_results_clean)

# Prediction function
predict_metric_gdd <- function(mod, yname, gdd_values) {
  new_data <- data.frame(total_gdd = gdd_values)
  
  # Construct model matrix for fixed effects only
  mm <- model.matrix(~ total_gdd, new_data)
  
  # Get fixed effect estimates
  fit <- mm %*% fixef(mod)
  
  # Calculate standard errors
  se <- sqrt(diag(mm %*% vcov(mod) %*% t(mm)))
  
  # Return a dataframe with fit and confidence intervals
  new_data %>%
    mutate(
      fit = as.vector(fit),
      se = se,
      lower = fit - 1.96 * se,
      upper = fit + 1.96 * se,
      metric = yname
    )
}

# Generate prediction data over a reasonable GDD range
gdd_range <- seq(min(cor.dats$total_gdd, na.rm = TRUE),
                 max(cor.dats$total_gdd, na.rm = TRUE),
                 length.out = 100)

predictions_gdd <- lapply(ys, function(y) {
  predict_metric_gdd(mods.annual_gdd[[y]], y, gdd_range)
})

pred_gdd_df <- bind_rows(predictions_gdd)

# Plotting
metric_labels <- c(
  "FunRedundancy.Pols" = "Pollinator Redundancy",
  "FunRedundancy.Plants" = "Plant Redundancy",
  "functional.complementarity.HL" = "Pollinator Complementarity",
  "functional.complementarity.LL" = "Plant Complementarity",
  "mean.number.of.links.HL" = "Pollinator Generalization",
  "mean.number.of.links.LL" = "Plant Generalization"
)

p_gdd <- ggplot(pred_gdd_df, aes(x = total_gdd, y = fit)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80", alpha = 0.4) +
  geom_point(
    data = cor.dats %>%
      pivot_longer(cols = all_of(ys), names_to = "metric", values_to = "value"),
    aes(x = total_gdd, y = value, color = Site, shape = as.factor(Year)),
    inherit.aes = FALSE,
    size = 2,
    alpha = 0.6
  ) +
  facet_wrap(~ metric, scales = "free_y", labeller = labeller(metric = metric_labels)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Total Growing Degree Days",
    y = "Predicted Network Metric",
    color = "Site",
    shape = "Year"
  ) +
  theme(strip.text = element_text(size = 12))

# Save and print plot
ggsave("figures/NetworkMetricsByTotalGDD.pdf", plot = p_gdd, width = 10, height = 7)
print(p_gdd)