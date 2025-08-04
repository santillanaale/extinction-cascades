rm(list=ls())
setwd("C:/Users/ale_s/University of Oregon Dropbox/Alejandro Santillana Fernandez/extinction-cascades/analysis/network")
source('src/initialize.R')
source("src/misc.R")
# source("plotting/src/predictIntervals.R")
# source("plotting/src/CIplotting.R")
# source("plotting/src/diagnostics.R")
# source("plotting/src/plotNetworkMets.R")
library(dplyr)
library(ggplot2)
library(tidyr)
library(broom.mixed)  # for tidy() and augment() with lme4 models
load('saved/mods/metrics.Rdata')

## ************************************************************
## network metrics by year
## ************************************************************
ys <- c("FunRedundancy.Pols",
        "FunRedundancy.Plants",
        "functional.complementarity.HL",
        "functional.complementarity.LL",
        "mean.number.of.links.HL",
        "mean.number.of.links.LL")

# Ensure Year is numeric for prediction grid
cor.dats$Year <- as.numeric(as.character(cor.dats$Year))

# Create prediction data frame
predict_metric <- function(mod, yname, years) {
  new_data <- data.frame(Year = years)
  
  # Design matrix
  mm <- model.matrix(~ Year, new_data)
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

# Get unique years from data
unique_years <- sort(unique(cor.dats$Year))

# Run prediction for each model
predictions <- lapply(ys, function(y) {
  predict_metric(mods.year[[y]], y, unique_years)
})

pred_df <- bind_rows(predictions)


# Define your labels
metric_labels <- c(
  "FunRedundancy.Pols" = "Pollinator Redundancy",
  "FunRedundancy.Plants" = "Plant Redundancy",
  "functional.complementarity.HL" = "Pollinator Complementarity",
  "functional.complementarity.LL" = "Plant Complementarity",
  "mean.number.of.links.HL" = "Pollinator Generalization",
  "mean.number.of.links.LL" = "Plant Generalization"
)

# Build the plot
p <- ggplot(pred_df, aes(x = Year, y = fit)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80", alpha = 0.4) +
  geom_point(data = cor.dats %>%
               pivot_longer(cols = all_of(ys), names_to = "metric", values_to = "value"),
             aes(x = Year, y = value),
             inherit.aes = FALSE,
             alpha = 0.4, shape = 1) +
  facet_wrap(~ metric, scales = "free_y", labeller = labeller(metric = metric_labels)) +
  theme_minimal(base_size = 14) +
  labs(x = "Year", y = "Predicted value") +
  theme(strip.text = element_text(size = 12)) +
  scale_x_continuous(
    breaks = seq(floor(min(pred_df$Year)), ceiling(max(pred_df$Year)), by = 3),
    labels = function(x) as.integer(x)
  )

# Save the plot as PDF
ggsave("figures/NetworkMetricsByYear.pdf", plot = p, width = 10, height = 7)

# Optional: display the plot in R console
print(p)






# ys <- names(mods.year)
# ylabs <- c("Pollinator \n redundancy",
#            "Plant \n redundancy",
#            "Pollinator \n complementarity",
#            "Plant \n complementarity",
#            "Pollinator \n generalization",
#            "Plant \n generalization")
# names(ylabs) <- ys

# makePanels(mods = mods.year,
#            xvars = "Year",
#            xlabel = "Year",
#            ys = ys)
# 
# 
# ## ************************************************************
# ## network metrics by Richness
# ## ************************************************************
# 
# 
# ys <- names(mods.rich)
# ylabs <- c("Pollinator \n redundancy",
#            "Plant \n redundancy",
#            "Pollinator \n complementarity",
#            "Plant \n complementarity",
#            "Plant \n generalization",
#            "Pollinator \n generalization")
# 
# names(ylabs) <- ys
# 
# makePanels(mods=mods.rich,
#            xvars= c("Richness", "FloralRichness"),
#            xlabel= "Richness",
#            ys=ys)