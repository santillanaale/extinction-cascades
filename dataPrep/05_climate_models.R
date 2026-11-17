rm(list=ls())
setwd("~/")
source("lab_paths.R")

library(dplyr)
library(lme4)

dir.SI <- file.path(local.path, "skyIslands")
dir.EC <- file.path(local.path, "extinction-cascades")

keys <- c("Site", "Year", "SampleRound")

# 1. Read community summaries ----------------------------------------

community <- read.csv(
  file.path(dir.SI, "data", "sitestats.csv"),
  colClasses = c(SampleRound = "character")
)

# 2. Read Q2, Q3, and climate outputs --------------------------------
# Replace these paths after confirming the existing exports.

# Temporary Q2 input: replace after tracing the current calculations.
network_metrics <- readRDS(
  file.path(
    dir.EC,
    "analysis", "network", "saved",
    "q2_metrics_table_expanded_61_N100.rds"
  )
)

# robustness_scores <- read.csv(
#   "PATH_TO_ROBUSTNESS_SCORES.csv",
#   colClasses = c(SampleRound = "character")
# )

# PCA input choice pending confirmation from Lauren.
climate <- read.csv(
  file.path(dir.EC, "data", "exposure.csv"),
  colClasses = c(SampleRound = "character")
) %>%
  select(all_of(keys), PC1, PC2)

# 3. Prepare and check event keys ------------------------------------

# Recover event keys if the metrics table contains only Network.
if (!all(keys %in% names(network_metrics))) {
  
  network_metrics <- tidyr::separate(
    network_metrics,
    Network,
    into = keys,
    sep = "\\.",
    remove = FALSE
  )
}

network_metrics <- network_metrics %>%
  mutate(
    Site = as.character(Site),
    Year = as.numeric(Year),
    SampleRound = as.character(SampleRound)
  )

stopifnot(
  !anyDuplicated(community[keys]),
  !anyDuplicated(network_metrics[keys]),
  !anyDuplicated(climate[keys])
)

# 4. Assemble Q1/Q2 analysis table -----------------------------------

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

cat("\n--- Missing PCA scores ---\n")
print(colSums(is.na(analysis[c("PC1", "PC2")])))

cat("\n--- RP 2018 round 1 ---\n")
print(
  analysis %>%
    filter(Site == "RP", Year == "2018", SampleRound == "1") %>%
    select(all_of(keys), PollAbundance, Network, PC1, PC2)
)

stopifnot(
  nrow(analysis) == nrow(community),
  all(is.finite(analysis$PC1)),
  all(is.finite(analysis$PC2))
)

# 5. Assemble Q3 analysis table --------------------------------------

robustness_analysis <- robustness_scores %>%
  left_join(
    analysis,
    by = keys,
    relationship = "many-to-one"
  )

# 6. Fit initial continuous-response models --------------------------

ys <- c(
  "BeeDiversity",
  "SyrphidDiversity",
  "MeanFloralDiversity"
)

model.dir <- file.path(dir.SI, "results", "climate_models")
dir.create(model.dir, recursive = TRUE, showWarnings = FALSE)

models <- list()

for (y in ys) {
  
  print(y)
  
  # Keep complete observations for this particular model.
  needed <- c(y, "PC1", "PC2", "Site", "Year")
  dat <- analysis[complete.cases(analysis[, needed]), ]
  
  print(nrow(dat))
  
  equation <- as.formula(
    paste(y, "~ PC1 + PC2 + (1 | Site) + (1 | Year)")
  )
  
  model <- lmer(equation, data = dat)
  
  print(summary(model))
  
  models[[y]] <- model
  
  saveRDS(
    model,
    file.path(model.dir, paste0(y, "_model.rds"))
  )
}