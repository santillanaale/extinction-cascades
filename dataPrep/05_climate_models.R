rm(list=ls())
setwd("~/")
source("lab_paths.R")

library(dplyr)
library(lme4)
library(performance)
library(see)

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

load(file.path(
  dir.SI, "analysis", "network", "saved",
  "corMets_PlantPollinator_YearSR_BeesSyrphidsButterflyHawk_N100.Rdata" 
))

network_metrics <- cor.dats

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

# 3. Assemble Q1/Q2 analysis table -----------------------------------

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
    select(all_of(keys), PollAbundance, PC1, PC2)
)

stopifnot(
  nrow(analysis) == nrow(community),
  all(is.finite(analysis$PC1)),
  all(is.finite(analysis$PC2))
)

# 4. Assemble Q3 analysis table --------------------------------------

# robustness_analysis <- robustness_scores %>%
#   left_join(
#     analysis,
#     by = keys,
#     relationship = "many-to-one"
#   )

# 5. Fit initial continuous-response models --------------------------

# temporary Q1
q1_ys <- c(
  "BeeDiversity",
  "SyrphidDiversity",
  "MeanFloralDiversity"
)

q2_ys <- c(
  # "number.of.species.HL",
  # "number.of.species.LL",
  "functional.complementarity.HL",
  "functional.complementarity.LL",
  "FunRedundancy.Pols",
  "FunRedundancy.Plants",
  "mean.number.of.links.HL",
  "mean.number.of.links.LL"
)

# ys <- q2_ys
#
# model.dir <- file.path(
#   dir.EC, "results", "climate_models", guild, "Q2"
# )
# dir.create(model.dir, recursive = TRUE, showWarnings = FALSE)

ys <- q1_ys

model.dir <- file.path(
  dir.EC, "results", "climate_models", "Q1_continuous"
)
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
  
  print(coef(summary(model)))       # Fixed-effect estimates and standard errors
  print(as.data.frame(VarCorr(model)))  # Random-effect and residual variances
  
  models[[y]] <- model
  
  saveRDS(
    model,
    file.path(model.dir, paste0(y, "_model.rds"))
  )
}

# 6. Check model diagnostics -----------------------------------------

pdf(file.path(model.dir, "model_diagnostics.pdf"), width = 11, height = 8.5)
tryCatch({
  for (y in names(models)) {
    plot(performance::check_model(models[[y]], panel = TRUE))
  }
}, finally = dev.off())

model_checks <- data.frame(
  response = names(models),
  n = vapply(models, nobs, numeric(1)),
  singular = vapply(models, lme4::isSingular, logical(1))
)

print(model_checks)
write.csv(
  model_checks,
  file.path(model.dir, "model_checks.csv"),
  row.names = FALSE
)