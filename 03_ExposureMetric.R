rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggfortify)

setwd("../skyIslands_saved/")
climate <- read_csv("data/relational/original/climate.csv")

setwd("../extinction-cascades/")

## ===============================================================
## Create multivariate climate exposure metric
##
## Uses sample-round climate anomalies to summarize the magnitude
## and direction of climate departures from the 30-year daily
## climatology.
##
## PCA axes represent the "Exposure" component of the climate
## vulnerability framework.
## ===============================================================


## ===============================================================
## Select climate anomaly variables
## ===============================================================
climate_matrix <- climate %>%
  select(
    Tmean_anom,
    Tmax_anom,
    Tmin_anom,
    Precip_anom,
    VPD_anom
  )

## ===============================================================
## Principal Components Analysis
##
## Standardize variables so each contributes equally regardless of
## measurement units.
## ===============================================================
pca <- prcomp(
  climate_matrix,
  center = TRUE,
  scale. = TRUE
)

summary(pca)
pca$rotation

## ===============================================================
## Save PCA object
##
## Preserves rotation matrix, scaling, and scores for later use.
## ===============================================================
saveRDS(
  pca,
  "Exposure_PCA.rds"
)

## ===============================================================
## Save PCA loadings
##
## Used to interpret the ecological meaning of each principal
## component.
## ===============================================================
loadings <- as.data.frame(pca$rotation)

loadings$Variable <- rownames(loadings)

write_csv(
  loadings,
  "Exposure_Loadings.csv"
)

## ===============================================================
## Save variance explained
## ===============================================================
variance <- data.frame(
  PC = paste0("PC", 1:5),
  Variance = summary(pca)$importance[2,]
)

write_csv(
  variance,
  "Exposure_Variance.csv"
)

## ===============================================================
## Create exposure dataset
##
## Append PCA scores to the original climate data.
## ===============================================================
scores <- as.data.frame(pca$x)

exposure <- bind_cols(
  climate,
  scores
)

write_csv(
  exposure,
  "Exposure.csv"
)

write.csv(exposure,"exposure.csv",row.names = FALSE)

## ===============================================================
## Visualization
## ===============================================================

### Site ordination

ggplot(
  exposure,
  aes(
    PC1,
    PC2,
    color = Site
  )
) +
  geom_point(size = 3) +
  stat_ellipse() +
  theme_classic()


### PCA biplot

autoplot(
  pca,
  data = exposure,
  colour = "Site",
  loadings = TRUE,
  loadings.label = TRUE
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  theme_classic()


### Ordination coloured by sampling year

ggplot(
  exposure,
  aes(
    PC1,
    PC2,
    color = factor(Year)
  )
) +
  geom_point(size = 3) +
  theme_classic()


### Ordination coloured by sample round

ggplot(
  exposure,
  aes(
    PC1,
    PC2,
    color = factor(SampleRound)
  )
) +
  geom_point(size = 3) +
  theme_classic()


### Scree plot

fv <- summary(pca)$importance

ggplot(
  data.frame(
    PC = factor(
      paste0("PC", 1:5),
      levels = paste0("PC", 1:5)
    ),
    Variance = fv[2,]
  ),
  aes(
    PC,
    Variance
  )
) +
  geom_col() +
  theme_classic() +
  ylab("Proportion of variance explained")


## ===============================================================
## Interpretation (update as variables or loadings change)
##
## PC1: Warm–dry versus cool–wet climate departures.
##
## PC2: Secondary temperature–precipitation gradient.
##
## Interpretation should always be confirmed using the loading
## matrix (Exposure_Loadings.csv).
## ===============================================================


## ===============================================================
## Future work
##
## - Join Exposure.csv with network robustness results.
## - Model:
##
##     Robustness ~ PC1 + PC2 + (1 | Site)
##
##  When you do the robustness regression,
##  do a separate model for each extinction scenario (vs all in one model with the fixed effect).
## - Compare PCA-based exposure with APi and GDD models.
## ===============================================================


# ===============================================================
# Optional:
#
# Quantify the magnitude of climate departure using Euclidean
# distance in standardized anomaly space.
#
# ClimateDistance <- sqrt(...)
#
# This would represent the magnitude of exposure, whereas the PCA
# axes represent the direction of exposure.
# ===============================================================