rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggfortify)
setwd("~/")
source("lab_paths.R")

climate <- read_csv(
  file.path(local.path, "skyIslands", "data", "sitestats.csv"),
  col_types = cols(SampleRound = col_character())
) %>%
  select(
    Site, Year, SampleRound,
    Tmean_anom, Tmax_anom, Tmin_anom,
    Precip_anom, VPD_anom
  )

setwd(file.path(local.path, "extinction-cascades"))

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
  "data/Exposure_PCA.rds"
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
  "data/Exposure_Loadings.csv"
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
  "data/Exposure_Variance.csv"
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

write.csv(exposure,"data/exposure.csv",row.names = FALSE)

## ===============================================================
## Save exposure axis interpretation
## ===============================================================

exposure_interpretation <- tibble(
  Axis = c("PC1", "PC2"),
  MovingRightMeans = c(
    if (pca$rotation["Tmean_anom", "PC1"] < 0) {
      paste(
        "Cooler conditions and lower vapor pressure deficit (VPD),",
        "with precipitation tending to increase"
      )
    } else {
      paste(
        "Warmer conditions and higher vapor pressure deficit (VPD),",
        "with precipitation tending to decrease"
      )
    },
    if (pca$rotation["Precip_anom", "PC2"] > 0) {
      paste(
        "More precipitation, warmer minimum temperatures,",
        "cooler maximum temperatures, and lower VPD"
      )
    } else {
      paste(
        "Less precipitation, cooler minimum temperatures,",
        "warmer maximum temperatures, and higher VPD"
      )
    }
  ),
  ZeroRepresents = "Pooled center across the rows used to fit the exposure PCA",
  VarianceExplainedPct = round(100 * variance$Variance[1:2], 1)
)

write_csv(exposure_interpretation, "data/Exposure_Interpretation.csv")

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