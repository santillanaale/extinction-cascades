library(tidyverse)
library(lubridate)
library(ggplot2)
season_climate_z <- read.csv("season_climate_z.csv")

season_climate_z <- read_csv("season_climate_z.csv") %>%
  filter(!Site %in% c("SS", "UK", "VC"))

climate_matrix <-
  season_climate_z %>%
  select(
    MeanTmean_z,
    MeanTmax_z,
    MeanTmin_z,
    CumPrecip_z,
    MeanVPD_z
  )

head(climate_matrix)

# PCA: What kind of climate departure occurred?
pca <-
  prcomp(
    climate_matrix,
    center = FALSE,
    scale. = FALSE
  )

summary(pca)
pca$rotation

scores <- as.data.frame(pca$x)

exposure <- bind_cols(
  season_climate_z,
  scores
)

write_csv(exposure,"Exposure.csv")

ggplot(exposure,
       aes(PC1,
           PC2,
           color = Site)) +
  geom_point(size = 3) +
  stat_ellipse()

# Euclidean distance: How large was the departure?
exposure$ClimateDistance <-
  sqrt(
    exposure$MeanTmean_z^2 +
      exposure$MeanTmax_z^2 +
      exposure$MeanTmin_z^2 +
      exposure$CumPrecip_z^2 +
      exposure$MeanVPD_z^2
  )

