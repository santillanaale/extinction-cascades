library(tidyverse)
library(lubridate)
library(ggplot2)
season_climate <- read_csv("season_climate.csv") %>%
  filter(!Site %in% c("SS", "UK", "VC"))

climate_matrix <- season_climate %>%
  select(
    MeanTmean_anom,
    MeanTmax_anom,
    MeanTmin_anom,
    CumPrecip_anom,
    MeanVPD_anom
  )

head(climate_matrix)

# PCA: What kind of climate departure occurred?
pca <-
  prcomp(
    climate_matrix,
    center = TRUE,
    scale. = TRUE
  )

summary(pca)
pca$rotation

scores <- as.data.frame(pca$x)

exposure <- bind_cols(
  season_climate,
  scores
)

write_csv(exposure,"Exposure.csv")

ggplot(exposure,
       aes(PC1,
           PC2,
           color = Site)) +
  geom_point(size = 3) +
  stat_ellipse()

# # Euclidean distance: How large was the departure?
# exposure$ClimateDistance <-
#   sqrt(
#     exposure$MeanTmean_z^2 +
#       exposure$MeanTmax_z^2 +
#       exposure$MeanTmin_z^2 +
#       exposure$CumPrecip_z^2 +
#       exposure$MeanVPD_z^2
#   )

