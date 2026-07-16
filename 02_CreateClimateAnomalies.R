library(tidyverse)
library(lubridate)
setwd("../skyIslands_saved/")
# 30-year daily normals path
daily_normals_file <- "data/PRISM_data/30-year_daily_normals/PRISM_ppt_tmin_tmean_tmax_vpdmax_30yr_normal_800m_daily_normals.csv"

daily_normals <- read_csv(daily_normals_file,
    skip = 10
  )

setwd("../extinction-cascades/")
season_climate <- read.csv("season_climate.csv")
season_windows <- read_csv("season_windows.csv")

## *******************************************************************
## Calculate seasonal climatology
## *******************************************************************
daily_normals <- daily_normals %>%
  mutate(
    NormalDate = as.Date(
      paste0(Date, "-2000"),
      format = "%B-%d-%Y"
    )
  )

season_windows <- season_windows %>%
  mutate(
    NormalStart =
      as.Date("2000-05-01"),
    NormalEnd =
      as.Date(
        paste0(
          "2000-",
          format(SeasonEnd, "%m-%d")
        )
      )
  )

climatology <-
  daily_normals %>%
  rename(
    Site = Name
  ) %>%
  inner_join(
    season_windows,
    by="Site"
  ) %>%
  filter(
    NormalDate >= NormalStart,
    NormalDate <= NormalEnd
  ) %>%
  group_by(
    Site,
    Year
  ) %>%
  summarise(
    MeanTmean =
      mean(`tmean (degrees C)`),
    MeanTmax =
      mean(`tmax (degrees C)`),
    MeanTmin =
      mean(`tmin (degrees C)`),
    CumPrecip =
      sum(`ppt (mm)`),
    MeanVPD =
      mean(`vpdmax (hPa)`),
    .groups="drop"
  )

## *******************************************************************
## Join with observed data
## *******************************************************************
season_climate <-
  season_climate %>%
  left_join(
    climatology,
    by=c("Site","Year"),
    suffix=c("","_normal")
  )

## *******************************************************************
## Compute Anomalies
## *******************************************************************
season_climate <-
  season_climate %>%
  group_by(Site) %>%
  mutate(
    MeanTmean_anom = MeanTmean - MeanTmean_normal,
    MeanTmax_anom  = MeanTmax  - MeanTmax_normal,
    MeanTmin_anom  = MeanTmin  - MeanTmin_normal,
    CumPrecip_anom = CumPrecip - CumPrecip_normal,
    MeanVPD_anom   = MeanVPD   - MeanVPD_normal
  ) %>%
  ungroup()

write_csv(season_climate,"season_climate.csv")
