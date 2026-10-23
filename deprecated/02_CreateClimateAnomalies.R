library(tidyverse)
library(lubridate)
setwd("../skyIslands_saved/")
# 30-year daily normals path
daily_normals_file <- "data/PRISM_data/30-year_daily_normals/PRISM_ppt_tmin_tmean_tmax_vpdmax_30yr_normal_800m_daily_normals.csv"

daily_normals <- read_csv(daily_normals_file,
    skip = 10
  )

setwd("../extinction-cascades/")
sample_climate <- read.csv("sample_climate.csv")
sample_windows <- read_csv("sample_windows.csv")

## *******************************************************************
## Calculate sample climatology
## *******************************************************************
daily_normals <- daily_normals %>%
  mutate(
    NormalDate = as.Date(
      paste0(Date, "-2000"),
      format = "%B-%d-%Y"
    )
  )

sample_windows <- sample_windows %>%
  mutate(
    NormalStart = as.Date(
      paste0(
        "2000-",
        format(WindowStart, "%m-%d")
      )
    ),
    NormalEnd = as.Date(
      paste0(
        "2000-",
        format(EndDate, "%m-%d")
      )
    )
  )

climatology <-
  daily_normals %>%
  rename(Site = Name) %>%
  inner_join(
    sample_windows,
    by = "Site",
    relationship = "many-to-many"
  ) %>%
  filter(
    NormalDate >= NormalStart,
    NormalDate <= NormalEnd
  ) %>%
  group_by(
    Site,
    Year,
    SampleRound
  ) %>%
  summarise(
    MeanTmean = mean(`tmean (degrees C)`),
    MeanTmax  = mean(`tmax (degrees C)`),
    MeanTmin  = mean(`tmin (degrees C)`),
    CumPrecip = sum(`ppt (mm)`),
    MeanVPD   = mean(`vpdmax (hPa)`),
    .groups = "drop"
  )

## *******************************************************************
## Join with observed data
## *******************************************************************
sample_climate <-
  sample_climate %>%
  left_join(
    climatology,
    by = c("Site","Year","SampleRound"),
    suffix = c("", "_normal")
  )

## *******************************************************************
## Compute Anomalies
## *******************************************************************
sample_climate <-
  sample_climate %>%
  mutate(
    MeanTmean_anom = MeanTmean - MeanTmean_normal,
    MeanTmax_anom  = MeanTmax  - MeanTmax_normal,
    MeanTmin_anom  = MeanTmin  - MeanTmin_normal,
    CumPrecip_anom = CumPrecip - CumPrecip_normal,
    MeanVPD_anom   = MeanVPD   - MeanVPD_normal
  )

write_csv(sample_climate, "sample_climate.csv")
