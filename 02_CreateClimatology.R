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

## *******************************************************************
## Calculate seasonal climatology
## *******************************************************************
daily_normals <- daily_normals %>%
  mutate(
    Month = month(as.Date(paste0(Date, "-2001"),
                          format = "%B-%d-%Y"))
  ) %>%
  filter(Month >= 5)

climatology <-
  daily_normals %>%
  group_by(Name) %>%
  summarise(
    
    MeanTmean =
      mean(`tmean (degrees C)`),
    
    MeanTmax =
      mean(`tmax (degrees C)`),
    
    MeanTmin =
      mean(`tmin (degrees C)`),
    
    MeanVPD =
      mean(`vpdmax (hPa)`),
    
    CumPrecip =
      sum(`ppt (mm)`),
    
    # CumSolsloped =
    #   sum(`solsloped (MJ/m^2)`),
    
    .groups="drop"
  )

## *******************************************************************
## Join with observed data
## *******************************************************************
season_climate <-
  season_climate %>%
  left_join(
    climatology,
    by=c("Site"="Name"),
    suffix=c("","_normal")
  )

## *******************************************************************
## Standardize
## *******************************************************************
season_climate <-
  season_climate %>%
  group_by(Site) %>%
  mutate(
    
    MeanTmean_z =
      (MeanTmean-MeanTmean_normal)/
      sd(MeanTmean),
    
    MeanTmax_z =
      (MeanTmax-MeanTmax_normal)/
      sd(MeanTmax),
    
    MeanTmin_z =
      (MeanTmin-MeanTmin_normal)/
      sd(MeanTmin),
    
    CumPrecip_z =
      (CumPrecip-CumPrecip_normal)/
      sd(CumPrecip),
    
    MeanVPD_z =
      (MeanVPD-MeanVPD_normal)/
      sd(MeanVPD)
    
  ) %>%
  ungroup()

write_csv(season_climate,"season_climate_z.csv")
