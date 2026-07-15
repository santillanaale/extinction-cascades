library(tidyverse)
library(lubridate)
setwd("../skyIslands_saved/")

sites_shp_path <- file.path("spatial", "sites.shp")
weather_csv <- "data/relational/original/weather.csv"

# Precipitation paths
daily_precip_dir <- "data/PRISM_data/PRISM_daily_precip"

# Temperature paths
daily_temp_dir <- "data/PRISM_data/PRISM_daily_temp"

# Max vapor pressure paths
daily_vpdmax_dir <- "data/PRISM_data/PRISM_daily_vpdmax"

# # Solar irradiance paths
# daily_solslope_dir <- "data/PRISM_data/PRISM_daily_solsloped"

## *******************************************************************
## Create a round level summary of weather (rounds extend over
## multiple days)
## *******************************************************************

weather <- read_csv(weather_csv) %>%
  mutate(
    Year = as.numeric(Year)
  )  %>%
  filter(SampleRound > 0)

# Summarize seasons
season_windows <- weather %>%
  group_by(Site, Year) %>%
  summarize(
    SeasonStart = as.Date(paste0(first(Year), "-05-01")),
    SeasonEnd   = max(Date),
    .groups="drop"
  )

## *******************************************************************
# Read in precipitation data and summarize
## *******************************************************************
daily_precip_files  <- list.files(daily_precip_dir, pattern = "\\.csv$",
                          full.names = TRUE)

daily_precip_data <- daily_precip_files  %>%
  set_names() %>%
  map_dfr(~ read_csv(.x, skip = 10) %>%
            select(Name, `ppt (mm)`, Date) %>%
            rename(Site = Name, Precip = `ppt (mm)`) %>%
            mutate(Date = ymd(Date),
                   Year = year(Date)))

season_precip <-
  daily_precip_data %>%
  inner_join(season_windows,
             by=c("Site","Year")) %>%
  filter(Date >= SeasonStart,
         Date <= SeasonEnd) %>%
  group_by(Site,Year) %>%
  summarise(
    
    CumPrecip=sum(Precip),
    
    .groups="drop")

## *******************************************************************
## Read in temperature data and summarize
## *******************************************************************
daily_temp_files  <- list.files(daily_temp_dir, pattern = "\\.csv$", full.names = TRUE)


daily_temp_data <- daily_temp_files %>%
  set_names() %>%
  map_dfr(~ read_csv(.x, skip = 10) %>%   # adjust skip if needed
            rename(
              Site = Name,
              tmin = `tmin (degrees C)`,
              tmax = `tmax (degrees C)`,
              tmean = `tmean (degrees C)`,
              Date = Date
            ) %>%
            mutate(
              Date = ymd(Date),
              Year = year(Date),
              Month = month(Date)
            ))

season_temp <-
  daily_temp_data %>%
  inner_join(season_windows,
             by=c("Site","Year")) %>%
  filter(Date>=SeasonStart,
         Date<=SeasonEnd) %>%
  group_by(Site,Year) %>%
  summarise(
    
    MeanTmean=mean(tmean),
    
    MeanTmax=mean(tmax),
    
    MeanTmin=mean(tmin),
    
    .groups="drop")

## *******************************************************************
## Read in vapor pressure data and summarize
## *******************************************************************
daily_vpdmax_files  <- list.files(daily_vpdmax_dir, pattern = "\\.csv$", full.names = TRUE)


daily_vpdmax_data <- daily_vpdmax_files %>%
  set_names() %>%
  map_dfr(~ read_csv(.x, skip = 10) %>%
            select(Name, `vpdmax (hPa)`, Date) %>%
            rename(Site = Name, vpdmax = `vpdmax (hPa)`) %>%
            mutate(Date = ymd(Date),
                   Year = year(Date)))

season_vpd <-
  daily_vpdmax_data %>%
  inner_join(season_windows,
             by=c("Site","Year")) %>%
  filter(Date >= SeasonStart,
         Date <= SeasonEnd) %>%
  group_by(Site,Year) %>%
  summarise(
    
    MeanVPD=mean(vpdmax),
    
    .groups="drop")

# ## *******************************************************************
# ## Read in solar irradiance data and summarize
# ## *******************************************************************
# daily_solslope_files  <- list.files(daily_solslope_dir, pattern = "\\.csv$", full.names = TRUE)
# 
# 
# daily_solslope_data <- daily_solslope_files %>%
#   set_names() %>%
#   map_dfr(~ read_csv(.x, skip = 10) %>%
#             select(Name, `solslope (MJ/m^2/day)`, Date) %>%
#             rename(Site = Name, solslope = `solslope (MJ/m^2/day)`) %>%
#             mutate(Date = ymd(Date),
#                    Year = year(Date)))
# 
# season_solar <-
#   daily_solslope_data %>%
#   inner_join(season_windows,
#              by=c("Site","Year")) %>%
#   filter(Date >= SeasonStart,
#          Date <= SeasonEnd) %>%
#   group_by(Site,Year) %>%
#   summarise(
#     
#     CumSolar=sum(solslope),
#     
#     .groups="drop")

## *******************************************************************
## Join seasonal climate datasets
## *******************************************************************
season_climate <-
  
  season_temp %>%
  
  left_join(season_precip) %>%
  
  left_join(season_vpd)
  
  # left_join(season_solar)

setwd("../extinction-cascades/")
write.csv(season_climate, file =
                              "season_climate.csv", row.names=FALSE)