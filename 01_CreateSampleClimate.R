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
  filter(Method == "Net")

# Summarize sample rounds
sample_windows <- weather %>%
  group_by(Site, Year, SampleRound) %>%
  summarize(
    StartDate = min(StartDate),
    EndDate   = max(Date),
    .groups = "drop"
  ) %>%
  group_by(Site, Year) %>%
  mutate(
    WindowStart = as.Date(paste0(Year, "-05-01"))
  ) %>%
  ungroup()

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

sample_precip <-
  daily_precip_data %>%
  inner_join(sample_windows,
             by=c("Site","Year")) %>%
  filter(Date >= WindowStart,
         Date <= EndDate) %>%
  group_by(Site,Year,SampleRound) %>%
  summarise(
    CumPrecip=sum(Precip),
    .groups="drop")


## *******************************************************************
## Read in temperature data and summarize
## *******************************************************************
daily_temp_files  <- list.files(daily_temp_dir, pattern = "\\.csv$", full.names = TRUE)


daily_temp_data <- daily_temp_files %>%
  set_names() %>%
  map_dfr(~ read_csv(.x, skip = 10) %>%   
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

sample_temp <-
  daily_temp_data %>%
  inner_join(sample_windows,
             by=c("Site","Year")) %>%
  filter(Date >= WindowStart,
         Date <= EndDate) %>%
  group_by(Site,Year,SampleRound) %>%
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

sample_vpd <-
  daily_vpdmax_data %>%
  inner_join(sample_windows,
             by=c("Site","Year")) %>%
  filter(Date >= WindowStart,
         Date <= EndDate) %>%
  group_by(Site,Year,SampleRound) %>%
  summarise(
    MeanVPD=mean(vpdmax),
    .groups="drop")

## *******************************************************************
## Join seasonal climate datasets
## *******************************************************************
sample_climate <-
  
  sample_temp %>%
  
  left_join(sample_precip,
    by = c("Site","Year","SampleRound")
  ) %>%
  
  left_join(sample_vpd,
    by = c("Site","Year","SampleRound")
  )

setwd("../extinction-cascades/")
write.csv(sample_climate, file = "sample_climate.csv", row.names=FALSE)
write.csv(sample_windows,file = "sample_windows.csv", row.names=FALSE)