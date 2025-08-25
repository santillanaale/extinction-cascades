rm(list = ls())
library(tidyverse)
library(terra)
library(sf)
setwd("~/")
source("lab_paths.R")
local.path
dir.bombus <- file.path(local.path, "extinction-cascades")

# This script was copied over from skyislands on 8/25/2025. It needs to be updated to match the paths in this repository.

# Precipitation and temperature data used in this analysis were downloaded
# from the PRISM Climate Group’s Explorer portal:
# https://prism.oregonstate.edu/explorer/
#
# - Date of download: August 17, 2025
# - Variables: Precipitation and Temperature
# - Temporal resolution: Monthly
# - Spatial resolution: 800 meters
# - File format: BIL (.bil)
# - Dataset used: PRISM Climate Group, Oregon State University


# ---- Load Site Shapefile ----
setwd(file.path(local.path, "skyIslands_saved"))

# sites_shp <- vect("spatial/sites.shp")
# print(crs(sites_shp))
# 
# # ---- Define path to PRISM data directory ----
# prism_dir <- "data/PRISM_data"
# 
# # ---- Import Precipitation Data ----
# 
# ## Summer Monsoon Precipitation (July, Aug, Sep) --------------------------
# years <- c("2011", "2016", "2017", "2020", "2021")
# months <- c("07", "08", "09")
# base_path_summer <- file.path(prism_dir, "summer_precip")
# 
# summer_rasters <- list()
# 
# for (year in years) {
#   monthly_rasters <- list()
#   for (month in months) {
#     file_path <- file.path(base_path_summer, year, paste0("PRISM_ppt_stable_4kmM3_", year, month, "_bil.bil"))
#     monthly_rasters[[month]] <- rast(file_path)
#   }
#   stacked_rasters <- rast(monthly_rasters)
#   summer_rasters[[year]] <- sum(stacked_rasters, na.rm = TRUE)
#   print(paste("Processed summer precipitation for year:", year))
# }
# 
# ## ---- Winter Precipitation (Dec previous year + Jan current year) ----------
# winter_years <- c("2012", "2017", "2018", "2021", "2022")
# base_path_winter <- file.path(prism_dir, "winter_precip")
# 
# winter_rasters <- list()
# 
# for (year in winter_years) {
#   dec_year <- as.character(as.numeric(year) - 1)
#   dec_path <- file.path(base_path_winter, dec_year, paste0("PRISM_ppt_stable_4kmM3_", dec_year, "12_bil.bil"))
#   jan_path <- file.path(base_path_winter, year, paste0("PRISM_ppt_stable_4kmM3_", year, "01_bil.bil"))
# 
#   dec_raster <- rast(dec_path)
#   jan_raster <- rast(jan_path)
# 
#   winter_rasters[[year]] <- sum(c(dec_raster, jan_raster), na.rm = TRUE)
#   print(paste("Processed winter precipitation for year:", year))
# }
# 
# # ---- Reproject sites to match raster CRS ------------------------------
# sites_shp <- project(sites_shp, crs(summer_rasters[[1]]))
# 
# # ---- Extract Precipitation Data -------------
# 
# ## ---- Extract Summer Monsoon Precipitation values at sites -----------------
# monsoon_results <- list()
# 
# for (year in years) {
#   monsoon_raster <- summer_rasters[[year]]
#   precip_values <- extract(monsoon_raster, sites_shp)
# 
#   site_df <- as.data.frame(sites_shp)
#   site_df$Monsoon_Precipitation <- precip_values[, 2]  # Extracted values
#   site_df$Year <- as.numeric(year)
# 
#   site_summary <- site_df %>%
#     group_by(Site, Year) %>%
#     summarize(Mean_Monsoon_Precip = mean(Monsoon_Precipitation, na.rm = TRUE), .groups = "drop")
# 
#   monsoon_results[[year]] <- site_summary
# }
# 
# monsoon_precip_data <- bind_rows(monsoon_results)
# 
# ## ---- Extract Winter Precipitation values at sites --------------------------
# winter_results <- list()
# 
# for (year in winter_years) {
#   winter_raster <- winter_rasters[[year]]
#   precip_values <- extract(winter_raster, sites_shp)
# 
#   site_df <- as.data.frame(sites_shp)
#   site_df$Winter_Precipitation <- precip_values[, 2]
#   site_df$Winter_Year <- as.numeric(year)
# 
#   site_summary <- site_df %>%
#     group_by(Site, Winter_Year) %>%
#     summarize(Mean_Winter_Precip = mean(Winter_Precipitation, na.rm = TRUE), .groups = "drop")
# 
#   winter_results[[year]] <- site_summary
# }
# 
# winter_precip_data <- bind_rows(winter_results)
# 
# # ---- Final Prep ----------
# 
# # Rename Winter_Year to Year for consistency
# winter_precip_data <- winter_precip_data %>%
#   rename(Year = Winter_Year)
# 
# # Shift monsoon data by one year (since it's antecedent year monsoon precip)
# monsoon_precip_data <- monsoon_precip_data %>%
#   mutate(Year = Year + 1)
# 
# # ---- Save processed precipitation data for downstream analyses -------------
# save(monsoon_precip_data, winter_precip_data, file = file.path('data/PRISM_data/precipitation_site_data.Rdata'))

# ---- Join precipitation data to spec_net --------

## ---- Load specimen table ----
setwd(dir.bombus)
load("data/spec_net.Rdata")

## ---- Load the processed climate data ----
setwd(file.path(local.path, "skyIslands_saved"))
load("data/PRISM_data/precipitation_site_data.Rdata")

## ---- Join climate datasets to spec.net ---------
spec.net <- spec.net %>%
  left_join(monsoon_precip_data, by = c("Site", "Year")) %>%
  left_join(winter_precip_data,  by = c("Site", "Year"))
