rm(list=ls())
setwd("C:/Users/ale_s/University of Oregon Dropbox/Alejandro Santillana Fernandez/extinction-cascades")
library(bipartite)
library(fossil)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(networkD3)
library(terra)
library(lme4)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(viridis)
library(webshot)
library(patchwork)
library(ggeffects)
library(lmerTest)

  # ---- Load Networks ----
load("../skyIslands/data/networks/Year_PlantPollinator_Bees.RData") #yearly networks for bees

# ---- Import and Visualize (Sankey) Networks ----

# Select the specific network "CH.2012"
example_network <- nets[["CH.2012"]]
network_name <- "CH.2012" # Set the name explicitly for reference

# Transform the network data for the Sankey plot
data_long <- as.data.frame(example_network) %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname) %>%
  filter(value > 0)

# Rename columns for clarity
colnames(data_long) <- c("source", "target", "value")
data_long$target <- paste(data_long$target, " ", sep = "")

# Create nodes data frame for the plot
nodes <- data.frame(name = c(as.character(data_long$source),
                             as.character(data_long$target)) %>% unique())

# Reformat the data for networkD3
data_long$IDsource <- match(data_long$source, nodes$name) - 1
data_long$IDtarget <- match(data_long$target, nodes$name) - 1

# Generate the Sankey plot
sankey_plot <- sankeyNetwork(Links = data_long, Nodes = nodes,
                             Source = "IDsource", Target = "IDtarget",
                             Value = "value", NodeID = "name",
                             fontSize = 18)

# Display the Sankey plot
sankey_plot

# ---- Import Site Shapefile ----

#The following script imports a shapefile for study sites.
sites_shp <- vect("Sites/sites.shp")

# ---- Geospatial Visualization ----
## ---- Sky Islands Map ----
#The map contextualizes the geographic positioning of the Sky Islands, a biodiversity hotspot spanning the southwestern U.S. and northern Mexico. Sites sampled are in the states of New Mexico and Arizona.
# Load the country or region boundary
world <- ne_countries(scale = "medium", returnclass = "sf")

# (Optional) Subset for specific countries if relevant
world <- world[world$name %in% c("United States of America"), ]

# Get U.S. state boundaries
states <- ne_states(country = "United States of America", returnclass = "sf")


sites <- st_as_sf(sites_shp)

# Plot
ggplot() +
  # Add basemap
  geom_sf(data = world, fill = "lightgray", color = "white") +
  
  # Add your sites
  geom_sf(data = sites, color = "red", size = 3, shape = 21, fill = "yellow") +
  
  # Add state boundaries
  geom_sf(data = states, fill = NA, color = "black", linetype = "dashed") +
  
  # Add optional scalebar and north arrow
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         style = north_arrow_fancy_orienteering) +
  
  # Define map extent for the southwestern U.S.
  coord_sf(xlim = c(-115, -102), ylim = c(30, 38), expand = FALSE) +
  
  # Customize the appearance
  theme_minimal() +
  labs(title = "Sky Island Sites",
       subtitle = "Geographic distribution in the Southwestern U.S.",
       x = "Longitude", y = "Latitude")


# ---- Precipitation Data ----

## ---- Import Monsoonal and Winter Rasters ----
### Import monthly precipitation data for each field season ###

### ---- Antecedent Year Monsoon Precipitation ----
# precip_2011_06_raster <- rast("PRISM_summer_precip/2011/PRISM_ppt_stable_4kmM3_201106_bil.bil")
precip_2011_07_raster <- rast("PRISM_summer_precip/2011/PRISM_ppt_stable_4kmM3_201107_bil.bil")
precip_2011_08_raster <- rast("PRISM_summer_precip/2011/PRISM_ppt_stable_4kmM3_201108_bil.bil")
precip_2011_09_raster <- rast("PRISM_summer_precip/2011/PRISM_ppt_stable_4kmM3_201109_bil.bil")

# precip_2016_06_raster <- rast("PRISM_summer_precip/2016/PRISM_ppt_stable_4kmM3_201606_bil.bil")
precip_2016_07_raster <- rast("PRISM_summer_precip/2016/PRISM_ppt_stable_4kmM3_201607_bil.bil")
precip_2016_08_raster <- rast("PRISM_summer_precip/2016/PRISM_ppt_stable_4kmM3_201608_bil.bil")
precip_2016_09_raster <- rast("PRISM_summer_precip/2016/PRISM_ppt_stable_4kmM3_201609_bil.bil")

# precip_2017_06_raster <- rast("PRISM_summer_precip/2017/PRISM_ppt_stable_4kmM3_201706_bil.bil")
precip_2017_07_raster <- rast("PRISM_summer_precip/2017/PRISM_ppt_stable_4kmM3_201707_bil.bil")
precip_2017_08_raster <- rast("PRISM_summer_precip/2017/PRISM_ppt_stable_4kmM3_201708_bil.bil")
precip_2017_09_raster <- rast("PRISM_summer_precip/2017/PRISM_ppt_stable_4kmM3_201709_bil.bil")

# precip_2020_06_raster <- rast("PRISM_summer_precip/2020/PRISM_ppt_stable_4kmM3_202006_bil.bil")
precip_2020_07_raster <- rast("PRISM_summer_precip/2020/PRISM_ppt_stable_4kmM3_202007_bil.bil")
precip_2020_08_raster <- rast("PRISM_summer_precip/2020/PRISM_ppt_stable_4kmM3_202008_bil.bil")
precip_2020_09_raster <- rast("PRISM_summer_precip/2020/PRISM_ppt_stable_4kmM3_202009_bil.bil")

# precip_2021_06_raster <- rast("PRISM_summer_precip/2021/PRISM_ppt_stable_4kmM3_202106_bil.bil")
precip_2021_07_raster <- rast("PRISM_summer_precip/2021/PRISM_ppt_stable_4kmM3_202107_bil.bil")
precip_2021_08_raster <- rast("PRISM_summer_precip/2021/PRISM_ppt_stable_4kmM3_202108_bil.bil")
precip_2021_09_raster <- rast("PRISM_summer_precip/2021/PRISM_ppt_stable_4kmM3_202109_bil.bil")

# Define the years and file paths
years <- c("2011", "2016", "2017", "2020", "2021")
months <- c("07", "08", "09")
base_path <- "PRISM_summer_precip"

# Initialize a list to store the summed rasters for each year
summer_rasters <- list()

# Loop through each year
for (year in years) {
  
  # Initialize an empty list to hold monthly rasters
  monthly_rasters <- list()
  
  # Loop through each month
  for (month in months) {
    # Construct the file path
    file_path <- file.path(base_path, year, paste0("PRISM_ppt_stable_4kmM3_", year, month, "_bil.bil"))
    
    # Load the raster and append to the list
    monthly_rasters[[month]] <- rast(file_path)
  }
  
  # Stack and sum the monthly rasters
  stacked_rasters <- rast(monthly_rasters)
  summer_raster <- sum(stacked_rasters, na.rm = TRUE)
  
  # Save the raster to the list
  summer_rasters[[year]] <- summer_raster
  
  
  # Print a message for progress tracking
  print(paste("Processed year:", year))
}

antecedent_monsoon_precip_2011_raster <- summer_rasters[["2011"]]
antecedent_monsoon_precip_2016_raster <- summer_rasters[["2016"]]
antecedent_monsoon_precip_2017_raster <- summer_rasters[["2017"]]
antecedent_monsoon_precip_2020_raster <- summer_rasters[["2020"]]
antecedent_monsoon_precip_2021_raster <- summer_rasters[["2021"]]

### ---- Winter Precipitation ----
# Define the years for which you want to calculate winter precipitation.
# Each year represents the January of the winter (e.g., winter 2011-2012 is represented by 2012)
winter_years <- c("2012", "2017", "2018", "2021", "2022")

# Base directory where your PRISM rasters are stored
base_path <- "PRISM_winter_precip"

# Initialize list to store the summed winter precipitation rasters
winter_rasters <- list()

# Loop through each winter year
for (year in winter_years) {
  
  # Determine the December year (previous calendar year)
  dec_year <- as.character(as.numeric(year) - 1)
  jan_year <- year
  
  # Construct file paths for December and January rasters
  dec_path <- file.path(base_path, dec_year, paste0("PRISM_ppt_stable_4kmM3_", dec_year, "12_bil.bil"))
  jan_path <- file.path(base_path, jan_year, paste0("PRISM_ppt_stable_4kmM3_", jan_year, "01_bil.bil"))
  
  # Load the rasters
  dec_raster <- rast(dec_path)
  jan_raster <- rast(jan_path)
  
  # Stack the rasters and calculate total winter precipitation
  winter_stack <- c(dec_raster, jan_raster)
  winter_total <- sum(winter_stack, na.rm = TRUE)
  
  # Save to the list
  winter_rasters[[year]] <- winter_total
  
  # Progress message
  print(paste("Processed winter ending in:", year))
}

# Optionally assign individual years to named variables
winter_precip_2012_raster <- winter_rasters[["2012"]]
winter_precip_2017_raster <- winter_rasters[["2017"]]
winter_precip_2018_raster <- winter_rasters[["2018"]]
winter_precip_2021_raster <- winter_rasters[["2021"]]
winter_precip_2022_raster <- winter_rasters[["2022"]]

# Identify variables with the pattern "precip_****_**_raster"
monthly_rasters <- ls(pattern = "^precip_\\d{4}_\\d{2}_raster$")
rm(list = monthly_rasters)

## ---- Extract Precipitation Data for Sites ----

### ---- Transform Site Spatial Data ----
#The following script imports a shapefile for study sites, ensuring consistency in coordinate reference systems (CRS) with PRISM rasters.
print(crs(sites_shp))
print(crs(winter_precip_2012_raster))
sites_shp <- project(sites_shp, crs(winter_precip_2012_raster))

### ---- Antecedent Monsoon ----
# Define the years you're interested in
years <- c(2011, 2016, 2017, 2020, 2021)

# Path to directory with monsoon rasters
base_path <- "PRISM_summer_precip"

# Initialize list to store results
monsoon_results <- list()

for (year in years) {
  
  # Define file paths for July, August, September
  months <- c("07", "08", "09")
  file_paths <- file.path(base_path, year, paste0("PRISM_ppt_stable_4kmM3_", year, months, "_bil.bil"))
  
  # Load and sum the 3 monthly rasters
  monthly_rasters <- lapply(file_paths, rast)
  monsoon_raster <- sum(rast(monthly_rasters), na.rm = TRUE)
  
  # Extract monsoon precip values at site locations
  precip_values <- extract(monsoon_raster, sites_shp)
  
  # Join with site shapefile data
  subsite_data <- as.data.frame(sites_shp)
  subsite_data$Monsoon_Precipitation <- precip_values[, 2]  # Column 2 has the values
  subsite_data$Year <- year
  
  # Group by Site and Year to get mean per site
  site_data <- subsite_data %>%
    group_by(Site, Year) %>%
    summarize(Mean_Monsoon_Precip = mean(Monsoon_Precipitation, na.rm = TRUE), .groups = "drop")
  
  # Store this year's data
  monsoon_results[[as.character(year)]] <- site_data
}

# Combine all years
monsoon_precip_data <- bind_rows(monsoon_results)

### ---- Winter ----
# Define winters by the year of the January raster
winter_years <- c(2012, 2017, 2018, 2021, 2022)

# Base directory for winter rasters
base_path <- "PRISM_winter_precip"

# Initialize list to store results
winter_results <- list()

for (year in winter_years) {
  
  # Get December of previous year and January of current year
  dec_year <- as.character(as.numeric(year) - 1)
  jan_year <- as.character(year)
  
  # Construct file paths
  dec_path <- file.path(base_path, dec_year, paste0("PRISM_ppt_stable_4kmM3_", dec_year, "12_bil.bil"))
  jan_path <- file.path(base_path, jan_year, paste0("PRISM_ppt_stable_4kmM3_", jan_year, "01_bil.bil"))
  
  # Load rasters
  dec_raster <- rast(dec_path)
  jan_raster <- rast(jan_path)
  
  # Sum December and January to get winter precipitation
  winter_raster <- sum(c(dec_raster, jan_raster), na.rm = TRUE)
  
  # Extract values at site locations
  precip_values <- extract(winter_raster, sites_shp)
  
  # Combine with shapefile data
  subsite_data <- as.data.frame(sites_shp)
  subsite_data$Winter_Precipitation <- precip_values[, 2]
  subsite_data$Winter_Year <- as.numeric(year)
  
  # Group by site and year, then average
  site_data <- subsite_data %>%
    group_by(Site, Winter_Year) %>%
    summarize(Mean_Winter_Precip = mean(Winter_Precipitation, na.rm = TRUE), .groups = "drop")
  
  # Store in results list
  winter_results[[as.character(year)]] <- site_data
}

# Combine all years into a single data frame
winter_precip_data <- bind_rows(winter_results)

## ---- Antecedent Monsoon Precipitation Map ----
# This code produces a faceted map that displays the spatial distribution of summer precipitation for the years 2011, 2016, 2017, 2020, and 2021.
# Define a list of the rasters
rasters <- list(
  precip_2011 = antecedent_monsoon_precip_2011_raster,
  precip_2016 = antecedent_monsoon_precip_2016_raster,
  precip_2017 = antecedent_monsoon_precip_2017_raster,
  precip_2020 = antecedent_monsoon_precip_2020_raster,
  precip_2021 = antecedent_monsoon_precip_2021_raster
)

# Define a list of corresponding years for labeling
years <- c("2011", "2016", "2017", "2020", "2021")

# Calculate global minimum and maximum precipitation values across all rasters
min_precip <- min(sapply(rasters, function(r) min(values(r), na.rm = TRUE)))
max_precip <- max(sapply(rasters, function(r) max(values(r), na.rm = TRUE)))

# Create an empty list to store all the data for ggplot
plot_data <- list()

# Loop through each raster and prepare data for ggplot
for (i in 1:length(rasters)) {
  
  # Get the current raster and year
  current_raster <- rasters[[i]]
  current_year <- years[i]
  
  # Convert the current raster to a data.frame for ggplot
  precip_df <- as.data.frame(current_raster, xy = TRUE, na.rm = TRUE)  # Remove NAs
  
  # Rename the third column to "precip"
  colnames(precip_df)[3] <- "precip"
  
  # Apply log transformation to precipitation values (log(x + 1) to avoid log(0))
  precip_df$precip_log <- log(precip_df$precip + 1)
  
  # Add the year to the data for faceting
  precip_df$year <- current_year
  
  # Store the data in the list
  plot_data[[i]] <- precip_df
}

# Combine all the data into a single data frame
plot_data_df <- do.call(rbind, plot_data)

# Convert SpatVector to sf object
sites <- st_as_sf(sites_shp)

# Create the plot with facet wrap
plot <- ggplot(plot_data_df) +
  # Add the precipitation raster layer
  geom_raster(aes(x = x, y = y, fill = precip_log)) +
  
  # Add a color scale for the precipitation with fixed limits for consistency
  scale_fill_viridis(option = "D", name = "Precipitation",
                     direction = 1) +
  # limits = c(log(min_precip + 1), log(max_precip + 1))) +
  
  # Add the site locations
  geom_sf(data = sites, color = "black", size = 3, shape = 21, fill = "grey") +
  
  # Customize map extent (adjust as needed)
  coord_sf(xlim = c(-111, -105), ylim = c(31, 37), expand = FALSE) +
  
  # Customize the appearance
  theme_minimal() +
  labs(
    title = "Sky Island Sites and Precipitation",
    subtitle = "Summer precipitation layer for different years (log-transformed)",
    x = "Longitude", 
    y = "Latitude"
  ) +
  # Round axis labels
  scale_x_continuous(breaks = seq(-110, -106, by = 2), 
                     labels = function(x) round(x, 2)) +  # Adjust longitude rounding
  scale_y_continuous(breaks = seq(31, 37, by = 1), 
                     labels = function(y) round(y, 2)) +  # Adjust latitude rounding
  
  # Facet wrap by year
  facet_wrap(~ year, nrow=1)

# Print the plot
print(plot)

# # Extract Average Annual Precipitation Data for Sites
# # This code extracts precipitation data from raster files for specified years and calculates average precipitation for each study site.
# years <- c(2011, 2016, 2017, 2020, 2021)
# results_list <- list()
# 
# for (year in years) {
#   # Load the raster for the year
#   raster <- rast(paste0("PRISM_precip/PRISM_ppt_stable_4kmM3_", year, "_bil.bil")) 
#   
#   # Extract precipitation values for subsites
#   precip_values <- extract(raster, sites_shp)
#   
#   # Combine extracted values with site shapefile data
#   subsite_data <- as.data.frame(sites_shp) # Convert shapefile to a data frame
#   subsite_data$Precipitation <- precip_values[, 2] # Add extracted values
#   subsite_data$Year <- year # Add the year
#   
#   # Group by Site and Year, then calculate mean precipitation
#   site_data <- subsite_data %>%
#     group_by(Site, Year) %>%
#     summarize(Mean_Precipitation = mean(Precipitation, na.rm = TRUE), .groups = "drop")
#   
#   # Store the yearly averaged data in the results list
#   results_list[[as.character(year)]] <- site_data
# }
# 
# # Combine all yearly site-level results into a single data frame
# annual_precip_data <- do.call(rbind, results_list)

# ---- Temperature Data ----
# Extract site coords for temperature data download
coords <- crds(sites_shp, df = TRUE)
coords$Site <- sites_shp$Site  # Adjust to match your actual site name column
# Average lat/lon per Site
sites_avg <- coords %>%
  group_by(Site) %>%
  summarize(
    latitude = mean(y),
    longitude = mean(x),
    .groups = "drop"
  )

## ---- Import Daily min max Temp Data ----
temp_2012 <- read.csv("PRISM_temp/PRISM_tmin_tmax_stable_800m_20120101_20121231.csv", skip = 10, stringsAsFactors = FALSE)
temp_2017 <- read.csv("PRISM_temp/PRISM_tmin_tmax_stable_800m_20170101_20171231.csv", skip = 10, stringsAsFactors = FALSE)
temp_2018 <- read.csv("PRISM_temp/PRISM_tmin_tmax_stable_800m_20180101_20181231.csv", skip = 10, stringsAsFactors = FALSE)
temp_2021 <- read.csv("PRISM_temp/PRISM_tmin_tmax_stable_800m_20210101_20211231.csv", skip = 10, stringsAsFactors = FALSE)
temp_2022 <- read.csv("PRISM_temp/PRISM_tmin_tmax_stable_800m_20220101_20220131.csv", skip = 10, stringsAsFactors = FALSE)

# Process each year's temp df
process_temp_data <- function(df) {
  df %>%
    rename(
      site_name = Name,
      longitude = Longitude,
      latitude = Latitude,
      elevation = Elevation..m.,
      date = Date,
      tmin = tmin..degrees.C.,
      tmax = tmax..degrees.C.
    ) %>%
    mutate(
      date = as.Date(date),
      year = year(date),
      month = month(date),
      gdd = calculate_gdd(tmin, tmax, base_temp = 10)
    )
}

## ---- Calculate Growing Degree Days ----
calculate_gdd <- function(tmin, tmax, base_temp = 10, upper_temp = NA) {
  tmin <- pmax(tmin, base_temp)
  tmax <- pmax(tmax, base_temp)
  
  if (!is.na(upper_temp)) {
    tmin <- pmin(tmin, upper_temp)
    tmax <- pmin(tmax, upper_temp)
  }
  
  gdd <- (tmax + tmin) / 2 - base_temp
  gdd[gdd < 0] <- 0
  return(gdd)
}

# Process each year
temp_2012_clean <- process_temp_data(temp_2012)
temp_2017_clean <- process_temp_data(temp_2017)
temp_2018_clean <- process_temp_data(temp_2018)
temp_2021_clean <- process_temp_data(temp_2021)
temp_2022_clean <- process_temp_data(temp_2022)

# Combine into one df
all_temp_data <- bind_rows(
  temp_2012_clean,
  temp_2017_clean,
  temp_2018_clean,
  temp_2021_clean,
  temp_2022_clean
)

## ---- Summarize by site and year ----
annual_gdd <- all_temp_data %>%
  group_by(site_name, year) %>%
  summarize(total_gdd = sum(gdd, na.rm = TRUE), .groups = "drop")

## ---- Monsoon season GDD ----
monsoon_gdd <- all_temp_data %>%
  filter(month %in% 7:9) %>%
  group_by(site_name, year) %>%
  summarize(monsoon_gdd = sum(gdd, na.rm = TRUE), .groups = "drop")

## ---- Winter GDD ----
# Add winter year column
all_temp_data <- all_temp_data %>%
  mutate(
    winter_year = case_when(
      month == 12 ~ year(date) + 1,
      month %in% 1:2 ~ year(date),
      TRUE ~ NA_real_
    )
  )

winter_gdd <- all_temp_data %>%
  filter(month %in% c(12, 1)) %>%
  filter(!is.na(winter_year)) %>%
  group_by(site_name, winter_year) %>%
  summarize(winter_gdd = sum(gdd, na.rm = TRUE), .groups = "drop") %>%
  rename(Year = winter_year)


# ---- Perform Extinction Simulations and Calculate Robustness ----

# This code simulates species extinctions within plant-pollinator networks, calculates their robustness (a measure of how well the network resists species loss), and 
# organizes the results into a data frame for analysis.


## ---- TCM Topological Coextinction Model ----

### ----  Increasing abundance  ----
# # Create an empty data frame to store the results
# TCM_increasing_abundance_robustness_results <- data.frame(Network = character(), Site = character(), Year = integer(), Robustness = numeric(), stringsAsFactors = FALSE)
# 
# # Loop through each network in the list `nets`
# for (net_name in names(nets)) {
#   # Extract the network
#   network <- nets[[net_name]]
#   
#   # Run extinction simulation
#   network_ext <- second.extinct(network, participant = "higher", method = "abun", nrep = 50, details = FALSE)
#   
#   # Calculate the robustness (area under second.extinct curve)
#   network_robustness <- robustness(network_ext)
#   
#   # Extract Site and Year from the network name
#   site_year <- unlist(strsplit(net_name, "\\.")) # Split by "."
#   site <- site_year[1]
#   year <- as.integer(site_year[2])
#   
#   # Append the results to the data frame
#   TCM_increasing_abundance_robustness_results <- rbind(TCM_increasing_abundance_robustness_results, data.frame(Network = net_name, Site = site, Year = year, Robustness = network_robustness))
# }
# 
# # Print the resulting data frame
# print(TCM_increasing_abundance_robustness_results)
# 
### ----  Decreasing abundance  ----
# # Create an empty data frame to store the results
# TCM_decreasing_abundance_robustness_results <- data.frame(Network = character(), Site = character(), Year = integer(), Robustness = numeric(), stringsAsFactors = FALSE)
# 
# source("modified_extinction.R")
# source("modified_second_extinct.R")
# 
# # Sanity Check
# modified.second.extinct(web = example_network, participant = "higher", method = "abundance", reverse = TRUE,     # 👈 removes most abundant species instead of least
#   nrep = 1, details = FALSE)
# 
# # Loop through each network in the list `nets`
# for (net_name in names(nets)) {
#   # Extract the network
#   network <- nets[[net_name]]
# 
#   # Run extinction simulation
#   network_ext <- modified.second.extinct(network, participant = "higher", method = "abun", reverse = TRUE, nrep = 50, details = FALSE)
# 
#   # Calculate the robustness (area under second.extinct curve)
#   network_robustness <- robustness(network_ext)
# 
#   # Extract Site and Year from the network name
#   site_year <- unlist(strsplit(net_name, "\\.")) # Split by "."
#   site <- site_year[1]
#   year <- as.integer(site_year[2])
# 
#   # Append the results to the data frame
#   TCM_decreasing_abundance_robustness_results <- rbind(TCM_decreasing_abundance_robustness_results, data.frame(Network = net_name, Site = site, Year = year, Robustness = network_robustness))
# }
# 
# # Print the resulting data frame
# print(TCM_decreasing_abundance_robustness_results)
# 
### ----  Decreasing Degree  ----
# # Create an empty data frame to store the results
# TCM_decreasing_degree_robustness_results <- data.frame(Network = character(), Site = character(), Year = integer(), Robustness = numeric(), stringsAsFactors = FALSE)
# 
# # Loop through each network in the list `nets`
# for (net_name in names(nets)) {
#   # Extract the network
#   network <- nets[[net_name]]
#   
#   # Run extinction simulation
#   network_ext <- second.extinct(network, participant = "higher", method = "degree", nrep = 50, details = FALSE)
#   
#   # Calculate the robustness (area under second.extinct curve)
#   network_robustness <- robustness(network_ext)
#   
#   # Extract Site and Year from the network name
#   site_year <- unlist(strsplit(net_name, "\\.")) # Split by "."
#   site <- site_year[1]
#   year <- as.integer(site_year[2])
#   
#   # Append the results to the data frame
#   TCM_decreasing_degree_robustness_results <- rbind(TCM_decreasing_degree_robustness_results, data.frame(Network = net_name, Site = site, Year = year, Robustness = network_robustness))
# }
# 
# # Print the resulting data frame
# print(TCM_decreasing_degree_robustness_results)
# 
### ----  Increasing Degree  ----
# # Create an empty data frame to store the results
# TCM_increasing_degree_robustness_results <- data.frame(Network = character(), Site = character(), Year = integer(), Robustness = numeric(), stringsAsFactors = FALSE)
# 
# # Loop through each network in the list `nets`
# for (net_name in names(nets)) {
#   # Extract the network
#   network <- nets[[net_name]]
#   
#   # Run extinction simulation
#   network_ext <- modified.second.extinct(network, participant = "higher", method = "degree", reverse = TRUE, nrep = 50, details = FALSE)
#   
#   # Calculate the robustness (area under second.extinct curve)
#   network_robustness <- robustness(network_ext)
#   
#   # Extract Site and Year from the network name
#   site_year <- unlist(strsplit(net_name, "\\.")) # Split by "."
#   site <- site_year[1]
#   year <- as.integer(site_year[2])
#   
#   # Append the results to the data frame
#   TCM_increasing_degree_robustness_results <- rbind(TCM_increasing_degree_robustness_results, data.frame(Network = net_name, Site = site, Year = year, Robustness = network_robustness))
# }
# 
# # Print the resulting data frame
# print(TCM_decreasing_abundance_robustness_results)

### ---- TCM: Combining all scenarios into one loop ----

# Load modified extinction functions
source("modified_extinction.R")
source("modified_second_extinct.R")

# Use only modified.second.extinct
scenarios <- list(
  abun = list(method = "abun", reverse = FALSE),
  abun_reverse = list(method = "abun", reverse = TRUE),
  degree = list(method = "degree", reverse = FALSE),
  degree_reverse = list(method = "degree", reverse = TRUE)
)


for (scenario_name in names(scenarios)) {
  scenario <- scenarios[[scenario_name]]
  
  results_list <- list()
  
  for (net_name in names(nets)) {
    network <- nets[[net_name]]
    
    # Run extinction simulation using modified function for all
    ext_result <- modified.second.extinct(
      web = network,
      participant = "higher",
      method = scenario$method,
      reverse = scenario$reverse,
      nrep = 50,
      details = FALSE
    )
    
    rob <- robustness(ext_result)
    
    site_year <- unlist(strsplit(net_name, "\\."))
    site <- site_year[1]
    year <- as.integer(site_year[2])
    
    results_list[[net_name]] <- data.frame(
      Network = net_name,
      Site = site,
      Year = year,
      Robustness = rob,
      stringsAsFactors = FALSE
    )
  }
  
  scenario_df <- do.call(rbind, results_list)
  
  # Save results
  save_path <- sprintf("analysis/network/robustness_metrics/%s_robustness_metrics.Rdata", scenario_name)
  save(scenario_df, file = save_path)
  
  assign(paste0("TCM_", scenario_name, "_robustness_results"), scenario_df)
  print(scenario_df)
}

## ---- SCM: Stochastic Coextinction Model ----

# #source('SCM/VieraAlmeida2015RFunctions/netcascade (April 2014).R')
# source("SCM/Dalsgaard2018Code/IterNodeDelMultiSim.R")
# 
# # # read in R values
# # rvals <- read.csv("SCM/Dalsgaard2018Code/data/rvalues.csv")
# # rvals$R <- trimws(rvals$R)
# 
# # Generate randomized rvals data frame
# set.seed(42)  # For reproducibility
# 
# # Define R categories as used in your bounds
# r_categories <- c("low", "med", "high", "low-med", "med-high", "low-high")
# 
# # Initialize list to collect rvals for each network
# rvals_list <- list()
# 
# # Loop through networks
# for (net_name in names(nets)) {
#   network <- nets[[net_name]]
#   
#   # Get plant species (assumes plants are rows)
#   plant_species <- rownames(network)
#   
#   # Assign each plant species a random R category
#   assigned_Rs <- sample(r_categories, length(plant_species), replace = TRUE)
#   
#   # Build data frame
#   df <- data.frame(
#     species = plant_species,
#     R = assigned_Rs,
#     file_name = net_name,
#     stringsAsFactors = FALSE
#   )
#   
#   rvals_list[[net_name]] <- df
# }
# 
# # Combine into one data frame
# rvals <- do.call(rbind, rvals_list)
# 
# # Create the r_bounds table matching your original logic
# r_bounds <- data.frame(
#   R = r_categories,
#   lb = NA_real_,
#   ub = NA_real_,
#   stringsAsFactors = FALSE
# )
# 
# # Define the bounds
# r_bounds[r_bounds$R == "low", c("lb", "ub")] <- c(0, 1/3)
# r_bounds[r_bounds$R == "med", c("lb", "ub")] <- c(1/3, 2/3)
# r_bounds[r_bounds$R == "high", c("lb", "ub")] <- c(2/3, 1)
# r_bounds[r_bounds$R == "low-med", c("lb", "ub")] <- c(0, 2/3)
# r_bounds[r_bounds$R == "med-high", c("lb", "ub")] <- c(1/3, 1)
# r_bounds[r_bounds$R == "low-high", c("lb", "ub")] <- c(0, 1)
# 
# # Optional: save to file
# write.csv(rvals, "rvals_randomized.csv", row.names = FALSE)
# 
# # Now rvals and r_bounds are ready to be passed into your simulation
# 
# # split rvals by web
# rvals <- split(rvals, rvals$file_name)
# 
# vulnerabilities <- NULL
# for(i in 1:length(nets)){
#   
#   # run simulations
#   most <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "bo", removal.taxa = "plant", remove.order = "from most connected", nsim = 10000)
#   strongest <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "bo", removal.taxa = "plant", remove.order = "from strongest", nsim = 10000)
#   random <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "bo", removal.taxa = "plant", remove.order = "random", nsim = 10000)
#   weakest <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "bo", removal.taxa = "plant", remove.order = "from weakest", nsim = 10000)
#   least <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "bo", removal.taxa = "plant", remove.order = "from least connected", nsim = 10000)
#   
#   r.most <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "random", removal.taxa = "plant", remove.order = "from most connected", nsim = 10000)
#   r.strongest <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "random", removal.taxa = "plant", remove.order = "from strongest", nsim = 10000)
#   r.random <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "random", removal.taxa = "plant", remove.order = "random", nsim = 10000)
#   r.weakest <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "random", removal.taxa = "plant", remove.order = "from weakest", nsim = 10000)
#   r.least <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "random", removal.taxa = "plant", remove.order = "from least connected", nsim = 10000)
#   
#   # gather results
#   most <- merge(most$mean_order, most$mean_vulnerability, by = "species")
#   most$web <- names(nets)[i]
#   most$removal.order <- "from most connected"
#   most$plantR <- "assigned"
#   
#   strongest <- merge(strongest$mean_order, strongest$mean_vulnerability, by = "species")
#   strongest$web <- names(nets)[i]
#   strongest$removal.order <- "from strongest"
#   strongest$plantR <- "assigned"
#   
#   random <- merge(random$mean_order, random$mean_vulnerability, by = "species")
#   random$web <- names(nets)[i]
#   random$removal.order <- "random"
#   random$plantR <- "assigned"
#   
#   weakest <- merge(weakest$mean_order, weakest$mean_vulnerability, by = "species")
#   weakest$web <- names(nets)[i]
#   weakest$removal.order <- "from weakest"
#   weakest$plantR <- "assigned"
#   
#   least <- merge(least$mean_order, least$mean_vulnerability, by = "species")
#   least$web <- names(nets)[i]
#   least$removal.order <- "from least connected"
#   least$plantR <- "assigned"
#   
#   r.most <- merge(r.most$mean_order, r.most$mean_vulnerability, by = "species")
#   r.most$web <- names(nets)[i]
#   r.most$removal.order <- "from most connected"
#   r.most$plantR <- "random"
#   
#   r.strongest <- merge(r.strongest$mean_order, r.strongest$mean_vulnerability, by = "species")
#   r.strongest$web <- names(nets)[i]
#   r.strongest$removal.order <- "from strongest"
#   r.strongest$plantR <- "random"
#   
#   r.random <- merge(r.random$mean_order, r.random$mean_vulnerability, by = "species")
#   r.random$web <- names(nets)[i]
#   r.random$removal.order <- "random"
#   r.random$plantR <- "random"
#   
#   r.weakest <- merge(r.weakest$mean_order, r.weakest$mean_vulnerability, by = "species")
#   r.weakest$web <- names(nets)[i]
#   r.weakest$removal.order <- "from weakest"
#   r.weakest$plantR <- "random"
#   
#   r.least <- merge(r.least$mean_order, r.least$mean_vulnerability, by = "species")
#   r.least$web <- names(nets)[i]
#   r.least$removal.order <- "from least connected"
#   r.least$plantR <- "random"
#   
#   vulnerabilities <- rbind(vulnerabilities,rbind(most,strongest,random,weakest,least,r.most,r.strongest,r.random,r.weakest,r.least))
#   
#   print(paste0(i,"/8"))
# }
# 
# vulnerabilities$index <- vulnerabilities$probability.of.extinction * (1-vulnerabilities$mean.order)
# # Key of how column names map onto vulnerability indices from the paper:
# # probability.of.extinction = PE
# # mean.order = SE
# # index = VE

## ---- Plotting Robustness ----
# Choose the robustness data (e.g., for 'abun' scenario)
robustness_df <- TCM_abun_robustness_results

# Correct for antecedent monsoon year
monsoon_precip_data_df <- monsoon_precip_data %>%
  mutate(Year = Year + 1)  # Shift year forward by one

# Change column names for winter precip
winter_precip_data_df <- winter_precip_data %>%
  rename(Year = Winter_Year)

# Change column names for GDD
annual_gdd_df <- annual_gdd %>%
  rename(Site = site_name, Year = year)

# Merge all climate data into robustness_df
robustness_climate <- robustness_df %>%
  left_join(monsoon_precip_data_df, by = c("Site", "Year")) %>%
  left_join(winter_precip_data_df, by = c("Site", "Year")) %>%
  left_join(annual_gdd_df, by = c("Site", "Year"))

# Plot
p1 <- ggplot(robustness_climate, aes(x = Mean_Monsoon_Precip, y = Robustness)) +
  geom_point(aes(color = Site)) + geom_smooth(method = "lm") + theme_minimal() +
  labs(title = "Monsoon Precip")

p2 <- ggplot(robustness_climate, aes(x = Mean_Winter_Precip, y = Robustness)) +
  geom_point(aes(color = Site)) + geom_smooth(method = "lm") + theme_minimal() +
  labs(title = "Winter Precip")

p3 <- ggplot(robustness_climate, aes(x = total_gdd, y = Robustness)) +
  geom_point(aes(color = Site)) + geom_smooth(method = "lm") + theme_minimal() +
  labs(title = "GDD")

# Combine into 1 row
p1 + p2 + p3 + plot_layout(ncol = 3)
ggsave("analysis/network/figures/RobustnessbyClimate.pdf", width = 12, height = 4)

## ---- Modeling Robustness ----
# Standardize predictors
robustness_climate <- robustness_climate %>%
  mutate(across(c(Mean_Monsoon_Precip, Mean_Winter_Precip, total_gdd),
                ~ scale(.)[,1],  # extract the numeric column from scale()
                .names = "z_{.col}"))
# 
# # Individual models
# mod_monsoon <- lmer(Robustness ~ z_Mean_Monsoon_Precip + (1 | Site), data = robustness_climate)
# mod_winter <- lmer(Robustness ~ z_Mean_Winter_Precip + (1 | Site), data = robustness_climate)
# mod_gdd    <- lmer(Robustness ~ z_total_gdd + (1 | Site), data = robustness_climate)

# Combined model
mod_all <- lmer(Robustness ~ z_Mean_Monsoon_Precip + z_Mean_Winter_Precip + z_total_gdd + (1 | Site),
                data = robustness_climate)
summary(mod_all)

# Get p-values
model_summary <- summary(mod_all)$coefficients

# Identify significance of each fixed effect (p < 0.05)
sig_monsoon <- model_summary["z_Mean_Monsoon_Precip", "Pr(>|t|)"] < 0.05
sig_winter  <- model_summary["z_Mean_Winter_Precip", "Pr(>|t|)"] < 0.05
sig_gdd     <- model_summary["z_total_gdd", "Pr(>|t|)"] < 0.05

# Predict effects for each variable
monsoon_eff <- ggpredict(mod_all, terms = "z_Mean_Monsoon_Precip")
winter_eff  <- ggpredict(mod_all, terms = "z_Mean_Winter_Precip")
gdd_eff     <- ggpredict(mod_all, terms = "z_total_gdd")

# Plot each effect separately (clean & labeled)
plot_effect <- function(pred_df, xlab, title, significant = FALSE) {
  line_type <- ifelse(significant, "solid", "dashed")
  
  ggplot(pred_df, aes(x = x, y = predicted)) +
    geom_line(size = 1, linetype = line_type) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
    labs(x = xlab, y = "Predicted Robustness", title = title) +
    theme_minimal(base_size = 14)
}

p1 <- plot_effect(monsoon_eff, "Standardized Monsoon Precipitation", "Effect of Monsoon Precipitation", significant = sig_monsoon)
p2 <- plot_effect(winter_eff,  "Standardized Winter Precipitation", "Effect of Winter Precipitation", significant = sig_winter)
p3 <- plot_effect(gdd_eff,     "Standardized Growing Degree Days (GDD)", "Effect of GDD", significant = sig_gdd)

combined_plot <- p1 + p2 + p3 + plot_layout(ncol = 1)
combined_plot
ggsave("analysis/network/figures/Robustness_ClimateEffects_Combined.pdf", plot = combined_plot, width = 8, height = 10)


# ---- Modeling Community Resistance ----
# Metrics: Network redundancy, complementarity, and generalization

setwd("C:/Users/ale_s/University of Oregon Dropbox/Alejandro Santillana Fernandez/extinction-cascades/analysis/network")
source('src/initialize.R')
source("src/misc.R")
library(dplyr)
library(ggplot2)
library(tidyr)
library(broom.mixed)  # for tidy() and augment() with lme4 models
library(lme4)
load('saved/mods/metrics.Rdata')


## ---- By Year ----
# ys <- c("FunRedundancy.Pols",
#         "FunRedundancy.Plants",
#         "functional.complementarity.HL",
#         "functional.complementarity.LL",
#         "mean.number.of.links.HL",
#         "mean.number.of.links.LL")
# 
# # Ensure Year is numeric for prediction grid
# cor.dats$Year <- as.numeric(as.character(cor.dats$Year))
# 
# # Create prediction data frame
# predict_metric <- function(mod, yname, years) {
#   new_data <- data.frame(Year = years)
#   
#   # Design matrix
#   mm <- model.matrix(~ Year, new_data)
#   fit <- mm %*% fixef(mod)
#   se <- sqrt(diag(mm %*% vcov(mod) %*% t(mm)))
#   
#   new_data %>%
#     mutate(
#       fit = fit,
#       se = se,
#       lower = fit - 1.96 * se,
#       upper = fit + 1.96 * se,
#       metric = yname
#     )
# }
# 
# # Get unique years from data
# unique_years <- sort(unique(cor.dats$Year))
# 
# # Run prediction for each model
# predictions <- lapply(ys, function(y) {
#   predict_metric(mods.year[[y]], y, unique_years)
# })
# 
# pred_df <- bind_rows(predictions)
# 
# 
# # Define your labels
# metric_labels <- c(
#   "FunRedundancy.Pols" = "Pollinator Redundancy",
#   "FunRedundancy.Plants" = "Plant Redundancy",
#   "functional.complementarity.HL" = "Pollinator Complementarity",
#   "functional.complementarity.LL" = "Plant Complementarity",
#   "mean.number.of.links.HL" = "Pollinator Generalization",
#   "mean.number.of.links.LL" = "Plant Generalization"
# )
# 
# ###  Visualization
# p <- ggplot(pred_df, aes(x = Year, y = fit)) +
#   geom_line(color = "black") +
#   geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80", alpha = 0.4) +
#   geom_point(data = cor.dats %>%
#                pivot_longer(cols = all_of(ys), names_to = "metric", values_to = "value"),
#              aes(x = Year, y = value),
#              inherit.aes = FALSE,
#              alpha = 0.4, shape = 1) +
#   facet_wrap(~ metric, scales = "free_y", labeller = labeller(metric = metric_labels)) +
#   theme_minimal(base_size = 14) +
#   labs(x = "Year", y = "Predicted value") +
#   theme(strip.text = element_text(size = 12)) +
#   scale_x_continuous(
#     breaks = seq(floor(min(pred_df$Year)), ceiling(max(pred_df$Year)), by = 3),
#     labels = function(x) as.integer(x)
#   )
# 
# # Save the plot as PDF
# ggsave("figures/NetworkMetricsByYear.pdf", plot = p, width = 10, height = 7)
# 
# # Optional: display the plot in R console
# print(p)

## ---- By Antecedent Monsoon Precipitation ----

## Join Monsoon Precip into cor.dats
# Shift monsoon year forward so 2011 precip matches 2012 network data
monsoon_precip_data_corrected <- monsoon_precip_data %>%
  mutate(Year = Year + 1)  # Shift year forward by one

# Make sure Site and Year are the same type in both data frames
monsoon_precip_data_corrected$Year <- as.numeric(monsoon_precip_data_corrected$Year)
cor.dats$Year <- as.numeric(cor.dats$Year)

# Merge monsoon precipitation into main dataset
cor.dats <- left_join(cor.dats, monsoon_precip_data_corrected, by = c("Site", "Year"))

# Define the metrics to analyze
ys <- c("FunRedundancy.Pols",
        "FunRedundancy.Plants",
        "functional.complementarity.HL",
        "functional.complementarity.LL",
        "mean.number.of.links.HL",
        "mean.number.of.links.LL")

## Fit models: metric ~ Mean_Monsoon_Precip
mods.monsoon <- lapply(ys, function(y) {
  formula <- as.formula(paste(y, "~ Mean_Monsoon_Precip + (1 | Site)"))
  lmer(formula, data = cor.dats, REML = FALSE)  
})
summary(mods.monsoon[[1]])
names(mods.monsoon) <- ys

## Prediction function using monsoon precipitation
predict_metric_precip <- function(mod, yname, precip_values) {
  new_data <- data.frame(Mean_Monsoon_Precip = precip_values)
  
  mm <- model.matrix(~ Mean_Monsoon_Precip, new_data)
  fit <- mm %*% fixef(mod)
  se <- sqrt(diag(mm %*% vcov(mod) %*% t(mm)))
  
  new_data %>%
    mutate(
      fit = fit,
      se = se,
      lower = fit - 1.96 * se,
      upper = fit + 1.96 * se,
      metric = yname
    )
}

##  Generate predictions across a reasonable precip range 
precip_range <- seq(min(cor.dats$Mean_Monsoon_Precip, na.rm = TRUE),
                    max(cor.dats$Mean_Monsoon_Precip, na.rm = TRUE),
                    length.out = 100)

predictions_precip <- lapply(ys, function(y) {
  predict_metric_precip(mods.monsoon[[y]], y, precip_range)
})
pred_precip_df <- bind_rows(predictions_precip)

## Plot: Metric ~ Monsoon Precip
metric_labels <- c(
  "FunRedundancy.Pols" = "Pollinator Redundancy",
  "FunRedundancy.Plants" = "Plant Redundancy",
  "functional.complementarity.HL" = "Pollinator Complementarity",
  "functional.complementarity.LL" = "Plant Complementarity",
  "mean.number.of.links.HL" = "Pollinator Generalization",
  "mean.number.of.links.LL" = "Plant Generalization"
)

p <- ggplot(pred_precip_df, aes(x = Mean_Monsoon_Precip, y = fit)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80", alpha = 0.4) +
  geom_point(
    data = cor.dats %>%
      pivot_longer(cols = all_of(ys), names_to = "metric", values_to = "value"),
    aes(x = Mean_Monsoon_Precip, y = value, color = Site, shape = as.factor(Year)),
    inherit.aes = FALSE,
    size = 2,
    alpha = 0.6
  ) +
  facet_wrap(~ metric, scales = "free_y", labeller = labeller(metric = metric_labels)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Monsoon Precipitation (mm)",
    y = "Predicted value",
    color = "Site",
    shape = "Year"
  ) +
  theme(strip.text = element_text(size = 12))

ggsave("figures/NetworkMetricsByMonsoonPrecip.pdf", plot = p, width = 10, height = 7)
print(p)

## ---- By Winter Precipitation ----

winter_precip_data <- winter_precip_data %>%
  rename(Year = Winter_Year)

# Ensure types match for join
winter_precip_data$Year <- as.numeric(winter_precip_data$Year)
cor.dats$Year <- as.numeric(cor.dats$Year)

# Join winter precip data to main dataset
cor.dats <- left_join(cor.dats, winter_precip_data, by = c("Site", "Year"))

# Define metrics to model
ys <- c("FunRedundancy.Pols",
        "FunRedundancy.Plants",
        "functional.complementarity.HL",
        "functional.complementarity.LL",
        "mean.number.of.links.HL",
        "mean.number.of.links.LL")

### Fit models: metric ~ Mean_Winter_Precip + (1 | Site)
mods.winter <- lapply(ys, function(y) {
  formula <- as.formula(paste(y, "~ Mean_Winter_Precip + (1 | Site)"))
  lmer(formula, data = cor.dats, REML = FALSE)
})
names(mods.winter) <- ys

### Prediction function for winter precip
predict_metric_winter <- function(mod, yname, precip_values) {
  new_data <- data.frame(Mean_Winter_Precip = precip_values)
  
  mm <- model.matrix(~ Mean_Winter_Precip, new_data)
  fit <- mm %*% fixef(mod)
  se <- sqrt(diag(mm %*% vcov(mod) %*% t(mm)))
  
  new_data %>%
    mutate(
      fit = fit,
      se = se,
      lower = fit - 1.96 * se,
      upper = fit + 1.96 * se,
      metric = yname
    )
}

### Generate predictions
precip_range_winter <- seq(min(cor.dats$Mean_Winter_Precip, na.rm = TRUE),
                           max(cor.dats$Mean_Winter_Precip, na.rm = TRUE),
                           length.out = 100)

predictions_winter <- lapply(ys, function(y) {
  predict_metric_winter(mods.winter[[y]], y, precip_range_winter)
})
pred_winter_df <- bind_rows(predictions_winter)

### Plot
metric_labels <- c(
  "FunRedundancy.Pols" = "Pollinator Redundancy",
  "FunRedundancy.Plants" = "Plant Redundancy",
  "functional.complementarity.HL" = "Pollinator Complementarity",
  "functional.complementarity.LL" = "Plant Complementarity",
  "mean.number.of.links.HL" = "Pollinator Generalization",
  "mean.number.of.links.LL" = "Plant Generalization"
)

p_winter <- ggplot(pred_winter_df, aes(x = Mean_Winter_Precip, y = fit)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80", alpha = 0.4) +
  geom_point(
    data = cor.dats %>%
      pivot_longer(cols = all_of(ys), names_to = "metric", values_to = "value"),
    aes(x = Mean_Winter_Precip, y = value, color = Site, shape = as.factor(Year)),
    inherit.aes = FALSE,
    size = 2,
    alpha = 0.6
  ) +
  facet_wrap(~ metric, scales = "free_y", labeller = labeller(metric = metric_labels)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Winter Precipitation (mm)",
    y = "Predicted value",
    color = "Site",
    shape = "Year"
  ) +
  theme(strip.text = element_text(size = 12))

ggsave("analysis/network/figures/NetworkMetricsByWinterPrecip.pdf", plot = p_winter, width = 10, height = 7)
print(p_winter)

## ---- By Degree Days ----
annual_gdd <- annual_gdd %>% rename(Site = site_name)
annual_gdd <- annual_gdd %>% rename(Year = year)

# Merge with your cor.dats
cor.dats <- left_join(cor.dats, annual_gdd, by = c("Site", "Year"))

ys <- c("FunRedundancy.Pols",
        "FunRedundancy.Plants",
        "functional.complementarity.HL",
        "functional.complementarity.LL",
        "mean.number.of.links.HL",
        "mean.number.of.links.LL")

# Model
mods.annual_gdd <- lapply(ys, function(y) {
  formula <- as.formula(paste(y, "~ total_gdd + (1 | Site)"))
  lmer(formula, data = cor.dats, REML = FALSE)
})
names(mods.annual_gdd) <- ys

# Prediction function
predict_metric_gdd <- function(mod, yname, gdd_values) {
  new_data <- data.frame(total_gdd = gdd_values)
  
  # Construct model matrix for fixed effects only
  mm <- model.matrix(~ total_gdd, new_data)
  
  # Get fixed effect estimates
  fit <- mm %*% fixef(mod)
  
  # Calculate standard errors
  se <- sqrt(diag(mm %*% vcov(mod) %*% t(mm)))
  
  # Return a dataframe with fit and confidence intervals
  new_data %>%
    mutate(
      fit = as.vector(fit),
      se = se,
      lower = fit - 1.96 * se,
      upper = fit + 1.96 * se,
      metric = yname
    )
}

# Generate prediction data over a reasonable GDD range
gdd_range <- seq(min(cor.dats$total_gdd, na.rm = TRUE),
                 max(cor.dats$total_gdd, na.rm = TRUE),
                 length.out = 100)

predictions_gdd <- lapply(ys, function(y) {
  predict_metric_gdd(mods.annual_gdd[[y]], y, gdd_range)
})

pred_gdd_df <- bind_rows(predictions_gdd)

# Plotting
metric_labels <- c(
  "FunRedundancy.Pols" = "Pollinator Redundancy",
  "FunRedundancy.Plants" = "Plant Redundancy",
  "functional.complementarity.HL" = "Pollinator Complementarity",
  "functional.complementarity.LL" = "Plant Complementarity",
  "mean.number.of.links.HL" = "Pollinator Generalization",
  "mean.number.of.links.LL" = "Plant Generalization"
)

p_gdd <- ggplot(pred_gdd_df, aes(x = total_gdd, y = fit)) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80", alpha = 0.4) +
  geom_point(
    data = cor.dats %>%
      pivot_longer(cols = all_of(ys), names_to = "metric", values_to = "value"),
    aes(x = total_gdd, y = value, color = Site, shape = as.factor(Year)),
    inherit.aes = FALSE,
    size = 2,
    alpha = 0.6
  ) +
  facet_wrap(~ metric, scales = "free_y", labeller = labeller(metric = metric_labels)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Total Growing Degree Days",
    y = "Predicted Network Metric",
    color = "Site",
    shape = "Year"
  ) +
  theme(strip.text = element_text(size = 12))

# Save and print plot
ggsave("analysis/network/figures/NetworkMetricsByTotalGDD.pdf", plot = p_gdd, width = 10, height = 7)
print(p_gdd)

















































# ### Coextinction simulations based on abundance ###
# 
# # Create an empty data frame to store the results
# robustness_results <- data.frame(Network = character(), Site = character(), Year = integer(), Robustness = numeric(), stringsAsFactors = FALSE)
# 
# # Loop through each network in the list `nets`
# for (net_name in names(nets)) {
#   # Extract the network
#   network <- nets[[net_name]]
#   
#   # Run extinction simulation
#   network_ext <- second.extinct(network, participant = "higher", method = "abun", nrep = 50, details = FALSE)
#   
#   # Calculate the robustness (area under second.extinct curve)
#   network_robustness <- robustness(network_ext)
#   
#   # Extract Site and Year from the network name
#   site_year <- unlist(strsplit(net_name, "\\.")) # Split by "."
#   site <- site_year[1]
#   year <- as.integer(site_year[2])
#   
#   # Append the results to the data frame
#   robustness_results <- rbind(robustness_results, data.frame(Network = net_name, Site = site, Year = year, Robustness = network_robustness))
# }
# 
# # Print the resulting data frame
# print(robustness_results)

















# 
# # second.extinct calculates the consequences of removing a species from a bipartite network. With plotting and slope-estimation functionality.
# # 
# # Details: The procedure of this function is simple. For example imagine the web to represent a pollination web, in which pollinators die one by one. Set all entries of a column to zero, see how may rows are now also all-zero (i.e. species that are now not pollinated any more), and count these as secondary extinctions of the primary exterminated pollinator.
# # 
# # Internally, each extermination is achieved by a call to extinction, followed by a call to empty, which counts the number of all-zero columns and rows.
# # 
# # Although written for pollination webs (hence the nomenclature of pollinator and plant), it can be similarly applied to other types of bipartite networks. It is called by networklevel.
# # Usage: second.extinct(web, participant = "higher", method = "abun", nrep = 10, details = FALSE, ext.row=NULL, ext.col=NULL)
# # 
# # Arguments:
# #   method: Random deletion of a species (random); according to its abundance, with least abundant going extinct first (abundance; default) or "degree" to build a sequence from the best-to-least connected species. This is the most extreme case, where the most generalist species goes extinct first (see Memmott et al. 1998). (Note that method="abundance" does not work with participant="both"; in this case a random sequence of species from both trophic levels is generated.) external will use the externally provided vector to determine extinction sequence.
# # ext.row: Optional vector giving the sequence in which lower-level species are to be deleted.
# # ext.col: Optional vector giving the sequence in which higher-level species are to be deleted.
# # 
# # Note: networklevel calls second.extinct; extinction and empty are internal helper functions, and slope.bipartite calculates extinction slopes from the output of second.extinct.
# # 
# # ```{r second.extinct}
# # # ?second.extinct
# # # ?slope.bipartite
# # # ?robustness
# # ```
# # 
# # bipartite function "networklevel" gives outputs:
# #   $‘extinction slope lower trophic level‘ [1] 2.286152 
# # $‘extinction slope higher trophic level‘ [1] 2.9211
# # 
# # Extinction slopes are hyperbolic fits to a simulated extinction sequence of the network, which causes secondary extinctions in the other trophic level (from Memmott 2004).
# # 
# # slope.bipartite fits a hyperbolic function to the extinction simulation of "second.extinct."
# # This function scales extinction sequences to values between 0 and 1 for each participant. 
# # The x-axis of the graph features the proportion of exterminated participants, while the y-axis depicts the proportion of secondary extinctions. Since these curves usually follow a hyperbolic function (see examples in Memmott et al. 2004), this is fitted to the data.
# # 
# 
# 
# 
# 
# 
# 
# #skyIslands/data/networks/Yr_PlantPollinator_Bees.RData #yearly networks for bees
# #Lili knows that sample rounds are not equal across all years
# #Years with less rounds would look a lot less diverse
# #view(nets$SM.2018)
# 
# 
# load("../skyIslands/data/networks/Year_PlantPollinator_Bees.RData") #yearly networks for bees
# networklevel(nets[["CH.2018"]])
# plotweb(nets[["CH.2018"]], method="normal", text.rot=90)
# 
# ## ***********************************************************
# ## sankey plots
# ## ***********************************************************
# 
# makeSankeyPlot<- function(i){
#   print(i)
#   print(names(nets)[i])
#   data_long <- as.data.frame(nets[[i]]) %>%
#     rownames_to_column %>%
#     gather(key = 'key', value = 'value', -rowname) %>%
#     filter(value > 0)
#   colnames(data_long) <- c("source", "target", "value")
#   data_long$target <- paste(data_long$target, " ", sep="")
#   
#   ## From these flows we need to create a node data frame: it lists
#   ## every entities involved in the flow
#   nodes <- data.frame(name=c(as.character(data_long$source),
#                              as.character(data_long$target)) %>% unique())
#   
#   ## With networkD3, connection must be provided using id, not using
#   ## real name like in the links dataframe.. So we need to reformat
#   ## it.
#   data_long$IDsource=match(data_long$source, nodes$name)-1 
#   data_long$IDtarget=match(data_long$target, nodes$name)-1  
#   net.p <- sankeyNetwork(Links = data_long, Nodes = nodes,
#                          Source = "IDsource", Target = "IDtarget",
#                          Value = "value", NodeID = "name", 
#                          fontSize=18)
#   ## colourScale=rainbow(
#   ##   prod(dim(bee.nets[[i]]))))
#   saveNetwork(net.p, sprintf("figures/sankey_%s.html",
#                              gsub(" ", "", names(nets)[i])))
#   webshot(sprintf("figures/sankey_%s.html",
#                   gsub(" ", "", names(nets)[i])),
#           file=sprintf("figures/sankey_%s.jpeg",
#                        gsub(" ", "",
#                             names(nets)[i])),
#           vwidth = 650, vheight = 1400)
#   
# }
# 
# lapply(1:length(nets), makeSankeyPlot)
# 
# # Visualize networks using sankey plots
# makeSankeyPlot<- function(i){
#   print(i)
#   print(names(nets)[i])
#   data_long <- as.data.frame(nets[[i]]) %>%
#     rownames_to_column %>%
#     gather(key = 'key', value = 'value', -rowname) %>%
#     filter(value > 0)
#   colnames(data_long) <- c("source", "target", "value")
#   data_long$target <- paste(data_long$target, " ", sep="")
#   
#   ## From these flows we need to create a node data frame: it lists
#   ## every entities involved in the flow
#   nodes <- data.frame(name=c(as.character(data_long$source),
#                              as.character(data_long$target)) %>% unique())
#   
#   ## With networkD3, connection must be provided using id, not using
#   ## real name like in the links dataframe.. So we need to reformat
#   ## it.
#   data_long$IDsource=match(data_long$source, nodes$name)-1
#   data_long$IDtarget=match(data_long$target, nodes$name)-1
#   net.p <- sankeyNetwork(Links = data_long, Nodes = nodes,
#                          Source = "IDsource", Target = "IDtarget",
#                          Value = "value", NodeID = "name",
#                          fontSize=18)
#   ## colourScale=rainbow(
#   ##   prod(dim(bee.nets[[i]]))))
#   saveNetwork(net.p, sprintf("figures/sankey_%s.html",
#                              gsub(" ", "", names(nets)[i])))
#   webshot(sprintf("figures/sankey_%s.html",
#                   gsub(" ", "", names(nets)[i])),
#           file=sprintf("figures/sankey_%s.jpeg",
#                        gsub(" ", "",
#                             names(nets)[i])),
#           vwidth = 650, vheight = 1400)
#   
# }
# lapply(1:length(nets), makeSankeyPlot)
# 
# 
# 
# # # Define a list of the rasters
# # rasters <- list(
# #   precip_2011 = precip_2011_raster,
# #   precip_2016 = precip_2016_raster,
# #   precip_2017 = precip_2017_raster,
# #   precip_2020 = precip_2020_raster,
# #   precip_2021 = precip_2021_raster
# # )
# # 
# # # Define a list of corresponding years for labeling
# # years <- c("2011", "2016", "2017", "2020", "2021")
# # 
# # # Calculate global minimum and maximum precipitation values across all rasters
# # min_precip <- min(sapply(rasters, function(r) min(values(r), na.rm = TRUE)))
# # max_precip <- max(sapply(rasters, function(r) max(values(r), na.rm = TRUE)))
# # 
# # # Loop through each raster and create a plot for each year
# # for (i in 1:length(rasters)) {
# #   
# #   # Get the current raster and year
# #   current_raster <- rasters[[i]]
# #   current_year <- years[i]
# #   
# #   # Convert the current raster to a data.frame for ggplot
# #   precip_df <- as.data.frame(current_raster, xy = TRUE, na.rm = TRUE)
# #   
# #   # Rename the third column to "precip"
# #   colnames(precip_df)[3] <- "precip"
# #   
# #   # Apply log transformation to precipitation values (log(x + 1) to avoid log(0))
# #   precip_df$precip_log <- log(precip_df$precip + 1)
# #   
# #   # Create the plot for the current year
# #   plot <- ggplot() +
# #     # Add the precipitation raster layer
# #     geom_raster(data = precip_df, aes(x = x, y = y, fill = precip_log)) +
# #     
# #     # Add a color scale for the precipitation
# #    scale_fill_viridis(option = "D", name = paste("Precipitation", year, "(log scale)"), direction = 1, 
# #                       limits = c(log(min_precip + 1), log(max_precip + 1))) +
# #     
# #     # Add the site locations
# #     geom_sf(data = sites, color = "red", size = 3, shape = 21, fill = "yellow") +
# #     
# #     # Customize map extent (adjust as needed)
# #     coord_sf(xlim = c(-111, -105), ylim = c(31, 37), expand = FALSE) +
# #     
# #     # Customize the appearance
# #     theme_minimal() +
# #     labs(
# #       title = paste("Sky Island Sites and Precipitation in", current_year, "(log-transformed)"),
# #       subtitle = paste("Precipitation layer for the year", current_year),
# #       x = "Longitude", 
# #       y = "Latitude"
# #     ) +
# #     theme(legend.position = "right")
# #   
# #   # Print the plot for the current year
# #   print(plot)
# # }
# 
# 
# # # Plot Robustness vs Mean_Precipitation by site
# # ggplot(final_data, aes(x = Mean_Precipitation, y = Robustness)) +
# #   geom_point(aes(color = Site), size = 3) +       # Points colored by Site
# #   geom_smooth(method = "lm", se = TRUE, color = "blue") + # Add linear regression line
# #   theme_minimal() +
# #   labs(title = "Robustness vs. Summer Precipitation",
# #        x = "Mean Precipitation (previous summer) (mm)",
# #        y = "Network Robustness")
# 
# # # Plot Robustness vs. Mean_Precipitation by year
# # ggplot(final_data, aes(x = Mean_Precipitation, y = Robustness, color = as.factor(Year))) +
# #   geom_point(size = 3) +
# #   geom_smooth(method = "lm", se = FALSE) +
# #   labs(title = "Network Robustness vs. Summer Precipitation by Year",
# #        x = "Mean Precipitation (previous summer) (mm)",
# #        y = "Network Robustness",
# #        color = "Year") +
# #   theme_minimal()
# 
# 
# # networklevel(nets[["CH.2018"]])
# # 
# # trial <- second.extinct(nets[["CH.2018"]], participant = "higher", method = "abun", nrep = 500, details = FALSE)
# # trial_slope <- slope.bipartite(trial, plot.it = TRUE)
# 
# 
# 
# 
