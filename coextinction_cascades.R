rm(list=ls())
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

# ---- Import Precipitation Rasters ----
### Import monthly precipitation data for each field season ###

precip_2011_06_raster <- rast("PRISM_summer_precip/2011/PRISM_ppt_stable_4kmM3_201106_bil.bil")
precip_2011_07_raster <- rast("PRISM_summer_precip/2011/PRISM_ppt_stable_4kmM3_201107_bil.bil")
precip_2011_08_raster <- rast("PRISM_summer_precip/2011/PRISM_ppt_stable_4kmM3_201108_bil.bil")
precip_2011_09_raster <- rast("PRISM_summer_precip/2011/PRISM_ppt_stable_4kmM3_201109_bil.bil")

precip_2016_06_raster <- rast("PRISM_summer_precip/2016/PRISM_ppt_stable_4kmM3_201606_bil.bil")
precip_2016_07_raster <- rast("PRISM_summer_precip/2016/PRISM_ppt_stable_4kmM3_201607_bil.bil")
precip_2016_08_raster <- rast("PRISM_summer_precip/2016/PRISM_ppt_stable_4kmM3_201608_bil.bil")
precip_2016_09_raster <- rast("PRISM_summer_precip/2016/PRISM_ppt_stable_4kmM3_201609_bil.bil")

precip_2017_06_raster <- rast("PRISM_summer_precip/2017/PRISM_ppt_stable_4kmM3_201706_bil.bil")
precip_2017_07_raster <- rast("PRISM_summer_precip/2017/PRISM_ppt_stable_4kmM3_201707_bil.bil")
precip_2017_08_raster <- rast("PRISM_summer_precip/2017/PRISM_ppt_stable_4kmM3_201708_bil.bil")
precip_2017_09_raster <- rast("PRISM_summer_precip/2017/PRISM_ppt_stable_4kmM3_201709_bil.bil")

precip_2020_06_raster <- rast("PRISM_summer_precip/2020/PRISM_ppt_stable_4kmM3_202006_bil.bil")
precip_2020_07_raster <- rast("PRISM_summer_precip/2020/PRISM_ppt_stable_4kmM3_202007_bil.bil")
precip_2020_08_raster <- rast("PRISM_summer_precip/2020/PRISM_ppt_stable_4kmM3_202008_bil.bil")
precip_2020_09_raster <- rast("PRISM_summer_precip/2020/PRISM_ppt_stable_4kmM3_202009_bil.bil")

precip_2021_06_raster <- rast("PRISM_summer_precip/2021/PRISM_ppt_stable_4kmM3_202106_bil.bil")
precip_2021_07_raster <- rast("PRISM_summer_precip/2021/PRISM_ppt_stable_4kmM3_202107_bil.bil")
precip_2021_08_raster <- rast("PRISM_summer_precip/2021/PRISM_ppt_stable_4kmM3_202108_bil.bil")
precip_2021_09_raster <- rast("PRISM_summer_precip/2021/PRISM_ppt_stable_4kmM3_202109_bil.bil")

# Define the years and file paths
years <- c("2011", "2016", "2017", "2020", "2021")
months <- c("06", "07", "08", "09")
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

precip_2011_raster <- summer_rasters[["2011"]]
precip_2016_raster <- summer_rasters[["2016"]]
precip_2017_raster <- summer_rasters[["2017"]]
precip_2020_raster <- summer_rasters[["2020"]]
precip_2021_raster <- summer_rasters[["2021"]]

# Identify variables with the pattern "precip_****_**_raster"
monthly_rasters <- ls(pattern = "^precip_\\d{4}_\\d{2}_raster$")
rm(list = monthly_rasters)


# ---- Import and Transform Spatial Data ----

#The following script imports a shapefile for study sites, ensuring consistency in coordinate reference systems (CRS) with PRISM rasters.
sites_shp <- vect("Sites/sites.shp")
print(crs(sites_shp))
print(crs(precip_2011_raster))
sites_shp <- project(sites_shp, crs(precip_2011_raster))



# ---- Geospatial Visualization ----

#Sky Islands Map
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



# ---- Import and Visualize (Sankey) Networks ----

# A Sankey plot is a type of flow diagram that visualizes the connections between two sets of elements, often used to show weighted interactions. In this context, the Sankey plot for the CH.2012 network represents a plant-pollinator interaction network, illustrating how different pollinators interacted with various plant species in the Chiricahua Mountains site in the survey year 2012.
load("../skyIslands/data/networks/Year_PlantPollinator_Bees.RData") #yearly networks for bees

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


# ---- Extract Precipitation Data for Sites ----

# This code extracts precipitation data from raster files for specified years and calculates average precipitation for each study site.
years <- c(2011, 2016, 2017, 2020, 2021)
results_list <- list()

for (year in years) {
  # Load the raster for the year
  raster <- rast(paste0("PRISM_precip/PRISM_ppt_stable_4kmM3_", year, "_bil.bil")) 
  
  # Extract precipitation values for subsites
  precip_values <- extract(raster, sites_shp)
  
  # Combine extracted values with site shapefile data
  subsite_data <- as.data.frame(sites_shp) # Convert shapefile to a data frame
  subsite_data$Precipitation <- precip_values[, 2] # Add extracted values
  subsite_data$Year <- year # Add the year
  
  # Group by Site and Year, then calculate mean precipitation
  site_data <- subsite_data %>%
    group_by(Site, Year) %>%
    summarize(Mean_Precipitation = mean(Precipitation, na.rm = TRUE), .groups = "drop")
  
  # Store the yearly averaged data in the results list
  results_list[[as.character(year)]] <- site_data
}

# Combine all yearly site-level results into a single data frame
all_precip_data <- do.call(rbind, results_list)


# ---- Perform Extinction Simulations and Calculate Robustness ----

# This code simulates species extinctions within plant-pollinator networks, calculates their robustness (a measure of how well the network resists species loss), and 
# organizes the results into a data frame for analysis.


## ---- TCM Topological Coextinction Model ----

### ----  Increasing abundance  ----
# Create an empty data frame to store the results
TCM_increasing_abundance_robustness_results <- data.frame(Network = character(), Site = character(), Year = integer(), Robustness = numeric(), stringsAsFactors = FALSE)

# Loop through each network in the list `nets`
for (net_name in names(nets)) {
  # Extract the network
  network <- nets[[net_name]]
  
  # Run extinction simulation
  network_ext <- second.extinct(network, participant = "higher", method = "abun", nrep = 50, details = FALSE)
  
  # Calculate the robustness (area under second.extinct curve)
  network_robustness <- robustness(network_ext)
  
  # Extract Site and Year from the network name
  site_year <- unlist(strsplit(net_name, "\\.")) # Split by "."
  site <- site_year[1]
  year <- as.integer(site_year[2])
  
  # Append the results to the data frame
  TCM_increasing_abundance_robustness_results <- rbind(TCM_increasing_abundance_robustness_results, data.frame(Network = net_name, Site = site, Year = year, Robustness = network_robustness))
}

# Print the resulting data frame
print(TCM_increasing_abundance_robustness_results)

### ----  Decreasing abundance  ----
# Create an empty data frame to store the results
TCM_decreasing_abundance_robustness_results <- data.frame(Network = character(), Site = character(), Year = integer(), Robustness = numeric(), stringsAsFactors = FALSE)

source("modified_extinction.R")
source("modified_second_extinct.R")

# Sanity Check
modified.second.extinct(web = example_network, participant = "higher", method = "abundance", reverse = TRUE,     # 👈 removes most abundant species instead of least
  nrep = 1, details = FALSE)

# Loop through each network in the list `nets`
for (net_name in names(nets)) {
  # Extract the network
  network <- nets[[net_name]]

  # Run extinction simulation
  network_ext <- modified.second.extinct(network, participant = "higher", method = "abun", reverse = TRUE, nrep = 50, details = FALSE)

  # Calculate the robustness (area under second.extinct curve)
  network_robustness <- robustness(network_ext)

  # Extract Site and Year from the network name
  site_year <- unlist(strsplit(net_name, "\\.")) # Split by "."
  site <- site_year[1]
  year <- as.integer(site_year[2])

  # Append the results to the data frame
  TCM_decreasing_abundance_robustness_results <- rbind(decreasing_abundance_robustness_results, data.frame(Network = net_name, Site = site, Year = year, Robustness = network_robustness))
}

# Print the resulting data frame
print(TCM_decreasing_abundance_robustness_results)

### ----  Decreasing Degree  ----
# Create an empty data frame to store the results
TCM_decreasing_degree_robustness_results <- data.frame(Network = character(), Site = character(), Year = integer(), Robustness = numeric(), stringsAsFactors = FALSE)

# Loop through each network in the list `nets`
for (net_name in names(nets)) {
  # Extract the network
  network <- nets[[net_name]]
  
  # Run extinction simulation
  network_ext <- second.extinct(network, participant = "higher", method = "degree", nrep = 50, details = FALSE)
  
  # Calculate the robustness (area under second.extinct curve)
  network_robustness <- robustness(network_ext)
  
  # Extract Site and Year from the network name
  site_year <- unlist(strsplit(net_name, "\\.")) # Split by "."
  site <- site_year[1]
  year <- as.integer(site_year[2])
  
  # Append the results to the data frame
  TCM_decreasing_degree_robustness_results <- rbind(TCM_decreasing_degree_robustness_results, data.frame(Network = net_name, Site = site, Year = year, Robustness = network_robustness))
}

# Print the resulting data frame
print(TCM_decreasing_degree_robustness_results)

### ----  Increasing Degree  ----
# Create an empty data frame to store the results
TCM_increasing_degree_robustness_results <- data.frame(Network = character(), Site = character(), Year = integer(), Robustness = numeric(), stringsAsFactors = FALSE)

# Loop through each network in the list `nets`
for (net_name in names(nets)) {
  # Extract the network
  network <- nets[[net_name]]
  
  # Run extinction simulation
  network_ext <- modified.second.extinct(network, participant = "higher", method = "degree", reverse = TRUE, nrep = 50, details = FALSE)
  
  # Calculate the robustness (area under second.extinct curve)
  network_robustness <- robustness(network_ext)
  
  # Extract Site and Year from the network name
  site_year <- unlist(strsplit(net_name, "\\.")) # Split by "."
  site <- site_year[1]
  year <- as.integer(site_year[2])
  
  # Append the results to the data frame
  TCM_increasing_degree_robustness_results <- rbind(TCM_increasing_degree_robustness_results, data.frame(Network = net_name, Site = site, Year = year, Robustness = network_robustness))
}

# Print the resulting data frame
print(TCM_decreasing_abundance_robustness_results)

## ---- SCM: Stochastic Coextinction Model ----

#source('SCM/VieraAlmeida2015RFunctions/netcascade (April 2014).R')
source("SCM/Dalsgaard2018Code/IterNodeDelMultiSim.R")

# read in R values
rvals <- read.csv("SCM/Dalsgaard2018Code/data/rvalues.csv")
rvals$R <- trimws(rvals$R)

# create table with upper and lower limits for R bounds
r_bounds <- as.data.frame(matrix(ncol = 3, nrow = length(unique(rvals$R)), dimnames = list(NULL,c("R","lb","ub"))))
r_bounds$R <- unique(rvals$R)
r_bounds[r_bounds$R == "low",][,c("lb","ub")] <- c(0,(1/3))
r_bounds[r_bounds$R == "high",][,c("lb","ub")] <- c((2/3),1)
r_bounds[r_bounds$R == "med-high",][,c("lb","ub")] <- c((1/3),1)
r_bounds[r_bounds$R == "low-med",][,c("lb","ub")] <- c(0,(2/3))
r_bounds[r_bounds$R == "med",][,c("lb","ub")] <- c((1/3),(2/3))
r_bounds[r_bounds$R == "low-high",][,c("lb","ub")] <- c(0,1)

# split rvals by web
rvals <- split(rvals, rvals$file_name)

vulnerabilities <- NULL
for(i in 1:length(nets)){
  
  # run simulations
  most <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "bo", removal.taxa = "plant", remove.order = "from most connected", nsim = 10000)
  strongest <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "bo", removal.taxa = "plant", remove.order = "from strongest", nsim = 10000)
  random <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "bo", removal.taxa = "plant", remove.order = "random", nsim = 10000)
  weakest <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "bo", removal.taxa = "plant", remove.order = "from weakest", nsim = 10000)
  least <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "bo", removal.taxa = "plant", remove.order = "from least connected", nsim = 10000)
  
  r.most <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "random", removal.taxa = "plant", remove.order = "from most connected", nsim = 10000)
  r.strongest <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "random", removal.taxa = "plant", remove.order = "from strongest", nsim = 10000)
  r.random <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "random", removal.taxa = "plant", remove.order = "random", nsim = 10000)
  r.weakest <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "random", removal.taxa = "plant", remove.order = "from weakest", nsim = 10000)
  r.least <- IterNodeDelMultiSim(network = nets[[names(nets[i])]], plant_r_assignments = rvals[[names(nets[i])]], r_ints = r_bounds, bo_or_random = "random", removal.taxa = "plant", remove.order = "from least connected", nsim = 10000)
  
  # gather results
  most <- merge(most$mean_order, most$mean_vulnerability, by = "species")
  most$web <- names(nets)[i]
  most$removal.order <- "from most connected"
  most$plantR <- "assigned"
  
  strongest <- merge(strongest$mean_order, strongest$mean_vulnerability, by = "species")
  strongest$web <- names(nets)[i]
  strongest$removal.order <- "from strongest"
  strongest$plantR <- "assigned"
  
  random <- merge(random$mean_order, random$mean_vulnerability, by = "species")
  random$web <- names(nets)[i]
  random$removal.order <- "random"
  random$plantR <- "assigned"
  
  weakest <- merge(weakest$mean_order, weakest$mean_vulnerability, by = "species")
  weakest$web <- names(nets)[i]
  weakest$removal.order <- "from weakest"
  weakest$plantR <- "assigned"
  
  least <- merge(least$mean_order, least$mean_vulnerability, by = "species")
  least$web <- names(nets)[i]
  least$removal.order <- "from least connected"
  least$plantR <- "assigned"
  
  r.most <- merge(r.most$mean_order, r.most$mean_vulnerability, by = "species")
  r.most$web <- names(nets)[i]
  r.most$removal.order <- "from most connected"
  r.most$plantR <- "random"
  
  r.strongest <- merge(r.strongest$mean_order, r.strongest$mean_vulnerability, by = "species")
  r.strongest$web <- names(nets)[i]
  r.strongest$removal.order <- "from strongest"
  r.strongest$plantR <- "random"
  
  r.random <- merge(r.random$mean_order, r.random$mean_vulnerability, by = "species")
  r.random$web <- names(nets)[i]
  r.random$removal.order <- "random"
  r.random$plantR <- "random"
  
  r.weakest <- merge(r.weakest$mean_order, r.weakest$mean_vulnerability, by = "species")
  r.weakest$web <- names(nets)[i]
  r.weakest$removal.order <- "from weakest"
  r.weakest$plantR <- "random"
  
  r.least <- merge(r.least$mean_order, r.least$mean_vulnerability, by = "species")
  r.least$web <- names(nets)[i]
  r.least$removal.order <- "from least connected"
  r.least$plantR <- "random"
  
  vulnerabilities <- rbind(vulnerabilities,rbind(most,strongest,random,weakest,least,r.most,r.strongest,r.random,r.weakest,r.least))
  
  print(paste0(i,"/8"))
}

vulnerabilities$index <- vulnerabilities$probability.of.extinction * (1-vulnerabilities$mean.order)
# Key of how column names map onto vulnerability indices from the paper:
# probability.of.extinction = PE
# mean.order = SE
# index = VE











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
