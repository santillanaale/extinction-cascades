rm(list=ls())
library(bipartite)
library(fossil)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(networkD3)


# second.extinct calculates the consequences of removing a species from a bipartite network. With plotting and slope-estimation functionality.
# 
# Details: The procedure of this function is simple. For example imagine the web to represent a pollination web, in which pollinators die one by one. Set all entries of a column to zero, see how may rows are now also all-zero (i.e. species that are now not pollinated any more), and count these as secondary extinctions of the primary exterminated pollinator.
# 
# Internally, each extermination is achieved by a call to extinction, followed by a call to empty, which counts the number of all-zero columns and rows.
# 
# Although written for pollination webs (hence the nomenclature of pollinator and plant), it can be similarly applied to other types of bipartite networks. It is called by networklevel.
# Usage: second.extinct(web, participant = "higher", method = "abun", nrep = 10, details = FALSE, ext.row=NULL, ext.col=NULL)
# 
# Arguments:
#   method: Random deletion of a species (random); according to its abundance, with least abundant going extinct first (abundance; default) or "degree" to build a sequence from the best-to-least connected species. This is the most extreme case, where the most generalist species goes extinct first (see Memmott et al. 1998). (Note that method="abundance" does not work with participant="both"; in this case a random sequence of species from both trophic levels is generated.) external will use the externally provided vector to determine extinction sequence.
# ext.row: Optional vector giving the sequence in which lower-level species are to be deleted.
# ext.col: Optional vector giving the sequence in which higher-level species are to be deleted.
# 
# Note: networklevel calls second.extinct; extinction and empty are internal helper functions, and slope.bipartite calculates extinction slopes from the output of second.extinct.
# 
# ```{r second.extinct}
# # ?second.extinct
# # ?slope.bipartite
# # ?robustness
# ```
# 
# bipartite function "networklevel" gives outputs:
#   $‘extinction slope lower trophic level‘ [1] 2.286152 
# $‘extinction slope higher trophic level‘ [1] 2.9211
# 
# Extinction slopes are hyperbolic fits to a simulated extinction sequence of the network, which causes secondary extinctions in the other trophic level (from Memmott 2004).
# 
# slope.bipartite fits a hyperbolic function to the extinction simulation of "second.extinct."
# This function scales extinction sequences to values between 0 and 1 for each participant. 
# The x-axis of the graph features the proportion of exterminated participants, while the y-axis depicts the proportion of secondary extinctions. Since these curves usually follow a hyperbolic function (see examples in Memmott et al. 2004), this is fitted to the data.
# 






#skyIslands/data/networks/Yr_PlantPollinator_Bees.RData #yearly networks for bees
#Lili knows that sample rounds are not equal across all years
#Years with less rounds would look a lot less diverse
#view(nets$SM.2018)


load("../skyIslands/data/networks/Yr_PlantPollinator_Bees.RData") #yearly networks for bees
networklevel(nets[["CH.2018"]])
plotweb(nets[["CH.2018"]], method="normal", text.rot=90)

## ***********************************************************
## sankey plots
## ***********************************************************

makeSankeyPlot<- function(i){
  print(i)
  print(names(nets)[i])
  data_long <- as.data.frame(nets[[i]]) %>%
    rownames_to_column %>%
    gather(key = 'key', value = 'value', -rowname) %>%
    filter(value > 0)
  colnames(data_long) <- c("source", "target", "value")
  data_long$target <- paste(data_long$target, " ", sep="")
  
  ## From these flows we need to create a node data frame: it lists
  ## every entities involved in the flow
  nodes <- data.frame(name=c(as.character(data_long$source),
                             as.character(data_long$target)) %>% unique())
  
  ## With networkD3, connection must be provided using id, not using
  ## real name like in the links dataframe.. So we need to reformat
  ## it.
  data_long$IDsource=match(data_long$source, nodes$name)-1 
  data_long$IDtarget=match(data_long$target, nodes$name)-1  
  net.p <- sankeyNetwork(Links = data_long, Nodes = nodes,
                         Source = "IDsource", Target = "IDtarget",
                         Value = "value", NodeID = "name", 
                         fontSize=18)
  ## colourScale=rainbow(
  ##   prod(dim(bee.nets[[i]]))))
  saveNetwork(net.p, sprintf("figures/sankey_%s.html",
                             gsub(" ", "", names(nets)[i])))
  webshot(sprintf("figures/sankey_%s.html",
                  gsub(" ", "", names(nets)[i])),
          file=sprintf("figures/sankey_%s.jpeg",
                       gsub(" ", "",
                            names(nets)[i])),
          vwidth = 650, vheight = 1400)
  
}

lapply(1:length(nets), makeSankeyPlot)

# Visualize networks using sankey plots
makeSankeyPlot<- function(i){
  print(i)
  print(names(nets)[i])
  data_long <- as.data.frame(nets[[i]]) %>%
    rownames_to_column %>%
    gather(key = 'key', value = 'value', -rowname) %>%
    filter(value > 0)
  colnames(data_long) <- c("source", "target", "value")
  data_long$target <- paste(data_long$target, " ", sep="")
  
  ## From these flows we need to create a node data frame: it lists
  ## every entities involved in the flow
  nodes <- data.frame(name=c(as.character(data_long$source),
                             as.character(data_long$target)) %>% unique())
  
  ## With networkD3, connection must be provided using id, not using
  ## real name like in the links dataframe.. So we need to reformat
  ## it.
  data_long$IDsource=match(data_long$source, nodes$name)-1
  data_long$IDtarget=match(data_long$target, nodes$name)-1
  net.p <- sankeyNetwork(Links = data_long, Nodes = nodes,
                         Source = "IDsource", Target = "IDtarget",
                         Value = "value", NodeID = "name",
                         fontSize=18)
  ## colourScale=rainbow(
  ##   prod(dim(bee.nets[[i]]))))
  saveNetwork(net.p, sprintf("figures/sankey_%s.html",
                             gsub(" ", "", names(nets)[i])))
  webshot(sprintf("figures/sankey_%s.html",
                  gsub(" ", "", names(nets)[i])),
          file=sprintf("figures/sankey_%s.jpeg",
                       gsub(" ", "",
                            names(nets)[i])),
          vwidth = 650, vheight = 1400)
  
}
lapply(1:length(nets), makeSankeyPlot)



# # Define a list of the rasters
# rasters <- list(
#   precip_2011 = precip_2011_raster,
#   precip_2016 = precip_2016_raster,
#   precip_2017 = precip_2017_raster,
#   precip_2020 = precip_2020_raster,
#   precip_2021 = precip_2021_raster
# )
# 
# # Define a list of corresponding years for labeling
# years <- c("2011", "2016", "2017", "2020", "2021")
# 
# # Calculate global minimum and maximum precipitation values across all rasters
# min_precip <- min(sapply(rasters, function(r) min(values(r), na.rm = TRUE)))
# max_precip <- max(sapply(rasters, function(r) max(values(r), na.rm = TRUE)))
# 
# # Loop through each raster and create a plot for each year
# for (i in 1:length(rasters)) {
#   
#   # Get the current raster and year
#   current_raster <- rasters[[i]]
#   current_year <- years[i]
#   
#   # Convert the current raster to a data.frame for ggplot
#   precip_df <- as.data.frame(current_raster, xy = TRUE, na.rm = TRUE)
#   
#   # Rename the third column to "precip"
#   colnames(precip_df)[3] <- "precip"
#   
#   # Apply log transformation to precipitation values (log(x + 1) to avoid log(0))
#   precip_df$precip_log <- log(precip_df$precip + 1)
#   
#   # Create the plot for the current year
#   plot <- ggplot() +
#     # Add the precipitation raster layer
#     geom_raster(data = precip_df, aes(x = x, y = y, fill = precip_log)) +
#     
#     # Add a color scale for the precipitation
#    scale_fill_viridis(option = "D", name = paste("Precipitation", year, "(log scale)"), direction = 1, 
#                       limits = c(log(min_precip + 1), log(max_precip + 1))) +
#     
#     # Add the site locations
#     geom_sf(data = sites, color = "red", size = 3, shape = 21, fill = "yellow") +
#     
#     # Customize map extent (adjust as needed)
#     coord_sf(xlim = c(-111, -105), ylim = c(31, 37), expand = FALSE) +
#     
#     # Customize the appearance
#     theme_minimal() +
#     labs(
#       title = paste("Sky Island Sites and Precipitation in", current_year, "(log-transformed)"),
#       subtitle = paste("Precipitation layer for the year", current_year),
#       x = "Longitude", 
#       y = "Latitude"
#     ) +
#     theme(legend.position = "right")
#   
#   # Print the plot for the current year
#   print(plot)
# }


# # Plot Robustness vs Mean_Precipitation by site
# ggplot(final_data, aes(x = Mean_Precipitation, y = Robustness)) +
#   geom_point(aes(color = Site), size = 3) +       # Points colored by Site
#   geom_smooth(method = "lm", se = TRUE, color = "blue") + # Add linear regression line
#   theme_minimal() +
#   labs(title = "Robustness vs. Summer Precipitation",
#        x = "Mean Precipitation (previous summer) (mm)",
#        y = "Network Robustness")

# # Plot Robustness vs. Mean_Precipitation by year
# ggplot(final_data, aes(x = Mean_Precipitation, y = Robustness, color = as.factor(Year))) +
#   geom_point(size = 3) +
#   geom_smooth(method = "lm", se = FALSE) +
#   labs(title = "Network Robustness vs. Summer Precipitation by Year",
#        x = "Mean Precipitation (previous summer) (mm)",
#        y = "Network Robustness",
#        color = "Year") +
#   theme_minimal()


# networklevel(nets[["CH.2018"]])
# 
# trial <- second.extinct(nets[["CH.2018"]], participant = "higher", method = "abun", nrep = 500, details = FALSE)
# trial_slope <- slope.bipartite(trial, plot.it = TRUE)




