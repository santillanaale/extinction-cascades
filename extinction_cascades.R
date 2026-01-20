rm(list=ls())
library(bipartite)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(lme4)
set.seed(123)
setwd("~/")
source("lab_paths.R")
local.path

# ---- Load Networks ----
setwd(file.path(local.path, "skyIslands"))
load("data/networks/YearSR_PlantPollinator_Bees.Rdata") #Site x Year x Survey Round

## --- Make net key ----
net_meta <- tibble(net_id = names(nets)) %>%
  separate(net_id, into = c("Site", "Year", "SampleRound"), sep = "\\.") %>%
  mutate(
    Year = as.numeric(Year),
    SampleRound = as.numeric(SampleRound)
  )

# ---- TCM: Combining all scenarios into one loop ----
## ---- Load modified extinction functions ----
setwd("~/")
setwd(file.path(local.path, "extinction-cascades"))

source("modified_extinction.R")
source("modified_second_extinct.R")

## ---- Unified helper function ----
run_robustness <- function(web, scenario, nrep = 100) {
  
  stopifnot(scenario$method %in% c("abundance", "degree"))
  
  ext <- modified.second.extinct(
    web = web,
    participant = "higher",      # bees
    method = scenario$method,
    reverse = scenario$reverse,
    nrep = nrep,
    details = FALSE
  )
  
  robustness(ext)
}

## ---- Define scenarios ----
scenarios <- tibble(
  scenario = c(
    "abundance_low", # remove rare first
    "abundance_high", # remove dominant first
    "degree_low", # remove specialists first
    "degree_high" # remove generalists first
  ),
  method = c(
    "abundance",
    "abundance",
    "degree",
    "degree"
  ),
  reverse = c(
    FALSE,
    TRUE,
    FALSE,
    TRUE
  )
)

## ---- Robustness Loop (Site × Year × SampleRound × scenario × Robustness----
robustness_results <- map_dfr(seq_len(nrow(scenarios)), function(s) {
  
  scenario <- scenarios[s, ]
  
  message("Running scenario: ", scenario$scenario)
  
  map_dfr(seq_along(nets), function(i) {
    
    web <- nets[[i]]
    
    rob <- run_robustness(
      web = web,
      scenario = scenario,
      nrep = 100
    )
    
    tibble(
      Site = net_meta$Site[i],
      Year = net_meta$Year[i],
      SampleRound = net_meta$SampleRound[i],
      scenario = scenario$scenario,
      Robustness = rob
    )
  })
})

# ---- Merge Climate ----
setwd("~/")
setwd(file.path(local.path, "skyIslands_saved"))

climate <- read.csv("data/relational/original/climate.csv")

robust_climate <- robustness_results %>%
  left_join(
    climate,
    by = c("Site", "Year", "SampleRound")
  )

setwd("~/")
setwd(file.path(local.path, "extinction-cascades"))
save(
  robustness_results,
  robust_climate,
  file = "analysis/network/robustness_metrics/YearSR_robustness_climate.Rdata"
)

# ---- Modeling Robustness vs Climate ----
## ---- Standardize predictors ----
robust_climate <- robust_climate %>%
  mutate(across(
    c(SpringPrecip, CumulativePrecip, RoundPrecip,
      SpringTmean, CumulativeTmean, RoundTmean, 
      SpringTmeanAnom, RoundTmeanAnom, CumulativeTmeanAnom),
    ~ scale(.)[,1],
    .names = "z_{.col}"
  ))

## ---- Core models ----
### ---- Precipitation only ----
mod_precip <- lmer(
  Robustness ~ z_SpringPrecip +
    z_CumulativePrecip +
    z_RoundPrecip +
    scenario +
    (1 | Site) +
    (1 | Year),
  data = robust_climate
)

summary(mod_precip)

### ---- Temperature only (absolute, TMean) ----
mod_tmean <- lmer(
  Robustness ~ z_SpringTmean +
    z_CumulativeTmean +
    z_RoundTmean +
    scenario +
    (1 | Site) +
    (1 | Year),
  data = robust_climate
)

summary(mod_tmean)

### ---- Temperature anomalies ----
mod_tmean_anom <- lmer(
  Robustness ~ z_SpringTmeanAnom +
    z_CumulativeTmeanAnom +
    z_RoundTmeanAnom +
    scenario +
    (1 | Site) +
    (1 | Year),
  data = robust_climate
)

summary(mod_tmean_anom)

## ---- Full Climate Model ----
mod_all <- lmer(
  Robustness ~ z_SpringPrecip +
    z_CumulativePrecip +
    z_RoundPrecip +
    z_CumulativeTmeanAnom +
    scenario +
    (1 | Site) +
    (1 | Year),
  data = robust_climate
)

summary(mod_all)

# Interpretation framework (this is paper-ready)
# 
# Spring terms → legacy / phenological setup
# 
# Cumulative terms → seasonal resource tracking
# 
# Round terms → short-term sensitivity
# 
# Scenario effects → robustness mechanism dependence
# 
# Anomalies vs absolute → physiological vs climatological drivers


# ---- Plotting robustness vs climate ----
## ---- TCM ----
### ---- Precipitation ----
p1 <- ggplot(robust_climate,
             aes(x = SpringPrecip, y = Robustness)) +
  geom_point(aes(color = Site), alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Spring Precipitation",
       x = "Spring precipitation",
       y = "Robustness")

p2 <- ggplot(robust_climate,
             aes(x = CumulativePrecip, y = Robustness)) +
  geom_point(aes(color = Site), alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Cumulative Precipitation",
       x = "Cumulative precipitation")

p3 <- ggplot(robust_climate,
             aes(x = RoundPrecip, y = Robustness)) +
  geom_point(aes(color = Site), alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Round Precipitation",
       x = "Round precipitation")

(p1 | p2 | p3)

ggsave(
  "analysis/network/figures/Robustness_by_Precip_RoundScale.pdf",
  width = 12, height = 4
)

### ---- Mean Temperature (absolute) ----
p1 <- ggplot(robust_climate,
             aes(x = SpringTmean, y = Robustness)) +
  geom_point(aes(color = Site), alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Spring Mean Temperature")

p2 <- ggplot(robust_climate,
             aes(x = CumulativeTmean, y = Robustness)) +
  geom_point(aes(color = Site), alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Cumulative Mean Temperature")

p3 <- ggplot(robust_climate,
             aes(x = RoundTmean, y = Robustness)) +
  geom_point(aes(color = Site), alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Round Mean Temperature")

(p1 | p2 | p3)

ggsave(
  "analysis/network/figures/Robustness_by_Temperature_RoundScale.pdf",
  width = 12, height = 4
)

### ---- Mean Temperature (Anomaly) ----
p1 <- ggplot(robust_climate,
             aes(x = SpringTmeanAnom, y = Robustness)) +
  geom_point(aes(color = Site), alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Spring Mean Temperature")

p2 <- ggplot(robust_climate,
             aes(x = CumulativeTmeanAnom, y = Robustness)) +
  geom_point(aes(color = Site), alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Cumulative Mean Temperature")

p3 <- ggplot(robust_climate,
             aes(x = RoundTmeanAnom, y = Robustness)) +
  geom_point(aes(color = Site), alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Round Mean Temperature")

(p1 | p2 | p3)

ggsave(
  "analysis/network/figures/Robustness_by_TemperatureAnomaly_RoundScale.pdf",
  width = 12, height = 4
)

### ---- By extinction scenario to show whether extinction order changes climate sensitivity ----
#### ---- Precip ----
ggplot(robust_climate,
       aes(x = SpringPrecip, y = Robustness)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ scenario) +
  theme_minimal()

ggplot(robust_climate,
       aes(x = CumulativePrecip, y = Robustness)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ scenario) +
  theme_minimal()

ggplot(robust_climate,
       aes(x = RoundPrecip, y = Robustness)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ scenario) +
  theme_minimal()

#### ---- Temp ----
ggplot(robust_climate,
       aes(x = SpringTmean, y = Robustness)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ scenario) +
  theme_minimal()

ggplot(robust_climate,
       aes(x = CumulativeTmean, y = Robustness)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ scenario) +
  theme_minimal()

ggplot(robust_climate,
       aes(x = RoundTmean, y = Robustness)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ scenario) +
  theme_minimal()

#### ---- TempAnom ----
ggplot(robust_climate,
       aes(x = SpringTmeanAnom, y = Robustness)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ scenario) +
  theme_minimal()

ggplot(robust_climate,
       aes(x = CumulativeTmeanAnom, y = Robustness)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ scenario) +
  theme_minimal()

ggplot(robust_climate,
       aes(x = RoundTmeanAnom, y = Robustness)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ scenario) +
  theme_minimal()



## ---- SCM ----
source("SCM/netcascade_mod.R")
source("SCM/scm_robustness_helpers.R")

### ---- Define SCM scenarios (parallel to TCM) ----
scm_scenarios <- tibble(
  scenario = c(
    "scm_random_plant",
    "scm_dominant_plant"
  ),
  TargetGuild = c("rows", "rows"),
  TargetSpecies = c("random_binary", "random_abundance")
)

### ---- Run SCM robustness in the same Site × Year x Sample Round loop ----
scm_results <- map_dfr(seq_len(nrow(scm_scenarios)), function(s) {
  
  scenario <- scm_scenarios[s, ]
  
  message("Running SCM scenario: ", scenario$scenario)
  
  map_dfr(seq_along(nets), function(i) {
    
    web <- as.matrix(nets[[i]])
    
    R_rows <- rep(0.5, nrow(web))  # or sample like VA2015
    R_cols <- rep(0.5, ncol(web))
    
    out <- run_scm_robustness(
      web = web,
      R_val_rows = R_rows,
      R_val_cols = R_cols,
      TargetGuild = scenario$TargetGuild,
      TargetSpecies = scenario$TargetSpecies,
      nrep = 100
    )
    
    tibble(
      Site = net_meta$Site[i],
      Year = net_meta$Year[i],
      SampleRound = net_meta$SampleRound[i],
      scenario = scenario$scenario,
      prop_cascade = out$prop_cascade,
      SCM_Robustness = out$SCM_robustness
    )
  })
})

### ---- Merge with Climate ----
scm_climate <- scm_results %>%
  left_join(climate, by = c("Site", "Year", "SampleRound"))










##########################################################################################################################

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


#### ---- Plotting some scenarios ----
# 
# # Example network
# example_network <- nets[["CH.2012"]]
# 
# # Extract extinction trajectories from the matrix output
# specialist_df <- data.frame(
#   step = specialist_ext[, "no"],  # number of pollinators removed
#   plants_remaining = 1 - (specialist_ext[, "ext.lower"] / max(specialist_ext[, "ext.lower"])),
#   scenario = "Specialist-first"
# )
# 
# generalist_df <- data.frame(
#   step = generalist_ext[, "no"],
#   plants_remaining = 1 - (generalist_ext[, "ext.lower"] / max(generalist_ext[, "ext.lower"])),
#   scenario = "Generalist-first"
# )
# 
# # Combine and plot
# extinction_df <- rbind(specialist_df, generalist_df)
# 
# ggplot(extinction_df, aes(x = step, y = plants_remaining, color = scenario)) +
#   geom_line(linewidth = 1.2) +
#   scale_color_manual(values = c("Specialist-first" = "#2196F3", "Generalist-first" = "#E64A19")) +
#   labs(
#     title = "Simulated Co-extinction Cascades in a Pollination Network",
#     x = "Number of Pollinators Removed",
#     y = "Fraction of Plants Remaining",
#     color = "Extinction Order"
#   ) +
#   theme_minimal(base_size = 14)

# trial <- second.extinct(nets[["CH.2018"]], participant = "higher", method = "abun", nrep = 500, details = FALSE)
# trial_slope <- slope.bipartite(trial, plot.it = TRUE)

#### ---- Abundance ----
# Run extinction simulation (abundance-based removal) (least to most abundant removed)
# abun <- second.extinct(example_network,
#                         participant = "higher",
#                         method = "abun",
#                         nrep = 500,
#                         details = FALSE)
# 
# n_plants <- nrow(example_network)
# 
# # Summarize extinction curves (mean across reps)
# # Each column is one replicate, rows correspond to % pollinators removed
# abun_mean <- rowMeans(abun[, -1])  # exclude 'no' column if present
# steps <- seq(0, 1, length.out = length(abun_mean))  # fraction of pollinators removed
# 
# # Build data frame for ggplot
# extinction_df <- data.frame(
#   step = steps,
#   plants_remaining = 1 - (abun_mean / n_plants),  # fraction of plants remaining
#   scenario = "Abundance-based"
# )
# 
# # Optional: Get robustness value (area under the curve)
# robustness_value <- robustness(abun)
# 
# # Plot with ggplot2
# ggplot(extinction_df, aes(x = step, y = plants_remaining, color = scenario)) +
#   geom_line(linewidth = 1.3) +
#   scale_color_manual(values = c("Abundance-based" = "#2196F3")) +
#   labs(
#     title = "Simulated Co-extinction Cascade",
#     subtitle = paste("Robustness =", round(robustness_value, 3)),
#     x = "Fraction of Pollinators Removed",
#     y = "Fraction of Plants Remaining",
#     color = "Extinction Order"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "bottom")
# Load modified extinction functions
# source("modified_extinction.R")
# source("modified_second_extinct.R")
# 
# 
# abun <- second.extinct(nets[["CH.2018"]], participant = "higher", method = "abun", nrep = 500, details = FALSE)
# slope.bipartite(abun, plot.it = TRUE)

#### ---- Reverse Abundance ----
# # Run extinction simulation (abundance-based removal) (most to least abundant removed)
# reverse_abun <- modified.second.extinct(example_network,
#                        participant = "higher",
#                        method = "abun",
#                        nrep = 500,
#                        details = FALSE,
#                        reverse = TRUE)
# 
# # Summarize extinction curves (mean across reps)
# # Each column is one replicate, rows correspond to % pollinators removed
# reverse_abun_mean <- rowMeans(reverse_abun[, -1])  # exclude 'no' column if present
# steps <- seq(0, 1, length.out = length(reverse_abun_mean))  # fraction of pollinators removed
# 
# # Build data frame for ggplot
# extinction_df <- data.frame(
#   step = steps,
#   plants_remaining = 1 - (reverse_abun_mean / n_plants),  # fraction of plants remaining
#   scenario = "Reverse Abundance-based"
# )
# 
# # Optional: Get robustness value (area under the curve)
# robustness_value <- robustness(reverse_abun)
# 
# # Plot with ggplot2
# ggplot(extinction_df, aes(x = step, y = plants_remaining, color = scenario)) +
#   geom_line(linewidth = 1.3) +
#   scale_color_manual(values = c("Reverse Abundance-based" = "#2196F3")) +
#   labs(
#     title = "Simulated Co-extinction Cascade",
#     subtitle = paste("Robustness =", round(robustness_value, 3)),
#     x = "Fraction of Pollinators Removed",
#     y = "Fraction of Plants Remaining",
#     color = "Extinction Order"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "bottom")

# reverse_abun <- modified.second.extinct(nets[["CH.2018"]], participant = "higher", method = "abun", nrep = 500, reverse = TRUE, details = FALSE)
# slope.bipartite(reverse_abun, plot.it = TRUE)
# 
# 
#### ---- Degree ----
# # Run extinction simulation (abundance-based removal) (most connected to least connected)
# degree_ext <- second.extinct(nets[["CH.2018"]], participant = "higher", method = "degree", nrep = 500, details = FALSE)
# slope.bipartite(degree_ext, plot.it = TRUE)
# 
#### ---- Reverse Abundance ----
# # Run extinction simulation (abundance-based removal) (least connected to most connected)
# reverse_degree_ext <- modified.second.extinct(nets[["CH.2018"]], participant = "higher", method = "degree", nrep = 500, details = FALSE, reverse = TRUE)
# slope.bipartite(reverse_degree_ext, plot.it = TRUE)



# # Use only modified.second.extinct
# scenarios <- list(
#   abun = list(method = "abun", reverse = FALSE),
#   abun_reverse = list(method = "abun", reverse = TRUE),
#   degree = list(method = "degree", reverse = FALSE),
#   degree_reverse = list(method = "degree", reverse = TRUE)
# )
# 
# 
# for (scenario_name in names(scenarios)) {
#   scenario <- scenarios[[scenario_name]]
#   
#   results_list <- list()
#   
#   for (net_name in names(nets)) {
#     network <- nets[[net_name]]
#     
#     # Run extinction simulation using modified function for all
#     ext_result <- modified.second.extinct(
#       web = network,
#       participant = "higher",
#       method = scenario$method,
#       reverse = scenario$reverse,
#       nrep = 50,
#       details = FALSE
#     )
#     
#     rob <- robustness(ext_result)
#     
#     site_year <- unlist(strsplit(net_name, "\\."))
#     site <- site_year[1]
#     year <- as.integer(site_year[2])
#     
#     results_list[[net_name]] <- data.frame(
#       Network = net_name,
#       Site = site,
#       Year = year,
#       Robustness = rob,
#       stringsAsFactors = FALSE
#     )
#   }
#   
#   scenario_df <- do.call(rbind, results_list)
#   
#   # Save results
#   save_path <- sprintf("analysis/network/robustness_metrics/%s_robustness_metrics.Rdata", scenario_name)
#   save(scenario_df, file = save_path)
#   
#   assign(paste0("TCM_", scenario_name, "_robustness_results"), scenario_df)
#   print(scenario_df)
# }

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

# ## ---- Plotting Robustness ----
# # Choose the robustness data (e.g., for 'abun' scenario)
# robustness_df <- TCM_abun_robustness_results
# 
# # Correct for antecedent monsoon year
# monsoon_precip_data_df <- monsoon_precip_data %>%
#   mutate(Year = Year + 1)  # Shift year forward by one
# 
# # Change column names for winter precip
# winter_precip_data_df <- winter_precip_data %>%
#   rename(Year = Winter_Year)
# 
# # Change column names for GDD
# annual_gdd_df <- annual_gdd %>%
#   rename(Site = site_name, Year = year)
# 
# # Merge all climate data into robustness_df
# robustness_climate <- robustness_df %>%
#   left_join(monsoon_precip_data_df, by = c("Site", "Year")) %>%
#   left_join(winter_precip_data_df, by = c("Site", "Year")) %>%
#   left_join(annual_gdd_df, by = c("Site", "Year"))
# 
# # Plot
# p1 <- ggplot(robustness_climate, aes(x = Mean_Monsoon_Precip, y = Robustness)) +
#   geom_point(aes(color = Site)) + geom_smooth(method = "lm") + theme_minimal() +
#   labs(title = "Monsoon Precip")
# 
# p2 <- ggplot(robustness_climate, aes(x = Mean_Winter_Precip, y = Robustness)) +
#   geom_point(aes(color = Site)) + geom_smooth(method = "lm") + theme_minimal() +
#   labs(title = "Winter Precip")
# 
# p3 <- ggplot(robustness_climate, aes(x = total_gdd, y = Robustness)) +
#   geom_point(aes(color = Site)) + geom_smooth(method = "lm") + theme_minimal() +
#   labs(title = "GDD")
# 
# # Combine into 1 row
# p1 + p2 + p3 + plot_layout(ncol = 3)
# ggsave("analysis/network/figures/RobustnessbyClimate.pdf", width = 12, height = 4)
# 
# ## ---- Modeling Robustness ----
# # Standardize predictors
# robustness_climate <- robustness_climate %>%
#   mutate(across(c(Mean_Monsoon_Precip, Mean_Winter_Precip, total_gdd),
#                 ~ scale(.)[,1],  # extract the numeric column from scale()
#                 .names = "z_{.col}"))
# # 
# # # Individual models
# # mod_monsoon <- lmer(Robustness ~ z_Mean_Monsoon_Precip + (1 | Site), data = robustness_climate)
# # mod_winter <- lmer(Robustness ~ z_Mean_Winter_Precip + (1 | Site), data = robustness_climate)
# # mod_gdd    <- lmer(Robustness ~ z_total_gdd + (1 | Site), data = robustness_climate)
# 
# # Combined model
# mod_all <- lmer(Robustness ~ z_Mean_Monsoon_Precip + z_Mean_Winter_Precip + z_total_gdd + (1 | Site),
#                 data = robustness_climate)
# summary(mod_all)
# 
# # Get p-values
# model_summary <- summary(mod_all)$coefficients
# 
# # Identify significance of each fixed effect (p < 0.05)
# sig_monsoon <- model_summary["z_Mean_Monsoon_Precip", "Pr(>|t|)"] < 0.05
# sig_winter  <- model_summary["z_Mean_Winter_Precip", "Pr(>|t|)"] < 0.05
# sig_gdd     <- model_summary["z_total_gdd", "Pr(>|t|)"] < 0.05
# 
# # Predict effects for each variable
# monsoon_eff <- ggpredict(mod_all, terms = "z_Mean_Monsoon_Precip")
# winter_eff  <- ggpredict(mod_all, terms = "z_Mean_Winter_Precip")
# gdd_eff     <- ggpredict(mod_all, terms = "z_total_gdd")
# 
# # Plot each effect separately (clean & labeled)
# plot_effect <- function(pred_df, xlab, title, significant = FALSE) {
#   line_type <- ifelse(significant, "solid", "dashed")
#   
#   ggplot(pred_df, aes(x = x, y = predicted)) +
#     geom_line(size = 1, linetype = line_type) +
#     geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
#     labs(x = xlab, y = "Predicted Robustness", title = title) +
#     theme_minimal(base_size = 14)
# }
# 
# p1 <- plot_effect(monsoon_eff, "Standardized Monsoon Precipitation", "Effect of Monsoon Precipitation", significant = sig_monsoon)
# p2 <- plot_effect(winter_eff,  "Standardized Winter Precipitation", "Effect of Winter Precipitation", significant = sig_winter)
# p3 <- plot_effect(gdd_eff,     "Standardized Growing Degree Days (GDD)", "Effect of GDD", significant = sig_gdd)
# 
# combined_plot <- p1 + p2 + p3 + plot_layout(ncol = 1)
# combined_plot
# ggsave("analysis/network/figures/Robustness_ClimateEffects_Combined.pdf", plot = combined_plot, width = 8, height = 10)
# 

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

library(broom.mixed)

# Extract tidy results for all models
monsoon_results <- lapply(names(mods.monsoon), function(name) {
  broom.mixed::tidy(mods.monsoon[[name]], effects = "fixed") %>%
    mutate(Metric = name)
}) %>%
  bind_rows()

# Keep only the climate predictor rows
monsoon_results_clean <- monsoon_results %>%
  filter(term == "Mean_Monsoon_Precip") %>%
  select(Metric, estimate, std.error, statistic, p.value)

# Print a neat summary table
print(monsoon_results_clean)

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

### Filter predictions and data for Plant Generalization
pred_monsoon_LL <- pred_precip_df %>%
  filter(metric == "mean.number.of.links.LL")

data_LL <- cor.dats %>%
  select(Site, Year, Mean_Monsoon_Precip, mean.number.of.links.LL) %>%
  rename(value = mean.number.of.links.LL)

# Extract from the combined results table
coef_LL <- monsoon_results_clean %>%
  filter(Metric == "mean.number.of.links.LL")

slope_LL <- coef_LL$estimate
p_LL <- coef_LL$p.value

label_text <- paste0(
  "Slope = ", round(slope_LL, 3),
  "\nP = ", signif(p_LL, 2)
)

### Plot only Plant Generalization
p_monsoon_LL <- ggplot(pred_monsoon_LL, aes(x = Mean_Monsoon_Precip, y = fit)) +
  geom_line(color = "black", linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "gray80", alpha = 0.4) +
  geom_point(data = data_LL,
             aes(x = Mean_Monsoon_Precip, y = value, color = Site, shape = as.factor(Year)),
             size = 2, alpha = 0.6) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Effect of Monsoon Precipitation on Plant Generalization",
    x = "Monsoon Precipitation (mm)",
    y = "Mean number of links (Plant Generalization)",
    color = "Site",
    shape = "Year"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 12)
  ) +
  annotate("label",
           x = max(pred_monsoon_LL$Mean_Monsoon_Precip) - 
             0.05 * diff(range(pred_monsoon_LL$Mean_Monsoon_Precip)),
           y = max(pred_monsoon_LL$fit) - 
             0.05 * diff(range(pred_monsoon_LL$fit)),
           label = label_text,
           hjust = 1, vjust = 1,
           size = 5,
           fill = "white",
           color = "black",
           alpha = 0.8)

### Print and save plot
print(p_monsoon_LL)

ggsave("figures/PlantGeneralizationByMonsoonPrecip.pdf", plot = p_monsoon_LL, width = 10, height = 7)

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

# Extract tidy results for all models
winter_results <- lapply(names(mods.winter), function(name) {
  broom.mixed::tidy(mods.winter[[name]], effects = "fixed") %>%
    mutate(Metric = name)
}) %>%
  bind_rows()

# Keep only the climate predictor rows
winter_results_clean <- winter_results %>%
  filter(term == "Mean_Winter_Precip") %>%
  select(Metric, estimate, std.error, statistic, p.value)

# Print a neat summary table
print(winter_results_clean)

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

ggsave("figures/NetworkMetricsByWinterPrecip.pdf", plot = p_winter, width = 10, height = 7)
print(p_winter)



### Filter predictions and data for Plant Generalization
pred_winter_LL <- pred_winter_df %>%
  filter(metric == "mean.number.of.links.LL")

data_LL <- cor.dats %>%
  select(Site, Year, Mean_Winter_Precip, mean.number.of.links.LL) %>%
  rename(value = mean.number.of.links.LL)

# Extract from the combined results table
coef_LL <- winter_results_clean %>%
  filter(Metric == "mean.number.of.links.LL")

slope_LL <- coef_LL$estimate
p_LL <- coef_LL$p.value

label_text <- paste0(
  "Slope = ", round(slope_LL, 3),
  "\nP = ", signif(p_LL, 2)
)

### Plot only Plant Generalization
p_winter_LL <- ggplot(pred_winter_LL, aes(x = Mean_Winter_Precip, y = fit)) +
  geom_line(color = "black", linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "gray80", alpha = 0.4) +
  geom_point(data = data_LL,
             aes(x = Mean_Winter_Precip, y = value, color = Site, shape = as.factor(Year)),
             size = 2, alpha = 0.6) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Effect of Winter Precipitation on Plant Generalization",
    x = "Winter Precipitation (mm)",
    y = "Mean number of links (Plant Generalization)",
    color = "Site",
    shape = "Year"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 12)
  )

p_winter_LL <- p_winter_LL +
  annotate("label",
           x = min(pred_winter_LL$Mean_Winter_Precip) + 
               0.05 * diff(range(pred_winter_LL$Mean_Winter_Precip)),
           y = max(pred_winter_LL$fit) - 
               0.05 * diff(range(pred_winter_LL$fit)),
           label = label_text,
           hjust = 0, vjust = 1,
           size = 5,
           fill = "white",
           color = "black",
           alpha = 0.8)


### Print plot
p_winter_LL
ggsave("figures/PlantGeneralizationByWinterPrecip.pdf", plot = p_winter_LL, width = 10, height = 7)

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

# Extract tidy results for all models
gdd_results <- lapply(names(mods.annual_gdd), function(name) {
  broom.mixed::tidy(mods.annual_gdd[[name]], effects = "fixed") %>%
    mutate(Metric = name)
}) %>%
  bind_rows()

# Keep only the climate predictor rows
gdd_results_clean <- gdd_results %>%
  filter(term == "mean_annual_gdd") %>%
  select(Metric, estimate, std.error, statistic, p.value)

# Print a neat summary table
print(gdd_results_clean)

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
ggsave("figures/NetworkMetricsByTotalGDD.pdf", plot = p_gdd, width = 10, height = 7)
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
