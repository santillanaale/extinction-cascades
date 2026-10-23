# R version 3.2.3 (2015-12-10)
# bipartite_2.08
# 
# This script:
#   
# Loads multiple hummingbird-plant networks and their species-level R classifications.
# 
# Assigns numerical ranges to R categories (e.g., low, medium, high).
# 
# Runs extinction simulations on plant species in networks using various removal strategies.
# 
# Runs each scenario under both real (bo) and random (random) R-values.
# 
# Merges and stores results for later analysis of species vulnerability patterns.

rm(list = ls())
library(bipartite)
source("IterNodeDelMultiSim.R")

# ---- read in hummingbird networks ----
# Loads all hummingbird-plant networks from a folder and stores each one in a list (web.store) with cleaned file names as keys.
web.store <- list() # Initialize list to hold networks
files <- list.files("./data/webs/", full.names = FALSE) # List all network files (assumed .txt)
for(i in 1:length(files)){
  web <- t(empty(as.matrix(read.table(paste0("./data/webs/",files[i]), header = TRUE, row.names = 1)))) # Read the matrix and transpose it
  web.store[[gsub("\\.txt","",files[i])]] <- web #  Store with name minus ".txt"
}

# ---- read in R values ----
rvals <- read.csv("./data/rvalues.csv") # Load table of R-values per plant
rvals$R <- trimws(rvals$R) # Remove any whitespace from R column

# ---- Define R Category Bounds ----
#Assigns numerical lower and upper bounds for R-value categories

# create table with upper and lower limits for R bounds
r_bounds <- as.data.frame(matrix(ncol = 3, nrow = length(unique(rvals$R)), dimnames = list(NULL,c("R","lb","ub")))) # Create empty bounds table
r_bounds$R <- unique(rvals$R) # Fill in unique R categories

# Set bounds for each R category (e.g., low = [0, 1/3])
r_bounds[r_bounds$R == "low",][,c("lb","ub")] <- c(0,(1/3))
r_bounds[r_bounds$R == "high",][,c("lb","ub")] <- c((2/3),1)
r_bounds[r_bounds$R == "med-high",][,c("lb","ub")] <- c((1/3),1)
r_bounds[r_bounds$R == "low-med",][,c("lb","ub")] <- c(0,(2/3))
r_bounds[r_bounds$R == "med",][,c("lb","ub")] <- c((1/3),(2/3))
r_bounds[r_bounds$R == "low-high",][,c("lb","ub")] <- c(0,1)

# ---- split rvals by web ----
rvals <- split(rvals, rvals$file_name) 

# ---- Run simulations for each network ----
vulnerabilities <- NULL # Initialize output table

# In total: 10 simulations per network: 5 strategies × 2 R-value assignment modes (assigned vs. random)
# Simulation parameters:
  # remove.order: Strategy for plant removal
  # removal.taxa = "plant": Only plant nodes are removed
  # nsim = 10000: Number of simulation iterations

for(i in 1:length(web.store)){ 
  
  # (a) Run Simulations with Assigned R-values (bo)
  most <- IterNodeDelMultiSim(network = web.store[[names(web.store[i])]], plant_r_assignments = rvals[[names(web.store[i])]], r_ints = r_bounds, bo_or_random = "bo", removal.taxa = "plant", remove.order = "from most connected", nsim = 10000)
  strongest <- IterNodeDelMultiSim(network = web.store[[names(web.store[i])]], plant_r_assignments = rvals[[names(web.store[i])]], r_ints = r_bounds, bo_or_random = "bo", removal.taxa = "plant", remove.order = "from strongest", nsim = 10000)
  random <- IterNodeDelMultiSim(network = web.store[[names(web.store[i])]], plant_r_assignments = rvals[[names(web.store[i])]], r_ints = r_bounds, bo_or_random = "bo", removal.taxa = "plant", remove.order = "random", nsim = 10000)
  weakest <- IterNodeDelMultiSim(network = web.store[[names(web.store[i])]], plant_r_assignments = rvals[[names(web.store[i])]], r_ints = r_bounds, bo_or_random = "bo", removal.taxa = "plant", remove.order = "from weakest", nsim = 10000)
  least <- IterNodeDelMultiSim(network = web.store[[names(web.store[i])]], plant_r_assignments = rvals[[names(web.store[i])]], r_ints = r_bounds, bo_or_random = "bo", removal.taxa = "plant", remove.order = "from least connected", nsim = 10000)
  
  # (b) Run Simulations with Random R-values
  r.most <- IterNodeDelMultiSim(network = web.store[[names(web.store[i])]], plant_r_assignments = rvals[[names(web.store[i])]], r_ints = r_bounds, bo_or_random = "random", removal.taxa = "plant", remove.order = "from most connected", nsim = 10000)
  r.strongest <- IterNodeDelMultiSim(network = web.store[[names(web.store[i])]], plant_r_assignments = rvals[[names(web.store[i])]], r_ints = r_bounds, bo_or_random = "random", removal.taxa = "plant", remove.order = "from strongest", nsim = 10000)
  r.random <- IterNodeDelMultiSim(network = web.store[[names(web.store[i])]], plant_r_assignments = rvals[[names(web.store[i])]], r_ints = r_bounds, bo_or_random = "random", removal.taxa = "plant", remove.order = "random", nsim = 10000)
  r.weakest <- IterNodeDelMultiSim(network = web.store[[names(web.store[i])]], plant_r_assignments = rvals[[names(web.store[i])]], r_ints = r_bounds, bo_or_random = "random", removal.taxa = "plant", remove.order = "from weakest", nsim = 10000)
  r.least <- IterNodeDelMultiSim(network = web.store[[names(web.store[i])]], plant_r_assignments = rvals[[names(web.store[i])]], r_ints = r_bounds, bo_or_random = "random", removal.taxa = "plant", remove.order = "from least connected", nsim = 10000)
  
  ## ---- Results Extraction & Annotation ----
  # For each simulation result, two tables are merged: mean_order: 
    # average extinction order for species, 
    # mean_vulnerability: probability of extinction
  most <- merge(most$mean_order, most$mean_vulnerability, by = "species")
  most$web <- names(web.store)[i]
  most$removal.order <- "from most connected"
  most$plantR <- "assigned"
  
  strongest <- merge(strongest$mean_order, strongest$mean_vulnerability, by = "species")
  strongest$web <- names(web.store)[i]
  strongest$removal.order <- "from strongest"
  strongest$plantR <- "assigned"
  
  random <- merge(random$mean_order, random$mean_vulnerability, by = "species")
  random$web <- names(web.store)[i]
  random$removal.order <- "random"
  random$plantR <- "assigned"
  
  weakest <- merge(weakest$mean_order, weakest$mean_vulnerability, by = "species")
  weakest$web <- names(web.store)[i]
  weakest$removal.order <- "from weakest"
  weakest$plantR <- "assigned"
  
  least <- merge(least$mean_order, least$mean_vulnerability, by = "species")
  least$web <- names(web.store)[i]
  least$removal.order <- "from least connected"
  least$plantR <- "assigned"
  
  r.most <- merge(r.most$mean_order, r.most$mean_vulnerability, by = "species")
  r.most$web <- names(web.store)[i]
  r.most$removal.order <- "from most connected"
  r.most$plantR <- "random"
  
  r.strongest <- merge(r.strongest$mean_order, r.strongest$mean_vulnerability, by = "species")
  r.strongest$web <- names(web.store)[i]
  r.strongest$removal.order <- "from strongest"
  r.strongest$plantR <- "random"
  
  r.random <- merge(r.random$mean_order, r.random$mean_vulnerability, by = "species")
  r.random$web <- names(web.store)[i]
  r.random$removal.order <- "random"
  r.random$plantR <- "random"
  
  r.weakest <- merge(r.weakest$mean_order, r.weakest$mean_vulnerability, by = "species")
  r.weakest$web <- names(web.store)[i]
  r.weakest$removal.order <- "from weakest"
  r.weakest$plantR <- "random"
  
  r.least <- merge(r.least$mean_order, r.least$mean_vulnerability, by = "species")
  r.least$web <- names(web.store)[i]
  r.least$removal.order <- "from least connected"
  r.least$plantR <- "random"
  
  ## --- Aggregate Results ----
  # All results from the 10 simulations are appended into a master dataframe vulnerabilities.
  vulnerabilities <- rbind(vulnerabilities,rbind(most,strongest,random,weakest,least,r.most,r.strongest,r.random,r.weakest,r.least))
  
  print(paste0(i,"/8")) # Logging progress, assumes 8 networks in total
}

# ---- Final Step: Calculate Composite Vulnerability Index ----
  # probability.of.extinction = PE
  # mean.order = SE (smaller = earlier extinction)
  # index = PE × (1 - SE) = vulnerability index (VE)
#This emphasizes species that go extinct early and with high probability.

vulnerabilities$index <- vulnerabilities$probability.of.extinction * (1-vulnerabilities$mean.order)



# Key of how column names map onto vulnerability indices from the paper:
# probability.of.extinction = PE
# mean.order = SE
# index = VE