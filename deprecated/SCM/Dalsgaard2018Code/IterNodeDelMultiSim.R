# ---- Function Definition and Arguments ----
# This function takes:
#   network: bipartite interaction matrix (e.g. hummingbirds × plants).
#   plant_r_assignments: a data frame mapping plant species to R-categories (e.g. low, high).
#   r_ints: numerical lower and upper bounds for each R-category.
#   bo_or_random: string, either "bo" (use expert-assigned plant R-values) or "random".
#   removal.taxa: string, which taxa to remove (plants in this case).
#   remove.order: removal strategy (e.g. "from most connected").
#   nsim: number of simulation iterations.

IterNodeDelMultiSim <- function(network,plant_r_assignments,r_ints,bo_or_random,removal.taxa,remove.order,nsim){
  #--------- Loading required packages and functions ---------------------------
  require(bipartite)
  require(data.table)
  source("SCM/Dalsgaard2018Code/IterativeNodeDeletion.R", local = TRUE)
  
  #--------- Defining some variables and initialising some results containers ---------------------------
  sim.store <- list() # Will hold results from all simulations
  nplant <- ncol(network) # Number of pol species (assuming pols = columns)
  
  #--------- Run iterative removal simulations ---------------------------
  # Each iteration simulates one scenario of progressive plant removal and resulting coextinctions.
  for(i in 1:nsim){
    anim.r <- runif(1,(2/3),1) # sample hummingbird value in 'high' range - all hummingbirds given same value
    
    ## ---- Option 1: Use Expert-Assigned R categories ----
    if(bo_or_random == "bo"){ # if using expert-assigned R values for plants
      r_ints$rval <- sapply(seq(nrow(r_ints)), function(x){ # sample one random value for each plant R classification
        runif(1,r_ints[x,"lb"],r_ints[x,"ub"])
      })
      plant_r_assignments.m <- merge(plant_r_assignments, r_ints, by = "R") # merge chosen R values with plant names
      plant_r_assignments.m <- plant_r_assignments.m[order(match(plant_r_assignments.m$plant_species,colnames(network))),] # make sure in same order as network
      plant.r <- plant_r_assignments.m$rval # assign r values for simulation
    ## ---- Option 2: Random R-value for all plants ----
    } else {
      plant.r <- runif(1,0,1)
    }
    
    # Run One Simulation: Runs IterativeNodeDeletion() (assumed to simulate species extinctions after plant removals), stores the result
    sim.store[[length(sim.store)+1]] <- IterativeNodeDeletion(network = network, anim.r = anim.r, plant.r = plant.r, removal.taxa = removal.taxa, 
                                                              remove.order = remove.order, silent = TRUE) # run simulation
    # Progress Print Every 1000 Sims
    if(i %% 1000 == 0){
      print(i)
    }
  }
  
  #--------- Calculating species probability of extinction (PE) ---------------------------
  # Extracts response_fates tables from each simulation, which record if a species survived or went extinct in that run.
  raw.vulnerability <- lapply(seq(sim.store), function(x) cbind(sim.store[[x]]$response_fates,x)) # extract the 'response_fates' tables
  raw.vulnerability <- do.call("rbind",raw.vulnerability) # rbind all 'response_fates' tables
  colnames(raw.vulnerability)[4] <- "iteration"
  
  # Aggregates across simulations to compute the probability each species went extinct.
  # fate == 1 means survived → so extinction probability = 1 - (survival rate)
  average.vulnerability <- aggregate(raw.vulnerability$fate, by = list(raw.vulnerability$species), FUN = "sum") # calculate number of times survived vs extinct across all runs
  colnames(average.vulnerability) <- c("species","probability.of.extinction")
  average.vulnerability$probability.of.extinction <- 1 - (average.vulnerability$probability.of.extinction/nsim) # convert to probability of extinction (fates = 1 when extant, and 0 when extinct. Therefore sum of fates divided by nsim gives the proportion of runs it surived. Therefore 1- is the number of times it went extinct)
  
  #--------- Calculating Speed of Extinction (SE) ---------------------------
  # Calculating the average proportion of plant species which had to be removed for a given hummingbird species to go extinct; the speed of extinction (SE)
  # order tells the step (or proportion of plants removed) at which a species went extinct.
  # Normalize by total number of plants.
  raw.order <- raw.vulnerability[complete.cases(raw.vulnerability),] # order excluding NAs
  raw.order$order <- raw.order$order/nplant # convert order to proportion of plant species removed
  # Averages extinction order across simulations.
  average.order <- aggregate(raw.order$order, by = list(raw.order$species), FUN = "mean") # calculate mean proportion of plant species removed
  colnames(average.order) <- c("species","mean.order")
  
  #--------- Output ---------------------------
  return(list(mean_vulnerability = average.vulnerability,
              mean_order = average.order))
}
 # IterNodeDelMultiSim() performs:
# Repeated simulations (n = nsim) of plant species removal from an ecological network.
# In each simulation:
#   Plants are removed one-by-one based on a strategy.
#   Coextinctions of hummingbirds are tracked.
# R-values influence how species respond to removal (resilience).
# Outputs two tables:
#   mean_vulnerability: estimated probability of extinction (PE) for each species.
#   mean_order: estimated speed of extinction (SE) — how soon species go extinct during the plant removal process.