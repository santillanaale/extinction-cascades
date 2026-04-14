rm(list=ls())
library(bipartite)
library(tidyverse)
set.seed(123)
setwd("~/")
source("lab_paths.R")

# ---- Load Networks ----
setwd(file.path(local.path, "skyIslands"))
load("data/networks/YearSR_PlantPollinator_Bees.Rdata") #Site x Year x Survey Round

## ---- Make net key ----
net_meta <- tibble(net_id = names(nets)) %>%
  separate(net_id, into = c("Site", "Year", "SampleRound"), sep = "\\.") %>%
  mutate(
    Year = as.numeric(Year),
    SampleRound = as.numeric(SampleRound)
  )

## ---- Remove bad network ----
bad_net <- "PL.2017.1"
keep    <- names(nets) != bad_net
nets     <- nets[keep]
net_meta <- net_meta[keep, ]

# ---- Load extinction functions ----
setwd("~/")
setwd(file.path(local.path, "extinction-cascades"))

source("modified_extinction.R")
source("modified_second_extinct.R")

## *******************************************************************
# ---- TCM: Topological Cascade Model ----
## *******************************************************************

## ---- Helper function ----
run_robustness <- function(web, scenario, nrep = 100) {
  
  stopifnot(scenario$method %in% c("abundance", "degree"))
  
  ext <- modified.second.extinct(
    web         = web,
    participant = "lower",      # plants
    method      = scenario$method,
    reverse     = scenario$reverse,
    nrep        = nrep,
    details     = FALSE
  )
  
  robustness(ext)
}

## ---- Define scenarios ----
scenarios <- tibble(
  scenario = c(
    "abundance_low",   # remove rare first
    "abundance_high",  # remove dominant first
    "degree_low",      # remove specialists first
    "degree_high"      # remove generalists first
  ),
  method  = c("abundance", "abundance", "degree", "degree"),
  reverse = c(FALSE, TRUE, FALSE, TRUE)
)

## ---- Robustness loop (Site x Year x SampleRound x scenario) ----
robustness_results <- map_dfr(seq_len(nrow(scenarios)), function(s) {
  
  scenario <- scenarios[s, ]
  message("Running TCM scenario: ", scenario$scenario)
  
  map_dfr(seq_along(nets), function(i) {
    
    rob <- run_robustness(
      web      = nets[[i]],
      scenario = scenario,
      nrep     = 100
    )
    
    tibble(
      Site        = net_meta$Site[i],
      Year        = net_meta$Year[i],
      SampleRound = net_meta$SampleRound[i],
      scenario    = scenario$scenario,
      Robustness  = rob
    )
  })
})

## *******************************************************************
# ---- SCM: Stochastic Cascade Model ----
## *******************************************************************

source("SCM/netcascade_mod.R")
source("SCM/scm_robustness_helpers.R")

## ---- Define SCM scenarios ----
scm_scenarios <- tibble(
  scenario      = c("scm_random_plant", "scm_dominant_plant"),
  TargetGuild   = c("rows", "rows"),
  TargetSpecies = c("random_binary", "random_abundance")
)

## ---- SCM robustness loop (Site x Year x SampleRound x scenario) ----
scm_results <- map_dfr(seq_len(nrow(scm_scenarios)), function(s) {
  
  scenario <- scm_scenarios[s, ]
  message("Running SCM scenario: ", scenario$scenario)
  
  map_dfr(seq_along(nets), function(i) {
    
    web    <- as.matrix(nets[[i]])
    R_rows <- rep(0.5, nrow(web))
    R_cols <- rep(0.5, ncol(web))
    
    out <- run_scm_robustness(
      web           = web,
      R_val_rows    = R_rows,
      R_val_cols    = R_cols,
      TargetGuild   = scenario$TargetGuild,
      TargetSpecies = scenario$TargetSpecies,
      nrep          = 100
    )
    
    tibble(
      Site          = net_meta$Site[i],
      Year          = net_meta$Year[i],
      SampleRound   = net_meta$SampleRound[i],
      scenario      = scenario$scenario,
      prop_cascade  = out$prop_cascade,
      SCM_Robustness = out$SCM_robustness
    )
  })
})

## *******************************************************************
# ---- Merge climate and save ----
## *******************************************************************

setwd("~/")
setwd(file.path(local.path, "skyIslands_saved"))
climate <- read.csv("data/relational/original/climate.csv")

setwd("~/")
setwd(file.path(local.path, "extinction-cascades"))

tcm_climate <- robustness_results %>%
  left_join(climate, by = c("Site", "Year", "SampleRound"))

scm_climate <- scm_results %>%
  left_join(climate, by = c("Site", "Year", "SampleRound"))

save(
  robustness_results,
  scm_results,
  tcm_climate,
  scm_climate,
  file = "analysis/network/robustness_metrics/YearSR_robustness_climate.Rdata"
)