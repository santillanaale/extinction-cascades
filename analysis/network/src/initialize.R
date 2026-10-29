library(lme4)
library(lmerTest)
library(RColorBrewer)

source('src/CalcMetrics.R')
source('src/misc.R')

save.path <- 'saved'

load('../../data/spec_net.Rdata')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 0){
    net.type <- args[1]
    species <- args[2]
    species.group <- args[3]
} else{
  net.type <- "YearSR"
  species <- "Plant"
  species.group <- "BeesFliesLepidoptera"
}


plants <- unique(spec.net$PlantGenusSpecies)
pols <- unique(spec.net$GenusSpecies)
parasites <- c("AscosphaeraSpp", "ApicystisSpp",
               "CrithidiaExpoeki", "CrithidiaBombi", "NosemaBombi",
               "NosemaCeranae")

if(species == "Plant"){
    species <- c("Plant", "Pollinator")
    lower.level  <- plants
    higher.level <- pols
} else if(species == "Parasite"){
    species <- c("Pollinator", "Parasite")
    lower.level  <- pols
    higher.level <- parasites
}
if(net.type == "YearSR"){
    nets.by.SR  <- TRUE
} else {
    nets.by.SR  <- FALSE
}

# Load expanded networks directly from skyIslands
network_file <- file.path(
  local.path,
  "skyIslands",
  "data",
  "networks",
  sprintf(
    "%s_%s_%s.Rdata",
    net.type,
    paste(species, collapse = ""),
    species.group
  )
)

stopifnot(file.exists(network_file))

load(network_file)

# Confirm expected networks
stopifnot(length(nets) == 61)

cat("Loaded", length(nets), "networks from:\n", network_file, "\n")


