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
    species.goup <- args[3]
} else{
    net.type <- "Year" ## Year or YearSR
    species <- "Plant" ## Plant or Parasite
    species.group <- "Bees" ## Bees or BeesSyrphids
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


load(file=sprintf("../../data/networks/%s_%s_%s.Rdata", net.type,
                  paste(species, collapse=""), species.group
                  ))


