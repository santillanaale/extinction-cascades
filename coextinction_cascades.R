rm(list=ls())
library(bipartite)
library(fossil)
library(tidyverse)
library(ggplot2)
library(dplyr)


#skyIslands/data/networks/Yr_PlantPollinator_Bees.RData #yearly networks for bees
load("../skyIslands/data/networks/Yr_PlantPollinator_Bees.RData") #sample round networks for bees
#Lili knows that sample rounds are not equal across all years
#Years with less rounds would look a lot less diverse
view(nets$SM.2018.1)

networklevel(nets[["CH.2018.1"]])
plotweb(nets[["CH.2018.1", method="normal", text.rot=90, labsize]])

#bipartite has a robustness index


load("../../../ca_survey_saved/data_original/networks/nets.Rdata")