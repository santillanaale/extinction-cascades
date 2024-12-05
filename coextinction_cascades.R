rm(list=ls())
library(bipartite)
library(fossil)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(networkD3)


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













