setwd("C:/Users/ale_s/University of Oregon Dropbox/Alejandro Santillana Fernandez/")
rm(list=ls())
setwd('extinction-cascades/analysis/network')
source('src/initialize.R')
source('src/vaznull2.R')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 0){
    N <- as.numeric(args[4])
} else{
    N <- 2
}

## ************************************************************
## calculate metrics and zscores ## beware this takes a while!
## ************************************************************

nets <- nets[apply(sapply(nets, dim) > 2, 2, all)]

mets <- lapply(nets, calcNetworkMetrics,
               N=N)

cor.dats <- prepDat(mets,  spec.net, net.type=net.type)

save(cor.dats, file=sprintf('saved/corMets_%s_%s.Rdata',
                             paste(species, collapse=""),
                             net.type))

## ************************************************************
load(file=sprintf('saved/corMets_%s_%s.Rdata',
                            paste(species, collapse=""), net.type))

cor.dats <- separate(cor.dats, col = Site, into = c("Site", "Year"), sep = "\\.")
cor.dats$Year <- factor(cor.dats$Year, levels = c("2012", "2017", "2018", "2021", "2022"))
cor.dats$Year <- as.numeric(as.character(cor.dats$Year))

ys <- c("FunRedundancy.Pols",
        "FunRedundancy.Plants",
        "functional.complementarity.HL",
        "functional.complementarity.LL",
        "mean.number.of.links.HL",
        "mean.number.of.links.LL")

xvar <- "Year"

## Year
formulas.year <- lapply(ys, function(x) {
  as.formula(paste(x, "~ Year + (1|Site)"))
})

## floral richness
formulas.floral.rich <-lapply(ys, function(x) {
  as.formula(paste(x, "~",
                   paste("scale(number.of.species.LL)",
                         "(1|Site)",
                         sep="+")))
})
## pollinator richness
formulas.pol.rich <-lapply(ys, function(x) {
  as.formula(paste(x, "~",
                   paste("scale(number.of.species.HL)",
                         "(1|Site)",
                         sep="+")))
})

## Models
mods.year <- lapply(formulas.year, function(fmla) {
  lmer(fmla, data = cor.dats)
})
mods.floral.rich <- lapply(formulas.floral.rich, function(x){
  lmer(x, data=cor.dats)
})
mods.pol.rich <- lapply(formulas.pol.rich, function(x){
  lmer(x, data=cor.dats)
})


names(mods.year) <- names(mods.floral.rich) <-
  names(mods.pol.rich) <- ys


## pollinator redundancy/complementarity/generalization should be regressed against
## the floral community, similarly plant redundancy/complementarity/generalization
## should be regressed against the pollinator

mods.rich <- c(mods.pol.rich["FunRedundancy.Pols"],
               mods.floral.rich["FunRedundancy.Plants"],
               mods.pol.rich["functional.complementarity.HL"],
               mods.floral.rich["functional.complementarity.LL"],
               mods.pol.rich["mean.number.of.links.LL"],
               mods.floral.rich["mean.number.of.links.HL"])

## results
lapply(mods.year, summary)
lapply(mods.rich, summary)

## check sig levels with method other than wald CI
lapply(mods.year, anova)
lapply(mods.rich, anova)

aics <- cbind(sapply(mods.year, AIC), sapply(mods.rich, AIC), NA)
aics[,3]  <- aics[,2]- aics[,1]
colnames(aics) <- c("year", "richness", "deltaAIC")
aics

save(mods.year, mods.rich, cor.dats,
     file=file.path(save.path, 'mods/metrics.Rdata'))