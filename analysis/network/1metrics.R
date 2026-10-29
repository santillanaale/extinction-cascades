rm(list=ls())
setwd("~/")
source("lab_paths.R")
setwd(file.path(local.path, "extinction-cascades/analysis/network"))
source('src/initialize.R')
source('src/vaznull2.R')

# Confirm expanded networks
stopifnot(length(nets) == 61)

q2_metric_names <- c(
  "functional.complementarity.HL",
  "functional.complementarity.LL",
  "FunRedundancy.Pols",
  "FunRedundancy.Plants",
  "mean.number.of.links.HL",
  "mean.number.of.links.LL",
  "number.of.species.HL",
  "number.of.species.LL"
)

q2_metrics <- lapply(nets, function(web) {
  calcMetric(
    web,
    index = c("functional complementarity", "number of species")
  )[q2_metric_names]
})

q2_metrics_df <- data.frame(
  Network = names(q2_metrics),
  do.call(rbind, q2_metrics),
  row.names = NULL,
  check.names = FALSE
)

stopifnot(
  nrow(q2_metrics_df) == 61,
  !anyDuplicated(q2_metrics_df$Network),
  all(is.finite(as.matrix(q2_metrics_df[q2_metric_names])))
)

saveRDS(
  q2_metrics_df,
  "saved/q2_observed_metrics_expanded_61.rds"
)

cat("Saved", nrow(q2_metrics_df),
    "networks with", length(q2_metric_names),
    "observed metrics each.\n")

## ************************************************************
## Calculate observed metrics and null-model z-scores ## beware this takes a while!
## ************************************************************
set.seed(123)
N <- 100

mets <- lapply(nets, calcNetworkMetrics, N=N)

cor.dats <- prepDat(mets,  spec.net, net.type=net.type)

# # Save full null-model outputs
saveRDS(
  mets,
  "saved/q2_null_metrics_expanded_61_N100.rds"
)

cat(
  "Saved", length(mets),
  "networks with", N,
  "null replicates each.\n"
)

# Save complete Q2 metrics for the expanded 61-network dataset
metrics_file <- "saved/corMets_PlantPollinator_YearSR_BeesFliesLepidoptera_N100.Rdata"
save(cor.dats, file = metrics_file)

## ************************************************************
## Prepare a Q2 analysis table with explicit network identifiers
## ************************************************************
load(metrics_file)

cor.dats <- cor.dats |>
  dplyr::rename(Network = Site) |>
  tidyr::separate(
    Network,
    into = c("Site", "Year", "SampleRound"),
    sep = "\\.",
    convert = TRUE,
    remove = FALSE
  )

# Validate the final table
stopifnot(
  nrow(cor.dats) == 61,
  ncol(cor.dats) == 52,
  !anyNA(cor.dats$Site),
  !anyNA(cor.dats$Year),
  !anyNA(cor.dats$SampleRound),
  !anyDuplicated(cor.dats$Network)
)

q2_standardized_names <- c(
  "functional.complementarity.HL",
  "functional.complementarity.LL",
  "FunRedundancy.Pols",
  "FunRedundancy.Plants",
  "mean.number.of.links.HL",
  "mean.number.of.links.LL"
)

q2_z_names <- paste0("z", q2_standardized_names)

stopifnot(
  all(q2_z_names %in% names(cor.dats)),
  all(is.finite(as.matrix(cor.dats[q2_z_names])))
)

# Save the analysis-ready table
saveRDS(
  cor.dats,
  "saved/q2_metrics_table_expanded_61_N100.rds"
)

cat(
  "Saved validated Q2 table:",
  nrow(cor.dats), "rows and",
  ncol(cor.dats), "columns.\n"
)

## ************************************************************
## Historical exploratory models (inactive; retained for reference)
## The original code assumed Site.Year identifiers, not Site.Year.Round.
## Updated climate models belong in 05_climate_models_plots.R.
## ************************************************************
# cor.dats$Year <- factor(cor.dats$Year, levels = c("2012", "2017", "2018", "2021", "2022"))
# cor.dats$Year <- as.numeric(as.character(cor.dats$Year))
# 
# ys <- c("FunRedundancy.Pols",
#         "FunRedundancy.Plants",
#         "functional.complementarity.HL",
#         "functional.complementarity.LL",
#         "mean.number.of.links.HL",
#         "mean.number.of.links.LL")
# 
# xvar <- "Year"
# 
# ## Year
# formulas.year <- lapply(ys, function(x) {
#   as.formula(paste(x, "~ Year + (1|Site)"))
# })
# 
# ## floral richness
# formulas.floral.rich <-lapply(ys, function(x) {
#   as.formula(paste(x, "~",
#                    paste("scale(number.of.species.LL)",
#                          "(1|Site)",
#                          sep="+")))
# })
# ## pollinator richness
# formulas.pol.rich <-lapply(ys, function(x) {
#   as.formula(paste(x, "~",
#                    paste("scale(number.of.species.HL)",
#                          "(1|Site)",
#                          sep="+")))
# })
# 
# ## Models
# mods.year <- lapply(formulas.year, function(fmla) {
#   lmer(fmla, data = cor.dats)
# })
# mods.floral.rich <- lapply(formulas.floral.rich, function(x){
#   lmer(x, data=cor.dats)
# })
# mods.pol.rich <- lapply(formulas.pol.rich, function(x){
#   lmer(x, data=cor.dats)
# })
# 
# 
# names(mods.year) <- names(mods.floral.rich) <-
#   names(mods.pol.rich) <- ys
# 
# 
# ## pollinator redundancy/complementarity/generalization should be regressed against
# ## the floral community, similarly plant redundancy/complementarity/generalization
# ## should be regressed against the pollinator
# 
# mods.rich <- c(mods.pol.rich["FunRedundancy.Pols"],
#                mods.floral.rich["FunRedundancy.Plants"],
#                mods.pol.rich["functional.complementarity.HL"],
#                mods.floral.rich["functional.complementarity.LL"],
#                mods.pol.rich["mean.number.of.links.LL"],
#                mods.floral.rich["mean.number.of.links.HL"])
# 
# ## results
# lapply(mods.year, summary)
# lapply(mods.rich, summary)
# 
# ## check sig levels with method other than wald CI
# lapply(mods.year, anova)
# lapply(mods.rich, anova)
# 
# aics <- cbind(sapply(mods.year, AIC), sapply(mods.rich, AIC), NA)
# aics[,3]  <- aics[,2]- aics[,1]
# colnames(aics) <- c("year", "richness", "deltaAIC")
# aics
# 
# save(mods.year, mods.rich, cor.dats,
#      file=file.path(save.path, 'mods/metrics.Rdata'))