# ================================================================
# 02_community_composition.R
# Q1: Community composition and diversity
#
# Unit: Site x Year x SampleRound
# Source: Prepared specimen records
#
# Outputs:
#   Q1_community_composition.csv
#   Q1_network_count_audit.csv
#
# Climate joins and statistical models belong in
# 05_climate_models_plots.R.
# ================================================================

rm(list = ls())

library(tidyverse)
library(vegan)

setwd("~/")
source("lab_paths.R")
setwd(file.path(local.path, "skyIslands"))

# ================================================================
# 1. Load prepared specimen records
# ================================================================

spec <- read_csv(
  "data/spec_all_methods.csv",
  col_types = cols(
    WindEnd = col_character(),
    .default = col_guess()
  ),
  show_col_types = FALSE
)

# Confirm no remaining parsing problems
stopifnot(nrow(problems(spec)) == 0)

stopifnot(
  all(c(
    "Site", "Year", "SampleRound", "Method",
    "Order", "Family", "GenusSpecies",
    "PlantGenusSpecies"
  ) %in% names(spec))
)

# ================================================================
# 2. Define canonical taxonomic universe
# ================================================================

bee.families <- c(
  "Andrenidae", "Apidae", "Colletidae",
  "Halictidae", "Megachilidae"
)

pollinator.families <- c(
  bee.families,
  "Syrphidae",
  "Bombyliidae"
)

# Select Net-collected records belonging to the network taxa.
# Require both pollinator and plant identifications.
#
# This selection was reconciled with the audited inventory:
# 14,308 records across 61 networks (November 2026).

spec.net.selected <- spec %>%
  filter(
    Method == "Net",
    !is.na(SampleRound),
    SampleRound != 0,
    Family %in% pollinator.families |
      Order == "Lepidoptera",
    !is.na(GenusSpecies),
    str_trim(GenusSpecies) != "",
    !is.na(PlantGenusSpecies),
    str_trim(PlantGenusSpecies) != ""
  ) %>%
  mutate(
    Network = paste(
      Site, Year, SampleRound, sep = "."
    )
  )

cat("\nSelected records:", nrow(spec.net.selected))
cat("\nNetworks:", n_distinct(spec.net.selected$Network))
cat("\nSites:", n_distinct(spec.net.selected$Site), "\n")

stopifnot(
  n_distinct(spec.net.selected$Network) == 61
)

# ================================================================
# 3. Audit against existing network inventory
# ================================================================
setwd(file.path(local.path, "extinction-cascades"))
inventory <- read_csv(
  "data/audit/master_sampling_inventory.csv",
  show_col_types = FALSE
)

q1_counts <- spec.net.selected %>%
  count(Network, name = "n_specimens")

count_audit <- inventory %>%
  select(Network, Interactions, N_Selected) %>%
  full_join(q1_counts, by = "Network") %>%
  mutate(
    difference = n_specimens - Interactions
  )

cat("\nNetwork-count audit:\n")

print(
  count_audit %>%
    summarise(
      n_networks = n(),
      n_missing = sum(is.na(n_specimens)),
      n_mismatches = sum(
        difference != 0, na.rm = TRUE
      ),
      total_difference = sum(
        difference, na.rm = TRUE
      )
    )
)

cat("\nNetworks with discrepancies:\n")

print(
  count_audit %>%
    filter(is.na(difference) | difference != 0)
)

# Preserve audit results for further reconciliation.
write_csv(
  count_audit,
  "data/audit/Q1_network_count_audit.csv"
)

# Fail if records or network coverage no longer match the inventory.
stopifnot(
  nrow(count_audit) == 61,
  !anyNA(count_audit$difference),
  all(count_audit$difference == 0),
  nrow(spec.net.selected) == 14308
)

# ================================================================
# 4. Assign mutually exclusive pollinator guilds
# ================================================================

spec.q1 <- spec.net.selected %>%
  mutate(
    Guild = case_when(
      GenusSpecies == "Apis mellifera" ~ "Honeybee",
      Family %in% bee.families ~ "OtherBees",
      Family == "Syrphidae" ~ "Syrphidae",
      Family == "Bombyliidae" ~ "Bombyliidae",
      Order == "Lepidoptera" ~ "Lepidoptera",
      TRUE ~ NA_character_
    )
  )

stopifnot(!anyNA(spec.q1$Guild))

cat("\nGuild record counts:\n")

print(
  spec.q1 %>%
    count(Guild, sort = TRUE)
)

# ================================================================
# 5. Guild-specific abundance, richness and Shannon diversity
# ================================================================

# Shannon diversity calculated from specimen/interaction
# record frequencies per pollinator taxon.

shannon_diversity <- function(x) {
  counts <- table(x)
  if (length(counts) == 0) return(NA_real_)
  vegan::diversity(counts, index = "shannon")
}

guild_metrics <- spec.q1 %>%
  group_by(
    Network, Site, Year, SampleRound, Guild
  ) %>%
  summarise(
    Abundance = n(),
    Richness = n_distinct(GenusSpecies),
    Shannon = shannon_diversity(GenusSpecies),
    .groups = "drop"
  )

# Create all five guild combinations for every network.
# Absent guilds have zero abundance/richness.
# Shannon is NA for absent guilds (undefined).

network_key <- spec.q1 %>%
  distinct(Network, Site, Year, SampleRound)

guild_levels <- c(
  "OtherBees", "Honeybee", "Syrphidae",
  "Bombyliidae", "Lepidoptera"
)

guild_metrics <- network_key %>%
  tidyr::crossing(Guild = guild_levels) %>%
  left_join(
    guild_metrics,
    by = c(
      "Network", "Site", "Year",
      "SampleRound", "Guild"
    )
  ) %>%
  mutate(
    Abundance = replace_na(Abundance, 0L),
    Richness = replace_na(Richness, 0L)
  )

guild_wide <- guild_metrics %>%
  pivot_wider(
    id_cols = c(Network, Site, Year, SampleRound),
    names_from = Guild,
    values_from = c(Abundance, Richness, Shannon),
    names_glue = "{.value}_{Guild}"
  )

# ================================================================
# 6. Visited-plant richness and diversity
# ================================================================

plant_metrics <- spec.q1 %>%
  group_by(Network, Site, Year, SampleRound) %>%
  summarise(
    PlantRichness = n_distinct(PlantGenusSpecies),
    PlantShannon = shannon_diversity(
      PlantGenusSpecies
    ),
    TotalRecords = n(),
    PollinatorRichness = n_distinct(
      GenusSpecies
    ),
    PollinatorShannon = shannon_diversity(
      GenusSpecies
    ),
    .groups = "drop"
  )

# ================================================================
# 7. Create Q1 response table
# ================================================================

q1_community <- guild_wide %>%
  left_join(
    plant_metrics,
    by = c(
      "Network", "Site", "Year", "SampleRound"
    )
  )

# ================================================================
# 8. Validate Q1 metrics
# ================================================================

stopifnot(
  nrow(q1_community) == 61,
  !anyDuplicated(q1_community$Network)
)

abundance_columns <- paste0(
  "Abundance_", guild_levels
)

q1_community <- q1_community %>%
  mutate(
    GuildAbundanceTotal = rowSums(
      across(all_of(abundance_columns))
    ),
    GuildCountDifference =
      GuildAbundanceTotal - TotalRecords
  )

stopifnot(
  all(q1_community$GuildCountDifference == 0)
)

cat("\nQ1 table dimensions:\n")
print(dim(q1_community))

cat("\nGuild presence across networks:\n")

print(
  guild_metrics %>%
    group_by(Guild) %>%
    summarise(
      networks_present = sum(Abundance > 0),
      networks_absent = sum(Abundance == 0),
      .groups = "drop"
    )
)

cat("\nMissing values in Q1 table:\n")

print(
  colSums(is.na(q1_community))
)

# ================================================================
# 9. Export Q1 response table
# ================================================================

write_csv(
  q1_community,
  "data/Q1_community_composition.csv"
)

cat("\nQ1 community composition table saved.\n")
cat("Network-count reconciliation: PASSED.\n")
cat("Guild and plant metrics: review remaining audit checks before final modeling.\n")