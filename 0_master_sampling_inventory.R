rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
set.seed(123)
setwd("~/")
source("lab_paths.R")

# ---- Load Networks ----
setwd(file.path(local.path, "skyIslands"))

# 1. Load networks into a separate environment
new.nets <- new.env()

load("data/networks/YearSR_PlantPollinator_BeesFliesLepidoptera.Rdata",
  envir = new.nets
)

nets <- new.nets$nets

# 2. Construct master network inventory
master_sampling_inventory <- tibble(
  Network = names(nets),
  Plants = vapply(nets, nrow, integer(1)),
  Pollinators = vapply(nets, ncol, integer(1)),
  Interactions = vapply(
    nets, function(x) sum(x), numeric(1)
  ),
  UniqueLinks = vapply(
    nets, function(x) sum(x > 0), integer(1)
  )
) %>%
  separate(
    Network,
    into = c("Site", "Year", "SampleRound"),
    sep = "\\.",
    remove = FALSE
  ) %>%
  mutate(
    Year = as.integer(Year),
    SampleRound = as.integer(SampleRound),
    YearSR = Network,
    
    # Previously excluded network
    Excluded = Network == "PL.2017.1",
    ExclusionReason = if_else(
      Excluded,
      "Previously excluded; reason to verify",
      NA_character_
    ),
    
    # Preliminary network eligibility
    NetworkAvailable = TRUE,
    NetworkIncluded = !Excluded,
    
    # Additional information to audit later
    StartDate = as.Date(NA),
    EndDate = as.Date(NA),
    ClimateAvailable = NA,
    PC1Available = NA
  ) %>%
  arrange(Site, Year, SampleRound)

# 3. Verify unique sampling units
stopifnot(
  !anyDuplicated(master_sampling_inventory$YearSR)
)

# 4. Summary by site and year
sampling_summary <- master_sampling_inventory %>%
  group_by(Site, Year) %>%
  summarise(
    N_Rounds = n(),
    N_Included = sum(NetworkIncluded),
    TotalInteractions = sum(Interactions),
    .groups = "drop"
  ) %>%
  arrange(Site, Year)

# 5. Overall audit summary
audit_summary <- master_sampling_inventory %>%
  summarise(
    N_Sites = n_distinct(Site),
    N_Years = n_distinct(Year),
    N_SiteYears = n_distinct(
      paste(Site, Year)
    ),
    N_Networks = n(),
    N_Included = sum(NetworkIncluded),
    N_Excluded = sum(Excluded),
    TotalInteractions = sum(Interactions)
  )

print(audit_summary)
print(sampling_summary, n = Inf)

# 6. Save reproducible outputs
dir.create("../extinction-cascades/data/audit", showWarnings = FALSE,
           recursive = TRUE)

saveRDS(
  master_sampling_inventory,
  "../extinction-cascades/data/audit/master_sampling_inventory.rds"
)

write_csv(
  master_sampling_inventory,
  "../extinction-cascades/data/audit/master_sampling_inventory.csv"
)

write_csv(
  sampling_summary,
  "../extinction-cascades/data/audit/sampling_summary.csv"
)

write_csv(
  audit_summary,
  "../extinction-cascades/data/audit/audit_summary.csv"
)

#################################################
# Add sampling dates to the master inventory
master_sampling_inventory <- master_sampling_inventory %>%
  select(-StartDate, -EndDate) %>%
  left_join(
    specimen_inventory %>%
      select(Site, Year, SampleRound, StartDate, EndDate),
    by = c("Site", "Year", "SampleRound")
  )

# Check for missing dates
master_sampling_inventory %>%
  filter(is.na(StartDate) | is.na(EndDate))

# Save updated inventory
saveRDS(
  master_sampling_inventory,
  "../data/audit/master_sampling_inventory.rds"
)

readr::write_csv(
  master_sampling_inventory,
  "../data/audit/master_sampling_inventory.csv"
)

#################################################
# Add the reconciliation result to master inventory
master_sampling_inventory <- master_sampling_inventory %>%
  left_join(
    selected_counts,
    by = c("Site", "Year", "SampleRound")
  ) %>%
  mutate(
    InteractionCountVerified =
      Interactions == N_Selected
  )

stopifnot(
  nrow(master_sampling_inventory) == 61,
  all(master_sampling_inventory$InteractionCountVerified)
)

saveRDS(
  master_sampling_inventory,
  "../data/audit/master_sampling_inventory.rds"
)

readr::write_csv(
  master_sampling_inventory,
  "../data/audit/master_sampling_inventory.csv"
)


#################################################
library(dplyr)
library(tidyr)

guild_inventory <- spec.net.selected %>%
  mutate(
    Guild = case_when(
      GenusSpecies == "Apis mellifera" ~ "Honeybee",
      Family %in% bee.families ~ "OtherBees",
      Family == "Syrphidae" ~ "Syrphids",
      Family == "Bombyliidae" ~ "BeeFlies",
      Order == "Lepidoptera" ~ "Lepidoptera",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(Site, Year, SampleRound, Guild) %>%
  summarise(
    Abundance = n(),
    Richness = n_distinct(GenusSpecies),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Guild,
    values_from = c(Abundance, Richness),
    values_fill = 0
  )

master_sampling_inventory <- master_sampling_inventory %>%
  left_join(
    guild_inventory,
    by = c("Site", "Year", "SampleRound")
  )

stopifnot(nrow(master_sampling_inventory) == 61)

saveRDS(
  master_sampling_inventory,
  "../../extinction-cascades/data/audit/master_sampling_inventory.rds"
)

readr::write_csv(
  master_sampling_inventory,
  "../../extinction-cascades/data/audit/master_sampling_inventory.csv"
)
