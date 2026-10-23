run_scm_robustness <- function(
    web,
    R_val_rows,
    R_val_cols,
    TargetGuild,
    TargetSpecies,
    nrep = 100
) {
  
  nrows <- nrow(web)
  ncols <- ncol(web)
  total_species <- nrows + ncols
  
  # Preallocate a results matrix: rows = nrep, cols = 2 (prop_cascade, robustness)
  results <- matrix(0, nrow = nrep, ncol = 2)
  colnames(results) <- c("prop_cascade", "robustness")
  
  for (rep_idx in seq_len(nrep)) {
    sim <- netcascade_mod(
      imatrix = web,
      R_val_rows = R_val_rows,
      R_val_cols = R_val_cols,
      TargetGuild = TargetGuild,
      TargetSpecies = TargetSpecies
    )
    
    n_extinct <- length(unique(sim$rows_species_data$lost_rows)) +
      length(unique(sim$cols_species_data$lost_cols))
    
    prop_cascade <- n_extinct / total_species
    robustness <- 1 - prop_cascade
    
    results[rep_idx, "prop_cascade"] <- prop_cascade
    results[rep_idx, "robustness"] <- robustness
  }
  
  # Compute mean values
  tibble(
    prop_cascade = mean(results[, "prop_cascade"]),
    SCM_robustness = mean(results[, "robustness"])
  )
}


# run_scm_robustness <- function(
#     web,
#     R_val_rows,
#     R_val_cols,
#     TargetGuild,
#     TargetSpecies,
#     nrep = 100
# ) {
#   
#   nrows <- nrow(web)
#   ncols <- ncol(web)
#   total_species <- nrows + ncols
#   
#   results <- replicate(nrep, {
#     
#     sim <- netcascade_mod(
#       imatrix = web,
#       R_val_rows = R_val_rows,
#       R_val_cols = R_val_cols,
#       TargetGuild = TargetGuild,
#       TargetSpecies = TargetSpecies
#     )
#     
#     n_extinct <- length(unique(sim$rows_species_data$lost_rows)) +
#       length(unique(sim$cols_species_data$lost_cols))
#     
#     prop_cascade <- n_extinct / total_species
#     robustness <- 1 - prop_cascade
#     
#     c(
#       prop_cascade = prop_cascade,
#       robustness = robustness
#     )
#   })
#   
#   tibble(
#     prop_cascade = mean(results["prop_cascade", ]),
#     SCM_robustness = mean(results["robustness", ])
#   )
# }
