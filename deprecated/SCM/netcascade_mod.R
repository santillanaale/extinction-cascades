#Pires et al 2020
#adapted from Vieira and Almeida-Neto 2015; doi: 10.1111/ele.12394
#Simulates extinction cascades by removing a target species
#Received from Susan Waters 1/20/2026

netcascade_mod <- function(
    imatrix,
    R_val_rows,
    R_val_cols,
    TargetGuild,
    TargetSpecies,
    extinct_cols = NULL,
    extinct_rows = NULL,
    return.matrix = FALSE
) {
  
  # --------- ARGUMENT CHECKS -----------------------------------
  if (!is.matrix(imatrix) || (nrow(imatrix) + ncol(imatrix)) < 3) {
    stop("'imatrix' must be a matrix with at least three species")
  }
  
  if (length(R_val_rows) != nrow(imatrix)) {
    stop("Length of R_val_rows must equal number of rows in imatrix")
  }
  
  if (length(R_val_cols) != ncol(imatrix)) {
    stop("Length of R_val_cols must equal number of columns in imatrix")
  }
  
  if (!(TargetGuild %in% c("rows", "cols", "random_richness", "random_binary"))) {
    stop("Invalid TargetGuild")
  }
  
  if (!is.numeric(TargetSpecies) && !(TargetSpecies %in% c("random_abundance", "random_binary"))) {
    stop("Invalid TargetSpecies")
  }
  
  if (!is.null(extinct_cols) && !is.integer(extinct_cols)) {
    stop("extinct_cols must be NULL or integer vector")
  }
  
  if (!is.null(extinct_rows) && !is.integer(extinct_rows)) {
    stop("extinct_rows must be NULL or integer vector")
  }
  
  # --------- DEFINING OBJECTS ----------------------------------
  nrows <- nrow(imatrix)
  ncols <- ncol(imatrix)
  rows <- 1:nrows
  cols <- 1:ncols
  rowsNA <- rows
  colsNA <- cols
  rowsNA[extinct_rows] <- NA
  colsNA[extinct_cols] <- NA
  degree_when_lost_rows <- c()
  degree_when_lost_cols <- c()
  
  # --------- DEPENDENCE MATRICES -------------------------------
  interaction_strength <- array(0, dim = c(nrows, ncols, 2))
  # columns’ dependence on rows
  interaction_strength[, , 1] <- t(t(imatrix) / rowSums(t(imatrix)))
  # rows’ dependence on columns
  interaction_strength[, , 2] <- imatrix / rowSums(imatrix)
  
  # --------- CHOOSING PRIMARY EXTINCTION -----------------------
  coext_rows <- c()
  coext_cols <- c()
  
  if (is.numeric(TargetSpecies)) {
    if (length(TargetSpecies) == 1) {
      if (TargetGuild == "rows") {
        if (TargetSpecies %in% extinct_rows) stop("Target species already extinct")
        coext_rows <- TargetSpecies
        degree_when_lost_rows <- 1
      }
      if (TargetGuild == "cols") {
        if (TargetSpecies %in% extinct_cols) stop("Target species already extinct")
        coext_cols <- TargetSpecies
        degree_when_lost_cols <- 1
      }
    } else {
      nspecies <- switch(TargetGuild, rows = nrows, cols = ncols)
      if (length(TargetSpecies) != nspecies) stop("TargetSpecies length mismatch")
      if (TargetGuild == "rows") {
        alive_rows <- rows[!is.na(rowsNA)]
        coext_rows <- sample(c(alive_rows, 0), 1, prob = c(TargetSpecies[alive_rows], 0))
        degree_when_lost_rows <- 1
      }
      if (TargetGuild == "cols") {
        alive_cols <- cols[!is.na(colsNA)]
        coext_cols <- sample(c(alive_cols, 0), 1, prob = c(TargetSpecies[alive_cols], 0))
        degree_when_lost_cols <- 1
      }
    }
  } else if (TargetSpecies == "random_abundance") {
    if (TargetGuild == "rows") {
      alive_rows <- rows[!is.na(rowsNA)]
      coext_rows <- sample(alive_rows, 1, prob = rowSums(imatrix[alive_rows, ]))
      degree_when_lost_rows <- 1
    }
    if (TargetGuild == "cols") {
      alive_cols <- cols[!is.na(colsNA)]
      coext_cols <- sample(alive_cols, 1, prob = colSums(imatrix[, alive_cols]))
      degree_when_lost_cols <- 1
    }
  } else if (TargetSpecies == "random_binary") {
    if (TargetGuild == "rows") {
      alive_rows <- rows[!is.na(rowsNA)]
      coext_rows <- sample(alive_rows, 1)
      degree_when_lost_rows <- 1
    }
    if (TargetGuild == "cols") {
      alive_cols <- cols[!is.na(colsNA)]
      coext_cols <- sample(alive_cols, 1)
      degree_when_lost_cols <- 1
    }
  }
  
  imatrix[coext_rows, ] <- 0
  imatrix[, coext_cols] <- 0
  lost_rows <- coext_rows
  lost_cols <- coext_cols
  
  # --------- CASCADE LOOP -------------------------------------
  equilibrium <- FALSE
  degree <- 1
  degree_table <- data.frame(
    degree = degree,
    guild = factor(TargetGuild, levels = c("rows", "cols")),
    n_extinctions = 1
  )
  
  while (!equilibrium) {
    extinct_rows <- coext_rows
    extinct_cols <- coext_cols
    rowsNA[extinct_rows] <- NA
    colsNA[extinct_cols] <- NA
    remaining_rows <- rows[!is.na(rowsNA)]
    remaining_cols <- cols[!is.na(colsNA)]
    
    if (length(extinct_rows) > 0) {
      for (i in extinct_rows) {
        Target <- R_val_cols[remaining_cols] * interaction_strength[i, remaining_cols, 1] > runif(length(remaining_cols))
        coext_cols <- c(coext_cols, remaining_cols[Target])
      }
      coext_cols <- unique(coext_cols)
      imatrix[, coext_cols] <- 0
      lost_cols <- c(lost_cols, coext_cols)
      if (length(coext_cols) > 0) {
        degree <- degree + 1
        degree_when_lost_cols <- c(degree_when_lost_cols, rep(degree, length(coext_cols)))
        degree_table[degree, ] <- data.frame(degree, "cols", length(coext_cols))
      }
      coext_rows <- c()
    } else {
      for (i in extinct_cols) {
        Target <- R_val_rows[remaining_rows] * interaction_strength[remaining_rows, i, 2] > runif(length(remaining_rows))
        coext_rows <- c(coext_rows, remaining_rows[Target])
      }
      coext_rows <- unique(coext_rows)
      imatrix[coext_rows, ] <- 0
      lost_rows <- c(lost_rows, coext_rows)
      if (length(coext_rows) > 0) {
        degree <- degree + 1
        degree_when_lost_rows <- c(degree_when_lost_rows, rep(degree, length(coext_rows)))
        degree_table[degree, ] <- data.frame(degree, "rows", length(coext_rows))
      }
      coext_cols <- c()
    }
    
    equilibrium <- (length(coext_rows) + length(coext_cols)) == 0
  }
  
  # --------- OUTPUT -------------------------------------------
  out <- list(
    cascade_data = degree_table,               # extinctions by degree
    rows_species_data = data.frame(            # which row species went extinct
      lost_rows = lost_rows,
      degree_of_extinction = degree_when_lost_rows
    ),
    cols_species_data = data.frame(            # which col species went extinct
      lost_cols = lost_cols,
      degree_of_extinction = degree_when_lost_cols
    )
  )
  
  if (return.matrix) {
    out$interaction_matrix <- imatrix
  }
  
  return(out)
}



# netcascade_mod <- function(
#     imatrix,
#     R_val_rows,
#     R_val_cols,
#     TargetGuild,
#     TargetSpecies,
#     extinct_cols = NULL,
#     extinct_rows = NULL,
#     return.matrix = FALSE
# ) {
#   
#   # --------- ARGUMENT CHECKS -----------------------------------
#   if (!is.matrix(imatrix) || (nrow(imatrix) + ncol(imatrix)) < 3) {
#     stop("'imatrix' must be a matrix with at least three species")
#   }
#   
#   if (length(R_val_rows) != nrow(imatrix)) stop("Length of R_val_rows must equal number of rows in imatrix")
#   if (length(R_val_cols) != ncol(imatrix)) stop("Length of R_val_cols must equal number of columns in imatrix")
#   if (!(TargetGuild %in% c("rows", "cols", "random_richness", "random_binary"))) stop("Invalid TargetGuild")
#   if (!is.numeric(TargetSpecies) && !(TargetSpecies %in% c("random_abundance", "random_binary"))) stop("Invalid TargetSpecies")
#   
#   if (!is.null(extinct_cols) && !is.integer(extinct_cols)) stop("extinct_cols must be NULL or integer vector")
#   if (!is.null(extinct_rows) && !is.integer(extinct_rows)) stop("extinct_rows must be NULL or integer vector")
#   
#   # --------- INITIALIZATION ------------------------------------
#   nrows <- nrow(imatrix)
#   ncols <- ncol(imatrix)
#   
#   rows <- 1:nrows
#   cols <- 1:ncols
#   
#   rowsNA <- rows
#   colsNA <- cols
#   
#   if (!is.null(extinct_rows)) rowsNA[extinct_rows] <- NA
#   if (!is.null(extinct_cols)) colsNA[extinct_cols] <- NA
#   
#   degree_when_lost_rows <- c()
#   degree_when_lost_cols <- c()
#   
#   # --------- DEPENDENCE MATRICES -------------------------------
#   interaction_strength <- array(0, dim = c(nrows, ncols, 2))
#   if (nrows > 0 && ncols > 0) {
#     interaction_strength[, , 1] <- t(t(imatrix) / rowSums(t(imatrix)))  # columns depend on rows
#     interaction_strength[, , 2] <- imatrix / rowSums(imatrix)           # rows depend on columns
#   }
#   
#   # --------- PRIMARY EXTINCTION -------------------------------
#   coext_rows <- integer(0)
#   coext_cols <- integer(0)
#   
#   alive_rows <- rows[!is.na(rowsNA)]
#   alive_cols <- cols[!is.na(colsNA)]
#   
#   if (is.numeric(TargetSpecies)) {
#     if (TargetGuild == "rows") coext_rows <- TargetSpecies
#     if (TargetGuild == "cols") coext_cols <- TargetSpecies
#   } else if (TargetSpecies == "random_abundance") {
#     if (TargetGuild == "rows") coext_rows <- sample(alive_rows, 1, prob = rowSums(imatrix[alive_rows, , drop = FALSE]))
#     if (TargetGuild == "cols") coext_cols <- sample(alive_cols, 1, prob = colSums(imatrix[, alive_cols, drop = FALSE]))
#   } else if (TargetSpecies == "random_binary") {
#     if (TargetGuild == "rows") coext_rows <- sample(alive_rows, 1)
#     if (TargetGuild == "cols") coext_cols <- sample(alive_cols, 1)
#   }
#   
#   # Remove invalid indices
#   coext_rows <- coext_rows[!is.na(coext_rows) & coext_rows != 0]
#   coext_cols <- coext_cols[!is.na(coext_cols) & coext_cols != 0]
#   
#   if (length(coext_rows) > 0) imatrix[coext_rows, , drop = FALSE] <- 0
#   if (length(coext_cols) > 0) imatrix[, coext_cols, drop = FALSE] <- 0
#   
#   lost_rows <- coext_rows
#   lost_cols <- coext_cols
#   
#   # --------- CASCADE LOOP -------------------------------------
#   equilibrium <- FALSE
#   degree <- 1
#   
#   degree_table <- data.frame(
#     degree = degree,
#     guild = factor(TargetGuild, levels = c("rows", "cols")),
#     n_extinctions = 1
#   )
#   
#   while (!equilibrium) {
#     
#     remaining_rows <- rows[!is.na(rowsNA)]
#     remaining_cols <- cols[!is.na(colsNA)]
#     
#     # Stop if one guild has no remaining interactions
#     if (length(remaining_rows) == 0 || length(remaining_cols) == 0) break
#     
#     extinct_rows <- coext_rows
#     extinct_cols <- coext_cols
#     
#     rowsNA[extinct_rows] <- NA
#     colsNA[extinct_cols] <- NA
#     
#     # ---- Rows go extinct, affect columns ----
#     if (length(extinct_rows) > 0) {
#       new_coext_cols <- integer(0)
#       for (i in extinct_rows) {
#         Target <- R_val_cols[remaining_cols] * interaction_strength[i, remaining_cols, 1] > runif(length(remaining_cols))
#         new_coext_cols <- c(new_coext_cols, remaining_cols[Target])
#       }
#       coext_cols <- unique(new_coext_cols)
#       coext_cols <- coext_cols[!is.na(coext_cols) & coext_cols != 0]
#       if (length(coext_cols) > 0) imatrix[, coext_cols, drop = FALSE] <- 0
#       lost_cols <- c(lost_cols, coext_cols)
#       if (length(coext_cols) > 0) {
#         degree <- degree + 1
#         degree_when_lost_cols <- c(degree_when_lost_cols, rep(degree, length(coext_cols)))
#         degree_table[degree, ] <- data.frame(degree, "cols", length(coext_cols))
#       }
#       coext_rows <- integer(0)
#     }
#     
#     # ---- Columns go extinct, affect rows ----
#     else if (length(extinct_cols) > 0) {
#       new_coext_rows <- integer(0)
#       for (i in extinct_cols) {
#         Target <- R_val_rows[remaining_rows] * interaction_strength[remaining_rows, i, 2] > runif(length(remaining_rows))
#         new_coext_rows <- c(new_coext_rows, remaining_rows[Target])
#       }
#       coext_rows <- unique(new_coext_rows)
#       coext_rows <- coext_rows[!is.na(coext_rows) & coext_rows != 0]
#       if (length(coext_rows) > 0) imatrix[coext_rows, , drop = FALSE] <- 0
#       lost_rows <- c(lost_rows, coext_rows)
#       if (length(coext_rows) > 0) {
#         degree <- degree + 1
#         degree_when_lost_rows <- c(degree_when_lost_rows, rep(degree, length(coext_rows)))
#         degree_table[degree, ] <- data.frame(degree, "rows", length(coext_rows))
#       }
#       coext_cols <- integer(0)
#     }
#     
#     equilibrium <- (length(coext_rows) + length(coext_cols)) == 0
#   }
#   
#   # --------- OUTPUT -------------------------------------------
#   out <- list(
#     cascade_data = degree_table,
#     rows_species_data = data.frame(
#       lost_rows = lost_rows,
#       degree_of_extinction = degree_when_lost_rows
#     ),
#     cols_species_data = data.frame(
#       lost_cols = lost_cols,
#       degree_of_extinction = degree_when_lost_cols
#     )
#   )
#   
#   if (return.matrix) out$interaction_matrix <- imatrix
#   
#   return(out)
# }
