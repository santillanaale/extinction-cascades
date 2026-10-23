# web: the bipartite network (usually a matrix of interactions).
# participant: "lower", "higher", or "both" — determines which trophic level the extinction is applied to.
# method: strategy used to pick which species to remove.
# ext.row / ext.col: custom extinction orders (only used with "external" method).

modified.extinction <- function(web, participant = "both", method = "random", 
                                ext.row = NULL, ext.col = NULL, reverse = FALSE) {
  
  ## ---- Input Matching and Validation ----
  partis <- c("lower", "higher", "both") # Matches participant argument to one of the expected values.
  partis.match <- pmatch(participant, partis) # pmatch() returns the index (1, 2, or 3).
  if (is.na(partis.match)) 
    stop("Choose participant: lower/higher/both.\n") # Throws error if no match found.
  
  meths <- c("random", "abundance", "degree", "external") # Matches method argument to one of the expected values.
  meths.match <- pmatch(method, meths) # pmatch() returns the index (1, 2, 3, or 4).
  if (is.na(meths.match)) 
    stop("Choose extinction method: random/abundance/degree.\n") # Throws error if no match found.
  
  # Get Matrix Dimensions
  nr <- NROW(web) # Number of rows (lower-level species)
  nc <- NCOL(web) # Number of columns (higher-level species)
  
  ## ---- METHOD: RANDOM ----
  # If participant = "both" and method = "random", pick randomly between lower (1) and higher (2).
  if (partis.match == 3 & meths.match == 1) 
    partis.match <- sample(2, 1)
  # Randomly selects a row (rexcl) or column (cexcl) and removes it (sets to 0), depending on participant.
  if (meths.match == 1) {
    rexcl <- sample(nr, 1)
    cexcl <- sample(nc, 1)
    if (partis.match == 1) 
      web[rexcl, ] <- 0
    if (partis.match == 2) 
      web[, cexcl] <- 0
  }
  
  ## ---- METHOD: ABUNDANCE (with reverse option) ----
  if (meths.match == 2) {
    # Randomly shuffles rows and columns to avoid bias from ties.
    web <- web[sample(1:nrow(web)), sample(1:ncol(web)), drop = FALSE]
    # Then calculates the order of rows and columns based on total abundance (sum of interactions).
    rseq <- order(rowSums(web))
    cseq <- order(colSums(web))
    # Reverse order if requested (remove most abundant instead of least)
    if (reverse) {
      rseq <- rev(rseq)
      cseq <- rev(cseq)
    }
    # Removes the least abundant species (row or column with the smallest sum).
    if (partis.match == 1) 
      web[rseq[1], ] <- 0
    if (partis.match == 2) 
      web[, cseq[1]] <- 0
    # f both trophic levels are allowed to go extinct ("both"), it compares minimum abundance between the two levels and 
    # removes the one with the lower abundance. If tied, choose randomly.
    if (partis.match == 3) {
      if (min(rowSums(web)) < min(colSums(web))) {
        web[rseq[1], ] <- 0
      }
      else {
        if (min(rowSums(web)) > min(colSums(web))) {
          web[, cseq[1]] <- 0
        }
        else {
          if (sample(2, 1) == 1) {
            web[rseq[1], ] <- 0
          }
          else {
            web[, cseq[1]] <- 0
          }
        }
      }
    }
  }
  
  ## ---- METHOD: DEGREE (with reverse option) ----
  # For lower level, compute degree (number of links), find the row with the most connections, and remove it.
  if (meths.match == 3) {
    if (partis.match == 1) {
      sequ <- rowSums(web > 0)
      order_seq <- order(sequ)
      if (reverse) order_seq <- rev(order_seq)
      which.ex <- which(sequ == sequ[order_seq[1]])
      ex <- if (length(which.ex) > 1) 
        sample(which.ex, 1) 
      else which.ex
      web[ex, ] <- 0
    }
  # For higher level: same logic, remove the most connected column.
    if (partis.match == 2) {
      sequ <- colSums(web > 0)
      order_seq <- order(sequ)
      if (reverse) order_seq <- rev(order_seq)
      which.ex <- which(sequ == sequ[order_seq[1]])
      ex <- if (length(which.ex) > 1) sample(which.ex, 1) else which.ex
      web[, ex] <- 0
    }
  }
  
  ## ---- METHOD: EXTERNAL ----
  # Uses the provided extinction sequence (ext.row or ext.col).
  # Always removes the first species in the sequence.
  if (meths.match == 4) {
    rseq <- ext.row
    cseq <- ext.col
    if (partis.match == 1) 
      web[rseq[1], ] <- 0
    if (partis.match == 2) 
      web[, cseq[1]] <- 0
  }
  
  # Returns the updated interaction matrix with one species removed.
  return(web)
}