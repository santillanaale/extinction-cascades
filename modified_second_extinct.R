# The function is defined as second.extinct, which takes several arguments:
#   web: A bipartite network object.
# participant: Specifies whether to extinct the higher trophic level ("higher"), lower trophic level ("lower"), or both ("both").
# method: The method used for extinction. Can be one of "abundance", "random", "degree", or "external".
# nrep: The number of iterations to perform if the method is "external". Default is 10.
# details: A logical indicating whether to return detailed output. Default is FALSE.
# ext.row and ext.col: Vectors providing external extinction sequences for the higher and lower trophic levels, respectively.

modified.second.extinct <- function (web, participant = "higher", method = "abun", nrep = 10, 
          details = FALSE, ext.row = NULL, ext.col = NULL, reverse = FALSE) 
{
  # ---- Input checks (unchanged) ----
  if (participant == "both" & method == "external") 
    stop("Sorry, that won't work. When you specify the sequence, you have to choose one of the two levels. 'both' won't work.")
  if (!is.null(ext.row) & length(ext.row) != NROW(web)) 
    stop("The length of the external row vector is different from the numbers of rows in the network!")
  if (!is.null(ext.col) & length(ext.col) != NCOL(web)) 
    stop("The length of the external col vector is different from the numbers of cols in the network!")
  if (participant == "higher" & method == "external" & is.null(ext.col)) 
    stop("You need to provide an external sequence of extinction for the higher trophic level!")
  if (participant == "lower" & method == "external" & is.null(ext.row)) 
    stop("You need to provide an external sequence of extinction for the lower trophic level!")
  
  # ---- Nested function: runs a single replicate ----
  # INTERNAL FUNCTION one.second.extinct: This internal function handles the actual extinction process.
  
  # This function performs a single extinction simulation. It's nested so the outer function can call it multiple times.
  one.second.extinct <- function(web, participant, method, ext.row, ext.col, reverse) {
    dead <- matrix(nrow = 0, ncol = 3) # It initializes an empty matrix dead to store extinct nodes.
    colnames(dead) <- c("no", "ext.lower", "ext.higher") # no: step number, ext.lower: surviving lower-level species, ext.higher: surviving higher-level species
    m2 <- web # It creates a working copy of the web
    i <- 1 # i is step counter.
    #The function enters a loop where it calculates extinct nodes based on the specified method:
    repeat {
      # For each iteration, it calls extinction on the current network with the given participant and method via the extinction() function (not defined here).
      # 👇 Pass reverse to extinction()
      n <- modified.extinction(m2, participant = participant, method = method, 
                               ext.row = ext.row, ext.col = ext.col, reverse = reverse)
      
      # Updates m2 with the modified network. Logs step number and survivors.
      dead <- rbind(dead, c(i, attributes(m2 <- empty(n, count = TRUE))$empty))
      
      # Stops if the remaining web is too small to continue meaningful extinction simulation.
      if (participant == "lower" & NROW(m2) < 2) # If the participant is "lower" and the number of nodes in the current network drops below 2, it breaks the loop.
        break
      if (participant == "higher" & NCOL(m2) < 2) # Similarly for the "higher" participant.
        break
      if (participant == "both" & min(dim(m2)) < 2) 
        break
      if (any(dim(n) == 1)) 
        break
      
      # If the method is "external", the function updates the custom extinction order (reduces indices due to species removals).
      if (method == "external") {
        ext.col[ext.col > ext.col[1]] <- ext.col[ext.col > 
                                                   ext.col[1]] - 1
        ext.row[ext.row > ext.row[1]] <- ext.row[ext.row > 
                                                   ext.row[1]] - 1
        ext.row <- ext.row[-1]
        ext.col <- ext.col[-1]
      }
      i <- i + 1 # next iteration
    }
    
    # Post-loop cleanup
    dead2 <- rbind(dead, c(NROW(dead) + 1, NROW(m2), NCOL(m2)))
    
    # Special case: fix inconsistency in the "degree"-based extinction pattern.
    if (participant == "lower" & method == "degree") {
      if (length(table(dead[, 2])) > 1) 
        dead2[, 2] <- 1
    }
    
    if (nrow(dead) + 1 != nrow(dead2)) 
      stop("PANIC! Something went wrong with the extinct sequence! Please contact the author to fix this!!")
    
    # Consistency check
    if (participant == "lower") supposed.length <- NROW(web)
    if (participant == "higher") supposed.length <- NCOL(web)
    if (participant == "both") supposed.length <- NROW(dead2)
    
    # If output is too short, pad with zeros
    if (NROW(dead2) != supposed.length) {
      missing <- supposed.length - NROW(dead2)
      addit1 <- (NROW(dead2) + 1):(NROW(dead2) + missing)
      addit2n3 <- rep(0, times = missing)
      dead2 <- rbind(dead2, as.matrix(data.frame(addit1, addit2n3, addit2n3)))
    }
    
    return(dead2)
  }
  
  # ---- Main logic ----
  if (is.vector(method)) 
    sequence <- method
  
  # If method is deterministic (abundance, degree, external), run once.
  if (pmatch(method, c("abundance", "random", "degree", "external")) %in% c(1, 3, 4)) {
    out <- one.second.extinct(web = web, participant = participant, 
                              method = method, ext.row = ext.row, ext.col = ext.col,
                              reverse = reverse)
  } else {
  # Otherwise (e.g., method = "random"), run multiple replicates.
  o <- replicate(nrep, one.second.extinct(web = web, participant = participant, 
                                          method = method, ext.row = ext.row, 
                                          ext.col = ext.col, reverse = reverse), 
                   simplify = FALSE)
  # If details = TRUE: return full list of results.
    if (details) {
      out <- o
    }
  # Else: return average extinction curve across replicates.
    else {
      lengths <- sapply(o, nrow)
      z <- o[[which.max(lengths)]]
      z[, 2:3] <- 0
      for (k in 1:length(o)) {
        nr <- nrow(o[[k]])
        z[1:nr, ] <- z[1:nr, ] + o[[k]]
        rm(nr)
      }
      out <- z/length(o)
      out[, 1] <- 1:max(lengths)
    }
  }
  
  # ---- Final Wrapping ----
  
  # Set class and attributes before returning the result.
  class(out) <- "bipartite"
  attr(out, "exterminated") <- c("both", "lower", "higher")[pmatch(participant, 
                                                                   c("both", "lower", "higher"))]
  out
}
# <bytecode: 0x0000025222977b48>
#   <environment: namespace:bipartite>