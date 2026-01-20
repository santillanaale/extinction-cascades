#Pires et al 2020
#adapted from Vieira and Almeida-Neto 2015; doi: 10.1111/ele.12394
#simulates extinction cascades by removing a target species
# The R arguments (R_rows, R_cols) can be a character string specifying an option (e.g. “normal”, “VAlow”, etc), or a numeric vector of length 1 or length == number of rows or columns. If a single number is given, this number is applied to all species in that guild.
# TargetGuild options: c(“rows”, “cols”, “random_richness”, “random_binary”) (VA2015 used random_richness)
# TargetSpecies options: numeric (1:nrow or 1:ncol), “random_abundance”, “random_binary”
#=========================================================================================================
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
  
  #---------ARGUMENT CHECKS-----------------------------------
  if(class(imatrix) != "matrix" || (nrow(imatrix) + ncol(imatrix)) < 3){
    stop("‘imatrix’ must be an object of class ‘matrix’, with at least three species")
  }
  
  if(length(R_val_rows) != nrow(imatrix)){
    stop("The length of vector ‘R_rows’ must be equal to number of rows (i.e. species in guild) in ‘imatrix’")
  }
  
  if(length(R_val_cols) != ncol(imatrix)){
    stop("The length of vector ‘R_cols’ must be equal to number of columns in ‘imatrix’")
  }
  
  if((TargetGuild %in% c("rows","cols","random_richness","random_binary")) == FALSE){
    stop("Invalid target guild for primary extinction. Valid targets guilds are “rows”, “cols”, “random_richness”, and “random_binary”")
  }
  
  if(is.numeric(TargetSpecies) == FALSE && !(TargetSpecies %in% c("random_abundance","random_binary"))){
    stop("Invalid value for the “TargetSpecies” argument. You may specify a single species by entering its row or column number or you may use a vector of relative probabilities for all species in the Target guild.")
  }
  
  if(!is.null(extinct_cols) && class(extinct_cols) != "integer"){
    stop("extinct_cols must be either NULL or an integer vector specifying the column numbers of species considered to be extinct on the original matrix")
  }
  
  if(!is.null(extinct_rows) && class(extinct_rows) != "integer"){
    stop("extinct_rows must be either NULL or an integer vector specifying the row numbers of species considered to be extinct on the original matrix")
  }
  
  #---------DEFINING SOME OBJECTS---------------------------
  nrows <- nrow(imatrix)
  ncols <- ncol(imatrix)
  cols <- 1:ncols
  rows <- 1:nrows
  colsNA <- 1:ncols
  rowsNA <- 1:nrows
  colsNA[extinct_cols] <- NA
  rowsNA[extinct_rows] <- NA
  
  degree_when_lost_cols <- c()
  degree_when_lost_rows <- c()
  
  #----------CALCULATING DEPENDENCE MATRICES-------------------
  interaction_strength <- array(0, dim = c(nrows, ncols, 2))
  
  #matrix of columns’ dependence on each row
  interaction_strength[,,1] <- t(t(imatrix)/rowSums(t(imatrix)))
  
  #matrix of rows’ dependence on each column
  interaction_strength[,,2] <- imatrix/rowSums(imatrix)
  
  #-----------CHOOSING TARGET SPECIES FOR PRIMARY EXTINCTION---
  coext_rows <- c()
  coext_cols <- c()
  if(is.numeric(TargetSpecies)){
    if(length(TargetSpecies) == 1){
      if(TargetGuild == "rows"){
        if(TargetSpecies %in% extinct_rows){
          stop("Specified target species for the primary extinction is already extinct")
        }
        coext_rows <- TargetSpecies
        degree_when_lost_rows <- 1
      }
      if(TargetGuild == "cols"){
        if(TargetSpecies %in% extinct_cols){
          stop("Specified target species for the primary extinction is already extinct")
        }
        coext_cols <- TargetSpecies
        degree_when_lost_cols <- 1
      }
    } else {
      nspecies <- switch(TargetGuild, rows = nrows, cols = ncols)
      if(length(TargetSpecies) == nspecies){
        if(TargetGuild == "rows"){
          alive_rows <- rows[!is.na(rowsNA)]
          coext_rows <- sample(c(alive_rows,0), 1, prob = c(TargetSpecies[!is.na(rowsNA)],0))
          degree_when_lost_rows <- 1
        }
        if(TargetGuild == "cols"){
          alive_cols <- cols[!is.na(colsNA)]
          coext_cols <- sample(c(alive_cols,0), 1, prob = c(TargetSpecies[!is.na(colsNA)],0))
          degree_when_lost_cols <- 1
        }
      } else {
        stop("Length of “TargetSpecies” must be 1 (specifying a single species within the Target guild) or else be equal to the number of species in the Target guild (specifying probabilities of primary extinction for each species in the Target guild)")
      }
    }
  } else if(TargetSpecies == "random_abundance"){
    if(TargetGuild == "rows"){
      alive_rows <- rows[!is.na(rowsNA)]
      coext_rows <- sample(x = alive_rows, size = 1, replace = FALSE,
                           prob = rowSums(imatrix[alive_rows,]))
      degree_when_lost_rows <- 1
    }
    if(TargetGuild == "cols"){
      alive_cols <- cols[!is.na(colsNA)]
      coext_cols <- sample(x = alive_cols, size = 1, replace = FALSE,
                           prob = colSums(imatrix[,alive_cols]))
      degree_when_lost_cols <- 1
    }
  } else if(TargetSpecies == "random_binary"){
    if(TargetGuild == "rows"){
      alive_rows <- rows[!is.na(rowsNA)]
      coext_rows <- sample(x = alive_rows, size = 1, replace = FALSE, prob = NULL)
      degree_when_lost_rows <- 1
    }
    if(TargetGuild == "cols"){
      alive_cols <- cols[!is.na(colsNA)]
      coext_cols <- sample(x = alive_cols, size = 1, replace = FALSE, prob = NULL)
      degree_when_lost_cols <- 1
    }
  } else {
    stop("Could not evaluate your non-numeric value for the “TargetSpecies” argument.")
  }
  
  imatrix[coext_rows, ] <- 0
  imatrix[, coext_cols] <- 0
  lost_rows <- coext_rows
  lost_cols <- coext_cols
  
  #-------------------CASCADE LOOP---------------------------
  equilibrium <- FALSE
  degree <- 1
  degree_table <- data.frame(degree, guild=factor(TargetGuild,levels=c("rows","cols")), n_extinctions=1)
  
  while(equilibrium == FALSE){
    extinct_rows <- coext_rows
    extinct_cols <- coext_cols
    colsNA[extinct_cols] <- NA
    rowsNA[extinct_rows] <- NA
    
    remaining_rows <- rows[!is.na(rowsNA)]
    remaining_cols <- cols[!is.na(colsNA)]
    
    if(length(extinct_rows) > 0){
      for(i in 1:length(extinct_rows)){
        Target <- R_val_cols[remaining_cols]*interaction_strength[extinct_rows[i],remaining_cols,1] > runif(length(remaining_cols))
        coext_cols <- c(coext_cols, remaining_cols[Target])
      }
      coext_rows <- c()
      coext_cols <- unique(coext_cols)
      colsNA[coext_cols] <- NA
      lost_cols <- c(lost_cols, coext_cols)
      imatrix[, coext_cols] <- 0
      
      for(i in 1:ncols){
        if(sum(imatrix[,i])==0){
          interaction_strength[,i,1] <- 0
        } else {
          interaction_strength[,i,1] <- imatrix[,i]/sum(imatrix[,i])
        }
      }
      
      if(length(coext_cols) > 0){
        degree <- degree + 1
        degree_when_lost_cols <- c(degree_when_lost_cols, rep(degree,length(coext_cols)))
        degree_table[degree,] <- data.frame(degree,"cols",length(coext_cols))
      }
      
    } else {
      for(i in 1:length(extinct_cols)){
        Target <- R_val_rows[remaining_rows]*interaction_strength[remaining_rows,extinct_cols[i],2] > runif(length(remaining_rows))
        coext_rows <- c(coext_rows, remaining_rows[Target])
      }
      coext_cols <- c()
      coext_rows <- unique(coext_rows)
      lost_rows <- c(lost_rows, coext_rows)
      rowsNA[coext_rows] <- NA
      imatrix[coext_rows,] <- 0
      
      for(i in 1:nrows){
        if(sum(imatrix[i,])==0){
          interaction_strength[i,,2] <- 0
        } else {
          interaction_strength[i,,2] <- imatrix[i,]/sum(imatrix[i,])
        }
      }
      
      if(length(coext_rows) > 0){
        degree <- degree + 1
        degree_when_lost_rows <- c(degree_when_lost_rows, rep(degree,length(coext_rows)))
        degree_table[degree,] <- data.frame(degree,"rows",length(coext_rows))
      }
    }
    
    equilibrium <- equilibrium + (length(coext_cols)+length(coext_rows))==0
  }
  
  #-------------------OUTPUT---------------------------
  if(length(lost_rows) > 0){
    spp_data_rows <- data.frame(lost_rows = lost_rows, degree_of_extinction = degree_when_lost_rows)
  } else {
    spp_data_rows <- "No rows were lost"
  }
  
  if(length(lost_cols) > 0){
    spp_data_cols <- data.frame(lost_cols = lost_cols, degree_of_extinction = degree_when_lost_cols)
  } else {
    spp_data_cols <- "No columns were lost"
  }
  
  if(return.matrix == TRUE){
    return(
      list(
        interaction_matrix = imatrix,
        lost_rows = lost_rows,
        lost_cols = lost_cols,
        cascade_data = degree_table,
        rows_species_data = spp_data_rows,
        cols_species_data = spp_data_cols
      )
    )
  } else {
    return(
      list(
        cascade_data = degree_table,
        rows_species_data = spp_data_rows,
        cols_species_data = spp_data_cols
      )
    )
  }
}
