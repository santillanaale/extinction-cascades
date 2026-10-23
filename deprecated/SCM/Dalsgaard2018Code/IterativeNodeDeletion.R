IterativeNodeDeletion <- function(network, anim.r, plant.r, removal.taxa, remove.order, silent = TRUE){
  #--------- LOADING IN REQUIRED PACKAGES AND FUNCTIONS ---------------------------
  require(bipartite)
  source("SCM/VieraAlmeida2015RFunctions/netcascade (April 2014).R", local = TRUE)
  
  #--------- ERROR CHECKS ---------------------------
  if(removal.taxa != "animal" & removal.taxa != "plant"){stop("'removal.taxa' must be either 'animal' or 'plant'")}
  if(any(rownames(network) %in% colnames(network))){stop("Species cannot occur in more than one trophic level")}
  if(!remove.order %in% c("random","from strongest","from weakest","from most connected","from least connected")){stop("'remove.order' must be set as either 'random', 'from strongest' or 'from weakest'")}
  
  #--------- DEFINING SOME VARIABLES AND SETTING UP RESULTS ARRAYS ---------------------------
  counter <- 0
  
  web <- network
  
  # if the web has no names, give them names:
  if(is.null(rownames(web)) == TRUE){rownames(web) <- paste0("L", seq.int(nrow(web)))}
  if(is.null(colnames(web)) == TRUE){colnames(web) <- paste0("H", seq.int(ncol(web)))}
  
  # number of animals and plants in the initial network
  nanim <- nrow(web)
  nplant <- ncol(web)
  
  # R values for animals - must be either a single value, applied to all species, or a vector of values of length equal to the number of species
  if(length(anim.r) == 1){
    ranim <- rep(anim.r, nrow(web))
  } else {
    if(length(anim.r) == nrow(web)){
      ranim <- anim.r
    } else {
      stop("anim.r must be either of length 1 or equal to the number of rows in the network")
    }
  }
  
  # R values for plants - must be either a single value, applied to all species, or a vector of values of length equal to the number of species
  if(length(plant.r) == 1){
    rplants <- rep(plant.r, ncol(web))
  } else {
    if(length(plant.r) == ncol(web)){
      rplants <- plant.r
    } else {
      stop("plant.r must be either of length 1 or equal to the number of rows in the network")
    }
  }
  
  # no dead animals or plants at start of simulation
  deadAnimals <- NULL
  deadPlants <- NULL
  
  # data.frame recording response fates
  response.fates <- as.data.frame(matrix(ncol = 3, nrow = nanim))
  colnames(response.fates) <- c("species","fate","order")
  response.fates$species <- rownames(web)
  response.fates$fate <- 1 # 1 = alive, 0 = extinct
  
  binary.web <- web
  binary.web[binary.web > 1] <- 1
  
  #--------- START OF ITERATIVE PLANT REMOVAL SIMULATION ---------------------------
  while(all(rowSums(web) == 0) == FALSE){
    counter <- counter + 1
    
    # run SCM simulation
    if(remove.order == "random"){
      sim <- netcascade(imatrix = web, ranim = ranim, rplants = rplants, deadPlants = deadPlants, deadAnimals = deadAnimals, targetGuild = "plant", target = rep(1, ncol(web)), return.matrix = TRUE) # run simulation
    } else if(remove.order == "from weakest"){
      weakest <- min(colSums(web)[colSums(web) != 0])
      weakest <- as.integer(which(colSums(web) == weakest))
      if(length(weakest) > 1){
        weakest <- sample(weakest, length(weakest), replace = FALSE)[1]
      }
      sim <- netcascade(imatrix = web, ranim = ranim, rplants = rplants, deadPlants = deadPlants, deadAnimals = deadAnimals, targetGuild = "plant", target = weakest, return.matrix = TRUE) # run simulation
    } else if(remove.order == "from strongest"){
      strongest <- as.integer(which(colSums(web) == max(colSums(web))))
      if(length(strongest) > 1){
        strongest <- sample(strongest, length(strongest), replace = FALSE)[1]
      }
      sim <- netcascade(imatrix = web, ranim = ranim, rplants = rplants, deadPlants = deadPlants, deadAnimals = deadAnimals, targetGuild = "plant", target = strongest, return.matrix = TRUE) # run simulation
    } else if(remove.order == "from most connected"){
      most.connected <- as.integer(which(colSums(binary.web) == max(colSums(binary.web))))
      if(length(most.connected) > 1){
        most.connected <- sample(most.connected, length(most.connected), replace = FALSE)[1]
      }
      sim <- netcascade(imatrix = web, ranim = ranim, rplants = rplants, deadPlants = deadPlants, deadAnimals = deadAnimals, targetGuild = "plant", target = most.connected, return.matrix = TRUE) # run simulation
    } else if(remove.order == "from least connected"){
      least.connected <- min(colSums(binary.web)[colSums(binary.web) != 0])
      least.connected <- as.integer(which(colSums(binary.web) == least.connected))
      if(length(least.connected) > 1){
        least.connected <- sample(least.connected, length(least.connected), replace = FALSE)[1]
      }
      sim <- netcascade(imatrix = web, ranim = ranim, rplants = rplants, deadPlants = deadPlants, deadAnimals = deadAnimals, targetGuild = "plant", target = least.connected, return.matrix = TRUE) # run simulation
    }
    # update web with post-simulation interaction matrix
    web <- sim$interaction_matrix
    binary.web <- web
    binary.web[binary.web > 1] <- 1
    
    # record disconnected species
    disconAnimals <- as.integer(which(rowSums(web) == 0))
    disconPlants <- as.integer(which(colSums(web) == 0))
    
    # feed extinct and disconnected species into the deadAnimals and deadPlants argument (disconnected species must be included to avoid infinite loops)
    deadAnimals <- as.integer(c(sim$lost_animals, disconAnimals))
    deadPlants <- as.integer(c(sim$lost_plants, disconPlants))
    
    # update response.fates
    if(length(sim$lost_animals) > 0){
      response.fates[response.fates$species %in% rownames(web)[sim$lost_animals],]$fate <- 0 
      response.fates[response.fates$species %in% rownames(web)[sim$lost_animals],]$order <- counter 
    }
    
    if(silent == FALSE){
      print(counter)
    }
  }
  
  #--------- OUTPUT ---------------------------
  return(list(response_fates = response.fates))
}
