#coextNumber performs simulations (n = nsims) of a single episode of primary extinction and its possible associated coextinction cascade and creates a frequency distribution for the total number of extinctions for each episode.
#primary extinctions are uniform random among both groups.'r' is sampled from a interval (rlow, rup)

coextNumber <- function(imatrix,rlow,rup,nsims){
  ext_counts <- c()
  for(sim in 1:nsims){
    rvalue <- runif(1,rlow,rup)
    ranim <- rep(rvalue, nrow(imatrix))
    rplants <- rep(rvalue, ncol(imatrix))
    guild <- sample(c('animal','plant'),1,F,c(nrow(imatrix),ncol(imatrix)))
    if(guild=="animal"){
      target <- rep(1,nrow(imatrix))
    }else{
      target <- rep(1,ncol(imatrix))
    }    
    sim_results <- netcascade(imatrix,ranim = ranim, rplants = rplants, targetGuild = guild, target = target)
    ext_counts[sim] <- sum(sim_results[[1]]$n_extinctions)
  }
  return(ext_counts)
}