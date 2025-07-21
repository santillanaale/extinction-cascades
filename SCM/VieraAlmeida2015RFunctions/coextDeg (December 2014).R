#coextMaxLvl performs simulations (n = nsims) of a single episode of primary extinction and its possible associated coextinction cascade and returns how often the extinction sequence ended in each extinction level (i.e. first level, second level, third level, etc.) 
#primary extinctions are uniform random among both groups.'r' is sampled from a interval (rlow, rup)

coextDeg <- function(imatrix,rlow,rup,nsims){
  degs <- c()
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
    profiles <- netcascade(imatrix=imatrix,ranim=ranim,rplants=rplants,targetGuild=guild,target=target)
    degs[sim] <- max(profiles[[1]]$degree)
  }
  return(degs)
}