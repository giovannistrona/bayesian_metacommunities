#Load R functions
source("get_matrix.R")

## Calculate the metapopulation capacity, the threshold condition above which
## a metapopulation persists. Defined as the leading eigenvalue of the Jacobian
## of the principal map, evaluated at 0. 
## Input:
## - spdat: species and landscape data. 
## - sp: species i. 
## Output: 
## - lambda: metapopulation capacity of species i. 
calc_metapopcapacity <- function(spdat, sp){
  M <- get_matrix(coords=cbind(spdat$x[spdat$species==sp],
                               spdat$y[spdat$species==sp]), 
                  kernel=spdat$kernel[spdat$species==sp][1],
                  xi=spdat$xi[spdat$species==sp][1]) ## sp's dispersal matrix 
  ## Jacobian of the principal map evaluated at 0:
  Eik_inverse <- ifelse(spdat$delta[spdat$species==sp]<=0, 1000,
                        -1/log(1-spdat$delta[spdat$species==sp]))
  Amat <- M*Eik_inverse
  lambda <- max(abs(eigen(Amat)$values)) ## metapopulation capacity
  return(lambda) ## return the metapopulation capacity
}
