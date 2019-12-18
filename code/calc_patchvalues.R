#Load R functions
source("get_matrix.R")

## Calculate patch values (relative patch values or betweenness centrality). 
## Input: 
## - spdat: species and landscape data
## - focal_sp: focal species for which patch values are calculated
## Output: 
## - pvec: vector with patch values for species i.   
calc_patchvalues <- function(spdat, focal_sp){
  M <- get_matrix(coords=cbind(spdat$x[spdat$species==focal_sp], 
                               spdat$y[spdat$species==focal_sp]), 
                  kernel=spdat$kernel[spdat$species==focal_sp][1],
                  xi=spdat$xi[spdat$species==focal_sp][1]) ## dispersal matrix
  ## Jacobian of the principal map evaluated at 0:
  Eik_inverse <- ifelse(spdat$delta[spdat$species==focal_sp]<=0, 1000,
                        -1/log(1-spdat$delta[spdat$species==focal_sp]))
  Amat <- M*Eik_inverse
  ## approximate relative patch values w/ leading eigenvector of the Jacobian 
  lev <- eigen(Amat)$vectors[,which.max(abs(eigen(Amat)$values))]
  pvec <- (lev*lev)/(as.numeric(lev%*%lev)) ## get relative patch values
  return(pvec) ## return vector with relative patch values
}