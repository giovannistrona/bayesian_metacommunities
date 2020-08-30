## Calculate spectral properties of principal map's Jacobian at 0, for given sp.
## Input:
## - spdat: species and landscape data
## - focal_sp: focal species for which patch values are calculated
## Output:
## - a list with elements:
##   * pvec: vector with patch values for species i.   
calc_spectral <- function(spdat, focal_sp){
  M <- get_matrix(coords=cbind(spdat$x[spdat$species==focal_sp], 
                               spdat$y[spdat$species==focal_sp]), 
                  kernel=spdat$kernel[spdat$species==focal_sp][1],
                  xi=spdat$xi[spdat$species==focal_sp][1]) ## dispersal matrix
  ## Jacobian of the principal map evaluated at 0:
  Eik <- -log(1-spdat$delta[spdat$species==focal_sp])
  Amat <- M/Eik
  ## approximate relative patch values w/ leading eigenvector of the Jacobian 
  esys <- eigen(Amat) ## eigensystem of Amat
  lew <- max(abs(esys$values)) ## dominant eigenvalue
  ind <- which.max(abs(esys$values)) ## index of dominant eigenvalue
  rev <- Re(esys$vectors[,ind]) ## dominant right eigenvector
  lev <- Re(eigen(t(Amat))$vectors[,ind]) ## dominant left eigenvector
  smat <- lev%o%rev/as.numeric(lev%*%rev) ## sensitivity matrix (diagonal = patch values)
  return(list(lambda=lew, smat=smat)) ## return dominant ew and sens. matrix
}