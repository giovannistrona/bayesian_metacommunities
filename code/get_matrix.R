## Create disperal matrix M. 
## Input: - coords: patch coordinatexs x and y, 
## - kernel: type of dispersal kernel (species-specific), 
## - xi: dispersal distance (species-specific). 
## Output:
## - dispersal matrix M. 
get_matrix <- function(coords, kernel, xi) {
  distances <- as.matrix(dist(coords))
  if (kernel=="Exponential") M <- exp(-distances/xi)
  if (kernel=="Gaussian")    M <- exp(-distances^2/(2*xi^2))
  if (kernel=="Rectangular") M <- 1*(distances<xi)
  diag(M) <- 0 # remove the diagonal from the matrix
  return(unname(M))
}
