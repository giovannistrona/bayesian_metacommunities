## Randomly draw patch coordinates. 
## Input:
## - N: number of patches in the landscape, 
## - dim: dimension of the landscape (1, 2 or 3). 
## Output: 
## - matrix with N rows and dim columns. 
generate_coordinates <- function(N, dim) {
  return(matrix(runif(N*dim, 0, 1), nrow=N, ncol=dim))
}