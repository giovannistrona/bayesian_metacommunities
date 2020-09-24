## Randomly draw patch coordinates. 
## Input:
## - N: number of patches in a 2D-landscape, 
## Output: 
## - matrix with N rows and 2 columns. 
generate_coordinates <- function(N) {
  return(matrix(runif(N*2, 0, 1), nrow=N, ncol=2))
}