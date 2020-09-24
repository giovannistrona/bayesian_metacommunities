
## Generate perturbed grid landscape
## Input
## - n: number of patches (approximate; corrected to fit onto grid)
## - pertsize: how much the position of each patch is perturbed
## Output
## - a dataframe with column 1 the x and column 2 the y-coordinates of the patches
lscframe <- function(n, pertsize) {
    grid_coords <- seq(0, 1, l=round(sqrt(n))) ## regularly spaced points in 1D
    grid_pts <- expand.grid(V1=grid_coords, V2=grid_coords) ## 2D grid
    pert <- as.data.frame(matrix(pertsize*runif(2*n, -1, 1), n, 2)) ## random part
    lsc <- (grid_pts+pert)%%1 ## grid + random, modulus 1 (to keep in unit square)
    colnames(lsc) <- c("x", "y") ## rename columns to x and y
    return(lsc)
}

