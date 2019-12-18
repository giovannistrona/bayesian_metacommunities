## Generate data file with species and landscape information. 
## Input: 
## - webfile: file containig the edgelist, i.e. who eats whom
## - datfile: species input file containing for each species sp its
##    - pi: baseline extinction probability,
##    - xi: dispersal distance,
##    - kernel: shape of the dispersal kernel ("Gaussian", "Exponential", or
##              "Rectangular").
## - landscapefile: file containing the patch coordinates.
## Output: 
## - spdat: data file combining species and landscape information. 
generate_spdat <- function(webfile, datfile, landscapefile){
  A <- read_csv(webfile, col_types=cols()) %>% ## read landscape file
    create_adj_matrix ## create adjacency matrix from links
  landscape <- read_csv(landscapefile, col_types=cols()) %>%
    as.matrix
  info <- read_csv(datfile, col_types=cols()) ## import species information
  N <- length(landscape[,1]) ## number of patches
  spdat <- tibble(patch=integer(),
                  x=numeric(),
                  y=numeric(),
                  species=character(),
                  pi=numeric(),
                  xi=numeric(),
                  kernel=character(), 
                  TL = numeric(),
                  OI = numeric(),
                  FT = character()) 
  dat <- tibble(patch=1:N, x=landscape[,1], y=landscape[,2])
  for(sp in rownames(A)){
    dat$species <- rep(sp, N)
    dat$pi <- rep(info$pi[info$species==sp], N)
    dat$xi <- rep(info$xi[info$species==sp], N)
    dat$kernel <- rep(info$kernel[info$species==sp], N)
    dat$TL <- rep(info$TL[info$species==sp], N)
    dat$OI <- rep(info$OI[info$species==sp], N)
    dat$FT <- rep(info$FT[info$species==sp], N)
    spdat <- bind_rows(spdat, dat)
  }
  return(spdat)
}