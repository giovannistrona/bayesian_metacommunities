## Monte Carlo evaluation of Bayesian network for a single species
## Input:
## - sp: name (as character string) of focal species
## - spdat: species information table
## - A: adjacency matrix
## - nreps: number of iterations for the Monte Carlo sim of the Bayesian network
## - alpha: first parameter of the Beta distribution
## - beta: second parameter of the Beta distribution
## Output:
## - a vector of length N (number of patches), each containing the marginal
##   extinction probability in patch k of species sp
bayesian_ext <- function(sp, spdat, A, nreps, alpha, beta) {
  N <- length(unique(spdat$patch)) ## no of patches
  marginal <- rep(0, N) ## marginal ext prob of sp in all patches
  for (k in 1:N) { ## for each patch:
    be <- spdat$pi[spdat$species==sp][k] ## baseline ext prob of sp in patch k
    deltas <- spdat %>% ## marginal ext probs of all species in kth patch
      filter(patch==unique(spdat$patch)[k]) %>% 
      mutate(d=1-(1-delta)*p) %>% 
      pull(d)
    for (i in 1:nreps) {
      extinct <- 1*(runif(nrow(A))<deltas)
      frac <- as.numeric(A[sp,] %*% extinct / sum(A[sp,]))
      Pext <- be + (1 - be) * pbeta(frac, alpha, beta) ## prob of extinction
      if (runif(1)<Pext) marginal[k] <- marginal[k] + 1
    }
  }
  return(marginal/nreps)
}
