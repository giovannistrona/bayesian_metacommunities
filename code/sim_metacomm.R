#Load R functions
source("get_matrix.R")
source("calc_spectral.R")

#Load C++ functions
sourceCpp("c_functions.cpp")

## Obtain metacommunity equilibrium state.
## Input:
## - A: adjacency matrix
## - spdat: species information table
## - nreps: number of iterations for the Monte Carlo sim of the Bayesian network
## - alpha: first parameter of the Beta distribution
## - beta: second parameter of the Beta distribution
## Output:
## - a species information table (spdat) updated with the correct equilibrium
##   marginal extinction probabilities delta_i^k and patch occupancies p_i^k
sim_metacomm <- function(A, spdat, nreps, alpha, beta) {
  N <- length(unique(spdat$patch)) ## no of patches
  marginal <- rep(0, N) ## marginal ext prob of sp in all patches
  for (sp in rownames(A)) { ## calculate deltas of each sp from Bayesian NW:
    if (sum(A[sp,])==0) { ## basal species
      marginal <- spdat$pi[spdat$species==sp]
    } else { ## consumer species
      for (k in 1:N) { ## for each patch:
        be <- spdat$pi[spdat$species==sp][k] ## baseline ext in patch k
        deltas <- spdat %>% ## marginal ext of all species in patch k
          filter(patch==unique(spdat$patch)[k]) %>% 
          mutate(d=1-(1-delta)*p) %>% 
          pull(d)
        marginal[k] <- getMarginal(A[sp,], deltas, be, alpha, beta, nreps)
      }
      marginal <- marginal/nreps
    }
    spdat$delta[spdat$species==sp] <- marginal ## assign ext prob to sp
    M <- get_matrix(coords=cbind(spdat$x[spdat$species==sp], 
                                 spdat$y[spdat$species==sp]), 
                    kernel=spdat$kernel[spdat$species==sp][1],
                    xi=spdat$xi[spdat$species==sp][1]) ## disp matrix of sp
    ## solve for equilibrium patch occupancies of sp numerically
    spdat$p[spdat$species==sp] <- itersol(spdat$p[spdat$species==sp], M,
                                          spdat$delta[spdat$species==sp])
    ## calculate sp's metapopultion capacity & patch values
    spectr <- calc_spectral(spdat, sp)
    spdat$lambda[spdat$species==sp] <- spectr$lambda
    spdat$pvalues[spdat$species==sp] <- diag(spectr$smat)
  }
  return(spdat)
}


