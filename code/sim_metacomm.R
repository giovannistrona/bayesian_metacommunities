#Load R functions
source("get_matrix.R")
source("bayesian_ext.R")
source("iteratemap.R")
source("calc_patchvalues.R")
source("calc_metapopcapacity.R")

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
  for (sp in rownames(A)) {
    M <- get_matrix(coords=cbind(spdat$x[spdat$species==sp], 
                                 spdat$y[spdat$species==sp]), 
                    kernel=spdat$kernel[spdat$species==sp][1],
                    xi=spdat$xi[spdat$species==sp][1]) ## dispersal matrix  
    if (sum(A[sp,])==0) { ## basal species
      spdat$delta[spdat$species==sp] <- spdat$pi[spdat$species == sp]
    } else { ## consumer species
      spdat$delta[spdat$species==sp] <- bayesian_ext(sp, spdat, A,
                                                     nreps, alpha, beta)
    }
    converge <- FALSE ## test for convergence: iterate 10 times, compare with
    while (!converge) { ## result before, and stop if they are close enough
      for (iter in 1:10) p_next <- iteratemap(
        p=spdat$p[spdat$species==sp],
        list(M=M, cap=nreps,
             delta=spdat$delta[spdat$species==sp]))
      if (sum((spdat$p[spdat$species==sp]-p_next)^2)<1e-8) converge <- TRUE
      spdat$p[spdat$species==sp] <- p_next
    }
    ## calculate sp's metapopultion capacity
    spdat$lambda[spdat$species==sp] <- calc_metapopcapacity(spdat, sp)
    ## calculate relative patch values
    spdat$pvalues[spdat$species==sp] <- calc_patchvalues(spdat, sp)
  }
  return(spdat)    
}
