#Load R functions
source("principalmap.R")

## Calculate the iteration map based on the principal map. 
## Input: -p: species i's perstistence probablity,  
## - pars: list of parameters
##         - M: dispersal matrix, 
##         - delta: extinction probability based on Bayesian network). 
##         - cap: maximum inverse extinction rate (useful to set to nreps)
## Output:
## - i's persistence prob. on each patch k. 
iteratemap <- function(p, pars) {
  return(principalmap(p, pars)/(principalmap(p, pars)+1))
}
