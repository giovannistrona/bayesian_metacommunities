## Calculate the principal map. 
## Input:
## - p: species i's persistence probability, 
## - pars: list of parameters
##         - M: dispersal matrix, 
##         - delta: extinction probability based on Bayesian network,
##         - cap: maximum inverse extinction rate (useful to set to nreps)
## Output:
## - i's principal map g_i,k = C_i,k / E_i,k, with C_i,k = sum_l M_i,kl * p_i,l
##   and E_i,k = -log(1-delta_i,k) 
principalmap <- function(p, pars) {
  Eik_inverse <- ifelse(pars$delta<=0, pars$cap, -1/log(1-pars$delta))
  return((pars$M%*%p)*Eik_inverse)
}
