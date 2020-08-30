############################################################
# Set up for simulations 
############################################################

#Load R functions
source("create_adj_matrix.R")
source("get_feeding_type.R")
source("generate_spdat.R")
source("sim_metacomm.R")

#Helper functions
'%ni%' <- function(x,y)!('%in%'(x,y)) 

## Run simulations. 
## Input: - web: name of web file (w/ path & extension), 
## - landsc: name of landscape file (w/ path & extension)
## - spinput: name of species input file (w/ path & extension)
## - spf: focal species for patch removal ("basal" or "top")
## - alpha: 1st parameter of beta distribution (Bayesian network)
## - beta: 2nd parameter of beta distribution (Bayesian network)
## - seed: seed for random numbers 
## - rem: number of patches removed per patch loss iteration
## - nreps:  number of iterations for Bayesian network. 
## Output:
## - spdat: simulation output 

run_instance <- function(web, landsc, spinput, spf, alpha, beta, seed, rem, nreps){
  # set.seed(seed) ## set a random seed (optional, included here for reproducibility) OR NOT BECAUSE OF PATCH SELECTION !!
  scenarios <- c("random", "best", "worst") #set of scenarios
  
  ###########################################################################################
  # 1: Simulate metacommunity once in the 'full' landscape without habitat loss. 
  ###########################################################################################
  A <- create_adj_matrix(read_csv(web, col_types=cols())) ## create adjacency matrix from links
  landscape <- as.matrix(read_csv(landsc, col_types=cols())) ## read landscape file
  N <- length(landscape[,1]) ## number of patches in the landscape
  
  spdat0 <- generate_spdat(web, spinput, landsc) %>% ## combine species and landscape information
    add_column(delta=0, lambda=0, p=0.5, pvalues=0)
    
  spdat0 <- sim_metacomm(A, spdat0, nreps, alpha, beta) %>% ## simulate metacommunity
    add_column(prem=0,
               web=web,
               spinput=spinput,
               landsc=landsc,
               spf=spf,
               seed=seed,
               alpha=alpha,
               beta=beta,
               scenario="all") 
  
  spdat <- spdat0
  
  ###########################################################################################
  # 2: Test whether one of the following conditions is met. If so, we stop here; 
  # else we continue and simulate the habitat loss scenarios. 
  ###########################################################################################
  
  flag <- FALSE 
  if ((all(spdat0 %>% filter(FT!="basal") %>% pull(lambda) %>% unique < 1)) || ## stop if no consumer can persist in the long-term
      (max(spdat0$p)<(1e-10)) || ## stop if persistence too low 
      ((spf!="basal") & (spdat0 %>% filter(species==rownames(A)[nrow(A)]) %>% pull(lambda) %>% unique<1))){ ## stop if top species is spf and extinct 
    flag <- TRUE
  }
  
  focal_sp <- ifelse(spf=="basal", rownames(A)[1], rownames(A)[nrow(A)]) ## name of focal species for patch removal 
  
  ###################################################################################################
  # 3: Simualte metacommunities with habitat loss (if none of the above listed conditions applies).
  ###################################################################################################
  
  if(!flag){
    for (scen in scenarios) {
      N0 <- N ## number of patches in landscape initially
      j <- 0
      rem0 <- rem
      pselec <- integer(0) ## initialize vector for patches selected for removal
      crashed <- FALSE ## test if populations are still extant
      spdat1 <- spdat0 %>% mutate(scenario=scen)
      while (!crashed) {
        j <- j + rem0
        if(scen %in% "best"){ ## best-case scenario; 
          ## we remove the lowest-value patches first 
          pselec <- spdat1 %>% ## vector of patches selected for removal
            filter(species==focal_sp) %>%
            top_n(-rem0, pvalues) %>% ## negative to select bottom rem0 rows
            pull(patch)
        }
        if(scen %in% "worst"){ ## worst-case scenario; 
          ## we remove the highest-value patches first
          pselec <- spdat1 %>%  ## vector of patches selected for removal
            filter(species==focal_sp) %>%  
            top_n(rem0, pvalues) %>%
            pull(patch) 
        }
        if(scen=="random"){ ## random scenario; 
          ## we remove patches at random
          pselec <- sample(unique(spdat1$patch), rem0) ## vector of patches selected for removal
        }
        pselec <- ifelse(length(pselec)<=rep(rem,length(pselec)), pselec, pselec[1:rem]) ## make sure pselec has a max. length of rem
        spdat_rem_patches <- spdat1 %>%
          filter(patch %ni% pselec) %>% ## remove selected patches
          mutate(prem=j)
        spdat1 <- sim_metacomm(A, spdat_rem_patches, nreps, alpha, beta)
        spdat <- rbind(spdat, spdat1)
        N0 <- length(unique(spdat1$patch)) ## number of patches in landscape
        if (N0<=rem0) rem0 <- 1
        if ((spf!="basal") & (spdat1 %>% filter(species==rownames(A)[nrow(A)]) %>% 
                              pull(lambda) %>% unique<1)){ ## stop if top species is spf and extinct
          crashed <- TRUE
        } 
        if((all((spdat1 %>% filter(FT!="basal") %>% pull(lambda) %>% unique)<1)) || ## stop if metapopulation capacity for all consumer species is below 1; 
           (max(spdat1$p)<=(1e-10)) || ## stop if overall persistence is too low; 
           (N0<=2)) {  ## or two or less patches remain in the landscape
          crashed <- TRUE
        }
      }
    }
  }
  return(spdat)
}

