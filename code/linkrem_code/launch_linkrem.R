
#Load packages
require(Rcpp, quietly=TRUE)
require(NetIndices, quietly=TRUE)
require(igraph, quietly=TRUE)
require(tidyverse, quietly=TRUE)

#Load R functions 
source("../generate_coordinates.R")
source("../create_adj_matrix.R")
source("../generate_spdat.R")
source("../get_matrix.R")
source("sim_metacomm_linkrem.R")
source("calc_spectral_linkrem.R")

#Load C++ functions
sourceCpp("../c_functions.cpp")

#Simulation options
linkremdir <- "../../data/linkrem_data/"
landscs <- Sys.glob("../../data/landscapes/landscape_*.csv")[1] #import landscapes (or else generate them, see below)
spinputs <- Sys.glob(paste0(linkremdir, "species_input/web*.csv"))
spfs <- "basal" #focal species for habitat loss; can either be 'basal' or 'top'
alpha <- c(1, 5) #1st parameter of beta distribution (Bayesian network) 
beta <- c(1, 5) #2nd parameter of beta distribution (Bayesian network) 
nreps <- 1000 #number of iterations for the Bayesian network

#Optional: generate landscapes if needed
# N <- 300 #number of patches
# landscape <- generate_coordinates(N) %>% as_tibble %>% rename(x=V1, y=V2) %>% #generate and
#     write_csv(paste0(linkremdir,"/landscapes/landscape_001.csv")) #save landscape
# landscs <-  Sys.glob(paste0(linkremdir,"/landscapes/landscape_*.csv"))

#Generate output directory
resdir <- "../../data/linkrem_data/results/"
if(file.exists(resdir)){
  print("Directory already exists.")
}else{
  dir.create(file.path(resdir)) #create direcotry for landscape files, if it does not exist yet
} 

#Simulation options
tab <- expand_grid(landsc=landscs, spinput=spinputs, spf=spfs,
                   outputdir=resdir, alpha=alpha, beta=beta) %>%
       mutate(web=paste0("../../data/webs/",
                    str_split_fixed(basename(spinput), "_", 2)[,1], ".csv")) %>%
      dplyr::select(web, landsc, spinput, spf, outputdir, alpha, beta)

seed <- 548261 #seed for random numbers
rem <- 1500 #number of links removed per link removal iteration
nreps <- 1000 #number of iterations for Bayesian network
scenarios <- c("best", "worst", "random") #set of scenarios

#Run the simultions with link removal
for(r in 1:nrow(tab)){
    web <- tab$web[r] #name of web file (w/ path & extension)
    landsc <- tab$landsc[r] #name of landscape file (w/ path & extension)
    spinput <- tab$spinput[r] #name of species input file (w/ path & extension)
    spf <- tab$spf[r] #(deprecated)
    outputdir <-  tab$outputdir[r] #ouptut directory
    alpha <- as.numeric(tab$alpha[r]) #1st parameter of beta distr (Bayesian network)
    beta <- as.numeric(tab$beta[r]) #2nd parameter of beta distr (Bayesian network)
  
    ## 1: simulate metacommunity in the full landscape (pre patch loss)
    A <- create_adj_matrix(read_csv(web, col_types=cols())) ## create adjacency matrix from links
    landscape <- as.matrix(read_csv(landsc, col_types=cols())) ## read landscape file
    N <- length(landscape[,1]) ## number of patches in the landscape
    LINKMAT <- matrix(1, N, N) ## matrix of spatial links (1=present, 0=removed)
    outfile <- paste0(outputdir, str_replace(basename(spinput), ".csv", ""), "_",
                      str_replace(basename(landsc), ".csv", ""), "_", spf, "_",
                      alpha, "_", beta, ".rds")
    
    ## combine species and landscape information
    spdat0 <- generate_spdat(web, spinput, landsc) %>%
      add_column(delta=0, lambda=0, p=0.5, pvalues=0)
    
    spdat0 <- sim_metacomm(A, spdat0, nreps, alpha, beta, LINKMAT) %>%
      add_column(lrem=0,
                 web=web,
                 spinput=spinput,
                 landsc=landsc,
                 spf=spf,
                 seed=seed,
                 alpha=alpha,
                 beta=beta,
                 scenario="all")
    
    spdat <- spdat0
    saveRDS(spdat, outfile)
    
    flag <- FALSE
    ## if one of the following conditions is TRUE we stop here else we continue and simulate the habitat loss scenarios:
    if ((all(spdat0 %>% filter(FT!="basal") %>% pull(lambda) %>% unique < 1)) || ## stop if no consumer can persist
        (max(spdat0$p)<(1e-10)) || ## stop if persistence too low
        ((spf!="basal") & (spdat0 %>% filter(species==rownames(A)[nrow(A)]) %>% pull(lambda) %>% unique<1))){ ## stop if top species is spf and extinct
      flag <- TRUE
    }
    
    focal_sp <- ifelse(spf=="basal", rownames(A)[1], rownames(A)[nrow(A)]) ## spf name
    
    if(!flag){
      ## 2: simualte metacommunities with link loss
      for (scen in scenarios) {
        j <- 0
        rem0 <- rem
        crashed <- FALSE ## test if populations are still extant
        LINKMAT <- matrix(1, N, N) ## reset matrix of spatial links
        lselec <- numeric()
        lselec0 <- 1:choose(N, 2)
        spdat1 <- spdat0 %>% mutate(scenario=scen)
        while (!crashed) {
          j <- j + rem0
          if(scen=="best"){ ## remove least important links
            Smat <- calc_spectral(spdat1, focal_sp, LINKMAT)$smat
            Smat <- ((Smat + t(Smat))/2)*LINKMAT
            diag(Smat) <- 0
            Smax <- max(Smat[upper.tri(Smat)])
            Smat[Smat==0] <- Smax + 1
            lselec <- c(lselec, order(Smat[upper.tri(Smat)])[1:rem0])
            LINKMAT[upper.tri(LINKMAT)][lselec] <- 0
            LINKMAT <- t(LINKMAT)
            LINKMAT[upper.tri(LINKMAT)][lselec] <- 0
          }
          if(scen=="worst"){ ## remove least important links
            Smat <- calc_spectral(spdat1, focal_sp, LINKMAT)$smat
            Smat <- ((Smat + t(Smat))/2)*LINKMAT
            diag(Smat) <- 0
            lselec <- c(lselec, rev(order(Smat[upper.tri(Smat)]))[1:rem0])
            LINKMAT[upper.tri(LINKMAT)][lselec] <- 0
            LINKMAT <- t(LINKMAT)
            LINKMAT[upper.tri(LINKMAT)][lselec] <- 0
          }
          if(scen=="random"){ ## remove links randomly
            set.seed(seed)
            lselec <- sample(choose(N, 2), j)
            LINKMAT[upper.tri(LINKMAT)][lselec] <- 0
            LINKMAT <- t(LINKMAT)
            LINKMAT[upper.tri(LINKMAT)][lselec] <- 0
          }
          spdat1 <- sim_metacomm(A,spdat1,nreps,alpha,beta,LINKMAT) %>% mutate(lrem=j)
          spdat <- rbind(spdat, spdat1)
          saveRDS(spdat, outfile)
          L0 <- sum(LINKMAT[upper.tri(LINKMAT)]) ## no. of links left (in upper tri)
          if (L0<=rem0) rem0 <- 150
          write(paste(L0, scen), stdout())
          # image(t(apply(LINKMAT,2,rev)))
          if((all((spdat1 %>% filter(FT!="basal") %>% pull(lambda) %>% unique)<1)) || ## stop if all consumers have crashed
             (max(spdat1$p)<=(1e-10)) || (L0<=100)) { ## stop if persistence too low or too few links
            crashed <- TRUE
          }
        }
      }
    }
}


##Summarize simulation output: 
#Load R function
source("generate_summary_linkrem.R")

resdir <- "../../data/linkrem_data/results/"  #Directory containing result files
sumdir <- "../../data/linkrem_data/summaries/"  #Output directory for summaries

if(file.exists(sumdir)){
  print("Directory already exists.")
}else{
  dir.create(file.path(sumdir)) #create direcotry for the summaries, if it does not exist yet
} 

webs <- tab %>% pull(web) %>% basename %>% str_extract(".*[^.csv]") %>% unique 

for(web in webs){ #generate one or multiple summary files 
  resultfiles <- paste0(resdir, web,"*.rds") %>% Sys.glob() #result files
  summaryfile <- paste0(sumdir, "summary_", web, ".rds") 
  #Note: Due to their size, we generated for the model food webs a separate summary file   
  #for each functional form of a consumer's response to resource loss. 
  summary <- generate_summary(resultfiles) #summarize simulation output
  summary %>% saveRDS(., summaryfile) #save summary in RDS format
}




