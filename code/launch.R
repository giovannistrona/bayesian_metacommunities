rm(list=ls())

#Load packages
require(NetIndices) #estimate network indices  
require(igraph) #network analysis 
require(tidyverse) #data handling 

##Generate random landscapes: 
#Load R functions 
source("generate_coordinates.R")
N <- 300 #number of patches
seed <- 54321 #random seed (optional, included here for reproducibility)

landscapedir <- "../data/landscapes/"
if(file.exists(landscapedir)){
  print("Directory already exists.")
}else{
  dir.create(file.path(landscapedir)) #create direcotry for landscape files, if it does not exist yet
} 

for(k in 1:5){
  set.seed(seed*k) 
  landscape <- generate_coordinates(N) %>% as_tibble %>% rename(x=V1, y=V2) %>% #generate and
    write_csv(sprintf("../data/landscapes/landscape_%03d.csv", k)) #save landscapes
}

##Generate species input: 
spinputdir <- "../data/specie_input/"
if(file.exists(spinputdir)){
  print("Directory already exists.")
}else{
  dir.create(file.path(spinputdir)) #create direcotry for species input files, if it does not exist yet
} 
#Source R script
source("generate_species_input.R") 

##Set up the simulations: 
#Simulation input 
landscs <- Sys.glob("../data/landscapes/landscape_*.csv") #landscapes
spinputs <- Sys.glob("../data/species_input/*.csv") #species input files

#Simulation options
spfs <- "basal" #focal species for habitat loss; can either be 'basal' or 'top'
alpha <- c(1, 5) #1st parameter of beta distribution (Bayesian network) 
beta <- c(5, 1)  #2nd parameter of beta distribution (Bayesian network)
seed <- 54321 #seed for random numbers ## did I forget to set this seed in the cluster simulations!?!
rem <- 10 #number of patches removed per patch loss iteration
nreps <- 1000 #number of iterations for the Bayesian network

tab <- expand.grid(landsc=landscs, spinput=spinputs, spf=spfs, alpha=alpha, beta=beta, 
                   seed=seed, rem=rem, nreps=nreps, stringsAsFactors=FALSE) %>% as_tibble 
tab <- tab %>% add_column(web = paste0("../data/webs/",  #extract web file 
                str_split_fixed(basename(tab %>% pull(spinput)), "_", 2)[,1], ".csv"), .before=1) #from species input file

#Output directory 
outputdir <- "../data/results/" 

if(file.exists(outputdir)){
  print("Directory already exists.")
}else{
  dir.create(file.path(outputdir)) #create direcotry for simulation output, if it does not exist yet
} 

#Load R functions
source("run_instance.R")
source("generate_summary.R")

#Simulate Bayesian metacommunities with habitat loss
for(r in 1:nrow(tab)){
  output <- run_instance(web=tab$web[r], landsc=tab$landsc[r], spinput=tab$spinput[r], spf=tab$spf[r],  
                alpha=tab$alpha[r], beta=tab$beta[r], seed=tab$seed[r], rem=tab$rem[r], nreps=tab$nreps[r])
  
  outfile <- paste0(outputdir, str_replace(basename(tab$spinput[r]), ".csv", ""), "_",
                  str_replace(basename(tab$landsc[r]), ".csv", ""), "_", tab$spf[r], "_",
                  tab$alpha[r], "_", tab$beta[r], ".rds") #output file (w/ path & extension)
  
  output %>% saveRDS(., outfile) #save the simulation output in RDS format
}

##Summarize simulation output: 
#Output directory for summaries
sumdir <- "../data/summaries/" 

if(file.exists(sumdir)){
  print("Directory already exists.")
}else{
  dir.create(file.path(sumdir)) #create direcotry for the summaries, if it does not exist yet
} 

resultfiles <- Sys.glob("../data/results/*.rds") #result files
summaryfile <- "../data/summaries/summary.rds" 

summary <- generate_summary(resultfiles) #summarize simulation output
summary %>% saveRDS(., summaryfile) #save summary in RDS format

