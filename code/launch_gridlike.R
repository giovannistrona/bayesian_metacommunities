rm(list=ls())

#Load packages
require(NetIndices, quietly=TRUE) #estimate network indices  
require(igraph, quietly=TRUE) #network analysis 
require(tidyverse, quietly=TRUE) #data handling 
require(Rcpp, quietly=TRUE) #interfacing with C++

##If needed, generate random landscapes: 
#Load R functions 
source("generate_coordinates.R")
N <- 300 #number of patches
seed <- 54321 #random seed (optional, included here for reproducibility)

##Generate or import food webs: 
#Generate or import model/empirical webs. Web files must be in csv format, 
#and their content must be in edgelist type with two columns: 'consumer' (column 1); 'resource' (column 2). 
#We generated webs (adjacency matrices) in C++ based on the allometric trophic network model in 
#Schneider et al. 2016 (c_code), and then we transformed the adj. matrices into edgelists. 
# source("make_edgelist.R") #transform adj. matrices into edgelists

##If needed, generate species input: 
spinputdir <- "../data/species_input/"
if(file.exists(spinputdir)){
  print("Directory already exists.")
}else{
  dir.create(file.path(spinputdir)) #create direcotry for species input files, if it does not exist yet
} 

#Source R script
# source("generate_species_input.R") 

##Set up the simulations: 
#Simulation input 
landscs <- Sys.glob("../data/landscapes/gridliek_*.csv") #landscapes
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
                  tab$alpha[r], "_", tab$beta[r], "_gridlike.rds") #output file (w/ path & extension)
  
  output %>% saveRDS(., outfile) #save the simulation output in RDS format
}

##Summarize simulation output: 
#Load R function
source("generate_summary.R")

sumdir <- "../data/summaries/"  #Output directory for summaries

if(file.exists(sumdir)){
  print("Directory already exists.")
}else{
  dir.create(file.path(sumdir)) #create direcotry for the summaries, if it does not exist yet
} 

webs <- "../data/webs/*" %>% Sys.glob %>% basename %>% str_split(.,".csv", simplify=TRUE) %>%
  as_tibble %>% pull(V1)

for(web in webs){ #generate one or multiple summary files 
  resultfiles <- paste0("../data/results/", web,"*_gridlike.rds") %>% Sys.glob() #result files
  summaryfile <- paste0("../data/summaries/summary_gridlike_", web, ".rds") 
  #Note: Due to their size, we generated for the model food webs a separate summary file   
  #for each functional form of a consumer's response to resource loss. 
  summary <- generate_summary(resultfiles) #summarize simulation output
  summary %>% saveRDS(., summaryfile) #save summary in RDS format
}



