rm(list=ls()) #clear R history

#Load packages
require(NetIndices, quietly=TRUE) #estimate network indices  
require(igraph, quietly=TRUE) #network analysis 
require(tidyverse, quietly=TRUE) #data handling 
require(Rcpp, quietly=TRUE) #interfacing with C++

#Load R functions
source("../create_adj_matrix.R")
source("../generate_spdat.R")
source("../generate_coordinates.R")
source("sim_metacomm_ryser.R")

#Input/output directory 
ryserdir <- "../../data/ryser_data/"
if(file.exists(ryserdir)){
  print("Directory already exists.")
}else{
  dir.create(file.path(ryserdir)) #create directory for input/output, if it does not exist yet
} 

#Generate landscapes 
no_patches <- seq(16, 70, 2) #vector with number of patches
for (n in no_patches) {
  lsc <- generate_coordinates(n) %>% as_tibble %>% rename(x=V1, y=V2)
  write_csv(lsc, paste0(ryserdir,"/landscapes/", n, ".csv"))
}

#Simulation options 
nreps <- 1000  #number of iterations for the Bayesian network
tab <- expand_grid(
  web=paste0(ryserdir,"/webs/web_1_1.csv"),
  landsc=Sys.glob(paste0(ryserdir,"/landscapes/*.csv")),
  spinputs=Sys.glob(paste0(ryserdir,"/species_input/*.csv")),
  spf="basal", #focal species for habitat loss; can either be 'basal' or 'top'
  outputdir=ryserdir, #output directory
  alpha=c(1, 5), #1st parameter of beta distribution (Bayesian network) 
  beta=c(1, 5),  #2nd parameter of beta distribution (Bayesian network)
  xiscale=10^(seq(-3, -1.2, l=31))) #vector with scaling factors for dispersal distances 

#Generate output file 
outfile <- paste0(ryserdir,"/results/tldat.csv") #name of the output file
fileconn <- file(outfile)
writeLines(paste0("spinput,xi,N,alpha,beta,tl"), fileconn) 
close(fileconn)

#Run the simulations
for (r in 1:nrow(tab)) {
  web <- tab$web[r]
  landsc <- tab$landsc[r]
  spinput <- tab$spinputs[r]
  spf <- tab$spf[r]
  outputdir <- tab$outputdir[r]
  alpha <- tab$alpha[r]
  beta <- tab$beta[r]
  xiscale <- tab$xiscale[r]
  A <- create_adj_matrix(read_csv(web, col_types=cols())) #adj matrix from links
  landscape <- as.matrix(read_csv(landsc, col_types=cols())) #landscape file
  spdat0 <- generate_spdat(web, spinput, landsc) %>% #generate species and landscape information
    add_column(delta=0, p=0.5) %>%
    mutate(xi=xi*xiscale)
  spdat <- sim_metacomm(A, spdat0, nreps, alpha, beta) %>% #simulate the metacommunity
    add_column(prem=0,
               web=web,
               spinput=spinput,
               landsc=landsc,
               spf=spf,
               seed=0,
               alpha=alpha,
               beta=beta,
               scenario="all")
  maxtl <- spdat %>% filter(p>1e-10) %>% pull(TL) %>% max #select maximum trophic level 
  write(paste(basename(spinput),xiscale,nrow(landscape),alpha,beta,maxtl,sep=","), #write output to file
             file=outfile, append=TRUE)
  cat(paste(r, "/", nrow(tab), "finished \n"))
}

