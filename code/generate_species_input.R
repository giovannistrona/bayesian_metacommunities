
############################################################
# Set up to generate species input
############################################################
#Load packages
require(NetIndices) #estimate network indices  
require(igraph) #network analysis 
require(tidyverse) #data handling 

#Load R functions
source("create_adj_matrix.R")
source("get_feeding_type.R")

#Required inputs
clargs <- commandArgs(trailingOnly=TRUE)
if (length(clargs) > 0) { #command-line arguments
  outputdir <- clargs[1] #output directory 
  webfiles <- clargs[2:length(clargs)] #name of web (w/ path & extension)
} else { #sample input parameters, if no command line arguments are given
  outputdir <-  "../data/species_input/"
  webfiles <- Sys.glob("../data/webs/*.csv")
}

#Parameter options
kshape <- "Exponential" #shape of dispersal kernel 

############################################################  
#Generate species input
############################################################
for (webfile in webfiles) {

  ## create species_input file 
  A <- read_csv(webfile) %>% create_adj_matrix ## create adjacency matrix from links
  w <- sub('\\.csv$', '', basename(webfile)) 
  n <- 1 ## offset for file ending
  
  ## From Eklöf et al. MEE 2013: baseline extinction probability 'pi' increases with trophic level 'TL'
  ## pi = 0.2 * TL/mean(TL)
  pi0 <- 0.2 
  xivec <- seq(0.01, 0.1, l=5) #dispesral distances
  
  ## dispersal distance 'xi' same for all species 
  for(xi0 in xivec){
    dat <- tibble(species = rownames(A)) %>% 
      add_column(pi = 0, xi = 0, kernel = kshape, TL = TrophInd(t(A))$TL, 
                 OI = TrophInd(t(A))$OI, FT = get_feeding_type(A))
    
    for(sp in 1:length(dat$species)){
      dat$pi[sp] <- pi0 * (dat$TL[sp]/mean(dat$TL))
      dat$xi[sp] <- xi0
    }
    dat$kernel <- kshape
    
    outfile <- paste0(outputdir, w, sprintf("_%03d.csv", n))
    dat %>% write_csv(outfile)     
    n <- n + 1
  }
  
  ## dispesral distance 'xi' increases with trophic level 
  dat <- tibble(species = rownames(A)) %>% 
    add_column(pi = 0, xi = 0, kernel = kshape, TL = TrophInd(t(A))$TL, 
               OI = TrophInd(t(A))$OI, FT =get_feeding_type(A))
  
  for(sp in 1:length(dat$species)){
    dat$pi[sp] <- pi0 * (dat$TL[sp]/mean(dat$TL))
    dat$xi[sp] <- mean(xivec) * (dat$TL[sp]/mean(dat$TL))
  }
  dat$kernel <- kshape
  
  outfile <- paste0(outputdir, w, sprintf("_%03d.csv", n))
  dat %>% write_csv(outfile)     
  n <- n + 1
  
  ## baseline extinction probability 'pi' same for all species (average over TL) 
  dat <- tibble(species = rownames(A)) %>% 
    add_column(pi = 0, xi = 0, kernel = kshape, TL = TrophInd(t(A))$TL, 
               OI = TrophInd(t(A))$OI, FT = get_feeding_type(A))
  
  for(sp in 1:length(dat$species))
    dat$xi[sp] <- mean(xivec) * (dat$TL[sp]/mean(dat$TL))
  
  dat$pi <- pi0
  dat$kernel <- kshape
  
  outfile <- paste0(outputdir, w, sprintf("_%03d.csv", n))
  dat %>% write_csv(outfile)     
  n <- n + 1
  
  ## dispersal distance 'xi' and baseline extinction probability 'pi' same for all species (average over TL) 
  dat <- tibble(species = rownames(A)) %>% 
    add_column(pi = 0, xi = 0, kernel = kshape, TL = TrophInd(t(A))$TL, 
               OI = TrophInd(t(A))$OI, FT = get_feeding_type(A))
  
  dat$pi <- pi0 
  dat$xi <- mean(xivec)
  dat$kernel <- kshape
  
  outfile <- paste0(outputdir, w, sprintf("_%03d.csv", n))
  dat %>% write_csv(outfile)     
  n <- n + 1
}

#Generate additional group based species input file for the serengeti web
serengeti_files <- Sys.glob(paste0(outputdir, "serengeti_*.csv")) 
groups <- read_csv("../data/serengeti/serengeti_groups.csv")

for(serengeti_file in serengeti_files)
  serengeti_file %>% read_csv %>% left_join(groups) %>% #add group as additional column
    write_csv(serengeti_file)

dat <- read_csv(paste0(outputdir, "serengeti_006.csv")) ## 'pi' and 'xi' TLB
for(sp in 1:length(dat$species)){
  dat$pi[sp] <- 0.2*(15-dat$group[sp])/mean(dat$group) #group based 'baseline extinction probability 'pi'
  dat$xi[sp] <- 0.055*(15-dat$group[sp])/mean(dat$group) #group based 'diserpersal distance 'xi'
}
dat %>% write_csv(paste0(outputdir, sprintf("serengeti_%03d.csv", n)))

