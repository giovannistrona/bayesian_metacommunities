
## Script to create edgelists using adjacency matrices from C-script as input 
require(igraph)
require(tidyverse)

## Make edge list from adjacency matrix. 
## Input: 
## - webfile: file containig the adjacency matrix
## Output: 
## - edgelist: edgelist, i.e who eats whom
make_edgelist <- function(webfile ){
  Amat <- read.table(webfile) %>% as.matrix ## make web into a matrix 
  colnames(Amat) <- rownames(Amat) <- paste0("species", sprintf("%03d", 1:nrow(Amat))) 
  g <- graph_from_adjacency_matrix(Amat, mode = "directed") ## generate graph
  if(!is_dag(g)) print("Web is cyclic.") ## test for cycles
  if(is_dag(g)){
    g2 <- graph_from_adjacency_matrix(Amat, mode = "directed") ## regeneraten directed graph 
    web <- g2 %>% get.edgelist %>% as_tibble %>% ## generate edgelist
      rename(consumer=V1, resource=V2) ## rename columns
  }
  return(web)
}

webs <- Sys.glob("../data/c_webs/*.csv") ## web files generated with C++
outputdir <- Sys.glob("../data/webs/") ## output directory for edgelists 

for(web in webs)
  web %>% make_edgelist %>% write_csv(., paste0(outputdir, web %>% basename))