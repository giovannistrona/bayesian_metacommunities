#Load R functions
source("create_adj_matrix.R")

## Extract feeding type from edgelist or adjacency matrix. 
## Input:
## - web: edgelist, i.e. who eats whom or adjacency matrix. 
## Output:
## - for each species, whether it is  basal ("basal"), primary consumer ("herb"),
##   higher trophic level than primary consumer but not top predator ("omni"),
##   or top predator ("pred")
get_feeding_type <- function(web){
  feeding_type <- list() ## list for feeding types
  if(!is.matrix(web) && names(web) == c("consumer", "resource")){
    Amat <- web %>% create_adj_matrix ## create adjacency matrix from links
  }else{
    Amat <- web
  }
  g <- graph_from_adjacency_matrix(Amat, mode = "directed") ## generate graph
  web2 <- g %>% get.edgelist %>% as_tibble %>% rename(consumer=V1, resource=V2) ## (re)generate edgelist
  S <- nrow(Amat) ## total number of species 
  Sb <- length(which(rowSums(Amat)==0)) ## number of basal species 
  Sc <- S - Sb ## number of consumer species
  prey <- c() ## initialize prey vector
  for (sp in rownames(Amat)){
    prey <- web2$resource[web2$consumer==sp] ## list of sp's prey
    if(length(prey)==0) feeding_type[[sp]] = "basal"
    if(sum(is.element(prey,rownames(Amat)[1:Sb])*1)>0 &&
       sum(is.element(prey,rownames(Amat)[(Sb+1):S])*1)==0) 
      feeding_type[[sp]] = "herbivore"
    if(sum(is.element(prey,rownames(Amat)[1:Sb])*1)>0 &&
       sum(is.element(prey,rownames(Amat)[(Sb+1):S])*1)>0) 
      feeding_type[[sp]] = "omnivore"
    if (sum(is.element(prey,rownames(Amat)[1:Sb])*1)==0 &&
        sum(is.element(prey,rownames(Amat)[(Sb+1):S])*1)>0)
      feeding_type[[sp]] = "carnivore"
  }
  return(unlist(feeding_type))
}