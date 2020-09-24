
## Randomly draw patch coordinates. 
## Input:
## - N: number of patches in the landscape, 
## - dim: dimension of the landscape (1, 2 or 3). 
## Output: 
## - matrix with N rows and dim columns. 
generate_coordinates <- function(N, dim) {
    return(matrix(runif(N*dim, 0, 1), nrow=N, ncol=dim))
}

## Create adjacency matrix from a data frame of edges.
## Input:
## - edgelist: data frame with two columns of taxon names, where consumers
##             are in column 1 and resources in column 2
## Output:
## - The adjacency matrix, with A[i,j]=1 if species i eats j and 0 otherwise.
##   Species are topologically sorted to make the matrix upper triangular.
##   The rows and columns of the matrix are labeled by taxon names.
create_adj_matrix <- function(edgelist) {
    ## Create adjacency matrix A; links point from resource to consumer, so
    ## the two columns of the edge list data frame are flipped
    A <- graph_from_data_frame(d=edgelist[,2:1]) %>% ## create igraph graph
        as_adjacency_matrix %>% ## convert to igraph adjacency matrix
        as.matrix %>% ## convert to a regular matrix
        t ## transpose result (so A[i,j] is 1 if i eats j and 0 otherwise)
    ## Find a sorting of A's rows and columns to make A lower triangular
    o <- graph_from_adjacency_matrix(A) %>% topo_sort(mode="in")
    ## Return sorted, lower triangular adjacency matrix A
    return(A[o,o])
}

## Create disperal matrix M. 
## Input: - coords: patch coordinatexs x and y, 
## - kernel: type of dispersal kernel (species-specific), 
## - xi: dispersal distance (species-specific). 
## Output:
## - dispersal matrix M. 
get_matrix <- function(coords, kernel, xi) {
    if(!exists("TORUS")) TORUS <- 0
    if(TORUS == 0) {
        distances <- as.matrix(dist(coords))
    } else if (TORUS == 1) {
        require("som.nn")
        distances <- as.matrix(dist.torus(coords))
        detach(package:som.nn, unload=TRUE)
    }
    if (kernel=="Exponential") M <- exp(-distances/xi)
    if (kernel=="Gaussian")    M <- exp(-distances^2/(2*xi^2))
    if (kernel=="Rectangular") M <- 1*(distances<xi)
    diag(M) <- 0 # remove the diagonal from the matrix
    return(unname(M))
}

## Generate data file with species and landscape information. 
## Input: 
## - webfile: file containig the edgelist, i.e. who eats whom
## - datfile: species input file containing for each species sp its
##    - pi: baseline extinction probability,
##    - xi: dispersal distance,
##    - kernel: shape of the dispersal kernel ("Gaussian", "Exponential", or
##              "Rectangular").
## - landscapefile: file containing the patch coordinates.
## Output: 
## - spdat: data file combining species and landscape information. 
generate_spdat <- function(webfile, datfile, landscapefile){
    A <- read_csv(webfile, col_types=cols()) %>% ## read landscape file
        create_adj_matrix ## create adjacency matrix from links
    landscape <- read_csv(landscapefile, col_types=cols()) %>%
        as.matrix
    info <- read_csv(datfile, col_types=cols()) ## import species information
    N <- length(landscape[,1]) ## number of patches
    spdat <- tibble(patch=integer(),
                    x=numeric(),
                    y=numeric(),
                    species=character(),
                    pi=numeric(),
                    xi=numeric(),
                    kernel=character(), 
                    TL = numeric(),
                    OI = numeric(),
                    FT = character()) 
    dat <- tibble(patch=1:N, x=landscape[,1], y=landscape[,2])
    for(sp in rownames(A)){
        dat$species <- rep(sp, N)
        dat$pi <- rep(info$pi[info$species==sp], N)
        dat$xi <- rep(info$xi[info$species==sp], N)
        dat$kernel <- rep(info$kernel[info$species==sp], N)
        dat$TL <- rep(info$TL[info$species==sp], N)
        dat$OI <- rep(info$OI[info$species==sp], N)
        dat$FT <- rep(info$FT[info$species==sp], N)
        spdat <- bind_rows(spdat, dat)
    }
    return(spdat)
}

## Calculate spectral properties of principal map's Jacobian at 0, for given sp.
## Input:
## - spdat: species and landscape data
## - focal_sp: focal species for which patch values are calculated
## - L: matrix of link removal
## Output:
## - a list with elements:
##   * pvec: vector with patch values for species i.   
calc_spectral <- function(spdat, focal_sp, L){
    M <- get_matrix(coords=cbind(spdat$x[spdat$species==focal_sp], 
                                 spdat$y[spdat$species==focal_sp]), 
                    kernel=spdat$kernel[spdat$species==focal_sp][1],
                    xi=spdat$xi[spdat$species==focal_sp][1]) ## dispersal matrix
    M <- M*L
    ## Jacobian of the principal map evaluated at 0:
    Eik <- -log(1-spdat$delta[spdat$species==focal_sp])
    Amat <- M/Eik
    ## approximate relative patch values w/ leading eigenvector of the Jacobian 
    esys <- eigen(Amat) ## eigensystem of Amat
    lew <- max(abs(esys$values)) ## dominant eigenvalue
    ind <- which.max(abs(esys$values)) ## index of dominant eigenvalue
    rev <- Re(esys$vectors[,ind]) ## dominant right eigenvector
    lev <- Re(eigen(t(Amat))$vectors[,ind]) ## dominant left eigenvector
    smat <- lev%o%rev/as.numeric(lev%*%rev) ## sensitivity matrix
    return(list(lambda=lew, smat=smat)) ## return dominant ew and sens. matrix
}

## Obtain metacommunity equilibrium state.
## Input:
## - A: adjacency matrix
## - spdat: species information table
## - nreps: number of iterations for the Monte Carlo sim of the Bayesian network
## - alpha: first parameter of the Beta distribution
## - beta: second parameter of the Beta distribution
## - L: matrix of link removal
## Output:
## - a species information table (spdat) updated with the correct equilibrium
##   marginal extinction probabilities delta_i^k and patch occupancies p_i^k
sim_metacomm <- function(A, spdat, nreps, alpha, beta, L) {
    N <- length(unique(spdat$patch)) ## no of patches
    marginal <- rep(0, N) ## marginal ext prob of sp in all patches
    for (sp in rownames(A)) { ## calculate deltas of each sp from Bayesian NW:
        if (sum(A[sp,])==0) { ## basal species
            marginal <- spdat$pi[spdat$species==sp]
        } else { ## consumer species
            for (k in 1:N) { ## for each patch:
                be <- spdat$pi[spdat$species==sp][k] ## baseline ext in patch k
                deltas <- spdat %>% ## marginal ext of all species in patch k
                    filter(patch==unique(spdat$patch)[k]) %>% 
                    mutate(d=1-(1-delta)*p) %>% 
                    pull(d)
                marginal[k] <- getMarginal(A[sp,], deltas, be, alpha, beta, nreps)
            }
            marginal <- marginal/nreps
        }
        spdat$delta[spdat$species==sp] <- marginal ## assign ext prob to sp
        M <- get_matrix(coords=cbind(spdat$x[spdat$species==sp], 
                                     spdat$y[spdat$species==sp]), 
                        kernel=spdat$kernel[spdat$species==sp][1],
                        xi=spdat$xi[spdat$species==sp][1]) ## disp matrix of sp
        M <- M*L ## remove links speciied in L
        ## solve for equilibrium patch occupancies of sp numerically
        spdat$p[spdat$species==sp] <- itersol(spdat$p[spdat$species==sp], M,
                                              spdat$delta[spdat$species==sp])
        ## calculate sp's metapopultion capacity
        spectr <- calc_spectral(spdat, sp, L)
        spdat$lambda[spdat$species==sp] <- spectr$lambda
    }
    return(spdat)
}

## Calculate food web properties. 
## Input:
## - A: adjacency matrix
## Output:
## - summary statistics on species' trophic levels
calc_NetInd <- function(A) {
    g <- graph_from_adjacency_matrix(t(A), mode="directed", weighted=NULL)
    ##g <- delete.vertices(simplify(g), degree(g) == 0)
    if (vcount(g)!=0) {
        fw <- list()
        fw$S <-vcount(g)
        fw$connectance <- sum(degree(g, mode ="in"))/fw$S^2
        fw$generality <- mean(degree(g, mode ="in"))
        fw$vulnerability <- mean(degree(g, mode = "out"))
        fw$sd.generality <- sd(degree(g, mode="in"))
        fw$sd.vulnerability <- sd(degree(g, mode="out"))
        fw$sd.linkedness <- sd(degree(g, mode="in") + degree(g, mode="out"))
        g2 <- graph_from_adjacency_matrix(t(A), mode="undirected", weighted=NULL)
        ##g2 <- delete.vertices(simplify(g2), degree(g2) == 0)
        fw$community <- cluster_louvain(g2)
        fw$modularity <- modularity(g2, membership(fw$community))
        adj.web <- get.adjacency(g, sparse = FALSE)
        if(sum(adj.web != t(A)*1) != 0) "g and t(A) are not the same."
        fw$TL <- TrophInd(t(A))$TL
        fw$TL.max <- max(fw$TL)
        fw$TL.mean <- mean(fw$TL)
        fw$TL.sd <- sd(fw$TL)
        fw$OI <- TrophInd(t(A))$OI
        fw$OI.max <- max(fw$OI)
        fw$OI.mean <- mean(fw$OI)
        fw$OI.sd <- sd(fw$OI)
    }
    return(fw)
}

## Extract feeding type. 
## Input:
## - webfile: file containig the edgelist, i.e. who eats whom
## Output:
## - for each species, whether it is  basal ("basal"), primary consumer ("herb"),
##   higher trophic level than primary consumer but not top predator ("omni"),
##   or top predator ("pred")
get_feeding_type <- function(webfile){
    FT <- list() ## list for feeding types
    Amat <- read_csv(webfile) %>% ## read edge list file
        create_adj_matrix ## create adjacency matrix from links
    g <- graph_from_adjacency_matrix(Amat, mode = "directed") ## generate graph
    web <- g %>% get.edgelist %>% as_tibble %>% rename(consumer=V1, resource=V2) 
    S = nrow(Amat)
    Sb = length(which(rowSums(Amat)==0)) 
    Sc = S - Sb
    prey <- 0
    for (sp in rownames(Amat)){
        prey <- web$resource[web$consumer==sp] ## list of sp's prey
        if(length(prey)==0) FT[[sp]] = "basal"
        if(sum(is.element(prey,rownames(Amat)[1:Sb])*1)>0 &&
           sum(is.element(prey,rownames(Amat)[(Sb+1):S])*1)==0) 
            FT[[sp]] = "herb"
        if(sum(is.element(prey,rownames(Amat)[1:Sb])*1)>0 &&
           sum(is.element(prey,rownames(Amat)[(Sb+1):S])*1)>0) 
            FT[[sp]] = "omni"
        if (sum(is.element(prey,rownames(Amat)[1:Sb])*1)==0 &&
            sum(is.element(prey,rownames(Amat)[(Sb+1):S])*1)>0)
            FT[[sp]] = "pred"
    }
    return(unlist(FT))
}


##----------------------- Helper functions --------------------##

'%ni%' <- function(x,y)!('%in%'(x,y))

numextract <- function(string){ 
    temp <- gsub("[^[:digit:].]", "", string)
    as.numeric(str_extract(temp, "\\-*\\d+\\.*\\d*"))
}

