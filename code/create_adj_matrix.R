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
