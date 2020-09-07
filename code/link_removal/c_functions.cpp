#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;


// calculate local extinction probs (deltas) of a focal species in a focal patch
// Input:
// - Arow: one row of the food web's adjacency matrix, corresponding to the
//         focal species
// - deltas: vector of species' extinction probs in focal patch
// - be: baseline extinction probability of focal species in focal patch
// - alpha: first parameter of the Beta distribution
// - beta: second parameter of the Beta distribution
// - reps: number of iterations for the Monte Carlo sim of the Bayesian network
// Output:
// - no of times the focal species in the focal patch went extinct (out of reps)
// [[Rcpp::export]]
int getMarginal(NumericVector Arow, NumericVector deltas, double be,
                double alpha, double beta, int reps) {
  int S=Arow.size(); // number of species
  int rep, i, marginal=0;
  double randnum, frac, pext, prodsum, tot;
  NumericVector extinct(S), rands_e(S), rands_r=Rcpp::runif(reps);
  for (rep=0; rep<reps; rep++) {
    prodsum=0.0; // cumulative sum of product of Arow[i] and extinct[i]
    tot=0.0; // sum of Arow's entries (row sum of A[sp,])
    rands_e=Rcpp::runif(S); // generate S random numbers, for future use
    for (i=0; i<S; i++) { // for each species:
      extinct[i]=(rands_e[i]<deltas[i] ? 1.0 : 0.0); // extinct or extant?
      prodsum+=Arow[i]*extinct[i]; // is ith species extinct?
      tot+=Arow[i]; // accumulate sum of entries of Arow
    }
    frac=prodsum/tot; // fraction of extinct prey
    pext=be+(1-be)*R::pbeta(frac, alpha, beta, true, false); // Bayesian NW eq
    marginal+=(rands_r[rep]<pext ? 1 : 0); // extinct if random number < pext
  }
  return(marginal);
}


// obtain equilibrium patch occupancies numerically
// Input:
// - p0: vector of initial patch occupancies of focal species
// - M: dispersal matrix of focal species
// - d: vector of extinction probabilities of focal species
// Output:
// - vector of updated patch occupancies for focal species in each patch
// [[Rcpp::export]]
NumericVector itersol(NumericVector p0, NumericMatrix M, NumericVector d) {
  int N=p0.size(); // number of patches
  int converge=0; // variable to check whether numerical iteration has converged
  int iter, k, l;
  NumericVector p_next=p0, result=p0, C(N), E(N), itermap(N);
  for (k=0; k<N; k++) E[k]=log(1.0/(1.0-d[k])); // extinction rates
  while (!converge) { // if solution hasn't converged yet:
    for (iter=0; iter<5; iter++) { // for ten iterations:
      for (k=0; k<N; k++) { // calculate colonization rate:
        C[k]=0.0; // initialize ith value of colonization
        for (l=0; l<N; l++) C[k]=C[k]+M(k,l)*p_next[l]; // colonization (M%*%p)
        p_next[k]=C[k]/(C[k]+E[k]); // iteration map
      }
    }
    if (sum(pow(result-p_next, 2))<1.0e-6) converge=1; // has result converged?
    result=p_next; // in either case, update result with current estimate of p
  }
  return(result);
}
