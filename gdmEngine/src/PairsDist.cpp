#include "Rcpp.h"

using namespace Rcpp;

//' @title Distance calculation for specified pairs
//' @description Calculate the distance between specified pairs. 
//' @param env_vars (matrix) A matrix of the standardised environment variables (cols) for each site (rows)
//' @param pair_rows (matrix) A matrix of the row index in 'env_vars' for site 1 (col 1) and site 2 (col 2) in each pair of sites (rows)
//' @return Vector, the Euclidean distance between the specified pairs.
//' @examples output = PairsDist(std.site.env.vars, ij.pairs)
//' @importFrom Rcpp cppFunction
//' @export
//[[Rcpp::export]]

NumericVector PairsDist(NumericMatrix env_vars, NumericMatrix pair_rows) {
  int n_vars = env_vars.ncol();
  int n_pairs = pair_rows.nrow();
  NumericVector out(n_pairs);
            
  for (int i = 0; i < n_pairs; i++) {
    double sum_sq_diff = 0;
    int RowOne = pair_rows(i,0) - 1;
    int RowTwo = pair_rows(i,1) - 1;
    for(int j = 0; j < n_vars; j++) {
      sum_sq_diff += pow((env_vars(RowOne,j) - env_vars(RowTwo,j)), 2);
      }    
    double pair_dist = sqrt(sum_sq_diff);
    out[i] = pair_dist;
    }

  return out;
  }

