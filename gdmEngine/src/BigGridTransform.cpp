#include <fstream>
#include <string>
#include <vector>
#include "Rcpp.h"
using namespace Rcpp;
//    Add header of the cpp code that does stuff here, e.g.:  #include "AdaptorFileUtils_V2.h"

//' @title Create GDM transformed grids for big grids
//' @description Calculate the dissimilarity between specified pairs of sites. 
//' @param site_spp (matrix) A matrix of the indices for species (col2) occurring in each site (col1), with each row an occurrence
//' @param pair_rows (matrix) A matrix of the site index for site 1 (col 1) and site 2 (col 2) in each pair of sites (rows)
//' @param site_rich (vector) A vector giving the total number of species in each site (element)
//' @param max_richness (integer)
//' @return Vector, the Sorensen dissimilarity between the specified pairs.
//' @examples output = PairsDissim(site.spp, ij.pairs, site.rich, max.rich)
//' @export
//[[Rcpp::export]]

int BigGridTransform(int nrows,
                     int ncols,
                     int nlayers,
                     Rcpp::StringVector InFilePaths,
                     Rcpp::StringVector OutFilePaths)
{
  int out;
  out = nrows + ncols; 
  return(out);
              
  // int n_sites = site_rich.size();
  // int n_pairs = pair_rows.nrow();
  // int n_records = site_spp.nrow();  
  // IntegerMatrix comp(ncols,nlayers);        
  // IntegerVector upto_index(n_sites);
  // NumericVector out(n_pairs);
  // 
  // 
  // for(int i_site=0; i_site < n_sites; i_site++) {
  //   upto_index[i_site] = 0;
  //   }
  // for(int i_rec = 0; i_rec < n_records; i_rec++) {
  //   int site_index = site_spp(i_rec,0) - 1;
  //   comp(site_index, upto_index[site_index]) = site_spp(i_rec,1);
  //   upto_index[site_index] += 1;
  //   }
  // 
  // for(int i = 0; i < n_pairs; i++) {
  //   int site_one_index = pair_rows(i,0) - 1;
  //   int site_two_index = pair_rows(i,1) - 1;
  //   float n_spp_common = 0;
  //   for(int i_spp_one=0; i_spp_one<site_rich[site_one_index]; i_spp_one++){
  //     for(int i_spp_two=0; i_spp_two<site_rich[site_two_index]; i_spp_two++){
  //       if(comp(site_one_index,i_spp_one) == comp(site_two_index,i_spp_two)){
  //         n_spp_common += 1;
  //         }
  //       }      
  //     }  
  //   float sum_rich = site_rich[site_one_index] + site_rich[site_two_index];
  //   out[i] = (1 - ((2 * n_spp_common) / (sum_rich)));
  //   }
  // return out;
  }
