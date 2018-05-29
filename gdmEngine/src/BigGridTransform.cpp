#include <fstream>
#include <string>
#include <vector>
#include "TransformGrids.h"
#include "DynamicArray.h"
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
                     int nlayers_in,
                     int geo,
                     int nSplinesTotal,
                     IntegerVector splines,
                     NumericVector knots,
                     NumericVector coefficients,
                     int extrap_code,
                     Rcpp::StringVector InFilePaths,
                     Rcpp::StringVector OutFilePaths)
{
  
  // For convenience, lets now shoot this straight to a stand-alone cpp function
  int i_spl;
  int i_prd;
  int out;
  int* Asplines;
  float* Aknots;
  float* Acoefficients;
  // hard code the maximum number of predictors for starters (n=50)
  std::string AInFilePaths[50];
  std::string AOutFilePaths[52];
  
  // Allocate arrays
  Asplines=AllocateInt1DArray(nlayers_in + geo);
  Aknots=AllocateFloat1DArray(nSplinesTotal);
  Acoefficients=AllocateFloat1DArray(nSplinesTotal);  
  
  // read the filepaths from R to the string arrays
  for(i_prd=0; i_prd<nlayers_in; i_prd++)
    {
    AInFilePaths[i_prd] = InFilePaths[i_prd];
    AOutFilePaths[i_prd] = OutFilePaths[i_prd];
    } // end for i_prd  
  if(geo>0)
    {
    AOutFilePaths[nlayers_in] = OutFilePaths[nlayers_in];
    AOutFilePaths[(nlayers_in+1)] = OutFilePaths[(nlayers_in+1)];
    }// end if geo>0
  
  // Read the data from R to the arrays
  for(i_prd=0; i_prd<(nlayers_in + geo); i_prd++)
    {
    Asplines[i_prd] = splines[i_prd];
    } // end for i_prd
  
  for(i_spl=0; i_spl<nSplinesTotal; i_spl++)
    {
    Aknots[i_spl] = knots[i_spl];
    Acoefficients[i_spl] = coefficients[i_spl]; 
    } // end for i_spl
  
  
  out = TransformGrids(nrows,
                       ncols,
                       nlayers_in,
                       geo,
                       nSplinesTotal,
                       Asplines,                         
                       Aknots,
                       Acoefficients,
                       extrap_code,
                       AInFilePaths,
                       AOutFilePaths);
  
  
  // Free the dynamic arrays
  FreeInt1DArray(Asplines);
  FreeFloat1DArray(Aknots);
  FreeFloat1DArray(Acoefficients);
  
  return(out);

  } // end BigGridTransform
