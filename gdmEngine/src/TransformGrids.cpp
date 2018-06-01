//---------------------------------------------------------------------------
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
//#include <FileCtrl.hpp>

#include "defines.h"
//#include "esri_headers.h"
//#include "Splash.h"
#include "CompactorFunctions.h"
//#include "BulkCompactor.h"
#include "Transpactor.h"
#include "TransformGrids.h"
//---------------------------------------------------------------------------
// This is the main cpp function that transforms the set of grids
//---------------------------------------------------------------------------
int TransformGrids(int arg_nrows,
                    int arg_ncols,
                    int arg_nlayers_in,
                    int arg_geo,
                    int arg_nSplinesTotal,
                    int* arg_splines,                         
                    float* arg_knots,
                    float* arg_coefficients,
                    int arg_extrap_code,
                    std::string arg_InFilePaths[],
                    std::string arg_OutFilePaths[])
{
///////////////////////////////////////////////////////////////////////////////  
  // Specify local variables
  int i_layer;
  int out=1;
  TTranspactor* Transpactor;
  
  
  // Fire up a Transpactor (how else would you do this?)
  Transpactor=new TTranspactor();
  Transpactor->Refresh();
  Transpactor->SetNLayers(arg_nlayers_in);
  Transpactor->SetNRows(arg_nrows);
  Transpactor->SetNCols(arg_ncols);
  Transpactor->SetExtrapolation(arg_extrap_code); 
  Transpactor->SetGeo(arg_geo);
  
  // Load the nsplines, knots & coefficients
  Transpactor->SetModelParameters(arg_geo,
                                  arg_nlayers_in,
                                  arg_splines,                         
                                  arg_knots,
                                  arg_coefficients);
  
  // Load the grid paths // BUG HERE
  Transpactor->SetGridPaths(arg_geo,
                            arg_nlayers_in,
                            arg_InFilePaths,
                            arg_OutFilePaths);   
    
  // Set up the extrapolation parameters
  for(i_layer=0; i_layer<arg_nlayers_in; i_layer++)
    {
    Transpactor->SetUpExtrapolation(i_layer);
    }

  // Now tranform the grids...
  out = Transpactor->TransBinaryGridsToFloat();
  out = 1;
  
  delete Transpactor;
  
  return(out);
  
} // end TransformGrids()
//---------------------------------------------------------------------------

