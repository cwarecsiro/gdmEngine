//---------------------------------------------------------------------------

#ifndef TransformGridsH
#define TransformGridsH
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
                   std::string arg_OutFilePaths[]);
//---------------------------------------------------------------------------
#endif
