//---------------------------------------------------------------------------

#ifndef esri_headersH
#define esri_headersH

//---------------------------------------------------------------------------

void SkipASCIIHeader(ifstream* argFileStream);

void LoadASCIIHeader(ifstream* argFileStream,                        long int* arg_i_ncols_ptr,						long int* arg_i_nrows_ptr,
						float* arg_f_xll_ptr,						float* arg_f_yll_ptr,						float* arg_f_cellsize_ptr,                        float* arg_f_nodata_ptr);
void OutputASCIIHeader(ofstream* argFileStream,					   const int arg_n_cols,					   const int arg_n_rows,					   const float arg_f_xll,					   const float arg_f_yll,					   const float arg_f_cellsize,					   const float arg_f_nodata);

void OutputFLTHeader(ofstream* argFileStream,
					   const int arg_n_cols,					   const int arg_n_rows,					   const float arg_f_xll,					   const float arg_f_yll,					   const float arg_f_cellsize,					   const float arg_f_nodata,					   const int arg_b_is_LSBFIRST);  // yes for Intel, no for UNIX
void OutputMCHHeader(ofstream* argFileStream,
					   const int arg_n_cols,                       const int arg_n_rows,					   const float arg_f_xll,
					   const float arg_f_yll,					   const float arg_f_cellsize,					   const float arg_f_nodata,					   const int arg_b_is_LSBFIRST,   // yes for Intel, no for UNIX					   const int arg_n_layers) ;
void OutputMCHHeader(ofstream* argFileStream,
					   const int arg_n_cols,                       const int arg_n_rows,					   const float arg_f_xll,
					   const float arg_f_yll,					   const float arg_f_cellsize,					   const float arg_f_nodata,					   const int arg_b_is_LSBFIRST,   // yes for Intel, no for UNIX					   const int arg_n_layers,
					   const int arg_n_cells) ;    // number of valid data cells
#endif
