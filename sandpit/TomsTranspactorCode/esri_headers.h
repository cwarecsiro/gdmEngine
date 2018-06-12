//---------------------------------------------------------------------------

#ifndef esri_headersH
#define esri_headersH

//---------------------------------------------------------------------------

void SkipASCIIHeader(ifstream* argFileStream);

void LoadASCIIHeader(ifstream* argFileStream,
						float* arg_f_xll_ptr,
void OutputASCIIHeader(ofstream* argFileStream,

void OutputFLTHeader(ofstream* argFileStream,
					   const int arg_n_cols,
void OutputMCHHeader(ofstream* argFileStream,
					   const int arg_n_cols,
					   const float arg_f_yll,
void OutputMCHHeader(ofstream* argFileStream,
					   const int arg_n_cols,
					   const float arg_f_yll,
					   const int arg_n_cells) ;    // number of valid data cells
#endif