//---------------------------------------------------------------------------
#include <string>
#include <iostream>
#include <fstream>
#include <vcl.h>
#pragma hdrstop

#include "esri_headers.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

using namespace std;

//*************************************************************************
// Skip ASCII header reads an ASCII header in from the supplied open file stream// and does nothing with it..just makes code in calling func more readable//*************************************************************************void SkipASCIIHeader(ifstream* argFileStream)
{string s_buffer;float f_buffer;long int i_buffer;
	// Read in input file headers	*argFileStream>>s_buffer;	*argFileStream>>i_buffer;	*argFileStream>>s_buffer;	*argFileStream>>i_buffer;	*argFileStream>>s_buffer;	*argFileStream>>f_buffer;	*argFileStream>>s_buffer;	*argFileStream>>f_buffer;	*argFileStream>>s_buffer;	*argFileStream>>f_buffer;	*argFileStream>>s_buffer;	*argFileStream>>i_buffer;}
//*************************************************************************// Load ASCII header reads an ASCII header in from the supplied open file stream// and fills calling variables..just makes code in calling func more readable//*************************************************************************
void LoadASCIIHeader(ifstream* argFileStream,						long int* arg_i_ncols_ptr,						long int* arg_i_nrows_ptr,						float* arg_f_xll_ptr,                        float* arg_f_yll_ptr,						float* arg_f_cellsize_ptr,
						float* arg_f_nodata_ptr){string s_buffer;float f_buffer;long int i_buffer;
	// Read in input file headers	*argFileStream>>s_buffer;	*argFileStream>>i_buffer;	*arg_i_ncols_ptr=i_buffer;	*argFileStream>>s_buffer;	*argFileStream>>i_buffer;	*arg_i_nrows_ptr=i_buffer;	*argFileStream>>s_buffer;
    *argFileStream>>f_buffer;	*arg_f_xll_ptr=f_buffer;
	*argFileStream>>s_buffer;
	*argFileStream>>f_buffer;	*arg_f_yll_ptr=f_buffer;
	*argFileStream>>s_buffer;	*argFileStream>>f_buffer;	*arg_f_cellsize_ptr=f_buffer;	*argFileStream>>s_buffer;	*argFileStream>>f_buffer;	*arg_f_nodata_ptr=f_buffer;}
//*************************************************************************// Output ASCII header writes an ASCII header in to the supplied open file stream//
//*************************************************************************
void OutputASCIIHeader(ofstream* argFileStream,					   const int arg_n_cols,					   const int arg_n_rows,                       const float arg_f_xll,					   const float arg_f_yll,
					   const float arg_f_cellsize,					   const float arg_f_nodata){  *argFileStream<<"ncols         ";  *argFileStream<<arg_n_cols<<endl;  *argFileStream<<"nrows         ";  *argFileStream<<arg_n_rows<<endl;  *argFileStream<<"xllcorner     ";  *argFileStream<<arg_f_xll<<endl;  *argFileStream<<"yllcorner    ";  *argFileStream<<arg_f_yll<<endl;  *argFileStream<<"cellsize      ";  *argFileStream<<arg_f_cellsize<<endl;  *argFileStream<<"NODATA_value  ";  *argFileStream<<arg_f_nodata<<endl;}

//*************************************************************************
// Output ASCII header writes an ASCII header in to the supplied open file stream////*************************************************************************
void OutputFLTHeader(ofstream* argFileStream,					   const int arg_n_cols,                       const int arg_n_rows,					   const float arg_f_xll,
					   const float arg_f_yll,					   const float arg_f_cellsize,					   const float arg_f_nodata,                       const int arg_b_is_LSBFIRST)  // yes for Intel, no for UNIX
{  *argFileStream<<"ncols         ";  *argFileStream<<arg_n_cols<<endl;  *argFileStream<<"nrows         ";  *argFileStream<<arg_n_rows<<endl;  *argFileStream<<"xllcorner     ";  *argFileStream<<arg_f_xll<<endl;  *argFileStream<<"yllcorner    ";  *argFileStream<<arg_f_yll<<endl;  *argFileStream<<"cellsize      ";  *argFileStream<<arg_f_cellsize<<endl;  *argFileStream<<"NODATA_value  ";  *argFileStream<<arg_f_nodata<<endl; // if(arg_b_is_LSBFIRST)	*argFileStream<<"BYTEORDER  LSBFIRST"<<endl;//  else//    *argFileStream<<"BYTEORDER  MSBFIRST"<<endl;
}//*************************************************************************
// Output MCH header writes an ASCII header in to the supplied open file stream// special extra data: NUM LAYERS////*************************************************************************
void OutputMCHHeader(ofstream* argFileStream,					   const int arg_n_cols,                       const int arg_n_rows,					   const float arg_f_xll,
					   const float arg_f_yll,					   const float arg_f_cellsize,					   const float arg_f_nodata,					   const int arg_b_is_LSBFIRST,   // yes for Intel, no for UNIX					   const int arg_n_layers)
{//Standard ESRI hdr  *argFileStream<<"ncols         ";  *argFileStream<<arg_n_cols<<endl;  *argFileStream<<"nrows         ";  *argFileStream<<arg_n_rows<<endl;  *argFileStream<<"xllcorner     ";  *argFileStream<<arg_f_xll<<endl;  *argFileStream<<"yllcorner    ";  *argFileStream<<arg_f_yll<<endl;  *argFileStream<<"cellsize      ";  *argFileStream<<arg_f_cellsize<<endl;  *argFileStream<<"NODATA_value  ";  *argFileStream<<arg_f_nodata<<endl;  *argFileStream<<"BYTEORDER  LSBFIRST"<<endl;    // keep byteordeer in so the file can be used as a hdr// Extra line  *argFileStream<<"nlayers       ";  *argFileStream<<arg_n_layers<<endl;
}
//*************************************************************************
// Output MCH header writes an ASCII header in to the supplied open file stream// special extra data: NUM LAYERS////*************************************************************************
void OutputMCHHeader(ofstream* argFileStream,					   const int arg_n_cols,                       const int arg_n_rows,					   const float arg_f_xll,
					   const float arg_f_yll,					   const float arg_f_cellsize,					   const float arg_f_nodata,					   const int arg_b_is_LSBFIRST,   // yes for Intel, no for UNIX					   const int arg_n_layers,					   const int arg_n_cells)
{//Standard ESRI hdr  *argFileStream<<"ncols         ";  *argFileStream<<arg_n_cols<<endl;  *argFileStream<<"nrows         ";  *argFileStream<<arg_n_rows<<endl;  *argFileStream<<"xllcorner     ";  *argFileStream<<arg_f_xll<<endl;  *argFileStream<<"yllcorner    ";  *argFileStream<<arg_f_yll<<endl;  *argFileStream<<"cellsize      ";  *argFileStream<<arg_f_cellsize<<endl;  *argFileStream<<"NODATA_value  ";  *argFileStream<<arg_f_nodata<<endl;  *argFileStream<<"BYTEORDER  LSBFIRST"<<endl;    // keep byteordeer in so the file can be used as a hdr// Extra line  *argFileStream<<"nlayers       ";  *argFileStream<<arg_n_layers<<endl;  *argFileStream<<"ndatacells       ";  *argFileStream<<arg_n_cells<<endl;
}
//*************************************************************************
// Load ASCII header reads an ASCII header in from the supplied open file stream// and fills calling variables..just makes code in calling func more readable//*************************************************************************
void LoadMCHHeader(ifstream* argFileStream,						long int* arg_i_ncols_ptr,						long int* arg_i_nrows_ptr,						float* arg_f_xll_ptr,                        float* arg_f_yll_ptr,						float* arg_f_cellsize_ptr,
						float* arg_f_nodata_ptr,						long int* arg_i_nlayers_ptr){string s_buffer;float f_buffer;long int i_buffer;
	// Read in input file headers	*argFileStream>>s_buffer;	*argFileStream>>i_buffer;	*arg_i_ncols_ptr=i_buffer;	*argFileStream>>s_buffer;	*argFileStream>>i_buffer;	*arg_i_nrows_ptr=i_buffer;	*argFileStream>>s_buffer;
    *argFileStream>>f_buffer;	*arg_f_xll_ptr=f_buffer;
	*argFileStream>>s_buffer;
	*argFileStream>>f_buffer;	*arg_f_yll_ptr=f_buffer;
	*argFileStream>>s_buffer;	*argFileStream>>f_buffer;	*arg_f_cellsize_ptr=f_buffer;	*argFileStream>>s_buffer;	*argFileStream>>f_buffer;	*arg_f_nodata_ptr=f_buffer;	//BYTE ORDER  2 STRINGS IGNORE	*argFileStream>>s_buffer;	*argFileStream>>s_buffer;	*argFileStream>>s_buffer; // nlayers	*arg_i_nlayers_ptr=i_buffer;}
