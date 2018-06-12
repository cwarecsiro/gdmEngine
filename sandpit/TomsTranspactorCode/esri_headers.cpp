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
// Skip ASCII header reads an ASCII header in from the supplied open file stream
{
	// Read in input file headers
//*************************************************************************
void LoadASCIIHeader(ifstream* argFileStream,
						float* arg_f_nodata_ptr)
	// Read in input file headers
    *argFileStream>>f_buffer;
	*argFileStream>>s_buffer;
	*argFileStream>>f_buffer;
	*argFileStream>>s_buffer;
//*************************************************************************
//*************************************************************************
void OutputASCIIHeader(ofstream* argFileStream,
					   const float arg_f_cellsize,

//*************************************************************************
// Output ASCII header writes an ASCII header in to the supplied open file stream
void OutputFLTHeader(ofstream* argFileStream,
					   const float arg_f_yll,
{
}
// Output MCH header writes an ASCII header in to the supplied open file stream
void OutputMCHHeader(ofstream* argFileStream,
					   const float arg_f_yll,
{
}
//*************************************************************************
// Output MCH header writes an ASCII header in to the supplied open file stream
void OutputMCHHeader(ofstream* argFileStream,
					   const float arg_f_yll,
{
}
//*************************************************************************
// Load ASCII header reads an ASCII header in from the supplied open file stream
void LoadMCHHeader(ifstream* argFileStream,
						float* arg_f_nodata_ptr,
	// Read in input file headers
    *argFileStream>>f_buffer;
	*argFileStream>>s_buffer;
	*argFileStream>>f_buffer;
	*argFileStream>>s_buffer;