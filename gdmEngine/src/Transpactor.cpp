//---------------------------------------------------------------------------
/*****************************************************************************
 * TRANSPACTOR.H
 * Object to manage compacting of grids, and transformation
 *
 * Written October 2017 by Tom Harwood.
 *
 ***************************************************************************/
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>

#include "defines.h"  // for PI
#include "esri_headers.h"
#include "DynamicArray.h"    // functions to ease allocation and destruction of dynamic arrays
#include "CompactorFunctions.h"
#include "Transpactor.h"

//---------------------------------------------------------------------------

//****************************************************************************
// Constructor 1 initialises
//****************************************************************************
  TTranspactor::TTranspactor(void)
	{
	int i_layer; //counter
	i_extrapolation=1; //0 don't extrapolate 1 do capped minimum extrapolation 2 10% extrapolation 3 linear extrapolation
       for(i_layer=0;i_layer<50;i_layer++)
	 {
	 b_active[i_layer]=0; // set all to inactive
//	 b_isbare[i_layer]=0; // none to bare
//	 b_iszero[i_layer]=0; // set none to zero (do this when loading coeffs
	 }
	}
//****************************************************************************
// Constructor 2  initialises, using the supplied extrapolation value
//****************************************************************************
TTranspactor::TTranspactor(const int arg_i_extrapolation)
{
int i_layer; //counter


	i_extrapolation=arg_i_extrapolation; //0 don't extrapolate 1 do capped minimum extrapolation

   for(i_layer=0;i_layer<50;i_layer++)
	 {
	 b_active[i_layer]=0; // set all to inactive
	// b_isbare[i_layer]=0; // none to bare
	// b_iszero[i_layer]=0; // set none to zero (do this when loading coeffs
	 }

}
/******************************************************************************
 * TTranspactor  DESTRUCTOR
 * Doesn't do anything
 *********************************************************************/
TTranspactor::~TTranspactor( void)
{

}
void TTranspactor::Refresh(void)
{
int i_layer; //counter
   // These may change between biomes
   for(i_layer=0;i_layer<50;i_layer++)
	 {
	 b_active[i_layer]=0; // set all to inactive
//	 b_isbare[i_layer]=0; // none to bare
//	 b_iszero[i_layer]=0; // set none to zero (do this when loading coeffs
	 }

}
//************************************************************
// TransCompact
// takes a lot of binary grid files, transforms them and squishes them a lot
// floats are converted to unsigned short ints
// the format is RUN LENGTH : datapoint 1 layer 1, datapoint1 layer 2 ..datapoint1 layer n:
// datapoint 2 layer 1, datapoint2 layer 2 ..datapoint2 layer n:  to datapoint n layer n
// NODATA runs are given as NEGATIVES, with no data afterwards.
// e.g. -5 5 34 23 24 25 21 23 25 12 25 25 -40 for a 2 layer grid.. is
// 5 nodata cells, 5 data cells, first cell, 34, 23 etc, followed by 40 nodatas
//****************************************************************
/*
int TTranspactor::TransCompactBinaryGrids(string arg_s_res_path)
{
int i_result_code;
long int ntrows,ntcols; //double checkers
float f_xllcorner,f_yllcorner, f_nodata,f_cellsize;
long int i_row, i_col;  // row and column counters
int i_layer;       //layer/file counter
float** f_data;
int i_pos; // position we have reached
ifstream* InFileStream_ptr;
ofstream OutFileStream;
//ofstream IndexFileStream;
ifstream HeaderInStream;
ofstream HeaderOutStream;
float* f_ptr;    // for reading data in
long int i_start_run, i_end_run;
short int si_run_length;
short* i_s_ptr;
short i_s;
float f_value; //general use and untransformed value
float f_transval; // GDM tranformed value
string s_path,s_index_path;
string s_name;
int i_short_size;
double d_spos; // total offset of the start of each new line.
int i_length;
int* i_ptr;
int i_float_size;
int i_data_format;
double* d_ptr;
long int* l_ptr;
int i_cell_count=0;

  i_result_code=1; // good unless we hear otherwise!
  i_short_size= sizeof(short);
  i_float_size=sizeof(float);

  // Allocate a set of arrays to hold the lines for all layers
  f_data= AllocateFloat2DArray(n_layers,
							   n_cols+1);

  // Open up all grid files to compact
  InFileStream_ptr = new ifstream[n_layers];
  for(i_layer=0;i_layer<n_layers;i_layer++)
	{
	InFileStream_ptr[i_layer].open(s_binary_in[i_layer].c_str(), ios::in |ios::binary);
	// If still missing, tell the user to get their act together!
	if(!InFileStream_ptr[i_layer].is_open())
	  return 0;
	}
  //Load a header file  first file  (BARE)
  s_path=s_binary_in[0];
  i_length=s_path.length();
  s_path.resize(i_length-3);     // take off "flt"
  s_path=s_path+"hdr";
  HeaderInStream.open(s_path.c_str());
  if(!HeaderInStream.is_open())
	return 0;
  LoadASCIIHeader(&HeaderInStream,
				  &ntcols,
				  &ntrows,
				  &f_xllcorner,
				  &f_yllcorner,
				  &f_cellsize,
				  &f_nodata);
  HeaderInStream.close();

  // Open big results file AS BINARY
  OutFileStream.open(arg_s_res_path.c_str(), ios::out |ios::binary);

  //Read whole lines into buffers
  for(i_row=0;i_row<n_rows;i_row++)
	{
    // read in whole line and process
    for(i_col=0;i_col<n_cols;i_col++)
	  {
      // Read in values for each cell
	  for(i_layer=0;i_layer<n_layers;i_layer++)
		{
		f_ptr=&(f_data[i_layer][i_col]);
        //Read to f_data array
		InFileStream_ptr[i_layer].read((char*)f_ptr,sizeof(float));
		}// end for i_layer
	  } // end for i_col

	// Assume all layers are aligned, so we need only
	// do fiddling in the topmost layer
	i_pos=0;
	i_col=0;
    //Work through the whole row
    while(i_pos<n_cols)
	  {
      // split into nodata or data runs
	  if(f_data[0][i_pos]<K_NO_DATA_ISH)
		{
        //NO_DATA
      //  i_col=i_pos;
		i_start_run=i_pos;
        si_run_length=0;
		while(f_data[0][i_pos]<K_NO_DATA_ISH&&i_pos<n_cols)
		  {
          si_run_length++;  // 1 first time through
          i_pos++;
		  if(si_run_length>=K_SHORT_MAX) // don't run off the end of a short for really big grids
            break;

		  }// end while still nodata
        // we have found the end of the run.
        // Since this is no_data.. convert to negative
		si_run_length=0-si_run_length;
        // Write out this only
		i_s_ptr=&(si_run_length);
		OutFileStream.write((char*)i_s_ptr,i_short_size);   // negative run length
	  //  d_spos+=i_short_size;
        } // end NODATA run
	  else
        {
		//DATA
	   // i_col=i_pos;
		i_start_run=i_pos;
        si_run_length=0;
		while(f_data[0][i_pos]>=K_NO_DATA_ISH&&i_pos<n_cols)
          {
		  si_run_length++;
		  i_pos++;
		  if(si_run_length>=K_SHORT_MAX)
			break;
          }// end while still data
		// we have found the end of the run.   record total amount of data
		i_cell_count+=si_run_length;
        //Go back to start of run and read out ALL LAYERS
        // Write out run length first
        i_s_ptr=&(si_run_length);
        OutFileStream.write((char*)i_s_ptr,i_short_size);   // positive run length
	 //	d_spos+=i_short_size;
		// Read out whole run, all layer data for each point
		// TRANSFORM GRIDS
		// Output data for each Grid cell for all layers
        // in order of the layers.. e.g. for 3 layers 1, 2,3, with run length 2
        // 2(runlength) 1 2 3 (firstcell) 1 2 3 (second cell)
		i_end_run=i_start_run+(int)si_run_length;
		for(i_col=i_start_run;i_col<i_end_run;i_col++)
		  {
          for(i_layer=0;i_layer<n_layers;i_layer++)
			{
			//Extract untransformed layer value
			f_value=f_data[i_layer][i_col];
			//Transform this value using the appropriate method
			if(i_layer==0) //BARE
			  f_transval=GetBareValue(f_value);
			else
			  {
			  if(IsZero(i_layer))
				f_transval=0;
			  else
				f_transval= GetExtraSplineValue(f_value,
											i_layer);
			  }
			f_ptr=&(f_transval);
            OutFileStream.write((char*)f_ptr,i_float_size);   // unshortened value for this layer
		  //	d_spos+=i_float_size;
            }// end for i_layer
          }// end for i_col start to end run
        } // end DATA run
      }
	// At the end of the row, output the offset from the start of the file
//    d_ptr=&d_spos;
//    IndexFileStream.write((char*)d_ptr,sizeof(double));

    }// end for i_row
  // Output an end of file token
  i_s=K_EOF;
  i_s_ptr=&i_s;
  OutFileStream.write((char*)i_s_ptr,i_short_size);   // eoftoken
  //Close files
  OutFileStream.close();
//  IndexFileStream.close();

  for(i_layer=0;i_layer<n_layers;i_layer++)
	 InFileStream_ptr[i_layer].close();

  // build header file name)
  s_path=arg_s_res_path;
  i_length=s_path.length();
  s_path.resize(i_length-4);
  s_path=s_path+".mdr";
  HeaderOutStream.open(s_path.c_str());
  OutputMCHHeader(&HeaderOutStream,
				  n_cols,
				  n_rows,
				  f_xllcorner,
				  f_yllcorner,
				  f_cellsize,
				  f_nodata,
				  1,   // yes for Intel, no for UNIX
				  n_layers,
				  i_cell_count);
  HeaderOutStream.close();
  // Free up arrays
  FreeFloat2DArray(f_data);
  delete [] InFileStream_ptr;

  return i_result_code;
}// end ConvertGridStackToCompactFLOATBinaryFile
*/


//************************************************************
// TransBinarygrods to Float.
// Exactly the same as transcompact binary grids, but exports to a list of float grids
// takes a lot of binary grid files, transforms them and squishes them a lot
// floats are converted to unsigned short ints
// the format is RUN LENGTH : datapoint 1 layer 1, datapoint1 layer 2 ..datapoint1 layer n:
// datapoint 2 layer 1, datapoint2 layer 2 ..datapoint2 layer n:  to datapoint n layer n
// NODATA runs are given as NEGATIVES, with no data afterwards.
// e.g. -5 5 34 23 24 25 21 23 25 12 25 25 -40 for a 2 layer grid.. is
// 5 nodata cells, 5 data cells, first cell, 34, 23 etc, followed by 40 nodatas
//****************************************************************
int TTranspactor::TransBinaryGridsToFloat(void) // (string arg_s_res_path)       // note names of transgrids are in s_binary_out
  {
  int i_result_code;
  int i_layer;       //layer/file counter
  int i_pos; // position we have reached
  int i_float_size;
  int i_data_format;
  int i_length;
  int i_short_size;
  long int ntrows,ntcols; //double checkers
  long int i_row, i_col;  // row and column counters
  long int i_start_run, i_end_run;
  short int si_run_length;
  short i_s;
  float f_xllcorner,f_yllcorner, f_nodata,f_cellsize;
  float f_value; //general use and untransformed value
  float f_transval; // GDM tranformed value
  float f_xdist, f_ydist;
  float f_xtransval, f_ytransval;
  double d_spos; // total offset of the start of each new line.
  
  // pointers for things
  int* i_ptr;
  long int* l_ptr;
  short* i_s_ptr;
  float* f_ptr;    // for reading data in
  float* f_ptr2;    // for reading dataout
  float* f_ptr3;
  float* f_ptr4;
  double* d_ptr;

  // filestreams
  ifstream* InFileStream_ptr;
  ifstream HeaderInStream;
  ofstream HeaderOutStream;
  ofstream* FLTOutStream_ptr;
  ofstream* FLTGeoOutStream_ptr;  
  //ofstream OutFileStream;
  //ofstream IndexFileStream;
  
  // Strings
  string s_path,s_index_path;
  string s_name;
  
  // establish values
  i_result_code=1; // good unless we hear otherwise!
  i_short_size= sizeof(short);
  i_float_size=sizeof(float);

  // Open up all grid files to transform
  InFileStream_ptr = new ifstream[n_layers];
  FLTOutStream_ptr = new ofstream[n_layers];
  for(i_layer=0;i_layer<n_layers;i_layer++)
	  {
	  InFileStream_ptr[i_layer].open(s_binary_in[i_layer].c_str(), ios::in |ios::binary);
	  // If still missing, tell the user to get their act together!
	  if(!InFileStream_ptr[i_layer].is_open())
	    return 0;
	  }

  //Load a header file - first file  
  s_path=s_binary_in[0];
  i_length=s_path.length();
  s_path.resize(i_length-3);     // take off "flt"
  s_path=s_path+"hdr";
  HeaderInStream.open(s_path.c_str());
  if(!HeaderInStream.is_open())
	  return 0;
  LoadASCIIHeader(&HeaderInStream,
                  &ntcols,
                  &ntrows,
                  &f_xllcorner,
                  &f_yllcorner,
                  &f_cellsize,
                  &f_nodata);
  HeaderInStream.close();

  // Create all the header files
  for(i_layer=0;i_layer<n_layers;i_layer++)
    {
    // build index file name ( to hold seekposes for the beginning of each line)
    s_path=s_binary_out[i_layer];  // with .flt at the end
    i_length=s_path.length();
    s_path.resize(i_length-3);     // take off "flt"
    s_path=s_path+"hdr";
    HeaderOutStream.open(s_path.c_str());
    OutputFLTHeader(&HeaderOutStream,
                    n_cols,
                    n_rows,
                    f_xllcorner,
                    f_yllcorner,
                    f_cellsize,
                    f_nodata,
                    1);   // yes for Intel, no for UNIX
    HeaderOutStream.close();
    } // end for output header loop

//~~//~~//~~//~~//~~//  
  // If we're generating trans grids for geo, open files to write to
  if(i_geo>0)
    {
    // Create all the header files
    for(i_layer=0;i_layer<2;i_layer++)
      {
      // build index file name ( to hold seekposes for the beginning of each line)
      s_path=s_geo_binary_out[i_layer];  // with .flt at the end
      i_length=s_path.length();
      s_path.resize(i_length-3);     // take off "flt"
      s_path=s_path+"hdr";
      HeaderOutStream.open(s_path.c_str());
      OutputFLTHeader(&HeaderOutStream,
                      n_cols,
                      n_rows,
                      f_xllcorner,
                      f_yllcorner,
                      f_cellsize,
                      f_nodata,
                      1);   // yes for Intel, no for UNIX
      HeaderOutStream.close();
      } // end for output header loop
    } // end if geo>0 
  //~~//~~//~~//~~//~~//
      
  // Open big results file AS BINARY
//  OutFileStream.open(arg_s_res_path.c_str(), ios::out |ios::binary);
  // build index file name ( to hold seekposes for the beginning of each line)
  //s_path=arg_s_res_path;
  //i_length=s_path.length();
  //s_path.resize(i_length-4);

  //~~//~~//~~//~~//~~//  
  // If we're generating trans grids for geo, open files to write to
  if(i_geo>0)
    {
    FLTGeoOutStream_ptr = new ofstream[2];
    // Open all the output files
    for(i_layer=0;i_layer<2;i_layer++)
      {
      // build index file name ( to hold seekposes for the beginning of each line)
      s_path=s_geo_binary_out[i_layer];  // with .flt at the end
      FLTGeoOutStream_ptr[i_layer].open(s_path.c_str(), ios::out |ios::binary);
      } // end for output fileloop
    } // end if geo   
  //~~//~~//~~//~~//~~//    
  
	// Open all the output files
  for(i_layer=0;i_layer<n_layers;i_layer++)
    {
    // build index file name ( to hold seekposes for the beginning of each line)
    s_path=s_binary_out[i_layer];  // with .flt at the end
    FLTOutStream_ptr[i_layer].open(s_path.c_str(), ios::out |ios::binary);
    } // end for output fileloop
  d_spos=0; // zero seekpos counter
  //Read in data
  f_ydist = 0;
  for(i_row=0;i_row<n_rows;i_row++)
	  {
    // read in whole line and process
    f_xdist = 0;
    for(i_col=0;i_col<n_cols;i_col++)
	    {
      // Read in values for each cell
	    for(i_layer=0;i_layer<n_layers;i_layer++)
		    {
		    f_ptr=&f_value;
		    f_ptr2=&(f_transval);
		    //Read data
		    InFileStream_ptr[i_layer].read((char*)f_ptr,sizeof(float));
		    if(f_value<-9000)
		      {
		      f_transval=(int)f_nodata;
		      }
		    else
		      {
		      f_transval= GetExtraSplineValue(f_value,
                                          i_layer);
		      } // end else valid value
		    FLTOutStream_ptr[i_layer].write((char*)f_ptr2,i_float_size);
		    }// end for i_layer
	    //~~//~~//~~//~~//~~//
	    if(i_geo>0)
	      {
	      // use the cells value read in from the last grid 
	      f_ptr3=&(f_xtransval);
	      f_ptr4=&(f_ytransval);
	      if(f_value<-9000)
	        {
	        f_xtransval=(int)f_nodata;
	        f_ytransval=(int)f_nodata;
	        }
	      else
	        {
	        f_xtransval= GetSplineValueGeo(f_xdist);
	        f_ytransval= GetSplineValueGeo(f_ydist);	        
	        } // end else valid value
	      FLTGeoOutStream_ptr[0].write((char*)f_ptr3,i_float_size);
	      FLTGeoOutStream_ptr[1].write((char*)f_ptr4,i_float_size);
	      }// end if geo >0
	    //~~//~~//~~//~~//~~//
	    f_xdist += f_cellsize;
	    } // end for i_col
    f_ydist += f_cellsize;
	  }// end for i_row
  
  // Output an end of file token
  i_s=K_EOF;
  i_s_ptr=&i_s;

  for(i_layer=0;i_layer<n_layers;i_layer++)
	  {
	  InFileStream_ptr[i_layer].close();
 //	 FLTOutStream_ptr[i_layer].write((char*)i_s_ptr,i_short_size);   // eoftoken
	  FLTOutStream_ptr[i_layer].close();
	  }

  //~~//~~//~~//~~//~~//
    if(i_geo>0)
    {   
    FLTGeoOutStream_ptr[0].close();
    FLTGeoOutStream_ptr[1].close();
    } // end if geo
  //~~//~~//~~//~~//~~// 

  
  delete [] InFileStream_ptr;
  delete [] FLTOutStream_ptr;
  return i_result_code;
}// end ConvertGridStackToTransGridsFLOATBinaryFiles
//************************************************************
// TransCompact
// takes a lot of binary grid files, transforms them and squishes them a lot
// floats are converted to unsigned short ints
// the format is RUN LENGTH : datapoint 1 layer 1, datapoint1 layer 2 ..datapoint1 layer n:
// datapoint 2 layer 1, datapoint2 layer 2 ..datapoint2 layer n:  to datapoint n layer n
// NODATA runs are given as NEGATIVES, with no data afterwards.
// e.g. -5 5 34 23 24 25 21 23 25 12 25 25 -40 for a 2 layer grid.. is
// 5 nodata cells, 5 data cells, first cell, 34, 23 etc, followed by 40 nodatas
//****************************************************************
/*
int TTranspactor::LoadGlobalGDMParameters(string arg_s_path)
{
//AnsiString as_floater;
ifstream* InFileStream_ptr;
ifstream InParamsFile;
int b_read,b_skip;
string s_buffer;
string s_PredSplValTest;
string s_PredCoefTest;
string s_PredSplVal= "PredSplVal";
string s_PredCoef= "PredCoef";
string s_left;
int n_Knots[18];
double d_sum_coeffs;
//double d_Coeff[18][4];
double d_test;
const char* c_string;
int i_char;
char c_equals= '=';
int i_file_layer; //what is the order of the layer in this file?
int i_alpha_layer=0;  // which layer are we when alphabetically ordered as in bulk compactor: currently hard coded for gdm in i_which_alpha_layer
int i_which_alpha_layer[18]={3,9,12,13,15,16,17,14,10,1,2,8,4,11,5,6,7,0 };//{2,8,11,12,14,15,16,13,9,0,1,7,3,10,4,5,6};//{9,10,0,12,14,15,16,11,1,8,13,2,3,7,4,5,6};  // conversion: put i_file_layer in as index and get out an alphabetically ordered layer
int i_knot=0;
int b_success=1; //unless otherwise!

  //Zero results arrays
  for(i_alpha_layer=0;i_alpha_layer<18;i_alpha_layer++)
	{
	n_splines[i_alpha_layer]=0;
	for (i_knot = 0; i_knot < 5; i_knot++)
	  {
	  d_quant[i_alpha_layer][i_knot]=0;
	  d_coeff[i_alpha_layer][i_knot]=0;
	  }
	} // end for i_aphpa_layer

  InParamsFile.open(arg_s_path.c_str());
  // Most of the parameters file is not useful here. We just want the model parameters:
  // For each of 3-4 splines we have a
  // PredSplVal1.1=0.3448585858858585858
  // PredSplVal1.2=0.67484748484848484848
  // PredCoef1.1=89;
  // Predcoef1.2=32.44
  //So we need to run down til we hit a PredSplVal then read a block.
  // Stop reading when when?
  i_file_layer=-1;// increment to 0 when we hit a block
  if(InParamsFile.is_open())
	{
	b_read=1;
	b_skip=0;
	while(b_read>0)
	  {
	  if(!b_skip)
		InParamsFile>>s_buffer;
	  b_skip=0;
	  if(IsQuant(s_buffer))
		{
		//*********start of block for this layer *************************
		i_file_layer++; // note initialised to -1 so first is zero
		i_alpha_layer=i_which_alpha_layer[i_file_layer]; // convert to alphabetically ordered layers
		i_knot=0;
		d_sum_coeffs=0;    // check for zero zum coeffs
		//Now work through the Quants..SplValues,(s_buffer counting knots until we hit a coeff
		while(IsQuant(s_buffer))
		  {
		  // will stay in her until end of PredSplVals
		  d_quant[i_alpha_layer][i_knot]=GetQuant(s_buffer);
		  i_knot++;
		  // read next line
		  InParamsFile>>s_buffer;
		  }  // end while is quant
		 //How many knots were there?
		 n_splines[i_alpha_layer]=i_knot;
		 b_active[i_alpha_layer]=1;
		i_knot=0;// reset for Coeffs
		//If all is well, we will have just read a nice PredCoeffsline
		while(IsCoeff(s_buffer))
		  {
		  // will stay in her until end of Coeffs.. a blank line
          // will stay in her until end of PredSplVals
		  d_coeff[i_alpha_layer][i_knot]=GetCoeff(s_buffer);
		  d_sum_coeffs+=d_coeff[i_alpha_layer][i_knot];
		  i_knot++;
		  // read next line
		  InParamsFile>>s_buffer;
		  }   // end while iscoeff
		//If all is well, we wil have just read the next PredSplVal!
		//End of block..flag
		b_skip=1;
		if(d_sum_coeffs<=0)
		  b_iszero[i_alpha_layer]=1;
		}// end if ISQuant (outer at start of block) so just back to the main read loop
	  else
		{
		// If not the start of a block
		if(IsTRANSPREDS(s_buffer))
		  b_read=0; //false..stop reading the file
        }

	  } // end while b_read

	}  // end if file open

  if(InParamsFile.is_open())
	InParamsFile.close();
  return b_success;
}
*/
/******************************************************************************
 * TTranspactor  SetSplines. Loads local spline values here to reduce arguments
 * herein and allow precalculation of extrapolation thang
 *********************************************************************/
void TTranspactor::SetSplines( const int arg_i_layer,           // layer to set splines for (default just use 0)
						   const int arg_n_splines,          // number of splines
						   const double *arg_d_coeffs,            //ptr to array of coeffs[arg_n_splines]
						   const double *arg_d_quants)            // ptr to array of quants[arg_n_splines]
{
int i_spline; //counter

   n_splines[arg_i_layer]=arg_n_splines;
   b_active[arg_i_layer]=1;
   for(i_spline=0;i_spline<n_splines[arg_i_layer];i_spline++)
	 {
	 d_quant[arg_i_layer][i_spline]=arg_d_quants[i_spline];  //Assume nobody uses 5 splines for overfitting reasons
	 d_coeff[arg_i_layer][i_spline]=arg_d_coeffs[i_spline];  //Assume nobody uses 5 splines for overfitting reasons

	 }// end for i_spline
}// end func set splines
/******************************************************************************
 * TTranspactor  Set Up Extrapolation
 * Works out extrapolations for this layer.
 * tests top and bottom separately.. is full extrapolation or edge 10% extrapolation weakest
 *
 * sets global variables d_m[] and d_c[]
 * these values take full account
 *********************************************************************/
void TTranspactor::SetUpExtrapolation(const int arg_i_layer) // set up linear parameters
{
double d_xrange; //bottom to top of fitted
double d_val_yrange; //bottom trans minus top trans
double d_bot_x; //local copies of spline lower quant [0]
double d_top_x; //local copies of spline upper quant  [nsplines-1]
double d_bot_val_y, d_top_val_y; //gdm tranformed values
double d_slope; //slope of the line
double d_intersect; // point at which line crosses y axis
double d_full_m;   // slope for full fit
double d_full_c;   //const for full fit
double d_ten_m;     //10% extrap slope
double d_ten_c;     //10% extrap constant
double d_bx,d_by,d_tx,d_ty;   //Values 10% inside and out for either top and bottom
double d_full_y; //test values lying outside spline for full model
int i;

  //First big line across full fitted range***********************
  //get x
  d_bot_x=d_quant[arg_i_layer][0];
  d_top_x=d_quant[arg_i_layer][n_splines[arg_i_layer]-1];
  //Now get GDM transformed values for top and bottom
  d_bot_val_y=GetSplineValue(d_bot_x,     // value to transform
							 arg_i_layer);
  d_top_val_y=GetSplineValue(d_top_x,     // value to transform
						   arg_i_layer);
  //Line fit calc.....................
  d_xrange=d_top_x-d_bot_x;
  d_val_yrange=d_top_val_y-d_bot_val_y;
  d_slope=d_val_yrange/d_xrange;
  d_intersect= d_bot_val_y-(d_bot_x*d_slope);     //Use point on lime (bot) to calculate offset
  //copy
  d_full_m=d_slope;
  d_full_c=d_intersect;

  // Now the top end************************************
  //get x
  d_bx=d_top_x-(d_xrange/10);  // Value 10% of range to calculate slope
  // get
  d_by=GetSplineValue(d_bx,
					 arg_i_layer);
  //Line fit calc..........................
  d_xrange=d_top_x-d_bx;
  d_val_yrange=d_top_val_y-d_by;
  d_slope=d_val_yrange/d_xrange;
  d_intersect= d_top_val_y-(d_top_x*d_slope);  //Use point on lime (top) to calculate offset
  //copy
  d_ten_m=d_slope;
  d_ten_c=d_intersect;
  //Calculate test point using 10% extrapolation
  d_tx=d_top_x+(d_xrange/10);  //Value 10% outside to calculate outer yval.
  d_ty=(d_ten_m*d_tx)+d_ten_c;  //10% extrap value
  d_full_y=(d_full_m*d_tx)+d_full_c;//full line value
  //Now do something with this info..check extrapolation setting
  switch(i_extrapolation)
	{
	case 0: // no extrapolation.. cap at ends
	  //Top end.[1]. cap at sum coefficients
	  d_m[arg_i_layer][1]=0; // leaving just c
	  d_c[arg_i_layer][1]=0; // set c to flat sum coefficients
	  for(i=0; i<n_splines[arg_i_layer]; i++)
		d_c[arg_i_layer][1]+=d_coeff[arg_i_layer][i];
	  break;
	case 1: // minimum linear extrapolation
	  //Top end: whichever gives the lowest
	  if(d_full_y<d_ty)
		{
		d_m[arg_i_layer][1]=d_full_m; // use full val
		d_c[arg_i_layer][1]=d_full_c; // set c to full val
		}    //end if full line is minimum
	  else
		{
		d_m[arg_i_layer][1]=d_ten_m; // use 10% val
		d_c[arg_i_layer][1]=d_ten_c; // set c to 10% val
		} // end if else 10% line is minimum
	  break;
	case 2:
	  //10%
		d_m[arg_i_layer][1]=d_ten_m; // use 10% val
		d_c[arg_i_layer][1]=d_ten_c; // set c to 10% val

	  break;
	case 3:
	default:
		//full
		d_m[arg_i_layer][1]=d_full_m; // use full val
		d_c[arg_i_layer][1]=d_full_c; // set c to full val
	  break;
	}  // end switch
  //**********************************************************


  // ** Now the bottom end************************************
  //get x
  d_tx=d_bot_x+(d_xrange/10);  // Value 10% of range above bottom to calculate slope   (tx becuase other one will be lower)
  // get
  d_ty=GetSplineValue(d_tx,
					  arg_i_layer);
  //Line fit calc..........................
  d_xrange=d_tx-d_bot_x;
  d_val_yrange=d_ty-d_bot_val_y;
  d_slope=d_val_yrange/d_xrange;
  d_intersect= d_bot_val_y-(d_bot_x*d_slope);  //Use point on lime (top) to calculate offset
  //copy
  d_ten_m=d_slope;
  d_ten_c=d_intersect;
  //Calculate test point using 10% extrapolation
  d_bx=d_bot_x-(d_xrange/10);  //Value 10% outside to calculate outer yval.
  d_by=(d_ten_m*d_bx)+d_ten_c;  //10% extrap value
  d_full_y=(d_full_m*d_bx)+d_full_c;//full line value
  //Now do something with this info..check extrapolation setting
  switch(i_extrapolation)
	{
	case 0: // no extrapolation.. cap at ends
	  //Top end.[1]. cap at sum coefficients
	  d_m[arg_i_layer][0]=0; // leaving just c
	  d_c[arg_i_layer][0]=0; // set c to flat zero
	  break;
	case 1: // minimum linear extrapolation
	  //bottom end: whichever gives the highest
	  if(d_full_y>d_by)
		{
		d_m[arg_i_layer][0]=d_full_m; // use full val
		d_c[arg_i_layer][0]=d_full_c; // set c to full val
		}    //end if full line is minimum
	  else
		{
		d_m[arg_i_layer][0]=d_ten_m; // use 10% val
		d_c[arg_i_layer][0]=d_ten_c; // set c to 10% val
		} // end if else 10% line is minimum
	  break;
	case 2:
	  //10%
		d_m[arg_i_layer][0]=d_ten_m; // use 10% val
		d_c[arg_i_layer][0]=d_ten_c; // set c to 10% val

	  break;
	case 3:
	default:
		//full
		d_m[arg_i_layer][0]=d_full_m;// Bug? changed from...  d_m[arg_i_layer][1]=d_full_m; // use full val           
	  d_c[arg_i_layer][0]=d_full_c;// Bug? changed from...  d_c[arg_i_layer][1]=d_full_c; // set c to full val
	  break;
	}  // end switch
}   // end func set up extrapolation
//***************************************************************************
// Set up file paths
// Deals with the multiple file path thign required for the global modelling
// Hardcoded for global
//**************************************************************************
/*
void TTranspactor::SetUpGlobFilePaths(const string arg_s_clim_path,
									  const string arg_s_land_path,
									  const string arg_s_subs_path,
									  const string arg_s_resid_path)
{
   s_binary_in[0]=arg_s_subs_path+"BARE.flt";
   s_binary_in[1]=arg_s_subs_path+"BD.flt";
   s_binary_in[2]=arg_s_subs_path+"CLAY.flt";
   s_binary_in[3]=arg_s_clim_path+"EAAS.flt";
   s_binary_in[4]=arg_s_clim_path+"EPI.flt";
   s_binary_in[5]=arg_s_resid_path+"MDSax1.flt";
   s_binary_in[6]=arg_s_resid_path+"MDSax2.flt";
   s_binary_in[7]=arg_s_resid_path+"MDSax3.flt";
   s_binary_in[8]=arg_s_subs_path+"PH.flt";
   s_binary_in[9]=arg_s_clim_path+"PTA.flt";
   s_binary_in[10]=arg_s_land_path+"RUG.flt";
   s_binary_in[11]=arg_s_subs_path+"SILT.flt";
   s_binary_in[12]=arg_s_clim_path+"TNI.flt";
   s_binary_in[13]=arg_s_clim_path+"TRX.flt";
   s_binary_in[14]=arg_s_land_path+"TWI_SAGE.flt";
   s_binary_in[15]=arg_s_clim_path+"TXX.flt";
   s_binary_in[16]=arg_s_clim_path+"WDI.flt";
   s_binary_in[17]=arg_s_clim_path+"WDX.flt";
}// end func set up file paths
 */
//***************************************************************************
// Set up output float file paths
//
// Hardcoded for global  for now
//**************************************************************************
/*
void TTranspactor::SetUpGlobOutFLTFilePaths(const string arg_s_res_path)
{
   s_binary_out[0]=arg_s_res_path+"BAREtran.flt";
   s_binary_out[1]=arg_s_res_path+"BDtrans.flt";
   s_binary_out[2]=arg_s_res_path+"CLAYtrans.flt";
   s_binary_out[3]=arg_s_res_path+"EAAStrans.flt";
   s_binary_out[4]=arg_s_res_path+"EPItrans.flt";
   s_binary_out[5]=arg_s_res_path+"MDSax1trans.flt";
   s_binary_out[6]=arg_s_res_path+"MDSax2trans.flt";
   s_binary_out[7]=arg_s_res_path+"MDSax3trans.flt";
   s_binary_out[8]=arg_s_res_path+"PHtrans.flt";
   s_binary_out[9]=arg_s_res_path+"PTAtrans.flt";
   s_binary_out[10]=arg_s_res_path+"RUGtrans.flt";
   s_binary_out[11]=arg_s_res_path+"SILTtrans.flt";
   s_binary_out[12]=arg_s_res_path+"TNItrans.flt";
   s_binary_out[13]=arg_s_res_path+"TRXtrans.flt";
   s_binary_out[14]=arg_s_res_path+"TWItrans.flt";
   s_binary_out[15]=arg_s_res_path+"TXXtrans.flt";
   s_binary_out[16]=arg_s_res_path+"WDItrans.flt";
   s_binary_out[17]=arg_s_res_path+"WDXtrans.flt";
}// end func set up file paths
 */
//***************************************************************************//
// Get Extra (polated) Spline Value
// Calculate the GDM transformed value for an environmental variable
// Applies any preset extrapolations assuming all nicely set up for this layer
double TTranspactor::GetExtraSplineValue(const double arg_d_value,     // value to transform
										 const int arg_i_layer)       // layer to use 0 by default
{
double d_result;
  //All layers and extrapolation types are the same
  // Check for underrun..use bot extrapolation [0]
  // Check for overrun.. use top extrap [1]
  // Otherwise just normal
  if(arg_d_value<d_quant[arg_i_layer][0])
	{
	//Underrun
	d_result=d_m[arg_i_layer][0]*arg_d_value+ d_c[arg_i_layer][0];
	} // end if underrun
  else
	{
	if(arg_d_value>d_quant[arg_i_layer][n_splines[arg_i_layer]-1])
	  {
	  //Overrun
	  d_result=d_m[arg_i_layer][1]*arg_d_value+ d_c[arg_i_layer][1];
	  }
	else
	  {
	  //Normal within bounds
	  d_result=GetSplineValue(arg_d_value,
							  arg_i_layer);
      }  // end if no extrapolation
	} // end if not underrun
  return d_result;
}   // end func get extrapolated Spline Value
//***************************************************************************//
// Get Extrapolation
// Calculate the amount of transformed extrapolation from the min-max fitted space
// Like the Manion ERR grids
// Applies any preset extrapolations assuming all nicely set up for this layer
double TTranspactor::GetExtrapolation(const double arg_d_value,     // value to transform
									  const int arg_i_layer)       // layer to use 0 by default
{
double d_result;
  //All layers and extrapolation types are the same
  // Check for underrun..use bot extrapolation [0]
  // Check for overrun.. use top extrap [1]
  // Otherwise just normal
  if(arg_d_value<d_quant[arg_i_layer][0])
	{
	//Underrun
	d_result=d_m[arg_i_layer][0]*arg_d_value+ d_c[arg_i_layer][0]-d_quant[arg_i_layer][0];
	} // end if underrun
  else
	{
	if(arg_d_value>d_quant[arg_i_layer][n_splines[arg_i_layer]-1])
	  {
	  //Overrun
	  d_result=d_m[arg_i_layer][1]*arg_d_value+ d_c[arg_i_layer][1]-d_quant[arg_i_layer][n_splines[arg_i_layer]-1];
	  }
	else
	  {
	  //Normal within bounds
	  d_result=0;
      }  // end if no extrapolation
	} // end if not underrun
  return d_result;
}   // end func get extrapolation
//***************************************************************************//
// Local var function for original GDM transforms
// Calculate the GDM transformed value for an environmental variable
double TTranspactor::GetSplineValue(const double arg_d_value,     // value to transform
                                    const int arg_i_layer)       // layer to use 0 by default
{
// declare local variables
int i;
double d_result = 0.0;
//double d_pQuants[3]= {1.509581, 7526.07642, 503404.0116};
//double d_pCoeffs[3]= {0.598526, 1.843695, 1.340533};
double* d_qptr;
double* d_cptr;

	   d_qptr=&(d_quant[arg_i_layer][0]);
	   d_cptr=&(d_coeff[arg_i_layer][0]);

	// loop over the splines
	// First spline
	for(i=0; i<n_splines[arg_i_layer]; i++)
	  {
	  if(i == 0) // first spline
		{
		d_result += DoSplineCalc( arg_d_value,
								  d_qptr[i],
								  d_qptr[i],
								  d_qptr[i+1] ) * d_cptr[i];
		} // end if i ==0
	  else if(i == n_splines[arg_i_layer]-1) // last spline
		{
		d_result += DoSplineCalc( arg_d_value,
								  d_qptr[i-1],
								  d_qptr[i],
								  d_qptr[i] ) * d_cptr[i];
		} // end else if
	  else  // a middle spline
		{
		d_result += DoSplineCalc( arg_d_value,
								  d_qptr[i-1],
								  d_qptr[i],
								  d_qptr[i+1] ) * d_cptr[i];
		} // end else
	  } // end for i
	return(d_result);
} // end GetSplineValue

//***************************************************************************//
// Local var function for original GDM transforms
// Calculate the GDM transformed value for an environmental variable
double TTranspactor::GetSplineValueGeo(const double arg_d_value)       // layer to use 0 by default
  {
    // declare local variables
    int i;
    double d_result = 0.0;
    double* d_qptr;
    double* d_cptr;
    
    d_qptr=&(d_geo_quant[0]);
    d_cptr=&(d_geo_coeff[0]);
    
    // loop over the splines
    // First spline
    for(i=0; i<i_geo_splines; i++)
      {
      if(i == 0) // first spline
        {
        d_result += DoSplineCalc(arg_d_value,
                                 d_qptr[i],
                                 d_qptr[i],
                                 d_qptr[i+1] ) * d_cptr[i];
        } // end if i ==0
      else if(i == i_geo_splines-1) // last spline
        {
        d_result += DoSplineCalc(arg_d_value,
                                 d_qptr[i-1],
                                 d_qptr[i],
                                 d_qptr[i] ) * d_cptr[i];
        } // end else if
      else  // a middle spline
        {
        d_result += DoSplineCalc(arg_d_value,
                                 d_qptr[i-1],
                                 d_qptr[i],
                                 d_qptr[i+1] ) * d_cptr[i];
      } // end else
    } // end for i
    return(d_result);
  } // end GetSplineValueGeo

//***************************************************************************//
// Calculate the GDM transformed value for an environmental variable
//  The value returned (d_result) is the gdm transformed value.
// Written by Glenn Manion.
//Note that for each knot you will need to reuse the respective quantiles
//Assuming 3 knots then for the 0% knot both q1 and q2 are the 0% quantile and q3 is the 50% quantile
// for the 50% knot q1 is the 0% quantile, q2 is the 50% quantile and q3 is the 100% quantile
// for the 100% knot q1 is the 50% quantile and q2 and q3 are the 100% quantile.
// Edited by Tom Harwood
//**Arguments*****************
//
//
double TTranspactor::CalculateSplineValue(const double arg_d_position,     // value to transform
									  const int arg_n_splines,          // number of splines
									  const double *arg_d_pCoeffs,            //ptr to array of coeffs[arg_n_splines]
									  const double *arg_d_pQuants)            // ptr to array of quants[arg_n_splines]
{
// declare local variables
int i;
double d_result = 0.0;
//double d_pQuants[3]= {1.509581, 7526.07642, 503404.0116};
//double d_pCoeffs[3]= {0.598526, 1.843695, 1.340533};

	// loop over the splines
	// First spline
	for(i=0; i<arg_n_splines; i++)
	  {
	  if(i == 0) // first spline
		{
		d_result += DoSplineCalc( arg_d_position, arg_d_pQuants[i], arg_d_pQuants[i], arg_d_pQuants[i+1] ) * arg_d_pCoeffs[i];
		} // end if i ==0
	  else if(i == arg_n_splines-1) // last spline
		{
		d_result += DoSplineCalc( arg_d_position, arg_d_pQuants[i-1], arg_d_pQuants[i], arg_d_pQuants[i] ) * arg_d_pCoeffs[i];
		} // end else if
	  else  // a middle spline
		{
		d_result += DoSplineCalc( arg_d_position, arg_d_pQuants[i-1], arg_d_pQuants[i], arg_d_pQuants[i+1] ) * arg_d_pCoeffs[i];
		} // end else
	  } // end for i
	return(d_result);
} // end CalculateTrans

//***************************************************************************//


//***************************************************************************//
// Calculate the I-Spline value for arg_d_val given quantiles q1, q2, q3
// Written by Glenn Manion
// Edited by Tom Harwood
//***************************************************************************//
double TTranspactor::DoSplineCalc(double arg_d_val,
							  double arg_d_q1,
							  double arg_d_q2,
							  double arg_d_q3)
{
	if ( arg_d_val <= arg_d_q1 )
	  return(0.0);
	else if ( arg_d_val >= arg_d_q3 ) return( 1.0 );
	else if ( ( arg_d_q1 < arg_d_val ) && ( arg_d_val < arg_d_q2 ) )
		return( ( ( ( arg_d_val - arg_d_q1 ) * ( arg_d_val - arg_d_q1 ) ) / ( ( arg_d_q2 - arg_d_q1 ) * ( arg_d_q3 - arg_d_q1 ) ) ) );

	else
		return( ( 1.0 - ( ( ( arg_d_q3 - arg_d_val ) * ( arg_d_q3 - arg_d_val ) ) / ( ( arg_d_q3 - arg_d_q2 ) * ( arg_d_q3 - arg_d_q1) ) ) ) );
}  // End func DoSplineCalc
//***************************************************************************//

/******************************************************************************
 * TTranspactor  SetModelParameters. Loads local spline, quatile and coefficient
 * values here to reduce arguments and allow precalculation of extrapolation 
 *********************************************************************/
void TTranspactor::SetModelParameters(int arg_geo,
                                      int arg_n_layers,
                                      int* arg_n_splines,                         
                                      float* arg_knots,
                                      float* arg_coefficients)            
{
  int i_spline, i_prd; 
  int i_upto=0;
  double d_geo_splines;
  
  // If geographic distance is included, move this to the end.
  if(arg_geo>0)
    {
    d_geo_splines=double(arg_n_splines[0]); 
    for(i_spline=0; i_spline<n_splines[i_prd]; i_spline++)
      {
      d_geo_quant[i_spline]=double(arg_knots[i_upto]);  
      d_geo_coeff[i_spline]=double(arg_coefficients[i_upto]);  
      i_upto += 1;
      }// end for i_spline
    for(i_prd=0; i_prd<arg_n_layers; i_prd++)
      {
      n_splines[i_prd]=double(arg_n_splines[(i_prd+1)]);
      for(i_spline=0; i_spline<n_splines[i_prd]; i_spline++)
        {
        d_quant[i_prd][i_spline]=double(arg_knots[i_upto]);  
        d_coeff[i_prd][i_spline]=double(arg_coefficients[i_upto]);  
        i_upto += 1;
        }// end for i_spline
      }// end for i_prd
    } // end if arg_geo>0
  else
    {
    for(i_prd=0; i_prd<arg_n_layers; i_prd++)
      {
      n_splines[i_prd]=double(arg_n_splines[i_prd]);
      for(i_spline=0; i_spline<n_splines[i_prd]; i_spline++)
        {
        d_quant[i_prd][i_spline]=double(arg_knots[i_upto]);  
        d_coeff[i_prd][i_spline]=double(arg_coefficients[i_upto]);  
        i_upto += 1;
        }// end for i_spline
      }// end for i_prd
    } // end else arg_geo>0
  
}// end func SetModelParameters
//******************************************************************************
  
  
  /******************************************************************************
   * TTranspactor  SetModelParameters. Loads local spline, quatile and coefficient
   * values here to reduce arguments and allow precalculation of extrapolation 
   *********************************************************************/
  void TTranspactor::SetGridPaths(int arg_geo,
                                  int arg_n_layers,
                                  std::string arg_InFilePaths[],
                                  std::string arg_OutFilePaths[])            
  {
    int i_prd; 

    for(i_prd=0; i_prd<arg_n_layers; i_prd++)
      {
      s_binary_in[i_prd] = arg_InFilePaths[i_prd];
      s_binary_out[i_prd] = arg_OutFilePaths[i_prd];
      } // end for i_prd
    
    // If geographic distance is included, add the filepaths for the x & y grids
    if(arg_geo>0)
      {
      s_geo_binary_out[0] = arg_OutFilePaths[arg_n_layers];
      s_geo_binary_out[1] = arg_OutFilePaths[(arg_n_layers+1)];      
      } // end if arg_geo>0
 
  }// end func SetGridPaths
  
  
    
//---------------------------------------------------------------------------

