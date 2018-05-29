//---------------------------------------------------------------------------
//*****************************************************************
// MURU COMPACTOR FUNCTIONS
// Set of routines to compact layers of transformed gdm grids
// or richness/condition grids to a floating point binary file
// with a form of run length encoding for nodata points
//*************
// Written by Tom Harwood CSIRO Ecosystem Sciences 2010
//*************
// Note that the input files for these functions are binary grids in
// the DIVA-GIS *.gri format, but any compatible floating point grid format
// can be used.
// Grid dimensions are passed to the functions, to avoid issues reading headers,
// and at present any value less than -500 is assumed to be a no_data point.
// In the GDM application this will not cause an issue, and works well with a
// standard nodata value of -9999.
//*******************************************************************
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include "defines.h"         // all #defines

#include "DynamicArray.h"    // functions to ease allocation and destruction of dynamic arrays
#include "CompactorFunctions.h"

using namespace std;
//************************************************************
// Convert Grid Stack To Compact FLOAT Binary File
// takes a lot of binary grid files, and squishes them a lot
// floats are converted to unsigned short ints
// the format is RUN LENGTH : datapoint 1 layer 1, datapoint1 layer 2 ..datapoint1 layer n:
// datapoint 2 layer 1, datapoint2 layer 2 ..datapoint2 layer n:  to datapoint n layer n
// NODATA runs are given as NEGATIVES, with no data afterwards.
// e.g. -5 5 34 23 24 25 21 23 25 12 25 25 -40 for a 2 layer grid.. is
// 5 nodata cells, 5 data cells, first cell, 34, 23 etc, followed by 40 nodatas
//****************************************************************
extern int ConvertGridStackToCompactFloatBinaryFileWithZeroGrids(const long int arg_i_num_rows,
                                                    const long int arg_i_num_cols,
													const int arg_i_num_layers,
													string arg_s_layer_files[],
													string arg_s_res_path)
{
int i_result_code;
long int n_rows,n_cols; //local copies of arguments
int n_layers;
long int i_row, i_col;  // row and column counters
int i_layer;       //layer/file counter
float** f_data;
int i_pos; // position we have reached
ifstream* InFileStream_ptr;
ofstream OutFileStream;
//ofstream IndexFileStream;
//ofstream HeaderFileStream;
float* f_ptr;    // for reading data in
long int i_start_run, i_end_run;
short int si_run_length;
short* i_s_ptr;
short i_s;
float f_value;
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

  i_result_code=1; // good unless we hear otherwise!
  i_short_size= sizeof(short);
  i_float_size=sizeof(float);
  n_rows= arg_i_num_rows;
  n_cols= arg_i_num_cols;
  n_layers= arg_i_num_layers;

  // Allocate a set of arrays to hold the lines for all layers
  f_data= AllocateFloat2DArray(n_layers,
                               n_cols+1);

  // Open up all grid files to compact
  InFileStream_ptr = new ifstream[n_layers];
  for(i_layer=0;i_layer<n_layers;i_layer++)
	{
	InFileStream_ptr[i_layer].open(arg_s_layer_files[i_layer].c_str(), ios::in |ios::binary);
	// Check file is present
	if(!InFileStream_ptr[i_layer].is_open())
	  {
	  //Problem.. assume that there is a Zero grid rather than a Trans grid
	  s_path=arg_s_layer_files[i_layer];
	  i_length=s_path.length();
	  s_path.resize(i_length-8);     // take off "Tran.flt"
	  s_path=s_path+"Zero.flt";     // put on
	  InFileStream_ptr[i_layer].clear();
	  InFileStream_ptr[i_layer].open(s_path.c_str(), ios::in |ios::binary);
	  }// end if not open
	// If still missing, tell the user to get their act together!
	if(!InFileStream_ptr[i_layer].is_open())
	  return 0;
    }
  // Open big results file AS BINARY
  OutFileStream.open(arg_s_res_path.c_str(), ios::out |ios::binary);
  // build index file name ( to hold seekposes for the beginning of each line)


  d_spos=0; // zero seekpos counter
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
        d_spos+=i_short_size;
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
        // we have found the end of the run.

        //Go back to start of run and read out ALL LAYERS
        // Write out run length first
        i_s_ptr=&(si_run_length);
        OutFileStream.write((char*)i_s_ptr,i_short_size);   // positive run length
        d_spos+=i_short_size;
        // Read out whoel run, all layer data for each point
        // Output data for each Grid cell for all layers
        // in order of the layers.. e.g. for 3 layers 1, 2,3, with run length 2
        // 2(runlength) 1 2 3 (firstcell) 1 2 3 (second cell)
        i_end_run=i_start_run+(int)si_run_length;
        for(i_col=i_start_run;i_col<i_end_run;i_col++)
          {
          for(i_layer=0;i_layer<n_layers;i_layer++)
            {
            f_value=f_data[i_layer][i_col];
            f_ptr=&(f_value);
            OutFileStream.write((char*)f_ptr,i_float_size);   // unshortened value for this layer
            d_spos+=i_float_size;
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

  // Free up arrays
  FreeFloat2DArray(f_data);
  delete [] InFileStream_ptr;

  return i_result_code;
}// end ConvertGridStackToCompactFLOATBinaryFile


//************************************************************
// Convert Grid Stack To Compact FLOAT Binary File
// takes a lot of binary grid files, and squishes them a lot
// floats are converted to unsigned short ints
// the format is RUN LENGTH : datapoint 1 layer 1, datapoint1 layer 2 ..datapoint1 layer n:
// datapoint 2 layer 1, datapoint2 layer 2 ..datapoint2 layer n:  to datapoint n layer n
// NODATA runs are given as NEGATIVES, with no data afterwards.
// e.g. -5 5 34 23 24 25 21 23 25 12 25 25 -40 for a 2 layer grid.. is
// 5 nodata cells, 5 data cells, first cell, 34, 23 etc, followed by 40 nodatas
//****************************************************************
extern int ConvertGridStackToCompactFloatBinaryFile(const long int arg_i_num_rows,
                                                    const long int arg_i_num_cols,
													const int arg_i_num_layers,
													string arg_s_layer_files[],
													string arg_s_res_path)
{
int i_result_code;
long int n_rows,n_cols; //local copies of arguments
int n_layers;
long int i_row, i_col;  // row and column counters
int i_layer;       //layer/file counter
float** f_data;
int i_pos; // position we have reached
ifstream* InFileStream_ptr;
ofstream OutFileStream;
ofstream IndexFileStream;
ofstream HeaderFileStream;
float* f_ptr;    // for reading data in
long int i_start_run, i_end_run;
short int si_run_length;
short* i_s_ptr;
short i_s;
float f_value;
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

  i_result_code=1; // good unless we hear otherwise!
  i_short_size= sizeof(short);
  i_float_size=sizeof(float);
  n_rows= arg_i_num_rows;
  n_cols= arg_i_num_cols;
  n_layers= arg_i_num_layers;

  // Allocate a set of arrays to hold the lines for all layers
  f_data= AllocateFloat2DArray(n_layers,
                               n_cols+1);

  // Open up all grid files to compact
  InFileStream_ptr = new ifstream[n_layers];
  for(i_layer=0;i_layer<n_layers;i_layer++)
    {
    InFileStream_ptr[i_layer].open(arg_s_layer_files[i_layer].c_str(), ios::in |ios::binary);
    }
  // Open big results file AS BINARY
  OutFileStream.open(arg_s_res_path.c_str(), ios::out |ios::binary);
  // build index file name ( to hold seekposes for the beginning of each line)
  s_path=arg_s_res_path;
  i_length=s_path.length();
  s_path.resize(i_length-4);
  s_name=".mci";
  s_index_path=s_path+s_name;
  IndexFileStream.open(s_index_path.c_str(), ios::out |ios::binary);

  s_path=s_path+".mch";
  HeaderFileStream.open(s_path.c_str(), ios::out |ios::binary);
  l_ptr=&n_cols;
  HeaderFileStream.write((char*)l_ptr,sizeof(long int));
  l_ptr=&n_rows;
  HeaderFileStream.write((char*)l_ptr,sizeof(long int));
  i_ptr=&n_layers;
  HeaderFileStream.write((char*)i_ptr,sizeof(int));
  i_data_format=K_FLOAT;
  i_ptr=&i_data_format;
  HeaderFileStream.write((char*)i_ptr,sizeof(int));
  HeaderFileStream.close();

  d_spos=0; // zero seekpos counter
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
        d_spos+=i_short_size;
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
        // we have found the end of the run.

        //Go back to start of run and read out ALL LAYERS
        // Write out run length first
        i_s_ptr=&(si_run_length);
        OutFileStream.write((char*)i_s_ptr,i_short_size);   // positive run length
        d_spos+=i_short_size;
        // Read out whoel run, all layer data for each point
        // Output data for each Grid cell for all layers
        // in order of the layers.. e.g. for 3 layers 1, 2,3, with run length 2
        // 2(runlength) 1 2 3 (firstcell) 1 2 3 (second cell)
        i_end_run=i_start_run+(int)si_run_length;
        for(i_col=i_start_run;i_col<i_end_run;i_col++)
          {
          for(i_layer=0;i_layer<n_layers;i_layer++)
            {
            f_value=f_data[i_layer][i_col];
            f_ptr=&(f_value);
            OutFileStream.write((char*)f_ptr,i_float_size);   // unshortened value for this layer
            d_spos+=i_float_size;
            }// end for i_layer
          }// end for i_col start to end run
        } // end DATA run
      }
    // At the end of the row, output the offset from the start of the file
    d_ptr=&d_spos;
    IndexFileStream.write((char*)d_ptr,sizeof(double));

    }// end for i_row
  // Output an end of file token
  i_s=K_EOF;
  i_s_ptr=&i_s;
  OutFileStream.write((char*)i_s_ptr,i_short_size);   // eoftoken
  //Close files
  OutFileStream.close();
  IndexFileStream.close();

  for(i_layer=0;i_layer<n_layers;i_layer++)
     InFileStream_ptr[i_layer].close();

  // Free up arrays
  FreeFloat2DArray(f_data);
  delete [] InFileStream_ptr;

  return i_result_code;
}// end ConvertGridStackToCompactFLOATBinaryFileWithZERO
//************************************************************
// Convert Condition Grid  To Compact FLOAT Binary File
// takes a binary grid file in DIVA GIS .gri format, and squishes them a lot
// floats are floats
// the format is RUN LENGTH : datapoint 1 layer 1, datapoint1 layer 2 ..datapoint1 layer n:
// datapoint 2 layer 1, datapoint2 layer 2 ..datapoint2 layer n:  to datapoint n layer n
// NODATA runs are given as NEGATIVES, with no data afterwards.
// e.g. -5 5 34 23 24 25 21 23 25 12 25 25 -40 for a 2 layer grid.. is
// 5 nodata cells, 5 data cells, first cell, 34, 23 etc, followed by 40 nodatas
//****************************************************************
extern int ConvertConditionGridToCompactFloatBinaryFile(const long int arg_i_num_rows,
                                                        const long int arg_i_num_cols,
                                                        string arg_s_condition_file,
                                                        string arg_s_res_path)
{
int i_result_code;
long int n_rows,n_cols; //local copies of arguments
int n_layers;
long int i_row, i_col;  // row and column counters
float* f_data;
//float** f_cell_pptr;     // Pointer to float pointer to make array of a pointer for each layer
long int i_pos; // position we have reached
ifstream InFileStream;
ofstream OutFileStream;
//ofstream IndexFileStream;
//ofstream HeaderFileStream;
float* f_ptr;    // for reading data in
long int i_start_run, i_end_run;
short int si_run_length;
short* i_s_ptr;
float f_value;
string s_path,s_index_path;
int i_short_size;
double d_spos; // total offset of the start of each new line.
int i_length;
int* i_ptr;
int i_float_size;
int i_data_format;
double* d_ptr;
long int* l_ptr;

  i_result_code=1; // good unless we hear otherwise!
  i_short_size= sizeof(short);
  i_float_size=sizeof(float);
  n_rows= arg_i_num_rows;
  n_cols= arg_i_num_cols;

  // Allocate an array to hold the lines
  f_data= AllocateFloat1DArray(n_cols+1);

  // Open up condition grid file to compact
  InFileStream.open(arg_s_condition_file.c_str(), ios::in |ios::binary);
  // Open results file AS BINARY
  OutFileStream.open(arg_s_res_path.c_str(), ios::out |ios::binary);
  // build index file name ( to hold seekposes for the beginning of each line)
//  s_path=arg_s_res_path;
//  i_length=s_path.length();
//  s_path.resize(i_length-4);
//  s_index_path=s_path+".mci";
//  IndexFileStream.open(s_index_path.c_str(), ios::out |ios::binary);

  // Quickly open and close and output header file details: cols. rows. layers.
//  s_path=s_path+".mch";
//  HeaderFileStream.open(s_path.c_str(), ios::out |ios::binary);
//  l_ptr=&n_cols;
//  HeaderFileStream.write((char*)l_ptr,sizeof(long int));
//  l_ptr=&n_rows;
//  HeaderFileStream.write((char*)l_ptr,sizeof(long int));
//  i_ptr=&n_layers;
//  n_layers=1;
//  HeaderFileStream.write((char*)i_ptr,sizeof(int));
//  i_data_format=K_FLOAT;
//  i_ptr=&i_data_format;
//  HeaderFileStream.write((char*)i_ptr,sizeof(int));
//  HeaderFileStream.close();

  d_spos=0; // zero seekpos counter
  //Read whole lines into buffers
  for(i_row=0;i_row<n_rows;i_row++)
    {
    // read in whole line and process
    for(i_col=0;i_col<n_cols;i_col++)
      {
      // Read in value for each cell
	  f_ptr=&(f_data[i_col]);
	  //Read to f_data array
	  InFileStream.read((char*)f_ptr,sizeof(float));
      } // end for i_col
    i_pos=0;
    i_col=0;
    //Work through the whole row
    while(i_pos<n_cols)
      {
      // split into nodata or data runs
      if(f_data[i_pos]<K_NO_DATA_ISH)
        {
        //NO_DATA
        //  i_col=i_pos;
        i_start_run=i_pos;
        si_run_length=0;
        while(f_data[i_pos]<K_NO_DATA_ISH&&i_pos<n_cols)
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
        d_spos+=i_short_size;
        } // end NODATA run
      else
        {
        //DATA
       // i_col=i_pos;
        i_start_run=i_pos;
        si_run_length=0;
        while(f_data[i_pos]>=K_NO_DATA_ISH&&i_pos<n_cols)
          {
          si_run_length++;
          i_pos++;
          if(si_run_length>=K_SHORT_MAX) // don't run off the end of a short for really big grids
            break;
          }// end while still data
        // we have found the end of the run.

        //Go back to start of run and read out ALL LAYERS
        // Write out run length first
        i_s_ptr=&(si_run_length);
        OutFileStream.write((char*)i_s_ptr,i_short_size);   // positive run length
        d_spos+=i_short_size;
        // Read out whoel run, all layer data for each point
        // Output data for each Grid cell for all layers
        // in order of the layers.. e.g. for 3 layers 1, 2,3, with run length 2
        // 2(runlength) 1 2 3 (firstcell) 1 2 3 (second cell)
        i_end_run=i_start_run+(int)si_run_length;
        for(i_col=i_start_run;i_col<i_end_run;i_col++)
          {
          f_value=f_data[i_col];
          f_ptr=&f_value;
          OutFileStream.write((char*)f_ptr,i_float_size);   // value for this layer
          d_spos+=i_float_size;
          }// end for i_col start to end run
        } // end DATA run
      } // end while
    // At the end of the row, output the offset from the start of the file
//    d_ptr=&d_spos;
//    IndexFileStream.write((char*)d_ptr,sizeof(double));
    }// end for i_row

  // Output an end of file token
  f_value=K_EOF;
  f_ptr=&f_value;
  OutFileStream.write((char*)f_ptr,i_float_size);   // eoftoken
  //Close files
  OutFileStream.close();
//  IndexFileStream.close();
  InFileStream.close();

  // Free up arrays
  FreeFloat1DArray(f_data);

  return i_result_code;
}// end ConvertConditionGridToCompactBinaryFile
//************************************************************
// Convert IntCondition Grid  To Compact FLOAT Binary File
// takes a binary grid file in DIVA GIS .gri format, and squishes them a lot
// floats are floats
// the format is RUN LENGTH : datapoint 1 layer 1, datapoint1 layer 2 ..datapoint1 layer n:
// datapoint 2 layer 1, datapoint2 layer 2 ..datapoint2 layer n:  to datapoint n layer n
// NODATA runs are given as NEGATIVES, with no data afterwards.
// e.g. -5 5 34 23 24 25 21 23 25 12 25 25 -40 for a 2 layer grid.. is
// 5 nodata cells, 5 data cells, first cell, 34, 23 etc, followed by 40 nodatas
//****************************************************************
extern int ConvertIntConditionGridToCompactFloatBinaryFile(const long int arg_i_num_rows,
                                                        const long int arg_i_num_cols,
                                                        string arg_s_condition_file,
                                                        string arg_s_res_path)
{
int i_result_code;
long int n_rows,n_cols; //local copies of arguments
int n_layers;
long int i_row, i_col;  // row and column counters
float* f_data;
//float** f_cell_pptr;     // Pointer to float pointer to make array of a pointer for each layer
long int i_pos; // position we have reached
ifstream InFileStream;
ofstream OutFileStream;
ofstream IndexFileStream;
ofstream HeaderFileStream;
float* f_ptr;    // for reading data in
long int i_start_run, i_end_run;
short int si_run_length;
short* i_s_ptr;
float f_value;
string s_path,s_index_path;
int i_short_size;
double d_spos; // total offset of the start of each new line.
int i_length;
int* i_ptr;
int i_float_size;
int i_data_format;
double* d_ptr;
long int* l_ptr;
int i_value;

  i_result_code=1; // good unless we hear otherwise!
  i_short_size= sizeof(short);
  i_float_size=sizeof(float);
  n_rows= arg_i_num_rows;
  n_cols= arg_i_num_cols;

  // Allocate an array to hold the lines
  f_data= AllocateFloat1DArray(n_cols+1);

  // Open up condition grid file to compact
  InFileStream.open(arg_s_condition_file.c_str(), ios::in |ios::binary);
  // Open results file AS BINARY
  OutFileStream.open(arg_s_res_path.c_str(), ios::out |ios::binary);
  // build index file name ( to hold seekposes for the beginning of each line)
  s_path=arg_s_res_path;
  i_length=s_path.length();
  s_path.resize(i_length-4);
  s_index_path=s_path+".mci";
  IndexFileStream.open(s_index_path.c_str(), ios::out |ios::binary);

  // Quickly open and close and output header file details: cols. rows. layers.
  s_path=s_path+".mch";
  HeaderFileStream.open(s_path.c_str(), ios::out |ios::binary);
  l_ptr=&n_cols;
  HeaderFileStream.write((char*)l_ptr,sizeof(long int));
  l_ptr=&n_rows;
  HeaderFileStream.write((char*)l_ptr,sizeof(long int));
  i_ptr=&n_layers;
  n_layers=1;
  HeaderFileStream.write((char*)i_ptr,sizeof(int));
  i_data_format=K_FLOAT;
  i_ptr=&i_data_format;
  HeaderFileStream.write((char*)i_ptr,sizeof(int));
  HeaderFileStream.close();

  d_spos=0; // zero seekpos counter
  i_ptr=&i_value; // to read in
  //Read whole lines into buffers
  for(i_row=0;i_row<n_rows;i_row++)
    {
    // read in whole line and process
    for(i_col=0;i_col<n_cols;i_col++)
      {
      // Read in value for each cell
      InFileStream.read((char*)i_ptr,sizeof(int));
      //Read to f_data array
      f_data[i_col]=(float)i_value;
      } // end for i_col
    i_pos=0;
    i_col=0;
    //Work through the whole row
    while(i_pos<n_cols)
      {
      // split into nodata or data runs
      if(f_data[i_pos]<K_NO_DATA_ISH)
        {
        //NO_DATA
        //  i_col=i_pos;
        i_start_run=i_pos;
        si_run_length=0;
        while(f_data[i_pos]<K_NO_DATA_ISH&&i_pos<n_cols)
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
        d_spos+=i_short_size;
        } // end NODATA run
      else
        {
        //DATA
       // i_col=i_pos;
        i_start_run=i_pos;
        si_run_length=0;
        while(f_data[i_pos]>=K_NO_DATA_ISH&&i_pos<n_cols)
          {
          si_run_length++;
          i_pos++;
          if(si_run_length>=K_SHORT_MAX) // don't run off the end of a short for really big grids
            break;
          }// end while still data
        // we have found the end of the run.

        //Go back to start of run and read out ALL LAYERS
        // Write out run length first
        i_s_ptr=&(si_run_length);
        OutFileStream.write((char*)i_s_ptr,i_short_size);   // positive run length
        d_spos+=i_short_size;
        // Read out whoel run, all layer data for each point
        // Output data for each Grid cell for all layers
        // in order of the layers.. e.g. for 3 layers 1, 2,3, with run length 2
        // 2(runlength) 1 2 3 (firstcell) 1 2 3 (second cell)
        i_end_run=i_start_run+(int)si_run_length;
        for(i_col=i_start_run;i_col<i_end_run;i_col++)
          {
          f_value=f_data[i_col];
          f_ptr=&f_value;
          OutFileStream.write((char*)f_ptr,i_float_size);   // value for this layer
          d_spos+=i_float_size;
          }// end for i_col start to end run
        } // end DATA run
      } // end while
    // At the end of the row, output the offset from the start of the file
    d_ptr=&d_spos;
    IndexFileStream.write((char*)d_ptr,sizeof(double));
    }// end for i_row

  // Output an end of file token
  f_value=K_EOF;
  f_ptr=&f_value;
  OutFileStream.write((char*)f_ptr,i_float_size);   // eoftoken
  //Close files
  OutFileStream.close();
  IndexFileStream.close();
  InFileStream.close();

  // Free up arrays
  FreeFloat1DArray(f_data);

  return i_result_code;
}// end ConvertIntConditionGridToCompactBinaryFile

//************************************************************
// Convert Richness Grid  To Compact FLOAT Binary File
// takes a binary grid file in DIVA GIS .gri format, and squishes them a lot
// floats are floats
// the format is RUN LENGTH : datapoint 1 layer 1, datapoint1 layer 2 ..datapoint1 layer n:
// datapoint 2 layer 1, datapoint2 layer 2 ..datapoint2 layer n:  to datapoint n layer n
// NODATA runs are given as NEGATIVES, with no data afterwards.
// e.g. -5 5 34 23 24 25 21 23 25 12 25 25 -40 for a 2 layer grid.. is
// 5 nodata cells, 5 data cells, first cell, 34, 23 etc, followed by 40 nodatas
//****************************************************************
extern int ConvertRichnessGridToCompactFloatBinaryFile(const long int arg_i_num_rows,
														const long int arg_i_num_cols,
														string arg_s_condition_file,
														string arg_s_res_path)
{
int i_result_code;
long int n_rows,n_cols; //local copies of arguments
int n_layers;
long int i_row, i_col;  // row and column counters
float* f_data;
//float** f_cell_pptr;     // Pointer to float pointer to make array of a pointer for each layer
long int i_pos; // position we have reached
ifstream InFileStream;
ofstream OutFileStream;
ofstream IndexFileStream;
ofstream HeaderFileStream;
float* f_ptr;    // for reading data in
long int i_start_run, i_end_run;
short int si_run_length;
short* i_s_ptr;
float f_value;
string s_path,s_index_path;
int i_short_size;
double d_spos; // total offset of the start of each new line.
int i_length;
int* i_ptr;
int i_float_size;
int i_data_format;
double* d_ptr;
long int* l_ptr;

  i_result_code=1; // good unless we hear otherwise!
  i_short_size= sizeof(short);
  i_float_size=sizeof(float);
  n_rows= arg_i_num_rows;
  n_cols= arg_i_num_cols;

  // Allocate an array to hold the lines
  f_data= AllocateFloat1DArray(n_cols+1);

  // Open up condition grid file to compact
  InFileStream.open(arg_s_condition_file.c_str(), ios::in |ios::binary);
  // Open results file AS BINARY
  OutFileStream.open(arg_s_res_path.c_str(), ios::out |ios::binary);
  // build index file name ( to hold seekposes for the beginning of each line)
  s_path=arg_s_res_path;
  i_length=s_path.length();
  s_path.resize(i_length-4);
  s_index_path=s_path+".mci";
  IndexFileStream.open(s_index_path.c_str(), ios::out |ios::binary);

  // Quickly open and close and output header file details: cols. rows. layers.
  s_path=s_path+".mch";
  HeaderFileStream.open(s_path.c_str(), ios::out |ios::binary);
  l_ptr=&n_cols;
  HeaderFileStream.write((char*)l_ptr,sizeof(long int));
  l_ptr=&n_rows;
  HeaderFileStream.write((char*)l_ptr,sizeof(long int));
  i_ptr=&n_layers;
  n_layers=1;
  HeaderFileStream.write((char*)i_ptr,sizeof(int));
  i_data_format=K_FLOAT;
  i_ptr=&i_data_format;
  HeaderFileStream.write((char*)i_ptr,sizeof(int));
  HeaderFileStream.close();

  d_spos=0; // zero seekpos counter
  //Read whole lines into buffers
  for(i_row=0;i_row<n_rows;i_row++)
    {
    // read in whole line and process
    for(i_col=0;i_col<n_cols;i_col++)
      {
      // Read in value for each cell
      f_ptr=&(f_data[i_col]);
      //Read to f_data array
      InFileStream.read((char*)f_ptr,sizeof(float));
      } // end for i_col
    i_pos=0;
    i_col=0;
    //Work through the whole row
    while(i_pos<n_cols)
      {
      // split into nodata or data runs
      if(f_data[i_pos]<K_NO_DATA_ISH)
        {
        //NO_DATA
        //  i_col=i_pos;
        i_start_run=i_pos;
        si_run_length=0;
        while(f_data[i_pos]<K_NO_DATA_ISH&&i_pos<n_cols)
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
        d_spos+=i_short_size;
        } // end NODATA run
      else
        {
        //DATA
       // i_col=i_pos;
        i_start_run=i_pos;
        si_run_length=0;
        while(f_data[i_pos]>=K_NO_DATA_ISH&&i_pos<n_cols)
          {
          si_run_length++;
          i_pos++;
          if(si_run_length>=K_SHORT_MAX) // don't run off the end of a short for really big grids
            break;
          }// end while still data
        // we have found the end of the run.

        //Go back to start of run and read out ALL LAYERS
        // Write out run length first
        i_s_ptr=&(si_run_length);
        OutFileStream.write((char*)i_s_ptr,i_short_size);   // positive run length
        d_spos+=i_short_size;
        // Read out whoel run, all layer data for each point
        // Output data for each Grid cell for all layers
        // in order of the layers.. e.g. for 3 layers 1, 2,3, with run length 2
        // 2(runlength) 1 2 3 (firstcell) 1 2 3 (second cell)
        i_end_run=i_start_run+(int)si_run_length;
        for(i_col=i_start_run;i_col<i_end_run;i_col++)
          {
          f_value=f_data[i_col];
          f_ptr=&f_value;
          OutFileStream.write((char*)f_ptr,i_float_size);   // value for this layer
          d_spos+=i_float_size;
          }// end for i_col start to end run
        } // end DATA run
      } // end while
    // At the end of the row, output the offset from the start of the file
    d_ptr=&d_spos;
    IndexFileStream.write((char*)d_ptr,sizeof(double));
    }// end for i_row

  // Output an end of file token
  f_value=K_EOF;
  f_ptr=&f_value;
  OutFileStream.write((char*)f_ptr,i_float_size);   // eoftoken
  //Close files
  OutFileStream.close();
  IndexFileStream.close();
  InFileStream.close();

  // Free up arrays
  FreeFloat1DArray(f_data);

  return i_result_code;
}// end ConvertRichnessGridToCompactBinaryFile

extern int IsQuant(string arg_s_buffer)
{
string s_PredSplValTest;
string s_PredSplVal= "PredSplVal";

s_PredSplValTest=arg_s_buffer.substr(0,10);
if(s_PredSplValTest==s_PredSplVal)
  return (1);
else
  return (0);
}
extern int IsCoeff(string arg_s_buffer)
{
string s_PredCoefTest;
string s_PredCoef= "PredCoef";

s_PredCoefTest=arg_s_buffer.substr(0,8);
if(s_PredCoefTest==s_PredCoef)
  return (1);
else
  return (0);
}

extern int IsTRANSPREDS(string arg_s_buffer)
{
string s_TRANS= "[TRANSPREDS]";

if(arg_s_buffer==s_TRANS)
  return (1);
else
  return (0);
}
extern void ReadMCHFile(string arg_s_mch_path,
					   long int* n_cols_ptr,
					   long int* n_rows_ptr,
					   int* n_layers_ptr,
					   int* i_data_format_ptr)
{
int i_result_code;
long int n_rows,n_cols; //local copies of arguments
int n_layers,i_data_format;
long int i_row, i_col;  // row and column counters
ifstream HeaderFileStream;
string s_path;
int* i_ptr;
long int* l_ptr;

 // open the little binary file and get the info out
 // do this bit with our local pointers..probably not necessary!
 s_path=arg_s_mch_path;
  HeaderFileStream.open(s_path.c_str(), ios::in |ios::binary);
  l_ptr=&n_cols;
  HeaderFileStream.read((char*)l_ptr,sizeof(long int));
  l_ptr=&n_rows;
  HeaderFileStream.read((char*)l_ptr,sizeof(long int));
  i_ptr=&n_layers;
  HeaderFileStream.read((char*)i_ptr,sizeof(int));
  i_ptr=&i_data_format;
  HeaderFileStream.read((char*)i_ptr,sizeof(int));
  HeaderFileStream.close();
 // copy into supplied pointers
 *n_cols_ptr=n_cols;
 *n_rows_ptr=n_rows;
 *n_layers_ptr=n_layers;
 *i_data_format_ptr=i_data_format;
}// end reade mch file
//---------------------------------------------------------------------------
 extern double GetQuant(string arg_s_buffer)
 {
//string s_buffer;
string s_PredSplValTest;
string s_PredCoefTest;
string s_PredSplVal= "PredSplVal";
string s_PredCoef= "PredCoef";
string s_left;
stringstream ss;
double d_result;
const char* c_string;
int i_char;
char c_equals= '=';
int i_layer=0;
int i_knot=0;

	//Found a predictor spline set
	//Extract the first spline value
	s_left= arg_s_buffer.substr(10,arg_s_buffer.length());
	c_string=s_left.c_str();
	i_char=0;
	while ( c_string[i_char]!=c_equals)
	  i_char++;
	s_left=s_left.substr(i_char+1,s_left.length() );
	ss<<s_left;
	ss>> d_result;
	return d_result;
 }
 //---------------------------------------------------------------------------
extern double GetCoeff(string arg_s_buffer)
{
//string s_buffer;
string s_PredSplValTest;
string s_PredCoefTest;
string s_PredSplVal= "PredSplVal";
string s_PredCoef= "PredCoef";
string s_left;
double d_result;
const char* c_string;
stringstream ss;
int i_char;
char c_equals= '=';
int i_layer=0;
int i_knot=0;

	//Found a predictor spline set
	//Extract the first spline value
	s_left= arg_s_buffer.substr(8,arg_s_buffer.length());
	c_string=s_left.c_str();
	i_char=0;
	while ( c_string[i_char]!=c_equals)
	  i_char++;
	s_left=s_left.substr(i_char+1,s_left.length() );
   ss<<s_left;
	ss>> d_result;
	return d_result;
}
//---------------------------------------------------------------------------

