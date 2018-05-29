//*****************************************************************************
// Dynamic Array . cpp
// Suite of functions for the allocation and deletion of dynamic arrays, 1& 2D
// Written Tom Harwood May 2010
//*****************************************************************************

#include "DynamicArray.h"

//*****************************************************************************
//*****************************************************************************
// ***** Allocate 1D Array ******** TEMPLATE VERSION **************************
// Creates a 1D array of the specified dimension
// TEMPLATE version of the function.. creates according to the template
// specification in the calling function
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use Free1DArray(TEMPLATE**)
// Arguments:
// const int arg_i_n_cols    number of columns
// Returns:
// Pointer to the TEMPLATE type pointing to the newly allocated 1D array
// ..which can be referenced as result[].
// usage example:
//      double* d_fish;
//      d_fish= Allocate1DArray<double>(100);
//      d_fish[21]=51.563;
//      Free1DArray(d_fish);
//*****************************************************************************
template < typename T >
T* Allocate1DArray(const int arg_i_n_cols)
{
  T* T_new_ptr = new T[arg_i_n_cols];
  return T_new_ptr;
} // End func Allocate1DArray
//*****************************************************************************
// ***** Free 1D Array ********
// Deletes the suuplied 1D array of the supplied TYPE
// Memory: Deletes argument
// Arguments:
// T* arg_T_Array    the array to delete  // any type
// Returns: none
//
//*****************************************************************************
template < typename T >
void Free1DArray(T* arg_T_Array)
{
  delete [] arg_T_Array;
}  // end func delete 1D array
//*****************************************************************************

//*****************************************************************************
//*****************************************************************************
// ***** Allocate Float 1D Array
// Creates a 1D array of the specified dimension
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use FreeFloat1DArray(float**)
// Arguments:
// const int arg_i_n_cols    number of columns
// Returns:
// Pointer  to the float* pointing to the newly allocated 1D array
//*****************************************************************************
float* AllocateFloat1DArray(const int arg_i_n_cols)
{
  float* f_new_ptr = new float[arg_i_n_cols];
  return f_new_ptr;
} // End func AllocateFloat1DArray
//*****************************************************************************
// ***** Free INT 1D Array ********
// Deletes the suuplied 1D array of the supplied TYPE
// Memory: Deletes argument
// Arguments:
// T* arg_T_Array    the array to delete  // any type
// Returns: none
//
//*****************************************************************************
void FreeFloat1DArray(float* arg_f_Array)
{
  delete [] arg_f_Array;
}  // end func delete float 1D array
//*****************************************************************************

//*****************************************************************************
//*****************************************************************************
// ***** Allocate Double 1D Array
// Creates a 1D array of the specified dimension
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use FreeDouble1DArray(double*)
// Arguments:
// const int arg_i_n_cols    number of columns
// Returns:
// Pointer  to the float* pointing to the newly allocated 1D array
//*****************************************************************************
double* AllocateDouble1DArray(const int arg_i_n_cols)
{
  double* d_new_ptr = new double[arg_i_n_cols];
  return d_new_ptr;
} // End func AllocateFloat1DArray
//*****************************************************************************
// ***** Free double 1D Array ********
// Deletes the suuplied 1D array of the supplied TYPE
// Memory: Deletes argument
// Arguments:
// T* arg_T_Array    the array to delete  // any type
// Returns: none
//
//*****************************************************************************
void FreeDouble1DArray(double* arg_d_Array)
{
  delete [] arg_d_Array;
}  // end func delete float 1D array
//*****************************************************************************

//*****************************************************************************
//*****************************************************************************
// ***** Allocate Int 1D Array
// Creates a 1D array of the specified dimension
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use FreeInt1DArray(int**)
// Arguments:
// const int arg_i_n_cols    number of columns
// Returns:
// Pointer  to the TEMPLATE type pointing to the newly allocated 1D array
//*****************************************************************************
int* AllocateInt1DArray(const int arg_i_n_cols)
{
  // Allocate memory for array
  int* i_new_ptr = new int[arg_i_n_cols];
  return i_new_ptr;
} // End func AllocateInt1DArray
//*****************************************************************************
// ***** Free INT 1D Array ********
// Deletes the suuplied 1D array of the supplied TYPE
// Memory: Deletes argument
// Arguments:
// T* arg_T_Array    the array to delete  // any type
// Returns: none
//
//*****************************************************************************
void FreeInt1DArray(int* arg_i_Array)
{
  // Free up array
  delete [] arg_i_Array;
}  // end func delete INT 1D array
//*****************************************************************************

//*****************************************************************************
//*****************************************************************************
// ***** Allocate Long 1D Array
// Creates a 1D array of the specified dimension
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use FreeInt1DArray(int**)
// Arguments:
// const long int arg_i_n_cols    number of columns
// Returns:
// Pointer  to the TEMPLATE type pointing to the newly allocated 1D array
//*****************************************************************************
long int* AllocateLong1DArray(const long int arg_i_n_cols)
{
  // Allocate memory for array
  long int* i_new_ptr = new long int[arg_i_n_cols];
  return i_new_ptr;
} // End func AllocateLong1DArray
//*****************************************************************************
// ***** Free Long 1D Array ********
// Deletes the suuplied 1D array of the supplied TYPE
// Memory: Deletes argument
// Arguments:
// T* arg_T_Array    the array to delete  // any type
// Returns: none
//
//*****************************************************************************
void FreeLong1DArray(long int* arg_i_Array)
{
  // Free up array
  delete [] arg_i_Array;
}  // end func delete LONG INT 1D array
//*****************************************************************************

//*****************************************************************************
//*****************************************************************************
// ***** Allocate Short 1D Array
// Creates a 1D array of the specified dimension
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use FreeShort1DArray(short**)
// Arguments:
// const int arg_i_n_cols    number of columns
// Returns:
// Poshorter  to the TEMPLATE type poshorting to the newly allocated 1D array
//*****************************************************************************
short* AllocateShort1DArray(const int arg_i_n_cols)
{
  // Allocate memory for array
  short* si_new_ptr = new short[arg_i_n_cols];
  return si_new_ptr;
} // End func AllocateShort1DArray
//*****************************************************************************
// ***** Free SHORT 1D Array ********
// Deletes the suuplied 1D array of the supplied TYPE
// Memory: Deletes argument
// Arguments:
// T* arg_T_Array    the array to delete  // any type
// Returns: none
//
//*****************************************************************************
void FreeShort1DArray(short* arg_si_Array)
{
  // Free up array
  delete [] arg_si_Array;
}  // end func delete SHORT 1D array
//**********************

//*****************************************************************************
//*****************************************************************************
// ***** Allocate Unsigned Short 1D Array
// Creates a 1D array of the specified dimension
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use FreeUnsignedShort1DArray(short**)
// Arguments:
// const int arg_i_n_cols    number of columns
// Returns:
// Poshorter  to the TEMPLATE type poshorting to the newly allocated 1D array
//*****************************************************************************
unsigned short* AllocateUnsignedShort1DArray(const int arg_i_n_cols)
{
  // Allocate memory for array
  unsigned short* si_new_ptr = new unsigned short[arg_i_n_cols];
  return si_new_ptr;
} // End func AllocateUnsigendShort1DArray
//*****************************************************************************
// ***** Free UNSIGNED SHORT 1D Array ********
// Deletes the suuplied 1D array of the supplied TYPE
// Memory: Deletes argument
// Arguments:
// T* arg_T_Array    the array to delete  // any type
// Returns: none
//
//*****************************************************************************
void FreeUnsignedShort1DArray(unsigned short* arg_si_Array)
{
  // Free up array
  delete [] arg_si_Array;
}  // end func delete SHORT 1D array
//**********************
//*****************************************************************************
//*****************************************************************************
// ***** Allocate 2D Array ******** TEMPLATE VERSION **************************
// Creates a 2D array of the specified dimensions
// TEMPLATE version of the function.. creates according to the template
// specification in the calling function
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use Free2DArray(TEMPLATE**)
// Arguments:
// const int arg_i_n_rows    number of rows
// const int arg_i_n_cols    number of columns
// Returns:
// Pointer to Pointer to the TEMPLATE type pointing to the newly allocated 2D array
// ..which can be referenced as result[][].
// usage example:
//      double** d_fish;
//      d_fish= Allocate2DArray<double>(100,177);
//      d_fish[21][54]=21.563;
//      Free2DArray(d_fish);
//*****************************************************************************
template < typename T >
T **Allocate2DArray(const int arg_i_n_rows,
                    const int arg_i_n_cols)
{
int i_row; // row counter
  // Allocate memory for array of elements of column ..
  // i.e. list of pointers to first element of each row
  T** T_new_pptr = new T*[arg_i_n_rows];
  // Allocate memory for array of elements of all the rows
  // Do this as one big block to speed things up
  T* T_row_ptr = new T[arg_i_n_rows*arg_i_n_cols];
  // Now point the pointers in the right place
  for(i_row=0; i_row<arg_i_n_rows; i_row++)
    {
    *(T_new_pptr+i_row)=T_row_ptr; // Could also be written as T_new_pptr[i_row]=
    T_row_ptr+=arg_i_n_cols;       // skip forward to the beginning of the next row
    }
  return T_new_pptr;
} // End func Allocate2DArray

//*****************************************************************************
// ***** Free 2D Array ********
// Deletes the suuplied 2D array of the supplied TYPE
// Memory: Deletes argument
// Arguments:
// T** arg_T_Array    the array to delete  // any type
// Returns: none
//
//*****************************************************************************
template < typename T >
void Free2DArray(T** arg_T_Array)
{
  // Free up both levels of the allocated array
  delete [] *arg_T_Array;
  delete [] arg_T_Array;
}  // end func delete 2D array
//*****************************************************************************



//*****************************************************************************
//*****************************************************************************
// ***** Allocate Float 2D Array ********
// Creates a 2D array of the specified dimensions
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use FreeFloat2DArray(float**)
// Arguments:
// const int arg_i_n_rows    number of rows
// const int arg_i_n_cols    number of columns
// Returns:
// Pointer to Pointer to Float pointing to the newly allocated 2D array
// ..which can be referenced as result[][].
// usage example:
//      float** f_fish;
//      f_fish= AllocateFloat2DArray(10,10);
//      f_fish[2][5]=21.563;
//      FreeFloat2DArray(f_fish);
//*****************************************************************************
float **AllocateFloat2DArray( const int arg_i_n_rows,
                              const int arg_i_n_cols)
{
int i_row;              // row counter
float** f_new_pptr;     // Pointer to pointer for the whole 2D buffer
float* f_row_ptr;       // Pointer to start of each row in the array

  // Allocate memory for array of elements of column ..
  // i.e. list of pointers to first element of each row
  f_new_pptr = new float*[arg_i_n_rows];
  // Allocate memory for array of elements of all the rows
  // Do this as one big block to speed things up
  f_row_ptr = new float [arg_i_n_rows* arg_i_n_cols];

  // Now point the pointers in the right place
  for(i_row=0; i_row<arg_i_n_rows; i_row++)
    {
    *(f_new_pptr+i_row) = f_row_ptr;
    f_row_ptr += arg_i_n_cols;
    }
  return f_new_pptr;
}  // end func Allocate Float 2D Array

//*****************************************************************************
// ***** Free Float 2D Array ********
// Deletes the suuplied 2D array
// Memory: Deletes argument
// Arguments:
// float** arg_f_Array    the array to delete
// Returns: none
//
//*****************************************************************************
void FreeFloat2DArray(float** arg_f_Array)
{
  // Free up both levels of the allocated array
  delete [] *arg_f_Array;
  delete [] arg_f_Array;
} // end func delete Float 2D array



//*****************************************************************************
//*****************************************************************************
// ***** Allocate Int 2D Array ********
// Creates a 2D array integers of the specified dimensions
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use FreeInt2DArray(int**)
// Arguments:
// const int arg_i_n_rows    number of rows
// const int arg_i_n_cols    number of columns
// Returns:
// Pointer to Pointer to Int pointing to the newly allocated 2D array
// ..which can be referenced as result[][].
// usage example:
//      int** i_fish;
//      i_fish= AllocateInt2DArray(60,300);
//      i_fish[27][35]=26.563;
//      FreeInt2DArray(i_fish);
//*****************************************************************************
int** AllocateInt2DArray( const int arg_i_n_rows,
                            const int arg_i_n_cols)
{
int i_row;              // row counter
int** i_new_pptr;     // Pointer to pointer for the whole 2D buffer
int* i_row_ptr;       // Pointer to start of each row in the array

  // Allocate memory for array of elements of column ..
  // i.e. list of pointers to first element of each row
  i_new_pptr = new int*[arg_i_n_rows];
  // Allocate memory for array of elements of all the rows
  // Do this as one big block to speed things up
  i_row_ptr = new int[arg_i_n_rows*arg_i_n_cols];

  // Now point the pointers in the right place
  for(i_row=0; i_row<arg_i_n_rows; i_row++)
    {
    *(i_new_pptr+i_row) = i_row_ptr;
    i_row_ptr += arg_i_n_cols;
    }
  return i_new_pptr;
}  // end func Allocate Int 2D Array
//*****************************************************************************
// ***** Free Int 2D Array ********
// Deletes the suuplied 2D array
// Memory: Deletes argument
// Arguments:
// int** arg_i_Array    the array to delete
// Returns: none
//
//*****************************************************************************
void FreeInt2DArray(int** arg_i_Array)
{
  // Free up both levels of the allocated array
  delete [] *arg_i_Array;
  delete [] arg_i_Array;
} // end func delete Int 2D array

//*****************************************************************************
//*****************************************************************************
// ***** Allocate Short 2D Array ********
// Creates a 2D array short integers of the specified dimensions
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use FreeShort2DArray(int**)
// Arguments:
// const int arg_i_n_rows    number of rows
// const int arg_i_n_cols    number of columns
// Returns:
// Pointer to Pointer to Int pointing to the newly allocated 2D array
// ..which can be referenced as result[][].
// usage example:
//      short int** i_fish;
//      i_fish= AllocateShort2DArray(60,300);
//      i_fish[27][35]=26.563;
//      FreeShort2DArray(i_fish);
//*****************************************************************************
short int** AllocateShort2DArray( const int arg_i_n_rows,
                                  const int arg_i_n_cols)
{
int i_row;              // row counter
short int** i_new_pptr;     // Pointer to pointer for the whole 2D buffer
short int* i_row_ptr;       // Pointer to start of each row in the array

  // Allocate memory for array of elements of column ..
  // i.e. list of pointers to first element of each row
  i_new_pptr = new short int*[arg_i_n_rows];
  // Allocate memory for array of elements of all the rows
  // Do this as one big block to speed things up
  i_row_ptr = new short int[arg_i_n_rows*arg_i_n_cols];

  // Now point the pointers in the right place
  for(i_row=0; i_row<arg_i_n_rows; i_row++)
    {
    *(i_new_pptr+i_row) = i_row_ptr;
    i_row_ptr += arg_i_n_cols;
    }
  return i_new_pptr;
}  // end func Allocate Int 2D Array
//*****************************************************************************
// ***** Free Short 2D Array ********
// Deletes the suuplied 2D array
// Memory: Deletes argument
// Arguments:
// int** arg_i_Array    the array to delete
// Returns: none
//
//*****************************************************************************
void FreeShort2DArray(short int** arg_i_Array)
{
  // Free up both levels of the allocated array
  delete [] *arg_i_Array;
  delete [] arg_i_Array;
} // end func delete Int 2D array


//---------------------------------------------------------------------------

 