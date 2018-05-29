//---------------------------------------------------------------------------
#ifndef CompactorFunctionsH
#define CompactorFunctionsH
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
using namespace std;
// Float versions
extern int ConvertGridStackToCompactFloatBinaryFile(const long int arg_i_num_rows,
                                                    const long int arg_i_num_cols,
                                                    const int arg_i_num_layers,
                                                    string arg_s_layer_files[],
													string arg_s_res_path);
extern int ConvertGridStackToCompactFloatBinaryFileWithZeroGrids(const long int arg_i_num_rows,
                                                    const long int arg_i_num_cols,
                                                    const int arg_i_num_layers,
                                                    string arg_s_layer_files[],
													string arg_s_res_path);

extern int ConvertConditionGridToCompactFloatBinaryFile(const long int arg_i_num_rows,
                                                        const long int arg_i_num_cols,
                                                        string arg_s_condition_file,
                                                        string arg_s_res_path);
extern int ConvertIntConditionGridToCompactFloatBinaryFile(const long int arg_i_num_rows,
                                                        const long int arg_i_num_cols,
                                                        string arg_s_condition_file,
                                                        string arg_s_res_path);
extern int ConvertRichnessGridToCompactFloatBinaryFile(const long int arg_i_num_rows,
                                                        const long int arg_i_num_cols,
                                                        string arg_s_condition_file,
                                                        string arg_s_res_path);

extern int IsQuant(string arg_s_buffer);
extern int IsCoeff(string arg_s_buffer);
extern int IsTRANSPREDS(string arg_s_buffer);
extern double GetQuant(string arg_s_buffer);
extern double GetCoeff(string arg_s_buffer);
//---------------------------------------------------------------------------
#endif
