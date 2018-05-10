//
// GdmUnconstrainedProbGrids.h
//
#ifndef __GDMUNCONSTRAINEDPROBGRIDS_H__
#define __GDMUNCONSTRAINEDPROBGRIDS_H__

#include "ANN.h"


int GetMaxTrainingClass(char *TrainingDataPath, FPTR fptr);


//
// Calculate the Sum Of Squares for an ANN dataset given 
// a particular number of nearest neighbours (nJ) to find
//
double GetSumOfSquaresANN( ANNkd_tree *theTree, 
	                       int nDataRows, int nDataCols, 
						   int nJ, double dR, 
						   ANNidxArray nn_idx,
						   ANNdistArray dist,
						   double *explu,
						   ANNpoint query_pnt,
						   ANNpointArray data_pnts,
						   int NumClasses, int *pClasses, double *pOutVec,
						   FPTR fptr);



//
// Determine the number of data rows from a comma delimited text file with header
// and determine the number of predictor columns after skipping Class,X and Y.
//
bool GetANNTableMetrics(char *TablePath, int *pRows, int *pCols);



//
// Extract the predictor columns from a comma delimited training data table 
// after skipping the header and after skipping the first three columns (class,X,Y)
//
ANNpointArray GetDataAsANNTable(char *TablePath, int nRows, int nCols, FPTR fptr);



//
// Extract the class field data [0] from a comma delimited training data table 
// after skipping the header and also determining the Maximum Class Index
//
int *GetClassesFromANN(char *TablePath, int nRows, int *MaxClassIndex, FPTR fptr);



//
// Extract a vector of class indices representing classes to skip in the probability grid creation
//
int *GetEmptyClasses(int *pNumEmptyClasses, char *pTrainingFile, int nMaxClass);



//
// Returns true if veg class in empty in the training data table
//
bool EmptyClass(int Class, int *pClassVector, int nSize);



//
// Returns true if we can extract a valid tansform grid path from a GDM parameter file
//
bool GetTransformGridAsDomain( char *pParams, char *pDomainPath );



#endif __GDMUNCONSTRAINEDPROBGRIDS_H__