//
// GdmGuiLib.h
//
#ifndef __GDMGUILIB_H__
#define __GDMGUILIB_H__

#include "stdafx.h"
#include "myCallback.h"


//
// Get a memory image of the composite data file
//
double *GetCompositeTableIntoMemory(char *pDataPath, GDM_INT *pTableRows, GDM_INT *pTableCols, bool BatchRun, FPTR fptr);


//
// Create a binary file image of the splined predictor data to pass to the matrix regression
// There will always be the same number of rows in the table data as the predictor matrix
// but the number of matrix cols may be less if not all predictors are 'In-use'.
//
bool CreatePredictorBinary(char *lpBinFile, char *pParams,
	                       double *pTableData, GDM_INT nTableRows, GDM_INT nTableCols, 
						   int *pSplineCounts, double *pQuantiles,
						   GDM_INT *pMatrixCols, bool BatchRun, FPTR fptr);


//
// Convert the binary matrix created in CreatePredictorBinary to a CSV file
//
void PrintBinary2CSV(char *lpBinFile, char *pParams, int *
	                 pSplineCounts, GDM_INT nTableRows, GDM_INT nMatrixCols, FPTR fptr);

//
// Dump the contents of the table memory back to file for debugging
//
void PrintTableData(char *pParams, double *pData, int nRows, int nCols, bool BatchRun);


//
// Dump the contents of the quantile vector memory back to file for debugging
//
void PrintQuantileVector(char *pParams,  double *pQuantiles, int nNumQuants, bool BatchRun);


//
// Dump the GDM results to text file for debugging
//
void PrintGDMResults(char *pParams, double dDevianceExplained, int nSplines, double *pCoeffs, bool BatchRun);



//
// Write the GDM results back to the parameter file
//
void UpdateParameterFile(char *pParams, 
		                 double dNullDeviance, double dGDMDeviance, double dDevianceExplained, 
						 double dIntercept, double *pCoefficients, 
						 int nTotalSplines, int *pSplineCounts, int nRows, FPTR fptr);

//
// Sum the GDM coefficients
//
double CalcCoefficientSum(double *pCoefficients, int nTotalSplines);


//
// Sum the number of active environmental predictors
//
int CalcNumberOfActivePredictors(char *pParams, double *pCoefficients, int *pSplineCounts );


//
// Get the number of active predictors in the current GDM model
//
int GetNumActivePredictors(char *pParams);


//
// Sum the number of coefficients > 0
//
int CalcNumberofCoeffGreaterThanZero(double *pCoefficients, int nTotalSplines);


//
// Get the spline counts for all the predictors including geographic distance 
// as a vector from the parameter file (ignore whether they are "in-use" or not)
//
int *GetSplineCountsFromParams(char *pParams, int *pItems);


//
// Get the quantiles for all the predictors including geographic distance 
// as a vector from the parameter file (ignore whether they are "in-use" or not)
//
double *GetQuantilesFromParams(char *pParams, int *pItems);


//
// Write the observed and predicted data to comma delimited textfile
//
void CreateObservedVsPredicted(char *pParams, double *pResponse, double *pPredData, double *pCoefficients, int nRows, int nCols, FPTR fptr);


//
// Get the spline counts for all the predictors including geographic distance
// as a vector from the parameter file if they are "in-use"
//
int *GetInUseSplineCountVectorFromParams(char *pParams, int *pItems);


//
// Get the quantiles for all the predictors including geographic distance as a vector if they are "in-use"
//
double *GetInUseQuantilesFromParams(char *pParams, int *pItems);


//
// Extract the coefficients from the "in-Use" predictors from a GDM Model param file
//
double *GetCoefficientsFromGDMParams(char *pModelParams, int *numCoefficients );


//
// Write the observed and predicted data to comma delimited textfile and to a pair of binary files for plotting
//
void CreatePredictionTable(char *pOutTablePath, double *pResponse, 
	                       double *pX0, double *pY0, double *pX1, double *pY1,
						   double *pPredData, double *pCoefficients, 
						   int nRows, int nCols, FPTR fptr);


#endif // __GDMGUILIB_H__