//
// NNLS_Double.h
//
#ifndef __NNLS_DOUBLE_H__
#define __NNLS_DOUBLE_H__

#include "myCallback.h"


#if defined _M_X64


double *WeightedNNLSRegression( char *lpTmpFile, 
							    double *pEnvDataMatrix, GDM_INT nRows, GDM_INT nCols, 
								double *pRespVector, double *pDeviance, double *pWeights, FPTR fptr );

double *WeightedNNLSResidualRegression( char *lpTmpFile, 
	                                    double *pResidualVector,
							            double *pEnvDataMatrix, GDM_INT nRows, GDM_INT nCols, 
								        double *pRespVector, double *pDeviance, double *pWeights, FPTR fptr );

double CalcGDMDevianceDouble( double *pY, double *pU, double *pW, GDM_INT nLen );

double GetWeightedNULLDeviance( GDM_INT nRows, double *pRespVector, double *pWeights );

double *CopyEnvMatrixDouble( double *pMatrix, GDM_INT nRows, GDM_INT nCols );

double *nnlsFITDouble( double *pEnvDataMatrix, GDM_INT nRows, GDM_INT nCols, double *pRespVector, double *pWeights );


#elif defined _WIN32


double *WeightedNNLSRegression( char *lpTmpFile, 
							    double *pEnvDataMatrix, int nRows, int nCols, 
								double *pRespVector, double *pDeviance, double *pWeights, FPTR fptr );

double *WeightedNNLSResidualRegression( char *lpTmpFile, 
	                                    double *pResidualVector,
							            double *pEnvDataMatrix, int nRows, int nCols, 
								        double *pRespVector, double *pDeviance, double *pWeights, FPTR fptr );


double CalcGDMDevianceDouble( double *pY, double *pU, double *pW, int nLen );

double GetWeightedNULLDeviance( int nRows, double *pRespVector, double *pWeights );

double GetTestWeightedNULLDeviance( double *pMatrix, int nCols, int nRows, double *pRespVector, double *pWeights );

double *CopyEnvMatrixDouble( double *pMatrix, int nRows, int nCols );

double *nnlsFITDouble( double *pEnvDataMatrix, int nRows, int nCols, double *pRespVector, double *pWeights );

#endif

#endif  // __NNLS_DOUBLE_H__