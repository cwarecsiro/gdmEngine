//
// GDMDensity.h
//
#ifndef __GDMDENSITY_H__
#define __GDMDENSITY_H__

#include "ANN.h"

//
//
// Count the Valid data cells in a Binary Export Grid 
//
int GetValidDataCellCount(char *GridPath, FPTR fptr);


//
// Populate ANN columns from transform grid
//
void DoESRIBinaryANNColumn( ANNpointArray data_pnts, int Index, char *GridPath, double **ppPoints, int nRecords );



//
// Get condition data from a GDM transform grid
//
void GetESRIBinaryConditionVals( char *lpDomainPath, float *fCondition, double **ppPoints, int nRecords );


//
// Creates a sample point file from Floating Point Domain Grid for ther Density Calculations
// 
bool CreateSamplePointMeshForDensity( char *lpParams, 
	                                  char *lpSamplePointPath, 
	                                  char *lpDomainPath, 
									  int nSamples, bool DoBatch, FPTR fptr );


//
// Populate DATA columns from transform grid
//
void DoESRIBinaryDataColumn( double **data_pnts, int Index, char *GridPath, double **ppPoints, int nRecords );


//
// get the environmental distance (ED)
//
double CalculateED(float **ppGridData, int ColIndex, double *pSample, int nTranGrids);



#endif // __GDMDENSITY_H__