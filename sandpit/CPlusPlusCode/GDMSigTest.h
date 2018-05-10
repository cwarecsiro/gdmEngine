#pragma once
//
// GDMSigTest.h
//
#ifndef __GDMSIGTEST_H__
#define __GDMSIGTEST_H__

double *GetCompositeTableIntoMemorySigTest(char *pDataPath, GDM_INT *pTableRows, GDM_INT *pTableCols);
int CountRows64SigTest(char *pPath);
int CountColumns64SigTest(char *pPath);

//
// Run the core GDM for a Significance Test
//
bool DoBaseSigGDM(char *pParams, bool DeleteBinaryFile,
	              double **pResponse, double **pWeights,
	              int *NumRows, int *NumCols,
	              FPTR fptr);


//
// Do a single GDM predictor significance test
//
double DoPredictorSigGDM(char *pParams, char *lpTmpFile,
	double *pResponse, double *pWeights, double baselineDevExpl,
	int NumRows, int NumCols,
	int Index, int NumIterations,
	int ColumnJump, int nPredSplines,
	double *pMin, double *pMax,
	FPTR fptr, UPTR uptr);


#endif // __GDMSIGTEST_H__


