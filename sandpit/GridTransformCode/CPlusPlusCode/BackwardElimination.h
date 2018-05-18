//
// BackwardElimination.h
//
#ifndef __BACKWARDELIMINATION_H__
#define __BACKWARDELIMINATION_H__

//
// Run the core GDM for a Backward Elimination
//
bool DoBaseGDM(char *pParams, bool BatchRun, bool DeleteBinaryFile, 
	           double **pResponse, double **pWeights, 
			   int *NumRows, int *NumCols,
			   FPTR fptr);


//
// Do an incremental GDM for the Backward Elimination Process
//
bool DoIncrementalGDM(char *pParams, bool BatchRun, bool DeleteBinaryFile, 
	                  char *lpTmpFile, char *myRunName,
	                  double *pResponse, double *pWeights, 
			          int NumRows, int NumCols, int Index,
			          FPTR fptr, UPTR uptr);

#endif // __BACKWARDELIMINATION_H__