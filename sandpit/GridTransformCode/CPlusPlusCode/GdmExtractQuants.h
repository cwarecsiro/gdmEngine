//
// GdmExtractQuants.h
//
#ifndef __GDMEXTRACTQUANTS_H__
#define __GDMEXTRACTQUANTS_H__


//
// Return the data table after extracting the quantiles and updating the parameter file
//
double *GetTableDataPostExtractAndUpdateQuantiles( char *pParams, 
	                                               GDM_INT *pRows, 
												   GDM_INT *pCols, 
												   bool BatchRun, 
												   FPTR fptr);


//
// Extract lower triangular data from lookup table into pPredDist and return number of items
//
int GetLowerTriFromLookup(char *pParams, int Index, double *pPredDist);


//
// sort a vector of doubles in ascending order
//
void SortVector(double *pDistVector, int NumItems);


//
// comparison routine for quicksort of doubles
//
int CompareDoubles( const void *, const void * );


//
// return the n% quantile of a sorted vector of doubles
//
double GetQuantAsPercentile( double *pVector, int nLen, int Percentile );


#endif // __GDMEXTRACTQUANTS_H__