//
// DistanceCalcs.h
//
#ifndef __DISTANCECALCS_H__
#define __DISTANCECALCS_H__


#define WT_NO_WEIGHTS 0
#define WT_WEIGHT_BY_SUM 1
#define WT_WEIGHT_BY_SQRTSUM 2


//
// Calculate the Bray-Curtis SIMILARITY Metric for a pair of vectors with presence absence values
//
double BrayCurtisSimilarityPresenceAbsence(int *p0, int *p1, int NumSpecies );

//
// Calculate the Bray-Curtis DISIMILARITY Metric for a pair of vectors with presence absence values
//
double BrayCurtisDisimilarityPresenceAbsence(int *p0, int *p1, int NumSpecies );

//
// Calculate the Bray-Curtis SIMILARITY Metric for a pair of vectors with abundance values
//
double BrayCurtisSimilarityAbundance(int *p0, int *p1, int NumSpecies );

//
// Calculate the Bray-Curtis DISIMILARITY Metric for a pair of vectors with abundance values
//
double BrayCurtisDissimilarityAbundance(int *p0, int *p1, int NumSpecies );

//
// Calculate the Bray-Curtis SIMILARITY Metric for a pair of vectors with abundance values as doubles
//
double BrayCurtisSimilarityAbundanceDouble(double *p0, double *p1, int NumSpecies );

//
// Calculate the Bray-Curtis DISIMILARITY Metric for a pair of vectors with abundance values as doubles
//
double BrayCurtisDissimilarityAbundanceDouble(double *p0, double *p1, int NumSpecies );

//
// Returns the Square Root of the sum of the two integer vectors
//
double GetWeightsSqrt(int *p0, int *p1, int NumSpecies);

//
// Returns the sum of the two integer vectors
//
double GetWeightsDefault(int *p0, int *p1, int NumSpecies);

//
// Returns the Square Root of the sum of the two double vectors
//
double GetWeightsDoubleSqrt(double *p0, double *p1, int NumSpecies);

//
// Returns the sum of the two double vectors
//
double GetWeightsDoubleDefault(double *p0, double *p1, int NumSpecies);

//
// Return the number of present records at a site in a SxS Table
//
int GetCountForSite(int *p, int NumSpecies);


#endif