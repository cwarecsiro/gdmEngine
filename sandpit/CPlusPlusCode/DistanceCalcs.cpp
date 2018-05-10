//
// DistanceCalcs.cpp
//
#include "stdafx.h"
#include "DistanceCalcs.h"


//
// Calculate the Bray-Curtis SIMILARITY Metric for a pair of vectors with presence absence values
//
double BrayCurtisSimilarityPresenceAbsence(int *p0, int *p1, int NumSpecies )
{
	double fA = 0.0;
	double fB = 0.0;
	double fC = 0.0;

	for ( int i=0; i<NumSpecies; i++ )
	{
		if (( p0[i] == 1 ) && (p1[i] == 1))
			fA += 1.0;

		else if (( p0[i] == 1 ) && (p1[i] == 0))
			fB += 1.0;

		else if (( p0[i] == 0 ) && (p1[i] == 1))
			fC += 1.0;
	}
	
	double dDenominator = (2.0 * fA) + fB + fC;
	if ( dDenominator == 0.0 )
		return(0.0);
	else
		return((2.0 * fA) / dDenominator);
}



//
// Calculate the Bray-Curtis DISIMILARITY Metric for a pair of vectors with presence absence values
//
double BrayCurtisDisimilarityPresenceAbsence(int *p0, int *p1, int NumSpecies )
{
	return(1.0 - BrayCurtisSimilarityPresenceAbsence(p0, p1, NumSpecies));
}



//
// Calculate the Bray-Curtis SIMILARITY Metric for a pair of vectors with abundance values
//
double BrayCurtisSimilarityAbundance(int *p0, int *p1, int NumSpecies )
{
	return(1.0 - BrayCurtisDissimilarityAbundance(p0, p1, NumSpecies));
}



//
// Calculate the Bray-Curtis DISIMILARITY Metric for a pair of vectors with abundance values
//
double BrayCurtisDissimilarityAbundance(int *p0, int *p1, int NumSpecies )
{
	double fA = 0.0;
	double fB = 0.0;

	for ( int i=0; i<NumSpecies; i++ )
	{
		fA += double(abs( p0[i] - p1[i]));
		fB += double(p0[i] + p1[i]);
	}

	if (fB == 0.0)
		return(0.0);
	else
		return(fA / fB);
}



//
// Calculate the Bray-Curtis SIMILARITY Metric for a pair of vectors with abundance values as doubles
//
double BrayCurtisSimilarityAbundanceDouble(double *p0, double *p1, int NumSpecies )
{
	return(1.0 - BrayCurtisDissimilarityAbundanceDouble(p0, p1, NumSpecies));
}



//
// Calculate the Bray-Curtis DISIMILARITY Metric for a pair of vectors with abundance values as doubles
//
double BrayCurtisDissimilarityAbundanceDouble(double *p0, double *p1, int NumSpecies )
{
	double fA = 0.0;
	double fB = 0.0;

	for ( int i=0; i<NumSpecies; i++ )
	{
		fA += fabs( p0[i] - p1[i]);
		fB += p0[i] + p1[i];
	}

	if (fB == 0.0)
		return(0.0);
	else
		return(fA / fB);
}



//
// Returns the Square Root of the sum of the two integer vectors
//
double GetWeightsSqrt(int *p0, int *p1, int NumSpecies)
{
	double dSum = 0.0;
	for ( int i=0; i<NumSpecies; i++ )
	{
		dSum += p0[i] + p1[i];
	}
	return(sqrt(dSum));
}


//
// Returns the Sum of the two integer vectors
//
double GetWeightsDefault(int *p0, int *p1, int NumSpecies)
{
	double dSum = 0.0;
	for ( int i=0; i<NumSpecies; i++ )
	{
		dSum += p0[i] + p1[i];
	}
	return(dSum);
}


//
// Returns the Square Root of the sum of the two double vectors
//
double GetWeightsDoubleSqrt(double *p0, double *p1, int NumSpecies)
{
	double dSum = 0.0;
	for ( int i=0; i<NumSpecies; i++ )
	{
		dSum += p0[i] + p1[i];
	}
	return(sqrt(dSum));
}


//
// Returns the sum of the two double vectors
//
double GetWeightsDoubleDefault(double *p0, double *p1, int NumSpecies)
{
	double dSum = 0.0;
	for ( int i=0; i<NumSpecies; i++ )
	{
		dSum += p0[i] + p1[i];
	}
	return(dSum);
}


//
// Return the number of present records at a site in a SxS Table
//
int GetCountForSite(int *p, int NumSpecies)
{
	int nSum = 0;
	for ( int i=0; i<NumSpecies; i++ )
	{
		nSum += p[i];
	}
	return(nSum);
}

