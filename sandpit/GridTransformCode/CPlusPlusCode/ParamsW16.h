//
// ParamsW16.h
//
#include "stdafx.h"

#ifndef __PARAMSW16_H__
#define __PARAMSW16_H__

#include "ProfileStrings.h"


///////////////////////////////////////////////////////////////////////////////////////
//
// functions for the [PREDICTORS] section
//
char *GetPredictorDomainPath(char *pParams);

int GetNumPredictors(char *pParams);

char *GetPredictorPathAt(char *pParams, int nIndex);

int GetPredictorInUseAt(char *pParams, int nIndex);

int GetPredictorTypeAt(char *pParams, int nIndex);

int GetPredictorDroppedAt(char *pParams, int nIndex);

bool GetPredictorNameAt(char *pParams, char *pBuff, int nIndex);

int GetEuclideanPredDroppedState(char *pParams);

int GetPredictorSplinesAt(char *pParams, int nIndex);

void SetPredictorSplineAt(char *pParams, int nPredIndex, int nSplineIndex, double dVal);

void SetPredictorCoeffAt(char *pParams, int nPredIndex, int nSplineIndex, double dVal);

int GetEuclideanSplines(char *pParams);

int GetEuclideanQuantType(char *pParams);

void SetEuclideanQuantType(char *pParams, int nType);

double GetEuclideanSplineAt(char *pParams, int nIndex);

void SetEuclideanSplineAt(char *pParams, int nIndex, double dVal);

void SetEuclideanCoeffAt(char *pParams, int nIndex, double dVal);

double GetEuclideanCoeffAt(char *pParams, int nIndex);

double GetPredictorSplineAt(char *pParams, int nPredIndex, int nSplineIndex);

double GetPredictorCoeffAt(char *pParams, int nPredIndex, int nSplineIndex);

char *GetEnvTableData(char *pParams);

void SetNumPredictors(char *pParams, int NumPreds);



///////////////////////////////////////////////////////////////////////////////////////
//
// functions for the [GDMODEL] section
//
char *GetWorkspacePath(char *pParams);

char *GetTmpRichnessPath(char *pParams);

char *GetResponseType(char *pParams);

char *GetPredictorType(char *pParams);

char *GetQuantileType(char *pParams);

bool UseEuclidean(char *pParams);

bool UseSubSample(char *pParams);

int GetNumSamples(char *pParams);

double GetIntercept(char *pParams);

double GetNullDeviance(char *pParams);

double GetGDMDeviance(char *pParams);

double GetDevianceExplained(char *pParams);

double GetSumOfCoefficients(char *pParams);

int GetNumDroppedPreds(char *pParams);

void SetIntercept(char *pParams, double dVal);

void SetNullDeviance(char *pParams, double dVal);

void SetGDMDeviance(char *pParams, double dVal);

void SetDevianceExplained(char *pParams, double dVal);

void SetSumOfCoefficients( char *pParams, double dVal );

void SetNumberOfActivePredictors( char *pParams, int nVal );

void SetNumberCoefficientsAboveZero( char *pParams, int nVal );

void SetNumberOfSitePairs( char *pParams, int nVal );

///////////////////////////////////////////////////////////////////////////////////////
//
// functions for the [RESPONSE] section
//
char *GetResponseData(char *pParams);
void SetFilteredResponseDataPath(char *pParams, char *path);
char *GetFilteredResponseData(char *pParams);
void SetFilteredSxSDataPath(char *pParams, char *path);
char *GetFilteredSxSData(char *pParams);

char *GetWeights(char *pParams);
bool DoAggregation(char *pParams);

char *GetResultsTable(char *pParams);

///////////////////////////////////////////////////////////////////////////////////////
//
// functions to test for return value validity
//
bool IsValidReturnString(char *lpVal);


//
// Clears geographic parameters in the [PREDICTORS] section
//
void ClearGeoPredictors(char *pParams);


//
// Append a blank row to a parameter file
//
void AppendRow(char *pParams);


//
// functions for the [GDMRESIDUALS] section
//
bool HoldEuclideanCoeffs(char *pParams);

bool GetHoldPredCoeffAt(char *pParams, int nIndex);

#endif //__PARAMSW16_H__

