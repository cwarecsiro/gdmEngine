//
// ParamsW16.cpp
//
#include "stdafx.h"
#include "ParamsW16.h"
#include "GDMBufferSizeDefs.h"


///////////////////////////////////////////////////////////////////////////////////////
//
// functions for the [PREDICTORS] section
//
char *GetPredictorDomainPath(char *pParams)
{
	char *lpBuffer = new char [FILEPATH_BUFFSIZE];
	GetProfileString( "PREDICTORS", "DomainGrid", lpBuffer, pParams );
	return(lpBuffer);
}



int GetNumPredictors(char *pParams)
{
	return(GetProfileInt( "PREDICTORS", "NumPredictors", pParams ));
}


char *GetPredictorPathAt(char *pParams, int nIndex)
{
	char tmp[32];
	sprintf( tmp, "EnvGrid%d", nIndex );

	char *lpBuffer = new char [FILEPATH_BUFFSIZE];
	GetProfileString( "PREDICTORS", tmp, lpBuffer, pParams );
	return(lpBuffer);
}


int GetPredictorInUseAt(char *pParams, int nIndex)
{
	char tmp[32];
	sprintf( tmp, "UseEnv%d", nIndex );
	return(GetProfileInt( "PREDICTORS", tmp, pParams ));
}


int GetPredictorTypeAt(char *pParams, int nIndex)
{
	char tmp[32];
	sprintf( tmp, "PredType%d", nIndex );
	return(GetProfileInt( "PREDICTORS", tmp, pParams ));
}


int GetPredictorDroppedAt(char *pParams, int nIndex)
{
	char tmp[32];
	sprintf( tmp, "DroppedPred%d", nIndex );
	return(GetProfileInt( "PREDICTORS", tmp, pParams ));
}


int GetEuclideanPredDroppedState(char *pParams)
{
	return(GetProfileInt( "PREDICTORS", "DroppedEuclPred", pParams ));
}


bool GetPredictorNameAt(char *pParams, char *pBuff, int nIndex)
{
	char tmp[32];
	sprintf( tmp, "EnvTab%d", nIndex );
	GetProfileString( "PREDICTORS", tmp, pBuff, pParams );
	return(true);
}


int GetPredictorSplinesAt(char *pParams, int nIndex)
{
	char tmp[32];
	sprintf( tmp, "PredSpl%d", nIndex );
	return(GetProfileInt( "PREDICTORS", tmp, pParams ));
}


void SetPredictorSplineAt(char *pParams, int nPredIndex, int nSplineIndex, double dVal)
{
	char tmp[32];
	sprintf( tmp, "PredSplVal%d.%d", nPredIndex, nSplineIndex );
	return(SetProfileDouble( "PREDICTORS", tmp, dVal, pParams ));
}


void SetPredictorCoeffAt(char *pParams, int nPredIndex, int nSplineIndex, double dVal)
{
	char tmp[32];
	sprintf( tmp, "PredCoef%d.%d", nPredIndex, nSplineIndex );
	return(SetProfileDouble( "PREDICTORS", tmp, dVal, pParams ));
}


int GetEuclideanSplines(char *pParams)
{
	return(GetProfileInt( "PREDICTORS", "EuclSpl", pParams ));
}


int GetEuclideanQuantType(char *pParams)
{
	return(GetProfileInt( "PREDICTORS", "EuclQuantType", pParams ));
}


void SetEuclideanQuantType(char *pParams, int nType)
{
	SetProfileInt( "PREDICTORS", "EuclQuantType", nType, pParams );
}


double GetEuclideanSplineAt(char *pParams, int nIndex)
{
	char tmp[32];
	sprintf( tmp, "EuclSplVal%d", nIndex );
	return(GetProfileDouble( "PREDICTORS", tmp, pParams ));
}


void SetEuclideanSplineAt(char *pParams, int nIndex, double dVal)
{
	char tmp[32];
	sprintf( tmp, "EuclSplVal%d", nIndex );
	SetProfileDouble( "PREDICTORS", tmp, dVal, pParams );
}


void SetEuclideanCoeffAt(char *pParams, int nIndex, double dVal)
{
	char tmp[32];
	sprintf( tmp, "EuclCoef%d", nIndex );
	SetProfileDouble( "PREDICTORS", tmp, dVal, pParams );
}


double GetEuclideanCoeffAt(char *pParams, int nIndex)
{
	char tmp[32];
	sprintf( tmp, "EuclCoef%d", nIndex );
	return(GetProfileDouble( "PREDICTORS", tmp, pParams ));
}


double GetPredictorSplineAt(char *pParams, int nPredIndex, int nSplineIndex)
{
	char tmp[32];
	sprintf( tmp, "PredSplVal%d.%d", nPredIndex, nSplineIndex );
	return(GetProfileDouble( "PREDICTORS", tmp, pParams ));
}


double GetPredictorCoeffAt(char *pParams, int nPredIndex, int nSplineIndex)
{
	char tmp[32];
	sprintf( tmp, "PredCoef%d.%d", nPredIndex, nSplineIndex );
	return(GetProfileDouble( "PREDICTORS", tmp, pParams ));
}


char *GetEnvTableData(char *pParams)
{
	char *lpBuffer = new char [FILEPATH_BUFFSIZE];
	GetProfileString( "PREDICTORS", "EnvTable", lpBuffer, pParams );
	return(lpBuffer);
}


void SetNumPredictors(char *pParams, int NumPreds)
{
	SetProfileInt( "PREDICTORS", "NumPredictors", NumPreds, pParams );
}

///////////////////////////////////////////////////////////////////////////////////////
//
// functions for the [GDMODEL] section
//
char *GetWorkspacePath(char *pParams)
{
	char *lpBuffer = new char [FILEPATH_BUFFSIZE];
	GetProfileString( "GDMODEL", "WorkspacePath", lpBuffer, pParams );
	return(lpBuffer);
}


char *GetTmpRichnessPath(char *pParams)
{
	char *lpBuffer = new char [FILEPATH_BUFFSIZE];
	GetProfileString( "GDMODEL", "WorkspacePath", lpBuffer, pParams );
	strcat(lpBuffer, "\\RichnessFilteredSXS.csv");
	return(lpBuffer);
}


char *GetResponseType(char *pParams)
{
	char *lpBuffer = new char [FILEPATH_BUFFSIZE];
	GetProfileString( "GDMODEL", "RespDataType", lpBuffer, pParams );
	return(lpBuffer);
}


char *GetPredictorType(char *pParams)
{
	char *lpBuffer = new char [FILEPATH_BUFFSIZE];
	GetProfileString( "GDMODEL", "PredDataType", lpBuffer, pParams );
	return(lpBuffer);
}


char *GetQuantileType(char *pParams)
{
	char *lpBuffer = new char [FILEPATH_BUFFSIZE];
	GetProfileString( "GDMODEL", "Quantiles", lpBuffer, pParams );
	return(lpBuffer);
}


bool UseEuclidean(char *pParams)
{
	return( 1 == GetProfileInt( "GDMODEL", "UseEuclidean", pParams ) );	
}


bool UseSubSample(char *pParams)
{
	return( 1 == GetProfileInt( "GDMODEL", "UseSubSample", pParams ) );	
}


int GetNumSamples(char *pParams)
{
	return( GetProfileInt( "GDMODEL", "NumSamples", pParams ) );	
}


double GetIntercept(char *pParams)
{
	return( GetProfileDouble( "GDMODEL", "Intercept", pParams ) );
}


double GetNullDeviance(char *pParams)
{
	return( GetProfileDouble( "GDMODEL", "NullDeviance", pParams ) );
}


double GetGDMDeviance(char *pParams)
{
	return( GetProfileDouble( "GDMODEL", "GDMDeviance", pParams ) );
}


double GetDevianceExplained(char *pParams)
{
	return( GetProfileDouble( "GDMODEL", "DevExplained", pParams ) );
}


double GetSumOfCoefficients(char *pParams)
{
	return( GetProfileDouble( "GDMODEL", "SumOfCoefficients", pParams ) );
}


int GetNumDroppedPreds(char *pParams)
{
	int nReturn = 0;
	int nPreds = GetNumPredictors(pParams);
	for ( int i=1; i<=nPreds; i++ )
	{
		if ( 1 == GetPredictorDroppedAt(pParams, i) )
			++nReturn;
	}
	return(nReturn);
}


void SetIntercept(char *pParams, double dVal)
{
	SetProfileDouble( "GDMODEL", "Intercept", dVal, pParams );
}


void SetNullDeviance(char *pParams, double dVal)
{
	SetProfileDouble( "GDMODEL", "NullDeviance", dVal, pParams );
}


void SetGDMDeviance(char *pParams, double dVal)
{
	SetProfileDouble( "GDMODEL", "GDMDeviance", dVal, pParams );
}


void SetDevianceExplained(char *pParams, double dVal)
{
	SetProfileDouble( "GDMODEL", "DevExplained", dVal, pParams );
}


void SetSumOfCoefficients( char *pParams, double dVal )
{
	SetProfileDouble( "GDMODEL", "SumOfCoefficients", dVal, pParams );
}


void SetNumberOfActivePredictors( char *pParams, int nVal )
{
	SetProfileInt( "GDMODEL", "NumberOfActivePredictors", nVal, pParams );
}


void SetNumberCoefficientsAboveZero( char *pParams, int nVal )
{
	SetProfileInt("GDMODEL", "NumberCoefficients>0", nVal, pParams );
}


void SetNumberOfSitePairs( char *pParams, int nVal )
{
	SetProfileInt("GDMODEL", "NumberSitePairs", nVal, pParams );
}



///////////////////////////////////////////////////////////////////////////////////////
//
// functions for the [RESPONSE] section
//
char *GetResponseData(char *pParams)
{
	char *lpBuffer = new char [FILEPATH_BUFFSIZE];
	GetProfileString( "RESPONSE", "InputData", lpBuffer, pParams );
	return(lpBuffer);
}


void SetFilteredResponseDataPath(char *pParams, char *path)
{
	SetProfileString( "RESPONSE", "FilteredInputData", path, pParams );
}


char *GetFilteredResponseData(char *pParams)
{
	char *lpBuffer = new char [FILEPATH_BUFFSIZE];
	GetProfileString( "RESPONSE", "FilteredInputData", lpBuffer, pParams );
	return(lpBuffer);
}


void SetFilteredSxSDataPath(char *pParams, char *path)
{
	SetProfileString( "RESPONSE", "FilteredSxSData", path, pParams );
}


char *GetFilteredSxSData(char *pParams)
{
	char *lpBuffer = new char [FILEPATH_BUFFSIZE];
	GetProfileString( "RESPONSE", "FilteredSxSData", lpBuffer, pParams );
	return(lpBuffer);
}


char *GetWeights(char *pParams)
{
	char *lpBuffer = new char [FILEPATH_BUFFSIZE];
	GetProfileString( "RESPONSE", "UseWeights", lpBuffer, pParams );
	return(lpBuffer);
}


bool DoAggregation(char *pParams)
{
	return( GetProfileInt( "RESPONSE", "Aggregate", pParams ) == 1 );
}




///////////////////////////////////////////////////////////////////////////////////////
//
// functions for the prediction results
//
char *GetResultsTable(char *pParams)
{
	char *lpBuffer = GetWorkspacePath(pParams);
	sprintf( lpBuffer, "%s\\ObservedVsPredicted.csv", GetWorkspacePath(pParams) );
	return(lpBuffer);
}


///////////////////////////////////////////////////////////////////////////////////////
//
// functions to test for return value validity
//
bool IsValidReturnString(char *lpVal)
{
	return(strlen(lpVal) > 0 );
}



//
// Clears the geographic parameters in the [PREDICTORS] section
//
void ClearGeoPredictors(char *pParams)
{
	WritePrivateProfileString( "PREDICTORS", "EuclSpl", NULL, pParams );
	WritePrivateProfileString( "PREDICTORS", "EuclMin", NULL, pParams );
	WritePrivateProfileString( "PREDICTORS", "EuclMid", NULL, pParams );
	WritePrivateProfileString( "PREDICTORS", "EuclMax", NULL, pParams );
	char lpKey[32];
	for ( int i=1; i<=3; i++ )
	{
		sprintf( lpKey, "EuclCoef%d", i );
		WritePrivateProfileString( "PREDICTORS", lpKey, NULL, pParams );
	}
	InsertNewLine( pParams );		
}


//
// Append a blank row to a parameter file
//
void AppendRow(char *pParams)
{
	FILE *fp = fopen(pParams, "a+t");
	fprintf(fp, "                                     \n" );
	fclose(fp);
}



///////////////////////////////////////////////////////////////////////////////////////
//
// functions for the [GDMRESIDUALS] section
//
bool HoldEuclideanCoeffs(char *pParams)
{
	return( 1 == GetProfileInt( "GDMRESIDUALS", "HoldEuclCoefficients", pParams ) );	
}


bool GetHoldPredCoeffAt(char *pParams, int nIndex)
{
	char tmp[32];
	sprintf( tmp, "HoldCoefficients%d", nIndex );
	return( 1 == GetProfileInt( "GDMRESIDUALS", tmp, pParams ) );	
}

