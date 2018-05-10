//
// GdmGuiLib.cpp
//
#include "stdafx.h"
#include "GdmBinLib.h"
#include "GdmGuiLib.h"
#include "GdmTabLib.h"
#include "NNLS_Double.h"
#include "Message.h"
#include "GDMBufferSizeDefs.h"
#include "ParamsW16.h"
#include "FilterSites.h"
#include "GdmExtractQuants.h"
#include "clsDoPath.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <io.h>
#include <fcntl.h>
#include <sys/stat.h>

#include <iostream>
#include <fstream>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Functions called from GUI /////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Show if 32 or 64 bit version of DLL
//
void TestDLLBitSize()
{
#if defined _M_X64

	Message("Compiled as 64 bit dll", "INFO");
	if (sizeof(int *) == 8)
	{
		Message("sizeof(int *) is 8", "INFO");
	}
	else if (sizeof(int *) == 4)
	{
		Message("sizeof(int *) is 4!!!!", "INFO");
	}

#elif defined _WIN32

	Message("Compiled as 32 bit dll", "INFO");
	if (sizeof(int *) == 8)
	{
		Message("sizeof(int *) is 8!!!!", "INFO");
	}
	else if (sizeof(int *) == 4)
	{
		Message("sizeof(int *) is 4", "INFO");
	}

#endif
}



//
// Fit a GDM model from residuals of a prior GDM model
//
bool GDM_FitFromResiduals( char *pParams, bool BatchRun, FPTR fptr )
{
	//
	// Get the data table path
	//
	char *DataTablePath = GetResponseData(pParams);
	//Message(DataTablePath, "DataTablePath");
	gmpath gmp;
	if (!gmp.FileExists(DataTablePath))
	{
		if (!BatchRun)
			Message("DataTablePath does NOT exist in GDM_FitFromResiduals", "ERROR");
		if (DataTablePath) delete[] DataTablePath;
		return(false);
	}


	//
	// Get a memory image of the composite data file
	//
	GDM_INT nTableRows = 0;
	GDM_INT nTableCols = 0;
	double *pTableData = GetTableDataPostExtractAndUpdateQuantiles( pParams, &nTableRows, &nTableCols, BatchRun, fptr);
	if (NULL == pTableData)
	{
		if (!BatchRun)
			Message("Cannot GetTableDataPostExtractAndUpdateQuantiles() in GDM_FitFromResiduals", "ERROR");
		if (DataTablePath) delete[] DataTablePath;
		return(false);
	}
	// Free the DataTablePath now that we don't need it anymore
	if (DataTablePath) delete[] DataTablePath;


	//
	// Get all the spline counts into a vector
	//
	int nSplineCounts = 0;
	int *theSplineCounts = GetSplineCountsFromParams(pParams, &nSplineCounts);


	//
	// Get all the quantiles into a vector 
	//
	int nQuantiles = 0;
	double *theQuantiles = GetQuantilesFromParams(pParams, &nQuantiles);

	
	//
	// Create a binary file image of the splined predictor data to pass to the matrix regression
	// There will always be the same number of rows in the table data as the predictor matrix
	// but the number of matrix cols may be less if not all predictors are 'In-use'.
	//
	GDM_INT nMatrixCols = 0;
	char *lpTmpFile = GetWorkspacePath(pParams);
	strcat(lpTmpFile, "\\gdmtmp.bin");

	if (!CreatePredictorBinary(lpTmpFile, pParams,
		                       pTableData, nTableRows, nTableCols, 
							   theSplineCounts, theQuantiles,
							   &nMatrixCols, BatchRun, fptr))
	{
		if (!BatchRun)
			Message("Cannot CreatePredictorBinary() in GDM_FitFromResiduals", "ERROR");
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		if (pTableData) delete[] pTableData;
		return(false);
	}
	PrintBinary2CSV(lpTmpFile, pParams, theSplineCounts, nTableRows, nMatrixCols, fptr);

	
	//
	// Extract the response column and the weights column before freeing the table data memory block
	//
	fptr("Extracting Response Vector...", 25);
	double *pResponse = new double [nTableRows];
	memcpy(pResponse, &pTableData[COL_RESPONSE * nTableRows], nTableRows * sizeof(double));
	if (NULL == pResponse)
	{
		if (!BatchRun)
			Message("Cannot extract Response column in GDM_FitFromResiduals", "ERROR");
		if (pTableData) delete[] pTableData;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}

	fptr("Extracting Weights Vector...", 50);
	double *pWeights = new double [nTableRows];
	memcpy(pWeights, &pTableData[COL_WEIGHTS * nTableRows], nTableRows * sizeof(double));
	if (NULL == pWeights)
	{
		if (!BatchRun)
			Message("Cannot extract Weights column in GDM_FitFromResiduals", "ERROR");
		if (pResponse) delete[] pResponse;
		if (pTableData) delete[] pTableData;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}

	// Free the Table data now that we don't need it anymore
	if (pTableData) delete[] pTableData;	

	
	//
	// Create the memory for the Predictor Matrix
	//
	fptr("Allocating Predictor Matrix...", 75);
	double *pPredictorMatrix = new double [nTableRows * nMatrixCols];
	if (NULL == pPredictorMatrix)
	{
		if (!BatchRun)
			Message("Cannot allocate pPredictorMatrix in GDM_FitFromResiduals", "ERROR");
		if (pResponse) delete[] pResponse;
		if (pWeights) delete[] pWeights;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}


	//
	// Populate the Predictor Matrix from the binary image created in CreatePredictorBinary()
	//
	int h = _open( lpTmpFile, _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE );
	if (h < 0)
	{
		if (!BatchRun)
			Message("Cannot open binary file for READ in GDM_FitFromResiduals", "ERROR");
		if (pPredictorMatrix) delete[] pPredictorMatrix;
		if (pResponse) delete[] pResponse;
		if (pWeights) delete[] pWeights;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}

	// incremental read
	double *pTmp = pPredictorMatrix;
	for ( int i=0; i<nMatrixCols; i++ )
	{
		_read(h, pTmp, (int)nTableRows * sizeof( double ));
		pTmp += nTableRows;
	}
	_close(h);


	//
	// Create the Coefficient Vector for the GDM Residual Model
	// Columns that do NOT have their coefficients help constant such as the intercept column ([0]), will have a value of -1.0. 
	// Values other than this will be used to reset the pNN vector, in the matrix regression to their parameter file coefficients.
	//
	double *pResidualCoefficients = new double [nMatrixCols];
	pResidualCoefficients[0] = -1.0;
	int nThisIndex = 1;
	if (UseEuclidean(pParams))
	{
		int nEuclSplines = GetEuclideanSplines(pParams);
		if (HoldEuclideanCoeffs(pParams))
		{
			for ( int i=0; i<nEuclSplines; i++ )
			{
				pResidualCoefficients[nThisIndex+i] = GetEuclideanCoeffAt(pParams, i+1);
			}
		}
		else
		{
			for ( int i=0; i<nEuclSplines; i++ )
			{
				pResidualCoefficients[nThisIndex+i] = -1.0;
			}
		}
		nThisIndex += nEuclSplines;
	}

	//
	// now do the environmental predictors
	//
	int nPreds = GetNumPredictors(pParams);
	for (int p=1; p<=nPreds; p++ )
	{
		int nSplines = GetPredictorSplinesAt(pParams, p);
		if (GetPredictorInUseAt(pParams, p)) 
		{
			if (GetHoldPredCoeffAt(pParams, p))
			{
				for ( int i=0; i<nSplines; i++ )
				{
					pResidualCoefficients[nThisIndex+i] = GetPredictorCoeffAt(pParams, p, i+1);
				}
			}
			else
			{
				for ( int i=0; i<nSplines; i++ )
				{
					pResidualCoefficients[nThisIndex+i] = -1.0;
				}
			}
			// increment the coefficient index
			nThisIndex += nSplines;			
		}
	}

	char qqq[256];
	sprintf(qqq, "%s\\ResidualVec.csv", GetWorkspacePath(pParams));
	FILE *fpResidual = fopen(qqq, "w+t");
	for ( int i=0; i<nMatrixCols; i++ )
	{
		fprintf(fpResidual, "%d, %lf\n", i, pResidualCoefficients[i]);
	}
	fclose(fpResidual);

	//
	// Do the matrix regression
	//
	//Message("Doing GDM Regression");
	fptr("Performing GDM regression...", 0);
	double dGDMDeviance;
	double *pCoefficients = WeightedNNLSResidualRegression( lpTmpFile, 
															pResidualCoefficients,									    
		                                                    pPredictorMatrix, 
											                nTableRows, 
											                nMatrixCols, 
											                pResponse, 
											                &dGDMDeviance, 
											                pWeights,
													        fptr);
	if ( NULL == pCoefficients )
	{
		if (!BatchRun)
			Message("pCoefficients are NULL", "ERROR in GDM_FitFromResiduals");
		if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );
		if (pPredictorMatrix) delete[] pPredictorMatrix;
		if (pResponse) delete[] pResponse;
		if (pWeights) delete[] pWeights;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		return(false);
		fptr("Ready:", 0);
		return(false);
	}


	//
	// remove the temporary matrix file if it exists (and it should!!) now that we don't need it
	//
	if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );
	if (lpTmpFile) delete[] lpTmpFile;


	//
	// create a NULL model and return the deviance
	//	
	//Message("Doing NULL Regression");
	double dNullDeviance = GetWeightedNULLDeviance( nTableRows, pResponse, pWeights );


	//
	// calculate deviance explained as a percentage
	//
	double dDevianceExplained = ( 1.0 - ( dGDMDeviance / dNullDeviance ) ) * 100;


	//
	// write relevent outputs back to the parameter file
	//	
	fptr("Updating parameter file...", 0);
	UpdateParameterFile(pParams, dNullDeviance, dGDMDeviance, dDevianceExplained,
						pCoefficients[0], &pCoefficients[1], (int)nMatrixCols-1, theSplineCounts, (int)nTableRows, fptr);

	//
	// write prediction table
	//
	CreateObservedVsPredicted(pParams, pResponse, pPredictorMatrix, pCoefficients, (int)nTableRows, (int)nMatrixCols, fptr);


	//
	// Clean up
	//
	if (pResidualCoefficients) delete[] pResidualCoefficients;
	if (pCoefficients) delete[] pCoefficients;
	if (pPredictorMatrix) delete[] pPredictorMatrix;
	if (pResponse) delete[] pResponse;
	if (pWeights) delete[] pWeights;
	if (theSplineCounts) delete[] theSplineCounts; 
	if (theQuantiles) delete[] theQuantiles;
	fptr("Ready:", 0);
	return(true);
}



//
// Fit a GDM model from a composite table supplied via a parameter file from a .NET interface
//
bool GDM_FitFromParamFile( char *pParams, bool BatchRun, FPTR fptr )
{
	//
	// Get the data table path
	//
	char *DataTablePath = GetResponseData(pParams);
	//Message(DataTablePath, "DataTablePath");
	gmpath gmp;
	if (!gmp.FileExists(DataTablePath))
	{
		if (!BatchRun)
			Message("DataTablePath does NOT exist in GDM_FitFromParamFile", "ERROR");
		if (DataTablePath) delete[] DataTablePath;
		return(false);
	}


	//
	// Get a memory image of the composite data file
	//
	GDM_INT nTableRows = 0;
	GDM_INT nTableCols = 0;
	double *pTableData = GetTableDataPostExtractAndUpdateQuantiles( pParams, &nTableRows, &nTableCols, BatchRun, fptr);
	if (NULL == pTableData)
	{
		if (!BatchRun)
			Message("Cannot GetTableDataPostExtractAndUpdateQuantiles() in GDM_FitFromParamFile", "ERROR");
		if (DataTablePath) delete[] DataTablePath;
		return(false);
	}

	// Free the DataTablePath now that we don't need it anymore
	if (DataTablePath) delete[] DataTablePath;


	//
	// Get all the spline counts into a vector
	//
	int nSplineCounts = 0;
	int *theSplineCounts = GetSplineCountsFromParams(pParams, &nSplineCounts);


	//
	// Get all the quantiles into a vector 
	//
	int nQuantiles = 0;
	double *theQuantiles = GetQuantilesFromParams(pParams, &nQuantiles);

	
	//
	// Create a binary file image of the splined predictor data to pass to the matrix regression
	// There will always be the same number of rows in the table data as the predictor matrix
	// but the number of matrix cols may be less if not all predictors are 'In-use'.
	//
	GDM_INT nMatrixCols = 0;
	char *lpTmpFile = GetWorkspacePath(pParams);
	strcat(lpTmpFile, "\\gdmtmp.bin");

	if (!CreatePredictorBinary(lpTmpFile, pParams,
		                       pTableData, nTableRows, nTableCols, 
							   theSplineCounts, theQuantiles,
							   &nMatrixCols, BatchRun, fptr))
	{
		if (!BatchRun)
			Message("Cannot CreatePredictorBinary() in GDM_FitFromParamFile", "ERROR");
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		if (pTableData) delete[] pTableData;
		return(false);
	}
	PrintBinary2CSV(lpTmpFile, pParams, theSplineCounts, nTableRows, nMatrixCols, fptr);

	
	//
	// Extract the response column and the weights column before freeing the table data memory block
	//
	fptr("Extracting Response Vector...", 25);
	double *pResponse = new double [nTableRows];
	memcpy(pResponse, &pTableData[COL_RESPONSE * nTableRows], nTableRows * sizeof(double));
	if (NULL == pResponse)
	{
		if (!BatchRun)
			Message("Cannot extract Response column in GDM_FitFromParamFile", "ERROR");
		if (pTableData) delete[] pTableData;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}

	fptr("Extracting Weights Vector...", 50);
	double *pWeights = new double [nTableRows];
	memcpy(pWeights, &pTableData[COL_WEIGHTS * nTableRows], nTableRows * sizeof(double));
	if (NULL == pWeights)
	{
		if (!BatchRun)
			Message("Cannot extract Weights column in GDM_FitFromParamFile", "ERROR");
		if (pResponse) delete[] pResponse;
		if (pTableData) delete[] pTableData;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}

	// Free the Table data now that we don't need it anymore
	if (pTableData) delete[] pTableData;	

	
	//
	// Create the memory for the Predictor Matrix
	//
	fptr("Allocating Predictor Matrix...", 75);
	double *pPredictorMatrix = new double [nTableRows * nMatrixCols];
	if (NULL == pPredictorMatrix)
	{
		if (!BatchRun)
			Message("Cannot allocate pPredictorMatrix in GDM_FitFromParamFile", "ERROR");
		if (pResponse) delete[] pResponse;
		if (pWeights) delete[] pWeights;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}


	//
	// Populate the Predictor Matrix from the binary image created in CreatePredictorBinary()
	//
	int h = _open( lpTmpFile, _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE );
	if (h < 0)
	{
		if (!BatchRun)
			Message("Cannot open binary file for READ in GDM_FitFromParamFile", "ERROR");
		if (pPredictorMatrix) delete[] pPredictorMatrix;
		if (pResponse) delete[] pResponse;
		if (pWeights) delete[] pWeights;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}

	// incremental read
	double *pTmp = pPredictorMatrix;
	for ( int i=0; i<nMatrixCols; i++ )
	{
		_read(h, pTmp, (int)nTableRows * sizeof( double ));
		pTmp += (int)nTableRows;
	}
	_close(h);


	//
	// Do the matrix regression
	//
	//Message("Doing GDM Regression");
	fptr("Performing GDM regression...", 0);
	double dGDMDeviance;
	double *pCoefficients = WeightedNNLSRegression( lpTmpFile, 
		                                            pPredictorMatrix, 
		(int)nTableRows,
											        nMatrixCols, 
											        pResponse, 
											        &dGDMDeviance, 
											        pWeights,
													fptr);
	if ( NULL == pCoefficients )
	{
		if (!BatchRun)
			Message("pCoefficients are NULL", "ERROR in GDM_FitFromParamFile");
		if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );
		if (pPredictorMatrix) delete[] pPredictorMatrix;
		if (pResponse) delete[] pResponse;
		if (pWeights) delete[] pWeights;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		return(false);
		fptr("Ready:", 0);
		return(false);
	}


	//
	// remove the temporary matrix file if it exists (and it should!!) now that we don't need it
	//
	if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );
	if (lpTmpFile) delete[] lpTmpFile;


	//
	// create a NULL model and return the deviance
	//	
	//Message("Doing NULL Regression");
	double dNullDeviance = GetWeightedNULLDeviance( nTableRows, pResponse, pWeights );


	//
	// calculate deviance explained as a percentage
	//
	double dDevianceExplained = ( 1.0 - ( dGDMDeviance / dNullDeviance ) ) * 100;


	//
	// write relevent outputs back to the parameter file
	//	
	fptr("Updating parameter file...", 0);
	UpdateParameterFile(pParams, dNullDeviance, dGDMDeviance, dDevianceExplained,
						pCoefficients[0], &pCoefficients[1], (int)nMatrixCols-1, theSplineCounts, (int)nTableRows, fptr);

	//
	// write prediction table
	//
	CreateObservedVsPredicted(pParams, pResponse, pPredictorMatrix, pCoefficients, (int)nTableRows, (int)nMatrixCols, fptr);


	//
	// Clean up
	//
	if (pCoefficients) delete[] pCoefficients;
	if (pPredictorMatrix) delete[] pPredictorMatrix;
	if (pResponse) delete[] pResponse;
	if (pWeights) delete[] pWeights;
	if (theSplineCounts) delete[] theSplineCounts; 
	if (theQuantiles) delete[] theQuantiles;
	fptr("Ready:", 0);
	return(true);
}


//
// Count the number of rows in a GDM Composite Table 
// Called from the .NET interface in the Run GDM panel
//
int GDM_CountRowsInTable(char *pPath, FPTR fptr)
{
	int nRows = CountRows(pPath, fptr);
	return(nRows);
}


//
// Extract the Minimum and Maximum values from the 
// geographic distance columns in a GDM Composite Table
//
void GDM_GetGeographicMinMax(char *pPath, double *pMin, double *pMax, FPTR fptr)
{
	char *pRow = new char [TABLE_ROW_BUFFSIZE];
	if ( NULL == pRow )
	{
		Message("Cannot allocate pRow in GDM_GetGeographicMinMax", "ERROR");
		return;
	}
	
	// open text file
	FILE *fp = fopen(pPath, "r+t");
	if ( NULL == fp)
	{
		Message("Cannot open pPath in GDM_GetGeographicMinMax", "ERROR");
		if (pRow) delete[] pRow;
		return;
	}

	// get header
	fgets( pRow, TABLE_ROW_BUFFSIZE, fp );

	// extract geographic distance
	int nCurrent = 0;
	fptr("Extracting Min and Max Geographic Distance From Data...", nCurrent);
	int nThis = 0;
	while(1)
	{
		++nThis;
		if (0 == nThis % 50000)  // increment progress bar every 50000 rows...
		{
			if (++nCurrent > 100) nCurrent = 1;
			fptr("Extracting Min and Max Geographic Distance From Data...", nCurrent);
		}
		
		if ( NULL == fgets( pRow, TABLE_ROW_BUFFSIZE, fp ) )
			break;

		// get the Response
		char *p = strtok(pRow, ",\n");

		// get the Weights
		p = strtok(NULL, ",\n");

		// get the first site 
		p = strtok(NULL, ",\n");
		double dX0 = atof(p);
		p = strtok(NULL, ",\n");
		double dY0 = atof(p);

		// get the second site 
		p = strtok(NULL, ",\n");
		double dX1 = atof(p);
		p = strtok(NULL, ",\n");
		double dY1 = atof(p);

		double distX = fabs(dX1 - dX0);
        double distY = fabs(dY1 - dY0);

        // calculate geographic distance
        double dTmp = sqrt((distX * distX) + (distY * distY));

        if (dTmp < *pMin) *pMin = dTmp;
        if (dTmp > *pMax) *pMax = dTmp;
	}
	if (pRow) delete[] pRow;
	if (fp) fclose(fp);
}


//
// Extract the Minimum and Maximum values from the predictor
// columns at the ONE-BASED nIndex and nIndex + nPreds in a GDM Composite Table
//
void GDM_GetPredictorMinMax(char *pPath, int nIndex, int nPreds, double *pMin, double *pMax, FPTR fptr)
{
	char *pRow = new char [TABLE_ROW_BUFFSIZE];
	if ( NULL == pRow )
	{
		Message("Cannot allocate p in GDM_GetPredictorMinMax", "ERROR");
		return;
	}
	
	// open text file
	FILE *fp = fopen(pPath, "r+t");
	if ( NULL == fp)
	{
		Message("Cannot open pPath in GDM_GetPredictorMinMax", "ERROR");
		if (pRow) delete[] pRow;
		return;
	}

	// get header
	fgets( pRow, TABLE_ROW_BUFFSIZE, fp );

	// extract predictor values
	double dTmp;
	int nCurrent = 0;
	fptr("Extracting Min and Max Values from Predictor...", nCurrent);
	int nThis = 0;
	while(1)
	{
		++nThis;
		if (0 == nThis % 50000)  // increment progress bar every 50000 rows...
		{
			if (++nCurrent > 100) nCurrent = 1;
			fptr("Extracting Min and Max Values from Predictor...", nCurrent);
		}
		
		if ( NULL == fgets( pRow, TABLE_ROW_BUFFSIZE, fp ) )
			break;

		// get the Response
		char *p = strtok(pRow, ",\n");

		// get the Weights
		p = strtok(NULL, ",\n");

		// get the first site 
		p = strtok(NULL, ",\n");
		p = strtok(NULL, ",\n");

		// get the second site 
		p = strtok(NULL, ",\n");
		p = strtok(NULL, ",\n");

		// get to the predictor for the first site at nIndex
		for (int i=0; i<nIndex; i++ )
		{
			p = strtok(NULL, ",\n");
		}
		dTmp = atof(p); 
        if (dTmp < *pMin) *pMin = dTmp;
        if (dTmp > *pMax) *pMax = dTmp;

		// get to the predictor for the last site at nIndex
		for (int i=0; i<nPreds; i++ )
		{
			p = strtok(NULL, ",\n");
		}
		dTmp = atof(p); 
        if (dTmp < *pMin) *pMin = dTmp;
        if (dTmp > *pMax) *pMax = dTmp;
	}
	if (pRow) delete[] pRow;
	if (fp) fclose(fp);
}



//
// Get a memory image of the composite data file
//
double *GetCompositeTableIntoMemory(char *pDataPath, GDM_INT *pTableRows, GDM_INT *pTableCols, bool BatchRun, FPTR fptr)
{
	//
	// Get the table metrics
	//
	int nRows = CountRows64(pDataPath, fptr);
	if ( nRows < 1 )
	{
		if (!BatchRun)
			Message("Cannot CountRows in GetCompositeTableIntoMemory()", "ERROR");
		*pTableRows = *pTableCols = 0;
		return(NULL);
	}

	//Message(nRows, "nRows");

	int nCols = CountColumns64(pDataPath, fptr);
	if ( nCols < 1 )
	{
		if (!BatchRun)
			Message("Cannot CountColumns in GetCompositeTableIntoMemory()", "ERROR");
		*pTableRows = *pTableCols = 0;
		return(NULL);
	}

	//Message(nCols, "nCols");
	//Message(nRows * nCols, "nRows * nCols");
	//Message(nRows * nCols * 8, "nRows * nCols * 8");

	if (sizeof(int *) == 4)  // 32 bit OS
	{
		if ((nRows * nCols * 8) <= 0)
		{
			if (!BatchRun)
			{
				Message(nRows, "NumRows");
				Message(nCols, "NumCols");
				Message("Cannot allocate pData in GetCompositeTableIntoMemory", "ERROR");
			}
			return(NULL);
		}
	}

	double *pData = new double [nRows * nCols];
	if (NULL == pData)
	{
		if (!BatchRun)
			Message("Cannot allocate pData in GetCompositeTableIntoMemory", "ERROR");
		return(NULL);
	}

	char *pRowData = new char [TABLE_ROW_BUFFSIZE];
	if (NULL == pRowData)
	{
		if (!BatchRun)
			Message("Cannot allocate pRowData in GetCompositeTableIntoMemory", "ERROR");
		if (pData) delete[] pData;
		return(NULL);
	}


	// open text file
	ifstream* pmyFile = new ifstream; // On the heap
	pmyFile->open( pDataPath );

	if (!pmyFile->is_open())
	{
		if (!BatchRun)
			Message("Cannot open pDataPath in GetCompositeTableIntoMemory", "ERROR");
		if (pData) delete[] pData;
		if (pRowData) delete[] pRowData;
		return(NULL);
	}

	// get header
	pmyFile->getline(pRowData,TABLE_ROW_BUFFSIZE);

	int nCurrent = 0;
	char pBuff[128];
	fptr("Getting composite table into memory", nCurrent);
	for ( int i=0; i<nRows; i++ )
	{
		if ( (long long)(i) * 100 / nRows > nCurrent )
		{
			nCurrent = (long long)(i) * 100 / nRows;
			sprintf(pBuff, "Getting composite table into memory (Row %d of %d)", i+1, nRows);
			//fptr("Getting composite table into memory", nCurrent);
			fptr(pBuff, nCurrent);
		}

		// get current row
		pmyFile->getline(pRowData,TABLE_ROW_BUFFSIZE);

		// extract data from current row
		char *p = strtok(pRowData, ",\n");
		if (NULL == p)
		{
			char qqq[64];
			sprintf(qqq, "nRows: %d   nCols: %d", nRows, nCols );
			Message(qqq);
			sprintf(qqq, "Got a NULL at Row: %d and Col: %d", i,0);
			Message(qqq);

			pmyFile->close();
			if (pRowData) delete[] pRowData;
			*pTableRows = 0;
			*pTableCols = 0;
			return(NULL);
		}
		pData[i] = atof(p);

		for ( int j=1; j<nCols; j++ )
		{
			p = strtok(NULL, ",\n");

			if (NULL == p)
			{
				char qqq[64];
				sprintf(qqq, "nRows: %d   nCols: %d", nRows, nCols );
				Message(qqq);
				sprintf(qqq, "Got a NULL at Row: %d and Col: %d", i,j);
				Message(qqq);

				pmyFile->close();
				if (pRowData) delete[] pRowData;
				*pTableRows = 0;
				*pTableCols = 0;
				return(NULL);
			}
			pData[(j*nRows)+i] = atof(p);
		}
	}

	fptr("Getting composite table into memory", 0);
	pmyFile->close();
	if (pRowData) delete[] pRowData;
	*pTableRows = (GDM_INT)nRows;
	*pTableCols = (GDM_INT)nCols;
	return(pData);
}


//
// Create a binary file image of the splined predictor data to pass to the matrix regression
// There will always be the same number of rows in the table data as the predictor matrix
// but the number of matrix cols may be less if not all predictors are 'In-use'.
//
bool CreatePredictorBinary(char *lpBinFile, char *pParams,
	                       double *pTableData, GDM_INT nTableRows, GDM_INT nTableCols, 
						   int *pSplineCounts, double *pQuantiles,
						   GDM_INT *pMatrixCols, bool BatchRun, FPTR fptr)
{
	int nCurrent = 0;
	if (fptr)
		fptr("Creating Predictor Binary...", nCurrent);

	//
	// Create the GDM binary data file
	//
	ofstream *out = new ofstream(lpBinFile, ios::binary | ios::trunc);
	if (!out->is_open())
	{
		if (!BatchRun)
			Message("Cannot create binary file for in CreatePredictorBinary()", "ERROR");
		return(false);
	}

	//
	// allocate temporary column memory
	//
	double *pTmpCol = new double [nTableRows];
	if (NULL == pTmpCol)
	{
		if (!BatchRun)
			Message("Cannot allocate pTmpCol for in CreatePredictorBinary()", "ERROR");
		out->close();
		return(false);
	}

	//
	// write the intercept columns of ones...
	//
	nCurrent = 0;
	if (fptr)
		fptr("Creating Predictor Binary...", nCurrent);
	for ( int i=0; i<nTableRows; i++ ) 
	{
		if ((long long)(i) * 100 / nTableRows > nCurrent)
		{
			nCurrent = (long long)(i) * 100 / (int)nTableRows;
			if (fptr)
				fptr("Creating Predictor Binary...", nCurrent);
		}

		pTmpCol[i] = 1.0;
	}
	out->write((char *)pTmpCol, (int)nTableRows * sizeof(double)); // write the column
	*pMatrixCols = 1;	                                      // increment the column count 

	//
	// Setup some pointers to the spline count vector and the quantile vector
	//
	int *pSplines = pSplineCounts;
	double *pQuants = pQuantiles;
	double dMin,dMid,dMax,dVal0,dVal1;

	//
	// Do the geographic distance predictor if defined
	//
	int nSplines = *pSplines;
	if (UseEuclidean(pParams))
	{		
		for ( int j=0; j<nSplines; j++ )
		{
			if (j==0)		          // first spline
			{
				dMin = pQuants[j+0];
				dMid = pQuants[j+0];
				dMax = pQuants[j+1];
			}

			else if (j==nSplines-1)   // last spline
			{
				dMin = pQuants[j-1];
				dMid = pQuants[j+0];
				dMax = pQuants[j+0];
			}

			else					  // an interior spline
			{
				dMin = pQuants[j-1];
				dMid = pQuants[j+0];
				dMax = pQuants[j+1];
			}

			nCurrent = 0;
			if (fptr)
				fptr("Creating Predictor Binary using Geographic Distance...", nCurrent);
			for ( int i=0; i<nTableRows; i++ ) 
			{
				if ((long long)(i) * 100 / nTableRows > nCurrent)
				{
					nCurrent = (long long)(i) * 100 / (int)nTableRows;
					if (fptr)
						fptr("Creating Predictor Binary using Geographic Distance...", nCurrent);
				}

				double distX = fabs(pTableData[(COL_SITE1_X0 * nTableRows)+i] - pTableData[(COL_SITE2_X1 * nTableRows)+i]);
				double distY = fabs(pTableData[(COL_SITE1_Y0 * nTableRows)+i] - pTableData[(COL_SITE2_Y1 * nTableRows)+i]);
				dVal0 = 0.0;
				dVal1 = sqrt((distX * distX) + (distY * distY));
			
				// calculate the I-Splines values for each site in the pair
				double d0 = DoSplineCalc(dVal0, dMin, dMid, dMax);
				double d1 = DoSplineCalc(dVal1, dMin, dMid, dMax);
				pTmpCol[i] = fabs(d0 - d1);
			} // for ( int i=0; i<nTableRows; i++ ) 

			out->write((char *)pTmpCol, nTableRows * sizeof(double)); // write the column
			*pMatrixCols += 1;						                  // increment the column count 
		} // for ( int j=0; j<nSplines; j++ )

	} // if (UseEuclidean(pParams))

	//
	// reset pointers ready for the environmental predictors to follow...
	//
	pQuants += nSplines;
	pSplines += 1;

	//
	// Do the environmental predictors
	//
	nCurrent = 0;
	char pBuff[128];
	int nPreds = GetNumPredictors(pParams);
	sprintf(pBuff, "Creating Predictor Binary From Pred %d of %d", 1, nPreds);
	if (fptr)
		fptr(pBuff, nCurrent);
	for ( int p=0; p<nPreds; p++ )
	{
		if (p * 100 / nPreds > nCurrent)
		{
			nCurrent = p * 100 / nPreds;
			sprintf(pBuff, "Creating Predictor Binary From Pred %d of %d", p+1, nPreds);
			if (fptr)
				fptr(pBuff, nCurrent);
		}

		int nSplines = *pSplines;

		if (GetPredictorInUseAt(pParams, p+1))
		{
			int Offset0 = (LEADING_COLS+p)*(int)nTableRows;
			int Offset1 = (LEADING_COLS+p+nPreds)*(int)nTableRows;

			for ( int j=0; j<nSplines; j++ )
			{
				if (j==0)		          // first spline
				{
					dMin = pQuants[j+0];
					dMid = pQuants[j+0];
					dMax = pQuants[j+1];
				}

				else if (j==nSplines-1)   // last spline
				{
					dMin = pQuants[j-1];
					dMid = pQuants[j+0];
					dMax = pQuants[j+0];
				}

				else					  // an interior spline
				{
					dMin = pQuants[j-1];
					dMid = pQuants[j+0];
					dMax = pQuants[j+1];
				}

				for ( int i=0; i<nTableRows; i++ ) 
				{
					dVal0 = pTableData[Offset0+i];
					dVal1 = pTableData[Offset1+i];
					
					// calculate the I-Splines values for each site in the pair
					double d0 = DoSplineCalc(dVal0, dMin, dMid, dMax);
					double d1 = DoSplineCalc(dVal1, dMin, dMid, dMax);
					pTmpCol[i] = fabs(d0 - d1);

				} // for ( int i=0; i<nTableRows; i++ ) 

				out->write((char *)pTmpCol, nTableRows * sizeof(double)); // write the column
				*pMatrixCols += 1;						          // increment the column count 

			} // for ( int j=0; j<nSplines; j++ )

		} // if (GetPredictorInUseAt(pParams, p+1))

		// increment pointers
		pQuants += nSplines;
		pSplines += 1;
	} // for ( int p=0; p<nPreds; p++ )
	if (fptr)
		fptr("Creating Predictor Binary...", nCurrent);

	//
	// cleanup
	//
	out->close();
	if (pTmpCol) delete[] pTmpCol;
	return(true);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Local Functions...
//
//
//
// Convert the binary matrix created in CreatePredictorBinary to a CSV file
//
void PrintBinary2CSV(char *lpBinFile, char *pParams, 
	                 int *pSplineCounts, GDM_INT nTableRows, GDM_INT nMatrixCols, FPTR fptr)
{
	char *pRowData = new char [TABLE_ROW_BUFFSIZE];
	char *pOutPath = new char [FILEPATH_BUFFSIZE];
	sprintf(pOutPath, "%s\\SplineData2CSV.csv", GetWorkspacePath(pParams));
	//Message(pOutPath, "pOutPath");
	int *pSplines = pSplineCounts;
	
	//
	// print the header
	//
	FILE *fp = fopen(pOutPath, "w+t");
	fprintf(fp, "Intercept");	
	if (UseEuclidean(pParams))
	{
		int nSplines = *pSplines;
		for ( int i=1; i<=nSplines; i++ )
		{
			fprintf(fp, ",Euclidean_%d", i);
		}
	}
	pSplines += 1;

	int nPreds = GetNumPredictors(pParams);
	for ( int i=1; i<=nPreds; i++ )
	{
		int nSplines = *pSplines;

		if (GetPredictorInUseAt(pParams,i))
		{
			char lpKey[64];
			sprintf(lpKey, "EnvTab%d", i);
			GetProfileString( "PREDICTORS", lpKey, pRowData, pParams );

			for ( int j=1; j<=nSplines; j++ )
			{
				
				fprintf(fp, ",%s_%d", pRowData, j);
			}
		}
		pSplines += 1;
	}
	fprintf(fp, "\n");

	//Message(nTableRows, "nTableRows");
	//Message(nMatrixCols, "nMatrixCols");

	//
	// print the model data
	//
	int h = _open( lpBinFile, _O_BINARY | _O_RDWR, S_IREAD | S_IWRITE  );
	if ( h < 0 )
	{
		Message( "Cannot open binary predictor file for READ in PrintBinary2CSV", "ERROR" );
		if (fp) fclose(fp);
		if (pRowData) delete[] pRowData;
		if (pOutPath) delete[] pOutPath;
		return;
	}

	double dVal;
	for (int i=0; i<nTableRows; i++ )
	{
		for ( int j=0; j<nMatrixCols; j++ )
		{
			_lseek(h,(long)((j*nTableRows)+i)*sizeof(double), SEEK_SET);
			_read(h, &dVal, sizeof(double));
			fprintf(fp, "%lf", dVal);
			if ( j <nMatrixCols-1 )
				fprintf(fp, ",");
			else
				fprintf(fp, "\n");
		}
	}

	//
	// cleanup
	//
	_close(h);
	if (fp) fclose(fp);
	if (pRowData) delete[] pRowData;
	if (pOutPath) delete[] pOutPath;
}


//
// Dump the contents of the table memory back to file for debugging
//
void PrintTableData(char *pParams, double *pData, int nRows, int nCols, bool BatchRun)
{
	char *pHeader = new char [TABLE_ROW_BUFFSIZE];
	if (NULL == pHeader)
	{
		if (!BatchRun)
			Message("Cannot allocate pRowData in PrintTableData", "ERROR");
		return;
	}
	
	sprintf( pHeader, "%s\\MemoryTable.csv", GetWorkspacePath(pParams));
	FILE *fpOut = fopen(pHeader, "w+t");
	if (NULL == fpOut)
	{
		if (!BatchRun)
			Message("Cannot open fpOut in PrintTableData", "ERROR");
		if (pHeader) delete[] pHeader;
		return;
	}

	FILE *fpIn = fopen(GetResponseData(pParams), "r+t");
	if (NULL == fpIn)
	{
		if (!BatchRun)
			Message("Cannot open fpIn in PrintTableData", "ERROR");
		if (fpOut) fclose(fpOut);
		if (pHeader) delete[] pHeader;
		return;
	}

	fgets(pHeader, TABLE_ROW_BUFFSIZE, fpIn);
	fprintf(fpOut, "%s", pHeader);
	if (fpIn) fclose(fpIn);

	for ( int i=0; i<nRows; i++ )
	{
		for ( int j=0; j<nCols; j++ )
		{
			fprintf(fpOut, "%lf", pData[(j*nRows)+i]);
			if (j < nCols-1)
				fprintf(fpOut, ",");
			else
				fprintf(fpOut, "\n");
		}		
	}

	if (fpOut) fclose(fpOut);
	if (pHeader) delete[] pHeader;
}


//
// Dump the contents of the quantile vector memory back to file for debugging
//
void PrintQuantileVector(char *pParams,  double *pQuantiles, int nNumQuants, bool BatchRun)
{
	char *p = new char [TABLE_ROW_BUFFSIZE];
	if (NULL == p)
	{
		if (!BatchRun)
			Message("Cannot allocate p in PrintQuantileVector", "ERROR");
		return;
	}
	
	sprintf( p, "%s\\QuantileVector.csv", GetWorkspacePath(pParams));
	FILE *fp = fopen(p, "w+t");
	if (NULL == fp)
	{
		if (!BatchRun)
			Message("Cannot open fp in PrintQuantileVector", "ERROR");
		if (p) delete[] p;
		return;
	}

	fprintf( fp, "Index,Quantile\n");
	for ( int i=0; i<nNumQuants; i++ )
	{
		fprintf( fp, "%d,%lf\n", i+1, pQuantiles[i]);
	}

	if ( fp) fclose(fp);
	if (p) delete[] p;
}


//
// Dump the GDM results to text file for debugging
//
void PrintGDMResults(char *pParams, double dDevianceExplained, int nSplines, double *pCoeffs, bool BatchRun)
{
	char *p = new char [TABLE_ROW_BUFFSIZE];
	if (NULL == p)
	{
		if (!BatchRun)
			Message("Cannot allocate p in PrintQuantileVector", "ERROR");
		return;
	}
	
	sprintf( p, "%s\\GDMResults.csv", GetWorkspacePath(pParams));
	FILE *fp = fopen(p, "w+t");
	if (NULL == fp)
	{
		if (!BatchRun)
			Message("Cannot open fp in PrintGDMResults", "ERROR");
		if (p) delete[] p;
		return;
	}

	fprintf( fp, "DevExpl.,%lf\n", dDevianceExplained);
	fprintf( fp, "Index,Coefficient\n");
	for ( int i=0; i<nSplines; i++ )
	{
		if ( i==0 )
			fprintf( fp, "Intercept,%lf\n", pCoeffs[i]);
		else
			fprintf( fp, "%d,%lf\n", i, pCoeffs[i]);
	}

	if ( fp) fclose(fp);
	if (p) delete[] p;
}


 

//
// Sum the GDM coefficients
//
double CalcCoefficientSum(double *pCoefficients, int nTotalSplines)
{
	double dVal = 0.0;
	for ( int i=0; i<nTotalSplines; i++ )
	{
		dVal += pCoefficients[i];
	}
	return(dVal);
}


//
// Sum the number of active predictors
//
int CalcNumberOfActivePredictors(char *pParams, double *pCoefficients, int *pSplineCounts )
{
	if ( 0 == GetNumActivePredictors(pParams) )
	{
		return(0);
	}

	int nItems = 0;
	int nPreds = GetNumPredictors(pParams);
	double *pCoeff = pCoefficients;
	for ( int i=0; i<nPreds; i++ )
	{
		if (1 == GetPredictorInUseAt(pParams, i+1))
		{
			double dSum = 0.0;
			int nSplines = pSplineCounts[i];
			for ( int j=0; j<nSplines; j++ )
			{
				dSum += pCoeff[j];
			}
			if ( dSum > 0.0 ) ++nItems;

			// increment pointer
			pCoeff += nSplines;
		}
	}
	if (UseEuclidean(pParams)) ++nItems;
	return(nItems);
}


//
// Get the number of active predictors in the current GDM model
//
int GetNumActivePredictors(char *pParams)
{
	int nActive = 0;
	int nPreds = GetNumPredictors(pParams);
	for ( int i=1; i<=nPreds; i++ )
	{
		nActive += GetPredictorInUseAt(pParams, i);
	}
	if (UseEuclidean(pParams)) ++nActive;
	return(nActive);
}


//
// Sum the number of coefficients > 0
//
int CalcNumberofCoeffGreaterThanZero(double *pCoefficients, int nTotalSplines)
{
	int nItems = 0;
	for ( int i=0; i<nTotalSplines; i++ )
	{
		if ( pCoefficients[i] > 0.0 ) ++nItems;
	}
	return(nItems);
}


//
// Get the spline counts for all the predictors including geographic distance 
// as a vector from the parameter file (ignore whether they are "in-use" or not)
//
int *GetSplineCountsFromParams(char *pParams, int *pItems)
{
	int nItems = GetNumPredictors(pParams) + 1;
	int *pSplineCounts = new int [nItems];

	pSplineCounts[0] = GetEuclideanSplines(pParams);
	for (int i=1; i<nItems; i++ ) 
	{
		pSplineCounts[i] = GetPredictorSplinesAt(pParams, i);
	}

	*pItems = nItems;
	return(pSplineCounts);
}


//
// Get the quantiles for all the predictors including geographic distance 
// as a vector from the parameter file (ignore whether they are "in-use" or not)
//
double *GetQuantilesFromParams(char *pParams, int *pItems)
{
	// count the total number of possible splines in the model
	int nSplines = GetEuclideanSplines(pParams);
	int nPreds = GetNumPredictors(pParams);
	for ( int i=1; i<=nPreds; i++ )
	{
		nSplines += GetPredictorSplinesAt(pParams, i);
	}

	// allocate and populate the Quantiles vector
	double *pQuantiles = new double [nSplines];

	int nThis = 0;
	for (int i=1; i<=GetEuclideanSplines(pParams); i++ )
	{
		pQuantiles[nThis++] = GetEuclideanSplineAt(pParams, i);
	}

	for ( int i=1; i<=nPreds; i++ )
	{
		for ( int j=1; j<=GetPredictorSplinesAt(pParams, i); j++ )
		{
			pQuantiles[nThis++] = GetPredictorSplineAt(pParams, i, j);
		}
	}

	*pItems = nSplines;
	return(pQuantiles);
}


//
// Write the observed and predicted data to comma delimited textfile and to a pair of binary files for plotting
//
void CreateObservedVsPredicted(char *pParams, double *pResponse, double *pPredData, double *pCoefficients, int nRows, int nCols, FPTR fptr)
{
	int nCurrent = 0;
	fptr("Writing predictions to file...", nCurrent );

	// open the binary files
	char *buff = new char [FILEPATH_BUFFSIZE];
	sprintf( buff, "%s\\Observed.bin", GetWorkspacePath(pParams) );
	ofstream *ObsOut = new ofstream(buff, ios::binary | ios::trunc);
	if (!ObsOut->is_open())
	{
		Message("Cannot create Observed.bin in CreateObservedVsPredicted()", "ERROR");
		return;
	}
	
	sprintf( buff, "%s\\Predicted.bin", GetWorkspacePath(pParams) );
	ofstream *PredOut = new ofstream(buff, ios::binary | ios::trunc);
	if (!PredOut->is_open())
	{
		Message("Cannot create Predicted.bin in CreateObservedVsPredicted()", "ERROR");
		ObsOut->close();
		return;
	}

	sprintf( buff, "%s\\GDM_SixCols.bin", GetWorkspacePath(pParams) );
	int hIn = _open( buff, _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE  );
	double myDouble;
	
	sprintf( buff, "%s\\ObservedVsPredicted.csv", GetWorkspacePath(pParams) );
	FILE *fp = fopen( buff, "w+t" );

	//fprintf( fp, "Observed,Predicted,Linked,Residual,Pearsons\n");
	fprintf( fp, "Observed,Predicted,Linked,Residual,Pearsons,Response,Weights,X0,Y0,X1,Y1\n");
	for ( int i=0; i<nRows; i++ )
	{
		if ( (long long)(i) * 100 / nRows > nCurrent )
		{
			nCurrent = (long long)(i) * 100 / nRows;
			fptr("Writing predictions to file...", nCurrent );
		}

		// write Observed response
		double dObserved = pResponse[i];
		fprintf( fp, "%lf,", dObserved);

		// write Predicted response
		double dPredicted = CalcDissimilarity( pPredData, pCoefficients, nRows, nCols, i );
		fprintf( fp, "%lf,", dPredicted);

		// write Linked predicted response
		double dLinked = 1.0 - exp(-dPredicted);
		fprintf( fp, "%lf,", dLinked);

		// write standard residual
		double dResidual = dObserved - dLinked;
		fprintf( fp, "%lf,", dResidual );

		// write Pearson residial
		//fprintf( fp, "%lf\n", dResidual / sqrt(dLinked * (1.0 - dLinked)) );
		fprintf( fp, "%lf,", dResidual / sqrt(dLinked * (1.0 - dLinked)) );

		// write response
		_read(hIn, &myDouble, sizeof(double));
		fprintf( fp, "%lf,", myDouble);
		
		// write weights
		_read(hIn, &myDouble, sizeof(double));
		fprintf( fp, "%lf,", myDouble);

		// write X0
		_read(hIn, &myDouble, sizeof(double));
		fprintf( fp, "%lf,", myDouble);

		// write Y0
		_read(hIn, &myDouble, sizeof(double));
		fprintf( fp, "%lf,", myDouble);

		// write X1
		_read(hIn, &myDouble, sizeof(double));
		fprintf( fp, "%lf,", myDouble);

		// write Y1
		_read(hIn, &myDouble, sizeof(double));
		fprintf( fp, "%lf\n", myDouble);

		// update the binary files
		ObsOut->write((char *)&dObserved, sizeof(double));
		PredOut->write((char *)&dPredicted, sizeof(double));
	}

	_close(hIn);
	fclose(fp);
	ObsOut->close();
	PredOut->close();
	if (buff) delete[] buff;
	nCurrent = 0;
	fptr("Writing predictions to file...", nCurrent );
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//														For GDM Predict
//
//
// Do predictions from a GDM model from a composite table supplied via a parameter file from a .NET interface
//
bool GDM_PredictFromParamFile( char *pModelParams, char *pInTablePath, char *pOutPath, FPTR fptr )
{
	//
	// Get all the spline counts into a vector
	//
	int nSplineCounts = 0;
	int *theSplineCounts = GetSplineCountsFromParams(pModelParams, &nSplineCounts);


	//
	// Get all the quantiles into a vector 
	//
	int nQuantiles = 0;
	double *theQuantiles = GetQuantilesFromParams(pModelParams, &nQuantiles);


	//
	// Get a memory image of the all the columns in the composite data file
	//
	GDM_INT nTableRows = 0;
	GDM_INT nTableCols = 0;
	double *pTableData = GetCompositeTableIntoMemory(pInTablePath, &nTableRows, &nTableCols, true, fptr);
	if (NULL == pTableData)
	{
		Message("Cannot GetCompositeTableIntoMemory() in GetTableDataPostExtractAndUpdateQuantiles", "ERROR");
		if (theSplineCounts) delete[] theSplineCounts;
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}


	//
	// Create a binary file image of the splined predictor data to pass to the matrix regression
	// There will always be the same number of rows in the table data as the predictor matrix
	// but the number of matrix cols may be less if not all predictors are 'In-use'.
	//
	GDM_INT nMatrixCols = 0;
	char *lpTmpFile = GetWorkspacePath(pModelParams);
	strcat(lpTmpFile, "\\gdmtmp.bin");

	if (!CreatePredictorBinary(lpTmpFile, pModelParams,
		                       pTableData, nTableRows, nTableCols, 
							   theSplineCounts, theQuantiles,
							   &nMatrixCols, false, fptr))
	{
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		if (pTableData) delete[] pTableData;
		return(false);
	}
	//PrintBinary2CSV(lpTmpFile, pModelParams, theSplineCounts, nTableRows, nMatrixCols, fptr);

	//
	// Create the memory for the Predictor Matrix
	//
	double *pPredictorMatrix = new double [nTableRows * nMatrixCols];
	if (NULL == pPredictorMatrix)
	{
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}


	//
	// Populate the Predictor Matrix from the binary image created in CreatePredictorBinary()
	//
	int h = _open( lpTmpFile, _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE );
	if (h < 0)
	{
		if (pPredictorMatrix) delete[] pPredictorMatrix;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}

	// incremental read
	double *pTmp = pPredictorMatrix;
	for ( int i=0; i<nMatrixCols; i++ )
	{
		_read(h, pTmp, (int)nTableRows * sizeof( double ));
		pTmp += nTableRows;
	}
	_close(h);


	//
	// remove the temporary matrix file if it exists (and it should!!) now that we don't need it
	//
	if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );
	if (lpTmpFile) delete[] lpTmpFile;
	

	//
	// Extract the response (observed) column 
	//
	fptr("Extracting Response Vector...", 25);
	double *pResponse = new double [nTableRows];	
	if (NULL == pResponse)
	{
		if (pTableData) delete[] pTableData;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		if (pPredictorMatrix) delete[] pPredictorMatrix;
		return(false);
	}
	memcpy(pResponse, &pTableData[COL_RESPONSE * nTableRows], nTableRows * sizeof(double));


	//
	// Extract the X0 column 
	//
	fptr("Extracting X0 Vector...", 35);
	double *pX0 = new double [nTableRows];
	if (NULL == pX0)
	{
		if (pTableData) delete[] pTableData;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		if (pPredictorMatrix) delete[] pPredictorMatrix;
		if (pResponse) delete[] pResponse;
		return(false);
	}
	memcpy(pX0, &pTableData[COL_SITE1_X0 * nTableRows], nTableRows * sizeof(double));


	//
	// Extract the Y0 column 
	//
	fptr("Extracting Y0 Vector...", 45);
	double *pY0 = new double [nTableRows];
	if (NULL == pY0)
	{
		if (pTableData) delete[] pTableData;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		if (pPredictorMatrix) delete[] pPredictorMatrix;
		if (pResponse) delete[] pResponse;
		if (pX0) delete[] pX0;
		return(false);
	}
	memcpy(pY0, &pTableData[COL_SITE1_Y0 * nTableRows], nTableRows * sizeof(double));


	//
	// Extract the X1 column 
	//
	fptr("Extracting X1 Vector...", 55);
	double *pX1 = new double [nTableRows];
	if (NULL == pX1)
	{
		if (pTableData) delete[] pTableData;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		if (pPredictorMatrix) delete[] pPredictorMatrix;
		if (pResponse) delete[] pResponse;
		if (pX0) delete[] pX0;
		if (pY0) delete[] pY0;
		return(false);
	}
	memcpy(pX1, &pTableData[COL_SITE2_X1 * nTableRows], nTableRows * sizeof(double));


	//
	// Extract the Y1 column 
	//
	fptr("Extracting Y1 Vector...", 65);
	double *pY1 = new double [nTableRows];
	if (NULL == pY1)
	{
		if (pTableData) delete[] pTableData;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		if (pPredictorMatrix) delete[] pPredictorMatrix;
		if (pResponse) delete[] pResponse;
		if (pX0) delete[] pX0;
		if (pY0) delete[] pY0;
		if (pX1) delete[] pX1;
		return(false);
	}
	memcpy(pY1, &pTableData[COL_SITE2_Y1 * nTableRows], nTableRows * sizeof(double));


	//
	// write prediction table
	//
	int numCoefficients = 0;
	double *pCoefficients = GetCoefficientsFromGDMParams(pModelParams, &numCoefficients );
	CreatePredictionTable(pOutPath, pResponse, pX0, pY0, pX1, pY1, pPredictorMatrix, pCoefficients, (int)nTableRows, (int)nMatrixCols, fptr);

	
	//
	// Cleanup
	//
	if (pTableData) delete[] pTableData;
	if (theSplineCounts) delete[] theSplineCounts; 
	if (theQuantiles) delete[] theQuantiles;
	if (pPredictorMatrix) delete[] pPredictorMatrix;
	if (pCoefficients) delete[] pCoefficients;
	if (pResponse) delete[] pResponse;
	if (pX0) delete[] pX0;
	if (pY0) delete[] pY0;
	if (pX1) delete[] pX1;
	if (pY1) delete[] pY1;
	return(true);
}



//
// Get the spline counts for all the predictors including geographic distance
// as a vector from the parameter file if they are "in-use"
//
int *GetInUseSplineCountVectorFromParams(char *pParams, int *pItems)
{
	// count the number of splines in the in-use predictors
	int nItems = 0;
	if (UseEuclidean(pParams))
	{
		++nItems;
	}

	int nPreds = GetNumPredictors(pParams);
	for (int i=1; i<=nPreds; i++ ) 
	{
		nItems += GetPredictorInUseAt(pParams, i); 
	}

	int *pSplineCounts = new int [nItems];

	if (UseEuclidean(pParams))
	{
		pSplineCounts[0] = GetEuclideanSplines(pParams);
		int nThisPred = 1;
		for (int i=1; i<=nPreds; i++ ) 
		{
			if (1 == GetPredictorInUseAt(pParams, i))
			{
				pSplineCounts[nThisPred] = GetPredictorSplinesAt(pParams, i);
				++nThisPred;
			}
		}
	}
	else
	{
		int nThisPred = 0;
		for (int i=1; i<=nPreds; i++ ) 
		{
			if (1 == GetPredictorInUseAt(pParams, i))
			{
				pSplineCounts[nThisPred] = GetPredictorSplinesAt(pParams, i);
				++nThisPred;
			}
		}
	}

	*pItems = nItems;
	return(pSplineCounts);
}



//
// Get the quantiles for all the predictors including geographic distance as a vector if they are "in-use"
//
double *GetInUseQuantilesFromParams(char *pParams, int *pItems)
{
	// count the total number of splines for the "in-Use" predictors in the model
	int nSplines = 0;
	if (UseEuclidean(pParams))
	{
		nSplines += GetEuclideanSplines(pParams);
	}

	int nPreds = GetNumPredictors(pParams);
	for (int i=1; i<=nPreds; i++ ) 
	{
		if (1 == GetPredictorInUseAt(pParams, i) )
		{
			nSplines += GetPredictorSplinesAt(pParams, i); 
		}
	}
	
	// allocate and populate the Quantiles vector
	double *pQuantiles = new double [nSplines];

	if (UseEuclidean(pParams))
	{
		int nThis = 0;
		for (int i=1; i<=GetEuclideanSplines(pParams); i++ )
		{
			pQuantiles[nThis++] = GetEuclideanSplineAt(pParams, i);
		}

		for ( int i=1; i<=nPreds; i++ )
		{
			if (1 == GetPredictorInUseAt(pParams, i) )
			{
				for ( int j=1; j<=GetPredictorSplinesAt(pParams, i); j++ )
				{
					pQuantiles[nThis++] = GetPredictorSplineAt(pParams, i, j);
				}
			}
		}
	}
	else
	{
		int nThis = 0;
		for ( int i=1; i<=nPreds; i++ )
		{
			if (1 == GetPredictorInUseAt(pParams, i) )
			{
				for ( int j=1; j<=GetPredictorSplinesAt(pParams, i); j++ )
				{
					pQuantiles[nThis++] = GetPredictorSplineAt(pParams, i, j);
				}
			}
		}
	}

	*pItems = nSplines;
	return(pQuantiles);
}



//
// Extract the coefficients from the "in-Use" predictors from a GDM Model param file
//
double *GetCoefficientsFromGDMParams(char *pModelParams, int *numCoefficients )
{
	// count the total number of splines for the "in-Use" predictors in the model
	int nCoeffs = 0;

	if (UseEuclidean(pModelParams))
	{
		nCoeffs += GetEuclideanSplines(pModelParams);
	}

	int nPreds = GetNumPredictors(pModelParams);
	for (int i=1; i<=nPreds; i++ ) 
	{
		if (1 == GetPredictorInUseAt(pModelParams, i) )
		{
			nCoeffs += GetPredictorSplinesAt(pModelParams, i); 
		}
	}
	
	// add a column for the intercept [0] and allocate and populate the Coefficients vector
	++nCoeffs; 

	double *pCoefficients = new double [nCoeffs];
	pCoefficients[0] = GetIntercept(pModelParams);

	if (UseEuclidean(pModelParams))
	{
		int nThis = 1;
		for (int i=1; i<=GetEuclideanSplines(pModelParams); i++ )
		{
			pCoefficients[nThis++] = GetEuclideanCoeffAt(pModelParams, i);
		}

		for ( int i=1; i<=nPreds; i++ )
		{
			if (1 == GetPredictorInUseAt(pModelParams, i) )
			{
				for ( int j=1; j<=GetPredictorSplinesAt(pModelParams, i); j++ )
				{
					pCoefficients[nThis++] = GetPredictorCoeffAt(pModelParams, i, j);
				}
			}
		}
	}
	else
	{
		int nThis = 1;
		for ( int i=1; i<=nPreds; i++ )
		{
			if (1 == GetPredictorInUseAt(pModelParams, i) )
			{
				for ( int j=1; j<=GetPredictorSplinesAt(pModelParams, i); j++ )
				{
					pCoefficients[nThis++] = GetPredictorCoeffAt(pModelParams, i, j);
				}
			}
		}
	}

	*numCoefficients = nCoeffs;
	return(pCoefficients);
}



//
// Write the observed and predicted data to comma delimited textfile and to a pair of binary files for plotting
//
void CreatePredictionTable(char *pOutTablePath, double *pResponse, 
	                       double *pX0, double *pY0, double *pX1, double *pY1,
						   double *pPredData, double *pCoefficients, 
						   int nRows, int nCols, FPTR fptr)
{
	int nCurrent = 0;
	fptr("Writing predictions to file...", nCurrent );
	FILE *fp = fopen( pOutTablePath, "w+t" );

	fprintf( fp, "X0,Y0,X1,Y1,Observed,Predicted,Linked,Residual,Pearsons\n");
	for ( int i=0; i<nRows; i++ )
	{
		if ( (long long)(i) * 100 / nRows > nCurrent )
		{
			nCurrent = (long long)(i) * 100 / nRows;
			fptr("Writing predictions to file...", nCurrent );
		}

		// write site pair data
		fprintf( fp, "%lf,%lf,%lf,%lf,", pX0[i], pY0[i], pX1[i], pY1[i]);

		// write Observed response
		double dObserved = pResponse[i];
		fprintf( fp, "%lf,", dObserved);

		// write Predicted response
		double dPredicted = CalcDissimilarity( pPredData, pCoefficients, nRows, nCols, i );
		fprintf( fp, "%lf,", dPredicted);

		// write Linked predicted response
		double dLinked = 1.0 - exp(-dPredicted);
		fprintf( fp, "%lf,", dLinked);

		// write standard residual
		double dResidual = dObserved - dLinked;
		fprintf( fp, "%lf,", dResidual );

		// write Pearson residial
		fprintf( fp, "%lf\n", dResidual / sqrt(dLinked * (1.0 - dLinked)) );
	}

	fclose(fp);
	nCurrent = 0;
	fptr("Writing predictions to file...", nCurrent );
}


//
// Determine the maximum number of contiguous doubles (8 bytes each) that can be allocated
//
double GetMaxMemBlock(FPTR fptr)
{
	GDM_UINT MaxVal = 1000000;
	double *p = NULL;
	int nCurrent = 0;
	int NumIts = 0;
	while(1)
	{
		if (nCurrent < 99)
			++nCurrent;
		else
			nCurrent = 1;
		fptr("Determining Maximum Memory Block...", nCurrent);
		
		try
		{

			p = new double [MaxVal];
			if (p) 
			{
				delete[] p;
				MaxVal += 1000000;
			}

			else
			{
				MaxVal -= 1000000;
				break;
			}
		}

		catch(...)
		{
			//Message("Got an exception thrown here", "INFO");
			break;
		}
	} 
	fptr("Status: Ready", 0);
	return((double)MaxVal);
}


//
// Write the GDM results back to the parameter file
//
void UpdateParameterFile(char *pParams, 
		                double dNullDeviance, double dGDMDeviance, double dDevianceExplained, 
						double dIntercept, double *pCoefficients, 
						int nTotalSplines, int *pSplineCounts, int nRows, FPTR fptr)
{
	if (fptr != NULL)
		fptr("Updating parameter file...", 0);
	SetNullDeviance(pParams, dNullDeviance);
	SetGDMDeviance(pParams, dGDMDeviance);
	SetDevianceExplained(pParams, dDevianceExplained);
	SetIntercept(pParams, dIntercept);

	if (fptr != NULL)
		fptr("Updating parameter file...", 20);
	SetSumOfCoefficients( pParams, CalcCoefficientSum(pCoefficients, nTotalSplines) );	
	SetNumberOfActivePredictors( pParams, CalcNumberOfActivePredictors(pParams, pCoefficients, pSplineCounts) );
	SetNumberCoefficientsAboveZero( pParams, CalcNumberofCoeffGreaterThanZero(pCoefficients, nTotalSplines) );
	SetNumberOfSitePairs( pParams, nRows );

	double *pCoeff = pCoefficients;
	int    *pCount = pSplineCounts;

	if (fptr != NULL)
		fptr("Updating parameter file...", 40);
	int nSplines = *pCount;
	if (UseEuclidean(pParams))
	{
		for ( int i=0; i<nSplines; i++ )
		{
			SetEuclideanCoeffAt(pParams, i+1, pCoeff[i]);
		}

		// increment the pointers...
		pCoeff += nSplines;
		pCount += 1;
	}

	else // set the euclidean quantiles and coefficients to zero as they were NOT used in the GDM
	{
		for ( int i=0; i<nSplines; i++ )
		{
			SetEuclideanCoeffAt(pParams, i+1, 0.0);
		}
	}
	

	if (fptr != NULL)
		fptr("Updating parameter file...", 50);
	int nCurrent = 50;
	int nPreds = GetNumPredictors(pParams);
	// now do the environmental predictors
	for (int p=1; p<=nPreds; p++ )
	{
		if (((p-1) * 50 / nPreds)+50 > nCurrent)
		{
			nCurrent = ((p-1) * 50 / nPreds)+50;
			if (fptr != NULL)
				fptr("Updating parameter file...", nCurrent);
		}

		nSplines = *pCount;
		if (GetPredictorInUseAt(pParams, p)) 
		{
			for ( int i=0; i<nSplines; i++ )
			{
				SetPredictorCoeffAt(pParams, p, i+1, pCoeff[i]);
			}
			// increment the coefficient pointer...
			pCoeff += nSplines;			
		}
		else
		{
			for ( int i=0; i<nSplines; i++ )
			{
				SetPredictorCoeffAt(pParams, p, i+1, 0.0);
			}
		}		
		// increment the spline count pointer...
		pCount += 1;
	}
	if (fptr != NULL)
		fptr("Updating parameter file...", 0);
}
