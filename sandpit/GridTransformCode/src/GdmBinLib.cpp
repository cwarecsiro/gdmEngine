//
// GdmBinLib.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>

// for the Binary File IO
#include <io.h>
#include <fcntl.h>
#include <sys/stat.h>

// support fucntions
#include "GdmBinLib.h"
//#include "TableFuncs.h"
#include "Message.h"
#include "myCallback.h"
//#include "gmPath.h"
#include "NNLS_Double.h"
//#include "GdmModel.h"


//
// Local function Declarations
//
bool CreateBinaryFile(GdmModel *gdmModel, FPTR fptr);
bool RunGDMRegression(GdmModel *gdmModel, FPTR fptr);
double CalculatePredictedResponseDouble(double *pPredictorData, int nRows, int nCols, double *pCoeffs, int Index );
bool CreateResponseBinaryFile(GdmModel *gdmModel, FPTR fptr);

//
// The main function called from the External Interface
//
bool RunGDMFromTable(char *lpParamFilePath, FPTR fptr)
{
	//
	// The GdmModel class holds all the spline and coefficient information
	// and is updated and used by CreateBinaryFile() and RunGDMRegression()
	//
	GdmModel *gdmModel = new GdmModel(lpParamFilePath);

	
	//
	// Create the Binary File representation of the 
	// I-Splined Predictor Data from the composite data table.
	//
	fptr( "Creating the Binary File...", 0 );
	if ( false == CreateBinaryFile(gdmModel, fptr) )
	{
		Message("RunGDMFromTable() - Cannot CreateBinaryFile", "ERROR");
		if (gdmModel) delete gdmModel;
		return(false);
	}


	//
	// Create the binary representation of the reponse data
	//
	fptr( "Creating the Response Binary File...", 0 );
	if ( false == CreateResponseBinaryFile(gdmModel, fptr) )
	{
		Message("RunGDMFromTable() - Cannot CreateResponseBinaryFile", "ERROR");
		if (gdmModel) delete gdmModel;
		return(false);
	}


	//
	// Run the GDM Regression.
	//
	fptr("Running the GDM Regression...", 0);
	if ( false == RunGDMRegression(gdmModel, fptr) )
	{
		Message("RunGDMFromTable() - Cannot RunGDMRegression", "ERROR");
		if (gdmModel) delete gdmModel;
		return(false);
	}
	

	//
	// Update the Parameter file and write the site pair statistics
	//
	fptr("Writing the GDM Results...", 0);
	gdmModel->Write(lpParamFilePath);
	gdmModel->WriteGdmSiteStats();


	//
	// Clean up
	//
	if (gdmModel) delete gdmModel;
	return(true);
}



//
// Create the Core GDM Binary File.
// This contains the Intercept data of ONES, 
// the optional Geographic predictor splines
// and the I-Splined Predictor data from the composite data table 
//
bool CreateBinaryFile(GdmModel *gdmModel, FPTR fptr)
{
	//
	// Create the binary file
	//
	int h = _open( gdmModel->GetBinaryFilePath(), _O_BINARY | _O_CREAT | O_RDWR, S_IREAD | S_IWRITE  );
	if ( h < 0 )
	{
		Message( "<CreateBinaryFile()> Cannot open binary file for READ/WRITE", "ERROR" );
		return(false);
	}

	int nRows = GetNumRowsFromTableFile(gdmModel->GetInputDataPath());
	int nPreds = GetNumPredictorsFromHeader(gdmModel->GetInputDataPath());


	// write the intercept data of ONES
	double dIntercept = 1.0;
	for ( int i=0; i<nRows; i++ )
	{
		_write(h, &dIntercept, sizeof(dIntercept));
	}
	fptr( "Running GDM From Table...", 10 );


	// write the optional geographic predictor data
	if ( gdmModel->UseGeoPred() )
	{
		double *pGeoData = GetGeographicISplines(nRows, gdmModel, gdmModel->GetInputDataPath());
		_write(h, pGeoData, nRows * gdmModel->GetGeoPredictor()->GetSplines() * sizeof(double));
		if ( pGeoData ) delete[] pGeoData;
	}
	fptr( "Running GDM From Table...", 40 );


	// write the I-Splined Predictor Data
	int nCurrent = 40;
	for ( int i=0; i<nPreds; i++ )
	{
		if ( (i * 60 / nPreds) + 40 > nCurrent )
		{
			nCurrent = (i * 60 / nPreds) + 40;
			fptr( "Running GDM From Table...", nCurrent );
		}

		if ( gdmModel->GetPredictorAt(i)->IsInUse() )
		{
			// convert predictor index from zero-base to one-base and get I-Spline data for this predictor		
			double *pPredData = GetPredictorISplines(nRows, gdmModel, gdmModel->GetInputDataPath(), i);
			_write(h, pPredData, nRows * gdmModel->GetPredictorAt(i)->GetSplines() * sizeof(double));
			if ( pPredData ) delete[] pPredData;
		}
	}
	fptr( "Running GDM From Table...", 100 );


	//
	// cleanup
	//
	_close(h);
	return(true);
}



//
// Run the GDM matrix regression
// and update the coefficients and
// GDM model data in the GdmModel Class
//
bool RunGDMRegression(GdmModel *gdmModel, FPTR fptr)
{
	double fNULLDeviance = 0.0;
	double fGDMDeviance = 0.0;
	double fExplainedDeviance = 0.0;

	int nRows = gdmModel->GetNumRows();
	int nPreds = GetNumPredictorsFromHeader(gdmModel->GetInputDataPath());
	int nCols = gdmModel->GetNumCols();

	double *pWeights = new double [nRows];
	for ( int i=0; i<nRows; i++ ) pWeights[i] = 1.0;
	double *pCoeffs = new double [nCols];
	double *pResponse = GetColumnAt(gdmModel->GetInputDataPath(), nRows, 0);

	//
	// Do the regression
	//
	NNLS_Double( gdmModel->GetBinaryFilePath(), nRows, nCols, 
				 pResponse, pWeights, pCoeffs, 
				 &fNULLDeviance, &fGDMDeviance, &fExplainedDeviance,
				 fptr);

	//
	// Apply the coefficients to the spline data to set the Predicted Response Values
	//
	double *pData = new double [nRows * nCols];
	if ( NULL == pData )
	{
		Message("RunGDMRegression() - Cannot allocate pData", "ERROR");
		return(false);
	}
	int h = _open( gdmModel->GetBinaryFilePath(), _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE  );
	if ( h < 0 )
	{
		Message( "RunGDMRegression() - Cannot open Binary file for READ", "ERROR" );
		if (pData) delete[] pData;
		return(false);
	}
	_read( h, pData, nRows * nCols * sizeof( double ) );
	_close( h );
	for ( int i=0; i<nRows; i++ )
	{
		gdmModel->GetPredicted()[i] = CalculatePredictedResponseDouble(pData, nRows, nCols, pCoeffs, i );
	}
	if (pData) delete[] pData;


	//
	// Copy Model Coefficients into the GdmModel Class
	//
	gdmModel->SetIntercept(pCoeffs[0]);
	gdmModel->SetNULLDeviance(fNULLDeviance);
	gdmModel->SetGDMDeviance(fGDMDeviance);
	gdmModel->SetExplainedDeviance(fExplainedDeviance);


	//
	// Copy Geographic Coefficients into the GdmModel Class
	//
	int nThisPred = 1;
	if ( gdmModel->UseGeoPred() )
	{
		for ( int i=0; i<gdmModel->GetGeoPredictor()->GetSplines(); i++ )
		{
			gdmModel->GetGeoPredictor()->SetCoeffAt(pCoeffs[i+1],i);
			++nThisPred;
		}
	}


	//
	// Copy Predictor Coefficients into the GdmModel Class
	//	
	for ( int i=0; i<nPreds; i++ )
	{
		if ( gdmModel->GetPredictorAt(i)->IsInUse() )
		{
			int nQuants = gdmModel->GetPredictorAt(0)->GetSplines();
			for ( int j=0; j<nQuants; j++ )
			{
				gdmModel->GetPredictorAt(i)->SetCoeffAt(pCoeffs[nThisPred], j);
				++nThisPred;
			}
		}
	}


	//
	// clean up
	//
	if (pResponse) delete[] pResponse;
	if (pWeights) delete[] pWeights;
	if (pCoeffs) delete[] pCoeffs;
	return(true);
}



//
// Calculate predicted Response as a double
//
double CalculatePredictedResponseDouble(double *pPredictorData, int nRows, int nCols, double *pCoeffs, int Index )
{
	double fVal = 0.0;

	// sum the terms
	for ( int i=0; i<nCols; i++ )
	{
		fVal += pCoeffs[i] * pPredictorData[Index + (i * nRows)];
	}

	return( fVal );
}



//
// Write a binary file version of the response data
//
bool CreateResponseBinaryFile(GdmModel *gdmModel, FPTR fptr)
{
	int nRows = gdmModel->GetNumRows();
	int nCols = gdmModel->GetNumCols();
	double *pResponse = GetColumnAt(gdmModel->GetInputDataPath(), nRows, 0);

	Message(nRows, "nRows in CreateResponseBinaryFile");
	Message(nCols, "nCols in CreateResponseBinaryFile");
	
	//
	// Create the binary file
	//
	char buff[512];
	sprintf(buff, "%s\\response.bin", gdmModel->GetWorkspacePath() );
	int h = _open( buff, _O_BINARY | _O_CREAT | O_RDWR, S_IREAD | S_IWRITE  );
	if ( h < 0 )
	{
		Message( "<CreateResponseBinaryFile()> Cannot open binary file for READ/WRITE", "ERROR" );
		return(false);
	}

	_write(h,pResponse, nRows * sizeof(double));
	_close(h);

	if (pResponse) delete[] pResponse;
	return(true);
}


