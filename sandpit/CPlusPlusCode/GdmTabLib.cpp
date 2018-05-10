//
// GdmTabLib
//
#include "stdafx.h"
#include "GdmBinLib.h"
#include "GdmTabLib.h"
#include "NNLS_Double.h"
#include "Message.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <io.h>
#include <fcntl.h>
#include <sys/stat.h>


#if defined _M_X64

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Functions called from R -Stats ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Main GDM fitting functions called from R
//
//
// Note that pData represents a column major matrix
//
void GDM_FitFromTable(char **wspath, 
				      double *pData, 
	                  int *pDoGeo, int *pPreds, 
				      int *pRows, int *pCols, 
				      int *pSplines, double *pQuantiles,
				      double *pGDMDev, double *pNullDev, double *pExpDev, 
				      double *pIntercept, double *pCoeffs,
				      double *pY, // observed
				      double *pX, // predicted
				      double *pE) // ecological dist
{
	//
	// Everything passed from R needs to be dereferenced...
	//
	int nDoGeo = *pDoGeo;
	int nPreds = *pPreds;
	long long nRows = *pRows;
	int nCols = *pCols;


	//
	// Get the response column
	//
	double *pResponse = &pData[COL_RESPONSE * nRows];
	if ( NULL == pResponse ) 
	{
		Message("pResponse is NULL", "ERROR in GDM_FitFromTable");
		return;
	}


	//
	// Get the weights column
	//
	double *pWeights = &pData[COL_WEIGHTS * nRows];
	if ( NULL == pWeights ) 
	{
		Message("pWeights is NULL", "ERROR in GDM_FitFromTable");
		return;
	}


	//
	// Get the predictor matrix
	//
	double *pPredData = ConstructMatrix(nDoGeo, pData, pQuantiles, nPreds, pSplines, nRows);
	if ( NULL == pPredData )
	{
		Message("pPredData is NULL", "ERROR in GDM_FitFromTable");
		return;
	}
	// Sum the splines vector for the total number of splines
	int nTotalSplines = GetTotalSplineCount(pSplines, nPreds);
	// Message(nTotalSplines, "nTotalSplines");
	// pQuantiles will have a length of nPreds * nTotalSplines 
	// ShowQuantiles(pQuantiles, nPreds, pSplines);
	//DebugPredMatrix("PredMatrix.csv", pPredData, nRows, nPreds, pSplines, nTotalSplines+1);

	
	//
	// Write a binary file image of the predictor matrix
	//
	char lpTmpFile[256];		
	sprintf(lpTmpFile, "%s/%s", *wspath, "gdmtmp.bin" );
	int h = _open( lpTmpFile, _O_BINARY | _O_CREAT | _O_TRUNC | _O_RDWR, S_IREAD | S_IWRITE );
	if ( h < 0 )
	{
		Message("Cannot create binary image file", "ERROR in GDM_FitFromTable");
		if (pPredData) delete[] pPredData;
		return;
	}	
	long long nThis = 0;
	for ( int i=0; i<(nTotalSplines+1); i++ )
	{
		_write( h, &pPredData[nThis], int(nRows * sizeof( double )) );
		nThis += nRows;
	}
	_close( h );


		
	//
	// Do the matrix regression
	//
	double dGDMDeviance;
	double *pCoefficients = WeightedNNLSRegression( lpTmpFile, 
		                                            pPredData, 
											        nRows, 
											        nTotalSplines+1, 
											        pResponse, 
											        &dGDMDeviance, 
											        pWeights,
													NULL);

	//
	// remove the temporary matrix file if it exists (and it should!!) now that we don't need it
	//
	if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );


	//
	// create a NULL model and return the deviance
	//
	double dNullDeviance = GetWeightedNULLDeviance( nRows, pResponse, pWeights );


	//
	// calculate deviance explained as a percentage
	//
	double dDevianceExplained = ( 1.0 - ( dGDMDeviance / dNullDeviance ) ) * 100;


	//
	// display results
	//
	pGDMDev[0] = dGDMDeviance;
	pNullDev[0] = dNullDeviance;
	pExpDev[0] = dDevianceExplained;
	pIntercept[0] = pCoefficients[0];
	for ( int i=1; i<nTotalSplines+1; i++ )
	{
		pCoeffs[i-1] = pCoefficients[i];
	}


	//
	// copy to response data to a vector
	//
	for ( long long i=0; i<nRows; i++ )
	{
		pY[i] = pResponse[i];
		pE[i] = CalcDissimilarity( pPredData, pCoefficients, nRows, nTotalSplines+1, i );
		pX[i] = 1.0 - exp(-pE[i]);
	}


	//
	// Clean up
	//
	if (pPredData) delete[] pPredData;
	if ( pCoefficients ) delete[] pCoefficients;
}





//
// Main GDM predict function called from R
//
//
// Note that pData represents a column major matrix
//
void GDM_PredictFromTable(double *pData, 
		                  int *pDoGeo, int *pPreds, int *pRows, 
					      double *pQuantiles, int *pSplines, double *pCoeffs,
					      double *pX) // predicted
{
	//
	// Everything passed from R needs to be dereferenced...
	//
	int nDoGeo = *pDoGeo;
	int nPreds = *pPreds;
	long long nRows  = (long long)(*pRows);
		

	//
	// Get the predictor matrix
	//
	double *pPredData = ConstructMatrix(nDoGeo, pData, pQuantiles, nPreds, pSplines, nRows);
	if ( NULL == pPredData )
	{
		Message("pPredData is NULL", "ERROR");
		return;
	}

	int nTotalSplines = GetTotalSplineCount(pSplines, nPreds);
	//Message(nTotalSplines, "nTotalSplines");
	// pQuantiles will have a length of nPreds * nTotalSplines 
	//ShowQuantiles(pQuantiles, nPreds, pSplines);
	//DebugPredMatrix("PredMatrix.csv", pPredData, nRows, nPreds, pSplines, nTotalSplines+1);


	//
	// copy predicted data to pX
	//
	for ( long long i=0; i<nRows; i++ )
	{
		pX[i] = 1.0 - exp(-CalcDissimilarity( pPredData, pCoeffs, nRows, nTotalSplines+1, i ));
	}


	//
	// Clean up
	//
	if (pPredData) delete[] pPredData;
}




//
// Populate pPredData as a transformed GDM predictor to plot in R-Stats
//
void GetPredictorPlotData( double *pPredData, int *pLength,
						   double *pCoefficients,
						   double *pQuantiles,
						   int *pSplines )
{
	int nLen  = *pLength;
	int nSplines = *pSplines;

	double dInc = fabs(pQuantiles[nSplines-1]-pQuantiles[0]) / nLen;
	double dVal = pQuantiles[0];
	for ( int i=0; i<nLen; i++ )
	{
		pPredData[i] = 0;
		for ( int j=0; j<nSplines; j++ )
		{
			if ( j==0 )
			{
				pPredData[i] += pCoefficients[j] * DoSplineCalc( dVal, pQuantiles[j], pQuantiles[j], pQuantiles[j+1] );
			}

			else if ( j == nSplines-1 )
			{
				pPredData[i] += pCoefficients[j] * DoSplineCalc( dVal, pQuantiles[j-1], pQuantiles[j], pQuantiles[j] );
			}

			else
			{
				pPredData[i] += pCoefficients[j] * DoSplineCalc( dVal, pQuantiles[j-1], pQuantiles[j], pQuantiles[j+1] );
			}
		}
		dVal += dInc;
	}
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Local Support Functions ///////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Predict dissimilarity from spline data and coefficients
//
double CalcDissimilarity( double *pData, double *pCoeffs, long long nRows, int nCols, long long nIndex )
{
	double dVal = 0.0;

	double *pTmp = &pData[nIndex];
	for ( int i=0; i<nCols; i++ )
	{
		dVal += *pTmp * pCoeffs[i];
		pTmp += nRows;
	}
	return(dVal);
}


//
// calculate total number of splines
//
int GetTotalSplineCount(int *pSplines, int nPreds)
{
	int nCount = 0;
	for ( int i=0; i<nPreds; i++ )
	{
		nCount += pSplines[i];
	}
	return(nCount);
}


//
// Construct a regression matrix for predictors
//
double *ConstructMatrix(int nDoGeo, double *pData, double *pQuants, int nPreds, int *pSplines, long long nRows)
{
	int nTotalSplines = GetTotalSplineCount(pSplines, nPreds);
	//Message(nTotalSplines, "nTotalSplines");

	//
	// Construct the predictor matrix with an extra column for the intercept
	//
	long long nSize = nRows * (nTotalSplines + 1);
	double *pPredData = new double [nSize];
	if ( NULL == pPredData )
	{
		Message("Cannot allocate Predictor Data", "ERROR in ConstructMatrix");
		return(NULL);
	}
	for ( long long i=0; i<nSize; i++ ) pPredData[i] = 0.0;


	//
	// Initialise the intercept column
	//
	double *pLoc = &pPredData[0];
	for ( long long i=0; i<nRows; i++ ) pLoc[i] = 1.0;
	pLoc += nRows;   // start at the first column after the intercept


	//
	// Initialise for the geographic distance if we are using it
	//
	if (nDoGeo)
	{
		//
		// Initialise the matrix for the geographic followed by the environmental predictors
		//		
		double dMin,dMid,dMax,dVal1,dVal2;
		int CurrentSpline = 0;

		for ( int pred=0; pred<nPreds; pred++ )
		{
			int d1Offset = ((LEADING_COLS+pred-1)*(int)nRows);
			int d2Offset = ((LEADING_COLS+pred+nPreds-2)*(int)nRows);

			int nThis = 0;
			for ( long long i=0; i<nRows; i++ )
			{
				if ( pred == 0 )
				{
					double distX = fabs( pData[(COL_SITE1_X0 * nRows)+i] - pData[(COL_SITE2_X1 * nRows)+i] );
					double distY = fabs( pData[(COL_SITE1_Y0 * nRows)+i] - pData[(COL_SITE2_Y1 * nRows)+i] );

					dVal1 = 0.0;
					dVal2 = sqrt( ( distX * distX ) + ( distY * distY ) ); 							
				}
				else
				{
					dVal1 = pData[d1Offset+i];
					dVal2 = pData[d2Offset+i];
				}

				for ( int spl = 0; spl<pSplines[pred]; spl++ )
				{
					if ( spl == 0 )                              // first spline
					{
						dMin = pQuants[CurrentSpline+spl+0];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1]; 				
					}
		
					else if ( spl == (pSplines[pred]-1) )        // last spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+0];
					}

					else	                                     // a middle spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1];
					}

					// calculate the spline values for each site
					double d1a = DoSplineCalc( dVal1, dMin, dMid, dMax );
					double d2a = DoSplineCalc( dVal2, dMin, dMid, dMax );

					// set the distance between sites for this spline
					pLoc[ (spl * nRows) + nThis ] = fabs( d2a - d1a );	

				} // for ( int spl = 0; spl<pSplines[pred]; spl++ )

				// increment for the next row
				++nThis;

			} // for ( long long i=0; i<nRows; i++ )

			// increment pointer to the start of the next predictor's spline columns
			pLoc += pSplines[pred] * nRows;

			// increment the current quantile index for the next pred
			CurrentSpline += pSplines[pred];

		} // for ( int pred=0; pred<nPreds; pred++ )
	}


	else
	{
		//
		// Initialise the matrix for the environmental predictors only
		//
		double dMin,dMid,dMax;
		int CurrentSpline = 0;

		for ( int pred=0; pred<nPreds; pred++ )
		{
			int d1Offset = ((LEADING_COLS+pred)*(int)nRows);
			int d2Offset = ((LEADING_COLS+nPreds+pred)*(int)nRows);

			int nThis = 0;
			for ( long long i=0; i<nRows; i++ )
			{
				double dVal1 = pData[d1Offset+i];
				double dVal2 = pData[d2Offset+i];

				for ( int spl = 0; spl<pSplines[pred]; spl++ )
				{
					if ( spl == 0 )                              // first spline
					{
						dMin = pQuants[CurrentSpline+spl+0];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1]; 
					
					}
		
					else if ( spl == (pSplines[pred]-1) )        // last spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+0];
					}

					else	                                     // a middle spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1];
					}


					// calculate the spline values for each site
					double d1a = DoSplineCalc( dVal1, dMin, dMid, dMax );
					double d2a = DoSplineCalc( dVal2, dMin, dMid, dMax );

					// set the distance between sites for this spline
					pLoc[ (spl * nRows) + nThis ] = fabs( d2a - d1a );	

				} // for ( int spl = 0; spl<pSplines[pred]; spl++ )

				// increment for the next row
				++nThis;

			} // for ( int i=0; i<nRows; i++ )

			// increment pointer to the start of the next predictor's spline columns
			pLoc += pSplines[pred] * nRows;

			// increment the current quantile index for the next pred
			CurrentSpline += pSplines[pred];

		} // for ( int pred=0; pred<nPreds; pred++ )

	} // environmental predictors only

	return(pPredData);
}




//
// Calculate the I-Spline value for dVal given quantiles q1, q2, q3
//
double DoSplineCalc( double dVal, double q1, double q2, double q3 )
{
	if ( dVal <= q1 ) return(0.0);

	else if ( dVal >= q3 ) return( 1.0 );

	else if ( ( q1 < dVal ) && ( dVal < q2 ) )
		return( ( ( ( dVal - q1 ) * ( dVal - q1 ) ) / ( ( q2 - q1 ) * ( q3 - q1 ) ) ) );

	else
		return( ( 1.0 - ( ( ( q3 - dVal ) * ( q3 - dVal ) ) / ( ( q3 - q2 ) * ( q3 - q1) ) ) ) );
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Local Debugging Functions /////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Display the contents of the quantile vector
//
void ShowQuantiles(double *pQuants, int nPreds, int *pSplines)
{
	char buff[256];
	double *pTmp = &pQuants[0];
	for ( int i=0; i<nPreds; i++ )
	{
		sprintf( buff, "Quant %d: ", i+1 );

		for ( int j=0; j<pSplines[i]; j++ )
		{
			sprintf( buff, "%s %lf ", buff, *pTmp);			
			++pTmp;
		}
		Message(buff);
	}
}




//
// Write the Predictor Matrix to a comma delimited file
//
void DebugPredMatrix(char *pPath, double *pPredData, long long nRows, int nPreds, int *pSplines, int nCols)
{
	FILE *fp = fopen(pPath, "w+t");

	// write the header
	fprintf( fp, "Intercept,");			
	for ( int p=0; p<nPreds; p++ )
	{
		for ( int s=0; s<pSplines[p]; s++ )
		{
			fprintf( fp, "Pred%dSpline%d", p+1,s+1);			
			if ( s < pSplines[p]-1 )
				fprintf( fp, "," );
		}
		if ( p < nPreds-1 )
			fprintf( fp, "," );
		else
			fprintf( fp, "\n" );

	}
	
	for ( long long i=0; i<nRows; i++ )
	{
		for ( int j=0; j<nCols; j++ )
		{
			fprintf( fp, "%lf", pPredData[(j*nRows)+i]);

			if ( j < nCols-1 ) 
				fprintf(fp, "," );
			else
				fprintf(fp, "\n" );
		}
	}

	if (fp) fclose(fp);
}




#elif defined _WIN32


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Functions called from R -Stats ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Main GDM fitting functions called from R
//
//
// Note that pData represents a column major matrix
//
void GDM_FitFromTable(char **wspath, 
				      double *pData, 
	                  int *pDoGeo, int *pPreds, 
				      int *pRows, int *pCols, 
				      int *pSplines, double *pQuantiles,
				      double *pGDMDev, double *pNullDev, double *pExpDev, 
				      double *pIntercept, double *pCoeffs,
				      double *pY, // observed
				      double *pX, // predicted
				      double *pE) // ecological dist
{
	//Message("At Entry to GDM_FitFromTable", "INFO");


	//
	// Everything passed from R needs to be dereferenced...
	//
	int nDoGeo = *pDoGeo;
	int nPreds = *pPreds;
	int nRows = *pRows;
	int nCols = *pCols;


	//
	// Get the response column
	//
	double *pResponse = &pData[COL_RESPONSE * nRows];
	if ( NULL == pResponse ) 
	{
		Message("pResponse is NULL", "ERROR in GDM_FitFromTable");
		return;
	}


	//
	// Get the weights column
	//
	double *pWeights = &pData[COL_WEIGHTS * nRows];
	if ( NULL == pWeights ) 
	{
		Message("pWeights is NULL", "ERROR in GDM_FitFromTable");
		return;
	}


	//
	// Get the predictor matrix
	//
	double *pPredData = ConstructMatrix(nDoGeo, pData, pQuantiles, nPreds, pSplines, nRows);
	if ( NULL == pPredData )
	{
		Message("pPredData is NULL", "ERROR in GDM_FitFromTable");
		return;
	}
	// Sum the splines vector for the total number of splines
	int nTotalSplines = GetTotalSplineCount(pSplines, nPreds);
	// Message(nTotalSplines, "nTotalSplines");
	// pQuantiles will have a length of nPreds * nTotalSplines 
	//ShowQuantiles(pQuantiles, nPreds, pSplines);
	DebugPredMatrix("D:\\GDM_2014\\GDM_Table_version_For_R\\GDM_BUG_11092014\\PredMatrix.csv", 
		            pPredData, nRows, nPreds, pSplines, nTotalSplines+1);

	
	//
	// Write a binary file image of the predictor matrix
	//
	char lpTmpFile[256];		
	sprintf(lpTmpFile, "%s/%s", *wspath, "gdmtmp.bin" );
	int h = _open( lpTmpFile, _O_BINARY | _O_CREAT | _O_TRUNC | _O_RDWR, S_IREAD | S_IWRITE );
	if ( h < 0 )
	{
		Message("Cannot create binary image file", "ERROR in GDM_FitFromTable");
		if (pPredData) delete[] pPredData;
		return;
	}	
	int nSize = nRows * (nTotalSplines+1);
	_write( h, pPredData, nSize * sizeof( double ) );
	_close( h );


		
	//
	// Do the matrix regression
	//
	double dGDMDeviance;
	double *pCoefficients = WeightedNNLSRegression( lpTmpFile, 
		                                            pPredData, 
											        nRows, 
											        nTotalSplines+1, 
											        pResponse, 
											        &dGDMDeviance, 
											        pWeights,
													NULL);
	if (NULL == pCoefficients)
	{
		Message("Got a NULL result from WeightedNNLSRegression()", "ERROR");
	}

	//
	// remove the temporary matrix file if it exists (and it should!!) now that we don't need it
	//
	if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );


	//
	// create a NULL model and return the deviance
	//
	double dNullDeviance = GetWeightedNULLDeviance( nRows, pResponse, pWeights );


	//
	// calculate deviance explained as a percentage
	//
	double dDevianceExplained = ( 1.0 - ( dGDMDeviance / dNullDeviance ) ) * 100;


	//
	// display results
	//
	pGDMDev[0] = dGDMDeviance;
	pNullDev[0] = dNullDeviance;
	pExpDev[0] = dDevianceExplained;
	pIntercept[0] = pCoefficients[0];
	for ( int i=1; i<nTotalSplines+1; i++ )
	{
		pCoeffs[i-1] = pCoefficients[i];
	}


	//
	// copy to response data to a vector
	//
	for ( int i=0; i<nRows; i++ )
	{
		pY[i] = pResponse[i];
		pE[i] = CalcDissimilarity( pPredData, pCoefficients, nRows, nTotalSplines+1, i );
		pX[i] = 1.0 - exp(-pE[i]);
	}


	//
	// Clean up
	//
	if (pPredData) delete[] pPredData;
	if ( pCoefficients ) delete[] pCoefficients;
}





//
// Main GDM predict function called from R
//
//
// Note that pData represents a column major matrix
//
void GDM_PredictFromTable(double *pData, 
		                  int *pDoGeo, int *pPreds, int *pRows, 
					      double *pQuantiles, int *pSplines, double *pCoeffs,
					      double *pX) // predicted
{
	//
	// Everything passed from R needs to be dereferenced...
	//
	int nDoGeo = *pDoGeo;
	int nPreds = *pPreds;
	int nRows  = *pRows;
		

	//
	// Get the predictor matrix
	//
	double *pPredData = ConstructMatrix(nDoGeo, pData, pQuantiles, nPreds, pSplines, nRows);
	if ( NULL == pPredData )
	{
		Message("pPredData is NULL", "ERROR");
		return;
	}

	int nTotalSplines = GetTotalSplineCount(pSplines, nPreds);
	//Message(nTotalSplines, "nTotalSplines");
	// pQuantiles will have a length of nPreds * nTotalSplines 
	//ShowQuantiles(pQuantiles, nPreds, pSplines);
	//DebugPredMatrix("PredMatrix.csv", pPredData, nRows, nPreds, pSplines, nTotalSplines+1);


	//
	// copy predicted data to pX
	//
	for ( int i=0; i<nRows; i++ )
	{
		pX[i] = 1.0 - exp(-CalcDissimilarity( pPredData, pCoeffs, nRows, nTotalSplines+1, i ));
	}


	//
	// Clean up
	//
	if (pPredData) delete[] pPredData;
}




//
// Populate pPredData as a transformed GDM predictor to plot in R-Stats
//
void GetPredictorPlotData( double *pPredData, int *pLength,
						   double *pCoefficients,
						   double *pQuantiles,
						   int *pSplines )
{
	int nLen  = *pLength;
	int nSplines = *pSplines;

	double dInc = fabs(pQuantiles[nSplines-1]-pQuantiles[0]) / nLen;
	double dVal = pQuantiles[0];
	for ( int i=0; i<nLen; i++ )
	{
		pPredData[i] = 0;
		for ( int j=0; j<nSplines; j++ )
		{
			if ( j==0 )
			{
				pPredData[i] += pCoefficients[j] * DoSplineCalc( dVal, pQuantiles[j], pQuantiles[j], pQuantiles[j+1] );
			}

			else if ( j == nSplines-1 )
			{
				pPredData[i] += pCoefficients[j] * DoSplineCalc( dVal, pQuantiles[j-1], pQuantiles[j], pQuantiles[j] );
			}

			else
			{
				pPredData[i] += pCoefficients[j] * DoSplineCalc( dVal, pQuantiles[j-1], pQuantiles[j], pQuantiles[j+1] );
			}
		}
		dVal += dInc;
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Local Support Functions ///////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Predict dissimilarity from spline data and coefficients
//
double CalcDissimilarity( double *pData, double *pCoeffs, int nRows, int nCols, int nIndex )
{
	double dVal = 0.0;

	double *pTmp = &pData[nIndex];
	for ( int i=0; i<nCols; i++ )
	{
		dVal += *pTmp * pCoeffs[i];
		pTmp += nRows;
	}
	return(dVal);
}


//
// calculate total number of splines
//
int GetTotalSplineCount(int *pSplines, int nPreds)
{
	int nCount = 0;
	for ( int i=0; i<nPreds; i++ )
	{
		nCount += pSplines[i];
	}
	return(nCount);
}


//
// Construct a regression matrix for predictors
//
double *ConstructMatrix(int nDoGeo, double *pData, double *pQuants, int nPreds, int *pSplines, int nRows)
{
	int nTotalSplines = GetTotalSplineCount(pSplines, nPreds);
	//Message(nTotalSplines, "nTotalSplines");

	//
	// Construct the predictor matrix with an extra column for the intercept
	//
	int nSize = nRows * (nTotalSplines + 1);
	double *pPredData = new double [nSize];
	if ( NULL == pPredData )
	{
		Message("Cannot allocate Predictor Data", "ERROR in ConstructMatrix");
		return(NULL);
	}
	for ( int i=0; i<nSize; i++ ) pPredData[i] = 0.0;


	//
	// Initialise the intercept column
	//
	double *pLoc = &pPredData[0];
	for ( int i=0; i<nRows; i++ ) pLoc[i] = 1.0;
	pLoc += nRows;   // start at the first column after the intercept


	//
	// Initialise for the geographic distance if we are using it
	//
	if (nDoGeo)
	{
		//
		// Initialise the matrix for the geographic followed by the environmental predictors
		//		
		double dMin,dMid,dMax,dVal1,dVal2;
		int CurrentSpline = 0;
		for ( int pred=0; pred<nPreds; pred++ )
		{
			int d1Offset = ((LEADING_COLS+pred-1)*nRows);
			int d2Offset = ((LEADING_COLS+pred+nPreds-2)*nRows);

			int nThis = 0;
			for ( int i=0; i<nRows; i++ )
			{
				if ( pred == 0 )
				{
					double distX = fabs( pData[(COL_SITE1_X0 * nRows)+i] - pData[(COL_SITE2_X1 * nRows)+i] );
					double distY = fabs( pData[(COL_SITE1_Y0 * nRows)+i] - pData[(COL_SITE2_Y1 * nRows)+i] );

					dVal1 = 0.0;
					dVal2 = sqrt( ( distX * distX ) + ( distY * distY ) ); 							
				}
				else
				{
					dVal1 = pData[d1Offset+i];
					dVal2 = pData[d2Offset+i];
				}

				for ( int spl = 0; spl<pSplines[pred]; spl++ )
				{
					if ( spl == 0 )                              // first spline
					{
						dMin = pQuants[CurrentSpline+spl+0];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1]; 				
					}
		
					else if ( spl == (pSplines[pred]-1) )        // last spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+0];
					}

					else	                                     // a middle spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1];
					}

					// calculate the spline values for each site
					double d1a = DoSplineCalc( dVal1, dMin, dMid, dMax );
					double d2a = DoSplineCalc( dVal2, dMin, dMid, dMax );

					// set the distance between sites for this spline
					pLoc[ (spl * nRows) + nThis ] = fabs( d2a - d1a );	

				} // for ( int spl = 0; spl<pSplines[pred]; spl++ )

				// increment for the next row
				++nThis;

			} // for ( int i=0; i<nRows; i++ )

			// increment pointer to the start of the next predictor's spline columns
			pLoc += pSplines[pred] * nRows;

			// increment the current quantile index for the next pred
			CurrentSpline += pSplines[pred];

		} // for ( int pred=0; pred<nPreds; pred++ )
	}


	else
	{
		//
		// Initialise the matrix for the environmental predictors only
		//
		double dMin,dMid,dMax;
		int CurrentSpline = 0;

		for ( int pred=0; pred<nPreds; pred++ )
		{
			int d1Offset = ((LEADING_COLS+pred)*nRows);
			int d2Offset = ((LEADING_COLS+nPreds+pred)*nRows);

			int nThis = 0;
			for ( int i=0; i<nRows; i++ )
			{
				double dVal1 = pData[d1Offset+i];
				double dVal2 = pData[d2Offset+i];

				for ( int spl = 0; spl<pSplines[pred]; spl++ )
				{
					if ( spl == 0 )                              // first spline
					{
						dMin = pQuants[CurrentSpline+spl+0];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1]; 
					
					}
		
					else if ( spl == (pSplines[pred]-1) )        // last spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+0];
					}

					else	                                     // a middle spline
					{
						dMin = pQuants[CurrentSpline+spl-1];
						dMid = pQuants[CurrentSpline+spl+0];
						dMax = pQuants[CurrentSpline+spl+1];
					}


					// calculate the spline values for each site
					double d1a = DoSplineCalc( dVal1, dMin, dMid, dMax );
					double d2a = DoSplineCalc( dVal2, dMin, dMid, dMax );

					// set the distance between sites for this spline
					pLoc[ (spl * nRows) + nThis ] = fabs( d2a - d1a );	

				} // for ( int spl = 0; spl<pSplines[pred]; spl++ )

				// increment for the next row
				++nThis;

			} // for ( int i=0; i<nRows; i++ )

			// increment pointer to the start of the next predictor's spline columns
			pLoc += pSplines[pred] * nRows;

			// increment the current quantile index for the next pred
			CurrentSpline += pSplines[pred];

		} // for ( int pred=0; pred<nPreds; pred++ )

	} // environmental predictors only

	return(pPredData);
}




//
// Calculate the I-Spline value for dVal given quantiles q1, q2, q3
//
double DoSplineCalc( double dVal, double q1, double q2, double q3 )
{
	if ( dVal <= q1 ) return(0.0);

	else if ( dVal >= q3 ) return( 1.0 );

	else if ( ( q1 < dVal ) && ( dVal < q2 ) )
		return( ( ( ( dVal - q1 ) * ( dVal - q1 ) ) / ( ( q2 - q1 ) * ( q3 - q1 ) ) ) );

	else
		return( ( 1.0 - ( ( ( q3 - dVal ) * ( q3 - dVal ) ) / ( ( q3 - q2 ) * ( q3 - q1) ) ) ) );
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Local Debugging Functions /////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Display the contents of the quantile vector
//
void ShowQuantiles(double *pQuants, int nPreds, int *pSplines)
{
	char buff[256];
	double *pTmp = &pQuants[0];
	for ( int i=0; i<nPreds; i++ )
	{
		sprintf( buff, "Quant %d: ", i+1 );

		for ( int j=0; j<pSplines[i]; j++ )
		{
			sprintf( buff, "%s %lf ", buff, *pTmp);			
			++pTmp;
		}
		Message(buff);
	}
}




//
// Write the Predictor Matrix to a comma delimited file
//
void DebugPredMatrix(char *pPath, double *pPredData, int nRows, int nPreds, int *pSplines, int nCols)
{
	FILE *fp = fopen(pPath, "w+t");

	// write the header
	fprintf( fp, "Intercept,");			
	for ( int p=0; p<nPreds; p++ )
	{
		for ( int s=0; s<pSplines[p]; s++ )
		{
			fprintf( fp, "Pred%dSpline%d", p+1,s+1);			
			if ( s < pSplines[p]-1 )
				fprintf( fp, "," );
		}
		if ( p < nPreds-1 )
			fprintf( fp, "," );
		else
			fprintf( fp, "\n" );

	}
	
	for ( int i=0; i<nRows; i++ )
	{
		for ( int j=0; j<nCols; j++ )
		{
			fprintf( fp, "%lf", pPredData[(j*nRows)+i]);

			if ( j < nCols-1 ) 
				fprintf(fp, "," );
			else
				fprintf(fp, "\n" );
		}
	}

	if (fp) fclose(fp);
}


#endif

