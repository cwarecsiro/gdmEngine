//
// GdmExtractQuants.cpp
//
#include "stdafx.h"
#include "GdmBinLib.h"
#include "GdmExtractQuants.h"
#include "GdmGuiLib.h"
#include "Message.h"
#include "ParamsW16.h"
#include "clsDoPath.h"
#include "GdmTabLib.h"
#include "GDM_PredictorTypes.h"
#include "GDMBufferSizeDefs.h"
#include "Gdm_GreatCircle_Calcs.h"

//enum predState { psData=0, psGrid=1, psUser=2 };
//if (predState::psData == GetEuclideanQuantType(pParams)) {}
//else if (predState::psGrid == GetEuclideanQuantType(pParams)) {}
//else if (predState::psUser == GetEuclideanQuantType(pParams)) {}

//
// Return the data table after extracting the quantiles and updating the parameter file
//
double *GetTableDataPostExtractAndUpdateQuantiles( char *pParams, 
	                                               GDM_INT *pRows, 
												   GDM_INT *pCols, 
												   bool BatchRun, 
												   FPTR fptr)
{
	//Message("GetTableDataPostExtractAndUpdateQuantiles");

	//
	// Get the data table path
	//
	char *DataTablePath = GetResponseData(pParams);
	gmpath gmp;
	if (!gmp.FileExists(DataTablePath))
	{
		if (!BatchRun)
		{
			Message(DataTablePath, "Cannot Locate in GetTableDataPostExtractAndUpdateQuantiles");
		}
		if (DataTablePath) delete[] DataTablePath;
		return(NULL);
	}

	
	//
	// Get a memory image of the composite data file
	//
	GDM_INT nTableRows = 0;
	GDM_INT nTableCols = 0;
	double *pTableData = GetCompositeTableIntoMemory(DataTablePath, &nTableRows, &nTableCols, BatchRun, fptr);
	if (NULL == pTableData)
	{
		if (!BatchRun)
			Message("Cannot GetCompositeTableIntoMemory() in GetTableDataPostExtractAndUpdateQuantiles", "ERROR");
		if (DataTablePath) delete[] DataTablePath;
		return(NULL);
	}
	// Free the DataTablePath now that we don't need it anymore
	if (DataTablePath) delete[] DataTablePath;


	//
	// Extract the geographic distance data and write a binary file for Response,Weights,X0,Y0,X1,Y1 on the way thru...
	//
	char *lpColFile = GetWorkspacePath(pParams);
	strcat(lpColFile, "\\GDM_SixCols.bin");
	ofstream *ColOut = new ofstream(lpColFile, ios::binary | ios::trunc);
	if (!ColOut->is_open())
	{
		Message("Cannot create GDM_SixCols.bin in GetTableDataPostExtractAndUpdateQuantiles()", "ERROR");
		if (pTableData) delete[] pTableData;
		if (DataTablePath) delete[] DataTablePath;
		return(NULL);
	}
	
	int nCurrent = 0;
	if (fptr) fptr("Extracting geographic distance...", nCurrent );
	double *pGeoDist = new double [nTableRows];
	for ( int i=0; i<nTableRows; i++ )
	{
		if (fptr)
		{
			if ( (long long)(i) * 100 / nTableRows > nCurrent )
			{
				nCurrent = (long long)(i) * 100 / (int)nTableRows;
				fptr("Extracting geographic distance...", nCurrent );
			}
		}

		// OLD
		//double distX = pTableData[(COL_SITE1_X0*nTableRows)+i] - pTableData[(COL_SITE2_X1*nTableRows)+i];
		//double distY = pTableData[(COL_SITE1_Y0*nTableRows)+i] - pTableData[(COL_SITE2_Y1*nTableRows)+i];
		//pGeoDist[i] = sqrt((distX * distX) + (distY * distY));
		//
		//
		// use the great circle distance calc instead
		//
		pGeoDist[i] = CalcGreatCircleDist((COL_SITE1_X0*(int)nTableRows)+i,
			                              (COL_SITE2_X1*(int)nTableRows)+i,
										  (COL_SITE1_Y0*(int)nTableRows)+i,
										  (COL_SITE2_Y1*(int)nTableRows)+i);

		//
		// create a binary file of the Response,Weights,X0,Y0,X1,Y1 as 8 byte doubles
		//
		ColOut->write((char *)&pTableData[(COL_RESPONSE*nTableRows)+i], sizeof(double));
		ColOut->write((char *)&pTableData[(COL_WEIGHTS*nTableRows)+i], sizeof(double));

		ColOut->write((char *)&pTableData[(COL_SITE1_X0*nTableRows)+i], sizeof(double));
		ColOut->write((char *)&pTableData[(COL_SITE1_Y0*nTableRows)+i], sizeof(double));
		ColOut->write((char *)&pTableData[(COL_SITE2_X1*nTableRows)+i], sizeof(double));
		ColOut->write((char *)&pTableData[(COL_SITE2_Y1*nTableRows)+i], sizeof(double));
	}
	ColOut->close();

	// sort the vector
	SortVector(pGeoDist, (int)nTableRows);

	// extract the geographic distance quantiles and update the parameter file
	int nQuants = GetEuclideanSplines(pParams);
	char lpKey[64];
	for ( int i=0; i<nQuants; i++ )
	{
		if ((i>0) && (i<nQuants-1))
		{
			int nPercentile = (int)(i * 100.0 / (double)(nQuants-1));
			double dQuant = GetQuantAsPercentile( pGeoDist, (int)nTableRows, nPercentile );

			sprintf( lpKey, "EuclSplVal%d", i+1 );
			SetProfileDouble( "PREDICTORS", lpKey, dQuant, pParams );
		}
	}
	nCurrent = 0;
	if (fptr) fptr("Extracting geographic distance...", nCurrent );

	// clean up
	if (pGeoDist) delete[] pGeoDist; 


	//
	// Extract the quantiles from the 'In-Use" predictors
	//
	double *pPredDist = new double [nTableRows*2];
	int nPreds = GetNumPredictors(pParams);
	nCurrent = 0;
	if (fptr) fptr("Extracting predictor data(1)...", nCurrent );
	for ( int p=0; p<nPreds; p++ )
	{
		if (fptr)
		{
			if ( (p+1) * 100 / nPreds > nCurrent )
			{
				nCurrent = (p+1) * 100 / nPreds;
				fptr("Extracting predictor data(1)...", nCurrent );
			}
		}
		
		if (GetPredictorInUseAt(pParams, p+1))
		{	
			int nThis = 0;
			for ( int i=0; i<nTableRows; i++ )
			{
				if (PRED_TYPE_GRID == GetPredictorTypeAt(pParams, p+1)) // grid predictor
				{
					pPredDist[nThis++] = pTableData[((LEADING_COLS+p)*nTableRows)+i];
					pPredDist[nThis++] = pTableData[((LEADING_COLS+p+nPreds)*nTableRows)+i];
				}

				else if (PRED_TYPE_COVARIANT == GetPredictorTypeAt(pParams, p+1)) // covariant predictor
				{
					pPredDist[nThis++] = pTableData[((LEADING_COLS+p+nPreds)*nTableRows)+i];
				}
			}

			// sort the vector
			if (PRED_TYPE_GRID == GetPredictorTypeAt(pParams, p+1)) // grid predictor
				SortVector(pPredDist, (int)nTableRows*2);

			else if (PRED_TYPE_COVARIANT == GetPredictorTypeAt(pParams, p+1)) // covariant predictor
				SortVector(pPredDist, (int)nTableRows);


			//
			// if the predictor type is CATEGORICAL, we need to extract the 
			// data from the lower triangular part of the lookup table
			//
			int nCatTableLength = 0;
			if (PRED_TYPE_CATEGORICAL == GetPredictorTypeAt(pParams, p+1)) // categorical predictor
			{
				nCatTableLength = GetLowerTriFromLookup(pParams, p+1, pPredDist);
				//Message(nCatTableLength, "nCatTableLength");
				SortVector(pPredDist, nCatTableLength);
			}


			// extract the predictor quantiles and update the parameter file
			nQuants = GetPredictorSplinesAt(pParams, p+1);
			for ( int i=0; i<nQuants; i++ )
			{
				int nPercentile = (int)(i * 100.0 / (double)(nQuants-1));
				double dQuant = 0.0;				

				if (PRED_TYPE_CATEGORICAL == GetPredictorTypeAt(pParams, p+1)) // categorical predictor
				{
					dQuant = GetQuantAsPercentile( pPredDist, nCatTableLength, nPercentile );
					SetPredictorSplineAt(pParams, p+1, i+1, dQuant);
				}

				if ((i>0) && (i<nQuants-1))
				{
					if (PRED_TYPE_GRID == GetPredictorTypeAt(pParams, p+1)) // grid predictor
					{
						dQuant = GetQuantAsPercentile( pPredDist, (int)nTableRows*2, nPercentile );
						SetPredictorSplineAt(pParams, p+1, i+1, dQuant);
					}

					else if (PRED_TYPE_COVARIANT == GetPredictorTypeAt(pParams, p+1)) // covariant predictor
					{
						dQuant = GetQuantAsPercentile( pPredDist, (int)nTableRows, nPercentile );
						SetPredictorSplineAt(pParams, p+1, i+1, dQuant);
					}
				} // if ((i>0) && (i<nQuants-1))
			} // for ( int i=0; i<nQuants; i++ )
		} // if (GetPredictorInUseAt(pParams, p+1))
	} // for ( int p=0; p<nPreds; p++ )

	nCurrent = 0;
	fptr("Extracting predictor data(1)...", nCurrent );
	if ( pPredDist) delete[] pPredDist;
	*pRows = nTableRows;
	*pCols = nTableCols;
	return(pTableData);
}



//
// Extract lower triangular data from lookup table into pPredDist and return number of items
//
int GetLowerTriFromLookup(char *pParams, int Index, double *pPredDist)
{
	// extract the lookup table filename
	char catstring[FILEPATH_BUFFSIZE];	
	strcpy(catstring, GetPredictorPathAt(pParams, Index));
	char *p = strtok(catstring, "+"); // get the first half of the string
	p = strtok(NULL, "\n"); // get the last half of the string


	// determine size
	char *buff = new char [TABLE_ROW_BUFFSIZE];
	FILE *fp = fopen(p, "r+t");
	fgets(buff, TABLE_ROW_BUFFSIZE, fp); // skip header
	int nRows = 0;
	while(1)
	{
		if (NULL == fgets(buff, TABLE_ROW_BUFFSIZE, fp))
			break;
		++nRows;
	}

	// copy the lower triangular items
	int NumItems = 0;
	rewind(fp);
	fgets(buff, TABLE_ROW_BUFFSIZE, fp); // skip header
	for ( int i=0; i<nRows; i++ )
	{
		fgets(buff, TABLE_ROW_BUFFSIZE, fp);

		char *p = strtok(buff, ",");  // skip first label column
		for ( int j=0; j<i; j++ )
		{
			p = strtok(NULL, ",");
			pPredDist[NumItems] = atof(p);
			++NumItems;
		}
	}
	fclose(fp);
	if (buff) delete[] buff;
	return(NumItems);
}



//
// Extract Quantiles from Table Data and update parameter file
//
bool ExtractAndUpdateQuantiles( char *pParams, bool BatchRun, FPTR fptr)
{
	//Message("ExtractAndUpdateQuantiles");

	//
	// Get the data table path
	//
	char *DataTablePath = GetResponseData(pParams);
	//gmpath gmp;
	/*if (!gmp.FileExists(DataTablePath))
	{
		if (!BatchRun)
		{
			Message(DataTablePath, "Cannot Locate in ExtractAndUpdateQuantiles");
		}
		if (DataTablePath) delete[] DataTablePath;
		return(false);
	}*/


	//
	// Get a memory image of the composite data file
	//
	GDM_INT nTableRows = 0;
	GDM_INT nTableCols = 0;
	double *pTableData = GetCompositeTableIntoMemory(DataTablePath, &nTableRows, &nTableCols, BatchRun, fptr);
	if (NULL == pTableData)
	{
		if (!BatchRun)
			Message("Cannot GetCompositeTableIntoMemory() in ExtractAndUpdateQuantiles", "ERROR");
		if (DataTablePath) delete[] DataTablePath;
		return(false);
	}


	// Free the DataTablePath now that we don't need it anymore
	if (DataTablePath) delete[] DataTablePath;


	//
	// Extract the geographic distance data
	//
	int nCurrent = 0;
	if (fptr) fptr("Extracting geographic distance...", nCurrent );
	double *pGeoDist = new double [nTableRows];
	for ( int i=0; i<nTableRows; i++ )
	{
		if (fptr)
		{
			if ( i * 100 / nTableRows > nCurrent )
			{
				nCurrent = i * 100 / (int)nTableRows;
				fptr("Extracting geographic distance...", nCurrent );
			}
		}


		//
		// OLD
		//double distX = pTableData[(COL_SITE1_X0*nTableRows)+i] - pTableData[(COL_SITE2_X1*nTableRows)+i];
		//double distY = pTableData[(COL_SITE1_Y0*nTableRows)+i] - pTableData[(COL_SITE2_Y1*nTableRows)+i];
		//pGeoDist[i] = sqrt((distX * distX) + (distY * distY));

		pGeoDist[i] = CalcGreatCircleDist((COL_SITE1_X0*(int)nTableRows)+i, 
			                              (COL_SITE2_X1*(int)nTableRows)+i,
										  (COL_SITE1_Y0*(int)nTableRows)+i,
										  (COL_SITE2_Y1*(int)nTableRows)+i);
	}

	// sort the vector
	SortVector(pGeoDist, (int)nTableRows);

	// extract the geographic distance quantiles and update the parameter file
	int nQuants = GetEuclideanSplines(pParams);
	char lpKey[64];
	for ( int i=0; i<nQuants; i++ )
	{
		if ((i>0) && (i<nQuants-1))
		{
			int nPercentile = (int)(i * 100.0 / (double)(nQuants-1));
			double dQuant = GetQuantAsPercentile( pGeoDist, (int)nTableRows, nPercentile );

			sprintf( lpKey, "EuclSplVal%d", i+1 );
			SetProfileDouble( "PREDICTORS", lpKey, dQuant, pParams );
		}
	}
	nCurrent = 0;
	if (fptr) fptr("Extracting geographic distance...", nCurrent );

	// clean up
	if (pGeoDist) delete[] pGeoDist; 


	//
	// Extract the quantiles from the 'In-Use" predictors
	//
	double *pPredDist = new double [nTableRows*2];
	int nPreds = GetNumPredictors(pParams);
	nCurrent = 0;
	if (fptr) fptr("Extracting predictor data(2)...", nCurrent );
	for ( int p=0; p<nPreds; p++ )
	{
		if (fptr)
		{
			if ( (p+1) * 100 / nPreds > nCurrent )
			{
				nCurrent = (p+1) * 100 / nPreds;
				fptr("Extracting predictor data(2)...", nCurrent );
			}
		}

		if (GetPredictorInUseAt(pParams, p+1))
		{
			int nThis = 0;
			for ( int i=0; i<nTableRows; i++ )
			{
				pPredDist[nThis++] = pTableData[((LEADING_COLS+p)*nTableRows)+i];
				pPredDist[nThis++] = pTableData[((LEADING_COLS+p+nPreds)*nTableRows)+i];
			}

			// sort the vector
			SortVector(pPredDist, (int)nTableRows*2);

			// extract the predictor quantiles and update the parameter file
			nQuants = GetPredictorSplinesAt(pParams, p+1);
			for ( int i=0; i<nQuants; i++ )
			{
				if ((i>0) && (i<nQuants-1))
				{
					int nPercentile = (int)(i * 100.0 / (double)(nQuants-1));
					double dQuant = GetQuantAsPercentile( pPredDist, (int)nTableRows*2, nPercentile );
					SetPredictorSplineAt(pParams, p+1, i+1, dQuant);
				} // if ((i>0) && (i<nQuants-1))
			} // for ( int i=0; i<nQuants; i++ )
		} // if (GetPredictorInUseAt(pParams, p+1))
	} // for ( int p=0; p<nPreds; p++ )

	nCurrent = 0;
	fptr("Extracting predictor data(2)...", nCurrent );

	if ( pPredDist) delete[] pPredDist;
	if (pTableData) delete[] pTableData;
	return(true);
}



//
// Extract ALL the Predictor Quantiles from Table Data and update parameter file
//
bool ExtractAndUpdateAllQuantiles( char *pParams, bool BatchRun, FPTR fptr)
{
	//Message("ExtractAndUpdateAllQuantiles");

	//
	// Get the data table path
	//
	char *DataTablePath = GetResponseData(pParams);
	gmpath gmp;
	if (!gmp.FileExists(DataTablePath))
	{
		if (!BatchRun)
		{
			Message(DataTablePath, "Cannot Locate in ExtractAndUpdateAllQuantiles");
		}
		if (DataTablePath) delete[] DataTablePath;
		return(false);
	}

	
	//
	// Get a memory image of the composite data file
	//
	GDM_INT nTableRows = 0;
	GDM_INT nTableCols = 0;
	double *pTableData = GetCompositeTableIntoMemory(DataTablePath, &nTableRows, &nTableCols, BatchRun, fptr);
	if (NULL == pTableData)
	{
		if (!BatchRun)
			Message("Cannot GetCompositeTableIntoMemory() in ExtractAndUpdateAllQuantiles", "ERROR");
		if (DataTablePath) delete[] DataTablePath;
		return(false);
	}
	// Free the DataTablePath now that we don't need it anymore
	if (DataTablePath) delete[] DataTablePath;


	//
	// Extract the geographic distance data
	//
	int nCurrent = 0;
	if (fptr) fptr("Extracting geographic distance...", nCurrent );
	double *pGeoDist = new double [nTableRows];
	for ( int i=0; i<nTableRows; i++ )
	{
		if (fptr)
		{
			if ( i * 100 / nTableRows > nCurrent )
			{
				nCurrent = i * 100 / (int)nTableRows;
				fptr("Extracting geographic distance...", nCurrent );
			}
		}

		//
		// OLD
		//double distX = pTableData[(COL_SITE1_X0*nTableRows)+i] - pTableData[(COL_SITE2_X1*nTableRows)+i];
		//double distY = pTableData[(COL_SITE1_Y0*nTableRows)+i] - pTableData[(COL_SITE2_Y1*nTableRows)+i];
		//pGeoDist[i] = sqrt((distX * distX) + (distY * distY));

		pGeoDist[i] = CalcGreatCircleDist((COL_SITE1_X0*(int)nTableRows)+i,
			                              (COL_SITE2_X1*(int)nTableRows)+i,
										  (COL_SITE1_Y0*(int)nTableRows)+i,
										  (COL_SITE2_Y1*(int)nTableRows)+i);

	}

	// sort the vector
	SortVector(pGeoDist, (int)nTableRows);

	// extract the geographic distance quantiles and update the parameter file
	int nQuants = GetEuclideanSplines(pParams);
	char lpKey[64];
	for ( int i=0; i<nQuants; i++ )
	{
		if (i==0)
		{
			sprintf( lpKey, "EuclSplVal%d", i+1 );
			SetProfileDouble( "PREDICTORS", lpKey, pGeoDist[0], pParams );
		}
		else if (i==nQuants-1)
		{
			sprintf( lpKey, "EuclSplVal%d", i+1 );
			SetProfileDouble( "PREDICTORS", lpKey, pGeoDist[nTableRows-1], pParams );
		}
		else //if ((i>0) && (i<nQuants-1))
		{
			int nPercentile = (int)(i * 100.0 / (double)(nQuants-1));
			double dQuant = GetQuantAsPercentile( pGeoDist, (int)nTableRows, nPercentile );

			sprintf( lpKey, "EuclSplVal%d", i+1 );
			SetProfileDouble( "PREDICTORS", lpKey, dQuant, pParams );
		}
	}
	nCurrent = 0;
	if (fptr) fptr("Extracting geographic distance...", nCurrent );

	// clean up
	if (pGeoDist) delete[] pGeoDist; 


	//
	// Extract the quantiles from the 'In-Use" predictors
	//
	double *pPredDist = new double [nTableRows*2];
	int nPreds = GetNumPredictors(pParams);
	nCurrent = 0;
	if (fptr) fptr("Extracting predictor data(3)...", nCurrent );
	for ( int p=0; p<nPreds; p++ )
	{
		if (fptr)
		{
			if ( (p+1) * 100 / nPreds > nCurrent )
			{
				nCurrent = (p+1) * 100 / nPreds;
				fptr("Extracting predictor data(3)...", nCurrent );
			}
		}

		int nThis = 0;
		for ( int i=0; i<nTableRows; i++ )
		{
			pPredDist[nThis++] = pTableData[((LEADING_COLS+p)*nTableRows)+i];
			pPredDist[nThis++] = pTableData[((LEADING_COLS+p+nPreds)*nTableRows)+i];
		}

		// sort the vector
		SortVector(pPredDist, (int)nTableRows*2);

		// extract the predictor quantiles and update the parameter file
		nQuants = GetPredictorSplinesAt(pParams, p+1);
		for ( int i=0; i<nQuants; i++ )
		{
			if (i==0)
			{
				SetPredictorSplineAt(pParams, p+1, i+1, pPredDist[0]);
			}

			else if (i==nQuants-1)
			{
				SetPredictorSplineAt(pParams, p+1, i+1, pPredDist[(nTableRows*2)-1]);
			}

			else
			{
				int nPercentile = (int)(i * 100.0 / (double)(nQuants-1));
				double dQuant = GetQuantAsPercentile( pPredDist, (int)nTableRows*2, nPercentile );
				SetPredictorSplineAt(pParams, p+1, i+1, dQuant);
			} // if ((i>0) && (i<nQuants-1))
		} // for ( int i=0; i<nQuants; i++ )
	} // for ( int p=0; p<nPreds; p++ )

	nCurrent = 0;
	fptr("Extracting predictor data(3)...", nCurrent );

	if ( pPredDist) delete[] pPredDist;
	if (pTableData) delete[] pTableData;
	return(true);
}



//
// sort a vector of doubles in ascending order
//
void SortVector(double *pDistVector, int NumItems)
{
	qsort( (void *)pDistVector, (size_t)NumItems, sizeof(double), CompareDoubles );
}



//
// comparison routine for quicksort of doubles
//
int CompareDoubles( const void *arg1, const void *arg2 )
{
	double *p1 = (double *)arg1;
	double *p2 = (double *)arg2;
	return( *p1 > *p2 ? 1 : *p1 < *p2 ? -1 : 0 );
}



//
// return the n% quantile of a sorted vector of doubles
//
double GetQuantAsPercentile( double *pVector, int nLen, int Percentile )
{
	// minimum
	if ( Percentile == 0 )
	{
		return(pVector[0]);
	}

	// maximum
	else if ( Percentile == 100 )
	{
		return(pVector[nLen-1]);
	}	

	// somewhere in between
	else
	{
		double dtmp = Percentile / 100.0;
		if ( fmod( nLen, dtmp ) < 0.5 )
		{
			return( pVector[ (int)floor(nLen * dtmp) ] );
		}
		else
		{
			return( pVector[ (int)ceil(nLen * dtmp) ] );
		}
	}
}

