//
// TableFuncs.cpp
//
#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "TableFuncs.h"
#include "Message.h"

#pragma warning( disable : 4996 )

#define MAX_HEADER 100000


//
// Return the number of predictors in the composite file.
// These columns start at index [5] and have 2 columns for each predictor
//
int GetNumPredictorsFromHeader(char *lpPath)
{
	int nPreds = 0;

	FILE *fp = fopen(lpPath, "r+t");
	if ( NULL == fp )
	{
		Message("GetNumPredictorsFromHeader()-Cannot open lpPath for Read", "ERROR" );
		return(nPreds);
	}
	
	char *Header = new char [MAX_HEADER];
	if ( NULL == Header )
	{
		fclose(fp);
		Message("GetNumPredictorsFromHeader()-Cannot allocate Header", "ERROR" );
		return(nPreds);
	}

	if ( NULL == fgets(Header,MAX_HEADER,fp) )
	{
		if (Header) delete[] Header;
		fclose(fp);
		Message("GetNumPredictorsFromHeader()-Cannot read Header", "ERROR" );
		return(nPreds);
	}

	char *p = strtok(Header,",\n");  // get the response column
	p = strtok(NULL,",\n");          // get the X0 column
	p = strtok(NULL,",\n");          // get the Y0 column
	p = strtok(NULL,",\n");          // get the X1 column
	p = strtok(NULL,",\n");          // get the Y1 column
	while(1)
	{
		p = strtok(NULL,",\n");
		if ( NULL == p )
			break;

		++nPreds;
	}

	// sanity check that there has been an EVEN number of predictor columns counted
	if ( nPreds % 2 )
	{
		if (Header) delete[] Header;
		fclose(fp);
		Message("GetNumPredictorsFromHeader()- ODD number of Predictor Columns", "ERROR" );
		return(nPreds);		
	}

	// divide by 2 to get the correct number of predictors
	nPreds /= 2;
	
	if (Header) delete[] Header;
	fclose(fp);
	return(nPreds);
}


//
// Return the number of data rows from a GDM Composite Table excluding the header
//
int GetNumRowsFromTableFile(char *lpPath)
{
	int nRows = 0;

	FILE *fp = fopen(lpPath, "r+t");
	if ( NULL == fp )
	{
		Message("GetNumRowsFromTableFile()-Cannot open lpPath for Read", "ERROR" );
		return(nRows);
	}
	
	char *Header = new char [MAX_HEADER];
	if ( NULL == Header )
	{
		fclose(fp);
		Message("GetNumRowsFromTableFile()-Cannot allocate Header", "ERROR" );
		return(nRows);
	}

	if ( NULL == fgets(Header,MAX_HEADER,fp) )
	{
		if (Header) delete[] Header;
		fclose(fp);
		Message("GetNumRowsFromTableFile()-Cannot read Header", "ERROR" );
		return(nRows);
	}

	while(1)
	{
		if ( NULL == fgets(Header,MAX_HEADER,fp) ) 
			break;

		++nRows;
	}

	if (Header) delete[] Header;
	fclose(fp);
	return(nRows);
}


//
// Return the ZERO-BASED Column as a vector of doubles
//
double *GetColumnAt(char *lpPath, int nRows, int Index)
{
	double *pData = new double [nRows];
	if ( NULL == pData )
	{
		Message("GetColumnAt() - Cannot Allocate pData", "ERROR");
		return(NULL);
	}

	// now get the data from the first column
	FILE *fp = fopen(lpPath, "r+t");
	if ( NULL == fp )
	{
		Message("GetColumnAt()-Cannot open lpPath for Read", "ERROR" );
		return(NULL);
	}
	
	char *RowData = new char [MAX_HEADER];
	if ( NULL == RowData )
	{
		fclose(fp);
		Message("GetColumnAt()-Cannot allocate RowData", "ERROR" );
		return(NULL);
	}

	if ( NULL == fgets(RowData,MAX_HEADER,fp) )
	{
		if (RowData) delete[] RowData;
		fclose(fp);
		Message("GetColumnAt()-Cannot read Header", "ERROR" );
		return(NULL);
	}

	
	for ( int i=0; i<nRows; i++ )
	{
		fgets(RowData,MAX_HEADER,fp);     // get the data row
		char *p = strtok(RowData,",\n");  // get the response column

		int j=0; 
		while( j < Index ) 
		{ 
			p = strtok(NULL,",\n");
			++j;
		}

		pData[i] = atof(p);               // convert to double
	}

	if (RowData) delete[] RowData;
	fclose(fp);
	return(pData);
}


//
// Return the ZERO-BASED Column as a vector of doubles
//
double *GetColumnAt(char *lpPath, int Index)
{
	int nRows = GetNumRowsFromTableFile(lpPath);
	if ( nRows < 1 )
	{
		Message("GetColumnAt() - Cannot GetNumRowsFromTableFile", "ERROR");
		return(NULL);
	}

	double *pData = new double [nRows];
	if ( NULL == pData )
	{
		Message("GetColumnAt() - Cannot Allocate pData", "ERROR");
		return(NULL);
	}

	// now get the data from the first column
	FILE *fp = fopen(lpPath, "r+t");
	if ( NULL == fp )
	{
		Message("GetColumnAt()-Cannot open lpPath for Read", "ERROR" );
		return(NULL);
	}
	
	char *RowData = new char [MAX_HEADER];
	if ( NULL == RowData )
	{
		fclose(fp);
		Message("GetColumnAt()-Cannot allocate RowData", "ERROR" );
		return(NULL);
	}

	if ( NULL == fgets(RowData,MAX_HEADER,fp) )
	{
		if (RowData) delete[] RowData;
		fclose(fp);
		Message("GetColumnAt()-Cannot read Header", "ERROR" );
		return(NULL);
	}

	
	for ( int i=0; i<nRows; i++ )
	{
		fgets(RowData,MAX_HEADER,fp);     // get the data row
		char *p = strtok(RowData,",\n");  // get the response column

		int j=0; 
		while( j < Index ) 
		{ 
			p = strtok(NULL,",\n");
			++j;
		}

		pData[i] = atof(p);               // convert to double
	}

	if (RowData) delete[] RowData;
	fclose(fp);
	return(pData);
}


//
// Get the data for the first site column for ZERO-BASED predictor
//
double *GetSiteADataForPredictorAt(char *lpPath, int nRows, int Index0)
{
	return(GetColumnAt(lpPath, nRows, Index0 + 5));
}


//
// Get the data for the first site column for ZERO-BASED predictor
//
double *GetSiteADataForPredictorAt(char *lpPath, int Index0)
{
	return(GetColumnAt(lpPath, Index0 + 5));
}


//
// Get the data for the second site column for ZERO-BASED predictor
//
double *GetSiteBDataForPredictorAt(char *lpPath, int nRows, int Index0)
{
	int nPreds = GetNumPredictorsFromHeader(lpPath);
	return(GetColumnAt(lpPath, nRows, nPreds + Index0 + 5));
}


//
// Get the data for the second site column for ZERO-BASED predictor
//
double *GetSiteBDataForPredictorAt(char *lpPath, int Index0)
{
	int nPreds = GetNumPredictorsFromHeader(lpPath);
	return(GetColumnAt(lpPath, nPreds + Index0 + 5));
}


//
// Dump first n items in the vector of doubles
//
void DumpFirstNItems(int n, double *pData, char *Label)
{
	char tmp[16];
	for ( int i=0; i<n; i++ )
	{
		sprintf(tmp, "Item_%d", i+1);
		Message(tmp, pData[i], Label);
	}
}


//
// Return an I-Spline representation of the Geographic Predictor Data
//
double *GetGeographicISplines(int nRows, int nQuants, char *lpTable)
{
	double *pData = new double [nRows * nQuants];
	if ( NULL == pData )
	{
		Message("GetGeographicISplines() - Cannot allocate pData", "ERROR");
		return(NULL);
	}

	// get the Geographic Distances as a vector
	double *pDistances = GetGeographicDistanceVector(nRows, lpTable);
	if ( NULL == pDistances )
	{
		Message("GetGeographicISplines() - Cannot allocate pDistances", "ERROR");
		if (pData) delete[] pData;
		return(NULL);
	}

	// extract the quantiles from the distance vector
	double *pQuants = ExtractDataQuantiles(nRows, nQuants, pDistances);	

	// finally, populate the I-Spline vector using the values in pDistances
	for ( int i=0; i<nRows; ++i )
	{
		for ( int j=0; j<nQuants; j++ )
		{
			double dMin, dMid, dMax;
			double *pLoc = &pData[j*nRows];

			if ( j == 0 ) 
			{
				// first slot, use min, min, mid for spline calc
				dMin = dMid = pQuants[j];
				dMax = pQuants[j+1];
			}

			else if ( j == nQuants-1 )
			{
				// last slot, use mid, max, max for spline calc
				dMin = pQuants[j-1];
				dMid = dMax = pQuants[j];
			}

			else
			{
				// somewhere in between, use min, mid, max for spline calc
				dMin = pQuants[j-1];
				dMid = pQuants[j];
				dMax = pQuants[j+1];
			}

			// set the cuurent spline value for the geographic distance predictor
			pLoc[i] = doSplineCalc( pDistances[i], dMin, dMid, dMax );
		}
	}

	if (pQuants) delete[] pQuants;
	if (pDistances) delete[] pDistances;
	return(pData);
}


//
// Return an I-Spline representation of the Geographic Predictor Data
//
double *GetGeographicISplines(int nRows, GdmModel *gdmModel, char *lpTable)
{
	int nQuants = gdmModel->GetGeoPredictor()->GetSplines();

	double *pData = new double [nRows * nQuants];
	if ( NULL == pData )
	{
		Message("GetGeographicISplines() - Cannot allocate pData", "ERROR");
		return(NULL);
	}

	// get the Geographic Distances as a vector
	double *pDistances = GetGeographicDistanceVector(nRows, lpTable);
	if ( NULL == pDistances )
	{
		Message("GetGeographicISplines() - Cannot allocate pDistances", "ERROR");
		if (pData) delete[] pData;
		return(NULL);
	}

	// extract the quantiles from the distance vector
	double *pQuants = ExtractDataQuantiles(nRows, nQuants, pDistances);


	// copy quantiles to the GdmMdodel Class
	for ( int i=0; i<gdmModel->GetGeoPredictor()->GetSplines(); i++ )
	{
		gdmModel->GetGeoPredictor()->SetQuantAt(pQuants[i], i);
	}


	// finally, populate the I-Spline vector using the values in pDistances
	for ( int i=0; i<nRows; ++i )
	{
		for ( int j=0; j<nQuants; j++ )
		{
			double dMin, dMid, dMax;
			double *pLoc = &pData[j*nRows];

			if ( j == 0 ) 
			{
				// first slot, use min, min, mid for spline calc
				dMin = dMid = pQuants[j];
				dMax = pQuants[j+1];
			}

			else if ( j == nQuants-1 )
			{
				// last slot, use mid, max, max for spline calc
				dMin = pQuants[j-1];
				dMid = dMax = pQuants[j];
			}

			else
			{
				// somewhere in between, use min, mid, max for spline calc
				dMin = pQuants[j-1];
				dMid = pQuants[j];
				dMax = pQuants[j+1];
			}

			// set the cuurent spline value for the geographic distance predictor
			pLoc[i] = doSplineCalc( pDistances[i], dMin, dMid, dMax );
		}
	}

	if (pQuants) delete[] pQuants;
	if (pDistances) delete[] pDistances;
	return(pData);
}


//
// Return an I-Spline representation of the Geographic Predictor Data
// and update the quantile values for the geographic predictor
//
double *GetPredictorISplines(int nRows, GdmModel *gdmModel, char *lpTablePath, int Index0)
{
	int nQuants = gdmModel->GetPredictorAt(Index0)->GetSplines();

	double *pData = new double [nRows * nQuants];
	if ( NULL == pData )
	{
		Message("GetGeographicISplines() - Cannot allocate pData", "ERROR");
		return(NULL);
	}

	// get the Site Data as a pair of vectors
	double *pData0 = GetSiteADataForPredictorAt(lpTablePath, nRows, Index0);
	double *pData1 = GetSiteBDataForPredictorAt(lpTablePath, nRows, Index0);
	if ( (!pData0) || (!pData1) )
	{
		Message("GetPredictorISplines() - Cannot allocate Site Pair Data Vectors", "ERROR");
		if (pData) delete[] pData;
		if (pData0) delete[] pData0;
		if (pData1) delete[] pData1;
		return(NULL);
	}

	// extract the quantiles from the distance vector
	double *pQuants = ExtractDataQuantiles(nRows, nQuants, pData0, pData1);

	// copy quantile values into the GdmModel Class
	for ( int i=0; i<nQuants; i++ )
	{
		gdmModel->GetPredictorAt(Index0)->SetQuantAt(pQuants[i], i);
	}

	// finally, populate the I-Spline vector using the values in pDistances
	for ( int i=0; i<nRows; ++i )
	{
		for ( int j=0; j<nQuants; j++ )
		{
			double dMin, dMid, dMax;
			double *pLoc = &pData[j*nRows];

			if ( j == 0 ) 
			{
				// first slot, use min, min, mid for spline calc
				dMin = dMid = pQuants[j];
				dMax = pQuants[j+1];
			}

			else if ( j == nQuants-1 )
			{
				// last slot, use mid, max, max for spline calc
				dMin = pQuants[j-1];
				dMid = dMax = pQuants[j];
			}

			else
			{
				// somewhere in between, use min, mid, max for spline calc
				dMin = pQuants[j-1];
				dMid = pQuants[j];
				dMax = pQuants[j+1];
			}

			// set the current spline value for the predictor
			double d0 = doSplineCalc( pData0[i], dMin, dMid, dMax );
			double d1 = doSplineCalc( pData1[i], dMin, dMid, dMax );
			pLoc[i] = fabs(d0-d1);
		}
	}

	if (pQuants) delete[] pQuants;
	if (pData0) delete[] pData0;
	if (pData1) delete[] pData1;
	return(pData);
}


//
// Return an I-Spline representation of the Predictor Data at ZERO-BASE Index0
//
double *GetPredictorISplines(int nRows, int nQuants, char *lpTablePath, int Index0)
{
	double *pData = new double [nRows * nQuants];
	if ( NULL == pData )
	{
		Message("GetPredictorISplines() - Cannot allocate pData", "ERROR");
		return(NULL);
	}

	// get the Site Data as a pair of vectors
	double *pData0 = GetSiteADataForPredictorAt(lpTablePath, nRows, Index0);
	double *pData1 = GetSiteBDataForPredictorAt(lpTablePath, nRows, Index0);
	if ( (!pData0) || (!pData1) )
	{
		Message("GetPredictorISplines() - Cannot allocate Site Pair Data Vectors", "ERROR");
		if (pData) delete[] pData;
		if (pData0) delete[] pData0;
		if (pData1) delete[] pData1;
		return(NULL);
	}

	// extract the quantiles from the distance vector
	double *pQuants = ExtractDataQuantiles(nRows, nQuants, pData0, pData1);
	
	// finally, populate the I-Spline vector using the values in pDistances
	for ( int i=0; i<nRows; ++i )
	{
		for ( int j=0; j<nQuants; j++ )
		{
			double dMin, dMid, dMax;
			double *pLoc = &pData[j*nRows];

			if ( j == 0 ) 
			{
				// first slot, use min, min, mid for spline calc
				dMin = dMid = pQuants[j];
				dMax = pQuants[j+1];
			}

			else if ( j == nQuants-1 )
			{
				// last slot, use mid, max, max for spline calc
				dMin = pQuants[j-1];
				dMid = dMax = pQuants[j];
			}

			else
			{
				// somewhere in between, use min, mid, max for spline calc
				dMin = pQuants[j-1];
				dMid = pQuants[j];
				dMax = pQuants[j+1];
			}

			// set the current spline value for the predictor
			double d0 = doSplineCalc( pData0[i], dMin, dMid, dMax );
			double d1 = doSplineCalc( pData1[i], dMin, dMid, dMax );
			pLoc[i] = fabs(d0-d1);
		}
	}

	if (pQuants) delete[] pQuants;
	if (pData0) delete[] pData0;
	if (pData1) delete[] pData1;
	return(pData);
}


//
// Return a vector of Euclidean Distances between sites
//
double *GetGeographicDistanceVector(int nRows, char *lpTable)
{
	double *pData = new double [nRows];
	if ( NULL == pData )
	{
		Message("GetGeographicDistanceVector() - Cannot allocate pData", "ERROR");
		return(NULL);
	}

	double *pX0 = GetColumnAt(lpTable, nRows, C_X0);
	double *pY0 = GetColumnAt(lpTable, nRows, C_Y0);
	double *pX1 = GetColumnAt(lpTable, nRows, C_X1);
	double *pY1 = GetColumnAt(lpTable, nRows, C_Y1);

	for ( int i=0; i<nRows; i++ )
	{
		pData[i] = CalcGeographicDistance(pX0[i], pY0[i], pX1[i], pY1[i]);
	}

	if (pX0) delete[] pX0;
	if (pY0) delete[] pY0;
	if (pX1) delete[] pX1;
	if (pY1) delete[] pY1;

	return(pData);
}


//
// Calculate the Euclidean Distance of a site pair
//
double CalcGeographicDistance(double X0, double Y0, double X1, double Y1)
{
	return(sqrt(((X0-X1)*(X0-X1))+((Y0-Y1)*(Y0-Y1))));
}


//
// Extract nQuants Quantiles from pData
//
double *ExtractDataQuantiles(int nRows, int nQuants, double *pData)
{
	double *pQuantiles = new double [nQuants];

	// copy the data and sort it
	double *ptmp = new double [nRows];
	if ( NULL == ptmp )
	{
		Message("ExtractDataQuantiles() - Cannot allocate ptmp", "ERROR");
		return(NULL);
	}
	memcpy(ptmp, pData, nRows * sizeof(double));
	qsort(ptmp,nRows,sizeof(double),dcompare);

	// now extract quantiles at appropriate positions in the sorted data
	pQuantiles[0] = ptmp[0];
	pQuantiles[nQuants-1] = ptmp[nRows-1];
	// middle quantiles...
	for ( int i=2; i<nQuants; i++ )
	{
		pQuantiles[i-1] = ptmp[nRows/(nQuants-1)*(i-1)];
	}

	// cleanup and return
	if (ptmp) delete[] ptmp;
	return(pQuantiles);
}


//
// Extract nQuants Quantiles from a composite vector made from pData0 and pData1 (used for geographic distance)
//
double *ExtractDataQuantiles(int nRows, int nQuants, double *pData0, double *pData1)
{
	double *pQuantiles = new double [nQuants];

	// copy the data and sort it
	double *ptmp = new double [nRows*2];
	if ( NULL == ptmp )
	{
		Message("ExtractDataQuantiles() - Cannot allocate ptmp", "ERROR");
		return(NULL);
	}
	memcpy(ptmp, pData0, nRows * sizeof(double));
	memcpy(&ptmp[nRows], pData1, nRows * sizeof(double));
	qsort(ptmp,nRows*2,sizeof(double),dcompare);

	// now extract quantiles at appropriate positions in the sorted data
	pQuantiles[0] = ptmp[0];
	pQuantiles[nQuants-1] = ptmp[(nRows*2)-1];
	// middle quantiles...
	for ( int i=2; i<nQuants; i++ )
	{
		pQuantiles[i-1] = ptmp[(nRows*2)/(nQuants-1)*(i-1)];
	}

	// cleanup and return
	if (ptmp) delete[] ptmp;
	return(pQuantiles);
}


//
// comparison function for quick sorting doubles
//
int dcompare( const void *p0, const void *p1 )
{
	double d0 = *((double *)p0);
	double d1 = *((double *)p1);

	if ( d0 < d1 )
		return( -1 );
	else if ( d0 > d1 )
		return( 1 );
	else
		return( 0 );
}



//
// Calculate the I-Spline value for dVal given quantiles q1, q2, q3
//
double doSplineCalc( double dVal, double q1, double q2, double q3 )
{
	if ( dVal <= q1 ) return( 0.0 );

	else if ( dVal >= q3 ) return( 1.0 );

	else if ( ( q1 < dVal ) && ( dVal < q2 ) )
	{
		return( ( ( ( dVal - q1 ) * ( dVal - q1 ) ) / ( ( q2 - q1 ) * ( q3 - q1 ) ) ) );
	}

	else
	{
		return( ( 1.0 - ( ( ( q3 - dVal ) * ( q3 - dVal ) ) / ( ( q3 - q2 ) * ( q3 - q1 ) ) ) ) );
	}
}


