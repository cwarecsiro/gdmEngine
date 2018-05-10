//
// Gi2CsvLib.cpp
//
#include "stdafx.h"

#include "GdmBinLib.h"
#include "Message.h"
#include "clsDoPath.h"

#include "clsBinaryFileClass.h"
#include "clsEsriBinaryHeader.h"

/////////////////////////////////////////////////////////////////////////////////////////////
// get these in a header file sometime...
bool CheckGrids(char *ParcelPath, char *DataPath);
int GetParcelCount(char *InputPath);
bool CountCells(char *ParcelPath, int *pCounts, int NumParcels);
bool WriteOutputTable(char *InputPath, char *OutputPath, int *pCounts, int NumParcels);
bool WriteOutputTable2(char *InputPath, char *OutputPath, 
	                   int *pCounts, float *pMin, float *pMax, float *pMean, float *pMedian, float *pStdDev, 
					   int NumParcels);
int fCompare( const void *arg1, const void *arg2 );
/////////////////////////////////////////////////////////////////////////////////////////////


void Test( char *ParcelPath, char *DataPath, char *InputPath, char *OutputPath, FPTR fptr )
{
	Message(ParcelPath, "ParcelPath");
	Message(DataPath, "DataPath");
	Message(InputPath, "InPath");
	Message(OutputPath, "OutputPath");
	for ( int i=0; i<100; i++ )
	{
		fptr("Testing...", i);
		for ( int j=0; j<0XFFFFFFF; j++ );
	}
}



bool RunGi2Csv( char *ParcelPath, char *DataPath, char *InputPath, char *OutputPath, bool InitialStep, FPTR fptr )
{
	/*Message(ParcelPath, "ParcelPath in RunGi2Csv");
	Message(DataPath, "DataPath");
	Message(InputPath, "InPath");
	Message(OutputPath, "OutputPath");*/

	// check that the shapefile and data binary grids are coincident
	if (!CheckGrids(ParcelPath, DataPath))
	{
		Message("Input grids are NOT coincident", "ERROR");
		return(false);
	}

	// count number of parcel records in shapefile binary grid
	int NumParcels = GetParcelCount(InputPath);
	if (NumParcels < 1)
	{
		Message("Number of Parcels Less Than 1", "ERROR");
		return(false);
	}

	// allocate the cell counter memory block
	int *pCounts = new int [NumParcels];
	for (int i=0; i<NumParcels; i++ )
	{
		pCounts[i] = 0;
	}

	// count the number of cells in each parcel and fill in the cell count array
	if (!CountCells(ParcelPath, pCounts, NumParcels))
	{
		Message("Cannot count cells", "ERROR");
		if (pCounts) delete[] pCounts;
		return(false);
	}

	// allocate and initialise the support memory blocks
	int *pCurrentIndex = new int [NumParcels];
	float *pPointValues = new float [NumParcels];
	float *pMin = new float [NumParcels];
	float *pMax = new float [NumParcels];
	float *pMean = new float [NumParcels];
	float *pMedian = new float [NumParcels];
	float *pStdDev = new float [NumParcels];
	float **ppData = new float * [NumParcels];
	for ( int i=0; i<NumParcels; i++ )
	{
		if (pCounts[i] > 0) 
		{
			ppData[i] = new float [pCounts[i]];
			for ( int j=0; j<pCounts[i]; j++ )
				ppData[i][j] = 0.0F;
		}
		else
			ppData[i] = NULL;

		pCurrentIndex[i] = 0;
		pPointValues[i] = pMin[i] = pMax[i] = pMean[i] = pMedian[i] = pStdDev[i] = 0.0F;
	}


	/////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////// Extract data from the binary grid into the memory blocks ////////////////////
	//
	gmpath gmPath;
	EsriBinaryHeader *ParcelHeader = new EsriBinaryHeader(gmPath.ChangeExtension(ParcelPath, ".hdr"));
	int nRows = ParcelHeader->GetNumRows();
	int nCols = ParcelHeader->GetNumCols();
	float fNoData = ParcelHeader->GetNoDataValue();
	double dMinX = ParcelHeader->GetMinX();
	double dMaxY = ParcelHeader->GetMaxY();
	double dSize = ParcelHeader->GetCellSize();
	if (ParcelHeader) delete ParcelHeader;

	float *pRowData = new float [nCols];
	BinaryFileClass *bfc = new BinaryFileClass(ParcelPath);
	float *pGridData = new float [nCols];
	BinaryFileClass *bfcdata = new BinaryFileClass(DataPath);
	int nCurrent = 0;
	for ( int i=0; i<nRows; i++ )
	{
		if (i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			fptr("Extracting Data...", nCurrent);			
		}

		bfc->ReadFloat(pRowData, nCols);
		bfcdata->ReadFloat(pGridData, nCols);
		for ( int j=0; j<nCols; j++ )
		{
			if ((pRowData[j] != fNoData) && (pGridData[j] != fNoData))
			{
				int nParcel = int(pRowData[j]);

				// sanity check...
				if (nParcel > NumParcels)
				{
					Message("Index exceeds NumParcels", "ERROR");
					bfc->Close();
					if (pRowData) delete[] pRowData;
					if (pGridData) delete[] pGridData;
					return(false);
				}

				if (pCounts[nParcel-1] > 0)
				{
					// convert from one-based to zero-based index
					ppData[nParcel-1][pCurrentIndex[nParcel-1]] = pGridData[j];
					pCurrentIndex[nParcel-1] += 1; // increment data slot index for this parcel
				}
			}
		}
	}
	bfc->Close();	
	bfcdata->SeekToStart();
	if (pRowData) delete[] pRowData;
	if (pGridData) delete[] pGridData;


	/////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////// Extract data from the binary grid for the singular points that have a count /////////
	/////////////// of zero but have non-zero X and Y Coodinates in the seed table //////////////////////
	//
	//FILE *fpParcels = fopen(InputPath, "r+t");
	//char buff[1000];
	//fgets(buff, 1000, fpParcels);
	//for ( int i=0; i<NumParcels; i++ )
	//{
	//	fgets(buff, 1000, fpParcels);

	//	if (pCounts[i] == 0)
	//	{
	//		char *p = strtok(buff, ",\n");   // get objectID

	//		p = strtok(NULL, ",\n");         // get area

	//		p = strtok(NULL, ",\n");         // get X_Coord
	//		double dX = atof(p);

	//		p = strtok(NULL, ",\n");         // get Y_Coord
	//		double dY = atof(p);

	//		if ((dX != 0.0) && (dY != 0.0))  // we have a coordinate to extract data from
	//		{
	//			long lOffset = long((fabs(dMaxY - dY) / dSize) * nCols) + long(fabs(dX - dMinX) / dSize);
	//			bfcdata->SeekTo(lOffset * sizeof(float));

	//			float fVal;
	//			bfcdata->ReadFloat(&fVal, 1);
	//			pPointValues[i] = fVal;
	//			Message(fVal, "fVal");
	//		}
	//		else
	//		{
	//			pPointValues[i] = fNoData;
	//		}
	//	}
	//}
	//fclose(fpParcels);
	bfcdata->Close();


	/////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////// Extract data statistics from the memory blocks /////////////////////////
	//
	nCurrent = 0;
	for (int i=0; i<NumParcels; i++ )
	{
		if (i * 100 / NumParcels > nCurrent)
		{
			nCurrent = i * 100 / NumParcels;
			fptr("Extracting Statistics...", nCurrent);
		}

		if (pCounts[i] > 0)
		{
			qsort( (void *)ppData[i], (size_t)pCounts[i], sizeof(float), fCompare );

			// assign stats...
			pMin[i] = ppData[i][0];
			pMax[i] = ppData[i][pCounts[i]-1];
			pMedian[i] = ppData[i][pCounts[i]/2];

			float fTotal = 0.0F;
			for ( int j=0; j<pCounts[i]; j++ )
			{
				fTotal += ppData[i][j];
			}

			// calculate mean
			pMean[i] = fTotal / pCounts[i];

			// calculate standard deviation
			if (pCounts[i] > 1)
			{
				float fSum = 0.0F;
				for ( int j=0; j<pCounts[i]; j++ ) 
					fSum += ((pMean[i]-ppData[i][j])*(pMean[i]-ppData[i][j]));
				pStdDev[i] = sqrt(fSum/(pCounts[i]-1));
			}
		}

		else if (pPointValues[i] != fNoData) // we have a single cell to include in the stats
		{
			pMin[i] = pMax[i] = pMean[i] = pMedian[i] = pPointValues[i];
			pStdDev[i] = 0.0F;
		}
	}


	if (InitialStep)
	{
		//
		// write the short output table to get X and Y from tiny parcels
		//
		if (!WriteOutputTable(InputPath, OutputPath, pCounts, NumParcels))
		{
			Message("Cannot Write Output Table", "ERROR");
		}
	}
	else
	{
		//
		// write the FULL output table
		//
		if (!WriteOutputTable2(InputPath, 
							   OutputPath, 
							   pCounts, pMin, pMax, pMean, pMedian, pStdDev,
							   NumParcels))
		{
			Message("Cannot Write Output Table2", "ERROR");
		}
	}

	// clean up
	if (pPointValues) delete[] pPointValues;
	if (pCounts) delete[] pCounts;
	if (pCurrentIndex) delete[] pCurrentIndex;
	if (pMin) delete[] pMin;
	if (pMax) delete[] pMax;
	if (pMean) delete[] pMean;
	if (pMedian) delete[] pMedian;
	if (pStdDev) delete[] pStdDev;
	for ( int i=0; i<NumParcels; i++ )
	{
		if (ppData[i]) delete[] ppData[i];
	}
	if (ppData) delete[] ppData;

	return(true);
}


//
// Check binary input grids for coincidence
//
bool CheckGrids(char *ParcelPath, char *DataPath)
{
	gmpath gmPath;
	EsriBinaryHeader *ParcelHeader = new EsriBinaryHeader(gmPath.ChangeExtension(ParcelPath, ".hdr"));
	EsriBinaryHeader *DataHeader = new EsriBinaryHeader(gmPath.ChangeExtension(DataPath, ".hdr"));

	if (ParcelHeader->GetCellSize() != DataHeader->GetCellSize())
	{
		Message("Error in grid size", "ERROR");
		if (ParcelHeader) delete ParcelHeader;
		if (DataHeader) delete DataHeader;
		return(false);
	}

	if (ParcelHeader->GetNumCols() != DataHeader->GetNumCols())
	{
		Message("Error in column count", "ERROR");
		if (ParcelHeader) delete ParcelHeader;
		if (DataHeader) delete DataHeader;
		return(false);
	}

	if (ParcelHeader->GetNumRows() != DataHeader->GetNumRows())
	{
		Message("Error in row count", "ERROR");
		if (ParcelHeader) delete ParcelHeader;
		if (DataHeader) delete DataHeader;
		return(false);
	}

	if (ParcelHeader->GetXllCorner() != DataHeader->GetXllCorner())
	{
		Message("Error in Xll corner", "ERROR");
		if (ParcelHeader) delete ParcelHeader;
		if (DataHeader) delete DataHeader;
		return(false);
	}

	if (ParcelHeader->GetYllCorner() != DataHeader->GetYllCorner())
	{
		Message("Error in Yll corner", "ERROR");
		if (ParcelHeader) delete ParcelHeader;
		if (DataHeader) delete DataHeader;
		return(false);
	}

	return(true);
}


//
// Count Number of Parcels in Seed Table
//
int GetParcelCount(char *InputPath)
{
	int nParcels = 0;
	char buff[1000];
	FILE *fp = fopen(InputPath, "r+t");
	fgets(buff, 1000, fp);
	while( 1 ) 
	{
		if ( NULL == fgets( buff, 1000, fp ) )
			break;
		++nParcels;
	}
	return(nParcels);
}


//
// Get the cell counts in the parcel grid
//
bool CountCells(char *ParcelPath, int *pCounts, int NumParcels)
{
	gmpath gmPath;
	EsriBinaryHeader *ParcelHeader = new EsriBinaryHeader(gmPath.ChangeExtension(ParcelPath, ".hdr"));
	int nRows = ParcelHeader->GetNumRows();
	int nCols = ParcelHeader->GetNumCols();
	float fNoData = ParcelHeader->GetNoDataValue();
	if (ParcelHeader) delete ParcelHeader;

	float *pRowData = new float [nCols];
	BinaryFileClass *bfc = new BinaryFileClass(ParcelPath);
	for ( int i=0; i<nRows; i++ )
	{
		bfc->ReadFloat(pRowData, nCols);
		for ( int j=0; j<nCols; j++ )
		{
			if (pRowData[j] != fNoData)
			{
				int nParcel = int(pRowData[j]);

				// sanity check...
				if (nParcel > NumParcels)
				{
					Message("Index exceeds NumParcels", "ERROR");
					bfc->Close();
					if (pRowData) delete[] pRowData;
					return(false);
				}

				// convert from one-based to zero-based index
				++pCounts[nParcel-1];  
			}
		}
	}
	bfc->Close();
	if (pRowData) delete[] pRowData;
	return(true);
}


//
// write short output csv table
//
bool WriteOutputTable(char *InputPath, char *OutputPath, int *pCounts, int NumParcels)
{
	char buff[1000];
	FILE *fpIn = fopen(InputPath, "r+t");
	fgets(buff, 1000, fpIn); // read header

	FILE *fpOut = fopen(OutputPath, "w+t");
	fprintf(fpOut, "ObjectID,Area_Ha,Count,X_Coord,Y_Coord\n");
	for (int i=0; i<NumParcels; i++ )
	{
		fgets(buff, 1000, fpIn);
		char *p = strtok(buff, ",\n");
		fprintf(fpOut, "%s,", p);

		p = strtok(NULL, ",\n");
		fprintf(fpOut, "%lf,", atof(p));

		fprintf(fpOut, "%d,", pCounts[i]);

		fprintf(fpOut, "%lf,%lf\n", 0.0,0.0);
	}

	if (fpIn) fclose(fpIn);
	if (fpOut) fclose(fpOut);
	return(true);
}


//
// write full output csv table
//
bool WriteOutputTable2(char *InputPath, char *OutputPath, 
	                   int *pCounts, float *pMin, float *pMax, float *pMean, float *pMedian, float *pStdDev, 
					   int NumParcels)
{
	char buff[1000];
	FILE *fpIn = fopen(InputPath, "r+t");
	fgets(buff, 1000, fpIn); // read header

	FILE *fpOut = fopen(OutputPath, "w+t");
	fprintf(fpOut, "ObjectID,Area_Ha,Count,Minimum,Maximum,Mean,Median,StdDev\n");
	for (int i=0; i<NumParcels; i++ )
	{
		fgets(buff, 1000, fpIn);
		char *p = strtok(buff, ",\n");
		fprintf(fpOut, "%s,", p);

		p = strtok(NULL, ",\n");
		fprintf(fpOut, "%lf,", atof(p));

		fprintf(fpOut, "%d,", pCounts[i]);
		fprintf(fpOut, "%f,", pMin[i]);
		fprintf(fpOut, "%f,", pMax[i]);
		fprintf(fpOut, "%f,", pMean[i]);
		fprintf(fpOut, "%f,", pMedian[i]);
		fprintf(fpOut, "%f\n", pStdDev[i]);
	}

	if (fpIn) fclose(fpIn);
	if (fpOut) fclose(fpOut);
	return(true);
}




////////////////////////////////////////////////////////////////////////////////////
// float comparison routine for Sort method
//
int fCompare( const void *arg1, const void *arg2 )
{
	float *p1 = (float *)arg1;
	float *p2 = (float *)arg2;
	if ( *p1 < *p2 )
		return(-1);
	else if ( *p1 > *p2 )
		return(1);
	else
		return(0);
}
