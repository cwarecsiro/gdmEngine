//
// BinaryGridTest.cpp
//
#include "stdafx.h"
#include "GdmBinLib.h"
#include "Message.h"
#include "clsDoPath.h"
#include "GDMBufferSizeDefs.h"
#include "clsEsriBinaryHeader.h"
#include "clsBinaryFileClass.h"


//
// Test that data cells in the domain and test grid are ALL coincident
//
// Assumes that the grids have exactly the same extent and cell size
//
bool TestEsriGrids(char *pDomainPath, char *TestGridPath, char *pLogFilePath, 
				   int nRows, int nCols, 
				   float fDomainNoData, float fGridNoData, 
				   double dCellSize, double dXMin, double dYMax,
				   FPTR fptr)
{
	bool DataOK = true;
	FILE *fpLog;
	bool HaveLogFile;
	if ( 0 == strcmp( "\0", pLogFilePath ) )
	{
		//Message( "No LogFilePath" );
		HaveLogFile = false;
		fpLog = NULL;
	}
	else
	{
		//Message( pLogFilePath );
		HaveLogFile = true;
		fpLog = fopen(pLogFilePath, "a+t");
	}


	//
	// Setup file paths with extensions
	//
	gmpath gmPath;
	char pDomainFlt[FILEPATH_BUFFSIZE];
	strcpy( pDomainFlt, gmPath.ChangeExtension( pDomainPath, ".flt" ) );

	char pTestGridFlt[FILEPATH_BUFFSIZE];
	strcpy( pTestGridFlt, gmPath.ChangeExtension( TestGridPath, ".flt" ) );

	//Message( pDomainFlt, "pDomainFlt" );
	//Message( pTestGridFlt, "pTestGridFlt" );

	float *pDomainRow = new float [nCols];
	float *pDataRow = new float [nCols];

	BinaryFileClass *bfc_Domain = new BinaryFileClass( pDomainFlt );
	if (!bfc_Domain->IsValid())
	{
		Message( "Cannot open domain file for READ", "TestEsriGrids" );
		if ( pDataRow ) delete[] pDataRow;
		if ( pDomainRow ) delete[] pDomainRow;
		return(false);
	}


	BinaryFileClass *bfc_Data = new BinaryFileClass( pTestGridFlt );
	if (!bfc_Data->IsValid())
	{
		Message( "Cannot open data file for READ", "TestEsriGrids" );
		if ( pDataRow ) delete[] pDataRow;
		if ( pDomainRow ) delete[] pDomainRow;
		return(false);
	}


	//
	// Check that the grids are coincident.....
	//
	int nInvalidData = 0;
	int nErrors = 0;
	int nCurrent = 0;
	char ValBuff[64];
	char feedback[FILEPATH_BUFFSIZE];
	sprintf(feedback, "Testing: %s", TestGridPath);
	fptr(feedback, nCurrent);	
	for ( int y=0; y<nRows; y++ )
	{
		if ( y * 100 / nRows > nCurrent )
		{
			nCurrent = y * 100 / nRows;
			fptr(feedback, nCurrent);	
		}

		bfc_Domain->ReadFloat( pDomainRow, nCols );
		bfc_Data->ReadFloat( pDataRow, nCols );

		for ( int x=0; x<nCols; x++ )
		{
			sprintf(ValBuff, "%lf", pDomainRow[x]);
			if (0 == strcmp(ValBuff, "-1.#IND00"))
			{
				++nInvalidData;
				fprintf( fpLog, "Invalid Domain Data (-1.#IND00) at Row: %d   Col: %d\n", y, x);
				DataOK = false;	
			}

			sprintf(ValBuff, "%lf", pDataRow[x]);
			if (0 == strcmp(ValBuff, "-1.#IND00"))
			{
				++nInvalidData;
				fprintf( fpLog, "Invalid Rest Grid Data (-1.#IND00) at Row: %d   Col: %d\n", y, x);
				DataOK = false;	
			}

			// if we have valid domain data
			if ( pDomainRow[x] != fDomainNoData )
			{
				// if the cells are NOT data coincident then 
				// increment error count and set return value to false
				if ( pDataRow[x] == fGridNoData )
				{
					++nErrors;
					if ( HaveLogFile )
					{
						if ((nErrors >= 1) && (nErrors <= 1000))
						{
							fprintf( fpLog, 
								     "Error at X: %lf   Y: %f   Row: %d   Col: %d\n", 
									 dXMin + (x * dCellSize) + (dCellSize / 2), 
									 dYMax - (y * dCellSize) - (dCellSize / 2),
									 y, x);
						}
						if (nErrors == 1000)
						{
							fprintf( fpLog, "Error count exceeds 1000\n" );
						}
					}
					DataOK = false;	
				}
			}
			// else check that the test grid has no data also
			else
			{
				// if the cells are NOT data coincident then 
				// increment error count and set return value to false
				if ( pDataRow[x] != fGridNoData )
				{
					++nErrors;
					if ( HaveLogFile )
					{
						if ((nErrors >= 1) && (nErrors <= 1000))
						{
							fprintf( fpLog, 
								     "Error at X: %lf   Y: %f   Row: %d   Col: %d\n", 
									 dXMin + (x * dCellSize) + (dCellSize / 2), 
									 dYMax - (y * dCellSize) - (dCellSize / 2),
									 y, x);
						}
						if (nErrors == 1000)
						{
							fprintf( fpLog, "Error count exceeds 1000\n" );
						}
					}
					DataOK = false;	
				}
			}
		}
	}

	if ( HaveLogFile )
	{
		fprintf( fpLog, "%d data coincidence errors\n\n", nErrors );
		fprintf( fpLog, "%d Invalid data errors\n\n", nInvalidData );
	}

	//
	// clean up
	//
	bfc_Domain->Close();
	bfc_Data->Close();
	if ( fpLog ) fclose(fpLog);
	if ( pDataRow ) delete[] pDataRow;
	if ( pDomainRow ) delete[] pDomainRow;
	return(DataOK);
}


