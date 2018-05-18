//
// ExtractPreds.cpp
//
#include "stdafx.h"
#include "GdmBinLib.h"
#include "ExtractPreds.h"
#include "myCallback.h"
#include "Message.h"
#include "ParamsW16.h"
#include "GDMBufferSizeDefs.h"
#include "clsDoPath.h"
#include "FilterSites.h"
#include "clsBinaryFileClass.h"
#include "clsEsriBinaryHeader.h"
#include "GDM_PredictorTypes.h"


//
// Main method called from Add Predictor Panel to append
// predictor values extracted from grids to an input table
// formatted Response,Weights,X0,Y0,X1,Y1
// Predictor names are prepended with S1_ and S2_ depending on their
// site pair column are also appended to the first row in the input table.
//
bool ExtractPredictors( char *pParams, bool DoBatch, FPTR fptr)
{
	//Message(pParams, "ExtractPredictors()");
	gmpath GmPath;
	char *pBuffer = new char [FILEPATH_BUFFSIZE];
	fptr("Extracting Predictors...", 0);


	//
	// Setup Output File
	//
	GetProfileString( "OUTPUT", "OutputTable", pBuffer, pParams );
	FILE *fpOut = fopen(pBuffer, "w+t");
	if ( NULL == fpOut )
	{
		if (!DoBatch)
			Message("Cannot open fpOut in ExtractPredictors", "ERROR");
		if (pBuffer) delete[] pBuffer;
		return(false);
	}


	//
	// write new full header to out table
	//
	fprintf( fpOut, "Response,Weights,X0,Y0,X1,Y1,");
	int nPreds = GetProfileInt( "PREDICTORS", "NumPredictors", pParams );

	char lpKey[32];
	for ( int i=1; i<=nPreds; i++ )
	{
		sprintf(lpKey, "PredType%d", i);
		if (PRED_TYPE_TABLE == GetProfileInt( "PREDICTORS", lpKey, pParams ))
		{
			sprintf(lpKey, "Tab0Name%d", i);
			GetProfileString( "PREDICTORS", lpKey, pBuffer, pParams );
			fprintf( fpOut, "S0_%s,", pBuffer);
		}
		else
		{
			sprintf(lpKey, "EnvGrid%d", i);
			GetProfileString( "PREDICTORS", lpKey, pBuffer, pParams );
			fprintf( fpOut, "S0_%s,", GmPath.GetName(pBuffer));
		}
	}


	for ( int i=1; i<=nPreds; i++ )
	{
		sprintf(lpKey, "PredType%d", i);
		if (PRED_TYPE_TABLE == GetProfileInt( "PREDICTORS", lpKey, pParams ))
		{
			sprintf(lpKey, "Tab1Name%d", i);
			GetProfileString( "PREDICTORS", lpKey, pBuffer, pParams );
			fprintf( fpOut, "S1_%s", pBuffer);
		}
		else
		{
			sprintf(lpKey, "EnvGrid%d", i);
			GetProfileString( "PREDICTORS", lpKey, pBuffer, pParams );
			fprintf( fpOut, "S1_%s", GmPath.GetName(pBuffer));
		}

		if ( i < nPreds )
			fprintf( fpOut, ",");
		else
			fprintf(fpOut, "\n");
	}


	//
	// Setup Input File
	//
	GetProfileString( "INPUT", "InputTable", pBuffer, pParams );
	/*FILE *fpIn = fopen(pBuffer, "r+t");
	if ( NULL == fpIn )
	{
		if (!DoBatch)
			Message("Cannot open fpIn in ExtractPredictors", "ERROR");
		if (pBuffer) delete[] pBuffer;
		return(false);
	}*/


	// open text file
	ifstream* pmyFile = new ifstream; // On the heap
	pmyFile->open( pBuffer );

	if (!pmyFile->is_open())
	{
		if (!DoBatch)
			Message("Cannot open fpIn in ExtractPredictors", "ERROR");
		if (pBuffer) delete[] pBuffer;
		return(false);
	}


	//
	// Setup the input grids
	//
	EsriBinaryHeader **ppHeaders = new EsriBinaryHeader * [nPreds];
	BinaryFileClass  **ppBfc = new BinaryFileClass * [nPreds];
	float *pGridValues_0 = new float [nPreds];
	float *pGridValues_1 = new float [nPreds];
	int *pPredType = new int [nPreds];
	double *** pppLookupTables = new double ** [nPreds];
	int *pLookupSize = new int [nPreds];

	// Setup the input tables, if any
	char pTableRow[FILEPATH_BUFFSIZE];
	FILE **ppCSVs = new FILE * [nPreds];


	int nCurrent = 0;
	fptr("Setting up predictors...", nCurrent);
	for ( int i=0; i<nPreds; i++ )
	{
		if ( i * 100 / nPreds > nCurrent )
		{
			nCurrent = i * 100 / nPreds;
			fptr("Setting up predictors...", nCurrent);
		}

		//
		// Special case for Tabular Predictors
		//
		sprintf(lpKey, "PredType%d", i+1);
		if (PRED_TYPE_TABLE == GetProfileInt( "PREDICTORS", lpKey, pParams ))
		{
			pPredType[i] = GetProfileInt( "PREDICTORS", lpKey, pParams );
			sprintf(lpKey, "TabPred%d", i+1);
			GetProfileString( "PREDICTORS", lpKey, pBuffer, pParams );
			ppCSVs[i] = fopen(pBuffer, "r+t");
			fgets(pTableRow, FILEPATH_BUFFSIZE, ppCSVs[i]); // get the header
			
			pppLookupTables[i] = NULL;
			ppHeaders[i] = NULL;
			ppBfc[i] = NULL;
			continue;
		}
		else
		{
			ppCSVs[i] = NULL;
		}


		sprintf(lpKey, "EnvGrid%d", i+1);
		GetProfileString( "PREDICTORS", lpKey, pBuffer, pParams );
		ppHeaders[i] = new EsriBinaryHeader(GmPath.ChangeExtension(pBuffer, ".hdr"));
		ppBfc[i] = new BinaryFileClass(GmPath.ChangeExtension(pBuffer, ".flt"));

		sprintf(lpKey, "PredType%d", i+1);
		pPredType[i] = GetProfileInt( "PREDICTORS", lpKey, pParams );

		if (pPredType[i] == PRED_TYPE_CATEGORICAL)
		{
			sprintf(lpKey, "CatLookup%d", i+1);
			GetProfileString( "PREDICTORS", lpKey, pBuffer, pParams );
			pppLookupTables[i] = ExtractLookupTable(pBuffer, DoBatch, &pLookupSize[i]);
			if (NULL == pppLookupTables[i])
			{
				// cleanup
				for ( int ii=0; ii<nPreds; ii++ )
				{
					if (ppHeaders[ii]) delete ppHeaders[ii];
					if (ppBfc[ii]) delete ppBfc[ii];
					if (pppLookupTables[ii]) 
					{
						for ( int j=0; j<pLookupSize[ii]; j++ )
						{
							if (pppLookupTables[ii][j]) delete[] pppLookupTables[ii][j];
						}
						delete[] pppLookupTables[ii]; 
					}
				}
				if (ppHeaders) delete[] ppHeaders;
				if (ppBfc) delete[] ppBfc;
				if (pppLookupTables)  delete pppLookupTables;
				if (pGridValues_0) delete[] pGridValues_0;
				if (pGridValues_1) delete[] pGridValues_1;
				if (pPredType) delete[] pPredType;
				if (pLookupSize) delete[] pLookupSize;
				if (pBuffer) delete[] pBuffer;
				if (ppCSVs[i]) fclose(ppCSVs[i]);
				if (ppCSVs) delete[] ppCSVs;
				return(false);
			}
		}
		else
		{
			pppLookupTables[i] = NULL;
		}

	} // for ( int i=0; i<nPreds; i++ )


	//
	// Extract the predictor data
	//
	char *pRow = new char [TABLE_ROW_BUFFSIZE]; // get header
	//fgets( pRow, TABLE_ROW_BUFFSIZE, fpIn );

	// get header
	pmyFile->getline(pRow,TABLE_ROW_BUFFSIZE);

	nCurrent = 0;
	int nThis = 0;
	char banner[100];
	sprintf(banner, "%s%d%s", "Extracting Predictors...(", nThis, " records processed)"); 
	fptr(banner, nCurrent);
	while(1)
	{
		if ( ++nThis % 1000 == 0 )   // update every 1000 records
		{
			if (++nCurrent > 100 )
				nCurrent = 1;

			sprintf(banner, "%s%d%s", "Extracting Predictors...(", nThis, " records processed)"); 
			fptr(banner, nCurrent);
		}

		// get the current record
		//if (NULL == fgets( pRow, TABLE_ROW_BUFFSIZE, fpIn ))
		//	break;

		//if (NULL ==	
		pmyFile->getline(pRow, TABLE_ROW_BUFFSIZE);
		if (pmyFile->eof())
			break;

		// extract the sites
		char *p = strtok(pRow, ",\n");   // get the response
		double dResponse = atof(p);

		p = strtok(NULL, ",\n");         // get the weight
		double dWeight = atof(p);

		p = strtok(NULL, ",\n");         // get X0
		double dX0 = atof(p);

		p = strtok(NULL, ",\n");         // get Y0
		double dY0 = atof(p);

		p = strtok(NULL, ",\n");         // get X1
		double dX1 = atof(p);

		p = strtok(NULL, ",\n");         // get Y1
		double dY1 = atof(p);

		bool fDataIsOK = true;


		//
		// This will greatly reduce the extraction time but assumes the grids are all GDM prepared. ie coincident!!!
		//
		long lOffset0 = 0L;
		long lOffset1 = 0L;
		if (pPredType[0] != PRED_TYPE_TABLE)
		{
			if ( (!SiteInBoundingRectangle(dX0,dY0,ppHeaders[0])) || (!SiteInBoundingRectangle(dX1,dY1,ppHeaders[0])))
			{
				fDataIsOK = false;
				for ( int j=0; j<nPreds; j++ )
				{
					if (ppCSVs[j]) fgets(pTableRow, FILEPATH_BUFFSIZE, ppCSVs[j]); // skip this row record
				}
				continue;
			}
			
			float fVal0 = 0.0F;
			float fVal1 = 0.0F;
			if ((!SiteHasData(dX0,dY0,&fVal0,ppHeaders[0],ppBfc[0])) || (!SiteHasData(dX1,dY1,&fVal0,ppHeaders[0],ppBfc[0])))
			{
				fDataIsOK = false;
				for ( int j=0; j<nPreds; j++ )
				{
					if (ppCSVs[j]) fgets(pTableRow, FILEPATH_BUFFSIZE, ppCSVs[j]); // skip this row record
				}
				continue;
			}
		
			lOffset0 = ((long)(fabs(ppHeaders[0]->GetMaxY() - dY0) / 
										ppHeaders[0]->GetCellSize()) * ppHeaders[0]->GetNumCols()) + 
												((long)(fabs(dX0 - ppHeaders[0]->GetMinX()) / ppHeaders[0]->GetCellSize())); 

			lOffset1 = ((long)(fabs(ppHeaders[0]->GetMaxY() - dY1) / 
										ppHeaders[0]->GetCellSize()) * ppHeaders[0]->GetNumCols()) + 
												((long)(fabs(dX1 - ppHeaders[0]->GetMinX()) / ppHeaders[0]->GetCellSize())); 
		}


		for ( int j=0; j<nPreds; j++ )
		{
			if (pPredType[j] == PRED_TYPE_GRID)
			{
				ppBfc[j]->SeekTo(lOffset0 * sizeof(float));
				ppBfc[j]->ReadFloat(&pGridValues_0[j],1);
				ppBfc[j]->SeekTo(lOffset1 * sizeof(float));
				ppBfc[j]->ReadFloat(&pGridValues_1[j],1);
			}

			//
			// if the predictor type is covariance, then we need to set the value to 0.0 for 
			// the first site in the pair and the difference in covariance for the second site in the pair
			//
			else if (pPredType[j] == PRED_TYPE_COVARIANT)
			{
				ppBfc[j]->SeekTo(lOffset0 * sizeof(float));
				ppBfc[j]->ReadFloat(&pGridValues_0[j],1);
				ppBfc[j]->SeekTo(lOffset1 * sizeof(float));
				ppBfc[j]->ReadFloat(&pGridValues_1[j],1);

				float fTmp = (float)fabs(pGridValues_0[j] - pGridValues_1[j]);

				pGridValues_0[j] = 0.0F;
				pGridValues_1[j] = fTmp;
			}

			//
			// if the predictor type is categorical, then we need to use the values
			// as one-based indices into a lookup table provided from the parameter file
			//
			else if (pPredType[j] == PRED_TYPE_CATEGORICAL)
			{
				ppBfc[j]->SeekTo(lOffset0 * sizeof(float));
				ppBfc[j]->ReadFloat(&pGridValues_0[j],1);
				ppBfc[j]->SeekTo(lOffset1 * sizeof(float));
				ppBfc[j]->ReadFloat(&pGridValues_1[j],1);

				float fTmp = (float)pppLookupTables[j][(int)pGridValues_0[j]][(int)pGridValues_1[j]];
				
				// reset categorical values as 0.0 for first site and lookup value as second site value
				pGridValues_0[j] = 0.0F;
				pGridValues_1[j] = fTmp;
			}

			//
			// if the predictor type is Tabular, then read the correct record from the comma delimited data file
			//
			else if (pPredType[j] == PRED_TYPE_TABLE)
			{
				fgets(pTableRow, FILEPATH_BUFFSIZE, ppCSVs[j]); // use this row record
				sprintf(lpKey, "Tab0Idx%d", j+1);
				int nIndex0 = GetProfileInt("PREDICTORS", lpKey, pParams);
				sprintf(lpKey, "Tab1Idx%d", j+1);
				int nIndex1 = GetProfileInt("PREDICTORS", lpKey, pParams);
						
				char *p = strtok(pTableRow, ",\n");
				if (nIndex0 == 0) 
				{
					pGridValues_0[j] = (float)atof(p);
				}
				if (nIndex1 == 0)
				{
					pGridValues_1[j] = (float)atof(p);
				}
				for ( int k=1; k<=max(nIndex0,nIndex1); k++ )
				{
					p = strtok(NULL, ",\n");
					if (k == nIndex0)
						pGridValues_0[j] = (float)atof(p);
					if (k == nIndex1)
						pGridValues_1[j] = (float)atof(p);
				}
			}

		} // for ( int j=0; j<nPreds; j++ )

		//
		// if the data is ok then we can write the current row to the output table
		//
		if (fDataIsOK)
		{
			fprintf(fpOut, "%lf,%lf,%lf,%lf,%lf,%lf,", dResponse, dWeight, dX0, dY0, dX1, dY1);
			for ( int j=0; j<nPreds; j++ )
			{
				//fprintf(fpOut, "%lf,", pGridValues_0[j]);
				fprintf(fpOut, "%0.3f,", pGridValues_0[j]);
			}
			for ( int j=0; j<nPreds; j++ )
			{
				//fprintf(fpOut, "%lf", pGridValues_1[j]);
				fprintf(fpOut, "%0.3f", pGridValues_1[j]);
				
				if ( j < nPreds-1 )
					fprintf( fpOut, "," );
				else
					fprintf( fpOut, "\n" );
			} // for ( int j=0; j<nPreds; j++ )
		} // if (fDataIsOK)
	} // while(1)


	//
	// Cleanup
	//
	for ( int i=0; i<nPreds; i++ )
	{
		if (ppHeaders[i]) delete ppHeaders[i];
		if (ppBfc[i]) delete ppBfc[i];
		if (pppLookupTables[i]) 
		{
			for ( int j=0; j<pLookupSize[i]; j++ )
			{
				if (pppLookupTables[i][j]) delete[] pppLookupTables[i][j];
			}
			delete[] pppLookupTables[i]; 
		}
		if (ppCSVs[i]) fclose(ppCSVs[i]);
	}
	if (pppLookupTables)  delete pppLookupTables;
	if (pLookupSize) delete[] pLookupSize;
	if (pPredType) delete[] pPredType;
	if (pGridValues_0) delete[] pGridValues_0;
	if (pGridValues_1) delete[] pGridValues_1;
	if (ppHeaders) delete[] ppHeaders;
	if (ppBfc) delete[] ppBfc;
	if (pRow) delete[] pRow;
	if (pBuffer) delete[] pBuffer;	
	if (ppCSVs) delete[] ppCSVs;
	//if (fpIn) fclose(fpIn);
	pmyFile->close();
	if (pmyFile) delete pmyFile;
	if (fpOut) fclose(fpOut);
	fptr("Ready", 0);
	return(true);
}



//
// Extract a comma delimited lookup table into a lookup matrix of doubles
// Assumes that the input file has both a leading names row 
// and a leading names column before the data.
//
double **ExtractLookupTable(char *pPath, bool DoBatch, int *pSize)
{
	// validate and determine the size of the lookup table
	if (!ValidateLookupFile(pPath, DoBatch, pSize))
	{
		if (!DoBatch)
			Message("Invalid Lookup Table", "ERROR");
		*pSize = 0;
		return(NULL);
	}

	int nSize = *pSize;
	if ( nSize < 1 )
	{
		if (!DoBatch)
			Message("Invalid Lookup Table Size", "ERROR");
		*pSize = 0;
		return(NULL);
	}

	double **ppLookup = new double * [nSize];
	for ( int i=0; i<nSize; i++ )
	{
		ppLookup[i] = new double [nSize];
	}

	//
	// Populate from the comma delimited text file
	//
	char *pRow = new char [TABLE_ROW_BUFFSIZE];

	FILE *fp = fopen(pPath, "r+t");
	if (NULL == fp)
	{
		if (!DoBatch)
			Message("Cannot open pPath in ExtractLookupTable", "ERROR");
		if (pRow) delete[] pRow;
		*pSize = 0;
		return(NULL);
	}

	
	// get the header row and count the columns (subtract one for the leading column)
	if (NULL == fgets( pRow, TABLE_ROW_BUFFSIZE, fp ))
	{
		if (!DoBatch)
			Message("Cannot get header row in ExtractLookupTable", "ERROR");
		if (pRow) delete[] pRow;
		*pSize = 0;
		return(NULL);
	}

	// read data
	for ( int i=0; i<nSize; i++ )
	{
		fgets( pRow, TABLE_ROW_BUFFSIZE, fp );
		char *p = strtok(pRow, ",\n"); // get leading row name

		for ( int j=0; j<nSize; j++ )
		{
			p = strtok(NULL, ",\n");
			ppLookup[i][j] = atof(p);
		}
	}

	if (pRow) delete[] pRow;
	return(ppLookup);
}



//
// Make sure that the lookup table file has the same number of row and columns
// and set pSize to the row/col count less one for the leading row/column
//
bool ValidateLookupFile(char *pPath, bool DoBatch, int *pSize)
{
	char *pRow = new char [TABLE_ROW_BUFFSIZE];

	FILE *fp = fopen(pPath, "r+t");
	if (NULL == fp)
	{
		if (!DoBatch)
			Message("Cannot open pPath in ValidateLookupFile", "ERROR");
		if (pRow) delete[] pRow;
		*pSize = 0;
		return(false);
	}

	
	// get the header row and count the columns (subtract one for the leading column)
	if (NULL == fgets( pRow, TABLE_ROW_BUFFSIZE, fp ))
	{
		if (!DoBatch)
			Message("Cannot get header row in ValidateLookupFile", "ERROR");
		if (pRow) delete[] pRow;
		*pSize = 0;
		return(false);
	}

	char *p = strtok(pRow, ",\n"); // leading row name
	int nCols = 0;
	while (NULL != (p = strtok(NULL, ",\n")))
	{
		++nCols;   // count the lookup columns
	}
	//Message(nCols, "nCols");


	//
	// now make sure that we have the same number of rows
	// and that each row has the correct number of columns
	//
	rewind(fp);
	fgets( pRow, TABLE_ROW_BUFFSIZE, fp ); // get header
	int nRows = 0;
	while(1)
	{
		if (NULL == fgets( pRow, TABLE_ROW_BUFFSIZE, fp ))
			break;

		char *p = strtok(pRow, ",\n"); // leading row name
		int nTmpCols = 0;
		while (NULL != (p = strtok(NULL, ",\n")))
		{
			++nTmpCols;   // count the lookup columns
		}

		if ( nCols != nTmpCols)
		{
			if (!DoBatch)
				Message("Row count does not equal col count in ValidateLookupFile", "ERROR");
			*pSize = 0;
			if (pRow) delete[] pRow;
			return(false);
		}

		++nRows;
	}
	
	if (nRows != nCols)
	{
		if (!DoBatch)
			Message("Row count does not equal col count in ValidateLookupFile", "ERROR");
		*pSize = 0;
		if (pRow) delete[] pRow;
		return(false);
	}

	*pSize = nRows;
	if (pRow) delete[] pRow;
	return(true);
}

