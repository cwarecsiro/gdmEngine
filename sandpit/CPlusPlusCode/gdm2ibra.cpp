//
// gdm2ibra.cpp
//
#include "stdafx.h"
#include "GdmBinLib.h"
#include "gdm2ibra.h"
#include "Message.h"
#include "clsDoPath.h"

#define GDM_CLASSES 100

bool CreateReport(char *ModelPath, 
	              char *IbraPath, 
				  char *LLSPath,
	              char *OutPath, 
				  char *IbraLookupPath,
				  int LLSIndex, 
				  bool CountCells, 
				  FPTR fptr)
{
	//Message(ModelPath, "ModelPath");
	//Message(IbraPath, "IbraPath");
	//Message(LLSPath, "LLSPath");
	//Message(OutPath, "OutPath");
	//Message(IbraLookupPath, "IbraLookupPath");
	//Message(LLSIndex, "LLSIndex");
	//Message(CountCells, "CountCells");
	//return(true);
	
	//
	// setup the data table
	//
	int nMaxIbra = GetIbraLookupMaxIndex(IbraLookupPath);
	int nMaxCols = nMaxIbra+1; // add a column for the class codes ([0])

	int **ppReportData = new int * [GDM_CLASSES];
	for ( int i=0; i<GDM_CLASSES; i++ )
	{
		fptr("Creating Data Table", i);

		ppReportData[i] = new int [nMaxCols];  
		for ( int j=0; j<nMaxCols; j++ )
		{
			ppReportData[i][j] = 0;
		}
	}

	// initialise the class codes
	for ( int i=0; i<GDM_CLASSES; i++ )
	{
		ppReportData[i][0] = i+1;
	}


	//
	// get the IBRA lookup strings
	//
	char **ppLookup = GetLookupStrings(IbraLookupPath, nMaxIbra);
	if (NULL == ppLookup)
	{
		Message("Cannot GetLookupStrings", "CreateReport");
	}


	//
	// setup the input grids
	//
	gmpath gmPath;
	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(ModelPath, ".hdr"));
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	float fNoData = header->GetNoDataValue();

	BinaryFileClass *bfcModel = new BinaryFileClass(ModelPath, BFC_ReadOnly);
	if ( !bfcModel->IsValid() )
	{
		Message( "Cannot open Model grid for READ", "CreateReport" );
		return(false);
	}

	BinaryFileClass *bfcIbra = new BinaryFileClass(IbraPath, BFC_ReadOnly);
	if ( !bfcIbra->IsValid() )
	{
		Message( "Cannot open IBRA grid for READ", "CreateReport" );
		return(false);
	}

	BinaryFileClass *bfcLLS = new BinaryFileClass(LLSPath, BFC_ReadOnly);
	if ( !bfcLLS->IsValid() )
	{
		Message( "Cannot open LLS grid for READ", "CreateReport" );
		return(false);
	}


	//
	// allocate the grid row vectors
	//
	float *pModelRow = new float [nCols];
	float *pIbraRow = new float [nCols];
	float *pLLSRow = new float [nCols];

	int *pTallies = new int [nMaxIbra];
	for ( int i=0; i<nMaxIbra; i++ )
	{
		pTallies[i] = 0;
	}

	int nCurrent = 0;
	fptr("Cross Tabulating...", nCurrent);
	for ( int i=0; i<nRows; i++ )
	{
		if (i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			fptr("Cross Tabulating...", nCurrent);
		}

		// read a row of grid data
		bfcModel->ReadFloat(pModelRow, nCols);
		bfcIbra->ReadFloat(pIbraRow, nCols);
		bfcLLS->ReadFloat(pLLSRow, nCols);

		for ( int j=0; j<nCols; j++ )
		{
			if ((pModelRow[j] != fNoData) && (pIbraRow[j] != fNoData) && (pLLSRow[j] != fNoData))
			{
				if (pLLSRow[j] == LLSIndex)  // we have the Local Land Service we asked for
				{
					int nClass = (int)pModelRow[j];
					int nIbraClass = (int)(pIbraRow[j]);
					
					pTallies[nIbraClass-1] = 1;
					ppReportData[nClass-1][nIbraClass] += 1;
				}
			} // if (pModelRow[j] != fNoData)
		} // for ( int j=0; j<nCols; j++ )
	} // for ( int i=0; i<nRows; i++ )


	//
	// write the table
	//
	WriteOutput(OutPath, nMaxIbra+1, ppReportData, ppLookup, pTallies, CountCells, header->GetCellSize(), fptr);


	//
	// cleanup
	//
	for (int i=0; i<nMaxIbra; i++ )
	{
		if (ppLookup[i]) delete[] ppLookup[i];
	}
	if (ppLookup) delete[] ppLookup;
	bfcModel->Close();
	bfcIbra->Close();
	bfcLLS->Close();
	for ( int i=0; i<GDM_CLASSES; i++ )
	{
		if (ppReportData[i]) delete[] ppReportData[i];
	}
	if (ppReportData) delete[] ppReportData;
	if (header) delete header;
	if (pModelRow) delete[] pModelRow;
	if (pIbraRow) delete[] pIbraRow;
	if (pLLSRow) delete[] pLLSRow;
	if (pTallies) delete[] pTallies;

	return(true);
}


//
// return the maximum one based index from the IBRA sub-region lookup table
//
int GetIbraLookupMaxIndex(char *IbraPath)
{
	int nMaxVal = 0;
	FILE *fp = fopen(IbraPath, "r+t");
	if ( NULL == fp )
	{
		Message("Cannot open IBRA Lookup Table", "GetIbraLookupMaxIndex");
		return(nMaxVal);
	}

	char *pRow = new char [512];
	// get header
	fgets( pRow, 512, fp );

	while(1)
	{
		if (NULL == fgets( pRow, 512, fp ))
			break;

		++nMaxVal;
	}
	if (pRow) delete[] pRow;
	fclose(fp);
	return(nMaxVal);
}


//
// extract the IBRA lookup strings from the lookup table
//
char **GetLookupStrings(char *IbraPath, int nMaxIbra)
{
	FILE *fp = fopen(IbraPath, "r+t");
	if ( NULL == fp )
	{
		Message("Cannot open IBRA Lookup Table", "GetLookupStrings");
		return(NULL);
	}


	//char len[64];

	char buff[256];
	char **ppLookup = new char * [nMaxIbra];
	for (int i=0; i<nMaxIbra; i++)
	{
		ppLookup[i] = new char [256];
	}

	fgets( buff, 256, fp ); // get header
	while(1)
	{
		if (NULL == fgets( buff, 256, fp ))
			break;

		char *p = strtok(buff, ",\n");   // get index
		int nIndex = atoi(p)-1;
		p = strtok(NULL, ",\n");         // get string
		strcpy(ppLookup[nIndex], p);
	}
	fclose(fp);
	return(ppLookup);
}


//
// write the output table
//
void WriteOutput(char *OutPath, 
	             int nCols, 
				 int **ppData, 
				 char **ppLookup, 
				 int *pTallies, 
				 bool CountCells, 
				 double CellSize,
				 FPTR fptr)
{
	FILE *fp = fopen(OutPath, "w+t");
	if ( NULL == fp )
	{
		Message("Cannot create Output file", "CreateReport");
		return;
	}

	fptr("Writing output table...", 0);	

	// write the header
	fprintf(fp, "GDM_Class");
	for ( int i=0; i<nCols-1; i++ )
	{
		if (pTallies[i] > 0)
			fprintf(fp, ",%s", ppLookup[i]);
	}
	fprintf(fp, "\n");


	// we may be wrting in hectares so calculate the factor 
	double dFactor = 1.0;
	if (CellSize < 1.0)  // probably lat/longs in degrees
	{
		dFactor = (CellSize / 0.001) * (CellSize / 0.001);
	}
	else // probably UTM in meters
	{
		dFactor = (CellSize / 100.0) * (CellSize / 100.0);
	}

	// write the data
	for ( int i=0; i<GDM_CLASSES; i++ )
	{
		fptr("Writing output table...", i);

		fprintf(fp, "%d", ppData[i][0]);
		for ( int j=1; j<nCols; j++ )
		{
			if (pTallies[j-1] > 0)
			{
				if (CountCells)
				{
					fprintf(fp, ",%d", ppData[i][j]);
				}
				else  // convert to hectares
				{
					fprintf(fp, ",%1.2lf", ppData[i][j] * dFactor);
				}
			}

			if ( j==nCols-1 )
				fprintf(fp, "\n");
		}
	}
	fclose(fp);
}
