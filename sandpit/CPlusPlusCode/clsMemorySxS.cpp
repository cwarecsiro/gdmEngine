//
// clsMemorySxS.cpp
//
#include "stdafx.h"
#include "clsMemorySxS.h"
#include "GDMBufferSizeDefs.h"
#include "Message.h"


//
// Assumes that the SxS table has a header and has an Easting and Northing column preceeding the Site data
//
MemorySxS::MemorySxS(char *pPath)
{
	numRows = CountRows(pPath);
	if (numRows <= 0)
	{
		Valid = false;
		return;
	}

	numDataCols = CountDataColumns(pPath);
	if (numDataCols <= 0)
	{
		Valid = false;
		return;
	}

	ppSxS = new int * [numRows];
	for ( int i=0; i<numRows; i++ )
	{
		ppSxS[i] = new int [numDataCols];
	}

	pX = new double [numRows];
	pY = new double [numRows];

	PopulateFromSxS(pPath);
	Valid = true;
}


//
// Destructor
//
MemorySxS::~MemorySxS()
{
	if (IsValid())
	{
		for ( int i=0; i<numRows; i++ )
		{
			if (ppSxS[i]) delete[] ppSxS[i];
		}
		if (ppSxS) delete[] ppSxS;

		if (pX) delete[] pX;
		if (pY) delete[] pY;
	}
}


//
// Count the number of species data column in the SXS after skipping the leading X and Y columns
//
int MemorySxS::CountDataColumns(char *pPath)
{
	FILE *fp = fopen(pPath, "r+t");
	if ( NULL == fp)
	{
		Message("Cannot open SxS File", "MemorySxS::CountDataColumns");
		return(-1);
	}

	char *pHeader = new char [TABLE_ROW_BUFFSIZE];
	if (NULL == pHeader)
	{
		Message("Cannot allocate pHeader", "MemorySxS::CountDataColumns");
		if (fp) fclose(fp);
		return(-1);
	}

	int nCols = 0;
	fgets( pHeader, TABLE_ROW_BUFFSIZE, fp );

	char *p = strtok(pHeader,",\n");  // get the easting
	p = strtok(NULL,",\n");           // get the northing

	// now count the species data columns
	while(1)
	{
		p = strtok(NULL,",\n");       // get the species record
		if (NULL == p) 
			break;
		else
			++nCols;
	}
	
	if (pHeader) delete[] pHeader;
	if (fp) fclose(fp);
	return(nCols);
}


//
// Count the number of rows in the SXS after skipping the header row
//
int MemorySxS::CountRows(char *pPath)
{
	FILE *fp = fopen(pPath, "r+t");
	if ( NULL == fp)
	{
		Message("Cannot open SxS File", "MemorySxS::CountRows");
		return(-1);
	}

	char *pRow = new char [TABLE_ROW_BUFFSIZE];
	if (NULL == pRow)
	{
		Message("Cannot allocate pRow", "MemorySxS::CountRows");
		if (fp) fclose(fp);
		return(-1);
	}

	fgets(pRow, TABLE_ROW_BUFFSIZE, fp);  // skip the header
	int nRows = 0;
	while(1)
	{
		if (NULL == fgets(pRow, TABLE_ROW_BUFFSIZE, fp))
			break;

		++nRows;
	}
	
	if (pRow) delete[] pRow;
	if (fp) fclose(fp);
	return(nRows);
}


//
// Populate the species memory data from the SxS file
//
bool MemorySxS::PopulateFromSxS(char *pPath)
{
	FILE *fp = fopen(pPath, "r+t");
	if ( NULL == fp)
	{
		Message("Cannot open SxS File", "MemorySxS::PopulateSxS");
		return(false);
	}

	char *pRow = new char [TABLE_ROW_BUFFSIZE];
	if (NULL == pRow)
	{
		Message("Cannot allocate pRow", "MemorySxS::PopulateSxS");
		if (fp) fclose(fp);
		return(false);
	}

	fgets(pRow, TABLE_ROW_BUFFSIZE, fp);  // skip the header
	for ( int i=0; i<GetRowCount(); i++ )
	{
		fgets(pRow, TABLE_ROW_BUFFSIZE, fp);
		char *p = strtok(pRow, ",\n");    // get Easting
		pX[i] = atof(p);

		p = strtok(NULL, ",\n");          // get Northing
		pY[i] = atof(p);

		for ( int j=0; j<GetColCount(); j++ )
		{
			p = strtok(NULL, ",\n");
			ppSxS[i][j] = atoi(p);
		}
	}

	if (pRow) delete[] pRow;
	if (fp) fclose(fp);
	return(true);
}

