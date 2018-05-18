//
// gdmUtilities.cpp
//
#include "stdafx.h"
#include "gdmUtilities.h"
#include "GDMBufferSizeDefs.h"
#include "Message.h"
#include "ParamsW16.h"


int CountRowsInResponse(char *pParams, bool UseFiltered, FPTR fptr)
{
	int nRows = 0;

	FILE *fp = NULL;
	if ( UseFiltered )
	{
		fp = fopen(GetFilteredResponseData(pParams), "r+t");
		//Message(GetFilteredResponseData(pParams), "UseFiltered = True");
	}
	else
	{
		fp = fopen(GetResponseData(pParams), "r+t");
		//Message(GetResponseData(pParams), "UseFiltered = False");
	}

	if (NULL == fp)
	{
		Message("Cannot open fp in CountRows", "ERROR");
		return(nRows);
	}

	char *pRow = new char [TABLE_ROW_BUFFSIZE];
	if (NULL == pRow)
	{
		Message("Cannot allocate pRow in CountRows", "ERROR");
		if (fp) fclose(fp);
		return(nRows);
	}

	fgets(pRow, TABLE_ROW_BUFFSIZE, fp);  // skip the header
	int nCurrent = 0;
	fptr("Counting Rows...", nCurrent);
	while(1)
	{
		if (0 == nRows % 100)  // increment progress bar every 100 rows...
		{
			if (++nCurrent > 100) nCurrent = 1;
			fptr("Counting Rows...", nCurrent);
		}

		if (NULL == fgets(pRow, TABLE_ROW_BUFFSIZE, fp))
			break;

		++nRows;
	}
	
	if (pRow) delete[] pRow;
	if (fp) fclose(fp);
	return(nRows);
}


