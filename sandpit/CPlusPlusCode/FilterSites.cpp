//
//
// FilterSites.cpp
//
//
#include "stdafx.h"
#include "GdmBinLib.h"
#include "FilterSites.h"
#include "GDMBufferSizeDefs.h"
#include "myCallback.h"
#include "ParamsW16.h"
#include "clsDoPath.h"
#include "ProfileStrings.h"
#include "clsBinaryFileClass.h"
#include "clsEsriBinaryHeader.h"


//
// Filter a GDM Composite thru ESRI binary export grid batch function
//
bool FilterCompositeTables(char *pParams, bool DoBatch, FPTR fptr)
{
	//Message(pParams, "pParams");
	gmpath GmPath;

	//
	// Setup the input composite file path
	//
	char *pInCompPath = new char [FILEPATH_BUFFSIZE];
	GetProfileString( "INPUT", "CompositeTable", pInCompPath, pParams );
	//Message(pInCompPath, "pInCompPath");	


	//
	// Setup the input composite file path
	//
	char *pOutCompPath = new char [FILEPATH_BUFFSIZE];
	GetProfileString( "OUTPUT", "OutDirectory", pOutCompPath, pParams );
	strcat(pOutCompPath, "\\");
	strcat(pOutCompPath, GmPath.GetName(pInCompPath));
	strcat(pOutCompPath, ".csv");
	//Message(pOutCompPath, "pOutCompPath");


	//
	// Setup the binary grids
	//
	int nGrids = GmPath.GetProfileInt("INPUT", "NumGrids", pParams);
	BinaryFileClass **ppBfc = new BinaryFileClass * [nGrids];
	EsriBinaryHeader **ppBfh = new EsriBinaryHeader * [nGrids];
	char lpKey[64];
	for ( int i=0; i<nGrids; i++ )
	{
		sprintf(lpKey, "FilterGrid%d", i+1);
		ppBfc[i] = new BinaryFileClass(GmPath.ChangeExtension(GmPath.GetProfileString("INPUT", lpKey, pParams), ".flt"));
		if (!ppBfc[i]->IsValid())
		{
			char tmp[256];
			sprintf(tmp, "Cannot open Filter Grid %s", GmPath.ChangeExtension(GmPath.GetProfileString("INPUT", lpKey, pParams), ".flt"));
			Message(tmp, "ERROR");
			for ( int j=0; j<i; j++ )
			{
				ppBfc[j]->Close();
				if (ppBfc[j]) delete ppBfc[j];
				if (ppBfh[j]) delete ppBfh[j];
			}
			if (ppBfc) delete[] ppBfc;
			if (ppBfh) delete[] ppBfh;
			if (pInCompPath) delete pInCompPath;
			if (pOutCompPath) delete pOutCompPath;
			return(false);
		}
		ppBfh[i] = new EsriBinaryHeader(GmPath.ChangeExtension(GmPath.GetProfileString("INPUT", lpKey, pParams), ".hdr"));
	}

	//
	// Open the input composite table
	//
	FILE *fpIn = fopen(pInCompPath, "r+t");
	if ( NULL == fpIn)
	{
		char tmp[256];
		sprintf(tmp, "Cannot open Input Composite File %s", pInCompPath);
		if (pInCompPath) delete pInCompPath;
		if (pOutCompPath) delete pOutCompPath;
		for ( int i=0; i<nGrids; i++ )
		{
			ppBfc[i]->Close();
			if (ppBfc[i]) delete ppBfc[i];
			if (ppBfh[i]) delete ppBfh[i];
		}
		if (ppBfc) delete[] ppBfc;
		if (ppBfh) delete[] ppBfh;
		Message(tmp, "ERROR");
		return(false);
	}

	// allocate the row buffers
	char *pInBuffer = new char [TABLE_ROW_BUFFSIZE];
	fgets(pInBuffer, TABLE_ROW_BUFFSIZE, fpIn);
	char *pOutBuffer = new char [TABLE_ROW_BUFFSIZE];
	strcpy(pOutBuffer, pInBuffer);


	//
	// Create the output filtered composite table and copy the input file header to it
	//
	FILE *fpOut = fopen(pOutCompPath, "w+t");
	if ( NULL == fpOut)
	{
		char tmp[256];
		sprintf(tmp, "Cannot create Output Composite File %s", pOutCompPath);
		if (fpIn) fclose(fpIn);
		if (pInCompPath) delete pInCompPath;
		if (pOutCompPath) delete pOutCompPath;
		for ( int i=0; i<nGrids; i++ )
		{
			ppBfc[i]->Close();
			if (ppBfc[i]) delete ppBfc[i];
			if (ppBfh[i]) delete ppBfh[i];
		}
		if (ppBfc) delete[] ppBfc;
		if (ppBfh) delete[] ppBfh;
		if (pInBuffer) delete pInBuffer;
		if (pOutBuffer) delete pOutBuffer;
		Message(tmp, "ERROR");
		return(false);
	}
	fprintf(fpOut, "%s", pOutBuffer);

	//
	// Filter the input table through the grids
	// 
	while(1)
	{
		if (NULL == fgets( pInBuffer, TABLE_ROW_BUFFSIZE, fpIn ))
			break;

		// copy the row
		strcpy(pOutBuffer, pInBuffer);

		// extract the sites
		char *p = strtok(pInBuffer, ",\n");   // get the response
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

		for (int j=0; j<nGrids; j++ )
		{
			if ( (!SiteInBoundingRectangle(dX0,dY0,ppBfh[j])) || (!SiteInBoundingRectangle(dX1,dY1,ppBfh[j])))
			{
				fDataIsOK = false;
				break;
			}

			float fVal0;
			//float fVal1;
			if ((!SiteHasData(dX0,dY0,&fVal0,ppBfh[j],ppBfc[0])) || (!SiteHasData(dX1,dY1,&fVal0,ppBfh[j],ppBfc[j])))
			{
				fDataIsOK = false;
				break;
			}
		}

		//
		// if the data is ok then we can write the current row to the output table
		//
		if (fDataIsOK)
		{
			fprintf(fpOut, "%s", pOutBuffer);
		} // if (fDataIsOK)
	} // while(1)


	//
	// Cleanup
	//
	if (fpIn) fclose(fpIn);
	if (fpOut) fclose(fpOut);
	if (pInCompPath) delete pInCompPath;
	if (pOutCompPath) delete pOutCompPath;
	for ( int i=0; i<nGrids; i++ )
	{
		ppBfc[i]->Close();
		if (ppBfc[i]) delete ppBfc[i];
		if (ppBfh[i]) delete ppBfh[i];
	}
	if (ppBfc) delete[] ppBfc;
	if (ppBfh) delete[] ppBfh;
	if (pInBuffer) delete pInBuffer;
	if (pOutBuffer) delete pOutBuffer;
	return(true);
}


//
// Filter GDM Composite Table Thru ESRI Binary Export Grids
//
bool FilterTableThruGrids(char *pParams, FPTR fptr)
{
	//Message(pParams, "pParams");
	gmpath GmPath;

	char *pRow = new char [TABLE_ROW_BUFFSIZE];
	char *pTmp = new char [TABLE_ROW_BUFFSIZE];
	char *pInPath = new char [TABLE_ROW_BUFFSIZE];
	char *pOutPath = new char [TABLE_ROW_BUFFSIZE];
	
	// open input text file
	strcpy(pInPath,GmPath.GetProfileString("INPUT_TABLE", "intable", pParams));
	//Message(pInPath, "INPUT");
	FILE *fpIn = fopen( pInPath, "r+t");
	if ( NULL == fpIn)
	{
		Message("Cannot open input table in FilterTableThruGrids", "ERROR");
		if (pRow) delete[] pRow;
		if (pTmp) delete pTmp;
		if (pInPath) delete[] pInPath;
		if (pOutPath) delete[] pOutPath;
		return(false);
	}


	// open output text file
	strcpy(pOutPath, GmPath.GetProfileString("OUTPUT_TABLE", "outtable", pParams));
	//Message(pOutPath, "OUTPUT");
	FILE *fpOut = fopen( pOutPath, "w+t");
	if ( NULL == fpOut)
	{
		Message("Cannot open output table in FilterTableThruGrids", "ERROR");
		if (fpIn) fclose(fpIn);
		if (pRow) delete[] pRow;
		if (pTmp) delete pTmp;
		if (pInPath) delete[] pInPath;
		if (pOutPath) delete[] pOutPath;
		return(false);
	}


	// setup the input filter grids
	int nGrids = GmPath.GetProfileInt("INPUT_GRIDS", "numgrids", pParams);
	BinaryFileClass **ppBfc = new BinaryFileClass * [nGrids];
	EsriBinaryHeader **ppBfh = new EsriBinaryHeader * [nGrids];
	char lpKey[64];
	for ( int i=0; i<nGrids; i++ )
	{
		sprintf(lpKey, "grid%d", i+1);
		ppBfc[i] = new BinaryFileClass(GmPath.ChangeExtension(GmPath.GetProfileString("INPUT_GRIDS", lpKey, pParams), ".flt"));
		ppBfh[i] = new EsriBinaryHeader(GmPath.ChangeExtension(GmPath.GetProfileString("INPUT_GRIDS", lpKey, pParams), ".hdr"));
	}


	// get header
	fgets( pRow, TABLE_ROW_BUFFSIZE, fpIn );

	// write header
	fprintf(fpOut, "%s", pRow);


	// filter the rows
	int nCurrent = 0;
	fptr("Filtering Table Through Grids...", nCurrent);
	int nThis = 0;
	while(1)
	{
		++nThis;
		if (0 == nThis % 100)  // increment progress bar every 100 rows...
		{
			if (++nCurrent > 100) nCurrent = 1;
			fptr("Filtering Table Through Grids...", nCurrent);
		}
		
		if ( NULL == fgets( pRow, TABLE_ROW_BUFFSIZE, fpIn ) )
			break;

		// get a copy
		strcpy(pTmp, pRow);

		// get coordinates
		char *p = strtok(pRow, ",\n");  // get response
		p = strtok(NULL, ",\n");        // get weights
		p = strtok(NULL, ",\n");        // get X0
		double dX0 = atof(p);
		p = strtok(NULL, ",\n");        // get Y0
		double dY0 = atof(p);
		p = strtok(NULL, ",\n");        // get X1
		double dX1 = atof(p);
		p = strtok(NULL, ",\n");        // get Y1
		double dY1 = atof(p);

		//////////////////////////////////////////////////////////////////////////////////////////
		// write record if All the grids have data for two sites
		if (ValidCoordinate(dX0, dY0, dX1, dY1, nGrids, ppBfc, ppBfh)) fprintf(fpOut, "%s", pTmp);
		//////////////////////////////////////////////////////////////////////////////////////////
	}
	for ( int i=0; i<nGrids; i++ )
	{
		ppBfc[i]->Close();
		if (ppBfc[i]) delete ppBfc[i];
		if (ppBfh[i]) delete ppBfh[i];
	}
	if (ppBfc) delete[] ppBfc;
	if (ppBfh) delete[] ppBfh;
	if (pRow) delete[] pRow;
	if (pTmp) delete pTmp;
	if (fpIn) fclose(fpIn);
	if (fpOut) fclose(fpOut);
	if (pInPath) delete[] pInPath;
	if (pOutPath) delete[] pOutPath;
	nCurrent = 0;
	fptr("Ready", nCurrent);
	return(true);
}



//
// Filter input table thru domain and rewrite the filtered output to 
// the Workspace with the table filename prefixed with "Filtered_"
//
bool FilterInputTable(char *pParams, FPTR fptr)
{
	if ( 0 == _stricmp("RD_SiteBySpecies", GetResponseType(pParams) ) )
	{
		return(FilterSiteTable(pParams, fptr));
	}

	else if ((0 == _stricmp("RD_SitePairPlusDistanceTable", GetResponseType(pParams))) ||
		     (0 == _stricmp("RD_SitePairPlusDistanceSumTable", GetResponseType(pParams))))
	{
		return(FilterSitePairTable(pParams, fptr));
	}

	else if (0 == _stricmp("RD_SitePlusSymmetricalTable", GetResponseType(pParams)))
	{
		return(FilterSitePairDistanceTable(pParams, fptr));
	}

	return(false);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Local Functions...
//
//
//
// Filter a SiteBySpecies table with the first two fields representing X0, Y0 thru a domain and 
// rewriting the filtered output to the Workspace with the table filename prefixed with "Filtered_"
//
bool FilterSiteTable(char *pParams, FPTR fptr)
{
	char *pOrig = new char [TABLE_ROW_BUFFSIZE];
	if ( NULL == pOrig )
	{
		Message("Cannot allocate pOrig", "FilterSiteTable");
		return(false);
	}

	char *pCopy = new char [TABLE_ROW_BUFFSIZE];
	if ( NULL == pCopy )
	{
		Message("Cannot allocate pCopy", "FilterSiteTable");
		if (pOrig) delete[] pOrig; 
		return(false);
	}

	
	// get the row count for the response table
	int nRows = CountRows(GetResponseData(pParams), fptr);

	
	// get the domain header
	char headerPath[FILEPATH_BUFFSIZE];
	sprintf(headerPath, "%s.hdr", GetPredictorDomainPath(pParams));
	EsriBinaryHeader *hdr = new EsriBinaryHeader(headerPath);


	// get the binary domain file
	char binaryPath[FILEPATH_BUFFSIZE];
	sprintf(binaryPath, "%s.flt", GetPredictorDomainPath(pParams));
	BinaryFileClass *bfc = new BinaryFileClass(binaryPath, BFC_ReadOnly);
		

	// open the text files...
	FILE *fpIn = fopen(GetResponseData(pParams), "r+t");
	if ( NULL == fpIn )
	{
		Message("Cannot open Input Table", "FilterSiteTable");
		if (pOrig) delete[] pOrig; 
		if (pCopy) delete[] pCopy; 
		return(false);
	}

	FILE *fpOut = fopen(ConstructFilteredTablePath(pParams), "w+t");
	if ( NULL == fpOut )
	{
		Message("Cannot open Output Table", "FilterSiteTable");
		if (pOrig) delete[] pOrig; 
		if (pCopy) delete[] pCopy; 
		return(false);
	}


	// get header
	fgets( pOrig, TABLE_ROW_BUFFSIZE, fpIn );
	fprintf(fpOut, "%s", pOrig);


	// filter data
	int nWritten = 0;
	int nCurrent = 0;
	bool DoSiteAggregation = DoAggregation(pParams);
	fptr("Filtering Sites through Domain...", nCurrent);
	for ( int i=0; i<nRows; i++ )
	{
		if ( i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			fptr("Filtering Sites through Domain...", nCurrent);
		}

		// read the row 
		fgets( pOrig, TABLE_ROW_BUFFSIZE, fpIn );

		// get a copy
		memcpy(pCopy, pOrig, TABLE_ROW_BUFFSIZE);

		// extract the X field 
		char *p = strtok(pOrig, ",\n");
		int nLen = (int)strlen(p);
		double dX = atof(p);

		// extract the Y field 
		p = strtok(NULL, ",\n");
		nLen += (int)strlen(p);
		double dY = atof(p);

		// check that the site has data before re-writing
		if ( SiteInBoundingRectangle(dX, dY, hdr) )
		{
			if (SiteHasData(dX, dY, hdr, bfc))
			{
				if (DoSiteAggregation)
				{
					// get to the first species column (strlen of X and Y fields plus 2 commas)
					nLen += 2;
					fprintf(fpOut, "%lf,%lf,%s",  AggregateX(dX, hdr), AggregateY(dY, hdr), &pCopy[nLen]);
				}
				else
				{
					fprintf(fpOut, "%s", pCopy);
				}
				++nWritten;
			}
		}
	}
	//sprintf(pOrig, "Original: %d   Written: %d", nRows, nWritten);
	//Message(pOrig, "FilterSiteTable");
	fptr("Filtering Sites through Domain...", 0);
	if (hdr) delete hdr;
	if (bfc) delete bfc;
	if (pOrig) delete[] pOrig;
	if (pCopy) delete[] pCopy;
	if ( fpIn) fclose(fpIn);
	if ( fpOut) fclose(fpOut);
	SetFilteredResponseDataPath(pParams, ConstructFilteredTablePath(pParams));
	return(true);
}



//
// Filter a table with the first four fields representing X0,Y0,X1,Y1 thru a domain and
// rewriting the filtered output to the Workspace with the table filename prefixed with "Filtered_"
//
bool FilterSitePairTable(char *pParams, FPTR fptr)
{
	char *pOrig = new char [TABLE_ROW_BUFFSIZE];
	if ( NULL == pOrig )
	{
		Message("Cannot allocate pOrig", "FilterSitePairTable");
		return(false);
	}

	char *pCopy = new char [TABLE_ROW_BUFFSIZE];
	if ( NULL == pCopy )
	{
		Message("Cannot allocate pCopy", "FilterSitePairTable");
		if (pOrig) delete[] pOrig; 
		return(false);
	}

	
	// get the row count for the response table
	int nRows = CountRows(GetResponseData(pParams), fptr);

	
	// get the domain header
	char headerPath[FILEPATH_BUFFSIZE];
	sprintf(headerPath, "%s.hdr", GetPredictorDomainPath(pParams));
	EsriBinaryHeader *hdr = new EsriBinaryHeader(headerPath);


	// get the binary domain file
	char binaryPath[FILEPATH_BUFFSIZE];
	sprintf(binaryPath, "%s.flt", GetPredictorDomainPath(pParams));
	BinaryFileClass *bfc = new BinaryFileClass(binaryPath, BFC_ReadOnly);
		

	// open the text files...
	FILE *fpIn = fopen(GetResponseData(pParams), "r+t");
	if ( NULL == fpIn )
	{
		Message("Cannot open Input Table", "FilterSitePairTable");
		if (pOrig) delete[] pOrig; 
		if (pCopy) delete[] pCopy; 
		return(false);
	}

	FILE *fpOut = fopen(ConstructFilteredTablePath(pParams), "w+t");
	if ( NULL == fpOut )
	{
		Message("Cannot open Output Table", "FilterSitePairTable");
		if (pOrig) delete[] pOrig; 
		if (pCopy) delete[] pCopy; 
		return(false);
	}


	// get header
	fgets( pOrig, TABLE_ROW_BUFFSIZE, fpIn );
	fprintf(fpOut, "%s", pOrig);


	// filter data
	int nWritten = 0;
	int nCurrent = 0;
	bool DoSiteAggregation = DoAggregation(pParams);
	fptr("Filtering Sites through Domain...", nCurrent);
	for ( int i=0; i<nRows; i++ )
	{
		if ( i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			fptr("Filtering Sites through Domain...", nCurrent);
		}

		// read the row 
		fgets( pOrig, TABLE_ROW_BUFFSIZE, fpIn );

		// get a copy
		memcpy(pCopy, pOrig, TABLE_ROW_BUFFSIZE);

		// extract the X0 field 
		char *p = strtok(pOrig, ",\n");
		int nLen = (int)strlen(p);
		double dX0 = atof(p);

		// extract the Y0 field 
		p = strtok(NULL, ",\n");
		nLen += (int)strlen(p);
		double dY0 = atof(p);

		// extract the X1 field 
		p = strtok(NULL, ",\n");
		nLen += (int)strlen(p);
		double dX1 = atof(p);

		// extract the Y1 field 
		p = strtok(NULL, ",\n");
		nLen += (int)strlen(p);
		double dY1 = atof(p);

		// check that the site has data before re-writing
		if ((SiteInBoundingRectangle(dX0, dY0, hdr)) && (SiteInBoundingRectangle(dX1, dY1, hdr)))
		{
			if ((SiteHasData(dX0, dY0, hdr, bfc)) && (SiteHasData(dX1, dY1, hdr, bfc)))
			{
				if (DoSiteAggregation)
				{					
					// get to the first species column (strlen of X0, Y0, X1 and Y1 fields plus 4 commas)
					nLen += 4;
					fprintf(fpOut, 
						    "%lf,%lf,%lf,%lf,%s",  
							AggregateX(dX0, hdr), 
							AggregateY(dY0, hdr), 
							AggregateX(dX1, hdr), 
							AggregateY(dY1, hdr), 
							&pCopy[nLen]);
				}
				else
				{
					fprintf(fpOut, "%s", pCopy);
				}
				++nWritten;
			}
		}
	}
	//sprintf(pOrig, "Original: %d   Written: %d", nRows, nWritten);
	//Message(pOrig, "FilterSiteTable");
	fptr("Filtering Sites through Domain...", 0);
	if (hdr) delete hdr;
	if (bfc) delete bfc;
	if (pOrig) delete[] pOrig;
	if (pCopy) delete[] pCopy;
	if ( fpIn) fclose(fpIn);
	if ( fpOut) fclose(fpOut);
	SetFilteredResponseDataPath(pParams, ConstructFilteredTablePath(pParams));
	return(true);
}



//
// Filter a SitePairPlusDistanceTable table thru a domain and
// rewriting the filtered output to the Workspace with the table filename prefixed with "Filtered_"
// NOTE that given that there is a distance table attached to the site pair columns, if a site fails
// to have data then it's corresponding column must also be deleted.
//
bool FilterSitePairDistanceTable(char *pParams, FPTR fptr)
{
	char *pOrig = new char [TABLE_ROW_BUFFSIZE];
	if ( NULL == pOrig )
	{
		Message("Cannot allocate pOrig", "FilterSitePairTable");
		return(false);
	}

	char *pCopy = new char [TABLE_ROW_BUFFSIZE];
	if ( NULL == pCopy )
	{
		Message("Cannot allocate pCopy", "FilterSitePairTable");
		if (pOrig) delete[] pOrig; 
		return(false);
	}

	char *pHeader = new char [TABLE_ROW_BUFFSIZE];
	if ( NULL == pHeader )
	{
		Message("Cannot allocate pHeader", "FilterSitePairTable");
		if (pOrig) delete[] pOrig; 
		if (pCopy) delete[] pCopy; 
		return(false);
	}

		
	// get the row count for the response table
	int nRows = CountRows(GetResponseData(pParams), fptr);


	// get the ToUse vector setup
	bool *pToUse = new bool [nRows];

	// get the X and Y vectors setup
	double *pX = new double [nRows];
	double *pY = new double [nRows];

	// setup the distance matrix 
	double **ppMatrix = new double * [nRows];
	for ( int i=0; i<nRows; i++ )
	{
		ppMatrix[i] = new double [nRows];
	}

	
	// get the domain header
	char headerPath[FILEPATH_BUFFSIZE];
	sprintf(headerPath, "%s.hdr", GetPredictorDomainPath(pParams));
	EsriBinaryHeader *hdr = new EsriBinaryHeader(headerPath);


	// get the binary domain file
	char binaryPath[FILEPATH_BUFFSIZE];
	sprintf(binaryPath, "%s.flt", GetPredictorDomainPath(pParams));
	BinaryFileClass *bfc = new BinaryFileClass(binaryPath, BFC_ReadOnly);
		

	// open the text files...
	FILE *fpIn = fopen(GetResponseData(pParams), "r+t");
	if ( NULL == fpIn )
	{
		Message("Cannot open Input Table", "FilterSitePairTable");
		if (pOrig) delete[] pOrig; 
		if (pCopy) delete[] pCopy; 
		if (pHeader) delete[] pHeader;
		if (pToUse) delete[] pToUse;
		return(false);
	}

	FILE *fpOut = fopen(ConstructFilteredTablePath(pParams), "w+t");
	if ( NULL == fpOut )
	{
		Message("Cannot open Output Table", "FilterSitePairTable");
		if (pOrig) delete[] pOrig; 
		if (pCopy) delete[] pCopy; 
		if (pHeader) delete[] pHeader;
		if (pToUse) delete[] pToUse;
		return(false);
	}


	// get header
	fgets( pHeader, TABLE_ROW_BUFFSIZE, fpIn );


	// filter data
	int nWritten = 0;
	int nCurrent = 0;
	bool DoSiteAggregation = DoAggregation(pParams);
	fptr("Filtering Sites through Domain...", nCurrent);
	for ( int i=0; i<nRows; i++ )
	{
		if ( i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			fptr("Filtering Sites through Domain...", nCurrent);
		}

		// read the row 
		fgets( pOrig, TABLE_ROW_BUFFSIZE, fpIn );

		// get a copy
		memcpy(pCopy, pOrig, TABLE_ROW_BUFFSIZE);

		// extract the X field 
		char *p = strtok(pOrig, ",\n");
		pX[i] = atof(p);

		// extract the Y field 
		p = strtok(NULL, ",\n");
		pY[i] = atof(p);

		// check that the site has data
		pToUse[i] = SiteInBoundingRectangle(pX[i], pY[i], hdr) && SiteHasData(pX[i], pY[i], hdr, bfc);

		// get the current row into the matrix
		for ( int j=0; j<nRows; j++ )
		{
			p = strtok(NULL, ",\n");
			ppMatrix[i][j] = atof(p);
		}
	}


	// find the last index of a TRUE ToUse index
	int nLastIndex;
	for ( int i=0; i<nRows; i++ )
	{
		if ( pToUse[i] ) nLastIndex = i;
	}


	// rewrite the header
	char *p = strtok(pHeader, ",\n");
	fprintf(fpOut, "%s,", p);
	p = strtok(NULL, ",\n");
	fprintf(fpOut, "%s,", p);
	for ( int i=0; i<=nLastIndex; i++ )
	{
		p = strtok(NULL, ",\n");
		if ( pToUse[i] )
		{
			fprintf(fpOut, "%s", p);
			if (i < nLastIndex) 
				fprintf(fpOut, ",");
			else
				fprintf(fpOut, "\n");
		}
	}	


	// now use the ToUse elements to recreate the site plus distance file
	for ( int i=0; i<=nLastIndex; i++ )
	{
		if ( pToUse[i] )
		{
			fprintf(fpOut, "%lf,%lf,", pX[i], pY[i]);
			for ( int j=0; j<=nLastIndex; j++ )
			{
				if ( pToUse[j] ) 
				{
					fprintf(fpOut, "%lf", ppMatrix[i][j]);
					if ( j < nLastIndex )
						fprintf( fpOut, "," );
					else
						fprintf( fpOut, "\n" );
				}
			}
			++nWritten;
		}
	}


	// clean up
	//sprintf(pOrig, "Original: %d   Written: %d", nRows, nWritten);
	//Message(pOrig, "FilterSiteTable");
	fptr("Filtering Sites through Domain...", 0);
	if (hdr) delete hdr;
	if (bfc) delete bfc;
	if (pOrig) delete[] pOrig;
	if (pCopy) delete[] pCopy;
	if ( fpIn) fclose(fpIn);
	if ( fpOut) fclose(fpOut);
	for ( int i=0; i<nRows; i++ ) if (ppMatrix[i]) delete[] ppMatrix[i];
	if (ppMatrix) delete[] ppMatrix;
	if (pToUse) delete[] pToUse;
	if (pX) delete[] pX;
	if (pY) delete[] pY;
	SetFilteredResponseDataPath(pParams, ConstructFilteredTablePath(pParams));
	return(true);
}



//
// Create a filepath from concatenation of the Workspace Path + \\Filter_ + ResponseFilePath
//
char *ConstructFilteredTablePath(char *pParams)
{
	gmpath myPath;
	char *pOut = new char[FILEPATH_BUFFSIZE];
	strcpy(pOut, GetWorkspacePath(pParams));
	strcat(pOut, "\\Filtered_");
	strcat(pOut, myPath.GetName(GetResponseData(pParams)));
	strcat(pOut, myPath.GetExtension(GetResponseData(pParams)));
	return(pOut);
}


//
// Using the 64bit fstream class, count the number of data rows excluding the header in a text file
//
int CountRows64(char *pPath, FPTR fptr)
{
	int nRows = 0;
	char *p = new char [TABLE_ROW_BUFFSIZE];
	if ( NULL == p )
	{
		Message("Cannot allocate p in CountRows", "ERROR");
		return(nRows);
	}

	// open text file
	ifstream* pmyFile = new ifstream; // On the heap
	pmyFile->open( pPath );

	if (!pmyFile->is_open())
	{
		Message("Cannot open pPath in CountRows", "ERROR");
		if (p) delete[] p;
		return(nRows);
	}
	
	// get header
	pmyFile->getline(p,TABLE_ROW_BUFFSIZE);

	int nCurrent = 0;
	fptr("Counting Rows...", nCurrent);
	int nThis = 0;
	while(!pmyFile->eof())
	{
		++nThis;
		if (0 == nThis % 50000)  // increment progress bar every 50000 rows...
		{
			if (++nCurrent > 100) nCurrent = 1;
			fptr("Counting Rows...", nCurrent);
		}
		
		pmyFile->getline(p,TABLE_ROW_BUFFSIZE);

		// make sure that there is some data in the line
		if (strlen(p) > 4)
			++nRows;
	}

	pmyFile->close();
	if (p) delete[] p;	
	return(nRows);
}


//
// Count the number of data rows excluding the header in a text file
//
int CountRows(char *pPath, FPTR fptr)
{
	int nRows = 0;
	char *p = new char [TABLE_ROW_BUFFSIZE];
	if ( NULL == p )
	{
		Message("Cannot allocate p in CountRows", "ERROR");
		return(nRows);
	}
	
	// open text file
	FILE *fp = fopen(pPath, "r+t");
	if ( NULL == fp)
	{
		Message("Cannot open pPath in CountRows", "ERROR");
		if (p) delete[] p;
		return(nRows);
	}

	// get header
	fgets( p, TABLE_ROW_BUFFSIZE, fp );

	// count rows
	int nCurrent = 0;
	fptr("Counting Rows...", nCurrent);
	int nThis = 0;
	while(1)
	{
		++nThis;
		if (0 == nThis % 50000)  // increment progress bar every 50000 rows...
		{
			if (++nCurrent > 100) nCurrent = 1;
			fptr("Counting Rows...", nCurrent);
		}
		
		if ( NULL == fgets( p, TABLE_ROW_BUFFSIZE, fp ) )
			break;

		++nRows;
	}
	if (fp) fclose(fp);
	if (p) delete[] p;	
	return(nRows);
}



//
// Using the 64bit fstream class, count the number of comma delimited columns in a text file
//
int CountColumns64(char *pPath, FPTR fptr)
{
	int nColumns = 0;
	char *pHeader = new char [TABLE_ROW_BUFFSIZE];
	if ( NULL == pHeader )
	{
		Message("Cannot allocate pHeader in CountColumns", "ERROR");
		return(nColumns);
	}
	
	ifstream* pmyFile = new ifstream; // On the heap
	pmyFile->open( pPath );

	if (!pmyFile->is_open())
	{
		Message("Cannot open pPath in CountColumns64", "ERROR");
		if (pHeader) delete[] pHeader;
		return(nColumns);
	}

	// get header and close file
	pmyFile->getline(pHeader,TABLE_ROW_BUFFSIZE);
	pmyFile->close();

	// count columns
	char *p = strtok(pHeader, ",\n");
	while(p)
	{
		++nColumns;
		p = strtok(NULL, ",\n");
	}

	if (pHeader) delete[] pHeader;
	return(nColumns);
}



//
// Count the number of comma delimited columns in a text file
//
int CountColumns(char *pPath, FPTR fptr)
{
	int nColumns = 0;
	char *pHeader = new char [TABLE_ROW_BUFFSIZE];
	if ( NULL == pHeader )
	{
		Message("Cannot allocate pHeader in CountColumns", "ERROR");
		return(nColumns);
	}
	
	// open text file
	FILE *fp = fopen(pPath, "r+t");
	if ( NULL == fp)
	{
		Message("Cannot open pPath in CountColumns", "ERROR");
		if (pHeader) delete[] pHeader;
		return(nColumns);
	}

	// get header and close file
	fgets( pHeader, TABLE_ROW_BUFFSIZE, fp );
	if (fp) fclose(fp);

	// count columns
	char *p = strtok(pHeader, ",\n");
	while(p)
	{
		++nColumns;
		p = strtok(NULL, ",\n");
	}

	if (pHeader) delete[] pHeader;
	return(nColumns);
}



//
// Return TRUE if the Site inside the Grid has Data
//
bool SiteHasData(double dX, double dY, EsriBinaryHeader *hdr, BinaryFileClass *bfc)
{
	long lOffset = ((long)(fabs(hdr->GetMaxY() - dY) / hdr->GetCellSize()) * hdr->GetNumCols()) + 
		           ((long)(fabs(dX - hdr->GetMinX()) / hdr->GetCellSize())); 

	bfc->SeekTo(lOffset*sizeof(float));
	float fVal;
	bfc->ReadFloat(&fVal,1);
	return(fVal != hdr->GetNoDataValue());
}


//
// Return TRUE if the Site inside the Grid has Data and extract the value
//
bool SiteHasData(double dX, double dY, float *pVal, EsriBinaryHeader *hdr, BinaryFileClass *bfc)
{
	long lOffset = ((long)(fabs(hdr->GetMaxY() - dY) / hdr->GetCellSize()) * hdr->GetNumCols()) + 
		           ((long)(fabs(dX - hdr->GetMinX()) / hdr->GetCellSize())); 

	bfc->SeekTo(lOffset*sizeof(float));
	//float fVal;
	bfc->ReadFloat(pVal,1);
	return(*pVal != hdr->GetNoDataValue());
}



//
// Return TRUE is the site coordinate is inside the Binary Grids Extent
//
bool SiteInBoundingRectangle(double dX, double dY, EsriBinaryHeader *hdr)
{
	return( (dX >= hdr->GetMinX()) && (dX < hdr->GetMaxX()) && (dY >= hdr->GetMinY()) && (dY < hdr->GetMaxY()) );
}



//
// Aggregate dX to the center of the enclosing gridcell
//
double AggregateX(double dX, EsriBinaryHeader *hdr)
{
	long lOffset = long((dX - hdr->GetMinX()) / hdr->GetCellSize()); 
	return(hdr->GetMinX() + (lOffset * hdr->GetCellSize()) + (hdr->GetCellSize() / 2));
}


//
// Aggregate dY to the center of the enclosing gridcell
//
double AggregateY(double dY, EsriBinaryHeader *hdr)
{
	long lOffset = long((hdr->GetMaxY() - dY) / hdr->GetCellSize()); 
	return(hdr->GetMaxY() - (lOffset * hdr->GetCellSize()) - (hdr->GetCellSize() / 2));
}


//
// return true if the is valid data for BOTH sites in ALL the grids
//
bool ValidCoordinate(double dX0, double dY0, double dX1, double dY1, 
	                 int nGrids, BinaryFileClass **ppBfc, EsriBinaryHeader **ppBfh)
{
	long lOffset;
	float fVal;

	for ( int i=0; i<nGrids; i++ )
	{
		//
		// make sure that the coordinates are in the grid
		//
		if (dX0 < ppBfh[i]->GetMinX()) return(false);
		if (dX1 < ppBfh[i]->GetMinX()) return(false);
		if (dY0 < ppBfh[i]->GetMinY()) return(false);
		if (dY1 < ppBfh[i]->GetMinY()) return(false);
		if (dX0 > ppBfh[i]->GetMaxX()) return(false);
		if (dX1 > ppBfh[i]->GetMaxX()) return(false);
		if (dY0 > ppBfh[i]->GetMaxY()) return(false);
		if (dY1 > ppBfh[i]->GetMaxY()) return(false);

		//
		// make sure that the site coordinates have valid data
		//
		lOffset = lGetCoordinateOffset(ppBfh[i], dX0, dY0);
		ppBfc[i]->SeekTo(lOffset * sizeof(float));
		ppBfc[i]->ReadFloat(&fVal, 1);
		if (ppBfh[i]->GetNoDataValue() == fVal) return(false);

		lOffset = lGetCoordinateOffset(ppBfh[i], dX1, dY1);
		ppBfc[i]->SeekTo(lOffset * sizeof(float));
		ppBfc[i]->ReadFloat(&fVal, 1);
		if (ppBfh[i]->GetNoDataValue() == fVal) return(false);
	}
	return(true);
}


//
// return linear cellwise offset of coordinate in binary file
//
long lGetCoordinateOffset(EsriBinaryHeader *header, double dX, double dY)
{
	return(((int)((header->GetMaxY() - dY)/header->GetCellSize())*header->GetNumCols()) +  
													(int)((dX - header->GetMinX())/header->GetCellSize()));
}
