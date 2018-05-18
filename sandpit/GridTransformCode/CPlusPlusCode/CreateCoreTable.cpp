//
// CreateCoreTable.cpp
//
#include "stdafx.h"
#include "GdmBinLib.h"
#include "CreateCoreTable.h"
#include "Message.h"
#include "GDMBufferSizeDefs.h"
#include "clsMemorySxS.h"
#include "clsMemorySxSDouble.h"
#include "DistanceCalcs.h"
#include "ConvertToSxS.h"
#include "clsDoPath.h"
#include <mbstring.h>
#include <string.h>

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	WORKING WITH COMPOSITE TABLE DATA AS THE INPUT...
//
//
// Create a GDM composite file from an existing GDM composite file.
// (Only copies the first 6 columns).
//
bool CreateCompositeFromComposite(char *InputPath, char *OutputPath, FPTR fptr)
{
	// sanity check there there are at least six columns
	if (!GotEnoughColsInComposite(InputPath))
	{
		Message("Not enough columns in Input Table in CreateCompositeFromComposite", "ERROR");
		return(false);
	}

	int nCurrent = 0;
	fptr("Creating Composite Table From Composite Table...", nCurrent);

	FILE *fpIn = fopen( InputPath, "r+t");
	if ( NULL == fpIn)
	{
		Message("Cannot open input table in CreateCompositeFromComposite", "ERROR");
		return(false);
	}

	// open output text file
	FILE *fpOut = fopen( OutputPath, "w+t");
	if ( NULL == fpOut)
	{
		Message("Cannot open output table in CreateCompositeFromComposite", "ERROR");
		if (fpIn) fclose(fpIn);
		return(false);
	}

	// get header
	char *pRow = new char [TABLE_ROW_BUFFSIZE];
	fgets( pRow, TABLE_ROW_BUFFSIZE, fpIn );

	// write new header
	fprintf(fpOut, "Response,Weights,X0,Y0,X1,Y1\n");

	// copy first 6 fields of the composite table file...
	int nThis = 0;
	while(1)
	{
		++nThis;
		if (0 == nThis % 1000)  // increment progress bar every 1000 rows...
		{
			if (++nCurrent > 100) nCurrent = 1;
			fptr("Creating Composite Table From Composite Table...", nCurrent);
		}
		
		if ( NULL == fgets( pRow, TABLE_ROW_BUFFSIZE, fpIn ) )
			break;

		// copy Response
		char *p = strtok(pRow, ",\n");
		fprintf( fpOut, "%lf,", atof(p));

		// copy Weights
		p = strtok(NULL, ",\n");
		fprintf( fpOut, "%lf,", atof(p));

		// copy X0
		p = strtok(NULL, ",\n");
		fprintf( fpOut, "%lf,", atof(p));

		// copy Y0
		p = strtok(NULL, ",\n");
		fprintf( fpOut, "%lf,", atof(p));

		// copy X1
		p = strtok(NULL, ",\n");
		fprintf( fpOut, "%lf,", atof(p));

		// copy Y1
		p = strtok(NULL, ",\n");
		fprintf( fpOut, "%lf\n", atof(p));
	}

	if (pRow) delete[] pRow;
	if (fpIn) fclose(fpIn);
	if (fpOut) fclose(fpOut);
	fptr("Creating Composite Table From Composite Table...", 0);
	return(true);
}



//
// Create a GDM composite file from an existing GDM composite 
// file but masking with an ESRI Binary Export Grid.
// (Only copies the first 6 columns).
//
bool CreateMaskedCompositeFromComposite(char *InputPath, 
	                                    char *MaskPath, 
										char *OutputPath, 
										bool DoAggregation,
										FPTR fptr)
{
	// sanity check there there are at least six columns
	if (!GotEnoughColsInComposite(InputPath))
	{
		Message("Not enough columns in Input Table in CreateCompositeFromComposite", "ERROR");
		return(false);
	}

	int nCurrent = 0;
	fptr("Creating Composite Table From Composite Table...", nCurrent);

	FILE *fpIn = fopen( InputPath, "r+t");
	if ( NULL == fpIn)
	{
		Message("Cannot open input table in CreateCompositeFromComposite", "ERROR");
		return(false);
	}

	// open output text file
	FILE *fpOut = fopen( OutputPath, "w+t");
	if ( NULL == fpOut)
	{
		Message("Cannot open output table in CreateCompositeFromComposite", "ERROR");
		if (fpIn) fclose(fpIn);
		return(false);
	}

	// get header
	char *pRow = new char [TABLE_ROW_BUFFSIZE];
	fgets( pRow, TABLE_ROW_BUFFSIZE, fpIn );

	// write new header
	fprintf(fpOut, "Response,Weights,X0,Y0,X1,Y1\n");

	gmpath GmPath;	
	char *pBuffer = new char [FILEPATH_BUFFSIZE];
	strcpy(pBuffer, MaskPath);
	EsriBinaryHeader *pHeader = new EsriBinaryHeader(GmPath.ChangeExtension(pBuffer, ".hdr"));
	BinaryFileClass *pBinary = new BinaryFileClass(GmPath.ChangeExtension(pBuffer, ".flt"));

	// copy first 6 fields of the composite table file...
	int nThis = 0;
	while(1)
	{
		++nThis;
		if (0 == nThis % 1000)  // increment progress bar every 1000 rows...
		{
			if (++nCurrent > 100) nCurrent = 1;
			fptr("Creating Masked Composite Table From Composite Table...", nCurrent);
		}
		
		if ( NULL == fgets( pRow, TABLE_ROW_BUFFSIZE, fpIn ) )
			break;

		// do Response
		char *p = strtok(pRow, ",\n");
		double dResponse = atof(p);		

		// do Weights
		p = strtok(NULL, ",\n");
		double dWeight = atof(p);

		// do X0
		p = strtok(NULL, ",\n");
		double dX0 = atof(p);

		// do Y0
		p = strtok(NULL, ",\n");
		double dY0 = atof(p);

		// copy X1
		p = strtok(NULL, ",\n");
		double dX1 = atof(p);

		// copy Y1
		p = strtok(NULL, ",\n");
		double dY1 = atof(p);

		bool fDataIsOK = true;
		if ( (!SiteInBoundingRectangle(dX0, dY0, pHeader)) || (!SiteInBoundingRectangle(dX1, dY1, pHeader)))
		{
			fDataIsOK = false;
			continue;
		}

		if ((!SiteHasData(dX0, dY0, pHeader, pBinary)) || (!SiteHasData(dX1, dY1, pHeader, pBinary)))
		{
			fDataIsOK = false;
			continue;
		}

		if (fDataIsOK)  // data is to be used
		{
			if (DoAggregation)
			{
				// extract the centroids of the current site pair
				AggregateSitePair(&dX0, &dY0, &dX1, &dY1, pHeader );

				// write the record after the sites have been aggregated to the gridcell centroid
				fprintf( fpOut, "%lf,%lf,%lf,%lf,%lf,%lf\n", dResponse, dWeight, dX0, dY0, dX1, dY1);
			}
			else
			{
				// write the record in its current state
				fprintf( fpOut, "%lf,%lf,%lf,%lf,%lf,%lf\n", dResponse, dWeight, dX0, dY0, dX1, dY1);
			}
		}
	} // while(1)

	if (pRow) delete[] pRow;
	if (fpIn) fclose(fpIn);
	if (fpOut) fclose(fpOut);
	if (pBuffer) delete[] pBuffer;
	if (pHeader) delete pHeader;
	if (pBinary) delete pBinary;
	fptr("Creating Composite Table From Composite Table...", 0);
	return(true);
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	WORKING WITH SITE PLUS SPECIES TABLE DATA AS THE INPUT...
//
//
// Create a GDM composite file from a Site Plus Species file.
//
bool CreateCompositeFromSitePlusSpecies( char *InputPath, 
	                                     char *OutputPath, 
										 bool HaveAbundance, 
										 int WeightType,
										 FPTR fptr)
{
	if (!GotEnoughColsInSitePlusSpecies(InputPath, HaveAbundance))
	{
		Message("Not enough columns in Input Table in CreateCompositeFromSitePlusSpecies", "ERROR");
		return(false);
	}

	//
	// An extra site by species file will be produced for user analysis purposes.
	// This file will have the same name as the output table but with _SXS appended to the filename
	//
	gmpath GmPath;
	char *pTmpSxS = new char [FILEPATH_BUFFSIZE];
	sprintf(pTmpSxS, "%s\\%s_SXS.csv", GmPath.GetDirectoryPath(OutputPath), GmPath.GetName(OutputPath));
	
	if (HaveAbundance)
	{		
		if (!SitePlusSpeciesAbundanceToSxS(InputPath, pTmpSxS, fptr))
		{
			Message("Cannot SitePlusSpeciesAbundanceToSxS in CreateCompositeFromSitePlusSpecies", "ERROR");
			if (pTmpSxS) delete[] pTmpSxS;
			return(false);
		}

		if (!CreateCompositeFromSiteBySpeciesAbundance( pTmpSxS, OutputPath, WeightType, fptr))
		{
			Message("Cannot CreateCompositeFromSiteBySpeciesAbundance in CreateCompositeFromSitePlusSpecies", "ERROR");
			if (pTmpSxS) delete[] pTmpSxS;
			return(false);
		}
	}
	else
	{
		if (!SitePlusSpeciesToSxS(InputPath, pTmpSxS, fptr))
		{
			Message("Cannot SitePlusSpeciesToSxS in CreateCompositeFromSitePlusSpecies", "ERROR");
			if (pTmpSxS) delete[] pTmpSxS;
			return(false);
		}

		if (!CreateCompositeFromSiteBySpecies( pTmpSxS, OutputPath, 0, WeightType, fptr))
		{
			Message("Cannot CreateCompositeFromSiteBySpecies in CreateCompositeFromSitePlusSpecies", "ERROR");
			if (pTmpSxS) delete[] pTmpSxS;
			return(false);
		}
	}

	//if (GmPath.FileExists(pTmpSxS)) remove(pTmpSxS); This file now persists after composite table creation
	if (pTmpSxS) delete[] pTmpSxS;
	return(true);
}


//
// Create a GDM composite file from a Site Plus Species file.
//
bool CreateCompositeFromMaskedSitePlusSpecies( char *InputPath, 
	                                           char *MaskPath,
	                                           char *OutputPath, 
											   bool DoAggregate,
										       bool HaveAbundance, 
										       int WeightType,
										       FPTR fptr)
{
	if (!GotEnoughColsInSitePlusSpecies(InputPath, HaveAbundance))
	{
		Message("Not enough columns in Input Table in CreateCompositeFromMaskedSitePlusSpecies", "ERROR");
		return(false);
	}

	gmpath GmPath;

	//
	// We need to create a version of the Site Plus Species Table that is Masked 
	// by the mask grid and optionally aggregated to gridcell centroids
	//
	char *pTmpSPlusS = new char [FILEPATH_BUFFSIZE];
	sprintf(pTmpSPlusS, "%stmpspluss.csv", GmPath.GetDirectoryPath(OutputPath));
	if (!MaskSitePlusSpeciesTable(InputPath, pTmpSPlusS, MaskPath, DoAggregate, HaveAbundance, fptr))
	{
		Message("Unable to MaskSitePlusSpeciesTable in CreateCompositeFromMaskedSitePlusSpecies", "ERROR");
		return(false);
	}

	
	//
	// An extra site by species file will be produced for user analysis purposes.
	// This file will have the same name as the output table but with _SXS appended to the filename
	//
	char *pTmpSxS = new char [FILEPATH_BUFFSIZE];
	sprintf(pTmpSxS, "%s\\%s_SXS.csv", GmPath.GetDirectoryPath(OutputPath), GmPath.GetName(OutputPath));

	if (HaveAbundance)
	{		
		if (!SitePlusSpeciesAbundanceToSxS(pTmpSPlusS, pTmpSxS, fptr))
		{
			Message("Cannot SitePlusSpeciesAbundanceToSxS in CreateCompositeFromMaskedSitePlusSpecies", "ERROR");
			if (pTmpSxS) delete[] pTmpSxS;
			return(false);
		}

		if (!CreateCompositeFromSiteBySpeciesAbundance( pTmpSxS, OutputPath, WeightType, fptr))
		{
			Message("Cannot CreateCompositeFromSiteBySpeciesAbundance in CreateCompositeFromMaskedSitePlusSpecies", "ERROR");
			if (pTmpSxS) delete[] pTmpSxS;
			return(false);
		}
	}
	else
	{
		if (!SitePlusSpeciesToSxS(pTmpSPlusS, pTmpSxS, fptr))
		{
			Message("Cannot SitePlusSpeciesToSxS in CreateCompositeFromMaskedSitePlusSpecies", "ERROR");
			if (pTmpSxS) delete[] pTmpSxS;
			return(false);
		}

		if (!CreateCompositeFromSiteBySpecies( pTmpSxS, OutputPath, 0, WeightType, fptr))
		{
			Message("Cannot CreateCompositeFromSiteBySpecies in CreateCompositeFromMaskedSitePlusSpecies", "ERROR");
			if (pTmpSxS) delete[] pTmpSxS;
			return(false);
		}
	}

	//if (GmPath.FileExists(pTmpSxS)) remove(pTmpSxS); This file now persists after composite table creation
	if (pTmpSxS) delete[] pTmpSxS;
	if (pTmpSPlusS) delete[] pTmpSPlusS;
	return(true);
}


//
// Create a version of the Site Plus Species Table that is Masked 
// by the mask grid and optionally aggregated to gridcell centroids
//
bool MaskSitePlusSpeciesTable(char *InputPath, 
	                          char *OutputPath, 
							  char *MaskPath, 
							  bool DoAggregation, 
							  bool HaveAbundance, 
							  FPTR fptr)
{
	// sanity check there there are at least six columns
	if (!GotEnoughColsInSitePlusSpecies(InputPath, HaveAbundance))
	{
		Message("Not enough columns in Input Table in MaskSitePlusSpeciesTable", "ERROR");
		return(false);
	}

	int nCurrent = 0;
	fptr("Masking Site Plus Species Table...", nCurrent);

	FILE *fpIn = fopen( InputPath, "r+t");
	if ( NULL == fpIn)
	{
		Message("Cannot open input table in MaskSitePlusSpeciesTable", "ERROR");
		return(false);
	}

	// open output text file
	FILE *fpOut = fopen( OutputPath, "w+t");
	if ( NULL == fpOut)
	{
		Message("Cannot open output table in MaskSitePlusSpeciesTable", "ERROR");
		if (fpIn) fclose(fpIn);
		return(false);
	}

	// allocate row buffers
	char *pRow = new char [TABLE_ROW_BUFFSIZE];
	char *pCopy = new char [TABLE_ROW_BUFFSIZE];
	char *pSpecCode = new char [256];
	char *pAbundance = new char [256];
	
	// get header
	fgets( pRow, TABLE_ROW_BUFFSIZE, fpIn );

	// write new header
	fputs( pRow, fpOut);

	gmpath GmPath;	
	char *pBuffer = new char [FILEPATH_BUFFSIZE];
	strcpy(pBuffer, MaskPath);
	EsriBinaryHeader *pHeader = new EsriBinaryHeader(GmPath.ChangeExtension(pBuffer, ".hdr"));
	BinaryFileClass *pBinary = new BinaryFileClass(GmPath.ChangeExtension(pBuffer, ".flt"));

	// Mask the table
	int nThis = 0;
	while(1)
	{
		++nThis;
		if (0 == nThis % 1000)  // increment progress bar every 1000 rows...
		{
			if (++nCurrent > 100) nCurrent = 1;
			fptr("Masking Site Plus Species Table...", nCurrent);
		}
		
		if ( NULL == fgets( pRow, TABLE_ROW_BUFFSIZE, fpIn ) )
			break;


		// get a copy of the row in case we are to re-write it
		memcpy(pCopy, pRow, TABLE_ROW_BUFFSIZE);


		// get X
		char *p = strtok(pRow, ",\n");
		double dX = atof(p);		

		// get Y
		p = strtok(NULL, ",\n");
		double dY = atof(p);

		// test if in grid bounding rectangle
		bool fDataIsOK = true;
		if (!SiteInBoundingRectangle(dX, dY, pHeader))
		{
			fDataIsOK = false;
			continue;
		}

		// test if has data
		if (!SiteHasData(dX, dY, pHeader, pBinary))
		{
			fDataIsOK = false;
			continue;
		}

		if (fDataIsOK)  // data is to be used
		{
			if (DoAggregation)
			{
				// extract the centroids of the current site pair
				AggregateSite(&dX, &dY, pHeader );

				// write the record after the sites have been aggregated to the gridcell centroid
				if (HaveAbundance)
				{
					p = strtok(NULL, ",\n");
					strcpy(pSpecCode, p);

					p = strtok(NULL, ",\n");
					strcpy(pAbundance, p);

					fprintf( fpOut, "%lf,%lf,%s,%s\n", dX, dY, pSpecCode, pAbundance);
				}
				else
				{
					p = strtok(NULL, ",\n");
					strcpy(pSpecCode, p);

					fprintf( fpOut, "%lf,%lf,%s\n", dX, dY, pSpecCode);
				}
			}
			else
			{
				// write the record in its current state
				fputs( pCopy, fpOut);
			}
		}
	} // while(1)

	// cleanup
	if (pRow) delete[] pRow;
	if (pCopy) delete[] pCopy;
	if (pSpecCode) delete[] pSpecCode;
	if (pAbundance) delete[] pAbundance;
	if (fpIn) fclose(fpIn);
	if (fpOut) fclose(fpOut);
	if (pBuffer) delete[] pBuffer;
	if (pHeader) delete pHeader;
	if (pBinary) delete pBinary;
	fptr("Masking Site Plus Species Table...", 0);
	return(true);
}





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	WORKING WITH SITE BY SPECIES TABLE DATA AS THE INPUT...
//
//
// Create a GDM composite file from a Site By Species file.
//
bool CreateCompositeFromSiteBySpecies( char *InputPath, 
	                                   char *OutputPath, 
									   int Cutpoint, 
									   int WeightType,
									   FPTR fptr)
{
	fptr("Creating In-Memory SxS", 50);
	MemorySxS *clsSxS = new MemorySxS(InputPath);
	int **ppData = clsSxS->GetSxS();

	FILE *fpOut = fopen(OutputPath, "w+t");
	if (NULL == fpOut)
	{
		Message("Cannot open fpOut in CreateCompositeFromSiteBySpecies", "ERROR");
		if (clsSxS) delete clsSxS;
		return(false);
	}
	// write header
	fprintf(fpOut, "Response,Weights,X0,Y0,X1,Y1\n");

	double dX0, dY0, dX1, dY1, dResponse, dWeights;

	// create the composite file...
	int nRows = clsSxS->GetRowCount();
	int nCols = clsSxS->GetColCount();

	int nItems = nRows * (nRows-1) / 2;
	int nThisItem = 0;
	int nCurrent = 0;
	int nUsed = 0;
	int nSkipped = 0;
	fptr("Creating Composite GDM Table from Site By Species...", nCurrent);	
	for( int i=0; i<nRows-1; i++ )
	{
		// check if we have more than the richness cutpoint
		if (GetCountForSite(ppData[i], nCols) < Cutpoint)
		{
			++nSkipped;
			continue;
		}

		for ( int j=i+1; j<nRows; j++ )
		{
			if ( ++nThisItem * 100 / nItems > nCurrent )
			{
				nCurrent = nThisItem * 100 / nItems;
				fptr("Creating Composite GDM Table from Site By Species...", nCurrent);
			}

			// check if we have more than the richness cutpoint
			if (GetCountForSite(ppData[j], nCols) < Cutpoint)
			{
				++nSkipped;
				continue;
			}

			dResponse = BrayCurtisDisimilarityPresenceAbsence(ppData[i], ppData[j], nCols);

			// default 
			dWeights = 1.0; // 	if (WeightType == WT_NO_WEIGHTS)
			// or the following...
			if (WeightType == WT_WEIGHT_BY_SUM)
			{
				dWeights = GetWeightsDefault(ppData[i], ppData[j], nCols);
			}
			else if (WeightType == WT_WEIGHT_BY_SQRTSUM)
			{
				dWeights = GetWeightsSqrt(ppData[i], ppData[j], nCols);
			}

			dX0 = clsSxS->GetX()[i];
			dY0 = clsSxS->GetY()[i];
			dX1 = clsSxS->GetX()[j];
			dY1 = clsSxS->GetY()[j];
			fprintf(fpOut, "%lf,%lf,%lf,%lf,%lf,%lf\n", dResponse, dWeights, dX0, dY0, dX1, dY1);
			++nUsed;
		}
	}
	//Message(nUsed, "nUsed");
	//Message(nSkipped, "nSkipped");

	fptr("Creating Composite GDM Table from Site By Species...", 0);
	if (fpOut) fclose(fpOut);
	if (clsSxS) delete clsSxS;
	//Message(nThisItem, "Items Written");
	return(true);
}


//
// Create a GDM composite file from a Site By Species Abundance file.
//
bool CreateCompositeFromSiteBySpeciesAbundance( char *InputPath, char *OutputPath, int WeightType, FPTR fptr)
{
	fptr("Creating In-Memory SxS", 50);
	MemorySxSDouble *clsSxS = new MemorySxSDouble(InputPath);
	double **ppData = clsSxS->GetSxS();

	FILE *fpOut = fopen(OutputPath, "w+t");
	if (NULL == fpOut)
	{
		Message("Cannot open fpOut in CreateCompositeFromSiteBySpeciesAbundance", "ERROR");
		if (clsSxS) delete clsSxS;
		return(false);
	}
	// write header
	fprintf(fpOut, "Response,Weights,X0,Y0,X1,Y1\n");

	double dX0, dY0, dX1, dY1, dResponse, dWeights;

	// create the composite file...
	int nRows = clsSxS->GetRowCount();
	int nCols = clsSxS->GetColCount();
	int nItems = nRows * (nRows-1) / 2;
	int nThisItem = 0;
	int nCurrent = 0;
	fptr("Creating Composite GDM Table from Site By Species...", nCurrent);	
	for( int i=0; i<nRows-1; i++ )
	{
		for ( int j=i+1; j<nRows; j++ )
		{
			if ( ++nThisItem * 100 / nItems > nCurrent )
			{
				nCurrent = nThisItem * 100 / nItems;
				fptr("Creating Composite GDM Table from Site By Species...", nCurrent);
			}

			dResponse = BrayCurtisDissimilarityAbundanceDouble(ppData[i], ppData[j], nCols);

			// default 
			dWeights = 1.0; // 	if (WeightType == WT_NO_WEIGHTS)
			// or the following...
			if (WeightType == WT_WEIGHT_BY_SUM)
			{
				dWeights = GetWeightsDoubleDefault(ppData[i], ppData[j], nCols);
			}
			else if (WeightType == WT_WEIGHT_BY_SQRTSUM)
			{
				dWeights = GetWeightsDoubleSqrt(ppData[i], ppData[j], nCols);
			}
			
			dX0 = clsSxS->GetX()[i];
			dY0 = clsSxS->GetY()[i];
			dX1 = clsSxS->GetX()[j];
			dY1 = clsSxS->GetY()[j];
			fprintf(fpOut, "%lf,%lf,%lf,%lf,%lf,%lf\n", dResponse, dWeights, dX0, dY0, dX1, dY1);
		}
	}
	fptr("Creating Composite GDM Table from Site By Species...", 0);
	if (fpOut) fclose(fpOut);
	if (clsSxS) delete clsSxS;
	//Message(nThisItem, "Items Written");
	return(true);
}



//
// Create a version of the Site By Species Table that is Masked 
// by the mask grid and optionally aggregated to gridcell centroids
//
bool MaskSiteBySpeciesTable(char *InputPath, 
	                        char *OutputPath, 
							char *MaskPath, 
							bool DoAggregation, 
							FPTR fptr)
{
	int nCurrent = 0;
	fptr("Masking Site By Species Table...", nCurrent);

	FILE *fpIn = fopen( InputPath, "r+t");
	if ( NULL == fpIn)
	{
		Message("Cannot open input table in MaskSiteBySpeciesTable", "ERROR");
		return(false);
	}

	// open output text file
	FILE *fpOut = fopen( OutputPath, "w+t");
	if ( NULL == fpOut)
	{
		Message("Cannot open output table in MaskSiteBySpeciesTable", "ERROR");
		if (fpIn) fclose(fpIn);
		return(false);
	}

	// allocate row buffers
	char *pRow = new char [TABLE_ROW_BUFFSIZE];
	char *pCopy = new char [TABLE_ROW_BUFFSIZE];
	
	// get header
	fgets( pRow, TABLE_ROW_BUFFSIZE, fpIn );

	// write new header
	fputs( pRow, fpOut);

	gmpath GmPath;	
	char *pBuffer = new char [FILEPATH_BUFFSIZE];
	strcpy(pBuffer, MaskPath);
	EsriBinaryHeader *pHeader = new EsriBinaryHeader(GmPath.ChangeExtension(pBuffer, ".hdr"));
	BinaryFileClass *pBinary = new BinaryFileClass(GmPath.ChangeExtension(pBuffer, ".flt"));

	// Mask the table
	int nThis = 0;
	while(1)
	{
		++nThis;
		if (0 == nThis % 100)  // increment progress bar every 100 rows...
		{
			if (++nCurrent > 100) nCurrent = 1;
			fptr("Masking Site By Species Table...", nCurrent);
		}
		
		if ( NULL == fgets( pRow, TABLE_ROW_BUFFSIZE, fpIn ) )
			break;


		// get a copy of the row in case we are to re-write it
		memcpy(pCopy, pRow, TABLE_ROW_BUFFSIZE);


		// get X
		char *p = strtok(pRow, ",\n");
		double dX = atof(p);		

		// get Y
		p = strtok(NULL, ",\n");
		double dY = atof(p);


		// test if in grid bounding rectangle
		bool fDataIsOK = true;
		if (!SiteInBoundingRectangle(dX, dY, pHeader))
		{
			fDataIsOK = false;
			continue;
		}


		// test if has data
		if (!SiteHasData(dX, dY, pHeader, pBinary))
		{
			fDataIsOK = false;
			continue;
		}


		// data is to be used
		if (fDataIsOK)  
		{
			if (DoAggregation)
			{
				// extract the centroids of the current site pair
				AggregateSite(&dX, &dY, pHeader );

				char *p0 = strchr(pCopy, ',');  // get the first comma
				++p0;                           // point to the first character after the first comma
				char *p1 = strchr(p0, ',');     // get the second comma
				++p1;                           // point to the first character after the second comma

				// write the agregated X and Y followed by the species data
				fprintf( fpOut, "%lf,%lf,%s", dX, dY, p1);
			}
			else
			{
				// write the record in its current state
				fputs( pCopy, fpOut);
			}
		}
	} // while(1)

	// cleanup
	if (pRow) delete[] pRow;
	if (pCopy) delete[] pCopy;
	if (fpIn) fclose(fpIn);
	if (fpOut) fclose(fpOut);
	if (pBuffer) delete[] pBuffer;
	if (pHeader) delete pHeader;
	if (pBinary) delete pBinary;
	fptr("Masking Site By Species Table...", 0);
	return(true);
}





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	SUPPORT METHODS...
//
//
// Check if site by species table contains values other than 0 or 1 and return true if so
//
bool SxSDataHasFloatingPoint( char *InputPath )
{
	char *pData = new char [TABLE_ROW_BUFFSIZE];
	FILE *fp = fopen(InputPath, "r+t");

	// get header
	fgets(pData, TABLE_ROW_BUFFSIZE, fp);

	while(1)
	{
		if ( NULL == fgets(pData, TABLE_ROW_BUFFSIZE, fp))
			break;

		char *p = strtok(pData, ",\n"); // get X
		p = strtok(NULL, ",\n");        // get Y

		// now search for site by species values not equal to zero or one
		while(1)
		{
			p = strtok(NULL, ",\n");
			if (NULL == p) 
				break;

			double dVal = atof(p);
			if ((dVal != 0.0) && (dVal != 1.0))
			{
				if (fp) fclose(fp);
				if (pData) delete[] pData;
				return(true);
			}
		}
	}

	if (fp) fclose(fp);
	if (pData) delete[] pData;
	return(false);
}



//
// Return true if there are at least six comma delimited fields in the data
//
bool GotEnoughColsInComposite(char *InputPath)
{
	FILE *fpIn = fopen( InputPath, "r+t");
	if ( NULL == fpIn)
	{
		Message("Cannot open input table in GotEnoughColsInComposite", "ERROR");
		return(false);
	}

	// get header and close file
	char *pHeader = new char [TABLE_ROW_BUFFSIZE];
	fgets( pHeader, TABLE_ROW_BUFFSIZE, fpIn );
	if (fpIn) fclose(fpIn);
	
	// get Response
	char *p = strtok(pHeader, ",\n");
	if (!p) 
	{
		if (pHeader) delete[] pHeader;
		return(false);
	}

	// get Weights
	p = strtok(NULL, ",\n");
	if (!p) 
	{
		if (pHeader) delete[] pHeader;
		return(false);
	}

	// get X0
	p = strtok(NULL, ",\n");
	if (!p) 
	{
		if (pHeader) delete[] pHeader;
		return(false);
	}

	// get Y0
	p = strtok(NULL, ",\n");
	if (!p) 
	{
		if (pHeader) delete[] pHeader;
		return(false);
	}

	// get X1
	p = strtok(NULL, ",\n");
	if (!p) 
	{
		if (pHeader) delete[] pHeader;
		return(false);
	}

	// get Y1
	p = strtok(NULL, ",\n");
	if (!p) 
	{
		if (pHeader) delete[] pHeader;
		return(false);
	}

	if (pHeader) delete[] pHeader;
	return(true);
}



//
// Return true if there are three comma delimited fields in the data or four if HaveAbundance
//
bool GotEnoughColsInSitePlusSpecies(char *InputPath, bool HaveAbundance)
{
	FILE *fpIn = fopen( InputPath, "r+t");
	if ( NULL == fpIn)
	{
		Message("Cannot open input table in GotEnoughColsInComposite", "ERROR");
		return(false);
	}

	// get header and close file
	char *pHeader = new char [TABLE_ROW_BUFFSIZE];
	fgets( pHeader, TABLE_ROW_BUFFSIZE, fpIn );
	if (fpIn) fclose(fpIn);

	// get X
	char *p = strtok(pHeader, ",\n");
	if (!p) 
	{
		if (pHeader) delete[] pHeader;
		return(false);
	}

	// get Y
	p = strtok(NULL, ",\n");
	if (!p) 
	{
		if (pHeader) delete[] pHeader;
		return(false);
	}

	// get species code
	p = strtok(NULL, ",\n");
	if (!p) 
	{
		if (pHeader) delete[] pHeader;
		return(false);
	}

	// get abundance field if HaveAbundance is True
	if (HaveAbundance)
	{
		p = strtok(NULL, ",\n");
		if (!p) 
		{
			if (pHeader) delete[] pHeader;
			return(false);
		}
	}
	if (pHeader) delete[] pHeader;
	return(true);
}



//
// Adjust the sites to center them in their respective gridcells
//
void AggregateSitePair(double *pX0, double *pY0, double *pX1, double *pY1, EsriBinaryHeader *pHeader )
{
	// This allows sites that are in the same gridcell but with slightly 
	// differing coordinates to be coerced into the same site. 
	//
	double FullCell = pHeader->GetCellSize();
	double HalfCell = pHeader->GetCellSize() / 2.0;
	double tmpX0 = *pX0;
	double tmpY0 = *pY0;
	double tmpX1 = *pX1;
	double tmpY1 = *pY1;

	// calculate the linear memory block offset from the X and Y values
	int nX0 = int( floor( ( tmpX0 - pHeader->GetMinX() ) / FullCell ) );
	int nY0 = int( floor( ( pHeader->GetMaxY() - tmpY0 ) / FullCell ) );
	int nX1 = int( floor( ( tmpX1 - pHeader->GetMinX() ) / FullCell ) );
	int nY1 = int( floor( ( pHeader->GetMaxY() - tmpY1 ) / FullCell ) );
	
	*pX0 = pHeader->GetMinX() + (nX0 * FullCell) + HalfCell;
	*pY0 = pHeader->GetMaxY() - (nY0 * FullCell) - HalfCell;
	*pX1 = pHeader->GetMinX() + (nX1 * FullCell) + HalfCell;
	*pY1 = pHeader->GetMaxY() - (nY1 * FullCell) - HalfCell;
}


//
// Adjust the site to center in gridcell centroid
//
void AggregateSite(double *pX0, double *pY0, EsriBinaryHeader *pHeader )
{
	// This allows sites that are in the same gridcell but with slightly 
	// differing coordinates to be coerced into the same site. 
	//
	double FullCell = pHeader->GetCellSize();
	double HalfCell = pHeader->GetCellSize() / 2.0;
	double tmpX0 = *pX0;
	double tmpY0 = *pY0;

	// calculate the linear memory block offset from the X and Y values
	int nX0 = int( floor( ( tmpX0 - pHeader->GetMinX() ) / FullCell ) );
	int nY0 = int( floor( ( pHeader->GetMaxY() - tmpY0 ) / FullCell ) );
	
	*pX0 = pHeader->GetMinX() + (nX0 * FullCell) + HalfCell;
	*pY0 = pHeader->GetMaxY() - (nY0 * FullCell) - HalfCell;
}

