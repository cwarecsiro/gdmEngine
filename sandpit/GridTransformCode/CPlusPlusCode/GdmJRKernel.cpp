//
// GdmJRKernel.cpp
//
#include "stdafx.h"
#include "GdmBinLib.h"
#include "GdmJRKernel.h"
#include "GDMBufferSizeDefs.h"
#include "clsDoPath.h"
#include "ParamsW16.h"
#include "Message.h"
#include "clsEsriBinaryHeader.h"
#include "clsBinaryFileClass.h"

bool CreateJRTrainingFile(char *pParams, char *JRInputXY, char *JROutputPath, FPTR fptr)
{
	//Message(JRInputXY , "JRInputXY");
	//Message(JROutputPath , "JROutputPath");
	gmpath gmPath;
	//
	char *lpDomainPath = GetPredictorDomainPath(pParams);
	if (!gmPath.FileExists(lpDomainPath))
	{
		// try using the first environmental grid in the predictor list
		GetProfileString( "PREDICTORS", "EnvGrid1", lpDomainPath, pParams );

		// is this grid there?
		if (!gmPath.FileExists(gmPath.ChangeExtension( lpDomainPath, ".flt" )))
		{
			return(false);
		}
	}

	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(lpDomainPath, ".hdr"));
	//
	int nGridRows = header->GetNumRows();
	int nGridCols = header->GetNumCols();
	double ResX = header->GetCellSize();
	float fNoData = header->GetNoDataValue();
	double dResX = header->GetCellSize();
	double dMinX = header->GetXllCorner();
	double dMinY = header->GetYllCorner();
	double dMaxX = dMinX + (dResX * nGridCols);
	double dMaxY = dMinY + (dResX * nGridRows);
	if (header) delete header;

	//
	// are we using euclidean ?
	//
	bool fUseEuclidean = ( 1 == GetProfileInt( "GDMODEL", "UseEuclidean", pParams ) ) ? true : false;


	//
	// determine the number of columns ( the number of transform grids )
	//
	int numTranGrids = 0;
	if (fUseEuclidean) numTranGrids += 2;

	// total number of possible preds, need to check for presence to get file transform count
	int nPreds = GetProfileInt( "PREDICTORS", "NumPredictors", pParams );
	char lpKey[64];
	char lpPredPath[BUFFLEN];
	for ( int i=1; i<=nPreds; i++ )
	{
		sprintf( lpKey, "PredTran%d", i );
		GetProfileString( "TRANSPREDS", lpKey, lpPredPath, pParams );

		if ( strlen( lpPredPath ) > 0 )
		{
			++numTranGrids;
		}
	}


	//
	// allocate strings for transformed input files
	//
	char **ppTranGridPaths = new char * [numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		ppTranGridPaths[i] = new char [BUFFLEN];
	}


	//
	// Collect the Tran grid paths
	//
	int nThis = 0;
	char tmpPath[BUFFLEN];

	// check geographic path
	if (fUseEuclidean)
	{
		GetProfileString( "TRANSPREDS", "EuclXTran", lpPredPath, pParams );
		strcpy(tmpPath, gmPath.ChangeExtension(lpPredPath, ".flt"));

		if (!gmPath.FileExists(tmpPath))
		{
			Message("Cannot find EuclXTran Grid", "ERROR");
			return(false);
		}
		else
		{
			strcpy (ppTranGridPaths[nThis++], tmpPath);
		}

		GetProfileString( "TRANSPREDS", "EuclYTran", lpPredPath, pParams );
		strcpy(tmpPath, gmPath.ChangeExtension(lpPredPath, ".flt"));

		if (!gmPath.FileExists(tmpPath))
		{
			Message("Cannot find EuclYTran Grid", "ERROR");
			return(false);
		}
		else
		{
			strcpy (ppTranGridPaths[nThis++], tmpPath);
		}

	}


	// now collect the environmental grid transform paths
	for (int i=1; i<=nPreds; i++)
	{
		sprintf(lpKey, "PredTran%d", i);
		GetProfileString( "TRANSPREDS", lpKey, lpPredPath, pParams );
		strcpy(tmpPath, gmPath.ChangeExtension(lpPredPath, ".flt"));

		if (gmPath.FileExists(tmpPath))
		{
			strcpy (ppTranGridPaths[nThis++], tmpPath);
		}
	}


	//
	// now create the table
	//
	FILE *fpOut = fopen( JROutputPath, "w+t" );

	// write the header
	fprintf( fpOut, "Class,X,Y" );
	for ( int i=0; i<numTranGrids; i++ )
	{
		char Name[64];
		_splitpath( ppTranGridPaths[i], NULL, NULL, Name, NULL );
		fprintf( fpOut, ",%s", Name );
	}
	fprintf( fpOut, "\n" );


	//
	// determine the number of sites
	//
	FILE *fpIn = fopen( JRInputXY, "r+t");
	char myRow[BUFFLEN];
	fgets( myRow, BUFFLEN, fpIn ); // get header
	// count rows
	int nRows = 0;
	while(1)
	{
		if ( NULL == fgets( myRow, BUFFLEN, fpIn ) )
			break;

		++nRows;
	}
	

	//
	// Create the Binary File Classes
	//
	BinaryFileClass *bfcDomain = new BinaryFileClass(gmPath.ChangeExtension( lpDomainPath, ".flt" ));
	if (!bfcDomain->IsValid())
	{
		Message("Cannot create domain file", "CreateJRTrainingFile");
		return(false);
	}
	BinaryFileClass **bfcTranGrids = new BinaryFileClass * [numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		bfcTranGrids[i] = new BinaryFileClass(ppTranGridPaths[i]);
		if (!bfcTranGrids[i]->IsValid())
		{
			Message("Cannot create transform grid file", "CreateJRTrainingFile");
			return(false);
		}
	}


	//
	// read the input data and for each valid location, extract the transform data
	//
	rewind(fpIn);                  // get back to start
	fgets( myRow, BUFFLEN, fpIn ); // get header
	char seps[] = ",\n";
	int nCurrent = 0;
	nThis = 0;
	while(1)
	{
		if ( NULL == fgets( myRow, BUFFLEN, fpIn ) )
			break;

		if ( ++nThis * 100 / nRows > nCurrent )
		{
			nCurrent = nThis * 100 / nRows;
			if ( nCurrent <= 100 )
				fptr("Filtering training sites", nCurrent);	
		}

		char *p = strtok( myRow, seps ); // get Class ID
		int nClass = atoi(p);
		
		p = strtok( NULL, seps );        // get easting
		double dX = atof(p);

		p = strtok( NULL, seps );        // get northing
		double dY = atof(p);

		//
		// check that the location has data in the GDM domain
		//
		if ( ( dX > dMinX ) && ( dX < dMaxX ) && ( dY > dMinY ) && ( dY < dMaxY ) )
		{
			// calculate the linear memory block offset from the X and Y values
			int nX = int( floor( ( dX - dMinX ) / ResX ) );
			int nY = int( floor( ( dMaxY - dY ) / ResX ) );

			// calculate the seek offset for this site
			int SeekOffset = ((nY * nGridCols) + nX) * sizeof(float);

			float fVal;
			bfcDomain->SeekTo(SeekOffset);
			bfcDomain->ReadFloat(&fVal, 1);

			// check for valid data in the domain grid
			if ( fVal != fNoData )
			{
				// write the class,X,Y fields
				fprintf(fpOut, "%d,%lf,%lf", nClass, dX, dY);

				// extract the transform grid data and write to table
				for ( int i=0; i<numTranGrids; i++ )
				{
					bfcTranGrids[i]->SeekTo(SeekOffset);
					bfcTranGrids[i]->ReadFloat(&fVal, 1);
					fprintf(fpOut, ",%f", fVal);
				}
				fprintf(fpOut, "\n");
			}
		}
	}
	fclose(fpIn);
	fclose(fpOut);
	
	//
	// Clean up
	//
	if (bfcDomain)
	{
		bfcDomain->Close();
		delete bfcDomain;
	}
	for ( int i=0; i<numTranGrids; i++ )
	{
		bfcTranGrids[i]->Close();
		if (bfcTranGrids[i]) delete bfcTranGrids[i];
		if (ppTranGridPaths[i]) delete[] ppTranGridPaths[i];
	}
	if (bfcTranGrids) delete[] bfcTranGrids;
	if (ppTranGridPaths) delete[] ppTranGridPaths;
	return(true);
}


