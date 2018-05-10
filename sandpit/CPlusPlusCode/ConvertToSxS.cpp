//
//
// ConvertToSxS.cpp
//
//
#include "stdafx.h"

#include "clsHashTab.h"
#include "clsSpecRec.h"
#include "clsAbundanceRec.h"
#include "clsStringLookup.h"
#include "clsMakeSxS.h"

#include "ConvertToSxS.h"
#include "GDMBufferSizeDefs.h"
#include "Message.h"
#include "ParamsW16.h"


///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Convert a Comma delimited textfile in the form X,Y,Species_Code
// to a Site By Species File in the form X,Y,Sp1,Sp2,...,SpN where
// Sp1..SpN are presence/absence values (1/0)
//
bool SitePlusSpeciesToSxS(char *pInput, char *pOutput, FPTR fptr)
{
	char buff[SIZE_KEY];
	int nIndex;
	fptr("Convert S+S data to Site By Species.", 0);

	///////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Create a SpecRec class containing all the X,Y,Codes from the input file p
	//
	GDM_SpecRec *pSpec = new GDM_SpecRec( pInput );
	int nRecs = pSpec->GetRecordCount();
	//
	//
	// Create a hash table for the concatenated X and Y values
	//
	Hashtable *hTab_Sites = new Hashtable(nRecs);
	//
	// populate the hash table with concatenated X and Y values comma delimited as strings
	//
	int nCurrent = 0;
	for ( int i=0; i<nRecs; i++ )
	{
		if ( i * 100 / nRecs > nCurrent )
		{
			nCurrent = i * 100 / nRecs;
			fptr("Convert S+S data to Site By Species.", nCurrent);
		}

		sprintf( buff, "%1.10lf,%1.10lf", pSpec->GetXAt(i),  pSpec->GetYAt(i) );
		NODE *node = new NODE(buff, buff);
		if (!hTab_Sites->contains(buff))
		{
			hTab_Sites->put( node );
		}
	}
	//
	// get the Unique sites count
	//
	int nUniqueSites = hTab_Sites->getSize();
	//
	// allocate a StringLU class for the location data
	//
	GDM_StringLU *pSiteLU = new GDM_StringLU( nUniqueSites );
	//
	// Populate the table with Unique Site records
	//
	nIndex = 0;
	hTab_Sites->initIterator(); 
	while(hTab_Sites->hasNext()) 
	{ 
		char tmp[SIZE_KEY];
		hTab_Sites->getNextKey(tmp); 
		pSiteLU->InsertAt( nIndex++, tmp );
	} 
	//
	// the records need sorting for the binary search routines to work
	//
	pSiteLU->Sort();


	///////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Create a hash table for the species codes
	//
	Hashtable *hTab_Codes = new Hashtable(nRecs);
	//
	// populate the hash table with concatenated X and Y values comma delimited as strings
	//
	nCurrent = 0;
	for ( int i=0; i<nRecs; i++ )
	{
		if ( i * 100 / nRecs > nCurrent )
		{
			nCurrent = i * 100 / nRecs;
			fptr("Convert S+S data to Site By Species.", nCurrent);
		}

		sprintf( buff, "%s", pSpec->GetCodeAt(i) );
		NODE *node = new NODE(buff, buff);
		if (!hTab_Codes->contains(buff))
		{
			hTab_Codes->put( node );
		}
	}
	//
	// get the Unique sites count
	//
	int nUniqueCodes = hTab_Codes->getSize();
	//
	// allocate a StringLU class for the location data
	//
	GDM_StringLU *pCodeLU = new GDM_StringLU( nUniqueCodes );
	//
	// Populate the table with Unique Site records
	//
	nIndex = 0;
	hTab_Codes->initIterator(); 
	while(hTab_Codes->hasNext()) 
	{ 
		char tmp[SIZE_KEY];
		hTab_Codes->getNextKey(tmp); 
		pCodeLU->InsertAt( nIndex++, tmp );
	} 
	//
	// the records need sorting for the binary search routines to work
	//
	pCodeLU->Sort();


	////////////////////////////////////////////////////////////////////////////////////////////////	
	//
	// create the Site by Species table
	//
	GDM_SiteBySpecies *pSXS = new GDM_SiteBySpecies( nUniqueSites, nUniqueCodes );	
	int nSSS = 0;
	nCurrent = 0;
	for ( int i=0; i<nRecs; i++ )
	{
		if ( i * 100 / nRecs > nCurrent )
		{
			nCurrent = i * 100 / nRecs;
			fptr("Convert S+S data to Site By Species.", nCurrent);
		}

		// get the location record index 
		sprintf( buff, "%1.10lf,%1.10lf", pSpec->GetXAt( i ),  pSpec->GetYAt( i ) );
		int nLocIndex = pSiteLU->GetRecordIndex( buff );

		// get the species record index 
		int nSpecIndex = pCodeLU->GetRecordIndex( pSpec->GetCodeAt( i ) );

		// set the SxS data at this index to present (1)
		pSXS->SetToPresent( nLocIndex, nSpecIndex );
	}


	////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// rewrite the filtered table as a Site By Species Table
	//
	//FuncPtrGlobal("Step 2: Convert input data to Site By Species.", 50, ProgBarGlobal);
	FILE *fp = fopen(pOutput, "w+t");
	//
	// write header
	//
	fprintf( fp, "X,Y," );
	for ( int i=0; i<nUniqueCodes; i++ )
	{
		fprintf( fp, pCodeLU->GetStringAt(i) );
		if ( i < nUniqueCodes-1 )
		{
			fprintf( fp, "," );
		}
		else
		{
			fprintf( fp, "\n" );	
		}
	}
	
	// write site by species data
	nCurrent = 0;
	for ( int i=0; i<pSiteLU->Count(); i++ )
	{
		if ( i * 100 / pSiteLU->Count() > nCurrent )
		{
			nCurrent = i * 100 / pSiteLU->Count();
			fptr("Convert S+S data to Site By Species.", nCurrent);
		}

		// write the site data
		fprintf( fp, "%s,", pSiteLU->GetData()[i].rec );

		// write the species data
		for ( int j=0; j<nUniqueCodes; j++ )
		{
			// write the site by species record
			fprintf( fp, "%d", pSXS->GetSxSData()[i][j] );

			if ( j < nUniqueCodes-1 )
				fprintf( fp, "," );
			else
				fprintf( fp, "\n" );
		}
	}
	fclose(fp);
	fptr("Convert S+S data to Site By Species.", 0);

	//
	// Cleanup
	//
	if ( pSpec ) delete pSpec;
	if ( hTab_Sites ) delete hTab_Sites;
	if ( pSiteLU ) delete pSiteLU;
	if ( pCodeLU ) delete pCodeLU;
	if ( pSXS ) delete pSXS;
	return(true);
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Convert a Comma delimited textfile in the form X,Y,Species_Code
// to a Site By Species File in the form X,Y,Sp1,Sp2,...,SpN where
// Sp1..SpN are abundance values (0..255)
//
bool SitePlusSpeciesAbundanceToSxS(char *pInput, char *pOutput, FPTR fptr)
{
	char buff[SIZE_KEY];
	int nIndex;
	fptr("Convert S+S abundance data to Site By Species.", 0);

	///////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Create a SpecRec class containing all the X,Y,Codes from the input file p
	//
	GDM_AbundanceRec *pSpec = new GDM_AbundanceRec( pInput );
	int nRecs = pSpec->GetRecordCount();
	//
	//
	// Create a hash table for the concatenated X and Y values
	//
	Hashtable *hTab_Sites = new Hashtable(nRecs);
	//
	// populate the hash table with concatenated X and Y values comma delimited as strings
	//
	for ( int i=0; i<nRecs; i++ )
	{
		sprintf( buff, "%1.10lf,%1.10lf", pSpec->GetXAt(i),  pSpec->GetYAt(i) );
		NODE *node = new NODE(buff, buff);
		if (!hTab_Sites->contains(buff))
		{
			hTab_Sites->put( node );
		}
	}
	//
	// get the Unique sites count
	//
	int nUniqueSites = hTab_Sites->getSize();
	//
	// allocate a StringLU class for the location data
	//
	GDM_StringLU *pSiteLU = new GDM_StringLU( nUniqueSites );
	//
	// Populate the table with Unique Site records
	//
	nIndex = 0;
	hTab_Sites->initIterator(); 
	while(hTab_Sites->hasNext()) 
	{ 
		char tmp[SIZE_KEY];
		hTab_Sites->getNextKey(tmp); 
		pSiteLU->InsertAt( nIndex++, tmp );
	} 
	//
	// the records need sorting for the binary search routines to work
	//
	pSiteLU->Sort();


	///////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Create a hash table for the species codes
	//
	Hashtable *hTab_Codes = new Hashtable(nRecs);
	//
	// populate the hash table with concatenated X and Y values comma delimited as strings
	//
	for ( int i=0; i<nRecs; i++ )
	{
		sprintf( buff, "%s", pSpec->GetCodeAt(i) );
		NODE *node = new NODE(buff, buff);
		if (!hTab_Codes->contains(buff))
		{
			hTab_Codes->put( node );
		}
	}
	//
	// get the Unique sites count
	//
	int nUniqueCodes = hTab_Codes->getSize();
	//
	// allocate a StringLU class for the location data
	//
	GDM_StringLU *pCodeLU = new GDM_StringLU( nUniqueCodes );
	//
	// Populate the table with Unique Site records
	//
	nIndex = 0;
	hTab_Codes->initIterator(); 
	while(hTab_Codes->hasNext()) 
	{ 
		char tmp[SIZE_KEY];
		hTab_Codes->getNextKey(tmp); 
		pCodeLU->InsertAt( nIndex++, tmp );
	} 
	//
	// the records need sorting for the binary search routines to work
	//
	pCodeLU->Sort();


	////////////////////////////////////////////////////////////////////////////////////////////////	
	//
	// create the Site by Species table
	//
	GDM_SiteBySpecies *pSXS = new GDM_SiteBySpecies( nUniqueSites, nUniqueCodes );	
	int nSSS = 0;
	for ( int i=0; i<nRecs; i++ )
	{
		// get the location record index 
		sprintf( buff, "%1.10lf,%1.10lf", pSpec->GetXAt( i ),  pSpec->GetYAt( i ) );
		int nLocIndex = pSiteLU->GetRecordIndex( buff );

		// get the species record index 
		int nSpecIndex = pCodeLU->GetRecordIndex( pSpec->GetCodeAt( i ) );

		// get the species abundance value
		int nAbundance = pSpec->GetAbundanceAt( i );

		// add to the respective SxS data at this index
		pSXS->AddAbundance( nLocIndex, nSpecIndex, nAbundance );
	}


	////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// rewrite the filtered table as a Site By Species Table
	//
	//FuncPtrGlobal("Step 2: Convert input data to Site By Species.", 50, ProgBarGlobal);
	FILE *fp = fopen(pOutput, "w+t");
	//
	// write header
	//
	fprintf( fp, "X,Y," );
	for ( int i=0; i<nUniqueCodes; i++ )
	{
		fprintf( fp, pCodeLU->GetStringAt(i) );
		if ( i < nUniqueCodes-1 )
		{
			fprintf( fp, "," );
		}
		else
		{
			fprintf( fp, "\n" );	
		}
	}
	
	// write site by species data
	int nCurrent = 0;
	for ( int i=0; i<pSiteLU->Count(); i++ )
	{
		if ( i * 100 / pSiteLU->Count() > nCurrent )
		{
			nCurrent = i * 100 / pSiteLU->Count();
			fptr("Convert S+S abundance data to Site By Species.", nCurrent);
		}

		// write the site data
		fprintf( fp, "%s,", pSiteLU->GetData()[i].rec );

		// write the species data
		for ( int j=0; j<nUniqueCodes; j++ )
		{
			// write the site by species record
			fprintf( fp, "%d", pSXS->GetSxSData()[i][j] );

			if ( j < nUniqueCodes-1 )
				fprintf( fp, "," );
			else
				fprintf( fp, "\n" );
		}
	}
	fclose(fp);
	fptr("Convert S+S abundance data to Site By Species.", 0);

	//
	// Cleanup
	//
	if ( pSpec ) delete pSpec;
	if ( hTab_Sites ) delete hTab_Sites;
	if ( pSiteLU ) delete pSiteLU;
	if ( pCodeLU ) delete pCodeLU;
	if ( pSXS ) delete pSXS;
	return(true);
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Get the Minimum and Maximum counts from a Site By Species Table
//
bool GetSxSMinMax(char *pInput, int *pMin, int *pMax, FPTR fptr)
{
	char *buff = new char [TABLE_ROW_BUFFSIZE];
	if ( NULL == buff )
	{
		Message("Cannot allocate buff", "GetSxSMinMax");
		if (buff) delete[] buff;
		return(false);
	}

	FILE *fp = fopen(pInput, "r+t");
	if ( NULL == fp )
	{
		Message("Cannot open pInput", "GetSxSMinMax");
		if (buff) delete[] buff;
		return(false);
	}

	// get header
	fgets( buff, TABLE_ROW_BUFFSIZE, fp );

	// count number of rows
	int nRows = 0;
	while(1)
	{
		if ( NULL == fgets( buff, TABLE_ROW_BUFFSIZE, fp ) )
		{
			break;
		}
		++nRows;
	}
	rewind(fp);
	
	// get header again
	fgets( buff, TABLE_ROW_BUFFSIZE, fp );

	// count presence or abundance in each row
	int nMin = 1000000;
	int nMax = 0;
	int nCurrent = 0;
	fptr("Extracting Min and Max from Site By Species Table...", nCurrent);
	for ( int i=0; i<nRows; i++ )
	{
		if ( i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			fptr("Extracting Min and Max from Site By Species Table...", nCurrent);
		}

		// get the current row
		fgets( buff, TABLE_ROW_BUFFSIZE, fp );

		// skip the X and Y fields
		char *p = strtok(buff,",");
		p = strtok(NULL,",");

		// Get the sum
		int nSum = 0;
		while (1)
		{
			p = strtok(NULL,",\n");
			if ( NULL == p ) 
				break;
			else
				nSum += atoi(p);
		}

		// update the global max/min vals
		if (nSum < nMin) nMin = nSum;
		if (nSum > nMax) nMax = nSum;
	}
	fptr("Extracting Min and Max from Site By Species Table...", 0);
	*pMin = nMin;
	*pMax = nMax;
	if (fp) fclose(fp);
	if (buff) delete[] buff;
	return(true);
}

