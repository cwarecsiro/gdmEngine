//
// clsMakeSxS.cpp
//
#include "stdafx.h"
#include "clsMakeSxS.h"
#include "Message.h"

//
// Constructor
//
GDM_SiteBySpecies::GDM_SiteBySpecies( int nRows, int nCols )
{
	numRows = nRows;
	numCols = nCols;

	ppSxS = new int * [ nRows ];
	if ( NULL == ppSxS )
	{
		Message( "Cannot allocate ppSxS", ", GDM_SiteBySpecies::GDM_SiteBySpecies( int nRows, int nCols )" );
		return;
	}


	//printf( "GDM_SiteBySpecies()      nRows: %d   nols: %d\n", nRows, nCols );

	for ( int i=0; i<nRows; i++ )
	{
		ppSxS[i] = new int [nCols];
		if ( NULL == ppSxS[i] )
		{
			Message( "Cannot allocate ppSxS[i]", ", GDM_SiteBySpecies::GDM_SiteBySpecies( int nRows, int nCols )" );
			// unwind any previous allocations
			for ( int j=0; j<i; j++ )
			{
				delete[] ppSxS[j];
			}
			delete[] ppSxS;
			return;
		}

		// initialise the row to all zeros
		for ( int j=0; j<nCols; j++ ) ppSxS[i][j] = 0;
	}
}



////////////////////////////////////////////////////////////////////////////////////
// Set the cell to present (1)
//
void 
GDM_SiteBySpecies::SetToPresent( int r, int c )
{
	ppSxS[r][c] = 1;
}


////////////////////////////////////////////////////////////////////////////////////
// Add the abundance value to the cell at index [r,c]
//
bool 
GDM_SiteBySpecies::AddAbundance( int r, int c, int abundance)
{
	if ( abundance < 0 )
	{
		Message("Abundance values MUST be 1..255", "Range Exception");
		return(false);
	}
	else
	{
		ppSxS[r][c] += abundance;
		return(true);
	}
}


////////////////////////////////////////////////////////////////////////////////////
// Write contents of SxS to comma delimited text file without header
//
void
GDM_SiteBySpecies::WriteTableToFileNoHeader( char *pPath )
{
	FILE *fp = fopen( pPath, "w+t" );

	// write the data
	for ( int i=0; i<numRows; i++ )
	{
		// write the one based ID
		fprintf( fp, "%d,", i+1 );

		for ( int j=0; j<numCols; j++ )
		{
			if ( j < (numCols-1) )
			{
				fprintf( fp, "%d,", ppSxS[i][j] );
			}
			else
			{
				fprintf( fp, "%d\n", ppSxS[i][j] );
			}
		}
	}

	fclose( fp );
}




////////////////////////////////////////////////////////////////////////////////////
// Write contents of SxS to comma delimited text file with header
//
void
GDM_SiteBySpecies::WriteTableToFileWithHeader( char *pPath, GDM_StringLU *pSpeciesLU )
{
	FILE *fp = fopen( pPath, "w+t" );

	// write a header
	fprintf( fp, "%s", "_ID" );
	for ( int i=0; i<numCols; i++ )
	{
		fprintf( fp, ",%s", pSpeciesLU->GetStringAt(i) );
	}
	fprintf( fp, "\n" );

	// write the data
	for ( int i=0; i<numRows; i++ )
	{
		// write the one based ID
		fprintf( fp, "%d,", i+1 );

		for ( int j=0; j<numCols; j++ )
		{
			if ( j < (numCols-1) )
			{
				fprintf( fp, "%d,", ppSxS[i][j] );
			}
			else
			{
				fprintf( fp, "%d\n", ppSxS[i][j] );
			}
		}
	}

	fclose( fp );
}




////////////////////////////////////////////////////////////////////////////////////
// Write species presents weights to comma delimited text file 
//
void 
GDM_SiteBySpecies::WriteWeightsToFile( char *pPath )
{
	FILE *fp = fopen( pPath, "w+t" );

	// write a header
	fprintf( fp, "%s,%s\n", "_ID", "Present" );

	for ( int i=0; i<numRows; i++ )
	{
		// count present items for this row
		int nCount = 0;
		for ( int j=0; j<numCols; j++ ) 
		{
			nCount += ppSxS[i][j];
		}

		// write the record
		fprintf( fp, "%d,%d\n", i+1, nCount );
	}

	fclose( fp );
}




////////////////////////////////////////////////////////////////////////////////////
// Write a summary table to comma delimited file detailing 
//		number of unique sites
//		number of unique species
//		number of present records
//		number of absence records
void
GDM_SiteBySpecies::WriteSummaryTablePresAbs( char *pPath )
{
	
	// get a count of the present and absent items
	int nPresent = 0;
	int nAbsent = 0;

	for ( int i=0; i<numRows; i++ )
	{
		for ( int j=0; j<numCols; j++ ) 
		{
			( ppSxS[i][j] == 0 ) ? ++nAbsent : ++nPresent;
		}
	}

	FILE *fp = fopen( pPath, "w+t" );

	// write a header
	fprintf( fp, "_Sites,_Species,_Present,_Absent\n" );

	// write the details
	fprintf( fp, "%d,%d,%d,%d\n", numRows, numCols, nPresent, nAbsent );

	fclose( fp );
}


////////////////////////////////////////////////////////////////////////////////////
// Write a summary table to comma delimited file detailing 
//		number of unique sites
//		number of unique species
//		summed abundance for each unique site
void
GDM_SiteBySpecies::WriteSummaryTableAbundance( char *pPath )
{
	FILE *fp = fopen( pPath, "w+t" );

	// write a header
	fprintf( fp, "_Sites,_Species\n" );

	// write the details
	fprintf( fp, "%d,%d\n", numRows, numCols );

	// write the abundance header
	fprintf( fp, "Site_ID,Abundance\n" );
	
	// write the abundance sums
	for ( int i=0; i<numRows; i++ )
	{
		int nSum = 0;
		for ( int j=0; j<numCols; j++ ) 
		{
			nSum += ppSxS[i][j];				
		}

		// write the details
		fprintf( fp, "%d,%d\n", i+1, nSum );
	}
	
	fclose( fp );
}



////////////////////////////////////////////////////////////////////////////////////
// Destructor
//
GDM_SiteBySpecies::~GDM_SiteBySpecies()
{
	for ( int i=0; i<numRows; i++ )
	{
		delete[] ppSxS[i];
	}
	delete[] ppSxS;
}

