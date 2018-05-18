//
// StringLookup.cpp
//
#include "stdafx.h"
#include "clsStringLookup.h"
#include "Message.h"

int Compare( const void *, const void * );
////////////////////////////////////////////////////////////////////////////////////
// String comparison routine for Sort method
//
int Compare( const void *arg1, const void *arg2 )
{
	STRING_REC64 *p1 = (STRING_REC64 *)arg1;
	STRING_REC64 *p2 = (STRING_REC64 *)arg2;
	return( _stricmp( p1->rec, p2->rec ) );
}
////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////
// Constructor
//
GDM_StringLU::GDM_StringLU( int nRecs )
{
	nSize = nRecs;

	pTable = new STRING_REC64 [nRecs];

	if ( NULL == pTable )
	{
		Message( "Unable to allocate pTable", "GDM_StringLU::GDM_StringLU( int nRecs )" );
	}

	for ( int i=0; i<nRecs; i++ )
	{
		memset( pTable[i].rec, 0, REC_LENGTH );
	}
}


////////////////////////////////////////////////////////////////////////////////////
// Return the data
//
STRING_REC64 *
GDM_StringLU::GetData()
{
	return(pTable);
}



////////////////////////////////////////////////////////////////////////////////////
// Insert a string at a position
//
void
GDM_StringLU::InsertAt( int nPos, char *p )
{
	if ( nPos > (nSize-1) )
	{
		Message( "Insert index out of bounds", "GDM_StringLU::InsertAt( int nPos, char *p )" );
	}

	strcpy( pTable[nPos].rec, p );
}



////////////////////////////////////////////////////////////////////////////////////
// Sort the lookup table in ascending order
//
void
GDM_StringLU::Sort()
{
	qsort( (void *)pTable, (size_t)nSize, sizeof(STRING_REC64), Compare );
}



////////////////////////////////////////////////////////////////////////////////////
// Search for and return a record from the lookup table
//
STRING_REC64 *
GDM_StringLU::FindRecord( char *p )
{
	return( STRING_REC64 * )bsearch( p, pTable, (size_t)nSize, sizeof( STRING_REC64 ), Compare );
}



////////////////////////////////////////////////////////////////////////////////////
// Search for and return the zero based index for a record from the lookup table
//
int 
GDM_StringLU::GetRecordIndex( char *p )
{
	STRING_REC64 *pTest = FindRecord( p );
	if ( NULL == pTest )
	{
		// not found
		Message( "Record Not Found", "GDM_StringLU::GetRecordIndex( char *p )" );
		return( -1 );
	}
	else
	{
		//return( int((int)pTest - (int)pTable) / sizeof(STRING_REC64));
		return(int((pTest-pTable) / sizeof(STRING_REC64)));
	}
}





////////////////////////////////////////////////////////////////////////////////////
// Write the lookup table contents into a text file
//
void
GDM_StringLU::WriteToFile( char *p )
{
	FILE *fp = fopen( p, "w+t" );
	if ( NULL == fp )
	{
		Message( "Cannot create text file", "GDM_StringLU::DumpToFile( char *p )" );
	}

	// write the header
	fprintf( fp, "_ID,X,Y\n" );

	for ( int i=0; i<nSize; i++ )
	{
		fprintf( fp, "%d,%s\n", i+1, pTable[i].rec );
	}

	fclose( fp );
}



////////////////////////////////////////////////////////////////////////////////////
// Write the lookup table contents into a binary file
//
void
GDM_StringLU::WriteToBinaryFile( char *Path )
{
	int h = _open( Path, _O_BINARY | _O_CREAT | _O_TRUNC | _O_RDWR, S_IREAD | S_IWRITE );
	if ( h < 0 )
	{
		Message( "Cannot open LOC Binary File", "ERROR" );
		return;
	}

	for ( int i=0; i<nSize; i++ )
	{
		char *p = strtok(pTable[i].rec, "," );
		float fX = (float)atof(p);

		p = strtok(NULL, "," );
		float fY = (float)atof(p);

		_write( h, &fX, sizeof(fX) );
		_write( h, &fY, sizeof(fY) );
	}

	_close( h );
}



////////////////////////////////////////////////////////////////////////////////////
// Destructor
//
GDM_StringLU::~GDM_StringLU()
{
	delete[] pTable;
}

