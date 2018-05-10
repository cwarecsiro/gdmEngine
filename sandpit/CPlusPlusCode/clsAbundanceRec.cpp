//
// clsAbundanceRec.cpp
//
#include "stdafx.h"
#include "clsAbundanceRec.h"
#include "GDMBufferSizeDefs.h"
#include "Message.h"


GDM_AbundanceRec::GDM_AbundanceRec( char *pPath )
{
	nRecs = nGetNumRecords( pPath );
	if ( nRecs < 1 )
	{
		Message( "Got No Records", "GDM_SpecRec::GDM_SpecRec( char *pPath )" );
		return;
	}

	// allocate the memory
	pSpec = new ABUNDANCE_REC [ nRecs ];
	if ( NULL == pSpec )
	{
		Message( "GDM_AbundanceRec::GDM_AbundanceRec( char *pPath )", "ERROR" );
		return;
	}


	char *buff = new char [TABLE_ROW_BUFFSIZE];	
	FILE *fp = fopen( pPath, "r+t" );
	fgets( buff, TABLE_ROW_BUFFSIZE, fp );

	// presume X,Y,Code,Abundance format
	char seps[] = ",\n";
	for ( int i=0; i<nRecs; i++ )
	{
		fgets( buff, TABLE_ROW_BUFFSIZE, fp );
		pSpec[i].nID = i+1;

		char *p = strtok( buff, seps ); 
		pSpec[i].dX = atof(p);

		p = strtok( NULL, seps );
		pSpec[i].dY = atof(p);

		p = strtok( NULL, seps );
		strcpy( pSpec[i].pSpec, p );

		p = strtok( NULL, seps );
		pSpec[i].Abundance = atoi(p);
	}
	fclose( fp );
	if ( buff) delete[] buff;
}



int 
GDM_AbundanceRec::nGetNumRecords( char *p )
{
	char *buff = new char [TABLE_ROW_BUFFSIZE];	
	int n = 0;

	FILE *fp = fopen( p, "r+t" );
	fgets( buff, TABLE_ROW_BUFFSIZE, fp );

	while( 1 )
	{
		if ( NULL == fgets( buff, TABLE_ROW_BUFFSIZE, fp ))
		{
			break;
		}

		++n;
	}

	fclose( fp );
	if ( buff) delete[] buff;

	return( n );
}




void 
GDM_AbundanceRec::DumpTableToFile( char *pPath )
{
	FILE *fp = fopen( pPath, "w+t" );
	if ( NULL == fp )
	{
		Message( "Cannot open pPath", "INFO" );
		return;
	}
	
	fprintf( fp, "_ID,X,Y,Spec\n" );

	for ( int i=0; i<nRecs; i++ )
	{
		fprintf( fp, "%d,%1.10lf,%1.10lf,%s\n", pSpec[i].nID, pSpec[i].dX,  pSpec[i].dY, pSpec[i].pSpec );
	}

	fclose( fp );
}




GDM_AbundanceRec::~GDM_AbundanceRec()
{
	delete[] pSpec;
}

