//
// clsMakeSxS.h
//
#ifndef __CLSMAKESXS_H__
#define __CLSMAKESXS_H__

#include "clsStringLookup.h"

class GDM_SiteBySpecies
{
	int **ppSxS;
	int numRows;
	int numCols;

	public:

		GDM_SiteBySpecies( int, int );

		void WriteTableToFileNoHeader( char * );

		void WriteTableToFileWithHeader( char *, GDM_StringLU * );

		void WriteWeightsToFile( char * );

		void WriteSummaryTablePresAbs( char * );

		void WriteSummaryTableAbundance( char * );

		void SetToPresent( int, int );

		bool AddAbundance( int, int, int );

		int **GetSxSData() { return(ppSxS); }

		~GDM_SiteBySpecies();
};


#endif // __CLSMAKESXS_H__