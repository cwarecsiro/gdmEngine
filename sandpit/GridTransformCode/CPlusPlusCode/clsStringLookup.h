//
// StringLookup.h
//
#ifndef __CLSSTRINGLOOKUP_H__
#define __CLSSTRINGLOOKUP_H__

#include "stdafx.h"

//#define REC_LENGTH 64
#define REC_LENGTH 128

typedef struct tabSTRING_REC64 {
	char rec[REC_LENGTH];
} STRING_REC64;



class GDM_StringLU
{
	STRING_REC64 *pTable;

	int nSize;

	public:

		GDM_StringLU( int );

		STRING_REC64 * GDM_StringLU::GetData();

		int Count() { return( nSize ); }

		void InsertAt( int, char * );

		STRING_REC64 *FindRecord( char * );

		int GetRecordIndex( char * );

		char *GetStringAt( int n ) { return( pTable[n].rec ); }

		void Sort();

		void WriteToFile( char * );

		void WriteToBinaryFile( char * );

		~GDM_StringLU();


	private:

};


#endif // __CLSSTRINGLOOKUP_H__
