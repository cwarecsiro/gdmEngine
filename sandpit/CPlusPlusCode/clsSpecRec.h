//
// clsSpecRec.h
//
#ifndef __CLSSPECREC_H__
#define __CLSSPECREC_H__


typedef struct tagSPEC_REC {

	int nID;
	double dX;
	double dY;
	char pSpec[64];

} SPEC_REC;


class GDM_SpecRec
{
	SPEC_REC *pSpec;
	int nRecs;


	public:

		GDM_SpecRec( char * );

		int GetRecordCount() { return( nRecs ); }

		SPEC_REC *GetSpecRecs() { return( pSpec ); }

		double GetXAt( int n ) { return( pSpec[n].dX ); }

		double GetYAt( int n ) { return( pSpec[n].dY ); }

		char *GetCodeAt( int n ) { return( pSpec[n].pSpec ); }

		void DumpTableToFile( char * );

		~GDM_SpecRec();


	private:

		int nGetNumRecords( char * );

};


#endif // __CLSSPECREC_H__