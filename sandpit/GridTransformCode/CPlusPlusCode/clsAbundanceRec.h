//
// clsAbundanceRec.h
//
#ifndef __CLSABUNDANCEREC_H__
#define __CLSABUNDANCEREC_H__


typedef struct tagABUNDANCE_REC {

	int nID;
	double dX;
	double dY;
	char pSpec[64];
	int Abundance;

} ABUNDANCE_REC;


class GDM_AbundanceRec
{
	ABUNDANCE_REC *pSpec;
	int nRecs;


	public:

		GDM_AbundanceRec( char * );

		int GetRecordCount() { return( nRecs ); }

		ABUNDANCE_REC *GetSpecRecs() { return( pSpec ); }

		double GetXAt( int n ) { return( pSpec[n].dX ); }

		double GetYAt( int n ) { return( pSpec[n].dY ); }

		char *GetCodeAt( int n ) { return( pSpec[n].pSpec ); }

		int GetAbundanceAt( int n ) { return( pSpec[n].Abundance ); }

		void DumpTableToFile( char * );

		~GDM_AbundanceRec();


	private:

		int nGetNumRecords( char * );

};


#endif // __CLSABUNDANCEREC_H__