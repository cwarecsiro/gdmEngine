//
// GdmTransform.h
//
#ifndef __GDMTRANSFORM_H__
#define __GDMTRANSFORM_H__

bool TransposeGridsFromParamFile( char *pParams, FPTR fptr );

bool CreateEastingGrid( char *pParams, FPTR fptr );
bool CreateNorthingGrid( char *pParams, FPTR fptr );
bool CreateEnvPredictorGrids( char *pParams, FPTR fptr );

bool CreateSummedAbsERRGrid( char *pParams, FPTR fptr);
bool SumAbsoluteGridValues(char *pSumGridFLT, char *OutlierFLT, FPTR fptr);

bool TransformLookupTable(char *To, char *From, char *pParams, int nIndex);

bool bPredictorHasNonZeroCoeffs( char *pParams, int index );
bool GetGridAsDomain(char *pParams, char *pDomainPath );

bool InitBinaryGridFromDomain( char *pPath, char *pDomain, FPTR fptr );
bool InitBinaryGridFromGrid( char *pPath, char *pGrid, FPTR fptr );

bool ApplyEastWestGradient( char *pPath, char *pParams, FPTR fptr );
bool ApplyNorthSouthGradient( char *pPath, char *pParams, FPTR fptr );
bool ApplySplineGradient( char *pPath, char *pParams, int nIndex, FPTR fptr );

char *GetOutlierFilename(char *filename, char *ext);

double DoTranSplineCalc( double dVal, double q1, double q2, double q3 );


#endif // __GDMTRANSFORM_H__