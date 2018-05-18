//
// GDMClassification.h
//
#ifndef __GDMCLASSIFICATION_H__
#define __GDMCLASSIFICATION_H__


typedef struct tagDemandSiteRecord {
	double dDist;
	double dX;
	double dY;
} DemandSiteRecord;



bool CreateSamplePointMesh( char *lpParams, char *lpSamplePointPath, char *lpDomainPath, int nSamples, bool DoBatch, FPTR fptr );

int NRound( float f );

double **ExtractDataFromTransformGrids( char *lpParams, char *lpDomainPath, char *lpSamplePath, int *pRows, int *pCols, bool DoBatch );

void DumpANNFile( char *lpParams, int *pClusterID, 
	              char *lpSamplePath, double *pSampleX, double *pSampleY, 
				  double **ppData, int nRows, int nCols );

void WriteDendrogramTable(char *lpParams, int *pClusterID, char *lpSamplePath, double **ppData, int nRows, int nCols);

bool CreateKMeansClosestPointTable(char *pParams, char *lpClassPath,
	                               char *lpInTablePath, char *lpOutTablePath, 
	                               int nClasses, FPTR fptr);

bool CreateKMeansDemandPointTable(char *pParams, char *lpClassPath,
                                  char *lpInTablePath, char *lpOutTablePath,
                                  int nClasses, FPTR fptr);

int GetNumGDMTransformedPreds(char *pParams);

char **GetTransformPaths(char *pParams);

void GetSampleCoords(char *lpSamplePath, double *pSampleX, double *pSampleY, int nRows);

int *GetHierachicalClusterClasses( double **ppData, int nRows, int nCols, int nClasses, bool DoBatch, FPTR fptr );

bool DoGridClassification( char *lpParams, 
						   char *lpSubDomain,
						   int *pClusterID, double **ppData, 
						   int nRows, int nCols, 
						   char *lpOutPath,
						   bool DoBatch, 
						   FPTR fptr );

bool DoClassificationColoring( char *lpParams, char *lpSamplePath, int nClasses, char *GridName, bool DoBatch, FPTR fptr );

bool PopulateClassColor( char *lpSamplePath, double ***pppData, int **ppClass, int *pRows, int *pCols );

double CalcManhattanDistance( double **ppClassAverages, int nCols, int x1, int x2 );

void CreateRGBGridLegends( char *ClassPath, int **ppRGB, int Version, int nClasses, char *GridName );

int nCountANNFileRows( char *p );

int nCountANNFileColumns( char *p );

bool ExtractDataFromTable(char *TablePath, int nRows, int nCols, int *pClusterID, double **pData, FPTR fptr);

bool DoPostClassificationMatrix(char *pParams, 
	                            char *lpDomainPath, 
								char *lpOutName, 
								char *lpSamplePath, 
								bool DoBatch, 
								FPTR fptr);

int dCompare( const void *arg1, const void *arg2 );

int *GetClassIDs(double *pSampleX, double *pSampleY, int nRows, char *pParams, char *lpOutName);


bool DoDemandPoints_01(char *pParams, 
	                   char *lpDomainPath, 
	                   char *lpOutName, 
					   char *lpSamplePath, 
					   bool DoBatch, 
					   FPTR fptr);


bool DoDemandPoints_02(char *pParams, 
	                   char *lpDomainPath, 
	                   char *lpOutName, 
					   char *lpSamplePath, 
					   bool DoBatch, 
					   FPTR fptr);


void DumpDemandPointFile( char *lpParams, char *lpSamplePath, 
	                      double *pSampleX, double *pSampleY, 
				          double **ppData, int nRows, int nCols );


double CalculateGDMTransform( double dVal, int nSplines, double *pQuantiles, double *pCoeffs );


#endif // __GDMCLASSIFICATION_H__