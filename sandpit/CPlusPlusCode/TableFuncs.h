//
// TableFuncs.h
//
#ifndef	__TABLEFUNCS_H__
#define	__TABLEFUNCS_H__

#include "stdafx.h"
#include "GdmModel.h"

// zero based column indices for the site data
#define C_X0 1
#define C_Y0 2
#define C_X1 3
#define C_Y1 4

int GetNumPredictorsFromHeader(char *lpPath);

int GetNumRowsFromTableFile(char *lpPath);

double *GetColumnAt(char *lpPath, int nRows, int Index0);
double *GetColumnAt(char *lpPath, int Index0);

double *GetSiteADataForPredictorAt(char *lpPath, int nRows, int Index0);
double *GetSiteADataForPredictorAt(char *lpPath, int Index0);

double *GetSiteBDataForPredictorAt(char *lpPath, int nRows, int Index0);
double *GetSiteBDataForPredictorAt(char *lpPath, int Index0);

void DumpFirstNItems(int n, double *pData, char *Label);

double *GetGeographicISplines(int nRows, int nQuants, char *lpTable);
double *GetGeographicISplines(int nRows, GdmModel *gdmModel, char *lpTable);

double *GetPredictorISplines(int nRows, int nQuants, char *lpTablePath, int Index0);
double *GetPredictorISplines(int nRows, GdmModel *gdmModel, char *lpTablePath, int Index0);

double *GetGeographicDistanceVector(int nRows, char *lpTable);
double CalcGeographicDistance(double X0, double Y0, double X1, double Y1);

double *ExtractDataQuantiles(int nRows, int nQuants, double *pData);
double *ExtractDataQuantiles(int nRows, int nQuants, double *pData0, double *pData1);

int dcompare( const void *p0, const void *p1 );
double doSplineCalc( double dVal, double q1, double q2, double q3 );

#endif // __TABLEFUNCS_H__