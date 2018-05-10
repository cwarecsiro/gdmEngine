//
// FilterSites.h
//
#ifndef __FILTERSITES_H__
#define __FILTERSITES_H__

#include "clsEsriBinaryHeader.h"
#include "clsBinaryFileClass.h"
#include "myCallback.h"
#include "Message.h"

extern "C" 
{
	//
	// Filter input table thru domain and rewrite the filtered output to 
	// the Workspace with the table filename prefixed with "Filtered_"
	//
	_declspec(dllexport) bool FilterInputTable(char *pParams, FPTR fptr);

}


//
// Local Declarations
//
bool SiteHasData(double dX, double dY, EsriBinaryHeader *hdr, BinaryFileClass *bfc);
bool SiteHasData(double dX, double dY, float *pVal, EsriBinaryHeader *hdr, BinaryFileClass *bfc);
bool SiteInBoundingRectangle(double dX, double dY, EsriBinaryHeader *hdr);
bool FilterSiteTable(char *pParams, FPTR fptr);
bool FilterSitePairTable(char *pParams, FPTR fptr);
bool FilterSitePairDistanceTable(char *pParams, FPTR fptr);
char *ConstructFilteredTablePath(char *pParams);
int CountRows(char *pPath, FPTR fptr);
int CountColumns(char *pPath, FPTR fptr);
double AggregateX(double dX, EsriBinaryHeader *hdr);
double AggregateY(double dY, EsriBinaryHeader *hdr);
bool ValidCoordinate(double dX0, double dY0, double dX1, double dY1, 
	                 int nGrids, BinaryFileClass **ppBfc, EsriBinaryHeader **ppBfh);
long lGetCoordinateOffset(EsriBinaryHeader *header, double dX, double dY);

// for 64bit mode
int CountColumns64(char *pPath, FPTR fptr);
int CountRows64(char *pPath, FPTR fptr);

#endif // __FILTERSITES_H__