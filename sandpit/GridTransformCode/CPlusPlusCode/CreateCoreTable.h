//
// CreateCoreTable.h
//
#ifndef __CREATECORETABLE_H__
#define __CREATECORETABLE_H__

#include "FilterSites.h"
#include "clsBinaryFileClass.h"
#include "clsEsriBinaryHeader.h"


//
// Return true if there are at least six comma delimited fields in the data
//
bool GotEnoughColsInComposite(char *InputPath);


//
// Return true if there are three comma delimited fields in the data or four if HaveAbundance
//
bool GotEnoughColsInSitePlusSpecies(char *InputPath, bool HaveAbundance);


//
// Adjust the sites to center them in their respective gridcells
//
void AggregateSitePair(double *pX0, double *pY0, double *pX1, double *pY1, EsriBinaryHeader *pHeader );


//
// Adjust the site to center in gridcell centroid
//
void AggregateSite(double *pX0, double *pY0, EsriBinaryHeader *pHeader );


//
// Create a version of the Site Plus Species Table that is Masked 
// by the mask grid and optionally aggregated to gridcell centroids
//
bool MaskSitePlusSpeciesTable(char *InputPath, 
	                          char *OutputPath, 
							  char *MaskPath, 
							  bool DoAggregation, 
							  bool HaveAbundance, 
							  FPTR fptr);


#endif // __CREATECORETABLE_H__