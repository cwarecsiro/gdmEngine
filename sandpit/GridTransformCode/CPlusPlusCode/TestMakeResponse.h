//
// TestMakeResponse.h
//
#ifndef __TESTMAKERESPONSE_H__
#define __TESTMAKERESPONSE_H__

#include "DistanceCalcs.h"
#include "myCallback.h"


extern "C"
{
	_declspec(dllexport) bool SitePlusSpeciesToResponsePresAbs(char *pIn, char *pResponse, char *pWeights, FPTR fptr);

	_declspec(dllexport) bool SitePlusSpeciesToResponseAbundance(char *pIn, char *pResponse, char *pWeights, FPTR fptr);

	_declspec(dllexport) bool SitePlusSpeciesToSites(char *pIn, char *pSites, FPTR fptr);
}

#endif // __TESTMAKERESPONSE_H__