//
// gdmUtilities.h
//
#ifndef __GDMUTILITIES_H__
#define __GDMUTILITIES_H__

#include "myCallback.h"

extern "C"
{
	_declspec(dllexport) int CountRowsInResponse(char *pParams, bool UseFiltered, FPTR fptr);

}

#endif