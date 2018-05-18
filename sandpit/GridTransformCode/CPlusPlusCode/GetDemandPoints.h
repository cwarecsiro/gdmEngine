//
// GetDemandPoints.h
//

#ifndef __GETDEMANDPOINTS_H__
#define __GETDEMANDPOINTS_H__

#include "stdafx.h"
#include "myCallback.h"

extern "C" 
{

_declspec(dllexport) bool GetDemandPoints_01( char *pParams, char *OutPath, int NumSamples, int NumDemandPoints, FPTR fptr);

}




#endif // __GETDEMANDPOINTS_H__