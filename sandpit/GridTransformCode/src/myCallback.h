//
// myCallback.h
//
#ifndef __MYCALLBACK_H__
#define __MYCALLBACK_H__

#include "stdafx.h"

typedef bool (CALLBACK *FPTR)( char *p, int i );

typedef bool (CALLBACK *UPTR)( double dExp, double dInc, double dIncCoeffSum, double dSum, int i);


#endif // __MYCALLBACK_H__