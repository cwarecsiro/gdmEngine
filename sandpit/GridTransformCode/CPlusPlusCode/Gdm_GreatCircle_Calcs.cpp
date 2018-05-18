//
// Gdm_GreatCircle_Calcs.cpp
//
#include "stdafx.h"
#include "Gdm_GreatCircle_Calcs.h"
#include "Message.h"
#include "ParamsW16.h"
#include "GDMBufferSizeDefs.h"

#include <math.h>

#define AVE_DIAM 6371.0

//
// convert degrees to radians
//
double d2r(double d)
{
	return(M_PI / 180.0 * d);
}


//
// calculate great circle distance with applying the spherical law of cosines (basic formula)
//
//double CalcGreatCircleDist(double dX0, double dY0, double dX1, double dY1)
//{
//	double dXAbs = fabs(dX0 - dX1);
//	double term1 = sin(d2r(dY0)) * sin(d2r(dY1));
//    double term2 = cos(d2r(dY0)) * cos(d2r(dY1)) * cos(d2r(dXAbs));
//    double dVal = acos(term1 + term2);
//    return (AVE_DIAM * dVal);
//}


//
// calculate great circle distance via euclidean distance Cosine version above NOT working properly
//
double CalcGreatCircleDist(double dX0, double dY0, double dX1, double dY1)
{
	double distX = dX0 - dX1;
	double distY = dY0 - dY1;
	return(sqrt((distX * distX) + (distY * distY)));
}