//
// GDM_GreatCircle_Calcs.h
//
#ifndef __GDM_GREATCIRCLE_CALCS_H__
#define __GDM_GREATCIRCLE_CALCS_H__

#define M_PI       3.14159265358979323846

// convert degrees to radians
double d2r(double d);

// calculate great circle distance with applying the spherical law of cosines (basic formula)
double CalcGreatCircleDist(double dX0, double dY0, double dX1, double dY1);

#endif // __GDM_GREATCIRCLE_CALCS_H__