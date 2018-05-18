//
// profileStrings.h
//
#include "stdafx.h"

#ifndef __PROFILESTRINGS_H__
#define __PROFILESTRINGS_H__

#define BUFFLEN 256
#define KBUFFLEN 1024
#define MBUFFLEN 1000000
#define PROFILE_ERROR -1000000


//
// set various values in the parameter file
//
void SetProfileString( char *lpSection, char *lpKey, char *lpString, char *lpParamFile );

void SetProfileInt( char *lpSection, char *lpKey, int nVal, char *lpParamFile );

void SetProfileDouble( char *lpSection, char *lpKey, double dVal, char *lpParamFile );


//
// check if profile string exists
//
bool ProfileStringExists( char *lpSection, char *lpKey, char *lpParamFile );


//
// extract various values from the parameter file
//
void GetProfileString( char *lpSection, char *lpKey, char *lpBuffer, char *lpParamFile );

int GetProfileInt( char *lpSection, char *lpKey, char *lpParamFile );

double GetProfileDouble( char *lpSection, char *lpKey, char *lpParamFile );



//
// delete sections from parameter file
//
void DeletePredictorSection( char *lpParamFile );

//
// insert a newline in the parameter file at the end
//
void InsertNewLine( char *lpParamFile );


#endif // __PROFILESTRINGS_H_
