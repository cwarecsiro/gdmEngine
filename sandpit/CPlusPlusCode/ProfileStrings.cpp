//
// profileStrings.cpp
//
#include "stdafx.h"
#include "ProfileStrings.h"
#include "Message.h"

#pragma warning( disable : 4996 )

void 
GetProfileString( char *lpSection, char *lpKey, char *lpBuffer, char *lpParamFile )
{
	GetPrivateProfileString( lpSection, lpKey, NULL, lpBuffer, BUFFLEN, lpParamFile );
}


int
GetProfileInt( char *lpSection, char *lpKey, char *lpParamFile )
{
	return( GetPrivateProfileInt( lpSection, lpKey, -1, lpParamFile ) );
}



double
GetProfileDouble( char *lpSection, char *lpKey, char *lpParamFile )
{
	char lpBuffer[BUFFLEN];
	GetPrivateProfileString( lpSection, lpKey, NULL, lpBuffer, BUFFLEN, lpParamFile );
	return( atof( lpBuffer ) );
}



//
// These methods check if the profile string exists
//
bool
ProfileStringExists( char *lpSection, char *lpKey, char *lpParamFile )
{
	char lpBuffer[BUFFLEN];
	GetPrivateProfileString( lpSection, lpKey, "0", lpBuffer, BUFFLEN, lpParamFile );
	return(0 == strcmp(lpBuffer, "0") ? false : true);
}



//
// if the lpSection argument does not exist in the parameter file, then it is created
// if the lpKey argument does not exist in the parameter file, then it is created
// if the lpKey argument is NULL, then the entire section including all entries in the section are deleted
// if the lpString argument is NULL, then the key pointer to by the lpKey value is deleted
//
void 
SetProfileString( char *lpSection, char *lpKey, char *lpString, char *lpParamFile )
{
	WritePrivateProfileString( lpSection, lpKey, lpString, lpParamFile );
}



void 
SetProfileInt( char *lpSection, char *lpKey, int nVal, char *lpParamFile )
{
	char lpBuffer[BUFFLEN];
	sprintf( lpBuffer, "%d", nVal );
	SetProfileString( lpSection, lpKey, lpBuffer, lpParamFile );
}



void 
SetProfileDouble( char *lpSection, char *lpKey, double dVal, char *lpParamFile )
{
	char lpBuffer[BUFFLEN];
	sprintf( lpBuffer, "%lf", dVal );
	SetProfileString( lpSection, lpKey, lpBuffer, lpParamFile );
}



void 
DeletePredictorSection( char *lpParamFile )
{
	SetProfileString( "PREDICTORS", NULL, "", lpParamFile );
}



void 
InsertNewLine( char *lpParamFile )
{
	FILE *fp = fopen( lpParamFile, "a+t" );
	fprintf( fp, "\n" );
	fclose( fp );
}


