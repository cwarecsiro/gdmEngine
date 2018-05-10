//
// clsDoPath.h
//
#ifndef __CLSDOPATH_H__
#define __CLSDOPATH_H__


class gmpath
{
	public: 

		char *GetName( char *thePath );

		char *GetDirectoryPath( char *thePath );

		char *GetDrive( char *thePath );

		char *GetExtension(  char *thePath );

		char *ChangeExtension(  char *thePath, char *theNewExt );

		char *GetProfileString( char *App, char *Key, char *thePath );

		int  GetProfileInt( char *App, char *Key, char *thePath );

		float GetProfileFloat( char *App, char *Key, char *thePath );

		double GetProfileDouble( char *App, char *Key, char *thePath );

		bool FileExists( char *thePath );
};



#endif // __CLSDOPATH_H__