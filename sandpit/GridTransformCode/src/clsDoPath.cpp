//
// clsDoPath.cpp
//
#include "stdafx.h"
#include <sys/stat.h>
#include "clsDoPath.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This class is build on the _splitpath and _makepath functions defined in stdlib.h
//
//void _splitpath(
//   const char *path,
//   char *drive,
//   char *dir,
//   char *fname,
//   char *ext 
//);
//
//	Parameters:
//
//	path: Full path
//
//	drive : Optional drive letter, followed by a colon (:)
//
//	dir: Optional directory path, including trailing slash. Forward slashes ( / ), backslashes ( \ ), or both may be used.
//
//	fname: Base filename (no extension)
//
//	ext: Optional filename extension, including leading period (.)
//
//
//void _makepath(
//   char *path,
//   const char *drive,
//   const char *dir,
//   const char *fname,
//   const char *ext 
//);
//
//  Parameters - The following arguments point to buffers containing the path elements.
//
//	drive: Contains a letter (A, B, and so on) corresponding to the desired drive 
//		   and an optional trailing colon. _makepath inserts the colon automatically 
//		   in the composite path if it is missing. If drive is a null character or 
//	       an empty string, no drive letter and colon appear in the composite path string.
//
//	dir:   Contains the path of directories, not including the drive designator or 
//	       the actual file name. The trailing slash is optional, and either a forward slash 
//	       (/) or a backslash (\) or both might be used in a single dir argument. 
//	       If a trailing slash (/ or \) is not specified, it is inserted automatically. 
//	       If dir is a null character or an empty string, no slash is inserted in the composite path string.
//
//  fname: Contains the base file name without any file name extensions. 
//	       If fname is NULL or points to an empty string, no file name is inserted in the composite path string.
//
//  ext:   Contains the actual file name extension, with or without a leading period (.). 
//	       _makepath inserts the period automatically if it does not appear in ext. 
//	       If ext is a null character or an empty string, no period is inserted in the composite path string.
//
//		The path argument must point to an empty buffer large enough to hold the complete path. 
//		Although there are no size limits on any of the fields that constitute path, 
//		the composite path must be no larger than the _MAX_PATH constant, defined in Stdlib.h. 
//		_MAX_PATH might be larger than the current operating-system version can handle.
//
//
//  Buffer sizes....
//
//	#define _MAX_PATH   260 /* max. length of full pathname */
//	#define _MAX_DRIVE  3   /* max. length of drive component */
//	#define _MAX_DIR    256 /* max. length of path component */
//	#define _MAX_FNAME  256 /* max. length of file name component */
//	#define _MAX_EXT    256 /* max. length of extension component */
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


char Path[_MAX_PATH];
char Drive[_MAX_DRIVE];
char Dir[_MAX_DIR];
char Name[_MAX_FNAME];
char Ext[_MAX_EXT];


char *gmpath::GetName( char *thePath )
{
	_splitpath( thePath, Drive, Dir, Name, Ext );
	return( Name );
}


char *gmpath::GetDirectoryPath( char *thePath )
{
	_splitpath( thePath, Drive, Dir, Name, Ext );
	_makepath( Path, Drive, Dir, NULL, NULL );
	return( Path );
}


char *gmpath::GetDrive( char *thePath )
{
	_splitpath( thePath, Drive, Dir, Name, Ext );
	return( Drive );
}


char *gmpath::GetExtension(  char *thePath )
{
	_splitpath( thePath, Drive, Dir, Name, Ext );
	return( Ext );
}


char *gmpath::ChangeExtension(  char *thePath, char *theNewExt )
{
	_splitpath( thePath, Drive, Dir, Name, Ext );
	_makepath( Path, Drive, Dir, Name, theNewExt );
	return( Path );
}


char *gmpath::GetProfileString( char *App, char *Key, char *thePath )
{
	GetPrivateProfileString( App, Key, "X", Path, _MAX_PATH, thePath ); 
	return( Path );
}


int gmpath::GetProfileInt( char *App, char *Key, char *thePath )
{
	GetPrivateProfileString( App, Key, "X", Path, _MAX_PATH, thePath ); 
	return( atoi(Path) );
}


float gmpath::GetProfileFloat( char *App, char *Key, char *thePath )
{
	GetPrivateProfileString( App, Key, "X", Path, _MAX_PATH, thePath ); 
	return( (float)atof(Path) );
}


double gmpath::GetProfileDouble( char *App, char *Key, char *thePath )
{
	GetPrivateProfileString( App, Key, "X", Path, _MAX_PATH, thePath ); 
	return( atof(Path) );
}


bool gmpath::FileExists( char *thePath )
{
	// check that the file exists, if not then return false
#if defined _M_X64

	struct _stat64 statbuf;
	if ( 0 == _stat64( thePath, &statbuf ) )
	{
		return( true );
	}
	else
	{
		return( false );
	}

#elif defined _WIN32

	struct _stat statbuf;
	if ( 0 == _stat( thePath, &statbuf ) )
	{
		return( true );
	}
	else
	{
		return( false );
	}
#endif
}