//
// clsBinaryFileClass.cpp
//
#include "stdafx.h"

#include "clsBinaryFileClass.h"
#include "Message.h"



//
// Constructor for opening a binary file as Read-Only
//
BinaryFileClass::BinaryFileClass(char *p)
{
	h = _open( p, _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE );
	if ( h < 0 )
	{
		char q[512];
		sprintf( q, "BinaryFileClass() - Cannot open %s for read", p );
		Message( q, "ERROR" );
		return;
	}
	CanWrite = false;
	strcpy(FilePath, p);
}



//
// Constructor for opening a binary file according to a selected state
//
// The states are defined in BinaryFileClass.h
// BFC_ReadOnly, BFCReadWrite, BFC_CreateOrTrunc
//
BinaryFileClass::BinaryFileClass(char *p, int state)
{
	switch (state)
	{
		case BFC_ReadOnly:
			h = _open( p, _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE );
			if ( h < 0 )
			{
				char q[512];
				sprintf( q, "BinaryFileClass() - Cannot open %s for read", p );
				Message( q, "ERROR" );
			}
			CanWrite = false;
			break;


		case BFC_ReadWrite:
			h = _open( p, _O_BINARY | _O_RDWR, S_IREAD | S_IWRITE );
			if ( h < 0 )
			{
				char q[512];
				sprintf( q, "BinaryFileClass() - Cannot open %s for readwrite", p );
				Message( q, "ERROR" );
			}
			CanWrite = true;
			break;


		case BFC_CreateOrTrunc:
			h = _open( p, _O_BINARY | _O_CREAT | _O_TRUNC | _O_RDWR, S_IREAD | S_IWRITE );
			if ( h < 0 )
			{
				char q[512];
				sprintf( q, "BinaryFileClass() - Cannot create or truncate %s ", p );
				Message( q, "ERROR" );
			}
			CanWrite = true;
			break;


		default:
			h = -1;

	}

	strcpy(FilePath, p);
}



//
// Sanity Check. Do we have a valid file descriptor?
//
bool
BinaryFileClass::IsValid()
{
	return( h >= 0 );
}



//
// Returns true if the binary file is open for writing
//
bool
BinaryFileClass::IsWritable()
{
	return(CanWrite);
}



//
// Force file closure 
// 
// This is done in the destructor also if the file 
// is still open when the object leaves scope.
//
void 
BinaryFileClass::Close()
{
	if ( this->IsValid() )
	{
		_close(h);
	}
	h = -1;
}



//
// Get current byte offset from file start
//
long 
BinaryFileClass::CurrentOffset()
{
	if ( this->IsValid() )
		return(_lseek(h, 0L, SEEK_CUR));
	else
		return(-1);
}



//
// Seek to position n bytes from start
// and return current position
//
long
BinaryFileClass::SeekTo(long n)
{
	if ( this->IsValid() )
		return(_lseek(h, n, SEEK_SET));
	else
		return(-1);
}



//
// Seek forward n bytes from current position
// and return current position
//
long 
BinaryFileClass::SeekForward(long n)
{
	if ( this->IsValid() )
		return(_lseek(h, n, SEEK_CUR));
	else
		return(-1);
}



//
// Seek backward n bytes from current position
// and return current position
//
long 
BinaryFileClass::SeekBackward(long n)
{
	if ( this->IsValid() )
	{
		_lseek(h, -n, SEEK_CUR);
		return(_lseek(h, 0L, SEEK_CUR));
	}
	else
		return(-1);
}



// 
// Set file descriptor to the start of the binary file
//
long 
BinaryFileClass::SeekToStart()
{
	if ( this->IsValid() )
	{
		return(_lseek(h, 0L, SEEK_SET));
	}
	else
		return(-1);
}



// 
// Set file descriptor to the end of the binary file
//
long 
BinaryFileClass::SeekToEnd()
{
	if ( this->IsValid() )
	{
		return(_lseek(h, 0L, SEEK_END));
	}
	else
		return(-1);
}



//
// Read n bytes into buffer 
// and return number of bytes read
//
long 
BinaryFileClass::ReadByte(char *p, long n)
{
	if ( this->IsValid() )
	{
		return(_read(h, p, n));
	}
	else
		return(-1);
}



//
// Read n X 4 byte integer shorts into buffer 
// and return number of bytes read
//
long 
BinaryFileClass::ReadShort(short *p, long n)
{
	if ( this->IsValid() )
	{
		return(_read(h, p, n * sizeof(short)));
	}
	else
		return(-1);
}



//
// Read n X 4 byte integers into buffer 
// and return number of bytes read
//
long 
BinaryFileClass::ReadInt(int *p, long n)
{
	if ( this->IsValid() )
	{
		return(_read(h, p, n * sizeof(int)));
	}
	else
		return(-1);
}



//
// Read n X 4 byte floats into buffer 
// and return number of bytes read
//
long 
BinaryFileClass::ReadFloat(float *p, long n)
{
	if ( this->IsValid() )
	{
		return(_read(h, p, n * sizeof(float)));
	}
	else
		return(-1);
}



//
// Read n X 8 byte doubles into buffer 
// and return number of bytes read
//
long 
BinaryFileClass::ReadDouble(double *p, long n)
{
	if ( this->IsValid() )
	{
		return(_read(h, p, n * sizeof(double)));
	}
	else
		return(-1);
}



//
// Write n bytes from buffer 
// and return number of bytes written
//
long 
BinaryFileClass::WriteByte(char *p, long n)
{
	if ( !CanWrite )
	{
		Message( "Binary file is NOT writable", "ERROR" );
		return(-1);
	}

	if ( this->IsValid() )
	{
		return(_write(h, p, n));
	}
	else
		return(-1);
}



//
// Write n X 4 byte integer shorts from buffer 
// and return number of bytes written
//
long 
BinaryFileClass::WriteShort(short *p, long n)
{
	if ( !CanWrite )
	{
		Message( "Binary file is NOT writable", "ERROR" );
		return(-1);
	}

	if ( this->IsValid() )
	{
		return(_write(h, p, n * sizeof(short)));
	}
	else
		return(-1);
}



//
// Write n X 4 byte integers from buffer 
// and return number of bytes written
//
long 
BinaryFileClass::WriteInt(int *p, long n)
{
	if ( !CanWrite )
	{
		Message( "Binary file is NOT writable", "ERROR" );
		return(-1);
	}

	if ( this->IsValid() )
	{
		return(_write(h, p, n * sizeof(int)));
	}
	else
		return(-1);
}



//
// Write n X 4 byte floats from buffer 
// and return number of bytes written
//
long 
BinaryFileClass::WriteFloat(float *p, long n)
{
	if ( !CanWrite )
	{
		Message( "Binary file is NOT writable", "ERROR" );
		return(-1);
	}

	if ( this->IsValid() )
	{
		return(_write(h, p, n * sizeof(float)));
	}
	else
		return(-1);
}



//
// Write n X 8 byte doubles from buffer 
// and return number of bytes written
//
long 
BinaryFileClass::WriteDouble(double *p, long n)
{
	if ( !CanWrite )
	{
		Message( "Binary file is NOT writable", "ERROR" );
		return(-1);
	}

	if ( this->IsValid() )
	{
		return(_write(h, p, n * sizeof(double)));
	}
	else
		return(-1);
}



//
// Destructor
//
BinaryFileClass::~BinaryFileClass()
{
	if ( h >= 0 )
	{
		_close(h);
	}
}



//
// return true if gridcell at x,y has valid data (ie NOT NoData)
//
bool 
BinaryFileClass::CoordHasData(long lOffset, float fNoData)
{
	if (lOffset < 0L)
		return(false);

	this->SeekTo(lOffset * sizeof(float));
	float f;
	this->ReadFloat(&f, 1L);
	return(f != fNoData);
}

