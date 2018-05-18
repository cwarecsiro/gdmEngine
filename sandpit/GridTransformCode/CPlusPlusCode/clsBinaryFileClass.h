//
// clsBinaryFileClass.h
//
#ifndef __CLSBINARYFILECLASS_H__
#define __CLSBINARYFILECLASS_H__


#include <stdio.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <io.h>


#define BUFFMEG  1000000
#define BUFF1024 1024
#define BUFF512   512
#define BUFF256   256
#define BUFF64     64


//
// The states that a file can be opened in....
//
#define BFC_ReadOnly 0
#define BFC_ReadWrite 1
#define BFC_CreateOrTrunc 2


class BinaryFileClass
{
	int h;
	bool CanWrite;
	char FilePath[BUFF1024];

	public:

		BinaryFileClass(char *);

		BinaryFileClass(char *, int);

		char *GetFilePath() { return(FilePath); }

		bool IsValid();

		bool IsWritable();

		void Close();

		long CurrentOffset();

		long SeekTo(long n);

		long SeekForward(long n);

		long SeekBackward(long n);

		long SeekToStart();

		long SeekToEnd();

		long ReadByte(char *p, long n );

		long ReadShort(short *p, long n );	

		long ReadInt(int *p, long n );

		long ReadFloat(float *p, long n );

		long ReadDouble(double *p, long n );

		long WriteByte(char *p, long n );

		long WriteShort(short *p, long n );	

		long WriteInt(int *p, long n );

		long WriteFloat(float *p, long n );

		long WriteDouble(double *p, long n );

		bool CoordHasData(long lOffset, float fNoData);

		~BinaryFileClass();


	private:

};


#endif // __CLSBINARYFILECLASS_H__