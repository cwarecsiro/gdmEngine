#pragma once
//
// ESRIBinaryGridClass.h
//
#ifndef __ESRIBINARYGRIDCLASS_H__
#define __ESRIBINARYGRIDCLASS_H__

#include "clsBinaryFileClass.h"
#include "clsEsriBinaryHeader.h"

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

#pragma once
class ESRIBinaryGridClass
{
private:

	// used for the binary data file
	BinaryFileClass *binaryFile;

	// used for the grid metrics text file
	EsriBinaryHeader *metricsFile;

	// a grid row used for read data or data to write
	float *pRowData;

	// support methods
	void Convert2FLT(char *, char *);
	void Convert2HDR(char *, char *);


public:
	ESRIBinaryGridClass(char *);
	ESRIBinaryGridClass(char *, int);
	~ESRIBinaryGridClass(void);

	BinaryFileClass *GetBinary() { return(this->binaryFile); }
	EsriBinaryHeader *GetMetrics() { return(this->metricsFile); }

	int NumCols() { return(metricsFile->GetNumCols()); }
	int NumRows() { return(metricsFile->GetNumRows()); }
	float NoData() { return(metricsFile->GetNoDataValue()); }

	void ReadRow();
	void WriteRow();

	float *GetRowData() { return(pRowData); }
};

#endif // __ESRIBINARYGRIDCLASS_H__