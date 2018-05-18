//
// ESRIBinaryGridClass.cpp
//
#include "stdafx.h"
#include "ESRIBinaryGridClass.h"
#include "Message.h"


//
// Constructor for opening a binary file as Read-Only
//
ESRIBinaryGridClass::ESRIBinaryGridClass(char *p)
{
	//Init
	metricsFile = NULL;
	binaryFile = NULL;
	char Error_Message[1024];

	//
	// open the text file containing the grid metrics...
	//
	char hdrPath[1024];
	Convert2HDR(p, hdrPath);
	metricsFile = new EsriBinaryHeader(hdrPath);
	if (NULL == metricsFile)
	{
		sprintf(Error_Message, "ESRIBinaryGridClass() - Cannot open %s for read", hdrPath);
		Message(Error_Message, "ERROR");
		return;
	}

	//
	// open the binary data file...
	//
	char fltPath[1024];
	Convert2FLT(p, fltPath);
	binaryFile = new BinaryFileClass(fltPath);
	if (NULL == binaryFile)
	{
		sprintf(Error_Message, "ESRIBinaryGridClass() - Cannot open %s for read", fltPath);
		Message(Error_Message, "ERROR");
		return;
	}

	//
	// allocate buffers
	//
	pRowData = new float[metricsFile->GetNumCols()];
	if (NULL == pRowData)
	{
		sprintf(Error_Message, "ESRIBinaryGridClass() - Cannot allocate pRowData");
		Message(Error_Message, "ERROR");
		return;
	}
}


//
// Constructor for opening a binary file according to a selected state
//
// The states are defined in BinaryFileClass.h
// BFC_ReadOnly, BFCReadWrite, BFC_CreateOrTrunc
//
ESRIBinaryGridClass::ESRIBinaryGridClass(char *p, int state)
{
	//Init
	metricsFile = NULL;
	binaryFile = NULL;
	char Error_Message[1024];

	//
	// open the text file containing the grid metrics...
	//
	char hdrPath[1024];
	Convert2HDR(p, hdrPath);
	metricsFile = new EsriBinaryHeader(hdrPath);
	if (NULL == metricsFile)
	{
		sprintf(Error_Message, "ESRIBinaryGridClass() - Cannot open %s for read", hdrPath);
		Message(Error_Message, "ERROR");
		return;
	}

	//
	// open the binary data file...
	//
	char fltPath[1024];
	Convert2FLT(p, fltPath);
	binaryFile = new BinaryFileClass(fltPath, state);
	if (!binaryFile->IsValid())
	{
		switch (state)
		{
		case BFC_ReadOnly:
			sprintf(Error_Message, "ESRIBinaryGridClass() - Cannot open %s for read", p);
			Message(Error_Message, "ERROR");
			return;

		case BFC_ReadWrite:
			sprintf(Error_Message, "ESRIBinaryGridClass() - Cannot open %s for readwrite", p);
			Message(Error_Message, "ERROR");
			return;

		case BFC_CreateOrTrunc:
			sprintf(Error_Message, "ESRIBinaryGridClass() - Cannot create or truncate %s ", p);
			Message(Error_Message, "ERROR");
			return;
		}
	}

	//
	// allocate buffers
	//
	pRowData = new float[metricsFile->GetNumCols()];
	if (NULL == pRowData)
	{
		sprintf(Error_Message, "ESRIBinaryGridClass() - Cannot allocate pRowData");
		Message(Error_Message, "ERROR");
		return;
	}
}


//
// Destructor
//
ESRIBinaryGridClass::~ESRIBinaryGridClass(void)
{
	if (metricsFile) delete metricsFile;
	if (binaryFile) delete binaryFile;
	if (pRowData) delete[] pRowData;
}


//
// Convert a file path string to have a .FLT extension
//
void ESRIBinaryGridClass::Convert2FLT(char *thePath, char *theReturnedPath)
{
	char Drive[_MAX_DRIVE];
	char Dir[_MAX_DIR];
	char Name[_MAX_FNAME];
	char Ext[_MAX_EXT];

	_splitpath(thePath, Drive, Dir, Name, Ext);
	_makepath(theReturnedPath, Drive, Dir, Name, ".flt");
}


//
// Convert a file path string to have a .HDR extension
//
void ESRIBinaryGridClass::Convert2HDR(char *thePath, char *theReturnedPath)
{
	char Drive[_MAX_DRIVE];
	char Dir[_MAX_DIR];
	char Name[_MAX_FNAME];
	char Ext[_MAX_EXT];

	_splitpath(thePath, Drive, Dir, Name, Ext);
	_makepath(theReturnedPath, Drive, Dir, Name, ".hdr");
}


void ESRIBinaryGridClass::ReadRow()
{
	binaryFile->ReadFloat(pRowData, NumCols());
}


void ESRIBinaryGridClass::WriteRow()
{
	binaryFile->WriteFloat(pRowData, NumCols());
}


