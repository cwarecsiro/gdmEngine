//
// CompositeJPEG.cpp
//
#include "stdafx.h"
#include "GdmBinLib.h"
#include "Message.h"
#include "clsDoPath.h"
#include "clsBinaryFileClass.h"
#include "clsEsriBinaryHeader.h"



bool RGBGridsOK(char *RGrid, char *GGrid, char *BGrid, FPTR fptr)
{
	//Message("Entering RGBGridsOK", "INFO");
	gmpath gmPath;

	EsriBinaryHeader *headerRed = new EsriBinaryHeader(gmPath.ChangeExtension(RGrid, ".hdr"));
	EsriBinaryHeader *headerGreen = new EsriBinaryHeader(gmPath.ChangeExtension(GGrid, ".hdr"));
	EsriBinaryHeader *headerBlue = new EsriBinaryHeader(gmPath.ChangeExtension(BGrid, ".hdr"));

	//
	// Cell Size
	//
	if ((headerRed->GetCellSize() != headerGreen->GetCellSize()) ||
	    (headerRed->GetCellSize() != headerBlue->GetCellSize()))
	{
		if (headerRed) delete headerRed;  
		if (headerGreen) delete headerGreen;
		if (headerBlue) delete headerBlue;
		Message("Different Cell Sizes in RGBGridsOK", "ERROR");
		return(false);
	}
	
	//
	// Row Count
	//
	if ((headerRed->GetNumRows() != headerGreen->GetNumRows()) ||
	    (headerRed->GetNumRows() != headerBlue->GetNumRows()))
	{		
		if (headerRed) delete headerRed;  
		if (headerGreen) delete headerGreen;
		if (headerBlue) delete headerBlue;
		Message("Different Row Counts in RGBGridsOK", "ERROR");
		return(false);
	}

	//
	// Column Count
	//
	if ((headerRed->GetNumCols() != headerGreen->GetNumCols()) ||
	    (headerRed->GetNumCols() != headerBlue->GetNumCols()))
	{
		if (headerRed) delete headerRed;  
		if (headerGreen) delete headerGreen;
		if (headerBlue) delete headerBlue;
		Message("Different Column Counts in RGBGridsOK", "ERROR");		
		return(false);
	}

	//
	// Left Edge
	//
	if ((headerRed->GetMinX() != headerGreen->GetMinX()) ||
	    (headerRed->GetMinX() != headerBlue->GetMinX()))
	{		
		if (headerRed) delete headerRed;  
		if (headerGreen) delete headerGreen;
		if (headerBlue) delete headerBlue;
		Message("Different Left Edge Coordinates in RGBGridsOK", "ERROR");
		return(false);
	}

	//
	// Right Edge
	//
	if ((headerRed->GetMaxX() != headerGreen->GetMaxX()) ||
	    (headerRed->GetMaxX() != headerBlue->GetMaxX()))
	{		
		if (headerRed) delete headerRed;  
		if (headerGreen) delete headerGreen;
		if (headerBlue) delete headerBlue;
		Message("Different Right Edge Coordinates in RGBGridsOK", "ERROR");
		return(false);
	}

	//
	// Top Edge
	//
	if ((headerRed->GetMinY() != headerGreen->GetMinY()) ||
	    (headerRed->GetMinY() != headerBlue->GetMinY()))
	{		
		if (headerRed) delete headerRed;  
		if (headerGreen) delete headerGreen;
		if (headerBlue) delete headerBlue;
		Message("Different Top Edge Coordinates in RGBGridsOK", "ERROR");
		return(false);
	}

	//
	// Bottom Edge
	//
	if ((headerRed->GetMaxY() != headerGreen->GetMaxY()) ||
	    (headerRed->GetMaxY() != headerBlue->GetMaxY()))
	{		
		if (headerRed) delete headerRed;  
		if (headerGreen) delete headerGreen;
		if (headerBlue) delete headerBlue;
		Message("Different Bottom Edge Coordinates in RGBGridsOK", "ERROR");
		return(false);
	}

	if (headerRed) delete headerRed;  
	if (headerGreen) delete headerGreen;
	if (headerBlue) delete headerBlue;
	return(true);
}



//
// Extract the minimum and maximum valid data values
//
bool GetGridMetrics(char *pPath, double *pMin, double *pMax, FPTR fptr)
{
	//Message("Entering GetGridMetrics", "INFO");
	gmpath gmPath;

	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(pPath, ".hdr"));
	BinaryFileClass *bfc_In = new BinaryFileClass(gmPath.ChangeExtension(pPath, ".flt"), BFC_ReadOnly);	
	if (!bfc_In->IsValid())
	{
		Message("Invalid Binary File Class in GetGridMetrics", "ERROR");
		if (header) delete header;
		return(false);
	}

	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	float fNoData = header->GetNoDataValue();
	float *pData = new float [nCols];

	double dMin = *pMin;
	double dMax = *pMax;

	int nCurrent = 0;
	fptr("Getting Grid Min/Max...", nCurrent);	
	for ( int i=0; i<nRows; i++ )
	{
		if ( i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			fptr("Getting Grid Min/Max...", nCurrent);	
		}

		// read a row
		bfc_In->ReadFloat( pData, nCols );

		for ( int j=0; j<nCols; j++ )
		{
			if ( pData[j] != fNoData ) 
			{
				double dVal = double(pData[j] * 1.0);
				if (dVal < dMin) dMin = dVal;
				if (dVal > dMax) dMax = dVal;
			}
		}
	}

	// reset min/max vals
	*pMin = dMin;
	*pMax = dMax;

	if (pData) delete[] pData;
	if (header) delete header;
	return(true);
}



//
// Merge and range standardise grids into a 24bbp binary file
//
bool CreateCompositeFile(char *RedPath, double dRedMin, double dRedMax,
                         char *GreenPath, double dGreenMin, double dGreenMax,
                         char *BluePath, double dBlueMin, double dBlueMax, 
                         char *CompPath, FPTR fptr)
{
	//Message("Entering CreateCompositeFile", "INFO");
	gmpath gmPath;
	EsriBinaryHeader *headerRed = new EsriBinaryHeader(gmPath.ChangeExtension(RedPath, ".hdr"));
	EsriBinaryHeader *headerGreen = new EsriBinaryHeader(gmPath.ChangeExtension(GreenPath, ".hdr"));
	EsriBinaryHeader *headerBlue = new EsriBinaryHeader(gmPath.ChangeExtension(BluePath, ".hdr"));

	int nRows = headerRed->GetNumRows();
	int nCols = headerRed->GetNumCols();
	float fRedNoData = headerRed->GetNoDataValue();
	float fGreenNoData = headerGreen->GetNoDataValue();
	float fBlueNoData = headerBlue->GetNoDataValue();

	double dRedRange = dRedMax - dRedMin;
	double dGreenRange = dGreenMax - dGreenMin;
	double dBlueRange = dBlueMax - dBlueMin;

	float *pRedData = new float [nCols];
	float *pGreenData = new float [nCols];
	float *pBlueData = new float [nCols];
	char  *pOutData = new char [nCols * 3];

	BinaryFileClass *bfc_Red = new BinaryFileClass(gmPath.ChangeExtension(RedPath, ".flt"), BFC_ReadOnly);	
	BinaryFileClass *bfc_Green = new BinaryFileClass(gmPath.ChangeExtension(GreenPath, ".flt"), BFC_ReadOnly);	
	BinaryFileClass *bfc_Blue = new BinaryFileClass(gmPath.ChangeExtension(BluePath, ".flt"), BFC_ReadOnly);	
	BinaryFileClass *bfc_Out = new BinaryFileClass(gmPath.ChangeExtension(CompPath, ".rgb"), BFC_CreateOrTrunc);	

	int nCurrent = 0;
	fptr("Creating Composite...", nCurrent);	
	for ( int i=0; i<nRows; i++ )
	{
		if ( i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			fptr("Creating Composite...", nCurrent);	
		}

		bfc_Red->ReadFloat(pRedData, nCols);
		bfc_Green->ReadFloat(pGreenData, nCols);
		bfc_Blue->ReadFloat(pBlueData, nCols);

		for ( int j=0; j<nCols; j++ )
		{
			if ((pRedData[j] == fRedNoData)  || (dRedRange == 0.0))
			{
				pOutData[(j*3) + 2] = 0x00;
			}

			else  
			{
					pOutData[(j*3) + 2] = char(((double)pRedData[j] - dRedMin) / dRedRange * 255);
			}


			if ((pGreenData[j] == fGreenNoData) || (dGreenRange == 0.0))
			{
				pOutData[(j*3) + 1] = 0x00;
			}

			else
			{
				pOutData[(j*3) + 1] = char(((double)pGreenData[j] - dGreenMin) / dGreenRange * 255);
			}


			if ((pBlueData[j] == fBlueNoData) || (dBlueRange == 0.0))
			{
				pOutData[(j*3) + 0] = 0x00;
			}

			else
			{
				pOutData[(j*3) + 0] = char(((double)pBlueData[j] - dBlueMin) / dBlueRange * 255);
			}
		}

		bfc_Out->WriteByte((char*)pOutData, nCols*3);
	}
	
	bfc_Red->Close();
	bfc_Green->Close();
	bfc_Blue->Close();
	bfc_Out->Close();
	if (headerRed) delete headerRed;  
	if (headerGreen) delete headerGreen;
	if (headerBlue) delete headerBlue;

	if (pRedData) delete[] pRedData;
	if (pGreenData) delete[] pGreenData;
	if (pBlueData) delete[] pBlueData;
	if (pOutData) delete[] pOutData;
	return(true);
}



//
// Get Row and Column Count form an ESRI Binary File
//
bool GetRowColCount(char *pPath, int *pRows, int *pCols)
{
	//Message("Entering GetRowColCount", "INFO");
	gmpath gmPath;
	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(pPath, ".hdr"));
	*pRows = header->GetNumRows();
	*pCols = header->GetNumCols();
	if (header) delete header;
	return(true);
}

