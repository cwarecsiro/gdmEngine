//
// GridMinMax.cpp
//
#include "stdafx.h"
#include "GdmBinLib.h"
#include "GridMinMax.h"
#include "clsEsriBinaryHeader.h"
#include "clsBinaryFileClass.h"
#include "clsDoPath.h"
#include "Message.h"


//
// Extract the Minimum and Maximum values from an ESRI Binary Export File
//
bool GetGridMinMax(char *GridPath, double *pMin, double *pMax, FPTR fptr)
{
	gmpath myGmPath;

	EsriBinaryHeader *ebh = new EsriBinaryHeader(myGmPath.ChangeExtension(GridPath, ".hdr")); 
	int nRows = ebh->GetNumRows();
	int nCols = ebh->GetNumCols();
	float fNoData = ebh->GetNoDataValue();
	if (ebh) delete ebh;

	BinaryFileClass *bfc = new BinaryFileClass(myGmPath.ChangeExtension(GridPath, ".flt"));
	if (!bfc->IsValid())
	{
		Message("Invalid Binary File", "GetGridMinMax()");
		return(false);
	}

	float *pRow = new float [nCols];

	float fMin = 100000000.0F;
	float fMax = -100000000.0F;
	int nCurrent = 0;
	fptr("Extracting Min and Max quantiles from grid", nCurrent );
	for ( int i=0; i<nRows; i++ )
	{
		if ( i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			fptr("Extracting Min and Max quantiles from grid", nCurrent );
		}

		bfc->ReadFloat(pRow, nCols);

		for ( int j=0; j<nCols; j++ )
		{
			if (fNoData != pRow[j])
			{
				// adjust min max vals...
				if ( pRow[j] < fMin ) fMin = pRow[j];
				if ( pRow[j] > fMax ) fMax = pRow[j];
			}
		}
	}
	fptr("Extracting Min and Max quantiles from grid", 0 );
	bfc->Close();
	if (bfc) delete bfc;
	*pMin = (double)fMin;
	*pMax = (double)fMax;
	if (pRow) delete[] pRow;
	return(true);
}

