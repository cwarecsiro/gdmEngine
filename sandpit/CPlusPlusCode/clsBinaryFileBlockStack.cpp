//
// clsBinaryFileBlockStack.cpp
//
#include "stdafx.h"

#include "clsBinaryFileBlockStack.h"
#include "Message.h"
#include "ProfileStrings.h"
#include "clsDoPath.h"


//
// ctor
//
BinaryFileBlockStack::BinaryFileBlockStack(char *pGDMParamFile, int PredGridType, int BlockSize)
{
	//
	// Firstly, do a sanity check that the GDM model only has environmental grids (no categorical or covariance)
	//
	char lpKey[64];
	int n = GetProfileInt( "PREDICTORS", "NumPredictors", pGDMParamFile );
	for ( int nx=1; nx<=n; nx++ )
	{
		sprintf(lpKey, "PredType%d", nx);
		if (GetProfileInt( "PREDICTORS", lpKey, pGDMParamFile ))
		{
			Message("GDM has covariance or categorical grids, cannot use", "BinaryFileStack");
			return;
		}
	}
		
	//
	// determine the number of valid grids
	//
	char pPath[BUFF512]; 
	nGrids=0;
	for ( int nx=1; nx<=n; nx++ )
	{
		if (PredGridType == BFS_PredGrids)
		{
			sprintf(lpKey, "UseEnv%d", nx);
			if (1 == GetProfileInt("PREDICTORS", lpKey, pGDMParamFile))
			{
				++nGrids;
			}
		}

		else if (PredGridType == BFS_TranGrids)
		{
			sprintf(lpKey, "PredTran%d", nx);
			if (ProfileStringExists("TRANSPREDS", lpKey, pGDMParamFile))
			{
				++nGrids;
			}
		}

		else
		{
			Message("Invalid BFS Grid Type", "ERROR");
			return;
		}
	}
	if (nGrids < 1)
	{
		Message("Got NO grids to use", "ERROR");
		return;
	}

	//
	// Open the Binary grids
	//
	gmpath gmPath;
	pHeader = NULL;
	ppGrids = new BinaryFileClass * [nGrids];
	int nThis = 0;
	for ( int nx=1; nx<=n; nx++ )
	{
		if (PredGridType == BFS_PredGrids)
		{
			sprintf(lpKey, "UseEnv%d", nx);
			if (1 == GetProfileInt("PREDICTORS", lpKey, pGDMParamFile))
			{
				sprintf(lpKey, "EnvGrid%d", nx);
				GetProfileString("PREDICTORS", lpKey, pPath, pGDMParamFile);
				ppGrids[nThis] = new BinaryFileClass(gmPath.ChangeExtension(pPath, ".flt"));
				if (!ppGrids[nThis]->IsValid())
				{
					Message("Invalid Binary Grid File arg", "ERROR");
					return;
				}
				
				// use the first grid for the header
				if (0 == nThis)
				{
					pHeader= new EsriBinaryHeader(gmPath.ChangeExtension(pPath, ".hdr"));
					nRows = pHeader->GetNumRows();
					nCols = pHeader->GetNumCols();
				}
				++nThis;
			}
		}

		else if (PredGridType == BFS_TranGrids)
		{
			sprintf(lpKey, "PredTran%d", nx);
			if (ProfileStringExists("TRANSPREDS", lpKey, pGDMParamFile))
			{
				GetProfileString("TRANSPREDS", lpKey, pPath, pGDMParamFile);
				ppGrids[nThis] = new BinaryFileClass(gmPath.ChangeExtension(pPath, ".flt"));
				if (!ppGrids[nThis]->IsValid())
				{
					Message("Invalid Binary Grid File arg", "ERROR");
					return;
				}
				// use the first grid for the header
				if (0 == nThis)
				{
					pHeader= new EsriBinaryHeader(gmPath.ChangeExtension(pPath, ".hdr"));
					nRows = pHeader->GetNumRows();
					nCols = pHeader->GetNumCols();
				}
				++nThis;
			}
		}
	}

	//
	// set the read block size and calculate the block remainder
	//
	nBlockSize = BlockSize;
	nNumBlocks = nRows / nBlockSize;
	nBlockRemainder = nRows % nBlockSize;

	//
	// allocate the row vectors
	//
	ppRowData = new float * [nGrids];
	for ( int i=0; i<nGrids; i++ )
	{
		ppRowData[i] = new float [nCols * nBlockSize];
	}
	pCellData = new float [nGrids];
}


//
// Get the next row of nBlockSize nCols from the current file position
//
void BinaryFileBlockStack::GetNextRow()
{
	if (!InsideRemainder())
	{
		for (int i=0; i<nGrids; i++)
			ppGrids[i]->ReadFloat(ppRowData[i], nCols * nBlockSize);
	}
	else
	{
		for (int i=0; i<nGrids; i++)
			ppGrids[i]->ReadFloat(ppRowData[i], nCols * nBlockRemainder);
	}
}


//
// get the cell vector for the nIndex column
//
void BinaryFileBlockStack::GetCellAt(int nBlockOffset, int nIndex)
{
	for (int i=0; i<nGrids; i++)
	{
		pCellData[i] = ppRowData[i][(nBlockOffset*nCols)+nIndex];
	}
}


//
// get the GDM dissimilarity between this cell vector and another cell vector
//
double BinaryFileBlockStack::GetGDMDissimilarity(float *pCompVec)
{
	double dDist = 0.0;

	for ( int i=0; i<nGrids; i++ )
	{
		dDist += (double)fabs(pCellData[i] - pCompVec[i]);
	}
	return(1.0 - exp(-dDist));
}


//
// copy the header
//
void BinaryFileBlockStack::CopyHeaderTo(char *pPath)
{
	pHeader->CopyTo(pPath);
}


bool BinaryFileBlockStack::InsideRemainder() 
{
	return(ppGrids[0]->CurrentOffset() / sizeof(float) / nCols >= nRows - nBlockRemainder);
}


//
// dtor
//
BinaryFileBlockStack::~BinaryFileBlockStack()
{
	for ( int i=0; i<nGrids; i++ )
	{
		if (ppRowData[i]) delete[] ppRowData[i];
		if (ppGrids[i]->IsValid()) ppGrids[i]->Close();
	}

	if (ppRowData) delete[] ppRowData;
	if (pCellData) delete[] pCellData;
	if (pHeader) delete pHeader;
	if (ppGrids) delete[] ppGrids;
}

