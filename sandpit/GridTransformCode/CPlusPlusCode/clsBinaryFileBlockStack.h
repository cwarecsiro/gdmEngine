//
// clsBinaryFileBlockStack.h
//
#ifndef __CLSBINARYFILEBLOCKSTACK_H__
#define __CLSBINARYFILEBLOCKSTACK_H__

#include "clsEsriBinaryHeader.h"
#include "clsBinaryFileClass.h"

#define BFS_PredGrids 0		// [PREDICTORS]
#define BFS_TranGrids 1		// [TRANSPREDS]

class BinaryFileBlockStack
{
	int nRows;
	int nCols;
	int nGrids;
	int nBlockSize;
	int nNumBlocks;
	int nBlockRemainder;

	float **ppRowData;
	float *pCellData;
	EsriBinaryHeader *pHeader;
	BinaryFileClass **ppGrids;	

public: 

	BinaryFileBlockStack(char *pGDMParamFile, int PredGridType, int BlockSize);

	int GetNumRows() { return(nRows); }
	int GetNumCols() { return(nCols); }
	int GetNumGrids() { return(nGrids); }
	int GetBlockSize() { return(nBlockSize); }
	int GetNumBlocks() { return(nNumBlocks); }
	int GetBlockRemainder() { return(nBlockRemainder); }

	float GetNoDataVal() { return(pHeader->GetNoDataValue()); }
	float **GetRowData() { return(ppRowData); }
	float *GetCellData() { return(&pCellData[0]); }

	void CopyHeaderTo(char *pPath);

	void GetNextRow();
	void GetCellAt(int nBlockOffset, int nIndex);

	double GetGDMDissimilarity(float *pCompVec);

	bool InsideRemainder();

	~BinaryFileBlockStack();

};

#endif // __CLSBINARYFILEBLOCKSTACK_H__