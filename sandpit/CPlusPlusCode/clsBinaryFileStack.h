//
// clsBinaryFileStack.h
//
#ifndef __CLSBINARYFILESTACK_H__
#define __CLSBINARYFILESTACK_H__

#include "clsEsriBinaryHeader.h"
#include "clsBinaryFileClass.h"

#define BFS_PredGrids 0		// [PREDICTORS]
#define BFS_TranGrids 1		// [TRANSPREDS]

class BinaryFileStack
{
	int nRows;
	int nCols;
	int nGrids;
	float **ppRowData;
	float *pCellData;
	EsriBinaryHeader *pHeader;
	BinaryFileClass **ppGrids;

public: 

	BinaryFileStack(char *pGDMParamFile, int PredGridType);

	int GetNumRows() { return(nRows); }
	int GetNumCols() { return(nCols); }
	int GetNumGrids() { return(nGrids); }

	float GetNoDataVal() { return(pHeader->GetNoDataValue()); }
	float **GetRowData() { return(ppRowData); }
	float *GetCellData() { return(&pCellData[0]); }

	void CopyHeaderTo(char *pPath);

	void GetNextRow();
	void GetRowAt(int nIndex);
	void GetCellAt(int nIndex);

	double GetGDMDissimilarity(float *pCompVec);
	double GetManhattanDist(float *pCompVec);

	~BinaryFileStack();

};


#endif // __CLSBINARYFILESTACK_H__