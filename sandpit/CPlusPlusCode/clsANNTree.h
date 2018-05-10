//
// clsANNTree.h
//
#ifndef __CLSANNTREE_H__
#define __CLSANNTREE_H__

#include "ANN.h"

class ANNTree
{
	int nRows;
	int nCols;

	ANNpointArray pTableData;
	ANNkd_tree *theTree;
	ANNidxArray nn_idx;
	ANNdistArray dist;
	ANNpoint query_pt;

	int GetNumRows(char *pDataCSVFile);
	int GetNumCols(char *pDataCSVFile, int nColsToSkip);
	void PopulateTable(char *pDataCSVFile, ANNpointArray pTableData, int nColsToSkip, int nRows, int nCols);
	void PrintLookupTable(char *pOutFile);
	void PrintIndexAndDistVectors(char *pOutFile);

public: 

	ANNTree(char *pDataCSVFile);
	ANNTree(char *pDataCSVFile, int nColsToSkip);

	bool Search();

	int GetNumRows() { return(nRows); }
	int GetNumCols() { return(nCols); }

	ANNdistArray GetDist() { return(dist); }
	ANNidxArray GetIndices() { return(nn_idx); }
	ANNpoint GetQueryPt() { return(query_pt); }

	ANNkd_tree *GetTree() { return(theTree); }

	~ANNTree();

};



#endif // __CLSANNTREE_H__


