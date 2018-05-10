//
// GDM2RVC.cpp
//
#include "stdafx.h"
#include "GdmBinLib.h"
#include "myCallback.h"
#include "Message.h"
#include "ParamsW16.h"
#include "GDMBufferSizeDefs.h"
#include "clsDoPath.h"
#include "clsBinaryFileClass.h"
#include "clsEsriBinaryHeader.h"
#include "GdmTransform.h"
#include "GDM_PredictorTypes.h"
#include "clsBinaryFileStack.h"
#include "clsBinaryFileBlockStack.h"
#include "ann.h" 

bool ExtractDataFromRVCTable(char *TablePath, int nRows, int nPreds, int *pClusterID, double **ppData, FPTR fptr);
int GetNumRVCRecords(char *ClassPath);
int GetRVCPredCount(char *ClassPath);


bool CreateRVCSimilarityTable(char *pParams, FPTR fptr)
{
	//Message("At entry to CreateRVCSimilarityTable", "INFO");
	int nCurrent = 0;

	//
	// setup the GDM transform grid stack
	//
	BinaryFileStack *bfs0 = new BinaryFileStack(pParams, BFS_TranGrids);
	if (NULL == bfs0)
	{
		Message("Cannot construct BinaryFileStack bfs0", "ERROR");
		return(false);
	}
	int nPreds = bfs0->GetNumGrids();


	//
	// setup the RVC Class grid
	//
	gmpath gmPath;
	char InPath[512];
	GetProfileString( "Input", "RVCClass", InPath, pParams );
	BinaryFileClass *bfc_RVC = new BinaryFileClass(gmPath.ChangeExtension(InPath, ".flt"));	
	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(InPath, ".hdr"));	
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	float fNoData = header->GetNoDataValue();
	float *pRVCRow = new float [nCols];

	//
	// get the minimum and maximum class values
	//
	int nMinClass = 100000;
	int nMaxClass = -1;
	fptr("Status: Creating RVC Lookup..", nCurrent);
	for ( int i=0; i<nRows; i++ )
	{
		if ( i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			fptr("Status: Creating RVC Lookup...", nCurrent);
		}

		// read the current scanline
		bfc_RVC->ReadFloat(pRVCRow, nCols);

		for ( int j=0; j<nCols; j++ )
		{
			if (pRVCRow[j] != fNoData)
			{
				if (pRVCRow[j] < nMinClass) nMinClass = (int)pRVCRow[j];
				if (pRVCRow[j] > nMaxClass) nMaxClass = (int)pRVCRow[j];
			}
		}
	}

	if (nMinClass < 0 )
	{
		Message("Minimum Class Valie is Less Than Zero!!!", "ERROR");
		return(false);
	}
	
	//Message(nMinClass, "nMinClass");
	//Message(nMaxClass, "nMaxClass");

	//
	// create the data matrix to hold running predictor means and class cell counts
	//
	int nClasses = nMaxClass + 1;  // allow for a zero valued class (0 = cleared in the case of RVCs)
	int *pClassCounts = new int [nClasses];
	for ( int i=0; i<nClasses; i++ ) pClassCounts[i] = 0;
	float **ppPredVal = new float * [nClasses];
	for ( int i=0; i<nClasses; i++ )
	{
		ppPredVal[i] = new float [nPreds];
		for ( int j=0; j<nPreds; j++ ) ppPredVal[i][j] = 0.0F;
	}


	//
	// re-read the RVC grid to get the cell counts in each class
	//
	bfc_RVC->SeekToStart();
	nCurrent = 0;
	fptr("Status: Counting cells in RVCs...", nCurrent);
	for ( int i=0; i<nRows; i++ )
	{
		if ( i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			fptr("Status: Counting cells in RVCs...", nCurrent);
		}

		// read the current scanline
		bfc_RVC->ReadFloat(pRVCRow, nCols);

		for ( int j=0; j<nCols; j++ )
		{
			if (pRVCRow[j] != fNoData)
			{
				pClassCounts[(int)pRVCRow[j]] += 1;
			}
		}
	}


	//
	// Extract the k-means centroid for each RVC 
	//
	bfc_RVC->SeekToStart();
	nCurrent = 0;
	fptr("Status: Extracting k-means...", nCurrent);
	for ( int i=0; i<nRows; i++ )
	{
		if ( i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			fptr("Status: Extracting k-means...", nCurrent);
		}

		bfc_RVC->ReadFloat(pRVCRow, nCols);
		bfs0->GetNextRow();
		float **ppData = bfs0->GetRowData();
		
		for (int j=0; j<nCols; j++ )
		{			
			if (pRVCRow[j] != bfs0->GetNoDataVal())
			{
				int nRVC = (int)pRVCRow[j];
				if (nRVC > 0 )
				{
					for ( int k=0; k<nPreds; k++ )
					{
						ppPredVal[nRVC][k] += ppData[k][j];
					}
				}
			}
		}
	}


	//
	// Take the mean values of all the class predictor values
	//
	for ( int i=0; i<nClasses; i++ )
	{
		if (pClassCounts[i] > 0)
		{
			for (int j=0; j<nPreds; j++ )
			{
				ppPredVal[i][j] = ppPredVal[i][j] / pClassCounts[i];
			}
		}
	}


	//
	// Write results to file
	//
	FILE *fpCount = fopen("D:\\GDM_CreateRVCs\\Class_Lookup.csv", "w+t");
	fprintf(fpCount, "Class,Count," );
	for ( int i=0; i<nPreds; i++ )
	{
		fprintf(fpCount, "Pred_%d", i+1);
		if (i < nPreds-1)
			fprintf(fpCount, ",");
		else
			fprintf(fpCount, "\n");
	}

	// we only write records that are greater than zero AND have a cell count > 0
	for ( int i=1; i<nClasses; i++ )
	{		
		if (pClassCounts[i] > 0)
		{
			fprintf(fpCount, "%d,%d,", i, pClassCounts[i] );
			for ( int j=0; j<nPreds; j++ )
			{
				fprintf(fpCount, "%f", ppPredVal[i][j]);
				if (j < nPreds-1)
					fprintf(fpCount, ",");
				else
					fprintf(fpCount, "\n");
			}
		}
	}
	fclose(fpCount);


	//
	// create a similarity table
	//
	float **ppSimilarity = new float * [nClasses];
	for ( int i=0; i<nClasses; i++ ) 
	{
		ppSimilarity[i] = new float [nClasses];
		for ( int j=0; j<nClasses-1; j++ ) ppSimilarity[i][j] = 0.0F;

		ppSimilarity[i][i] = 1.0F; // initialise diagonal
	}


	//
	// Calculate RVC Class Similarities
	//
	for ( int i=1; i<nClasses-1; i++ )
	{
		for ( int j=i+1; j<nClasses; j++ )
		{
			if ((pClassCounts[i] > 0) && (pClassCounts[j] > 0))
			{
				float fVal = 0.0F;
				for ( int k=0; k<nPreds; k++ )
				{
					fVal += fabs(ppPredVal[i][k] - ppPredVal[j][k]);
				}

				ppSimilarity[i][j] = ppSimilarity[j][i] = exp(-fVal);
			}
		}
	}


	//
	// write the similarity table to file
	//
	FILE *fpSimilarity = fopen("D:\\GDM_CreateRVCs\\similarity.csv", "w+t");	
	for ( int i=0; i<nClasses; i++ )
	{
		for ( int j=0; j<nClasses; j++ )
		{
			fprintf(fpSimilarity, "%f", ppSimilarity[i][j]);
			if (j < nClasses-1)
				fprintf(fpSimilarity, ",");
			else
				fprintf(fpSimilarity, "\n");
		}
	}
	fclose(fpSimilarity);


	//
	// clean up
	//
	for ( int i=0; i<nClasses; i++ )  if (ppSimilarity[i]) delete[] ppSimilarity[i];
	if (ppSimilarity) delete[] ppSimilarity;
	if (pClassCounts) delete[] pClassCounts;
	for ( int i=0; i<nClasses; i++ )  if (ppPredVal[i]) delete[] ppPredVal[i];
	if (ppPredVal) delete[] ppPredVal;
	bfc_RVC->Close();
	if (header) delete header;
	if (pRVCRow) delete[] pRVCRow;
	if (bfs0) delete bfs0;
	nCurrent = 0;
	fptr("Status: Ready", nCurrent);
	return(true);
}



bool DoRVClassification(char *pParams, FPTR fptr)
{
	//Message("At entry to DoRVClassification", "INFO");
	int nCurrent = 0;

	//
	// setup the GDM transform grid stack
	//
	BinaryFileStack *bfs0 = new BinaryFileStack(pParams, BFS_TranGrids);
	if (NULL == bfs0)
	{
		Message("Cannot construct BinaryFileStack bfs0", "ERROR");
		return(false);
	}
	int nPreds = bfs0->GetNumGrids();


	//
	// setup the RVC Class grid
	//
	gmpath gmPath;
	char InPath[512];
	GetProfileString( "Input", "RVCClass", InPath, pParams );
	BinaryFileClass *bfc_RVC = new BinaryFileClass(gmPath.ChangeExtension(InPath, ".flt"));	
	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(InPath, ".hdr"));	
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	float fNoData = header->GetNoDataValue();
	float *pRVCRow = new float [nCols];


	//
	// setup the output classification grid
	//
	char OutPath[512];
	GetProfileString( "Output", "OutClass", OutPath, pParams );
	BinaryFileClass *bfc_Out = new BinaryFileClass(gmPath.ChangeExtension(OutPath, ".flt"), BFC_CreateOrTrunc);	
	float *pOutRow = new float [nCols];


	//
	// setup the ANN lookup table
	//
	char ClassPath[512];
	GetProfileString( "Input", "ClassLookup", ClassPath, pParams );
	int nRecords = GetNumRVCRecords(ClassPath);
	int nTablePreds = GetRVCPredCount(ClassPath);

	//Message(nRecords, "nRecords");
	//Message(nTablePreds, "nTablePreds");

	if (nTablePreds != nPreds)
	{
		Message("Got different predictor count in lookup and grid stack!!!", "ERROR");
		return(false);
	}


	//
	// extract classes and predictor data
	//
	int *pClusterID = new int [nRecords];
	double **ppData = new double * [nRecords];
	for (int i=0; i<nRecords; i++ )
	{
		ppData[i] = new double [nTablePreds];
	}
	ExtractDataFromRVCTable(ClassPath, nRecords, nTablePreds, pClusterID, ppData, fptr);

	
	//
	// construct the ANN tree and populate from ppData
	//
	ANNpointArray data_pnts = annAllocPts( nRecords, nTablePreds );
	for ( int i=0; i<nRecords; i++ )
	{
		for ( int j=0; j<nTablePreds; j++ )
		{
			data_pnts[i][j] = ppData[i][j];
		}
	}


	//
	// now build the kd tree search structure with 
	// the data points, number of data points and dimension of search space
	//
	ANNkd_tree *theTree = new ANNkd_tree ( data_pnts, nRecords, nTablePreds ); 


	//
	// find the nearest neighbour classes
	//
	// 
	// setup some structures for the queries to follow
	ANNidxArray nn_idx = new ANNidx[1];
	ANNdistArray dist = new ANNdist[1];
	ANNpoint query_pt = annAllocPt( nTablePreds );
		

	//
	// Do the classification
	//
	fptr("Status: Creating GDM To RVC Classification...", nCurrent);
	for ( int i=0; i<nRows; i++ )
	{
		if ( i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			fptr("Status: Creating GDM To RVC Classification...", nCurrent);
		}

		// read the current scanline
		bfc_RVC->ReadFloat(pRVCRow, nCols);
		bfs0->GetNextRow();
		float **ppPredData = bfs0->GetRowData();

		for ( int j=0; j<nCols; j++ )
		{
			if (pRVCRow[j] != fNoData)
			{
				if (pRVCRow[j] == 0) // classify cleared classes to nearest uncleared class
				{
					for ( int k=0; k<nPreds; k++ )
					{
						query_pt[k] = ppPredData[k][j];
					}

					// find the nearest neighbour class in the ANN tree
					theTree->annkSearch( query_pt, 1, nn_idx, dist, 0 );		
	
					// set the class in the output grid
					pOutRow[j] = (float)pClusterID[ nn_idx[0] ];
				}
				else  // ignore non-cleared classes
				{
					pOutRow[j] = 0.0F;
				}
			}
			else
			{
				pOutRow[j] = fNoData;
			}
		}

		bfc_Out->WriteFloat(pOutRow, nCols);
	}


	//
	// clean up
	//
	if (pClusterID) delete[] pClusterID;
	for ( int i=0; i<nRecords; i++ ) if ( ppData[i] ) delete ppData[i];
	if ( ppData ) delete[] ppData;
	if (theTree) delete theTree;
	header->CopyTo(gmPath.ChangeExtension(OutPath, ".hdr"));
	bfc_Out->Close();
	bfc_RVC->Close();
	if (header) delete header;
	if (pOutRow) delete[] pOutRow;
	if (pRVCRow) delete[] pRVCRow;
	if (bfs0) delete bfs0;
	nCurrent = 0;
	fptr("Status: Ready", nCurrent);
	return(true);
}



bool ExtractDataFromRVCTable(char *TablePath, int nRows, int nPreds, int *pClusterID, double **ppData, FPTR fptr)
{
	//
	// extract from the training table
	//
	char *seps = ",\n";
	char *Buff = new char [TABLE_ROW_BUFFSIZE];
	FILE *fp = fopen(TablePath, "r+t");
	fgets(Buff, TABLE_ROW_BUFFSIZE, fp); // get header
	int nCurrent = 0;
	for ( int i=0; i<nRows; i++ )
	{
		if (i * 100 / nRows > nCurrent)
		{
			fptr("Extracting Table Data...", nCurrent);	
		}

		fgets(Buff, TABLE_ROW_BUFFSIZE, fp);

		// get Class_id
		char *p = strtok(Buff, seps);
		pClusterID[i] = atoi(p);	  // set class
		
		// count
		p = strtok(NULL, seps);       // get count

		// preds...
		for ( int j=0; j<nPreds; j++ )
		{
			p = strtok(NULL, seps);
			ppData[i][j] = atof(p);
		}
	}
	fclose(fp);
	if (Buff) delete[] Buff;
	fptr("Extracting Table Data...", 0);	
	return(true);
}


//
// Count number of RVC Record in lookup file assuming it also has a header line
// 
int GetNumRVCRecords(char *ClassPath)
{
	int nCount = 0;
	char *Buff = new char [TABLE_ROW_BUFFSIZE];
	FILE *fp = fopen(ClassPath, "r+t");
	fgets(Buff, TABLE_ROW_BUFFSIZE, fp); // get header
	while(1)
	{
		if (NULL == fgets(Buff, TABLE_ROW_BUFFSIZE, fp))
			break;

		++nCount;
	}
	fclose(fp);
	return(nCount);
}


//
// Count number of RVC Preds in lookup file assuming it also has a class and countcolumn and a header line
// 
int GetRVCPredCount(char *ClassPath)
{	
	char *Buff = new char [TABLE_ROW_BUFFSIZE];
	FILE *fp = fopen(ClassPath, "r+t");
	fgets(Buff, TABLE_ROW_BUFFSIZE, fp); // get header
	fclose(fp);

	// skip first two cols
	char *seps = ",\n";
	
	// get Class_id
	char *p = strtok(Buff, seps);
	
	// get count
	p = strtok(NULL, seps);      

	// now count predictor columns
	int nPreds = 0;
	p = strtok(NULL, seps);      
	while(p)
	{
		++nPreds;
		p = strtok(NULL, seps);      
	}
	return(nPreds);
}

