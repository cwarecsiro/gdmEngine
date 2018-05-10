//
// GDMClassification.cpp
//
#include "stdafx.h"
#include "GdmBinLib.h"
#include "Message.h"
#include "GDMBufferSizeDefs.h"
#include "ParamsW16.h"
#include "FilterSites.h"
#include "GdmExtractQuants.h"
#include "clsDoPath.h"
#include "GDMClassification.h"
#include "Cluster.h"
#include "pcoMatrix.h"
#include "clsBinaryFileClass.h"
#include "clsEsriBinaryHeader.h"
#include "ESRIBinaryGridClass.h"
#include "GdmTransform.h"
#include "GdmUnconstrainedProbGrids.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// for the ANN Library
#include "ann.h" 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <io.h>
#include <fcntl.h>
#include <sys/stat.h>

#include <iostream>
#include <fstream>

static double EPSILON = 0.0000001;

//
// write a csv table _ID, X, Y, Val=0 with XYs being the centroid of each gridcell in the user defined extent
//
bool CreateRasterCentroids(double dCellSize, double dLeftEdge, double dTopEdge, int nNumRows, int nNumCols, char *pOutPath, FPTR fptr)
{
	double dHalfCell = dCellSize / 2.0;
	int nID = 0;

	FILE *fp = fopen(pOutPath, "w+t");
	fprintf(fp, "_ID,X,Y,Value\n");
	int nCurrent = 0;
	fptr("Creating Centroid Table...", nCurrent);
	for (int i = 0; i < nNumRows; i++)
	{
		if (i * 100 / nNumRows > nCurrent)
		{
			nCurrent = i * 100 / nNumRows;
			fptr("Creating Centroid Table...", nCurrent);
		}

		for (int j = 0; j < nNumCols; j++)
		{
			++nID;
			
			double xVal = (j * dCellSize) + dLeftEdge + dHalfCell;
			double yVal = dTopEdge - (i * dCellSize) - dHalfCell;

			fprintf(fp, "%d,%lf,%lf,%d\n", nID, xVal, yVal, 0);
		}

	}
	fptr("Ready...", 0);
	fclose(fp);
	return(true);
}





//
// write a csv table derived from the above for use in R for creating Dendrograms
// It will contain a record for each class, followed by the predictor data that is the 
// average of each predictor in the respective class in the table created in DumpANNFile()
//
// Note: values in pClusterID are ONE-BASED so subtract 1 if using for array indices
//
void WriteDendrogramTable(char *lpParams, int *pClusterID, char *lpOutPath, double **ppData, int nRows, int nCols)
{
	gmpath gmPath;

	FILE *fp = fopen(lpOutPath, "w+t");
	if (NULL == fp)
	{
		Message("Cannot create Class K-Means Path", "WriteDendrogramTable");
		return;
	}

	// write header
	if (1 == GetProfileInt("GDMODEL", "UseEuclidean", lpParams))
	{
		fprintf(fp, "Class,Count,XTran,YTran");
	}
	else
	{
		fprintf(fp, "Class,Count");
	}

	// extract the filenames fror the column header
	int nPreds = GetProfileInt("PREDICTORS", "NumPredictors", lpParams);
	char lpKey[64];
	char lpPredPath[BUFFLEN];
	for (int i = 1; i <= nPreds; i++)
	{
		sprintf(lpKey, "PredTran%d", i);
		GetProfileString("TRANSPREDS", lpKey, lpPredPath, lpParams);

		if (strlen(lpPredPath) > 0)
		{
			sprintf(lpPredPath, "%s", gmPath.ChangeExtension(lpPredPath, ".flt"));
			fprintf(fp, ",%s", gmPath.GetName(lpPredPath));
		}

		if (i == nPreds)
			fprintf(fp, "\n");
	}

	// get maximum class number
	int nMaxClass = 0;
	for (int i = 0; i < nRows; i++)
	{
		if (nMaxClass < pClusterID[i]) nMaxClass = pClusterID[i];
	}

	// allocate the count and running average vectors and initialise to zero
	int *pCounts = new int[nMaxClass];
	for (int i = 0; i < nMaxClass; i++) pCounts[i] = 0;
	double **ppValues = new double *[nRows];
	for (int i = 0; i < nRows; i++)
	{
		ppValues[i] = new double[nCols];
		for (int j = 0; j < nCols; j++) ppValues[i][j] = 0.0;
	}

	// extract data
	for (int i = 0; i < nRows; i++)
	{
		pCounts[pClusterID[i] - 1] += 1;

		for (int j = 0; j < nCols; j++)
		{
			ppValues[pClusterID[i] - 1][j] += ppData[i][j];
		}
	}

	// now write the K-Means of the transformed predictor fields
	for (int i = 0; i < nMaxClass; i++)
	{
		fprintf(fp, "%d,%d,", i + 1, pCounts[i]);
		for (int j = 0; j < nCols; j++)
		{
			fprintf(fp, "%lf", ppValues[i][j] / pCounts[i]);

			if (j < (nCols - 1))
				fprintf(fp, ",");
			else
				fprintf(fp, "\n");
		}
	}
	fflush(fp);
	if (fp) fclose(fp);
	if (pCounts) delete[] pCounts;
	for (int i = 0; i < nMaxClass; i++)
	{
		if (ppValues[i]) delete[] ppValues[i];
	}
	if (ppValues) delete[] ppValues;
}



//
// Create a csv table with the following fields: Class,Env_Distance,X,Y representing the closest 
// Manhattan distance between the K-Means centroid for each class and its' GDM classification cells
// Also include the GDM Transformed values for the selected sites.
//
bool CreateKMeansClosestPointTable(char *pParams,
								   char *lpClassPath,
	                               char *lpInTablePath,
	                               char *lpOutTablePath,
	                               int nClasses,
	                               FPTR fptr)
{
	int nCurrent = 0;
	gmpath gmPath;
	
	// setup the GDM Classification grid...
	ESRIBinaryGridClass *ebgClass = new ESRIBinaryGridClass(gmPath.ChangeExtension(lpClassPath, ".flt"));
	int nRows = ebgClass->GetMetrics()->GetNumRows();
	int nCols = ebgClass->GetMetrics()->GetNumCols();
	float fNoData = ebgClass->GetMetrics()->GetNoDataValue();

	// setup class counts and initialise to zero...
	int *pClassCounts = new int[nClasses];
	for (int i = 0; i < nClasses; i++) pClassCounts[i] = 0;

	// setup class index extents and init to index of opposite axes 
	int *pClassXMin = new int[nClasses];
	int *pClassXMax = new int[nClasses];
	int *pClassYMin = new int[nClasses];
	int *pClassYMax = new int[nClasses];
	for (int i = 0; i < nClasses; i++)
	{
		pClassXMin[i] = nCols;
		pClassXMax[i] = -1;
		pClassYMin[i] = nRows;
		pClassYMax[i] = -1;
	}


	//
	// count the gridcells in each of the GDM classes...
	// also determine the zero based indices of the class extents
	//
	for (int i = 0; i < ebgClass->GetMetrics()->GetNumRows(); i++)
	{
		ebgClass->ReadRow();
		for (int j = 0; j < nCols; j++)
		{
			float fVal = ebgClass->GetRowData()[j];
			if (fVal != fNoData)
			{
				// convert ONE-Based index to ZERO based index...
				int nIndex = ((int)fVal) - 1;
				
				//adjust count
				pClassCounts[nIndex] += 1;

				// adjust indices
				if (pClassXMin[nIndex] > j) pClassXMin[nIndex] = j;
				if (pClassXMax[nIndex] < j) pClassXMax[nIndex] = j;
				if (pClassYMin[nIndex] > i) pClassYMin[nIndex] = i;
				if (pClassYMax[nIndex] < i) pClassYMax[nIndex] = i;
			}
		}
	}


	//
	// cycle through each class, create an ANN array and extract the nearest site from the KMeans centroid table
	//
	int nTranPreds = GetNumGDMTransformedPreds(pParams);
	char **ppTranGridPaths = GetTransformPaths(pParams);
	ESRIBinaryGridClass **ebgTrans = new ESRIBinaryGridClass * [nTranPreds];
	for (int i = 0; i < nTranPreds; i++)
	{
		ebgTrans[i] = new ESRIBinaryGridClass(gmPath.ChangeExtension(ppTranGridPaths[i], ".flt"));
	}

	int nNearest = 1;
	ANNpoint query_pt = annAllocPt(nTranPreds);
	ANNidxArray nn_idx = new ANNidx[nNearest];
	ANNdistArray dist = new ANNdist[nNearest];

	char buffer[BUFF1024];
	FILE *fp = fopen(lpInTablePath, "r+t");
	fgets(buffer, BUFF1024, fp); // skip header

	FILE *fpOut = fopen(lpOutTablePath, "w+t");
	fprintf(fpOut, "Class,X,Y,Distance");
	for (int tp = 0; tp < nTranPreds; tp++)
	{
		fprintf(fpOut, ",%s", gmPath.GetName(ppTranGridPaths[tp]));
	}
	fprintf(fpOut, "\n");

	nCurrent = 0;
	for (int i = 0; i < nClasses; i++)
	{
		if (i * 100 / nClasses > nCurrent)
		{
			nCurrent = i * 100 / nClasses;
			fptr("Creating K-Means Closest Point Table...", nCurrent);
		}

		//
		// extract the transform values for sites in the current class
		//
		ANNpointArray data_pnts = annAllocPts(pClassCounts[i], nTranPreds);
		double *pXLoc = new double[pClassCounts[i]];
		double *pYLoc = new double[pClassCounts[i]];

		// position to first data row for current class
		ebgClass->GetBinary()->SeekTo(pClassYMin[i] * nCols * sizeof(float)); 
		for (int tp = 0; tp < nTranPreds; tp++)
		{
			ebgTrans[tp]->GetBinary()->SeekTo(pClassYMin[i] * nCols * sizeof(float));
		}

		int nThisRow = 0;
		for (int y = pClassYMin[i]; y <= pClassYMax[i]; y++)
		{
			ebgClass->ReadRow();
			for (int tp = 0; tp < nTranPreds; tp++)
			{
				ebgTrans[tp]->ReadRow();
			}

			for (int x = pClassXMin[i]; x <= pClassXMax[i]; x++)
			{
				float fVal = ebgClass->GetRowData()[x];
				if (fVal != fNoData)
				{
					if (fVal == (i + 1))  // class indices are stored as one-based in the classification grid
					{
						// extract transform values for this cell
						for (int tp = 0; tp < nTranPreds; tp++)
						{
							data_pnts[nThisRow][tp] = ebgTrans[tp]->GetRowData()[x];
						}

						pXLoc[nThisRow] = ebgClass->GetMetrics()->GetMinX() + 
							              (x * ebgClass->GetMetrics()->GetCellSize()) +
							              (ebgClass->GetMetrics()->GetCellSize() / 2);

						pYLoc[nThisRow] = ebgClass->GetMetrics()->GetMaxY() -
							              (y * ebgClass->GetMetrics()->GetCellSize()) -
							              (ebgClass->GetMetrics()->GetCellSize() / 2);

						++nThisRow;
					}
				}
			} 
		} 
		
		//
		// now build the kd tree search structure with 
		// the data points, number of data points and dimension of search space
		//
		ANNkd_tree *theTree = new ANNkd_tree(data_pnts, pClassCounts[i], nTranPreds);

		//
		// Initialise the query vector from the K-Means table file
		//
		fgets(buffer, BUFF1024, fp);
		char *p = strtok(buffer, ",\n"); // get the class field
		p = strtok(NULL, ",\n");         // get the count field
		// now get the transform values
		for (int tp = 0; tp < nTranPreds; tp++)
		{
			p = strtok(NULL, ",\n");
			query_pt[tp] = atof(p);
		}

		//
		// Search for the nearest neighbours in Transformed GDM space...
		//
		theTree->annkSearch(query_pt, nNearest, nn_idx, dist, 0);

		//
		// Write the best closest site X/Y's for this class
		//
		for (int d = 0; d < nNearest; d++)
		{
			fprintf(fpOut, "%d,%lf,%lf,%lf", i + 1, pXLoc[nn_idx[d]], pYLoc[nn_idx[d]], dist[d]);
			for (int tp = 0; tp < nTranPreds; tp++)
			{
				fprintf(fpOut, ",%lf", data_pnts[nn_idx[d]][tp]);
			}
			fprintf(fpOut, "\n");
		}

		// cleanup
		if (theTree) delete theTree;
		if (pXLoc) delete[] pXLoc;
		if (pYLoc) delete[] pYLoc;		
	} // for (int i = 0; i < nClasses; i++)

	if (fp) fclose(fp);
	if (fpOut) fclose(fpOut);
	if (ebgClass) delete ebgClass;
	for (int i = 0; i < nTranPreds; i++)
	{
		if (ebgTrans[i]) delete ebgTrans[i];
	}
	if (ebgTrans) delete ebgTrans;
	if (pClassCounts) delete[] pClassCounts;
	if (pClassXMin) delete[] pClassXMin;
	if (pClassXMax) delete[] pClassXMax;
	if (pClassYMin) delete[] pClassYMin;
	if (pClassYMax) delete[] pClassYMax;
	for (int i = 0; i < nTranPreds; i++)
	{
		if (ppTranGridPaths[i]) delete[] ppTranGridPaths[i];
	}
	if (ppTranGridPaths) delete[] ppTranGridPaths;
	if (nn_idx) delete[] nn_idx;
	if (dist) delete[] dist;
	fptr("Ready:", 0);
	return(true);
}


//
// Create a csv table with the following fields: Class,X,Y representing the closest 
// Manhattan distance between the K-Means centroid for each class and its' GDM classification cells
// Also include the GDM Transformed values for the selected sites.
//
bool CreateKMeansDemandPointTable(char *pParams,
	                              char *lpClassPath,
	                              char *lpInTablePath,
	                              char *lpOutTablePath,
	                              int nClasses,
	                              FPTR fptr)
{
	int nCurrent = 0;
	gmpath gmPath;

	// setup the GDM Classification grid...
	ESRIBinaryGridClass *ebgClass = new ESRIBinaryGridClass(gmPath.ChangeExtension(lpClassPath, ".flt"));
	int nRows = ebgClass->GetMetrics()->GetNumRows();
	int nCols = ebgClass->GetMetrics()->GetNumCols();
	float fNoData = ebgClass->GetMetrics()->GetNoDataValue();

	// setup class counts and initialise to zero...
	int *pClassCounts = new int[nClasses];
	for (int i = 0; i < nClasses; i++) pClassCounts[i] = 0;

	// setup class index extents and init to index of opposite axes 
	int *pClassXMin = new int[nClasses];
	int *pClassXMax = new int[nClasses];
	int *pClassYMin = new int[nClasses];
	int *pClassYMax = new int[nClasses];
	for (int i = 0; i < nClasses; i++)
	{
		pClassXMin[i] = nCols;
		pClassXMax[i] = -1;
		pClassYMin[i] = nRows;
		pClassYMax[i] = -1;
	}


	//
	// count the gridcells in each of the GDM classes...
	// also determine the zero based indices of the class extents
	//
	for (int i = 0; i < ebgClass->GetMetrics()->GetNumRows(); i++)
	{
		ebgClass->ReadRow();
		for (int j = 0; j < nCols; j++)
		{
			float fVal = ebgClass->GetRowData()[j];
			if (fVal != fNoData)
			{
				// convert ONE-Based index to ZERO based index...
				int nIndex = ((int)fVal) - 1;

				//adjust count
				pClassCounts[nIndex] += 1;

				// adjust indices
				if (pClassXMin[nIndex] > j) pClassXMin[nIndex] = j;
				if (pClassXMax[nIndex] < j) pClassXMax[nIndex] = j;
				if (pClassYMin[nIndex] > i) pClassYMin[nIndex] = i;
				if (pClassYMax[nIndex] < i) pClassYMax[nIndex] = i;
			}
		}
	}

	//
	// cycle through each class, create an ANN array and extract the nearest site from the KMeans centroid table
	//
	int nTranPreds = GetNumGDMTransformedPreds(pParams);
	char **ppTranGridPaths = GetTransformPaths(pParams);
	ESRIBinaryGridClass **ebgTrans = new ESRIBinaryGridClass *[nTranPreds];
	for (int i = 0; i < nTranPreds; i++)
	{
		ebgTrans[i] = new ESRIBinaryGridClass(gmPath.ChangeExtension(ppTranGridPaths[i], ".flt"));
	}

	int nNearest = 1;
	ANNpoint query_pt = annAllocPt(nTranPreds);
	ANNidxArray nn_idx = new ANNidx[nNearest];
	ANNdistArray dist = new ANNdist[nNearest];

	char buffer[BUFF1024];
	FILE *fp = fopen(lpInTablePath, "r+t");
	fgets(buffer, BUFF1024, fp); // skip header

	FILE *fpOut = fopen(lpOutTablePath, "w+t");
	fprintf(fpOut, "Class,X,Y");
	for (int tp = 0; tp < nTranPreds; tp++)
	{
		fprintf(fpOut, ",%s", gmPath.GetName(ppTranGridPaths[tp]));
	}
	fprintf(fpOut, "\n");

	nCurrent = 0;
	for (int i = 0; i < nClasses; i++)
	{
		if (i * 100 / nClasses > nCurrent)
		{
			nCurrent = i * 100 / nClasses;
			fptr("Creating K-Means DemandPointTable...", nCurrent);
		}

		//
		// extract the transform values for sites in the current class
		//
		ANNpointArray data_pnts = annAllocPts(pClassCounts[i], nTranPreds);
		double *pXLoc = new double[pClassCounts[i]];
		double *pYLoc = new double[pClassCounts[i]];

		// position to first data row for current class
		ebgClass->GetBinary()->SeekTo(pClassYMin[i] * nCols * sizeof(float));
		for (int tp = 0; tp < nTranPreds; tp++)
		{
			ebgTrans[tp]->GetBinary()->SeekTo(pClassYMin[i] * nCols * sizeof(float));
		}

		int nThisRow = 0;
		for (int y = pClassYMin[i]; y <= pClassYMax[i]; y++)
		{
			ebgClass->ReadRow();
			for (int tp = 0; tp < nTranPreds; tp++)
			{
				ebgTrans[tp]->ReadRow();
			}

			for (int x = pClassXMin[i]; x <= pClassXMax[i]; x++)
			{
				float fVal = ebgClass->GetRowData()[x];
				if (fVal != fNoData)
				{
					if (fVal == (i + 1))  // class indices are stored as one-based in the classification grid
					{
						// extract transform values for this cell
						for (int tp = 0; tp < nTranPreds; tp++)
						{
							data_pnts[nThisRow][tp] = ebgTrans[tp]->GetRowData()[x];
						}

						pXLoc[nThisRow] = ebgClass->GetMetrics()->GetMinX() +
							(x * ebgClass->GetMetrics()->GetCellSize()) +
							(ebgClass->GetMetrics()->GetCellSize() / 2);

						pYLoc[nThisRow] = ebgClass->GetMetrics()->GetMaxY() -
							(y * ebgClass->GetMetrics()->GetCellSize()) -
							(ebgClass->GetMetrics()->GetCellSize() / 2);

						++nThisRow;
					}
				}
			}
		}

		//
		// now build the kd tree search structure with 
		// the data points, number of data points and dimension of search space
		//
		ANNkd_tree *theTree = new ANNkd_tree(data_pnts, pClassCounts[i], nTranPreds);

		//
		// Initialise the query vector from the K-Means table file
		//
		fgets(buffer, BUFF1024, fp);
		char *p = strtok(buffer, ",\n"); // get the class field
		p = strtok(NULL, ",\n");         // get the count field
										 // now get the transform values
		for (int tp = 0; tp < nTranPreds; tp++)
		{
			p = strtok(NULL, ",\n");
			query_pt[tp] = atof(p);
		}

		//
		// Search for the nearest neighbours in Transformed GDM space...
		//
		theTree->annkSearch(query_pt, nNearest, nn_idx, dist, 0);

		//
		// Write the best closest site X/Y's for this class
		//
		for (int d = 0; d < nNearest; d++)
		{
			fprintf(fpOut, "%d,%lf,%lf", i + 1, pXLoc[nn_idx[d]], pYLoc[nn_idx[d]]);
			for (int tp = 0; tp < nTranPreds; tp++)
			{
				fprintf(fpOut, ",%lf", data_pnts[nn_idx[d]][tp]);
			}
			fprintf(fpOut, "\n");
		}

		// cleanup
		if (theTree) delete theTree;
		if (pXLoc) delete[] pXLoc;
		if (pYLoc) delete[] pYLoc;
	} // for (int i = 0; i < nClasses; i++)

	if (fp) fclose(fp);
	if (fpOut) fclose(fpOut);
	if (ebgClass) delete ebgClass;
	for (int i = 0; i < nTranPreds; i++)
	{
		if (ebgTrans[i]) delete ebgTrans[i];
	}
	if (ebgTrans) delete ebgTrans;
	if (pClassCounts) delete[] pClassCounts;
	if (pClassXMin) delete[] pClassXMin;
	if (pClassXMax) delete[] pClassXMax;
	if (pClassYMin) delete[] pClassYMin;
	if (pClassYMax) delete[] pClassYMax;
	for (int i = 0; i < nTranPreds; i++)
	{
		if (ppTranGridPaths[i]) delete[] ppTranGridPaths[i];
	}
	if (ppTranGridPaths) delete[] ppTranGridPaths;
	if (nn_idx) delete[] nn_idx;
	if (dist) delete[] dist;
	fptr("Ready:", 0);
	return(true);
}



//
// Get number of transformed predictors from a GDM model
//
int GetNumGDMTransformedPreds(char *pParams)
{
	//
	// are we using euclidean ?
	//
	bool fUseEuclidean = (1 == GetProfileInt("GDMODEL", "UseEuclidean", pParams)) ? true : false;

	//
	// determine the number of columns ( the number of transform grids )
	//
	int nCols = 0;
	if (fUseEuclidean) nCols += 2;

	// total number of possible preds, need to check for presence to get file transform count
	int nPreds = GetProfileInt("PREDICTORS", "NumPredictors", pParams);
	char lpKey[64];
	char lpPredPath[BUFFLEN];
	for (int i = 1; i <= nPreds; i++)
	{
		sprintf(lpKey, "PredTran%d", i);
		GetProfileString("TRANSPREDS", lpKey, lpPredPath, pParams);

		if (strlen(lpPredPath) > 0)
		{
			++nCols;
		}
	}
	return(nCols);
}


//
// Get an array of GDM Transform Grid Paths from the parameter file
//
char **GetTransformPaths(char *pParams)
{
	int nTransforms = GetNumGDMTransformedPreds(pParams);
	char **ppPaths = new char *[nTransforms];
	for (int i = 0; i < nTransforms; i++)
	{
		ppPaths[i] = new char[BUFFLEN];
	}

	//
	// are we using euclidean ?
	//
	bool fUseEuclidean = (1 == GetProfileInt("GDMODEL", "UseEuclidean", pParams)) ? true : false;

	// total number of possible preds, need to check for presence to get file transform count
	int nPreds = GetProfileInt("PREDICTORS", "NumPredictors", pParams);
	char lpKey[64];
	char lpPredPath[BUFFLEN];
	int nThisPath = 0;
	for (int i = 1; i <= nPreds; i++)
	{
		sprintf(lpKey, "PredTran%d", i);
		GetProfileString("TRANSPREDS", lpKey, lpPredPath, pParams);

		if (strlen(lpPredPath) > 0)
		{
			strcpy(ppPaths[nThisPath], lpPredPath);
			++nThisPath;
		}
	}

	return(ppPaths);
}



//
// Perform an unsupervised classification from a composite GDM table
//
bool ClassifyFromCompositeTableData( char *ParamPath,
	                                 int nClasses, 
									 int nRows,
                                     int nPreds, 
	                                 char *InputTablePath, 
                                     char *OutputTablePath, 
                                     FPTR fptr)
{
	//Message("Entering ClassifyFromCompositeTableData");
	//return(true);
	fptr("Doing GDM Classification...", 10);	


	//
	// get the GDM model splines, quantiles and coefficients
	//
	char lpKey[64];
	int *pSplines = new int [nPreds];
	for ( int i=0; i<nPreds; i++ )
	{
		sprintf(lpKey, "PredSpl%d", i+1);
		pSplines[i] = GetProfileInt("PREDICTORS", lpKey, ParamPath);
	}

	double **ppQuantiles = new double * [nPreds];
	for ( int i=0; i<nPreds; i++ )
	{
		ppQuantiles[i] = new double [pSplines[i]];

		for ( int j=0; j<pSplines[i]; j++ )
		{
			sprintf(lpKey, "PredSplVal%d.%d", i+1, j+1); 
			ppQuantiles[i][j] = GetProfileDouble("PREDICTORS", lpKey, ParamPath);
		}
	}

	double **ppCoeffs = new double * [nPreds];
	for ( int i=0; i<nPreds; i++ )
	{
		ppCoeffs[i] = new double [pSplines[i]];

		for ( int j=0; j<pSplines[i]; j++ )
		{
			sprintf(lpKey, "PredCoef%d.%d", i+1, j+1); 
			ppCoeffs[i][j] = GetProfileDouble("PREDICTORS", lpKey, ParamPath);
		}
	}


	//
	// Allocate and populate data memory
	//
    double **ppData = new double * [nRows];
	for ( int i=0; i<nRows; i++ )
	{
		ppData[i] = new double [nPreds];
	}
	char *buff = new char [TABLE_ROW_BUFFSIZE];
	FILE *fp = fopen(InputTablePath, "r+t");
	fgets(buff, TABLE_ROW_BUFFSIZE, fp);
	for ( int i=0; i<nRows; i++ )
	{
		fgets(buff, TABLE_ROW_BUFFSIZE, fp);
		char *p = strtok(buff, ",\n");			// get class
		p = strtok(NULL, ",\n");				// get X
		p = strtok(NULL, ",\n");				// get Y
		
		for ( int j=0; j<nPreds; j++ )
		{
			p = strtok(NULL, ",\n");			// get pred
			ppData[i][j] = CalculateGDMTransform( atof(p), pSplines[j], ppQuantiles[j], ppCoeffs[j] );
		}
	}
	fclose(fp);


	//
	// Do a grid based nearest neighbour classification using the transform grids generated 
	// by the GDM process and the comma delimited table written in CreateSamplePointFile().
	//
	fptr("Doing Hierachical Clustering...", 25);	
	int *pClusterID = GetHierachicalClusterClasses( ppData, nRows, nPreds, nClasses, false, fptr );


	//
	// Write the classification table
	//
	fptr("Writing GDM Classification Table...", 75);	
	FILE *fpIn = fopen(InputTablePath, "r+t");
	fgets(buff, TABLE_ROW_BUFFSIZE, fpIn);
	
	FILE *fpOut = fopen(OutputTablePath, "w+t");
	fputs(buff, fpOut);
	for ( int i=0; i<nRows; i++ )
	{
		fgets(buff, TABLE_ROW_BUFFSIZE, fpIn);
		char *p = strtok(buff, ","); 
		
		p = strtok(NULL, ",\n");
		double dX = atof(p);

		p = strtok(NULL, ",\n");
		double dY = atof(p);

		//p = strtok(NULL, "\n");
		//fprintf(fpOut, "%d,%lf,%lf,%s\n", pClusterID[i], dX, dY, p );

		fprintf(fpOut, "%d,%lf,%lf,", pClusterID[i], dX, dY );
		for ( int j=0; j<nPreds; j++ )
		{
			fprintf(fpOut, "%lf", ppData[i][j] );
			if (j < nPreds-1)
				fprintf(fpOut, ",");
			else
				fprintf(fpOut, "\n");
		}
	}

	fclose(fpIn);
	fclose(fpOut);


	//
	// Do the classification colouring to produce a suite of recoloured classifications
	//
	fptr("Doing GDM Classification Colouring...", 90);	
	if ( false == DoClassificationColoring( NULL, OutputTablePath, nClasses, "CompositeColor", false, fptr ) )
	{
		Message( "DoClassificationColoring returned false", "ClassifyFromCompositeTableData" );
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Get the maximum class index (total number of classes)
	//
	int maxIndex = 0;
	for ( int i=0; i<nRows; i++ )
	{
		if ( pClusterID[i] > maxIndex )
		{
			maxIndex = pClusterID[i];
		}
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Allocate class counters (count the number of samples in each class)
	//
	int *ClassCount = new int [maxIndex];
	for ( int i=0; i<maxIndex; i++ )
	{
		ClassCount[i] = 0;
	}

	//
	// Initialise class counters
	//
	for ( int i=0; i<nRows; i++ )
	{
		++ClassCount[pClusterID[i]-1];
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Allocate class averages [ClassIndex][Sum then divide by count for average]
	//
	double **ppClassAverages = new double * [maxIndex];
	for (int i=0; i<maxIndex; i++ )
	{
		ppClassAverages[i] = new double [nPreds];

		// init to zero
		for ( int j=0; j<nPreds; j++ ) ppClassAverages[i][j] = 0.0;
	}


	//
	// Get the sum of the transform values for each class
	//
	for (int i=0; i<nRows; i++ )
	{
		for ( int j=0; j<nPreds; j++ )
		{
			ppClassAverages[pClusterID[i]-1][j] += ppData[i][j];
		}
	}


	//
	// Divide by the class count to get the average
	//
	for (int i=0; i<maxIndex; i++ )
	{
		for ( int j=0; j<nPreds; j++ )
		{
			if (ClassCount[i] < 1) 
				ppClassAverages[i][j] = 0.0;
			else
				ppClassAverages[i][j] /= ClassCount[i];
		}
	}


	//
	// Create an inter-class distance matrix that uses the MEAN between classes as the distance metric
	//
	double **DistMatrix_Mean = new double * [maxIndex];
	for ( int i=0; i<maxIndex; i++ )
	{
		DistMatrix_Mean[i] = new double [maxIndex];
		for ( int j=0; j<maxIndex; j++ )
		{
			DistMatrix_Mean[i][j] = 1.0;

			if (ClassCount[i] < 1) 
				DistMatrix_Mean[i][j] = 0.0;
		}
	}


	//
	// Use the GDM link function to calculate similarities
	//
	for ( int i=0; i<maxIndex; i++)
	{
		for ( int j=0; j<maxIndex; j++)
		{
			if ((ClassCount[i] < 1) || (ClassCount[j] < 1))
			{
				DistMatrix_Mean[i][j] = DistMatrix_Mean[j][i] = 0.0;
			}
			
			else if ( i != j)
			{
				double dVal = 0.0;
				for ( int k=0; k<nPreds; k++ )
				{
					dVal += fabs(ppClassAverages[i][k] - ppClassAverages[j][k]);
				}

				double dDist = exp(-(dVal));
				DistMatrix_Mean[i][j] = DistMatrix_Mean[j][i] = dDist;
			}
		}
	}


	//
	// Write distance matrix to .csv
	//
	gmpath GmPath;
	char namestring [FILEPATH_BUFFSIZE];
	sprintf(namestring, "%s\\%s_Sim_Mean.csv", GmPath.GetDirectoryPath(OutputTablePath), GmPath.GetName(OutputTablePath));
	FILE *fpSim = fopen( namestring, "w+t" );
	for ( int i=0; i<maxIndex; i++ )
	{
		for ( int j=0; j<maxIndex; j++ )
		{
			fprintf(fpSim, "%lf", DistMatrix_Mean[i][j]);
			if ( j < maxIndex-1 ) 
				fprintf(fpSim, ",");
			else
				fprintf(fpSim, "\n");
		}
	}
	fclose(fpSim);


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Allocate class median calculation area and values area
	//
	//
	double **ppClassMedianCalc = new double * [maxIndex];
	double **ppClassMedianVals = new double * [maxIndex];
	for (int i=0; i<maxIndex; i++ )
	{
		// this data is for calculating the medians
		if (ClassCount[i] > 0)
		{
			ppClassMedianCalc[i] = new double [ClassCount[i]];
		}
		else
		{
			ppClassMedianCalc[i] = NULL;
		}

		// this data holds the characteristic transform vectors (based on the median)
		ppClassMedianVals[i] = new double [nPreds];
		for ( int j=0; j<nPreds; j++ ) ppClassMedianVals[i][j] = 0.0;
	}


	//
	// determine the median for each predictor within each class sample
	//
	int *pCurrent = new int [maxIndex];
	for (int thisPred=0; thisPred<nPreds; thisPred++)
	{
		// init current indices to zero
		for ( int i=0; i<maxIndex; i++ ) pCurrent[i] = 0;

		// load the data into the class vectors for this pred
		for ( int i=0; i<nRows; i++ )
		{
			int nClass = pClusterID[i];
			ppClassMedianCalc[nClass-1][pCurrent[nClass-1]] = ppData[i][thisPred];
			++pCurrent[nClass-1];
		}

		// sort the vectors
		for ( int i=0; i<maxIndex; i++ )
		{
			if (ClassCount[i] > 0)
				qsort( (void *)ppClassMedianCalc[i], (size_t)ClassCount[i], sizeof(double), dCompare );
		}

		// extract the median
		for ( int i=0; i<maxIndex; i++ )
		{
			if (ClassCount[i] < 1)
			{
				ppClassMedianVals[i][thisPred] = 0.0;
			}
			else if (1 == ClassCount[i] % 2)  // odd number of samples for this class
			{
				ppClassMedianVals[i][thisPred] = ppClassMedianCalc[i][ClassCount[i]/2];
			}
			else // even number of samples for this class
			{
				int lowerIndex = (ClassCount[i]/2) - 1;
				int upperIndex = lowerIndex + 1;
				ppClassMedianVals[i][thisPred] = ppClassMedianCalc[i][lowerIndex] + 
					                                         ((ppClassMedianCalc[i][upperIndex] - ppClassMedianCalc[i][lowerIndex]) / 2.0);
			}
		}
	}


	//
	// Create an inter-class distance matrix that uses the MEDIAN between classes as the distance metric
	//
	double **DistMatrix_Median = new double * [maxIndex];
	for ( int i=0; i<maxIndex; i++ )
	{
		DistMatrix_Median[i] = new double [maxIndex];

		for ( int j=0; j<maxIndex; j++ )
		{
			DistMatrix_Median[i][j] = 1.0;

			if (ClassCount[i] < 1) 
				DistMatrix_Median[i][j] = 0.0;
		}
	}


	//
	// Use the GDM link function to calculate similarities from medians
	//
	for ( int i=0; i<maxIndex; i++)
	{
		for ( int j=0; j<maxIndex; j++)
		{
			if ((ClassCount[i] < 1) || (ClassCount[j] < 1))
			{
				DistMatrix_Median[i][j] = DistMatrix_Median[j][i] = 0.0;				
			}
			
			else if ( i != j)
			{
				double dVal = 0.0;
				for ( int k=0; k<nPreds; k++ )
				{
					dVal += fabs(ppClassMedianVals[i][k] - ppClassMedianVals[j][k]);
				}

				double dDist = exp(-(dVal));
				DistMatrix_Median[i][j] = DistMatrix_Median[j][i] = dDist;
			}
		}
	}


	//
	// Write distance matrix to .csv
	//
	sprintf(namestring, "%s\\%s_Sim_Median.csv", GmPath.GetDirectoryPath(OutputTablePath), GmPath.GetName(OutputTablePath));
	FILE *fpSimMedian = fopen( namestring, "w+t" );
	for ( int i=0; i<maxIndex; i++ )
	{
		for ( int j=0; j<maxIndex; j++ )
		{
			fprintf(fpSimMedian, "%lf", DistMatrix_Median[i][j]);
			if ( j < maxIndex-1 ) 
				fprintf(fpSimMedian, ",");
			else
				fprintf(fpSimMedian, "\n");
		}
	}
	fclose(fpSimMedian);


	//
	// Extract some distance stats
	//
	int MinMeanClass_A = 0;
	int MinMeanClass_B = 0;
	int MaxMeanClass_A = 0;
	int MaxMeanClass_B = 0;
	int MinMedianClass_A = 0;
	int MinMedianClass_B = 0;
	int MaxMedianClass_A = 0;
	int MaxMedianClass_B = 0;
	double dMinMean = 100.0;
	double dMaxMean = 0.0;
	double dMinMedian = 100.0;
	double dMaxMedian = 0.0;
	
	for ( int i=0; i<maxIndex; i++ )
	{
		for ( int j=0; j<maxIndex; j++ )
		{
			if (i != j)
			{
				if (DistMatrix_Mean[i][j] < dMinMean)
				{
					dMinMean = DistMatrix_Mean[i][j];
					MinMeanClass_A = i+1;
					MinMeanClass_B = j+1;
				}

				if (DistMatrix_Mean[i][j] > dMaxMean)
				{
					dMaxMean = DistMatrix_Mean[i][j];
					MaxMeanClass_A = i+1;
					MaxMeanClass_B = j+1;
				}


				if (DistMatrix_Median[i][j] < dMinMedian)
				{
					dMinMedian = DistMatrix_Median[i][j];
					MinMedianClass_A = i+1;
					MinMedianClass_B = j+1;
				}

				if (DistMatrix_Median[i][j] > dMaxMedian)
				{
					dMaxMedian = DistMatrix_Median[i][j];
					MaxMedianClass_A = i+1;
					MaxMedianClass_B = j+1;
				}
			}
		}
	}


	sprintf(namestring, "%s\\%s_Sim_Stats.csv", GmPath.GetDirectoryPath(OutputTablePath), GmPath.GetName(OutputTablePath));
	FILE *fpSimStats = fopen( namestring, "w+t" );
	fprintf(fpSimStats, 
		    "############################################\nStats for %d Similarity Classes for %s\n############################################\n",
			maxIndex, GmPath.GetName(OutputTablePath));

	fprintf(fpSimStats, 
		    "Minimum Similarity between classes %d and %d for mean class transform values of %lf\n\n", 
			MinMeanClass_A, MinMeanClass_B, dMinMean );

	fprintf(fpSimStats, 
		    "Maximum Similarity between classes %d and %d for mean class transform values of %lf\n\n", 
			MaxMeanClass_A, MaxMeanClass_B, dMaxMean );

	fprintf(fpSimStats, 
		    "Minimum Similarity between classes %d and %d for median class transform values of %lf\n\n", 
			MinMedianClass_A, MinMedianClass_B, dMinMedian );

	fprintf(fpSimStats, 
		    "Maximum Similarity between classes %d and %d for median class transform values of %lf\n\n", 
			MaxMedianClass_A, MaxMedianClass_B, dMaxMedian );

	fclose(fpSimStats);


	fptr("Doing GDM Classification...", 0);	

	//
	// Clean up
	//
	if (pCurrent) delete[] pCurrent;
	for ( int i=0; i<maxIndex; i++ ) 
	{
		if (DistMatrix_Mean[i]) delete[] DistMatrix_Mean[i];
		if (DistMatrix_Median[i]) delete[] DistMatrix_Median[i];		
		if (ppClassMedianCalc[i]) delete[] ppClassMedianCalc[i];
		if (ppClassMedianVals[i]) delete[] ppClassMedianVals[i];
	}
	if (DistMatrix_Mean) delete[] DistMatrix_Mean;
	if (DistMatrix_Median) delete[] DistMatrix_Median;
	if (ppClassMedianCalc) delete[] ppClassMedianCalc;
	if (ppClassMedianVals) delete[] ppClassMedianVals;
	if (ClassCount) delete[] ClassCount;
	if (pSplines) delete[] pSplines;
	for ( int i=0; i<nPreds; i++ )
	{
		if (ppQuantiles[i]) delete[] ppQuantiles[i];
		if (ppCoeffs[i]) delete[] ppCoeffs[i];
	}
	if (ppQuantiles) delete[] ppQuantiles;
	if (ppCoeffs) delete[] ppCoeffs;
	if (pClusterID) delete[] pClusterID;
	for ( int i=0; i<nRows; i++ )
		if (ppData[i]) delete[] ppData[i];
	if (ppData) delete[] ppData;
	return(true);
}


//
// Perform an unsupervised classification from a user defined table
//
bool ClassifyFromTable(char *pParams, 
                       char *TablePath, 
                       int  nSamples, 
                       int  nPreds, 
                       char *lpDomainPath, 
                       char *lpOutName, 
                       FPTR fptr)
{
	int nRows = nSamples;
	int nCols = nPreds;

	//
	// allocate table and vector data
	//
	int *pClusterID = new int [nRows];
	double **ppData = new double * [nRows];
	for (int i=0; i<nRows; i++ )
	{
		ppData[i] = new double [nCols];
	}
	
	//
	// extract training data from table 
	//
	fptr("Extracting Table Data...", 0);	
	ExtractDataFromTable(TablePath, nRows, nCols, pClusterID, ppData, fptr);

	//
	// use the training data to create an unsupervised classification grid
	//
	fptr("Doing Classification...", 0);	
	if ( false == DoGridClassification( pParams, lpDomainPath, 
		                                pClusterID, ppData, 
										nRows, nCols, 
										lpOutName, false, 
										fptr ) )
	{
		Message( "ClassifyFromTable returned false", "ClassifyFromTable" );
		if ( pClusterID ) delete[] pClusterID;
		for ( int i=0; i<nRows; i++ ) if ( ppData[i] ) delete ppData[i];
		if ( ppData ) delete[] ppData;
		return( false );
	}


	//
	// Get the maximum class id
	//
	int MaxClass = 0;
	for (int i=0; i<nRows; i++ )
	{
		if (pClusterID[i] > MaxClass) MaxClass = pClusterID[i];
	}


	//
	// cleanup
	//
	if (pClusterID) delete[] pClusterID;
	for ( int i=0; i<nRows; i++ ) if ( ppData[i] ) delete ppData[i];
	if ( ppData ) delete[] ppData;
	char pWorkspace[BUFFLEN];
	GetProfileString( "GDMODEL", "WorkspacePath", pWorkspace, pParams );


	//
	// write the metadata for this classification
	//
	char myFilePath[BUFFLEN];
	sprintf(myFilePath, "%s_meta.txt", lpOutName);
	FILE *fpMeta = fopen(myFilePath, "w+t");
	fprintf(fpMeta, "/////////////////////////////////////////////////////////////////////////////////////////////////\n");
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "Metadata for unsupervised classification of a GDM model from a training data table.\n");
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "/////////////////////////////////////////////////////////////////////////////////////////////////\n");
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "GDM Parameter File: %s\n", pParams);
	fprintf(fpMeta, "Classification Domain: %s\n", lpDomainPath);
	fprintf(fpMeta, "Training Data Path: %s\n", TablePath);
	fprintf(fpMeta, "Output Path: %s\\%s\n", pWorkspace, lpOutName);
	fprintf(fpMeta, "Number of Classes: %d\n", MaxClass);
	fprintf(fpMeta, "Number of Training Sites: %d\n", nRows);
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "/////////////////////////////////////////////////////////////////////////////////////////////////\n");
	if (fpMeta) fclose(fpMeta);


	//
	// write the class lookup table used in the BFT for this number of classes
	//
	sprintf(myFilePath, "%s_BFT_Class_Table.csv", lpOutName);
	FILE *fpClassLU = fopen(myFilePath, "w+t");

	fprintf(fpClassLU, "CLASS,CLASS_NAME,FORMATION,FORM_CODE,FORM_ID\n");
	fprintf(fpClassLU, "0,Cleared,Cleared,0,0\n");
	for ( int i=1; i<= MaxClass; i++ )
	{
		fprintf(fpClassLU, "%d,Class_%d,Class_%d,%d,%d\n", i,i,i,i,i );
	}
	if (fpClassLU) fclose(fpClassLU);


	//
	// finally do the classification colouring to produce
	// a suite of recoloured classifications
	//
	if ( false == DoClassificationColoring( pParams, TablePath, MaxClass, lpOutName, false, fptr ) )
	{
		Message( "DoClassificationColoring returned false", "DoGDMClassification" );
	}

	return(true);
}



//
// extract training data from table 
//
bool ExtractDataFromTable(char *TablePath, int nRows, int nCols, int *pClusterID, double **ppData, FPTR fptr)
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
		
		// site 
		p = strtok(NULL, seps);       // get X 
		p = strtok(NULL, seps);       // get Y

		// preds...
		for ( int j=0; j<nCols; j++ )
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
// Create a Training Data Set for user re-definition of the CLASS field
// prior to running an unsupervised classification with this training table
//
bool CreateClassTrainingTable(char *pParams, char *lpDomainPath, char *lpPredDomainPath, char *lpOutPath, int NumSamples, FPTR fptr )
{
	//
	// Extract a set of samples points from the domain.
	//
	if (!CreateSamplePointMesh( pParams, lpOutPath, lpDomainPath, NumSamples, false, fptr ))
	{
		Message("Cannot CreateSamplePointFile()", "CreateClassTrainingTable");
		return(false);
	}
	

	//
	// Extract transform grid data at the sample point locations
	// into a matrix to pass into the classification routines
	//
	int nRows, nCols;
	double **ppData = ExtractDataFromTransformGrids( pParams, lpPredDomainPath, lpOutPath, &nRows, &nCols, false );
	if ( NULL == ppData )
	{
		Message( "ExtractDataFromTransformGrids failure", "CreateClassTrainingTable" );
		return( false );
	}
	
	
	int *pClusterID = new int [nRows];
	for ( int i=0; i<nRows; i++ )
	{
		pClusterID[i] = 0;
	}


	// get the sample point coordinates
	char pPointsFile[BUFFLEN];
	sprintf( pPointsFile, "%s", lpOutPath );
	double *pSampleX = new double [nRows];
	double *pSampleY = new double [nRows];
	GetSampleCoords(pPointsFile, pSampleX, pSampleY, nRows);


	// write the training file for debug and for use in deriving the colors
	DumpANNFile( pParams, pClusterID, pPointsFile, pSampleX, pSampleY, ppData, nRows, nCols );

	//
	// cleanup
	//
	if (pClusterID) delete[] pClusterID;
	for ( int i=0; i<nRows; i++ ) if ( ppData[i] ) delete ppData[i];
	if ( ppData ) delete[] ppData;
	if (pSampleX) delete[] pSampleX;
	if (pSampleY) delete[] pSampleY;
	return(true);
}



//
// Do a sanity check to make sure that the cell size of the transform grids and the domain are the same
// and that the domain extent is equal to or is a valid interior extent to the transform grid extents.
//
bool DomainIsValidSubGridForClassification(char *SubDomainPath, char *TransformPath)
{
	return (true);
	//gmpath gmPath;
	//EsriBinaryHeader *headerDomain = new EsriBinaryHeader(gmPath.ChangeExtension(SubDomainPath, ".hdr"));	
	//EsriBinaryHeader *headerTransf = new EsriBinaryHeader(gmPath.ChangeExtension(TransformPath, ".hdr"));	

	//// check gridcell size
	//if (headerDomain->GetCellSize() != headerTransf->GetCellSize())
	//{
	//	delete headerDomain;
	//	delete headerTransf;
	//	return(false);
	//}

	//// check that the subdomain has a valid sub-extent of the transform grid
	//if ((headerDomain->GetMinX() - headerTransf->GetMinX()) < 0.0)
	//{
	//	delete headerDomain;
	//	delete headerTransf;
	//	return(false);
	//}

	//if ((headerDomain->GetMaxX() - headerTransf->GetMaxX()) > 0.0)
	//{
	//	delete headerDomain;
	//	delete headerTransf;
	//	return(false);
	//}

	//if ((headerDomain->GetMinY() - headerTransf->GetMinY()) < 0.0)
	//{
	//	delete headerDomain;
	//	delete headerTransf;
	//	return(false);
	//}

	//if ((headerDomain->GetMaxY() - headerTransf->GetMaxY()) > 0.0)
	//{
	//	delete headerDomain;
	//	delete headerTransf;
	//	return(false);
	//}

	//delete headerDomain;
	//delete headerTransf;
 //   return (true);
}


//
// A Wrapper function to call DoGDMClassification via params provided via a parameter file
//
bool DoANNClassificationViaParamFile(char *ParamFilePath, FPTR fptr)
{
	char lpGDMParamPath[BUFF1024];
	GetProfileString("CONFIGURATION", "GDMPARAMS", lpGDMParamPath, ParamFilePath);

	char lpDateTime[BUFF256];
	GetProfileString("CONFIGURATION", "LASTCHANGE", lpDateTime, ParamFilePath);

	char lpDomainPath[BUFF1024];
	if (GetProfileInt("INPUTS", "HaveMask", ParamFilePath) == 1)
		GetProfileString("INPUTS", "MaskPath", lpDomainPath, ParamFilePath);
	else if (false == GetTransformGridAsDomain(ParamFilePath, lpDomainPath))
	{
		Message("Cannot extract transform domain path", "DoANNClassificationViaParamFile");
		return(false);
	}

	char lpOutName[BUFF256];
	GetProfileString("OUTPUTS", "OutName", lpOutName, ParamFilePath);

	char lpOutDirectory[BUFF1024];
	GetProfileString("OUTPUTS", "OutDirectory", lpOutDirectory, ParamFilePath);

	int nClasses = GetProfileInt("INPUTS", "NumClasses", ParamFilePath);
	int NumSamples = GetProfileInt("INPUTS", "NumSamples", ParamFilePath);

	return(DoGDMClassification(lpGDMParamPath, lpDomainPath, 
		                       lpOutName, lpOutDirectory, 
		                       nClasses, NumSamples, 
		                       false, lpDateTime, fptr));
}



//
// Main GDM Unsupervised Classification called from .NET interface
//
bool DoGDMClassification(char *pParams, 
						 char *lpDomainPath,
						 char *lpOutName,
						 char *lpOutPath,
   					     int nClasses, 
						 int NumSamples,
						 bool DoBatch,
						 char *TimeString,
						 FPTR fptr )
{
	int NumClasses = nClasses;

	//
	// Extract a set of samples points for the nearest neighbour classification and write to
	// a comma delimited text file in the Working directory called ClassificationSample.csv
	//
	char lpSamplePath[FILEPATH_BUFFSIZE];
	sprintf(lpSamplePath, "%s\\%s.csv", lpOutPath, lpOutName);
	if (!CreateSamplePointMesh( pParams, lpSamplePath, lpDomainPath, NumSamples, DoBatch, fptr ))
	{
		if (!DoBatch)
			Message("Cannot CreateSamplePointFile()", "DoGDMClassification");
		return(false);
	}


	//
	// Extract transform grid data at the sample point locations
	// into a matrix to pass into the classification routines
	//
	int nRows, nCols;
	double **ppData = ExtractDataFromTransformGrids( pParams, lpDomainPath, lpSamplePath, &nRows, &nCols, DoBatch );
	if ( NULL == ppData )
	{
		if (!DoBatch)
			Message( "ExtractDataFromTransformGrids failure", "DoClassification" );
		return( false );
	}

	//
	// The sample point code may NOT always return the exact number of points requested.
	// This could present a problem if the number of classes was GREATER than the number of sample points.
	// So if the number of classes is greater than the actual number of samples (rows) then adjust NumClasses.
	//
	if (NumClasses > nRows ) NumClasses = nRows;
	

	//
	// Do a grid based nearest neighbour classification using the transform grids generated 
	// by the GDM process and the comma delimited table written in CreateSamplePointFile().
	//
	fptr("Doing Hierachical Clustering...", 25);	
	int *pClusterID = GetHierachicalClusterClasses( ppData, nRows, nCols, NumClasses, DoBatch, fptr );


	// get the sample point coordinates
	double *pSampleX = new double [nRows];
	double *pSampleY = new double [nRows];
	GetSampleCoords(lpSamplePath, pSampleX, pSampleY, nRows);

	// write the training file for debug and for use in deriving the colors
	DumpANNFile( pParams, pClusterID, lpSamplePath, pSampleX, pSampleY, ppData, nRows, nCols );

	
	//
	// cleanup
	//
	if (pSampleX) delete[] pSampleX;
	if (pSampleY) delete[] pSampleY;


	//
	// use the training data to create an unsupervised classification grid
	//
	fptr("Doing Classification...", 75);	
	if ( false == DoGridClassification( pParams, lpDomainPath, 
		                                pClusterID, ppData, 
										nRows, nCols, 
										lpSamplePath, 
										DoBatch, 
										fptr ) )
	{
		if (!DoBatch)
			Message( "DoGridClassification returned false", "DoGDMClassification" );
		if ( pClusterID ) delete[] pClusterID;
		for ( int i=0; i<nRows; i++ ) if ( ppData[i] ) delete ppData[i];
		if ( ppData ) delete[] ppData;
		return( false );
	}


	//
	// write the metadata for this classification
	//
	char myFilePath[BUFFLEN];
	sprintf(myFilePath, "%s\\%s_meta.txt", lpOutPath, lpOutName);
	FILE *fpMeta = fopen(myFilePath, "w+t");
	fprintf(fpMeta, "/////////////////////////////////////////////////////////////////////////////////////////////////\n");
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "Metadata for unsupervised classification of a GDM model: %s\n", TimeString);
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "/////////////////////////////////////////////////////////////////////////////////////////////////\n");
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "GDM Parameter File: %s\n", pParams);
	fprintf(fpMeta, "Classification Domain: %s\n", lpDomainPath);
	fprintf(fpMeta, "Output Path: %s\\%s\n", lpOutPath, lpOutName);
	fprintf(fpMeta, "Number of Classes: %d\n", NumClasses);
	fprintf(fpMeta, "Number of Training Sites: %d\n", nRows);
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "/////////////////////////////////////////////////////////////////////////////////////////////////\n");
	if (fpMeta) fclose(fpMeta);


	////
	//// write the class lookup table used in the BFT for this number of classes
	////
	//sprintf(myFilePath, "%s\\%s_BFT_Class_Table.csv", lpOutPath, lpOutName);
	//FILE *fpClassLU = fopen(myFilePath, "w+t");

	//fprintf(fpClassLU, "CLASS,CLASS_NAME,FORMATION,FORM_CODE,FORM_ID\n");
	//fprintf(fpClassLU, "0,Cleared,Cleared,0,0\n");
	//for ( int i=1; i<= NumClasses; i++ )
	//{
	//	fprintf(fpClassLU, "%d,Class_%d,Class_%d,%d,%d\n", i,i,i,i,i );
	//}
	//if (fpClassLU) fclose(fpClassLU);


	//
	// Do the classification colouring to produce a suite of recoloured classifications
	//
	if ( false == DoClassificationColoring( pParams, lpSamplePath, NumClasses, lpOutName, DoBatch, fptr ) )
	{
		if (!DoBatch)
			Message( "DoClassificationColoring returned false", "DoGDMClassification" );
	}

	
	////
	//// Do post classification sample table for creation of inter-class distance matrix
	////
	//if ( false == DoPostClassificationMatrix(pParams, lpDomainPath, lpOutName, lpSamplePath, DoBatch, fptr))
	//{
	//	if (!DoBatch)
	//		Message( "DoPostClassificationMatrix returned false", "DoGDMClassification" );
	//}


	//
	// Test for Dan Faith's problem with clustering of demand points
	//
	/*if ( false == DoDemandPoints_02(pParams, lpDomainPath, lpOutName, lpSamplePath, DoBatch, fptr))
	{
		if (!DoBatch)
			Message( "DoPostClassificationMatrix returned false", "DoGDMClassification" );
	}*/

	
	//
	// write a csv table derived from the above for use in R for creating Dendrograms
	// It will contain a record for each class, followed by the predictor data that is the 
	// average of each predictor in the respective class in the table created in DumpANNFile()
	//
	char lpKMPath[FILEPATH_BUFFSIZE];
	sprintf(lpKMPath, "%s\\%s_KMeans4RDend.csv", lpOutPath, lpOutName);
	WriteDendrogramTable(pParams, pClusterID, lpKMPath, ppData, nRows, nCols);


	//
	// write a table with the closest (ED) sites in the GDM classes 
	// to the K-Median of each class written in %s\\%s_KMeans4RDend.csv
	//
	/*char lpClassPath[BUFFLEN];
	sprintf(lpClassPath, "%s\\%s", lpOutPath, lpOutName);
	char lpOutTablePath[BUFFLEN];
	sprintf(lpOutTablePath, "%s\\%s_Closest2KMeans.csv", lpOutPath, lpOutName);
	CreateKMeansClosestPointTable(pParams, lpClassPath, lpKMPath, lpOutTablePath, nClasses, fptr);

	sprintf(lpOutTablePath, "%s\\%s_KMeansDemandSites.csv", lpOutPath, lpOutName);
	CreateKMeansDemandPointTable(pParams, lpClassPath, lpKMPath, lpOutTablePath, nClasses, fptr);*/

	//
	// clean up
	//
	if (pClusterID) delete[] pClusterID;
	for (int i = 0; i<nRows; i++) if (ppData[i]) delete ppData[i];
	if (ppData) delete[] ppData;
	return(true);
}



//
// Do post classification sample table for creation of inter-class distance matrix
//
bool DoPostClassificationMatrix(char *pParams, 
	                            char *lpDomainPath, 
	                            char *lpOutName, 
								char *lpSamplePath, 
								bool DoBatch, 
								FPTR fptr)
{
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// This is a post classification step where we sample the classification grid 
	// for NUMSAMPLES samples and get the class and the transform values.
	// This data will be used to extract a charanteristic vector from the 
	// transform values where each class will be represented a predictor
	// value equal to the median of the transform values within each class.
	//
	//
	// Extract a set of samples points for the nearest neighbour classification and write to
	// a comma delimited text file in the Working directory called (Classification NAME)_Sample.csv
	//
	gmpath GmPath;
	char namestring [FILEPATH_BUFFSIZE];
	int NUMSAMPLES = 50000;
	sprintf(namestring, "%s\\%s_Sample.csv", GmPath.GetDirectoryPath(lpSamplePath), GmPath.GetName(lpSamplePath));
	if (!CreateSamplePointMesh( pParams, namestring, lpDomainPath, NUMSAMPLES, false, fptr ))
	{
		if (!DoBatch)
			Message("Cannot CreateSamplePointFile()", "DoGDMClassification");
		return(false);
	}

	//
	// Extract transform grid data at the sample point locations
	// into a matrix to pass into the classification routines
	//
	int nRows = 0;
	int nCols = 0;
	double **ppData = ExtractDataFromTransformGrids( pParams, lpDomainPath, namestring, &nRows, &nCols, DoBatch );
	if ( NULL == ppData )
	{
		if (!DoBatch)
			Message( "ExtractDataFromTransformGrids failure", "DoClassification" );
		return( false );
	}


	//
	// Get the sample X,Y's
	// 
	double *pSampleX = new double [nRows];
	double *pSampleY = new double [nRows];
	GetSampleCoords(namestring, pSampleX, pSampleY, nRows);

	//
	// delete the points file now we don't need it
	//
	if (GmPath.FileExists(namestring))
	{
		//Message(namestring, "namestring");
		remove(namestring);
	}


	//
	// Get the class ID's
	//
	int *pClassID = GetClassIDs(pSampleX, pSampleY, nRows, pParams, lpOutName);

	//
	// Now write the table of sites,class,transformed predictor values...
	//
	/*sprintf(namestring, "%s\\%s_XYClass.csv", GmPath.GetDirectoryPath(lpSamplePath), GmPath.GetName(lpSamplePath));
	FILE *fpDist = fopen( namestring, "w+t" );
	fprintf(fpDist, "X,Y,Class");
	int nPreds = GetProfileInt( "PREDICTORS", "NumPredictors", pParams );
	char lpKey[64];
	char lpPredPath[BUFFLEN];
	for (int i=1; i<=nPreds; i++ )
	{
		sprintf( lpKey, "PredTran%d", i );
		GetProfileString( "TRANSPREDS", lpKey, lpPredPath, pParams );

		if ( strlen( lpPredPath ) > 0 )
		{
			sprintf(lpPredPath, "%s", GmPath.ChangeExtension( lpPredPath, ".flt" ));
			fprintf(fpDist, ",%s", GmPath.GetName(lpPredPath));
		}

		if (i == nPreds)
			fprintf(fpDist, "\n");
	}
	for ( int i=0; i<nRows; i++ )
	{
		fprintf(fpDist, "%lf,%lf,%d", pSampleX[i], pSampleY[i], pClassID[i]);
		for ( int j=0; j<nCols; j++ )
		{
			fprintf(fpDist, ",%lf", ppData[i][j]);
			if ( j == nCols-1 )
				fprintf(fpDist, "\n");
		}
	}
	fclose(fpDist);*/

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Get the maximum class index (total number of classes)
	//
	int maxIndex = 0;
	for ( int i=0; i<nRows; i++ )
	{
		if ( pClassID[i] > maxIndex )
		{
			maxIndex = pClassID[i];
		}
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Allocate class counters (count the number of samples in each class)
	//
	int *ClassCount = new int [maxIndex];
	for ( int i=0; i<maxIndex; i++ )
	{
		ClassCount[i] = 0;
	}

	//
	// Initialise class counters
	//
	for ( int i=0; i<nRows; i++ )
	{
		++ClassCount[pClassID[i]-1];
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Allocate class averages [ClassIndex][Sum then divide by count for average]
	//
	double **ppClassAverages = new double * [maxIndex];
	for (int i=0; i<maxIndex; i++ )
	{
		ppClassAverages[i] = new double [nCols];

		// init to zero
		for ( int j=0; j<nCols; j++ ) ppClassAverages[i][j] = 0.0;
	}


	//
	// Get the sum of the transform values for each class
	//
	for (int i=0; i<nRows; i++ )
	{
		for ( int j=0; j<nCols; j++ )
		{
			ppClassAverages[pClassID[i]-1][j] += ppData[i][j];
		}
	}


	//
	// Divide by the class count to get the average
	//
	for (int i=0; i<maxIndex; i++ )
	{
		for ( int j=0; j<nCols; j++ )
		{
			if (ClassCount[i] < 1) 
				ppClassAverages[i][j] = 0.0;
			else
				ppClassAverages[i][j] /= ClassCount[i];
		}
	}


	//
	// Write class counts and averages
	//
	/*sprintf(namestring, "%s\\%s_ClassCount.csv", GmPath.GetDirectoryPath(lpSamplePath), GmPath.GetName(lpSamplePath));
	FILE *fpCount = fopen( namestring, "w+t" );
	fprintf(fpCount, "Class,Count");
	for ( int j=0; j<nCols; j++ )
	{
		fprintf(fpCount, ",Pred_0%d", j+1);
		if ( j == nCols-1 ) fprintf(fpCount, "\n");
	}
	
	for ( int i=0; i<maxIndex; i++ )
	{
		fprintf(fpCount, "%d,%d", i+1, ClassCount[i]);
		for ( int j=0; j<nCols; j++ )
		{
			fprintf(fpCount, ",%lf", ppClassAverages[i][j]);
			if ( j == nCols-1 ) fprintf(fpCount, "\n");
		}
	}
	fclose(fpCount);*/


	//
	// Create an inter-class distance matrix that uses the MEAN between classes as the distance metric
	//
	double **DistMatrix_Mean = new double * [maxIndex];
	for ( int i=0; i<maxIndex; i++ )
	{
		DistMatrix_Mean[i] = new double [maxIndex];
		for ( int j=0; j<maxIndex; j++ )
		{
			DistMatrix_Mean[i][j] = 1.0;

			if (ClassCount[i] < 1) 
				DistMatrix_Mean[i][j] = 0.0;
		}
	}


	//
	// Use the GDM link function to calculate similarities
	//
	for ( int i=0; i<maxIndex; i++)
	{
		for ( int j=0; j<maxIndex; j++)
		{
			if ((ClassCount[i] < 1) || (ClassCount[j] < 1))
			{
				DistMatrix_Mean[i][j] = DistMatrix_Mean[j][i] = 0.0;
			}
			
			else if ( i != j)
			{
				double dVal = 0.0;
				for ( int k=0; k<nCols; k++ )
				{
					dVal += fabs(ppClassAverages[i][k] - ppClassAverages[j][k]);
				}

				double dDist = exp(-(dVal));
				DistMatrix_Mean[i][j] = DistMatrix_Mean[j][i] = dDist;
			}
		}
	}


	//
	// Write distance matrix to .csv
	//
	sprintf(namestring, "%s\\%s_Sim_Mean.csv", GmPath.GetDirectoryPath(lpSamplePath), GmPath.GetName(lpSamplePath));
	FILE *fpSim = fopen( namestring, "w+t" );
	for ( int i=0; i<maxIndex; i++ )
	{
		for ( int j=0; j<maxIndex; j++ )
		{
			fprintf(fpSim, "%lf", DistMatrix_Mean[i][j]);
			if ( j < maxIndex-1 ) 
				fprintf(fpSim, ",");
			else
				fprintf(fpSim, "\n");
		}
	}
	fclose(fpSim);


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Allocate class median calculation area and values area
	//
	//
	double **ppClassMedianCalc = new double * [maxIndex];
	double **ppClassMedianVals = new double * [maxIndex];
	for (int i=0; i<maxIndex; i++ )
	{
		// this data is for calculating the medians
		if (ClassCount[i] > 0)
		{
			ppClassMedianCalc[i] = new double [ClassCount[i]];
		}
		else
		{
			ppClassMedianCalc[i] = NULL;
		}

		// this data holds the characteristic transform vectors (based on the median)
		ppClassMedianVals[i] = new double [nCols];
		for ( int j=0; j<nCols; j++ ) ppClassMedianVals[i][j] = 0.0;
	}


	//
	// determine the median for each predictor within each class sample
	//
	int *pCurrent = new int [maxIndex];
	for (int thisPred=0; thisPred<nCols; thisPred++)
	{
		// init current indices to zero
		for ( int i=0; i<maxIndex; i++ ) pCurrent[i] = 0;

		// load the data into the class vectors for this pred
		for ( int i=0; i<nRows; i++ )
		{
			int nClass = pClassID[i];
			ppClassMedianCalc[nClass-1][pCurrent[nClass-1]] = ppData[i][thisPred];
			++pCurrent[nClass-1];
		}

		// sort the vectors
		for ( int i=0; i<maxIndex; i++ )
		{
			if (ClassCount[i] > 0)
				qsort( (void *)ppClassMedianCalc[i], (size_t)ClassCount[i], sizeof(double), dCompare );
		}

		// extract the median
		for ( int i=0; i<maxIndex; i++ )
		{
			if (ClassCount[i] < 1)
			{
				ppClassMedianVals[i][thisPred] = 0.0;
			}
			else if (1 == ClassCount[i] % 2)  // odd number of samples for this class
			{
				ppClassMedianVals[i][thisPred] = ppClassMedianCalc[i][ClassCount[i]/2];
			}
			else // even number of samples for this class
			{
				int lowerIndex = (ClassCount[i]/2) - 1;
				int upperIndex = lowerIndex + 1;
				ppClassMedianVals[i][thisPred] = ppClassMedianCalc[i][lowerIndex] + 
					                                         ((ppClassMedianCalc[i][upperIndex] - ppClassMedianCalc[i][lowerIndex]) / 2.0);
			}
		}
	}


	//
	// Write class predictor medians
	//
	/*sprintf(namestring, "%s\\%s_ClassMedian.csv", GmPath.GetDirectoryPath(lpSamplePath), GmPath.GetName(lpSamplePath));
	FILE *fpMedian = fopen( namestring, "w+t" );
	fprintf(fpMedian, "Class,Count");
	for ( int j=0; j<nCols; j++ )
	{
		fprintf(fpMedian, ",Pred_0%d", j+1);
		if ( j == nCols-1 ) fprintf(fpMedian, "\n");
	}
	
	for ( int i=0; i<maxIndex; i++ )
	{
		fprintf(fpMedian, "%d,%d", i+1, ClassCount[i]);
		for ( int j=0; j<nCols; j++ )
		{
			fprintf(fpMedian, ",%lf", ppClassMedianVals[i][j]);
			if ( j == nCols-1 ) fprintf(fpMedian, "\n");
		}
	}
	fclose(fpMedian);*/


	//
	// Create an inter-class distance matrix that uses the MEDIAN between classes as the distance metric
	//
	double **DistMatrix_Median = new double * [maxIndex];
	for ( int i=0; i<maxIndex; i++ )
	{
		DistMatrix_Median[i] = new double [maxIndex];

		for ( int j=0; j<maxIndex; j++ )
		{
			DistMatrix_Median[i][j] = 1.0;

			if (ClassCount[i] < 1) 
				DistMatrix_Median[i][j] = 0.0;
		}
	}


	//
	// Use the GDM link function to calculate similarities from medians
	//
	for ( int i=0; i<maxIndex; i++)
	{
		for ( int j=0; j<maxIndex; j++)
		{
			if ((ClassCount[i] < 1) || (ClassCount[j] < 1))
			{
				DistMatrix_Median[i][j] = DistMatrix_Median[j][i] = 0.0;				
			}
			
			else if ( i != j)
			{
				double dVal = 0.0;
				for ( int k=0; k<nCols; k++ )
				{
					dVal += fabs(ppClassMedianVals[i][k] - ppClassMedianVals[j][k]);
				}

				double dDist = exp(-(dVal));
				DistMatrix_Median[i][j] = DistMatrix_Median[j][i] = dDist;
			}
		}
	}


	//
	// Write distance matrix to .csv
	//
	sprintf(namestring, "%s\\%s_Sim_Median.csv", GmPath.GetDirectoryPath(lpSamplePath), GmPath.GetName(lpSamplePath));
	FILE *fpSimMedian = fopen( namestring, "w+t" );
	for ( int i=0; i<maxIndex; i++ )
	{
		for ( int j=0; j<maxIndex; j++ )
		{
			fprintf(fpSimMedian, "%lf", DistMatrix_Median[i][j]);
			if ( j < maxIndex-1 ) 
				fprintf(fpSimMedian, ",");
			else
				fprintf(fpSimMedian, "\n");
		}
	}
	fclose(fpSimMedian);


	//
	// Extract some distance stats
	//
	int MinMeanClass_A = 0;
	int MinMeanClass_B = 0;
	int MaxMeanClass_A = 0;
	int MaxMeanClass_B = 0;
	int MinMedianClass_A = 0;
	int MinMedianClass_B = 0;
	int MaxMedianClass_A = 0;
	int MaxMedianClass_B = 0;
	double dMinMean = 100.0;
	double dMaxMean = 0.0;
	double dMinMedian = 100.0;
	double dMaxMedian = 0.0;
	
	for ( int i=0; i<maxIndex; i++ )
	{
		for ( int j=0; j<maxIndex; j++ )
		{
			if (i != j)
			{
				if (DistMatrix_Mean[i][j] < dMinMean)
				{
					dMinMean = DistMatrix_Mean[i][j];
					MinMeanClass_A = i+1;
					MinMeanClass_B = j+1;
				}

				if (DistMatrix_Mean[i][j] > dMaxMean)
				{
					dMaxMean = DistMatrix_Mean[i][j];
					MaxMeanClass_A = i+1;
					MaxMeanClass_B = j+1;
				}


				if (DistMatrix_Median[i][j] < dMinMedian)
				{
					dMinMedian = DistMatrix_Median[i][j];
					MinMedianClass_A = i+1;
					MinMedianClass_B = j+1;
				}

				if (DistMatrix_Median[i][j] > dMaxMedian)
				{
					dMaxMedian = DistMatrix_Median[i][j];
					MaxMedianClass_A = i+1;
					MaxMedianClass_B = j+1;
				}
			}
		}
	}


	sprintf(namestring, "%s\\%s_Sim_Stats.csv", GmPath.GetDirectoryPath(lpSamplePath), GmPath.GetName(lpSamplePath));
	FILE *fpSimStats = fopen( namestring, "w+t" );
	fprintf(fpSimStats, 
		    "############################################\nStats for %d Similarity Classes for %s\n############################################\n",
			maxIndex, GmPath.GetName(lpSamplePath));

	fprintf(fpSimStats, 
		    "Minimum Similarity between classes %d and %d for mean class transform values of %lf\n\n", 
			MinMeanClass_A, MinMeanClass_B, dMinMean );

	fprintf(fpSimStats, 
		    "Maximum Similarity between classes %d and %d for mean class transform values of %lf\n\n", 
			MaxMeanClass_A, MaxMeanClass_B, dMaxMean );

	fprintf(fpSimStats, 
		    "Minimum Similarity between classes %d and %d for median class transform values of %lf\n\n", 
			MinMedianClass_A, MinMedianClass_B, dMinMedian );

	fprintf(fpSimStats, 
		    "Maximum Similarity between classes %d and %d for median class transform values of %lf\n\n", 
			MaxMedianClass_A, MaxMedianClass_B, dMaxMedian );

	fclose(fpSimStats);


	//
	// Cleanup
	//
	if (pCurrent) delete[] pCurrent;
	for ( int i=0; i<maxIndex; i++ ) 
	{
		if (DistMatrix_Mean[i]) delete[] DistMatrix_Mean[i];
		if (DistMatrix_Median[i]) delete[] DistMatrix_Median[i];		
		if (ppClassMedianCalc[i]) delete[] ppClassMedianCalc[i];
		if (ppClassMedianVals[i]) delete[] ppClassMedianVals[i];
	}
	if (DistMatrix_Mean) delete[] DistMatrix_Mean;
	if (DistMatrix_Median) delete[] DistMatrix_Median;
	if (ppClassMedianCalc) delete[] ppClassMedianCalc;
	if (ppClassMedianVals) delete[] ppClassMedianVals;
	if (ClassCount) delete[] ClassCount;
	if (pClassID) delete[] pClassID;
	if (pSampleX) delete[] pSampleX;
	if (pSampleY) delete[] pSampleY;
	for ( int i=0; i<nRows; i++ ) if ( ppData[i] ) delete ppData[i];
	if ( ppData ) delete[] ppData;
	return(true);
}


////////////////////////////////////////////////////////////////////////////////////
// String comparison routine for Sort method
//
int dCompare( const void *arg1, const void *arg2 )
{
	double *p1 = (double *)arg1;
	double *p2 = (double *)arg2;
	if ( *p1 < *p2 )
		return(-1);
	else if ( *p1 > *p2 )
		return(1);
	else
		return(0);
}



//
// Extract Class IDs from the classification Grid
//
int *GetClassIDs(double *pSampleX, double *pSampleY, int nRows, char *pParams, char *lpOutName)
{
	int *pClasses = new int [nRows];
	char lpWorkspace[BUFFLEN];
	char lpOutBin[BUFFLEN];
	char lpOutHdr[BUFFLEN];

	GetProfileString( "GDMODEL", "WorkspacePath", lpWorkspace, pParams );
	sprintf( lpOutBin, "%s\\%s\\%s", lpWorkspace, lpOutName, lpOutName );
	gmpath GmPath;
	sprintf(lpOutBin, "%s", GmPath.ChangeExtension(lpOutBin, ".flt"));
	sprintf(lpOutHdr, "%s", GmPath.ChangeExtension(lpOutBin, ".hdr"));

	BinaryFileClass *bfc = new BinaryFileClass(lpOutBin);
	EsriBinaryHeader *bfh = new EsriBinaryHeader(lpOutHdr);
	
	for ( int i=0; i<nRows; i++ )
	{
		double dX = pSampleX[i];
		double dY = pSampleY[i];

		long lOffset = ((long)(fabs(bfh->GetMaxY() - dY) / bfh->GetCellSize()) * bfh->GetNumCols()) + 
		               ((long)(fabs(dX - bfh->GetMinX()) / bfh->GetCellSize())); 

		bfc->SeekTo(lOffset*sizeof(float));
		float fVal;
		bfc->ReadFloat(&fVal,1);

		pClasses[i] = (int)fVal;
	}

	if (bfc) delete bfc;
	if (bfh) delete bfh;
	return(pClasses);
}



//
// Creates a sample point file from Floating Point Domain Grid
// 
bool CreateSamplePointMesh( char *lpParams, char *lpSamplePointPath, char *lpDomainPath, int nSamples, bool DoBatch, FPTR fptr )
{
	gmpath gmPath;
	//
	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(lpDomainPath, ".hdr"));
	//
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	double ResX = header->GetCellSize();
	float fNoData = header->GetNoDataValue();
	double dResX = header->GetCellSize();
	double dMinX = header->GetXllCorner();
	double dMinY = header->GetYllCorner();
	double dMaxX = dMinX + (dResX * nCols);
	double dMaxY = dMinY + (dResX * nRows);

	//
	// open the file for reading
	//
	char lpFLTPath[BUFFLEN];
	sprintf(lpFLTPath, "%s", gmPath.ChangeExtension( lpDomainPath, ".flt" ) );

	BinaryFileClass *bfcDomain = new BinaryFileClass(lpFLTPath, BFC_ReadOnly);
	if ( !bfcDomain->IsValid() )
	{
		if (!DoBatch)
			Message( "Cannot open Domain file for READ", "CreateSamplePointMesh" );
		if (header) delete header;
		return(false);
	}


	// allocate a block for the domain
	float **ppDomain = new float * [ nRows ];
	for ( int i=0; i<nRows; i++ )
	{
		ppDomain[i] = new float [ nCols ];
	}


	// count the valid data cells in the domain that are in the ROI
	int nValidCells = 0;
	float *pRows = new float [nCols];
	int x,y;
	for ( y=0; y<nRows; y++ )
	{
		// seek to the current row
		bfcDomain->SeekTo(y*nCols*sizeof(float));

		// read the current row
		bfcDomain->ReadFloat(ppDomain[y], nCols);

		// count valid cells
		for ( x=0; x<nCols; x++ )
		{
			//if ( (ppDomain[y][x] != fNoData) && (ppDomain[y][x] > 0.0F) )		// Threshold Set Here
			if (ppDomain[y][x] != fNoData)
			{
				++nValidCells;
			}
		}
	}

	bfcDomain->Close();
	if (header) delete header;


	//
	// if the number of valid cells is less than nSamples then use the valid cell count
	//
	if ( nValidCells < nSamples )
	{
		// create the classification sample point file
		FILE *fp = fopen( lpSamplePointPath, "w+t" );

		// write a header
		fprintf( fp, "_ID,X,Y\n" );
		int nThisPoint = 1;
		int nCurrent = 0;
		fptr("Create sample points...", nCurrent);	
		for ( y=0; y<nRows; y++ )
		{
			if ( y * 100 / nRows > nCurrent )
			{
				nCurrent = y * 100 / nRows;
				fptr("Create sample points...", nCurrent);	
			}

			for ( x=0; x<nCols; x++ )
			{
				//if ( (ppDomain[y][x] != fNoData) && (ppDomain[y][x] > 0.0F) )	// Threshold Set Here
				if (ppDomain[y][x] != fNoData)
				{
					// write the new point record
					fprintf( fp, 
						     "%d,%f,%f\n", 
						     nThisPoint++, 
							 dMinX + (x * dResX) + (dResX / 2.0),
							 dMaxY - (y * dResX) - (dResX / 2.0)
							 );
				}
			}
		}

		fclose( fp );
	}
	else // use a subsample 
	{		
		int nDemandCount = nSamples;
		int nTotalCells = nRows * nCols;
		float fCellProp = (float)nValidCells / (float)nTotalCells;
		int nScaledDemandCount = NRound( (float)nDemandCount / fCellProp );
		float fAspectRatio = (float)nCols / (float)nRows;
		float fRoot = sqrtf( (float)nScaledDemandCount / fAspectRatio);
		int nDemandWidth = (int)floorf( fRoot * fAspectRatio );
		int nDemandHeight = nScaledDemandCount / nDemandWidth;		
		int nDemandTotal = nDemandWidth * nDemandHeight;

		//
		// generate a file of the points here
		//	
		double dXInc = ( dMaxX - dMinX ) / nDemandWidth;
		double dXOffset = dXInc / 2;
		double dXCurrent = dMinX + dXOffset;

		double dYInc = ( dMaxY - dMinY ) / nDemandHeight;
		double dYOffset = dYInc / 2;
		double dYCurrent = dMinY + dYOffset;		

		int nThisPoint = 1;

		// create the classification sample point file
		FILE *fp = fopen( lpSamplePointPath, "w+t" );

		// write a header
		fprintf( fp, "_ID,X,Y\n" );
		int nCurrent = 0;
		fptr("Create sample points...", nCurrent);	
		for ( int yy=0; yy<nDemandHeight; yy++ )
		{
			if ( yy * 100 / nDemandHeight > nCurrent )
			{
				nCurrent = yy * 100 / nDemandHeight;
				fptr("Create sample points...", nCurrent);	
			}

			// reset for the current row
			dXCurrent = dMinX + dXOffset;
			
			for ( int xx=0; xx<nDemandWidth; xx++ )
			{
				if ( dXCurrent < dMaxX )
				{
					//
					// only write the point if it has valid data in the domain grid
					//
					int nXIndex = int( floor ( fabs( dXCurrent-dMinX ) / dResX ) );
					int nYIndex = int( floor ( fabs( dYCurrent-dMaxY ) / dResX ) );

					// calculate the linear memory block offset from the X and Y values
					//if ( ppDomain[ nYIndex ][ nXIndex ] > 0.0F )  // Threshold Set Here
					if (ppDomain[nYIndex][nXIndex] != fNoData)
					{
						//
						// recalculate coordinate to place it in the center of its gridcell
						//
						dXCurrent = dMinX + (nXIndex * dResX) + (dResX / 2);
						dYCurrent = dMaxY - (nYIndex * dResX) - (dResX / 2 );

						// write the new point record
						fprintf( fp, "%d,%f,%f\n", nThisPoint, dXCurrent, dYCurrent );

						// update for successive points on this row
						++nThisPoint;
					}
				}
				dXCurrent += dXInc;
			}
			// reset for the next row
			dYCurrent += dYInc; 
		}
		fclose( fp );
	}

	//
	// clean up 
	//
	for ( int i=0; i<nRows; i++ ) if ( ppDomain[i] ) delete[] ppDomain[i];
	if ( ppDomain ) delete[] ppDomain;

	return( true );
}



//
// takes a float and if floating point remainder < 0.5 then return floor as an int
// otherwise return ceiling as an int
//
int NRound( float f )
{
	float fTmp = fmodf( f, floorf( f ) );
	if ( fTmp < (float)0.5 )
		return( int( floorf( f ) ) );
	else
		return( int( ceilf( f ) ) );
}



//
// extract data values from transformed predictor 
// grids into a matrix of doubles, setting the row
// and column counts into the supplied parameters
//
double **ExtractDataFromTransformGrids( char *lpParams, char *lpDomainPath, char *lpSamplePath, int *pRows, int *pCols, bool DoBatch )
{
	//
	// are we using euclidean ?
	//
	bool fUseEuclidean = ( 1 == GetProfileInt( "GDMODEL", "UseEuclidean", lpParams ) ) ? true : false;

	//
	// determine the number of columns ( the number of transform grids )
	//
	int nCols = 0;
	if ( fUseEuclidean ) nCols += 2;

	// total number of possible preds, need to check for presence to get file transform count
	int nPreds = GetProfileInt( "PREDICTORS", "NumPredictors", lpParams );
	char lpKey[64];
	char lpPredPath[BUFFLEN];
	for ( int i=1; i<=nPreds; i++ )
	{
		sprintf( lpKey, "PredTran%d", i );
		GetProfileString( "TRANSPREDS", lpKey, lpPredPath, lpParams );

		if ( strlen( lpPredPath ) > 0 )
		{
			++nCols;
		}
	}

	// 
	// determine the number of rows ( the number of records in the points file )
	//
	int nRows = 0;

	// open the classification sample point file
	char buff[BUFFLEN];
	FILE *fp = fopen( lpSamplePath, "r+t" );
	// get the header
	fgets( buff, BUFFLEN, fp );
	
	while( 1 )
	{
		if ( NULL == fgets( buff, BUFFLEN, fp ) )
			break;

		++nRows;
	}



	gmpath gmPath;
	//
	//char *lpDomainPath = GetPredictorDomainPath(lpParams);
	//if (!gmPath.FileExists(lpDomainPath))
	if (!gmPath.FileExists(gmPath.ChangeExtension( lpDomainPath, ".flt" )))
	{
		// try using the first environmental grid in the predictor list
		GetProfileString( "PREDICTORS", "EnvGrid1", lpDomainPath, lpParams );

		// is this grid there?
		if (!gmPath.FileExists(gmPath.ChangeExtension( lpDomainPath, ".flt" )))
		{
			if (!DoBatch)
				Message("Cannot locate a domain grid", "ExtractDataFromTransformGrids");
			return(NULL);
		}
	}
	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(lpDomainPath, ".hdr"));
	//
	int nGridRows = header->GetNumRows();
	int nGridCols = header->GetNumCols();
	double ResX = header->GetCellSize();
	float fNoData = header->GetNoDataValue();
	double dResX = header->GetCellSize();
	double dMinX = header->GetXllCorner();
	double dMinY = header->GetYllCorner();
	double dMaxX = dMinX + (dResX * nGridCols);
	double dMaxY = dMinY + (dResX * nGridRows);
	if (header) delete header;

	
	//
	// allocate the memory block for the OUTPUT data
	//
	char qqq[128];
	double **ppData = new double * [ nRows ];
	if ( NULL == ppData )
	{
		sprintf( qqq, "Cannot allocate ppData <%d> bytes", int(nRows * sizeof(double)) );
		if (!DoBatch)
			Message(qqq, "INFO");
		return( NULL );
	}

	for ( int i=0; i<nRows; i++ ) 
	{
		ppData[i] = new double [ nCols ];
		if ( NULL == ppData[i] )
		{					
			sprintf( qqq, "Cannot allocate ppData[%d] <%d> bytes", i, int(nCols * sizeof(double)) );
			if (!DoBatch)
				Message(qqq, "INFO");			
			return( NULL );
		}
	}


	float pGridData;
	//
	// extract the data
	//
	int nCurrentCol = 0;
	char *seps = ",\n";
	if ( fUseEuclidean )
	{
		// get the file path
		GetProfileString( "TRANSPREDS", "EuclXTran", buff, lpParams );
		sprintf(buff, "%s", gmPath.ChangeExtension( buff, ".flt" ) );

		//
		// open the file for reading
		//
		BinaryFileClass *bfcIn = new BinaryFileClass(buff, BFC_ReadOnly);
		if ( !bfcIn->IsValid() )
		{
			if (!DoBatch)
				Message( "Cannot open XTran file for READ", "ExtractDataFromTransformGrids" );
			return(false);
		}	


		// extract grid data for points and copy to memory block
		rewind( fp );
		fgets( buff, BUFFLEN, fp ); // get the header
		int nThisRow = 0;
		while( 1 )
		{
			if ( NULL == fgets( buff, BUFFLEN, fp ) )
				break;

			// skip past the class_id
			char *p = strtok( buff, seps );		

			p = strtok( NULL, seps );
			double dX = atof( p );	     // get the X

			p = strtok( NULL, seps );
			double dY = atof( p );       // get the Y

			int nOffset = int( floor( ( ( dX - dMinX ) / dResX ) ) ) + ( int( floor( ( dMaxY - dY ) / dResX ) ) * nGridCols );
			
			bfcIn->SeekTo(nOffset * sizeof( float ));
			bfcIn->ReadFloat(&pGridData, 1);
			ppData[nThisRow][nCurrentCol] = pGridData;
			++nThisRow;
		}

		// close the Euclidean X Grid
		bfcIn->Close();

		// increment column index in memory block
		++nCurrentCol;


		// get the file path
		GetProfileString( "TRANSPREDS", "EuclYTran", buff, lpParams );
		sprintf(buff, "%s", gmPath.ChangeExtension( buff, ".flt" ) );

		//
		// open the file for reading
		//
		bfcIn = new BinaryFileClass(buff, BFC_ReadOnly);
		if ( !bfcIn->IsValid() )
		{
			if (!DoBatch)
				Message( "Cannot open YTran file for READ", "ExtractDataFromTransformGrids" );
			return(false);
		}

		// extract grid data for points and copy to memory block
		rewind( fp );
		fgets( buff, BUFFLEN, fp ); // get the header
		nThisRow = 0;
		while( 1 )
		{
			if ( NULL == fgets( buff, BUFFLEN, fp ) )
				break;

			// skip past the class_id
			char *p = strtok( buff, seps );

			p = strtok( NULL, seps );
			double dX = atof( p );	     // get the X

			p = strtok( NULL, seps );
			double dY = atof( p );       // get the Y

			int nOffset = int( floor( ( ( dX - dMinX ) / dResX ) ) ) + ( int( floor( ( dMaxY - dY ) / dResX ) ) * nGridCols );

			bfcIn->SeekTo(nOffset * sizeof( float ));
			bfcIn->ReadFloat(&pGridData, 1);
			ppData[nThisRow][nCurrentCol] = pGridData;
			++nThisRow;
		}

		// close the Euclidean Y Grid
		bfcIn->Close();

		// increment column index in memory block
		++nCurrentCol;

	}


	//
	// now do the environmental transforms
	//
	for ( int i=1; i<=nPreds; i++ )
	{
		sprintf( lpKey, "PredTran%d", i );
		GetProfileString( "TRANSPREDS", lpKey, lpPredPath, lpParams );
		
		if ( strlen( lpPredPath ) > 0 )
		{
			sprintf(lpPredPath, "%s", gmPath.ChangeExtension( lpPredPath, ".flt" ) );

			BinaryFileClass *bfcIn = new BinaryFileClass(lpPredPath, BFC_ReadOnly);
			if ( !bfcIn->IsValid() )
			{
				if (!DoBatch)
					Message( "Cannot open PredTran file for READ", "ExtractDataFromTransformGrids" );
				return(false);
			}

			// extract grid data for points and copy to memory block
			rewind( fp );
			fgets( buff, BUFFLEN, fp ); // get the header
			int nThisRow = 0;
			while( 1 )
			{
				if ( NULL == fgets( buff, BUFFLEN, fp ) )
					break;

				// skip past the class_id
				char *p = strtok( buff, seps );

				p = strtok( NULL, seps );
				double dX = atof( p );	     // get the X

				p = strtok( NULL, seps );
				double dY = atof( p );       // get the Y

				int nOffset = int( floor( ( ( dX - dMinX ) / dResX ) ) ) + ( int( floor( ( dMaxY - dY ) / dResX ) ) * nGridCols );

				bfcIn->SeekTo(nOffset * sizeof( float ));
				bfcIn->ReadFloat(&pGridData, 1);
				ppData[nThisRow][nCurrentCol] = pGridData;
				++nThisRow;
			}

			bfcIn->Close();

			// increment column index in memory block
			++nCurrentCol;
			
		}
	}
	if ( fp ) fclose( fp );

	//
	// set the matrix metrics into the parameter pointers
	//
	*pRows = nRows;
	*pCols = nCols;

	return( ppData );
}


//
// Dump contents of Demand Point Sample file
//
void DumpDemandPointFile( char *lpParams, char *lpSamplePath, 
	                      double *pSampleX, double *pSampleY, 
				          double **ppData, int nRows, int nCols )
{
	gmpath gmPath;
	FILE *fp = fopen(lpSamplePath, "w+t" );
	if (NULL == fp)
	{
		Message("Cannot create lpSamplePath", "DumpANNFile");
		return;
	}

	// write header
	if (1 == GetProfileInt( "GDMODEL", "UseEuclidean", lpParams ) )
	{
		fprintf(fp, "X,Y,XTran,YTran");
	}
	else
	{
		fprintf(fp, "X,Y");
	}

	int nPreds = GetProfileInt( "PREDICTORS", "NumPredictors", lpParams );
	char lpKey[64];
	char lpPredPath[BUFFLEN];
	for (int i=1; i<=nPreds; i++ )
	{
		sprintf( lpKey, "PredTran%d", i );
		GetProfileString( "TRANSPREDS", lpKey, lpPredPath, lpParams );

		if ( strlen( lpPredPath ) > 0 )
		{
			sprintf(lpPredPath, "%s", gmPath.ChangeExtension( lpPredPath, ".flt" ) );
			fprintf(fp, ",%s", gmPath.GetName(lpPredPath));
		}

		if (i == nPreds)
			fprintf(fp, "\n");
	}

	// now write the transformed predictor fields
	for ( int i=0; i<nRows; i++ )
	{
		fprintf( fp, "%lf,%lf,",  pSampleX[i], pSampleY[i] );
		
		for ( int j=0; j<nCols; j++ )
		{
			if ( j < (nCols-1) )
				fprintf( fp, "%lf,", ppData[i][j] );

			else
				fprintf( fp, "%lf\n", ppData[i][j] );
		}
	}
	fflush(fp);
	if ( fp ) fclose( fp );
}



//
// Dump contents of classification training file for debug
//
// NOTE: This comma delimited text file has NO HEADER
//       It just contains the class column followed 
//       by the data values for each predictor
//       for each sample point in the [CLASSIFY] PointFile
//
void DumpANNFile( char *lpParams, int *pClusterID, 
	              char *lpSamplePath, double *pSampleX, double *pSampleY, 
				  double **ppData, int nRows, int nCols )
{
	gmpath gmPath;
	FILE *fp = fopen(lpSamplePath, "w+t" );
	if (NULL == fp)
	{
		Message("Cannot create lpSamplePath", "DumpANNFile");
		return;
	}

	// write header
	if (1 == GetProfileInt( "GDMODEL", "UseEuclidean", lpParams ) )
	{
		fprintf(fp, "Class,X,Y,XTran,YTran");
	}
	else
	{
		fprintf(fp, "Class,X,Y");
	}

	int nPreds = GetProfileInt( "PREDICTORS", "NumPredictors", lpParams );
	char lpKey[64];
	char lpPredPath[BUFFLEN];
	for (int i=1; i<=nPreds; i++ )
	{
		sprintf( lpKey, "PredTran%d", i );
		GetProfileString( "TRANSPREDS", lpKey, lpPredPath, lpParams );

		if ( strlen( lpPredPath ) > 0 )
		{
			sprintf(lpPredPath, "%s", gmPath.ChangeExtension( lpPredPath, ".flt" ) );
			fprintf(fp, ",%s", gmPath.GetName(lpPredPath));
		}

		if (i == nPreds)
			fprintf(fp, "\n");
	}

	// now write the transformed predictor fields
	for ( int i=0; i<nRows; i++ )
	{
		fprintf( fp, "%d,%lf,%lf,",  pClusterID[i], pSampleX[i], pSampleY[i] );
		
		for ( int j=0; j<nCols; j++ )
		{
			if ( j < (nCols-1) )
				fprintf( fp, "%lf,", ppData[i][j] );

			else
				fprintf( fp, "%lf\n", ppData[i][j] );
		}
	}
	fflush(fp);
	if ( fp ) fclose( fp );
}


//
// Extract sample point coordinates to and X and Y vector
//
void GetSampleCoords(char *lpSamplePath, double *pSampleX, double *pSampleY, int nRows)

{
	char buff[BUFFLEN];
	FILE *fp = fopen( lpSamplePath, "r+t" );
	if (NULL == fp)
	{
		Message("Cannot open lpSamplePath", "GetSampleCoords");
		return;
	}
	

	// get the header
	fgets( buff, BUFFLEN, fp );
	char *seps = ",\n";
	for ( int i=0; i<nRows; i++ )
	{
		fgets( buff, BUFFLEN, fp );

		char *p = strtok( buff, seps );  // get the ID field

		p = strtok( NULL, seps );
		pSampleX[i] = atof( p );	     // get the X

		p = strtok( NULL, seps );
		pSampleY[i] = atof( p );         // get the Y
	}
	fclose(fp);
}


//
// Do a hierachical classification of a matrix of doubles
// and return a vector of ints representing the one-based
// classes derived from clustering the matrix data to a
// user defined number of classes derived from the parameter file.
//
int *GetHierachicalClusterClasses( double **ppData, int nRows, int nCols, int nClasses, bool DoBatch, FPTR fptr )
{
	//
	// This section creates a mask matrix like the 
	// ppData Matrix of ONES for the classification code
	//
	int **ppMask = new int * [ nRows ];
	for ( int i=0; i<nRows; i++ ) 
	{
		ppMask[i] = new int [ nCols ];
		for ( int j=0; j<nCols; j++ ) ppMask[i][j] = 1;
	}
	
	// initialise a weights vector of ones...
	double *pWeights = new double [ nCols ];
	for ( int i=0; i<nCols; i++) pWeights[i] = 1.0;

	fptr("Doing Hierachical Clustering...", 35);
	
	// hierarchical clustering here
	for ( int i=0; i<nCols; i++) pWeights[i] = 1.0;
	int nnodes = nRows-1;

	fptr("Doing Hierachical Clustering...", 45);	

	Node *pTree = treecluster( nRows, nCols, ppData, ppMask, pWeights, 0, 'b', 'a', 0, fptr ); 
	if (NULL == pTree )
	{ 
		// Indication that the treecluster routine failed 
		if (!DoBatch)
			Message( "treecluster routine failed due to insufficient memory", "ExtractDataFromTransformGrids" );
		if ( pWeights ) delete[] pWeights;
		for ( int i=0; i<nRows; i++ ) 
		{
			if ( ppMask[i] ) delete[] ppMask[i];
		}
		if ( ppMask ) delete[] ppMask;
		return( NULL );
	}


	/*FILE *fpDB = fopen("tree01.csv", "w+t");
	fprintf(fpDB, "Element,Left,Right,Distance\n");
	for ( int i=0; i<nRows-1; i++ )
	{
		fprintf(fpDB, "%d,%d,%d,%lf\n", i+1, pTree[i].left, pTree[i].right, pTree[i].distance);
	}
	fclose(fpDB);*/

	fptr("Doing Hierachical Clustering...", 70);	

	// get the desired number of classes
	int *pClusterID = new int [ nRows ];
	cuttree ( nRows, pTree, nClasses, pClusterID );


	/*FILE *fpCT = fopen("treeCut01.csv", "w+t");
	fprintf(fpCT, "Row,Class\n");
	for ( int i=0; i<nRows; i++ )
	{
		fprintf(fpCT, "%d,%d\n", i+1, pClusterID[i]+1);
	}
	fclose(fpCT);*/



	//
	// The Cluster ID vector is currently ZERO-BASED,
	// therefore we need to increment all the values
	// by ONE to convert into ONE_BASED ID values
	// for the grid classification that follows.
	//
	for ( int i=0; i<nRows; i++ ) pClusterID[i] += 1;

	//
	// clean up
	//
	if ( pTree ) delete[] pTree;
	if ( pWeights ) delete[] pWeights;
	for ( int i=0; i<nRows; i++ ) if ( ppMask[i] ) delete[] ppMask[i];
	if ( ppMask ) delete[] ppMask;

	return( pClusterID );
}



//
// Calculate the relative RGB's based on the distance between classes
//
bool DoClassificationColoring( char *lpParams, char *lpSamplePath, int nClasses, char *GridName, bool DoBatch, FPTR fptr )
{
	int nCurrent = 10;
	fptr("Doing Classification Coloring...", nCurrent);	

	//
	// get ANNFile.csv into memory
	//
	double **ppData = NULL;
	int *pClass = NULL;
	int nRows, nCols;
	if ( false == PopulateClassColor( lpSamplePath, &ppData, &pClass, &nRows, &nCols ) )
	{
		if (!DoBatch)
			Message( "PopulateClassColor failed", "DoClassificationColoring()" );
		return( false );
	}

	//
	// setup an nClasses X nCols matrix
	//
	int *pClassCounts = new int [ nClasses ];
	for ( int i=0; i<nClasses; i++ ) pClassCounts[i] = 0;

	double **ppClassAverages = new double * [ nClasses ];
	for ( int i=0; i<nClasses; i++ )
	{
		ppClassAverages[i] = new double [ nCols ];
		for ( int j=0; j<nCols; j++ ) ppClassAverages[i][j] = 0.0;
	}


	//
	// Calculate the class averages 
	//
	// ( Note that they are stored as ONE_BASED
	// so we subtract one for the ZERO_BASED index )
	//
	for ( int i=0; i<nRows; i++ )
	{
		// get class index
		int nThis = pClass[i] - 1;

		// increment count for the current class
		pClassCounts[ nThis ] += 1;

		// increment the transformed grid values for this class
		for ( int j=0; j<nCols; j++ )
		{
			ppClassAverages[ nThis ][j] += ppData[i][j];		
		}
	}
	// finally calculate average
	for ( int i=0; i<nClasses; i++ )
	{
		for ( int j=0; j<nCols; j++ )
		{
			ppClassAverages[i][j] /= pClassCounts[i];
		}
	}

	nCurrent = 20;
	fptr("Doing Classification Coloring...", nCurrent);	


	//
	// clean up the ANNFile memory allocations
	//
	if ( pClass ) delete[] pClass;
	for ( int i=0; i<nRows; i++ ) if ( ppData[i] ) delete ppData[i];
	if ( ppData ) delete[] ppData;
	// and the ClassCounts
	if ( pClassCounts ) delete[] pClassCounts;

	
	//
	// create a Manhattan metric Similarity Matrix
	//
	double **ppSimilarMatrix = new double * [ nClasses ];
	// initialise to zero to get the diagonals set to 0.0
	for ( int i=0; i<nClasses; i++ ) 
	{
		ppSimilarMatrix[i] = new double [ nClasses ];

		for ( int j=0; j<nClasses; j++ )  ppSimilarMatrix[i][j] = 0.0;
	}


	//
	// calculate the manhattan distances between 
	// the rows in the class averages matrix
	//
	double dTestMin = 0.0;
	for ( int i=1; i<nClasses; i++ )
	{
		for ( int j=0; j<i; j++ )
		{
			double dDist = CalcManhattanDistance( ppClassAverages, nCols, i, j );

			// apply the transformation to create a similarity metric
			// ( after Tanaka Tarumi 1995, P188 )
			ppSimilarMatrix[i][j] = ppSimilarMatrix[j][i] = -( dDist * dDist ) / 2;

			// adjust minimum value
			if ( ppSimilarMatrix[i][j] < dTestMin ) dTestMin = ppSimilarMatrix[i][j];
		}
	}


	//
	// Write a comma delimited Similarity Matrix Table
	//
	gmpath myGmPath;
	//char mySimPath[FILEPATH_BUFFSIZE];
	//sprintf(mySimPath, "%s\%s_SIMILARITY.csv", myGmPath.GetDirectoryPath(lpSamplePath), myGmPath.GetName(lpSamplePath));
	//Message(mySimPath, "myClassPath");
	//Message(GridName, "GridName");



	//char lpSimPath[FILEPATH_BUFFSIZE];
	////sprintf( lpSimPath, "%s_SIMILARITY.csv", myClassPath );
	////sprintf( lpSimPath, "%s_SIMILARITY.csv", GridName );
	//sprintf(lpSimPath, "%s\\%s_BFT_SIMILARITY.csv", myGmPath.GetDirectoryPath(lpSamplePath), myGmPath.GetName(lpSamplePath));
	//FILE *fpSim = fopen(lpSimPath, "w+t");
	//for ( int i=0; i<nClasses; i++ )
	//{
	//	for ( int j=0; j<nClasses; j++ )
	//	{
	//		double dTempSim = 1.0 - (fabs(ppSimilarMatrix[i][j] / dTestMin));

	//		fprintf(fpSim, "%lf", dTempSim);
	//		//fprintf(fpSim, "%lf", ppSimilarMatrix[i][j]);
	//		if ( j < nClasses-1 )
	//			fprintf(fpSim, "," );
	//		else
	//			fprintf(fpSim, "\n" );
	//	}
	//}
	//fclose(fpSim);


	//
	// clean up the class averages matrix 
	//
	for ( int i=0; i<nClasses; i++ )
		if ( ppClassAverages[i] ) delete[] ppClassAverages[i];
	if ( ppClassAverages ) delete[] ppClassAverages;


	nCurrent = 40;
	fptr("Doing Classification Coloring...", nCurrent);	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Do the Principal Coordinate Analysis of the similarity matrix
	// and use the first three coordinate axes as RGB color indices
	//
	matrix sim = new_matrix( nClasses, nClasses );
	for ( int i=0; i<nClasses; i++ )
	{
		for ( int j=0; j<nClasses; j++ )
		{
			sim[i][j] = ppSimilarMatrix[i][j];
		}
	}


	// clean up similarity matrix
	for ( int i=0; i<nClasses; i++ ) if ( ppSimilarMatrix[i] ) delete[] ppSimilarMatrix[i];
	if ( ppSimilarMatrix ) delete[] ppSimilarMatrix;


	// prepare for double centering
	vector l_mean = new_vector( nClasses );
	init_vector( l_mean, nClasses );
	double t_mean = 0.0;
	for ( int i=0; i<nClasses; i++ )
	{
		for ( int j=0; j<nClasses; j++ )
		{
			l_mean[i] += sim[i][j];
		}
		t_mean += l_mean[i];
		l_mean[i] /= nClasses;
	}
	t_mean /= nClasses * nClasses;


	// double centering
	for ( int i=0; i<nClasses; i++ )
	{
		for ( int j=0; j<nClasses; j++ )
		{
			sim[i][j] = sim[i][j] - l_mean[i] - l_mean[j] + t_mean;
		}
	}


	// principal coordinate analysis calcs
	matrix evec = new_matrix( nClasses, nClasses );
	for ( int i=0; i<nClasses; i++ )
	{
		for ( int j=0; j<nClasses; j++ )
		{
			evec[i][j] = sim[i][j];
		}
	} 


	matrix evec_work = new_matrix( nClasses, nClasses );
	vector eval = new_vector( nClasses );
	vector eval_work = new_vector( nClasses );
	if ( eigen( nClasses, evec, eval, eval_work ) )
	{
		for ( int i=0; i<nClasses; i++ )
		{
			for ( int j=0; j<nClasses; j++ )
			{
				evec[i][j] = sim[i][j];
			}
		} 	


		if ( jacobi( nClasses, evec, eval, evec_work ) )
		{
			free_vector( eval );
			free_vector( eval_work );
			free_matrix( evec );
			free_matrix( sim );
			free_matrix( evec_work );
			if (!DoBatch)
				Message( "error in jacobi", "DoClassificationColoring()" );
			return( false );
		}
	}

		
	nCurrent = 60;
	fptr("Doing Classification Coloring...", nCurrent);	

	// do the eigen values
	vector cont = new_vector( nClasses );
	vector cumcont = new_vector( nClasses );
	double dSum1 = 0.0;
	int pd = 1;
	for ( int i=0; i<nClasses; i++ )
	{
		if ( eval[i] > 0.0 )
			dSum1 += eval[i];

		if ( eval[i] < - EPSILON )
		{
			pd = 0;
			break;
		}
	}
	double dSum2 = 0.0;
	for ( int i=0; i<nClasses; i++ )
	{
		if ( eval[i] > 0.0 )
			dSum2 += eval[i];

		if ( pd )
		{
			cont[i] = eval[i] / dSum1;
			cumcont[i] = dSum2 / dSum1;
		}
		else
		{
			cont[i] = 0.0;
			cumcont[i] = 0.0;
		}
	}

	// do the eigen vector
	for ( int i=0; i<nClasses; i++ )
	{
		for ( int j=0; j<nClasses; j++ )
		{
			evec_work[i][j] = evec[j][i];
		}
	}

	// get the PCO scores
	for ( int i=0; i<nClasses; i++ )
	{
		for ( int j=0; j<nClasses; j++ )
		{
			if ( eval[j] > 0.0 )
			{
				evec_work[i][j] = evec_work[i][j] * sqrt( eval[j] );
			}
			else
			{
				evec_work[i][j] = 0.0;
			}			
		}
	}


	nCurrent = 70;
	fptr("Doing Classification Coloring...", nCurrent);	

	// get some min and max stats from the resulting PCO scores
	double dMin = DBL_MAX;
	double dMax = DBL_MIN;
	for ( INT i=0; i<nClasses; i++ )
	{
		for ( int j=0; j<3; j++ )
		{
			if ( evec_work[i][j] < dMin )
			{
				dMin = evec_work[i][j];
			}
			if ( evec_work[i][j] > dMax )
			{
				dMax = evec_work[i][j];
			}
		}
	}

	
	double dRange = dMax - dMin;

	nCurrent = 90;
	fptr("Doing Classification Coloring...", nCurrent);	

	//
	// setup an RGB matrix to convert the PCO 
	// scores into RGB values for legend coloring
	//
	int **ppRGB = new int * [ nClasses ];
	for ( int i=0; i<nClasses; i++ ) ppRGB[i] = new int [ 3 ];
	for ( int i=0; i<nClasses; i++ )
	{
		for ( int j=0; j<3; j++ )
		{
			double dCol = ( ( evec_work[i][j] - dMin ) / dRange ) * 255.0;

			ppRGB[i][j] = int( floor( dCol ) );
		}
	}

	
	//
	// create the RGB legends
	//
	gmpath GmPath;
	char ClassPath[FILEPATH_BUFFSIZE];
	sprintf(ClassPath, "%s\\%s", GmPath.GetDirectoryPath(lpSamplePath), GmPath.GetName(lpSamplePath));
	//Message(ClassPath, "ClassPath");

	nCurrent = 0;
	fptr("Doing Classification Coloring...", nCurrent);	
	for ( int i=1; i<=6; i++ )
	{
		if ( i * 100 / 6 > nCurrent )
		{
			nCurrent = i * 100 / 6;
			fptr("Doing Classification Coloring...", nCurrent);	
		}

		CreateRGBGridLegends( ClassPath, ppRGB, i, nClasses, GridName );
	}
	

	for ( int i=0; i<nClasses; i++ ) if ( ppRGB[i] ) delete[] ppRGB[i];
	if ( ppRGB ) delete[] ppRGB;
	free_vector( l_mean );
	free_vector( eval );
	free_vector( eval_work );
	free_matrix( evec );
	free_matrix( sim );
	free_matrix( evec_work );
	free_vector( cont );
	free_vector( cumcont );
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	fptr("Status: Ready.", 0);	
	return( true );
}


bool PopulateClassColor( char *lpSamplePath, double ***pppData, int **ppClass, int *pRows, int *pCols )
{
	//
	// get the ANNFile into an in-memory vector for the class column[0]
	// and an in-memory array for the transpose grids data columns
	//
	int nRows = nCountANNFileRows( lpSamplePath );
	int nCols = nCountANNFileColumns( lpSamplePath );
	// subtract three to exclude the X, Y, and class columns
	nCols -= 3;

	int *pClass = new int [ nRows ];

	double **ppData = new double * [ nRows ];
	for ( int i=0; i<nRows; i++ ) ppData[i] = new double [ nCols ];

	char *pRowBuff = new char [ TABLE_ROW_BUFFSIZE ];
	FILE *fp = fopen( lpSamplePath, "r+t" );

	// get header
	fgets(pRowBuff, TABLE_ROW_BUFFSIZE, fp);

	char *seps=",\n";
	for ( int i=0; i<nRows; i++ )
	{
		if ( NULL == fgets( pRowBuff, TABLE_ROW_BUFFSIZE, fp ) )
			break;

		// get the classID
		char *p = strtok( pRowBuff, seps );
		pClass[i] = atoi(p);

		// get the X Coord
		p = strtok( NULL, seps );

		// get the Y Coord
		p = strtok( NULL, seps );

		// get the transposed grid data
		for ( int j=0; j<nCols; j++ )
		{
			p = strtok( NULL, seps );
			ppData[i][j] = atof(p);
		}
	}

	fclose( fp );
	if ( pRowBuff )  delete[] pRowBuff;

	//
	// setup the pointers
	//
	*pppData = ppData;
	*ppClass = pClass;
	*pRows = nRows;
	*pCols = nCols;

	return( true );
}


double CalcManhattanDistance( double **ppClassAverages, int nCols, int x1, int x2 )
{
	double dDist = 0.0;

	for ( int i=0; i<nCols; i++ )
	{
		dDist += fabs( ppClassAverages[x1][i] - ppClassAverages[x2][i] );
	}

	return( dDist );
}


//
// rewite the .GRD file for the classification grid using 
// the legend RGBs created via PCO over the number  of classes
//
void CreateRGBGridLegends( char *ClassPath, int **ppRGB, int Version, int nClasses, char *GridName )
{
	int i;

	//char lpWorkspace[BUFFLEN];
	//GetProfileString( "GDMODEL", "WorkspacePath", lpWorkspace, lpParams );

	//
	// dump out the variants of the RGB.csv file to the Workdirectory
	//
	char lpColorPathCSV[FILEPATH_BUFFSIZE];
	char lpColorPathCLR[FILEPATH_BUFFSIZE];

	//
	// 1 of 6 (RGB)
	//
	sprintf( lpColorPathCSV, "%s_RGB.csv", ClassPath );
	sprintf( lpColorPathCLR, "%s_RGB.clr", ClassPath );
	//sprintf( lpColorPathCSV, "%s_RGB.csv", GridName );
	//sprintf( lpColorPathCLR, "%s_RGB.clr", GridName );

	FILE *fpCSV = fopen( lpColorPathCSV, "w+t" );
	FILE *fpCLR = fopen( lpColorPathCLR, "w+t" );

	fprintf(fpCSV, "Red,Green,Blue\n");
	for ( i=0; i<nClasses; i++ )
	{
		fprintf( fpCSV, "%d,%d,%d\n", ppRGB[i][0], ppRGB[i][1], ppRGB[i][2] );
		fprintf( fpCLR, "%d %d %d %d\n", i+1, ppRGB[i][0], ppRGB[i][1], ppRGB[i][2] );
	}

	fclose(fpCSV);
	fclose(fpCLR);


	//
	// 2 of 6 (RBG)
	//
	sprintf( lpColorPathCSV, "%s_RBG.csv", ClassPath );
	sprintf( lpColorPathCLR, "%s_RBG.clr", ClassPath );
	//sprintf( lpColorPathCSV, "%s_RBG.csv", GridName );
	//sprintf( lpColorPathCLR, "%s_RBG.clr", GridName );

	fpCSV = fopen( lpColorPathCSV, "w+t" );
	fpCLR = fopen( lpColorPathCLR, "w+t" );

	fprintf(fpCSV, "Red,Green,Blue\n");
	for ( i=0; i<nClasses; i++ )
	{
		fprintf( fpCSV, "%d,%d,%d\n", ppRGB[i][0], ppRGB[i][2], ppRGB[i][1] );
		fprintf( fpCLR, "%d %d %d %d\n", i+1, ppRGB[i][0], ppRGB[i][2], ppRGB[i][1] );
	}

	fclose(fpCSV);
	fclose(fpCLR);


	//
	// 3 of 6 (BRG)
	//
	sprintf( lpColorPathCSV, "%s_BRG.csv", ClassPath );
	sprintf( lpColorPathCLR, "%s_BRG.clr", ClassPath );
	//sprintf( lpColorPathCSV, "%s_BRG.csv", GridName );
	//sprintf( lpColorPathCLR, "%s_BRG.clr", GridName );

	fpCSV = fopen( lpColorPathCSV, "w+t" );
	fpCLR = fopen( lpColorPathCLR, "w+t" );

	fprintf(fpCSV, "Red,Green,Blue\n");
	for ( i=0; i<nClasses; i++ )
	{
		fprintf( fpCSV, "%d,%d,%d\n", ppRGB[i][2], ppRGB[i][0], ppRGB[i][1] );
		fprintf( fpCLR, "%d %d %d %d\n", i+1, ppRGB[i][2], ppRGB[i][0], ppRGB[i][1] );
	}

	fclose(fpCSV);
	fclose(fpCLR);


	//
	// 4 of 6 (BGR)
	//
	sprintf( lpColorPathCSV, "%s_BGR.csv", ClassPath );
	sprintf( lpColorPathCLR, "%s_BGR.clr", ClassPath );
	//sprintf( lpColorPathCSV, "%s_BGR.csv", GridName );
	//sprintf( lpColorPathCLR, "%s_BGR.clr", GridName );

	fpCSV = fopen( lpColorPathCSV, "w+t" );
	fpCLR = fopen( lpColorPathCLR, "w+t" );

	fprintf(fpCSV, "Red,Green,Blue\n");
	for ( i=0; i<nClasses; i++ )
	{
		fprintf( fpCSV, "%d,%d,%d\n", ppRGB[i][2], ppRGB[i][1], ppRGB[i][0] );
		fprintf( fpCLR, "%d %d %d %d\n", i+1, ppRGB[i][2], ppRGB[i][1], ppRGB[i][0] );
	}

	fclose(fpCSV);
	fclose(fpCLR);


	//
	// 5 of 6 (GRB)
	//
	sprintf( lpColorPathCSV, "%s_GRB.csv", ClassPath );
	sprintf( lpColorPathCLR, "%s_GRB.clr", ClassPath );
	//sprintf( lpColorPathCSV, "%s_GRB.csv", GridName );
	//sprintf( lpColorPathCLR, "%s_GRB.clr", GridName );

	fpCSV = fopen( lpColorPathCSV, "w+t" );
	fpCLR = fopen( lpColorPathCLR, "w+t" );

	fprintf(fpCSV, "Red,Green,Blue\n");
	for ( i=0; i<nClasses; i++ )
	{
		fprintf( fpCSV, "%d,%d,%d\n", ppRGB[i][1], ppRGB[i][0], ppRGB[i][2] );
		fprintf( fpCLR, "%d %d %d %d\n", i+1, ppRGB[i][1], ppRGB[i][0], ppRGB[i][2] );
	}

	fclose(fpCSV);
	fclose(fpCLR);


	//
	// 6 of 6 (GBR)
	//
	sprintf( lpColorPathCSV, "%s_GBR.csv", ClassPath );
	sprintf( lpColorPathCLR, "%s_GBR.clr", ClassPath );
	//sprintf( lpColorPathCSV, "%s_GBR.csv", GridName );
	//sprintf( lpColorPathCLR, "%s_GBR.clr", GridName );

	fpCSV = fopen( lpColorPathCSV, "w+t" );
	fpCLR = fopen( lpColorPathCLR, "w+t" );

	fprintf(fpCSV, "Red,Green,Blue\n");
	for ( i=0; i<nClasses; i++ )
	{
		fprintf( fpCSV, "%d,%d,%d\n", ppRGB[i][1], ppRGB[i][2], ppRGB[i][0] );
		fprintf( fpCLR, "%d %d %d %d\n", i+1, ppRGB[i][1], ppRGB[i][2], ppRGB[i][0] );
	}

	fclose(fpCSV);
	fclose(fpCLR);
}


//
// return the number of effective rows in a comma delimited text file 
//
int nCountANNFileRows( char *p )
{
	int n = 0;

	char *pRow = new char [ TABLE_ROW_BUFFSIZE ];
	if ( pRow )
	{
		FILE *fp = fopen( p, "r+t" );
		if ( NULL == fp )
		{
			Message( "Cannot open text file argument", "nCountFileRowsNoHeader" );
			delete[] pRow;
			return( n );
		}

		// get header
		fgets( pRow, TABLE_ROW_BUFFSIZE, fp );
		
		while( 1 )
		{
			if ( NULL == fgets( pRow, TABLE_ROW_BUFFSIZE, fp ) )
				break;

			++n;
		}
		fclose( fp );
		delete[] pRow;
	}

	return( n );
}



//
// return the number of effective columns in a comma delimited text file 
//         This function assumes that there is AT LEAST ONE ROW
//
int nCountANNFileColumns( char *p )
{
	int n = 0;

	char *pRow = new char [ TABLE_ROW_BUFFSIZE ];
	if ( pRow )
	{
		FILE *fp = fopen( p, "r+t" );
		if ( NULL == fp )
		{
			Message( "Cannot open text file argument", "nCountFileColumns" );
			delete[] pRow;
			return( n );
		}
		
		char *seps = ",\n";

		// get the first row
		fgets( pRow, TABLE_ROW_BUFFSIZE, fp );

		// get the first comma delimited token
		char *t = strtok( pRow, seps );
		n = 1;

		while( 1 )
		{
			t = strtok( NULL, seps );

			// are we done?
			if ( NULL == t ) break;

			// else...
			++n;
		}
		fclose( fp );
		delete[] pRow;
	}

	return( n );
}



//
// use the class,tran1,..tranN table and the 
// GDM transformed grids to produce a single 
// classification grid using the ANN tree. 
//
bool DoGridClassification( char *lpParams, 
						   char *lpSubDomain,
						   int *pClusterID, double **ppData, 
						   int nRows, int nCols, 
						   char *lpOutPath,
						   bool DoBatch,
						   FPTR fptr )
{
	//
	// construct the ANN tree
	//
	ANNpointArray data_pnts = annAllocPts( nRows, nCols );
	// populate from ppData
	for ( int i=0; i<nRows; i++ )
	{
		for ( int j=0; j<nCols; j++ )
		{
			data_pnts[i][j] = ppData[i][j];
		}
	}


	//
	// now build the kd tree search structure with 
	// the data points, number of data points and dimension of search space
	//
	ANNkd_tree *theTree = new ANNkd_tree ( data_pnts, nRows, nCols ); 


	//
	// find the nearest neighbour classes
	//
	// 
	// setup some structures for the queries to follow
	ANNidxArray nn_idx = new ANNidx[1];
	ANNdistArray dist = new ANNdist[1];
	ANNpoint query_pt = annAllocPt( nCols );


	//
	// Extract the filepaths for the transform grids in the parameter file
	//
	bool fHaveEuclidean = GetProfileInt( "GDMODEL", "UseEuclidean", lpParams ) == 1 ? true : false;
	int nGrids = 0;
	if ( fHaveEuclidean ) nGrids += 2;

	int nEnvGrids = GetProfileInt( "PREDICTORS", "NumPredictors", lpParams );
	char lpKey[64];
	for ( int i=1; i<=nEnvGrids; i++ )
	{
		//
		// we ONLY transform grids that have NON-ZERO coefficients
		//
		char lpPredPath[256];
		sprintf( lpKey, "PredTran%d", i );
		GetProfileString( "TRANSPREDS", lpKey, lpPredPath, lpParams );
		if ( strlen( lpPredPath ) > 0 )      //if ( fPredictorHasNonZeroCoeffs( lpParams, i ) )
		{
			++nGrids;
		}
	}


	//
	// create an array for the file handles to the transform grids
	//
	int *pHandles = new int [ nGrids ];

	// open the grid binaries...
	char pGridPath[BUFFLEN];
	gmpath gmPath;
	if ( fHaveEuclidean )
	{
		GetProfileString( "TRANSPREDS", "EuclXTran", pGridPath, lpParams );
		sprintf(pGridPath, "%s", gmPath.ChangeExtension(pGridPath, ".flt"));

		pHandles[0] = _open( pGridPath, _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE  );
		if ( pHandles[0] < 0 )
		{
			if (!DoBatch)
				Message( "Cannot open EuclXTran.flt file", "DoGridClassification" );
			if ( pHandles ) delete[] pHandles;
			return( false );
		}


		GetProfileString( "TRANSPREDS", "EuclYTran", pGridPath, lpParams );
		sprintf(pGridPath, "%s", gmPath.ChangeExtension(pGridPath, ".flt"));

		pHandles[1] = _open( pGridPath, _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE  );
		if ( pHandles[1] < 0 )
		{
			if (!DoBatch)
				Message( "Cannot open EuclYTran.flt file", "DoGridClassification" );
			if ( pHandles ) delete[] pHandles;
			return( false );
		}

		// now do the environmental grids
		int nThis = 2; // starting slot for the environmental grid handles

		for ( int i=1; i<=nEnvGrids; i++ )
		{
			//
			// we ONLY transform grids that have NON-ZERO coefficients
			//
			char lpKey[64];
			char lpPredPath[256];
			sprintf( lpKey, "PredTran%d", i );
			GetProfileString( "TRANSPREDS", lpKey, lpPredPath, lpParams );
			if ( strlen( lpPredPath ) > 0 )  //if ( fPredictorHasNonZeroCoeffs( lpParams, i ) )
			{
				sprintf( lpKey, "PredTran%d", i );
				GetProfileString( "TRANSPREDS", lpKey, pGridPath, lpParams );
				sprintf(pGridPath, "%s", gmPath.ChangeExtension(pGridPath, ".flt"));

				pHandles[nThis] = _open( pGridPath, _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE  );
				if ( pHandles[nThis] < 0 )
				{
					if (!DoBatch)
						Message( "Cannot open EnvTransform.GRI file", "DoGridClassification" );
					for ( int j=0; j<i; j++ ) _close( pHandles[i] );
					if ( pHandles ) delete[] pHandles;
					return( false );
				}
				
				++nThis;
			}
		}
	}
	else
	{
		int nThis = 0; // starting slot for the environmental grid handles

		for ( int i=1; i<=nEnvGrids; i++ )
		{
			//
			// we ONLY transform grids that have NON-ZERO coefficients
			//
			char lpKey[64];
			char lpPredPath[256];
			sprintf( lpKey, "PredTran%d", i );
			GetProfileString( "TRANSPREDS", lpKey, lpPredPath, lpParams );
			if ( strlen( lpPredPath ) > 0 ) 
			{
				sprintf( lpKey, "PredTran%d", i );
				GetProfileString( "TRANSPREDS", lpKey, pGridPath, lpParams );
				sprintf(pGridPath, "%s", gmPath.ChangeExtension(pGridPath, ".flt"));

				pHandles[nThis] = _open( pGridPath, _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE  );
				if ( pHandles[nThis] < 0 )
				{
					if (!DoBatch)
						Message( "Cannot open EnvTransform.GRI file", "DoGridClassification" );
					for ( int j=0; j<i; j++ ) _close( pHandles[i] );
					if ( pHandles ) delete[] pHandles;
					return( false );
				}
				
				++nThis;
			}
		}
	}

	
	//
	// Open domain grid
	//
	char lpDomainPath[BUFFLEN];
	strcpy(lpDomainPath,lpSubDomain);
	sprintf(lpDomainPath, "%s", gmPath.ChangeExtension(lpDomainPath, ".flt"));
	int hDomain = _open( lpDomainPath, _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE  );

	//
	// get some grid metrics from the domain grid
	//	
	EsriBinaryHeader *DomainHeader = new EsriBinaryHeader(gmPath.ChangeExtension(lpDomainPath, ".hdr"));
	int nDomainGridRows = DomainHeader->GetNumRows();
	int nDomainGridCols = DomainHeader->GetNumCols();
	float fDomainNoData = DomainHeader->GetNoDataValue();
	double dDomainResX = DomainHeader->GetCellSize();
	double dDomainMinX = DomainHeader->GetXllCorner();
	double dDomainMinY = DomainHeader->GetYllCorner();
	double dDomainMaxX = dDomainMinX + (dDomainResX * nDomainGridCols);
	double dDomainMaxY = dDomainMinY + (dDomainResX * nDomainGridRows);

	//
	// get some grid metrics from the last transform grid
	//	
	EsriBinaryHeader *TranHeader = new EsriBinaryHeader(gmPath.ChangeExtension(pGridPath, ".hdr"));
	int nTranGridRows = TranHeader->GetNumRows();
	int nTranGridCols = TranHeader->GetNumCols();
	float fTranNoData = TranHeader->GetNoDataValue();
	double dTranResX = TranHeader->GetCellSize();
	double dTranMinX = TranHeader->GetXllCorner();
	double dTranMinY = TranHeader->GetYllCorner();
	double dTranMaxX = dTranMinX + (dTranResX * nTranGridCols);
	double dTranMaxY = dTranMinY + (dTranResX * nTranGridRows);
	if (TranHeader) delete TranHeader;


	//
	// Calculate starting row and column indices from Domain Extent VS Transform Grid Extent
	//
	int nStartRow = (int)((dTranMaxY - dDomainMaxY) / dDomainResX);
	int nStartCol = (int)((dDomainMinX - dTranMinX) / dDomainResX);

	if ((nStartRow < 0) || (nStartCol < 0))
	{
		Message("Domain Extends Beyond Transform Grid Range, \ncheck the .hdr file against the transform grid .hdr files", "ERROR");
		_close( hDomain);
		if (DomainHeader) delete DomainHeader;
		for ( int i=0; i<nGrids; i++ )
		{
			_close( pHandles[i] );
		}
		return(false);
	}

	//Message(nStartRow, "nStartRow");
	//Message(nStartCol, "nStartCol");

	//
	// Create output grid
	//
	//char lpWorkspace[BUFFLEN];
	//char lpOutPath[BUFFLEN];
	//GetProfileString( "GDMODEL", "WorkspacePath", lpWorkspace, lpParams );
	//sprintf( lpOutPath, "%s\\%s", lpWorkspace, lpOutName );
	//sprintf(lpOutPath, "%s", gmPath.ChangeExtension(lpOutPath, ".flt"));
	//int hOut = _open( lpOutPath, _O_BINARY | _O_RDWR | _O_CREAT | _O_TRUNC, S_IREAD | S_IWRITE  );
	int hOut = _open( gmPath.ChangeExtension(lpOutPath, ".flt"), _O_BINARY | _O_RDWR | _O_CREAT | _O_TRUNC, S_IREAD | S_IWRITE  );
	if ( hOut < 0 )
	{
		if (!DoBatch)
			Message( "Cannot open Output grid", lpOutPath );
		return(false);
	}


	// 
	// Allocate row memory blocks...
	//
	float *pOutRow = new float [nDomainGridCols];
	float **ppInputRows = new float * [nGrids];
	for ( int i=0; i<nGrids; i++ ) ppInputRows[i] = new float [nTranGridCols];
	int *pSubDomainRow = new int [nDomainGridCols];


	//
	// Do Classification
	//
	int nCurrent = 0;
	fptr("Doing nearest neighbour classification...", nCurrent);	
	for ( int i=0; i<nTranGridRows; i++ )
	{
		if ( i * 100 / nTranGridRows > nCurrent )
		{
			nCurrent = i * 100 / nTranGridRows;
			fptr("Doing nearest neighbour classification...", nCurrent);	
		}

		// get the row data for each transform grid
		for ( int g=0; g<nGrids; g++ )
		{
			_read( pHandles[g], ppInputRows[g], nTranGridCols * sizeof( float ) );
		}

		if ( (i >= nStartRow) && (i < (nStartRow + nDomainGridRows)))
		{

			// read the sub-domain row
			_read( hDomain, pSubDomainRow, nDomainGridCols * sizeof( float ) );

			for ( int j=0; j<nTranGridCols; j++ )
			{
				if ( (j >= nStartCol) && (j < (nStartCol + nDomainGridCols)))
				{
					// do the nearest neighbour search and allocate the resulting class ID			
					if ( pSubDomainRow[j-nStartCol] > fDomainNoData )
					{
						// populate the query_pt vector
						for ( int g=0; g<nGrids; g++ )
						{
							query_pt[g] = ppInputRows[g][j];
						}
			
						// find the nearest neighbour class in the ANN tree
						theTree->annkSearch( query_pt, 1, nn_idx, dist, 0 );		
	
						// set the class in the output grid
						pOutRow[j-nStartCol] = (float)pClusterID[ nn_idx[0] ];
					}
	
					// we have NO_DATA, then write no-data 
					else
					{
						pOutRow[j-nStartCol] = fDomainNoData;
					}

				} // if ( (j >= nStartCol) && (j < (nStartCol + nDomainGridCols)))

			} // for ( int j=0; j<nTranGridCols; j++ )

			// write this output row
			_write( hOut, pOutRow, nDomainGridCols * sizeof( float ) );

		} // if ( (i >= nStartRow) && (i < (nStartRow + nDomainGridRows)))
	} // for ( int i=0; i<nTranGridRows; i++ )


	//
	// close the binary grid files
	//
	_close( hOut );
	_close( hDomain);
	for ( int i=0; i<nGrids; i++ )
	{
		_close( pHandles[i] );
	}
	

	//
	// clean up
	//
	for ( int i=0; i<nGrids; i++ ) if ( ppInputRows[i] ) delete[] ppInputRows[i];
	if ( ppInputRows ) delete[] ppInputRows;
	if (pSubDomainRow) delete[] pSubDomainRow;
	if (pOutRow) delete[] pOutRow;
	if ( pHandles ) delete[] pHandles;
	if (theTree) delete theTree;


	//
	// finally, rewrite the legend file for the classification grid
	//
	//DomainHeader->CopyTo(gmPath.ChangeExtension(lpOutPath, ".hdr"));
	//if (DomainHeader) delete DomainHeader;


	//
	// finally, rewrite the legend file for the classification grid using the first transform grid header
	//
	bool haveWrittenOutHeader = false;
	for (int i = 1; i <= nEnvGrids; i++)
	{
		//
		// we ONLY transform grids that have NON-ZERO coefficients
		//
		char lpKey[64];
		char lpPredPath[1024];
		sprintf(lpKey, "PredTran%d", i);
		GetProfileString("TRANSPREDS", lpKey, lpPredPath, lpParams);
		if (strlen(lpPredPath) > 0)
		{
			sprintf(lpKey, "PredTran%d", i);
			GetProfileString("TRANSPREDS", lpKey, pGridPath, lpParams);
			EsriBinaryHeader *outHeader = new EsriBinaryHeader(gmPath.ChangeExtension(pGridPath, ".hdr"));
			outHeader->CopyTo(gmPath.ChangeExtension(lpOutPath, ".hdr"));
			if (outHeader) delete outHeader;
			break;
		}
	}

	if (!haveWrittenOutHeader) // no transform grids??? shouldn't happen!!!
	{
		DomainHeader->CopyTo(gmPath.ChangeExtension(lpOutPath, ".hdr"));
		if (DomainHeader) delete DomainHeader;
	}

	return( true );
}


//
// Check that all the valid cells in the domain grid have valid data in the first GDM transform grid
//
bool DomainAndTransformOK(char *pDomain, char *pGDMParamPath)
{
	int nEnvGrids = GetProfileInt("PREDICTORS", "NumPredictors", pGDMParamPath);
	if (nEnvGrids < 1)
	{
		Message("NumEnvGrids < 1", "ERROR");
		return(false);
	}

	char lpKey[64];
	char lpPredPath[1024];
	memset(lpPredPath, 0, 1024);
	for (int i = 1; i <= nEnvGrids; i++)
	{
		//
		// we ONLY transform grids that have NON-ZERO coefficients
		//
		sprintf(lpKey, "PredTran%d", i);
		GetProfileString("TRANSPREDS", lpKey, lpPredPath, pGDMParamPath);
		if (strlen(lpPredPath) > 0)
		{
			break;
		}
	}
	if (NULL == lpPredPath)
	{
		Message("Have No Valid GDM Transform Paths", "ERROR");
		return(false);
	}

	//
	// Do a cell by cell comparison of the doamin grid to the transform grid
	// 
	// If any valid domain cell does NOT contain vlid data in the transform grid then return false as the grids are NOT coincident
	//





	return(true);
}



//
// Calculate the GDM Transform value
//
double CalculateGDMTransform( double dVal, int nSplines, double *pQuantiles, double *pCoeffs )
{
	double dCalc = 0.0;

	for ( int i=0; i<nSplines; i++ )
	{
		if ( i == 0 )    // first spline
		{
			dCalc += DoTranSplineCalc( dVal, pQuantiles[0], pQuantiles[0], pQuantiles[1] ) * pCoeffs[i];
		}

		else if ( i == nSplines-1 )  // last spline
		{
			dCalc += DoTranSplineCalc( dVal, pQuantiles[nSplines-2], pQuantiles[nSplines-1], pQuantiles[nSplines-1] ) * pCoeffs[i];
		}

		else // a middle spline
		{
			dCalc += DoTranSplineCalc( dVal, pQuantiles[i-1], pQuantiles[i], pQuantiles[i+1] ) * pCoeffs[i];
		}
	}
	return(dCalc);
}



//////
////// Test for Dan Faith's problen with clustering of demand points
//////
////bool DoDemandPoints_02(char *pParams, 
////	                   char *lpDomainPath, 
////	                   char *lpOutName, 
////					   char *lpSamplePath, 
////					   bool DoBatch, 
////					   FPTR fptr)
////{
////	//
////	// Extract a set of samples points for the nearest neighbour classification and write to
////	// a comma delimited text file in the Working directory called (Classification NAME)_Sample.csv
////	//
////	gmpath GmPath;
////	char namestring [FILEPATH_BUFFSIZE];
////	int NUMSAMPLES = 100000;
////	//int NUMSAMPLES = 1000;
////	sprintf(namestring, "%s\\%s_Sample.csv", GmPath.GetDirectoryPath(lpSamplePath), GmPath.GetName(lpSamplePath));
////	if (!CreateSamplePointMesh( pParams, namestring, lpDomainPath, NUMSAMPLES, false, fptr ))
////	{
////		if (!DoBatch)
////			Message("Cannot CreateSamplePointFile()", "DoGDMClassification");
////		return(false);
////	}
////
////	//
////	// Extract transform grid data at the sample point locations
////	// into a matrix to pass into the classification routines
////	//
////	int nRows = 0;
////	int nCols = 0;
////	double **ppData = ExtractDataFromTransformGrids( pParams, lpDomainPath, namestring, &nRows, &nCols, DoBatch );
////	if ( NULL == ppData )
////	{
////		if (!DoBatch)
////			Message( "ExtractDataFromTransformGrids failure", "DoClassification" );
////		return( false );
////	}
////
////
////	//
////	// Get the sample X,Y's
////	// 
////	double *pSampleX = new double [nRows];
////	double *pSampleY = new double [nRows];
////	GetSampleCoords(namestring, pSampleX, pSampleY, nRows);
////
////	//
////	// delete the points file now we don't need it
////	//
////	if (GmPath.FileExists(namestring))
////	{
////		//Message(namestring, "namestring");
////		remove(namestring);
////	}
////
////
////	//
////	// write the training file for debug and for use in deriving the colors
////	//
////	int *pClusterID = new int [nRows];
////	for ( int i=0; i<nRows; i++ ) pClusterID[i] = 1;
////	DumpDemandPointFile( pParams, lpSamplePath, pSampleX, pSampleY, ppData, nRows, nCols );
////	
////
////	//
////	// Cleanup
////	//
////	if (pSampleX) delete[] pSampleX;
////	if (pSampleY) delete[] pSampleY;
////	for ( int i=0; i<nRows; i++ )
////	{
////		if (ppData[i]) delete ppData[i];
////	}
////	if (ppData) delete[] ppData;
////	if (pClusterID) delete[] pClusterID;
////	return(true);
////}
////
////
////
////
////bool DoDemandPoints_01(char *pParams, 
////	                   char *lpDomainPath, 
////	                   char *lpOutName, 
////					   char *lpSamplePath, 
////					   bool DoBatch, 
////					   FPTR fptr)
////{
////	//
////	// Extract a set of samples points for the nearest neighbour classification and write to
////	// a comma delimited text file in the Working directory called (Classification NAME)_Sample.csv
////	//
////	gmpath GmPath;
////	char namestring [FILEPATH_BUFFSIZE];
////	//int NUMSAMPLES = 10000;
////	int NUMSAMPLES = 1000;
////	sprintf(namestring, "%s\\%s_Sample.csv", GmPath.GetDirectoryPath(lpSamplePath), GmPath.GetName(lpSamplePath));
////	if (!CreateSamplePointMesh( pParams, namestring, lpDomainPath, NUMSAMPLES, false, fptr ))
////	{
////		if (!DoBatch)
////			Message("Cannot CreateSamplePointFile()", "DoGDMClassification");
////		return(false);
////	}
////
////	//
////	// Extract transform grid data at the sample point locations
////	// into a matrix to pass into the classification routines
////	//
////	int nRows = 0;
////	int nCols = 0;
////	double **ppData = ExtractDataFromTransformGrids( pParams, lpDomainPath, namestring, &nRows, &nCols, DoBatch );
////	if ( NULL == ppData )
////	{
////		if (!DoBatch)
////			Message( "ExtractDataFromTransformGrids failure", "DoClassification" );
////		return( false );
////	}
////
////
////	//
////	// Get the sample X,Y's
////	// 
////	double *pSampleX = new double [nRows];
////	double *pSampleY = new double [nRows];
////	GetSampleCoords(namestring, pSampleX, pSampleY, nRows);
////
////	//
////	// delete the points file now we don't need it
////	//
////	if (GmPath.FileExists(namestring))
////	{
////		//Message(namestring, "namestring");
////		remove(namestring);
////	}
////
////
////	//
////	// Now write the table of sites,transformed predictor values...
////	//
////	sprintf(namestring, "%s\\%s_DEMAND.csv", GmPath.GetDirectoryPath(lpSamplePath), GmPath.GetName(lpSamplePath));
////	FILE *fpDist = fopen( namestring, "w+t" );
////	fprintf(fpDist, "X,Y");
////	int nPreds = GetProfileInt( "PREDICTORS", "NumPredictors", pParams );
////	char lpKey[64];
////	char lpPredPath[BUFFLEN];
////	for (int i=1; i<=nPreds; i++ )
////	{
////		sprintf( lpKey, "PredTran%d", i );
////		GetProfileString( "TRANSPREDS", lpKey, lpPredPath, pParams );
////
////		if ( strlen( lpPredPath ) > 0 )
////		{
////			sprintf(lpPredPath, "%s", GmPath.ChangeExtension( lpPredPath, ".flt" ));
////			fprintf(fpDist, ",%s", GmPath.GetName(lpPredPath));
////		}
////
////		if (i == nPreds)
////			fprintf(fpDist, "\n");
////	}
////	for ( int i=0; i<nRows; i++ )
////	{
////		fprintf(fpDist, "%lf,%lf", pSampleX[i], pSampleY[i]);
////		for ( int j=0; j<nCols; j++ )
////		{
////			fprintf(fpDist, ",%lf", ppData[i][j]);
////			if ( j == nCols-1 )
////				fprintf(fpDist, "\n");
////		}
////	}
////	fclose(fpDist); 
////
////
////	//
////	// Allocate the distance matrix
////	//
////	double **ppDistMatrix = new double * [nRows];
////	for ( int i=0; i<nRows; i++ )
////	{
////		ppDistMatrix[i] = new double [nRows];
////		for ( int j=0; j<nRows; j++ )
////		{
////			ppDistMatrix[i][j] = 0.0;
////		}
////	}
////
////
////	//
////	// Calculate the pairwise dissimilarities
////	//
////	for ( int i=0; i<nRows-1; i++ )
////	{
////		for (int j=i+1; j<nRows; j++ )
////		{
////			double dManhattan = 0.0;
////			for ( int k=0; k<nCols; k++ )
////			{
////				dManhattan += fabs(ppData[i][k] - ppData[j][k]);
////			}
////			ppDistMatrix[i][j] = ppDistMatrix[j][i] = dManhattan;
////		}
////	}
////
////
////	//
////	// Write the distance matrix
////	//
////	/*sprintf(namestring, "%s\\%s_DEMAND_DIST.csv", GmPath.GetDirectoryPath(lpSamplePath), GmPath.GetName(lpSamplePath));
////	FILE *fpDemandDist = fopen(namestring, "w+t");
////	for ( int i=0; i<nRows; i++ )
////	{
////		for ( int j=0; j<nRows; j++ )
////		{
////			fprintf( fpDemandDist, "%lf", ppDistMatrix[i][j] );
////			if ( j < nRows-1 )
////				fprintf(fpDemandDist, ",");
////			else
////				fprintf(fpDemandDist, "\n");
////		}
////	}
////	fclose(fpDemandDist);*/
////
////
////	//
////	// Allocate the Dist,X,Y Table and initialise to sample site coords
////	//
////	DemandSiteRecord *pDemandSites = new DemandSiteRecord [nRows];
////	for ( int i=0; i<nRows; i++ )
////	{
////		pDemandSites[i].dDist = 0.0;
////		pDemandSites[i].dX = pSampleX[i];
////		pDemandSites[i].dY = pSampleY[i];
////	}
////
////
////	// 
////	// calculate the proportion as a percentage of the sites with a dissimilarity greater than the threshold
////	//
////	double dThreshhold = 0.7;  // equates to a dissimilarity of 0.503414696
////	for ( int i=0; i<nRows; i++ )
////	{
////		int nCount = 0;
////		for ( int j=0; j<nRows; j++ )
////		{
////			if (ppDistMatrix[i][j] > dThreshhold)
////			{
////				++nCount;
////			}
////		}
////		pDemandSites[i].dDist = nCount * 100 / nRows;
////	}
////
////
////	//
////	// write results to file
////	//
////	sprintf(namestring, "%s\\%s_DEMAND_RECS.csv", GmPath.GetDirectoryPath(lpSamplePath), GmPath.GetName(lpSamplePath));
////	FILE *fpDemandRecs = fopen(namestring, "w+t");
////	fprintf(fpDemandRecs, "X,Y,Percentage\n");
////	for (int i=0; i<nRows; i++ )
////	{
////		fprintf(fpDemandRecs, "%lf,%lf,%lf\n", pDemandSites[i].dX, pDemandSites[i].dY, pDemandSites[i].dDist);
////	}
////	fclose(fpDemandRecs);
////
////
////	//
////	// Cleanup
////	//
////	if (pSampleX) delete[] pSampleX;
////	if (pSampleY) delete[] pSampleY;
////	for ( int i=0; i<nRows; i++ )
////	{
////		if (ppData[i]) delete ppData[i];
////		if (ppDistMatrix[i]) delete ppDistMatrix[i];
////	}
////	if (ppData) delete[] ppData;
////	if (ppDistMatrix) delete[] ppDistMatrix;
////	if (pDemandSites) delete[] pDemandSites;
////	return(true);
////}