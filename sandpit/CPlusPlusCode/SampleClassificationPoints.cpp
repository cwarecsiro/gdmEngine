//
// SampleClassificationPoints.cpp
//
#include "stdafx.h"
#include "GdmBinLib.h"
#include "Message.h"
#include "GDMBufferSizeDefs.h"
#include "ParamsW16.h"
#include "FilterSites.h"
#include "GdmExtractQuants.h"
#include "clsDoPath.h"
#include "clsEsriBinaryHeader.h"
#include "ESRIBinaryGridClass.h"
#include "GdmTransform.h"

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

#include <time.h>       /* time */
#include <random>

inline double GetXCentroid(long lOffset, int nCols, double dSize, double dXMin);
inline double GetYCentroid(long lOffset, int nCols, double dSize, double dYMax);
int GetNumGDMTransformedPreds1(char *pParams);
char **GetTransformPaths1(char *pParams);

//
// Create a Spatially Balanced sample table of N classes extracted from a GDM Classification Grid
//
bool CreateSpatiallyBalancedSamplesFromGDMClassification(char *pParams, int NumSamples, char *pOutPath, FPTR fptr)
{
	Message("CreateSpatiallyBalancedSamplesFromGDMClassification", "INFO");
	Message(pParams, "pParams");

	gmpath gmPath;
	char GDMParams[BUFF1024];
	GetProfileString("CONFIGURATION", "GDMPARAMS", GDMParams, pParams);
	Message(GDMParams, "GDMParams");

	char ClassificationDirPath[BUFF1024];
	GetProfileString("OUTPUTS", "OutDirectory", ClassificationDirPath, pParams);
	Message(ClassificationDirPath, "ClassificationDirPath");

	char ClassificationName[BUFF1024];
	GetProfileString("OUTPUTS", "OutName", ClassificationName, pParams);
	Message(ClassificationName, "ClassificationName");

	char ClassificationPath[BUFF1024];
	sprintf(ClassificationPath, "%s\\%s", ClassificationDirPath, ClassificationName);
	Message(ClassificationPath, "ClassificationPath");

	char Classification_FLT[BUFF1024];
	strcpy(Classification_FLT, gmPath.ChangeExtension(ClassificationPath, ".flt"));
	Message(Classification_FLT, "Classification_FLT");

	char Classification_HDR[BUFF1024];
	strcpy(Classification_HDR, gmPath.ChangeExtension(ClassificationPath, ".hdr"));
	Message(Classification_HDR, "Classification_HDR");



	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int NumClasses = GetProfileInt("INPUTS", "NumClasses", pParams);

	int *pClassCounts = new int[NumClasses];
	bool *pCanUse = new bool[NumClasses];
	bool *pCanSample = new bool[NumClasses]; // Can sample if there are more cells in the class than in the TargetSamplesPerGroup

	for (int i = 0; i < NumClasses; i++)
	{
		pClassCounts[i] = 0;
		pCanUse[i] = false;
		pCanSample[i] = false;
	}
	Message(NumClasses, "NumClasses");
	Message(NumSamples, "NumSamples");

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	EsriBinaryHeader *ebhClass = new EsriBinaryHeader(Classification_HDR);
	int nRows = ebhClass->GetNumRows();
	int nCols = ebhClass->GetNumCols();
	float fNoData = ebhClass->GetNoDataValue();
	double dSize = ebhClass->GetCellSize();
	double dXMin = ebhClass->GetMinX();
	double dYMax = ebhClass->GetMaxY();

	ESRIBinaryGridClass *ebgClass = new ESRIBinaryGridClass(Classification_FLT);
	if (!ebgClass->GetBinary()->IsValid())
	{
		Message("Invalid Classification Binary Float File", "ERROR");
		return(false);
	}

	//Message(dSize, "Cellsize");
	//Message(nRows, "nRows");
	//Message(nCols, "nCols");
	//Message(fNoData, "NoData");

	float *pClassRow = new float[nCols];
	int nCurrent = 0;
	fptr("Getting Grid Class Statistics...", nCurrent);
	for (int i = 0; i < nRows; i++)
	{
		if (i * 100 / nRows > nCurrent)
		{
			nCurrent = i * 100 / nRows;
			fptr("Getting Grid Class Statistics...", nCurrent);
		}

		// read a row of data...
		ebgClass->GetBinary()->ReadFloat(pClassRow, nCols);

		for (int j = 0; j < nCols; j++)
		{
			if (pClassRow[j] != fNoData)
			{
				//++pClassCounts[int(pClassRow[j]) + 1];
				++pClassCounts[int(pClassRow[j])-1];
			}
		}
	}


	//
	// Count the number of cells in classes greater than zero...
	//
	int nMaxCount = 0;
	int nValidClasses = 0;
	for (int i = 0; i < NumClasses; i++)
	{
		if (pClassCounts[i] > 0) ++nValidClasses;

		if (pClassCounts[i] > nMaxCount) nMaxCount = pClassCounts[i];
	}
	//Message(nValidClasses, "nValidClasses");
	//Message(nMaxCount, "nMaxCount");

	int nTargetSamplesPerGroup = NumSamples / NumClasses;
	if (nValidClasses > 0)
	{
		nTargetSamplesPerGroup = NumSamples / NumClasses;
	}
	//Message(nTargetSamplesPerGroup, "nTargetSamplesPerGroup");


	//FILE *fp = fopen("D:\\SGA_August_2017_Tutorial\\WD_01\\Class_500\\ClassCounts.csv", "w+t");
	//fprintf(fp, "Class,Count,CanUse,CanSample\n");
	for (int i = 0; i < NumClasses; i++)
	{
		pCanUse[i] = (pClassCounts[i] > 0) ? true : false;
		pCanSample[i] = (pClassCounts[i] > nTargetSamplesPerGroup) ? true : false;
		//fprintf(fp, "%d,%d,%d,%d\n", i + 1, pClassCounts[i], pCanUse[i], pCanSample[i]);
	}
	//fclose(fp);


	//
	// Loop thru all the Valid Classes and extract samples
	//
	//
	// allocate a vector for the offsets
	long *pOffsets = new long[nMaxCount];

	int nTranPreds = GetNumGDMTransformedPreds1(GDMParams);
	//Message(nTranPreds, "nTranPreds");

	char **ppTranGridPaths = GetTransformPaths1(GDMParams);
	ESRIBinaryGridClass **ebgTrans = new ESRIBinaryGridClass *[nTranPreds];
	for (int i = 0; i < nTranPreds; i++)
	{
		//Message(ppTranGridPaths[i], "ppTranGridPaths[i]");
		ebgTrans[i] = new ESRIBinaryGridClass(gmPath.ChangeExtension(ppTranGridPaths[i], ".flt"));
	}

	FILE *fpOut = fopen(pOutPath, "w+t");
	fprintf(fpOut, "Class,X,Y");
	for (int i = 0; i < nTranPreds; i++)
	{
		fprintf(fpOut, ",%s", gmPath.GetName(ppTranGridPaths[i]));
		if (i == nTranPreds - 1) 
			fprintf(fpOut, "\n");
	}

	nCurrent = 0;
	fptr("Getting Grid Class Samples...", nCurrent);

	for (int cls = 0; cls < NumClasses; cls++)
	{
		if (cls * 100 / NumClasses > nCurrent)
		{
			nCurrent = cls * 100 / NumClasses;
			fptr("Getting Grid Class Samples...", nCurrent);
		}

		if (pCanUse[cls])
		{
			int nThis = 0;
			//for (int i = 0; i < nMaxCount; i++) pOffsets[i] = -1L;
			ebgClass->GetBinary()->SeekToStart();

			// count the number of cells in this class
			for (int i = 0; i < nRows; i++)
			{
				// read a row of data...
				ebgClass->GetBinary()->ReadFloat(pClassRow, nCols);

				for (int j = 0; j < nCols; j++)
				{
					if (pClassRow[j] == cls + 1)
					{
						pOffsets[nThis] = (i*nCols) + j;
						++nThis;
					}
				}
			}

			// get random sample or all cells from the Offsets vector...
			int NumOffsets = nThis;
			if (pCanSample[cls])
			{
				//
				// use a sample of size nTargetSamplesPerGroup from the cells in this class
				//
				std::srand(time(NULL));
				const int nMin = 0;
				const int nMax = NumOffsets - 1;
				std::default_random_engine generator;
				std::uniform_int_distribution<int> distribution(nMin, nMax);

				for (int i = 0; i < nTargetSamplesPerGroup; i++)
				{
					// get a random index
					int nIndex = distribution(generator);

					fprintf(fpOut,
						"%d,%lf,%lf",
						cls + 1,
						GetXCentroid(pOffsets[nIndex], nCols, dSize, dXMin),
						GetYCentroid(pOffsets[nIndex], nCols, dSize, dYMax));

					// get the predictor values
					for (int k = 0; k < nTranPreds; k++)
					{
						ebgTrans[k]->GetBinary()->SeekTo(pOffsets[nIndex]*sizeof(float));
						float val;
						ebgTrans[k]->GetBinary()->ReadFloat(&val, 1);
						fprintf(fpOut, ",%2.6f", val);

						if (k == nTranPreds - 1)
							fprintf(fpOut, "\n");
					}
				}
			}
			else
			{
				// use all the cells in this class
				for (int i = 0; i < pClassCounts[cls]; i++)
				{
					// just use all the cells in the class
					fprintf(fpOut,
						"%d,%lf,%lf",
						cls + 1,
						GetXCentroid(pOffsets[i], nCols, dSize, dXMin),
						GetYCentroid(pOffsets[i], nCols, dSize, dYMax)); 

					// get the predictor values
					for (int k = 0; k < nTranPreds; k++)
					{
						ebgTrans[k]->GetBinary()->SeekTo(pOffsets[i] * sizeof(float));
						float val;
						ebgTrans[k]->GetBinary()->ReadFloat(&val, 1);
						fprintf(fpOut, ",%2.6f", val);

						if (k == nTranPreds - 1)
							fprintf(fpOut, "\n");
					}
				}
			}
		}

	} // for (int cls = 0; cls < NumClasses; cls++)

	if (fpOut)
		fclose(fpOut);

	for (int i = 0; i < nTranPreds; i++)
	{
		if (ebgTrans[i]) delete ebgTrans[i];
	}
	if (ebgTrans) delete ebgTrans;

	if (ebhClass) delete ebhClass;
	if (ebgClass) delete ebgClass;
	if (pClassRow) delete[] pClassRow;
	if (pClassCounts) delete[] pClassCounts;
	if (pCanUse) delete[] pCanUse;
	if (pCanSample) delete[] pCanSample;
	if (pOffsets) delete[] pOffsets;
	fptr("Ready...", 0);
	Message("Done");
	return(true);
}


inline double GetXCentroid(long lOffset, int nCols, double dSize, double dXMin)
{
	return(dXMin + ((lOffset % nCols) * dSize + (dSize / 2.0)));
}


inline double GetYCentroid(long lOffset, int nCols, double dSize, double dYMax)
{
	return(dYMax - ((lOffset / nCols) * dSize - (dSize / 2.0)));
}


//
// Get number of transformed predictors from a GDM model
//
int GetNumGDMTransformedPreds1(char *pParams)
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
char **GetTransformPaths1(char *pParams)
{
	int nTransforms = GetNumGDMTransformedPreds1(pParams);
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

