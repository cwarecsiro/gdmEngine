//
// GDMSigTest.cpp
//
#include "stdafx.h"
#include "GdmBinLib.h"
#include "myCallback.h"
#include "Message.h"
#include "GDMSigTest.h"
#include "GdmGuiLib.h"
#include "GDMBufferSizeDefs.h"
#include "ParamsW16.h"
#include "clsDoPath.h"
#include "GdmTabLib.h"
#include "NNLS_Double.h"
#include "GdmExtractQuants.h"
#include <algorithm>
#include <vector>
#include <ctime>
#include <cstdlib>


//
// define an Epsilon value that allows the marginal change in dev explained for a random permutation
//
double d_EXCEED_EPSILON = 0.0;


//
// Save basic GDM parameters for file (if we have no quantiles etc)
//
void SaveGDMParams(char **WDPath, char **paramFileName, char **tablePath, int *NumPreds, bool UseGeoDist)
{
	//Message(*WDPath, "Workspace Path");
	//Message(*paramFileName, "Parameter File Name");

	gmpath gmp;
	char ParamFilePath[1024];
	sprintf(ParamFilePath, "%s%s%s", *WDPath, "\\", *paramFileName);
	
	//sprintf(ParamFilePath, "%s\\%s_PARAMS.txt", gmp.GetDirectoryPath(*WDPath), gmp.GetName(*paramFileName));
	//Message(ParamFilePath, "Parameter File Path");
	
	FILE *fp = fopen(ParamFilePath, "w+t");
	fprintf(fp, "################################################################################################\n");
	fprintf(fp, "###\n");
	fprintf(fp, "### GDM Parameter File generated from R\n");
	fprintf(fp, "###\n");
	fprintf(fp, "\n\n");
	fprintf(fp, "[GDMODEL]\n");
	fprintf(fp, "WorkspacePath=%s\n", *WDPath);
	fprintf(fp, "UseEuclidean=%d\n", UseGeoDist?1:0);
	fprintf(fp, "UseSubSample=0\n");
	fprintf(fp, "NumSamples=0\n");
	fprintf(fp, "Intercept=0\n");
	fprintf(fp, "NullDeviance=0\n");
	fprintf(fp, "GDMDeviance=0\n");
	fprintf(fp, "DevExplained=0\n");
	fprintf(fp, "NumberSitePairs=0\n");
	fprintf(fp, "SumOfCoefficients=0\n");
	fprintf(fp, "NumberOfActivePredictors=0\n");
	fprintf(fp, "NumberCoefficients>0=0\n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "[RESPONSE]\n");
	fprintf(fp, "InputData=%s\n", *tablePath);
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "[PREDICTORS]\n");
	int nPreds = *NumPreds;
	fprintf(fp, "NumPredictors=%d\n", nPreds);

	char tabheader[4096];
	FILE *fpTab = fopen(*tablePath, "r+t");
	fgets(tabheader, 4096, fpTab);
	fclose(fpTab);
	char seps[] = ",\n";
	char *p = strtok(tabheader, seps);
	for (int i = 0; i < 5; i++) 
		p = strtok(NULL, seps);
	for (int i = 1; i <= nPreds; i++)
	{
		p = strtok(NULL, seps);
		fprintf(fp, "EnvTab%d=%s\n", i, p);
	}
	fprintf(fp, "\n");
	fprintf(fp, "\n");

	for (int i = 1; i <= nPreds; i++)
	{
		fprintf(fp, "UseEnv%d=1\n", i);
	}
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	for (int i = 1; i <= nPreds; i++)
	{
		fprintf(fp, "PredSpl%d=3\n", i);
	}
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	for (int i = 1; i <= nPreds; i++)
	{
		fprintf(fp, "QuantType%d=3\n", i);
	}
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "EuclSpl=3\n");
	fprintf(fp, "EuclSplVal1=0\n");
	fprintf(fp, "EuclSplVal2=0\n");
	fprintf(fp, "EuclSplVal3=0\n");
	fprintf(fp, "EuclCoef1=0\n");
	fprintf(fp, "EuclCoef2=0\n");
	fprintf(fp, "EuclCoef3=0\n");
	fprintf(fp, "EuclQuantType=0\n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fclose(fp);
}


//
// Extract Quantiles from Table Data and update parameter file
//
bool ExtractAndUpdateQuantilesSigTest(char **ppParams)
{
	//Message("ExtractAndUpdateQuantilesSigTest");
	char *pParams = *ppParams;
	//Message(pParams, "pParams");

	//
	// Get the data table path
	//
	char *DataTablePath = GetResponseData(pParams);
	
	//
	// Get a memory image of the composite data file
	//
	GDM_INT nTableRows = 0;
	GDM_INT nTableCols = 0;
	double *pTableData = GetCompositeTableIntoMemorySigTest(DataTablePath, &nTableRows, &nTableCols);
	if (NULL == pTableData)
	{
		Message("Cannot GetCompositeTableIntoMemory() in ExtractAndUpdateQuantiles", "ERROR");
		if (DataTablePath) delete[] DataTablePath;
		return(false);
	}

	// Free the DataTablePath now that we don't need it anymore
	if (DataTablePath) delete[] DataTablePath;

	//
	// Extract the geographic distance data
	//
	int nCurrent = 0;
	double *pGeoDist = new double[nTableRows];
	for (int i = 0; i<nTableRows; i++)
	{
		double distX = pTableData[(COL_SITE1_X0*nTableRows) + i] - pTableData[(COL_SITE2_X1*nTableRows) + i];
		double distY = pTableData[(COL_SITE1_Y0*nTableRows) + i] - pTableData[(COL_SITE2_Y1*nTableRows) + i];
		pGeoDist[i] = sqrt((distX * distX) + (distY * distY));
	}

	// sort the vector
	SortVector(pGeoDist, (int)nTableRows);

	// extract the geographic distance quantiles and update the parameter file
	int nQuants = GetEuclideanSplines(pParams);
	char lpKey[64];
	for (int i = 0; i<nQuants; i++)
	{
		if (i == 0)
		{
			sprintf(lpKey, "EuclSplVal%d", i + 1);
			SetProfileDouble("PREDICTORS", lpKey, pGeoDist[0], pParams);
		}

		else if ((i>0) && (i<nQuants - 1))
		{
			int nPercentile = (int)(i * 100.0 / (double)(nQuants - 1));
			double dQuant = GetQuantAsPercentile(pGeoDist, (int)nTableRows, nPercentile);
			sprintf(lpKey, "EuclSplVal%d", i + 1);
			SetProfileDouble("PREDICTORS", lpKey, dQuant, pParams);
		}

		else
		{
			sprintf(lpKey, "EuclSplVal%d", i + 1);
			SetProfileDouble("PREDICTORS", lpKey, pGeoDist[(int)nTableRows-1], pParams);
		}
	}

	// clean up
	if (pGeoDist) delete[] pGeoDist;

	//
	// Extract the quantiles from the 'In-Use" predictors
	//
	double *pPredDist = new double[nTableRows * 2];
	int nPreds = GetNumPredictors(pParams);
	for (int p = 0; p<nPreds; p++)
	{
		if (GetPredictorInUseAt(pParams, p + 1))
		{
			int nThis = 0;
			for (int i = 0; i<nTableRows; i++)
			{
				pPredDist[nThis++] = pTableData[((LEADING_COLS + p)*nTableRows) + i];
				pPredDist[nThis++] = pTableData[((LEADING_COLS + p + nPreds)*nTableRows) + i];
			}

			// sort the vector
			SortVector(pPredDist, (int)nTableRows * 2);

			// extract the predictor quantiles and update the parameter file
			nQuants = GetPredictorSplinesAt(pParams, p + 1);
			for (int i = 0; i<nQuants; i++)
			{
				if (i == 0)
				{
					SetPredictorSplineAt(pParams, p + 1, i + 1, pPredDist[0]);
				}

				else if ((i>0) && (i<nQuants - 1))
				{
					int nPercentile = (int)(i * 100.0 / (double)(nQuants - 1));
					double dQuant = GetQuantAsPercentile(pPredDist, (int)nTableRows * 2, nPercentile);
					SetPredictorSplineAt(pParams, p + 1, i + 1, dQuant);
				} // if ((i>0) && (i<nQuants-1))

				else
				{
					SetPredictorSplineAt(pParams, p + 1, i + 1, pPredDist[(int)((nTableRows * 2)-1)]);
				}
			} // for ( int i=0; i<nQuants; i++ )

			AppendRow(pParams);

		} // if (GetPredictorInUseAt(pParams, p+1))
	} // for ( int p=0; p<nPreds; p++ )

	if (pPredDist) delete[] pPredDist;
	if (pTableData) delete[] pTableData;
	return(true);
}


//
// Get a memory image of the composite data file
//
double *GetCompositeTableIntoMemorySigTest(char *pDataPath, GDM_INT *pTableRows, GDM_INT *pTableCols)
{
	//
	// Get the table metrics
	//
	int nRows = CountRows64SigTest(pDataPath);
	if (nRows < 1)
	{
		Message("Cannot CountRows in GetCompositeTableIntoMemory()", "ERROR");
		*pTableRows = *pTableCols = 0;
		return(NULL);
	}

	//Message(nRows, "nRows");

	int nCols = CountColumns64SigTest(pDataPath);
	if (nCols < 1)
	{
		Message("Cannot CountColumns in GetCompositeTableIntoMemory()", "ERROR");
		*pTableRows = *pTableCols = 0;
		return(NULL);
	}

	//Message(nCols, "nCols");
	//Message(nRows * nCols, "nRows * nCols");
	//Message(nRows * nCols * 8, "nRows * nCols * 8");

	if (sizeof(int *) == 4)  // 32 bit OS
	{
		if ((nRows * nCols * 8) <= 0)
		{
			Message(nRows, "NumRows");
			Message(nCols, "NumCols");
			Message("Cannot allocate pData in GetCompositeTableIntoMemory", "ERROR");
			return(NULL);
		}
	}

	double *pData = new double[nRows * nCols];
	if (NULL == pData)
	{
		Message("Cannot allocate pData in GetCompositeTableIntoMemory", "ERROR");
		return(NULL);
	}

	char *pRowData = new char[TABLE_ROW_BUFFSIZE];
	if (NULL == pRowData)
	{
		Message("Cannot allocate pRowData in GetCompositeTableIntoMemory", "ERROR");
		if (pData) delete[] pData;
		return(NULL);
	}


	// open text file
	ifstream* pmyFile = new ifstream; // On the heap
	pmyFile->open(pDataPath);

	if (!pmyFile->is_open())
	{
		Message("Cannot open pDataPath in GetCompositeTableIntoMemory", "ERROR");
		if (pData) delete[] pData;
		if (pRowData) delete[] pRowData;
		return(NULL);
	}

	// get header
	pmyFile->getline(pRowData, TABLE_ROW_BUFFSIZE);

	int nCurrent = 0;
	char pBuff[128];
	for (int i = 0; i<nRows; i++)
	{
		if ((long long)(i) * 100 / nRows > nCurrent)
		{
			nCurrent = (long long)(i) * 100 / nRows;
			sprintf(pBuff, "Getting composite table into memory (Row %d of %d)", i + 1, nRows);
		}

		// get current row
		pmyFile->getline(pRowData, TABLE_ROW_BUFFSIZE);

		// extract data from current row
		char *p = strtok(pRowData, ",\n");
		if (NULL == p)
		{
			char qqq[64];
			sprintf(qqq, "nRows: %d   nCols: %d", nRows, nCols);
			Message(qqq);
			sprintf(qqq, "Got a NULL at Row: %d and Col: %d", i, 0);
			Message(qqq);

			pmyFile->close();
			if (pRowData) delete[] pRowData;
			*pTableRows = 0;
			*pTableCols = 0;
			return(NULL);
		}
		pData[i] = atof(p);

		for (int j = 1; j<nCols; j++)
		{
			p = strtok(NULL, ",\n");

			if (NULL == p)
			{
				char qqq[64];
				sprintf(qqq, "nRows: %d   nCols: %d", nRows, nCols);
				Message(qqq);
				sprintf(qqq, "Got a NULL at Row: %d and Col: %d", i, j);
				Message(qqq);

				pmyFile->close();
				if (pRowData) delete[] pRowData;
				*pTableRows = 0;
				*pTableCols = 0;
				return(NULL);
			}
			pData[(j*nRows) + i] = atof(p);
		}
	}

	pmyFile->close();
	if (pRowData) delete[] pRowData;
	*pTableRows = (GDM_INT)nRows;
	*pTableCols = (GDM_INT)nCols;
	return(pData);
}


//
// Using the 64bit fstream class, count the number of data rows excluding the header in a text file
//
int CountRows64SigTest(char *pPath)
{
	int nRows = 0;
	char *p = new char[TABLE_ROW_BUFFSIZE];
	if (NULL == p)
	{
		Message("Cannot allocate p in CountRows", "ERROR");
		return(nRows);
	}

	// open text file
	ifstream* pmyFile = new ifstream; // On the heap
	pmyFile->open(pPath);

	if (!pmyFile->is_open())
	{
		Message("Cannot open pPath in CountRows", "ERROR");
		if (p) delete[] p;
		return(nRows);
	}

	// get header
	pmyFile->getline(p, TABLE_ROW_BUFFSIZE);

	while (!pmyFile->eof())
	{
		pmyFile->getline(p, TABLE_ROW_BUFFSIZE);

		// make sure that there is some data in the line
		if (strlen(p) > 4)
			++nRows;
	}

	pmyFile->close();
	if (p) delete[] p;
	return(nRows);
}



//
// Using the 64bit fstream class, count the number of comma delimited columns in a text file
//
int CountColumns64SigTest(char *pPath)
{
	int nColumns = 0;
	char *pHeader = new char[TABLE_ROW_BUFFSIZE];
	if (NULL == pHeader)
	{
		Message("Cannot allocate pHeader in CountColumns", "ERROR");
		return(nColumns);
	}

	ifstream* pmyFile = new ifstream; // On the heap
	pmyFile->open(pPath);

	if (!pmyFile->is_open())
	{
		Message("Cannot open pPath in CountColumns64", "ERROR");
		if (pHeader) delete[] pHeader;
		return(nColumns);
	}

	// get header and close file
	pmyFile->getline(pHeader, TABLE_ROW_BUFFSIZE);
	pmyFile->close();

	// count columns
	char *p = strtok(pHeader, ",\n");
	while (p)
	{
		++nColumns;
		p = strtok(NULL, ",\n");
	}

	if (pHeader) delete[] pHeader;
	return(nColumns);
}



//
// Run a GDM progressively dropping each of the predictors
// from the analysis to determine its contribution to the fit.
//
bool DoSigTestGDM(char **ppParams, int *pIterations)
{
	char *pParams = *ppParams;
	int nIterations = *pIterations;

	gmpath gmp;
	char logpath[512];
	sprintf(logpath, "%s\\%s_LOG.txt", gmp.GetDirectoryPath(pParams), gmp.GetName(pParams));
	FILE *fpLog = NULL;
	fpLog = fopen(logpath, "w+t");
	fprintf(fpLog, "Logfile for GDM Significance Test: %s\n", pParams);
	fprintf(fpLog, "Iterations per predictor: %d\n", nIterations);
	int nEnvPreds = GetNumPredictors(pParams);
	if (UseEuclidean(pParams)) ++nEnvPreds;
	fprintf(fpLog, "Number of predictors: %d\n\n", nEnvPreds);
	fclose(fpLog);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Do the full GDM and retain the binary matrix file
	// for sub-setting to use in the significance test
	//
	double *pResponse = NULL;
	double *pWeights = NULL;
	int NumRows, NumCols;
	time_t now = time(0);
	char* dt = ctime(&now);
	fpLog = fopen(logpath, "a+t");
	fprintf(fpLog, "Start Baseline GDM at %s", dt);
	fclose(fpLog);
	if (!DoBaseSigGDM(pParams, false, &pResponse, &pWeights, &NumRows, &NumCols, NULL))
	{
		Message("Cannot DoBaseSigGDM in DoSigTestGDM", "ERROR");
		return(false);
	}
	now = time(0);
	dt = ctime(&now);
	fpLog = fopen(logpath, "a+t");
	fprintf(fpLog, " Done Baseline GDM at %s\n", dt);
	fclose(fpLog);

	if (pResponse == NULL)
	{
		Message("Got a NULL Response Vector", "ERROR");
		return(false);
	}

	if (pWeights == NULL)
	{
		Message("Got a NULL Weights Vector", "ERROR");
		return(false);
	}

	double DevExplained = GetProfileDouble("GDMODEL", "DevExplained", pParams);
	//Message(DevExplained, "DevExplained");
	//uptr(-1, -1, DevExplained, -1, -1, -1, -1);
	//
	////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// setup path to full GDM binary matrix file
	char *lpTmpFile = GetWorkspacePath(pParams);
	std::strcat(lpTmpFile, "\\gdmtmp.bin");

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Loop thru each predictor (including Geographic Dist if DoGeo is true)
	//
	//
	// do euclidean distance if required
	double dEuclSigProportion = -1.0;
	double dEuclMinExceed = 0.0;
	double dEuclMaxExceed = 0.0;
	int nEuclideanSplines = GetEuclideanSplines(pParams);
	int ColumnJump = 1; // skip the intercept column
	if (UseEuclidean(pParams))
	{
		now = time(0);
		dt = ctime(&now);
		fpLog = fopen(logpath, "a+t");
		fprintf(fpLog, "Start Significance Test for Geographic Distance at %s", dt);
		fclose(fpLog);

		// do nIterations of significance testing for the geographic distance predictor
		dEuclSigProportion = DoPredictorSigGDM(pParams, lpTmpFile, 
			                                   pResponse, pWeights, DevExplained, 
			                                   NumRows, NumCols,		                                       	                                   
			                                   0, nIterations,
			                                   ColumnJump, nEuclideanSplines,
											   &dEuclMinExceed, &dEuclMaxExceed,
			                                   NULL, NULL);

		//uptr(dEuclSigProportion, dEuclSigProportion, dEuclSigProportion, -1, -1, -1, 0);
		ColumnJump += nEuclideanSplines;

		now = time(0);
		dt = ctime(&now);
		fpLog = fopen(logpath, "a+t");
		fprintf(fpLog, " Done Significance Test for Geographic Distance at %s\n", dt);
		fclose(fpLog);
	}

	// do the environmental predictors
	int nPreds = GetNumPredictors(pParams);
	double *pPredSignificance = new double[nPreds];
	double *pPredMinExceed = new double[nPreds];
	double *pPredMaxExceed = new double[nPreds];

	// loop thru the environmental predictors...
	char pBuff[128];
	for (int i = 1; i <= nPreds; i++)
	{
		// initialise
		pPredSignificance[i-1] = -1.0;
		int nSplines = GetPredictorSplinesAt(pParams, i);
		
		if (GetPredictorInUseAt(pParams, i))
		{
			GetPredictorNameAt(pParams, pBuff, i);
			now = time(0);
			dt = ctime(&now);
			fpLog = fopen(logpath, "a+t");
			fprintf(fpLog, "Start Significance Test for %s at %s", pBuff, dt);
			fclose(fpLog);

			// do nIterations of significance testing for this environmental predictor
			pPredSignificance[i-1] = DoPredictorSigGDM(pParams, lpTmpFile,
				                                       pResponse, pWeights, DevExplained,
				                                       NumRows, NumCols,
				                                       i, nIterations,
				                                       ColumnJump, nSplines,
				                                       &pPredMinExceed[i-1], &pPredMaxExceed[i-1],
				                                       NULL, NULL);

			//uptr(pPredSignificance[i-1], pPredSignificance[i-1], pPredSignificance[i-1], -1, -1, -1, i);

			// update column jump for this predictor
			ColumnJump += nSplines;

			now = time(0);
			dt = ctime(&now);
			fpLog = fopen(logpath, "a+t");
			fprintf(fpLog, " Done Significance Test for %s at %s\n", pBuff, dt);
			fclose(fpLog);
		}

	} // for (int i = 1; i <= nPreds; i++)


	//
	// Write results to the GDM parameter file...
	//
	/*SetProfileDouble("SIGNIFICANCE", "BaseLineDevExplained", DevExplained, pParams);
	SetProfileInt("SIGNIFICANCE", "Iterations", nIterations, pParams);
	char buff1[64];
	char buff2[64];
	if (UseEuclidean(pParams))
	{
		sprintf(buff2, "%0.8lf", dEuclSigProportion);
		SetProfileString("SIGNIFICANCE", "PVal_Geo", buff2, pParams);

		sprintf(buff2, "%0.8lf", dEuclMinExceed);
		SetProfileString("SIGNIFICANCE", "MinX_Geo", buff2, pParams);

		sprintf(buff2, "%0.8lf", dEuclMaxExceed);
		SetProfileString("SIGNIFICANCE", "MaxX_Geo", buff2, pParams);
	}	
	for (int i = 1; i <= nPreds; i++)
	{
		if (GetPredictorInUseAt(pParams, i))
		{
			sprintf(buff1, "PVal_%02d", i);
			sprintf(buff2, "%0.8lf", pPredSignificance[i-1]);
			SetProfileString("SIGNIFICANCE", buff1, buff2, pParams);

			sprintf(buff1, "MinX_%02d", i);
			sprintf(buff2, "%0.8lf", pPredMinExceed[i-1]);
			SetProfileString("SIGNIFICANCE", buff1, buff2, pParams);

			sprintf(buff1, "MaxX_%02d", i);
			sprintf(buff2, "%0.8lf", pPredMaxExceed[i-1]);
			SetProfileString("SIGNIFICANCE", buff1, buff2, pParams);
		}
	}
*/
	//
	// write backup significance file like the BWE tool produces
	//
	//gmpath gmp;
	char tmpBuff[512];
	char resultspath[512];
	sprintf(resultspath, "%s\\%s_SIG.csv", gmp.GetDirectoryPath(pParams), gmp.GetName(pParams));
	FILE *fp = fopen(resultspath, "w+t");
	fprintf(fp, "Index,Predictor,Iterations,Baseline_%%DevExp.,P_Value,Min_Marginal_%%DevExp,Max_Marginal_%%DevExp\n");

	if (UseEuclidean(pParams))
	{
		fprintf(fp, 
			    "%d,%s,%d,%0.8lf,%0.8lf,%0.8lf,%0.8lf\n", 
			    0, "GeoDist", nIterations, DevExplained, dEuclSigProportion, dEuclMinExceed, dEuclMaxExceed);
	}

	for (int i = 1; i <= nPreds; i++)
	{
		GetPredictorNameAt(pParams, tmpBuff, i);

		if (GetPredictorInUseAt(pParams, i))
		{
			fprintf(fp, 
				    "%d,%s,%d,%0.8lf,%0.8lf,%0.8lf,%0.8lf\n", 
				    i, tmpBuff, nIterations, DevExplained, pPredSignificance[i - 1], pPredMinExceed[i - 1], pPredMaxExceed[i - 1]);
		}

		/*else
		{
			fprintf(fp, "%s,%d,%s,%s,%s,%s\n", tmpBuff, 0, "N/A", "N/A", "N/A", "N/A");
		}*/
	}

	fclose(fp);

	//SetProfileString("BASELINEMODEL", "GDMParamPath", pParams, newpath);
	//SetProfileDouble("BASELINEMODEL", "BaseLineDevExplained", DevExplained, newpath);
	//SetProfileInt("SIGNIFICANCE", "Iterations", nIterations, newpath);
	//if (UseEuclidean(pParams))
	//{
	//	sprintf(buff2, "%0.8lf", dEuclSigProportion);
	//	SetProfileString("SIGNIFICANCE", "PVal_Geo", buff2, newpath);
	//	sprintf(buff2, "%0.8lf", dEuclMinExceed);
	//	SetProfileString("SIGNIFICANCE", "MinX_Geo", buff2, newpath);
	//	sprintf(buff2, "%0.8lf", dEuclMaxExceed);
	//	SetProfileString("SIGNIFICANCE", "MaxX_Geo", buff2, newpath);
	//}
	//for (int i = 1; i <= nPreds; i++)
	//{
	//	if (GetPredictorInUseAt(pParams, i))
	//	{
	//		sprintf(buff1, "PVal_%02d", i);
	//		sprintf(buff2, "%0.8lf", pPredSignificance[i - 1]);
	//		SetProfileString("SIGNIFICANCE", buff1, buff2, newpath);
	//		sprintf(buff1, "MinX_%02d", i);
	//		sprintf(buff2, "%0.8lf", pPredMinExceed[i - 1]);
	//		SetProfileString("SIGNIFICANCE", buff1, buff2, newpath);
	//		sprintf(buff1, "MaxX_%02d", i);
	//		sprintf(buff2, "%0.8lf", pPredMaxExceed[i - 1]);
	//		SetProfileString("SIGNIFICANCE", buff1, buff2, newpath);
	//	}
	//}


	// remove the binary matrix file and clean up
	if ((_access(lpTmpFile, 0)) != -1) remove(lpTmpFile);
	if (lpTmpFile) delete[] lpTmpFile;
	if (pPredSignificance) delete[] pPredSignificance;
	if (pResponse) delete[] pResponse;
	if (pWeights) delete[] pWeights;

	now = time(0);
	dt = ctime(&now);
	fpLog = fopen(logpath, "a+t");
	fprintf(fpLog, "Done Significance Test at %s\n", dt);
	fclose(fpLog);

	return(true);
}




//
// Run the core GDM for a Significance Test
//
bool DoBaseSigGDM(char *pParams, bool DeleteBinaryFile,
	              double **pResponse, double **pWeights,
	              int *NumRows, int *NumCols,
	              FPTR fptr)
{
	//
	// Get all the spline counts into a vector
	//
	int nSplineCounts = 0;
	int *theSplineCounts = GetSplineCountsFromParams(pParams, &nSplineCounts);


	//
	// Get all the quantiles into a vector 
	//
	int nQuantiles = 0;
	double *theQuantiles = GetQuantilesFromParams(pParams, &nQuantiles);


	//
	// Get the data table path
	//
	//char *DataTablePath = GetFilteredResponseData(pParams);
	char *DataTablePath = GetResponseData(pParams);
	//Message(DataTablePath, "DataTablePath");
	gmpath gmp;
	if (!gmp.FileExists(DataTablePath))
	{
		Message("DataTablePath does NOT exist in DoBaseSigGDM", "ERROR");
		if (DataTablePath) delete[] DataTablePath;
		return(false);
	}


	//
	// Get a memory image of the composite data file
	//
	GDM_INT nTableRows = 0;
	GDM_INT nTableCols = 0;
	//double *pTableData = GetCompositeTableIntoMemory(DataTablePath, &nTableRows, &nTableCols, false, fptr);
	double *pTableData = GetCompositeTableIntoMemorySigTest(DataTablePath, &nTableRows, &nTableCols);
	if (NULL == pTableData)
	{
		Message("Cannot GetCompositeTableIntoMemory() in DoBaseSigGDM", "ERROR");
		if (DataTablePath) delete[] DataTablePath;
		if (theSplineCounts) delete[] theSplineCounts;
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}
	//PrintTableData(pParams, pTableData, nTableRows, nTableCols, false);
	// Free the DataTablePath now that we don't need it anymore
	if (DataTablePath) delete[] DataTablePath;


	//
	// Create a binary file image of the splined predictor data to pass to the matrix regression
	// There will always be the same number of rows in the table data as the predictor matrix
	// but the number of matrix cols may be less if not all predictors are 'In-use'.
	//
	GDM_INT nMatrixCols = 0;
	char *lpTmpFile = GetWorkspacePath(pParams);
	strcat(lpTmpFile, "\\gdmtmp.bin");

	if (!CreatePredictorBinary(lpTmpFile, pParams,
		                       pTableData, nTableRows, nTableCols,
		                       theSplineCounts, theQuantiles,
		                       &nMatrixCols, false, NULL))
	{
		Message("Cannot CreatePredictorBinary() in DoBaseSigGDM", "ERROR");
		if (theSplineCounts) delete[] theSplineCounts;
		if (theQuantiles) delete[] theQuantiles;
		if (pTableData) delete[] pTableData;
		return(false);
	}
	//PrintBinary2CSV(lpTmpFile, pParams, theSplineCounts, nTableRows, nMatrixCols, fptr);

	//
	// Extract the response column and the weights column before freeing the table data memory block
	//
	*pResponse = new double[nTableRows];
	memcpy(*pResponse, &pTableData[COL_RESPONSE * nTableRows], nTableRows * sizeof(double));
	if (NULL == *pResponse)
	{
		Message("Cannot extract Response column in DoBaseSigGDM", "ERROR");
		if (pTableData) delete[] pTableData;
		if (theSplineCounts) delete[] theSplineCounts;
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}

	*pWeights = new double[nTableRows];
	memcpy(*pWeights, &pTableData[COL_WEIGHTS * nTableRows], nTableRows * sizeof(double));
	if (NULL == *pWeights)
	{
		Message("Cannot extract Weights column in DoBaseSigGDM", "ERROR");
		if (*pResponse) delete[] * pResponse;
		if (pTableData) delete[] pTableData;
		if (theSplineCounts) delete[] theSplineCounts;
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}
	// Free the Table data now that we don't need it anymore
	if (pTableData) delete[] pTableData;


	//
	// Create the memory for the Predictor Matrix
	//
	double *pPredictorMatrix = new double[nTableRows * nMatrixCols];
	if (NULL == pPredictorMatrix)
	{
		Message("Cannot allocate pPredictorMatrix in DoBaseSigGDM", "ERROR");
		if (*pResponse) delete[] * pResponse;
		if (*pWeights) delete[] * pWeights;
		if (theSplineCounts) delete[] theSplineCounts;
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}


	//
	// Populate the Predictor Matrix from the binary image created in CreatePredictorBinary()
	//
	int h = _open(lpTmpFile, _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE);
	if (h < 0)
	{
		Message("Cannot open binary file for READ in DoBaseSigGDM", "ERROR");
		if (pPredictorMatrix) delete[] pPredictorMatrix;
		if (*pResponse) delete[] * pResponse;
		if (*pWeights) delete[] * pWeights;
		if (theSplineCounts) delete[] theSplineCounts;
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}
	//_read(h, pPredictorMatrix, nTableRows * nMatrixCols * sizeof( double ));
	double *pTmp = pPredictorMatrix;
	for (int i = 0; i<nMatrixCols; i++)
	{
		_read(h, pTmp, (unsigned)nTableRows * sizeof(double));
		pTmp += nTableRows;
	}
	_close(h);

	//int nSize = int(nMatrixCols * nTableRows * sizeof(double));
	//Message(nSize, "Matrix Size");

	//
	// Do the matrix regression
	//
	if (fptr)
		fptr("Performing GDM regression...", 5);
	double dGDMDeviance;
	double *pCoefficients = WeightedNNLSRegression(lpTmpFile,
		pPredictorMatrix,
		nTableRows,
		nMatrixCols,
		*pResponse,
		&dGDMDeviance,
		*pWeights,
		NULL);
	if (NULL == pCoefficients)
	{
		Message("pCoefficients are NULL", "ERROR in DoBaseSigGDM");
		if ((_access(lpTmpFile, 0)) != -1) remove(lpTmpFile);
		if (pPredictorMatrix) delete[] pPredictorMatrix;
		if (*pResponse) delete[] * pResponse;
		if (*pWeights) delete[] * pWeights;
		if (theSplineCounts) delete[] theSplineCounts;
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}


	//
	// remove the temporary matrix file if function parameters define removal
	//
	if (DeleteBinaryFile)
	{
		Message("About to delete GDM Binary File");
		if ((_access(lpTmpFile, 0)) != -1) remove(lpTmpFile);
	}


	//
	// create a NULL model and return the deviance
	//
	double dNullDeviance = GetWeightedNULLDeviance(nTableRows, *pResponse, *pWeights);


	//
	// calculate deviance explained as a percentage
	//
	double dDevianceExplained = (1.0 - (dGDMDeviance / dNullDeviance)) * 100;


	//
	// write relevent outputs back to the parameter file
	//	
	//fptr("Updating parameter file...", 0);
	PrintGDMResults(pParams, dDevianceExplained, nMatrixCols, pCoefficients, false);
	UpdateParameterFile(pParams, dNullDeviance, dGDMDeviance, dDevianceExplained,
		                pCoefficients[0], &pCoefficients[1], (int)nMatrixCols - 1, theSplineCounts, (int)nTableRows, NULL);

	//
	// assign matrix dimensions
	//
	*NumRows = (int)nTableRows;
	*NumCols = (int)nMatrixCols;


	//
	// Clean up
	//
	if (lpTmpFile) delete[] lpTmpFile;
	if (pCoefficients) delete[] pCoefficients;
	if (pPredictorMatrix) delete[] pPredictorMatrix;
	if (theSplineCounts) delete[] theSplineCounts;
	if (theQuantiles) delete[] theQuantiles;
	return(true);
}





//
// Do a single GDM predictor significance test
//
double DoPredictorSigGDM(char *pParams, char *lpTmpFile, 
	                     double *pResponse, double *pWeights, double baselineDevExpl,
	                     int NumRows, int NumCols, 
	                     int Index, int NumIterations,
	                     int ColumnJump, int nPredSplines,
						 double *pMin, double *pMax,
	                     FPTR fptr, UPTR uptr)
{
	if (fptr != NULL)
		fptr("Performing GDM significance test...", 5);

	// change colour to processing (orange)
	//if (uptr)
	//	uptr(-5, -5, -5, -5, -5, -5, Index);

	double *pPredictorMatrix = new double[NumRows * NumCols];
	if (NULL == pPredictorMatrix)
	{
		Message("Cannot allocate pPredictorMatrix in DoPredictorSigGDM()", "ERROR");
		return(false);
	}

	// setup path to partial GDM binary matrix file for use in backward elimination
	char *lpSigFile = GetWorkspacePath(pParams);
	std::strcat(lpSigFile, "\\tmpPredSig.bin");

	//
	// Permute the data columns in the Predictor Matrix for this predictor only
	//
	std::vector<int> myVector;
	for (int i = 0; i < NumRows; i++)
	{
		myVector.push_back(i);
	}

	//
	// Perform the matrix regression with permuted data for this predictor
	//
	double *pShuffleData = new double[nPredSplines * NumRows];
	double *pIterDevExpl = new double[NumIterations];
	int nCurrent = 5;
	for (int i = 0; i < NumIterations; i++)
	{
		if (i * 100 / NumIterations)
		{
			nCurrent = i * 100 / NumIterations;
			if (fptr != NULL)
				fptr("Performing GDM significance test...", nCurrent);
		}


		//
		// Populate the Predictor Matrix from the binary image created in DoBaseSigGDM()
		//
		int h = _open(lpTmpFile, _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE);
		if (h < 0)
		{
			Message("Cannot open binary file for READ in DoPredictorSigGDM()", "ERROR");
			if (pPredictorMatrix) delete[] pPredictorMatrix;
			return(false);
		}


		//
		// create new binary matrix using intercept column followed by predictor columns
		//
		int hOut = _open(lpSigFile, _O_BINARY | _O_CREAT | _O_TRUNC | _O_RDWR, S_IREAD | S_IWRITE);
		if (hOut < 0)
		{
			Message("Cannot create binary file for in DoPredictorSigGDM()", "ERROR");
			if (lpSigFile) delete[] lpSigFile;
			return(-1.0);
		}


		//
		// copy to temporary significance test file while reading into predictor matrix
		//
		lseek(hOut, 0L, SEEK_SET);
		double *pTmp = pPredictorMatrix;
		for (int j = 0; j < NumCols; j++)
		{
			_read(h, pTmp, NumRows * sizeof(double));
			_write(hOut, pTmp, NumRows * sizeof(double));
			pTmp += NumRows;
		}
		_close(h);

		
		//
		// Permute the data columns in the Predictor Matrix for this predictor only
		//
		std::random_shuffle(myVector.begin(), myVector.end());
		for (int j = 0; j < nPredSplines; j++)
		{
			int nOffset = (ColumnJump + j) * NumRows;
			for (int k = 0; k < NumRows; k++)
			{
				pShuffleData[(j * NumRows) + k] = pPredictorMatrix[nOffset + myVector[k]];
			}
		}


		//
		// Copy permuted sub matrix to lpSigFile binary for use in the NNLS regression
		//		
		pTmp = &pPredictorMatrix[ColumnJump * NumRows];
		memcpy(pTmp, pShuffleData, nPredSplines * NumRows * sizeof(double));

		_lseek(hOut, ColumnJump * NumRows * sizeof(double), SEEK_SET);
		_write(hOut, pTmp, nPredSplines * sizeof(double));
		_close(hOut);

  
		//
		// Do the matrix regression
		//
		double dGDMDeviance;
		double *pCoefficients = NULL;
		if (NumIterations <= 10)
		{
			if (fptr != NULL)
				fptr("Performing GDM regression...", 5);
			pCoefficients = WeightedNNLSRegression(lpSigFile,
				                                   pPredictorMatrix,
				                                   NumRows,
				                                   NumCols,
				                                   pResponse,
				                                   &dGDMDeviance,
				                                   pWeights,
				                                   NULL);
		}
		else
		{
			pCoefficients = WeightedNNLSRegression(lpSigFile,
				                                   pPredictorMatrix,
				                                   NumRows,
				                                   NumCols,
				                                   pResponse,
				                                   &dGDMDeviance,
				                                   pWeights,
				                                   NULL);
		}

		if (NULL == pCoefficients)
		{
			Message("pCoefficients are NULL", "ERROR in DoPredictorSigGDM");
			if ((_access(lpSigFile, 0)) != -1) remove(lpSigFile);
			if (pPredictorMatrix) delete[] pPredictorMatrix;
			return(false);
		}


		//
		// create a NULL model and return the deviance
		//
		double dNullDeviance = GetWeightedNULLDeviance(NumRows, pResponse, pWeights);


		//
		// calculate deviance explained as a percentage
		//
		double dDevianceExplained = (1.0 - (dGDMDeviance / dNullDeviance)) * 100;
		pIterDevExpl[i] = dDevianceExplained;

	} // for (int i = 0; i < NumIterations; i++)


	// determine the proportion of randomised GDM Deviance Explained that exceed the baseline GdM
	int nExceed = 0;
	double dMin = 100.0;
	double dMax = -100.0;
	double dSum = 0.0;
	for (int i = 0; i < NumIterations; i++)
	{
		double dGap = pIterDevExpl[i] - (baselineDevExpl + d_EXCEED_EPSILON);

		if (dGap > 0.0)
		{
			++nExceed;

			if (dGap < dMin)
			{
				dMin = dGap;
			}

			if (dGap > dMax)
			{
				dMax = dGap;
			}

			dSum += dGap;
		}
	}


	// these will be used to update the parameter file
	if (nExceed == 0)
	{
		*pMin = 0.0;
		*pMax = 0.0;
	}
	else
	{
		*pMin = dMin;
		*pMax = dMax;
	}

	// change colour back to ready (green)
	if (fptr != NULL)
		fptr("Performing GDM significance test...", 100);
	//if (uptr)
	//	uptr(-10, -10, -10, -10, -10, -10, Index);

	//
	// Cleanup
	//
	if ((_access(lpSigFile, 0)) != -1) remove(lpSigFile);
	if (lpSigFile) delete[] lpSigFile;
	if (pShuffleData) delete[] pShuffleData;
	if (pIterDevExpl) delete[] pIterDevExpl;
	if (pPredictorMatrix) delete[] pPredictorMatrix;

	//
	// return significance
	//
	return(((double)nExceed) / NumIterations);
}



