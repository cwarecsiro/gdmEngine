//
// FilterPCT.cpp
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
#include "GDM_PredictorTypes.h"
#include "FilterPCT.h"

////////////////////////////////////////////////////////////////////////////////////
// ComparePCT comparison routine for Sort method (decending)
//
int ComparePCT(const void *arg1, const void *arg2)
{
	PCT_VEC *p1 = (PCT_VEC *)arg1;
	PCT_VEC *p2 = (PCT_VEC *)arg2;

	if (p1->dProb > p2->dProb)
		return(-1);
	else if (p1->dProb < p2->dProb)
		return(1);
	else
		return(0);
}
////////////////////////////////////////////////////////////////////////////////////



bool FilterPCTsFromR(char **ppPCTInputTable, char **ppWorkspace, int *pNumPCTs, int *pPCTs)
{
	char InPath[BUFF1024];
	strcpy(InPath, *ppPCTInputTable);  //Message(InPath, "InPath");
	int NumPCTs = *pNumPCTs;	       //Message(NumPCTs, "NumPCTs");
	int nRows = GetNumRows(InPath);    //Message(nRows, "nRows");
	int nCols = GetNumCols(InPath);    //Message(nCols, "nCols");


	//
	// Setup the PCT envelope grids
	//
	BinaryFileClass **ppFLT = new BinaryFileClass *[NumPCTs];
	EsriBinaryHeader **ppHDR = new EsriBinaryHeader *[NumPCTs];
	gmpath gmp;
	char buffHDR[BUFF1024];
	char buffFLT[BUFF1024];
	for (int i = 0; i < NumPCTs; i++)
	{
		sprintf(buffHDR, "%s\\PCT%d.hdr", *ppWorkspace, pPCTs[i]);    //Message(buffHDR, "PCT Paths.HDR");
		sprintf(buffFLT, "%s\\PCT%d.flt", *ppWorkspace, pPCTs[i]);    //Message(buffFLT, "PCT Paths.FLT");

		if ((gmp.FileExists(buffHDR)) && (gmp.FileExists(buffFLT)))   //Message(pPCTs[i], "EXISTS");
		{
			ppFLT[i] = new BinaryFileClass(buffFLT);
			ppHDR[i] = new EsriBinaryHeader(buffHDR);
		}
		else
		{
			Message(pPCTs[i], "DOES NOT EXIST");
			return(false);
		}
	}


	//
	// setup the projected PCT probabilty lookup table
	//
	double *pPCTprobs = new double[NumPCTs];
	PCT_VEC *myPCTVEC = new PCT_VEC[NumPCTs];

	//
	// process the table row by row...
	//
	char outPath[BUFF1024];
	char RowData[BUFFMEG];
	sprintf(outPath, "%s\\%s_VVV.csv", gmp.GetDirectoryPath(InPath), gmp.GetName(InPath));
	FILE *fpIn = fopen(InPath, "r+t");
	fgets(RowData, BUFFMEG, fpIn);
	FILE *fpOut = fopen(outPath, "w+t");
	//fprintf(fpOut, "ID,SC,X,Y,PCT\n");
	fprintf(fpOut, "PCT\n");

	char *seps = ",\n";
	while (1)
	{
		if (NULL == fgets(RowData, BUFFMEG, fpIn))
			break;

		// ID
		char *p = strtok(RowData, seps);
		int nID = atoi(p);

		// SC
		p = strtok(NULL, seps);
		int nSC = atoi(p);

		// X
		p = strtok(NULL, seps);
		double dX = atof(p);

		// Y
		p = strtok(NULL, seps);
		double dY = atof(p);

		// optional
		p = strtok(NULL, seps);

		// PCT probabilities
		for (int i = 0; i < NumPCTs; i++)
		{
			p = strtok(NULL, seps);
			pPCTprobs[i] = atof(p);
		}

		// Max Probability
		p = strtok(NULL, seps);
		double MaxBaseProb = atof(p);

		// best pct
		p = strtok(NULL, seps);
		int nBestPCT = atoi(p);


		//
		// Determine if the site is within the valid data for the Best PCT and if not,
		// add 10000 to the BestPCT value to signify that it is out of the envelope.
		//
		int nBaseIndex = getPCTIndex(nBestPCT);  // holds the index for the projected maximum likelyhood PCT
		long lOffset = GetGridOffset(dX, dY, ppHDR[nBaseIndex]);

		// if the site is NOT in the Maximum pProjected PCT Envelope...
		if (!((ppFLT[nBaseIndex])->CoordHasData(lOffset, ppHDR[nBaseIndex]->GetNoDataValue())))
		{
			// set this PCT to old PCT + 10000
			nBestPCT += 10000;

			// initialise the PCTVEC
			for (int i = 0; i < NumPCTs; i++)
			{
				myPCTVEC[i].PCT = pPCTs[i];
				myPCTVEC[i].dProb = pPCTprobs[i];
			}
			// and sort it
			qsort((void *)myPCTVEC, (size_t)NumPCTs, sizeof(PCT_VEC), ComparePCT);

			//
			// now we need to find the PCT for this site that (Inside the PCT envelope) AND (Has the highest projected probability except for the BestPCT)
			// The best PCT will be at the highest so we start from the second highest and go to the lowest while site in PCT
			//
			for (int i = 1; i < NumPCTs; i++)  
			{
				// get the PCT from this slot in the SORTED PCT List
				int tmpPCT = getPCTIndex(myPCTVEC[i].PCT);

				// if this site is inside the envelope for this PCT
				long tmpOffset = GetGridOffset(dX, dY, ppHDR[tmpPCT]);

				if (((ppFLT[tmpPCT])->CoordHasData(tmpOffset, ppHDR[tmpPCT]->GetNoDataValue())))
				{
					nBestPCT = myPCTVEC[i].PCT;
					break;
				}

			} // for (int i = 1; i < NumPCTs; i++)  

		}  // if the site is NOT in the Maximum pProjected PCT Envelope...
		
		// write output row
		fprintf(fpOut, "%d\n", nBestPCT);

	}
	fclose(fpIn);
	fclose(fpOut);

	//
	// cleanup
	//
	if (pPCTprobs) delete[] pPCTprobs;
	if (myPCTVEC) delete[] myPCTVEC;
	for (int i = 0; i < NumPCTs; i++)
	{
		if (ppFLT[i]) delete ppFLT[i];
		if (ppHDR[i]) delete ppHDR[i];
	}
	if (ppFLT) delete[] ppFLT;
	if (ppHDR) delete[] ppHDR;

	return(true);
}



//
// get the zero based binary file slot for SC1200
//
//int getPCTIndex(int pct)
//{
//	if (pct == 23)
//		return(0);
//	else if (pct == 26)
//		return(1);
//	else if (pct == 77)
//		return(2);
//
//	else
//	{
//		//Message("Error in getPCTIndex");
//		return(0);
//	}
//}

//
// get the zero based binary file slot for SC500
//
//int getPCTIndex(int pct)
//{
//	if (pct == 10)
//		return(0);
//
//	else if (pct == 11)
//		return(1);
//
//	else if (pct == 13)
//		return(2);
//
//	else if (pct == 15)
//		return(3);
//
//	else if (pct == 16)
//		return(4);
//
//	else if (pct == 245)
//		return(5);
//
//	else if (pct == 40)
//		return(6);
//
//	else if (pct == 74)
//		return(7);
//
//	else if (pct == 9)
//		return(8);
//
//	else
//	{
//		//Message("Error in getPCTIndex");
//		return(0);
//	}
//}


//
// get the zero based binary file slot for SC1100
//
//int getPCTIndex(int pct)
//{
//	if (pct == 2)
//		return(0);
//
//	else if (pct == 249)
//		return(1);
//
//	else if (pct == 278)
//		return(2);
//
//	else if (pct == 299)
//		return(3);
//
//	else if (pct == 302)
//		return(4);
//
//	else if (pct == 5)
//		return(5);
//
//	else if (pct == 7)
//		return(6);
//
//	else if (pct == 79)
//		return(7);
//
//	else if (pct == 8)
//		return(8);
//
//	else if (pct == 85)
//		return(9);
//
//	else
//	{
//		//Message("Error in getPCTIndex");
//		return(0);
//	}
//}


////
//// get the zero based binary file slot for SC700
////
//int getPCTIndex(int pct)
//{
//	if (pct == 12)
//		return(0);
//
//	else if (pct == 181)
//		return(1);
//
//	else if (pct == 182)
//		return(2);
//
//	else if (pct == 238)
//		return(3);
//
//	else if (pct == 24)
//		return(4);
//
//	else if (pct == 240)
//		return(5);
//
//	else if (pct == 335)
//		return(6);
//
//	else if (pct == 336)
//		return(7);
//
//	else if (pct == 47)
//		return(8);
//
//	else if (pct == 53)
//		return(9);
//
//	else
//	{
//		//Message("Error in getPCTIndex");
//		return(0);
//	}
//}


//
// get the zero based binary file slot for SC300
//
int getPCTIndex(int pct)
{
	if (pct == 1100)
		return(0);

	else if (pct == 1196)
		return(1);

	else if (pct == 300)
		return(2);

	else if (pct == 638)
		return(3);

	else if (pct == 639)
		return(4);

	else if (pct == 753)
		return(5);

	else
	{
		//Message("Error in getPCTIndex");
		return(0);
	}
}


//
// get gridcell offset (the number of cells away from the start of the file) based on X and Y coordinate
//
long GetGridOffset(double dX, double dY, EsriBinaryHeader *ebh)
{
	//
	// get mask grid metrics
	//
	int nCols = ebh->GetNumCols();
	double dSize = ebh->GetCellSize();
	double dMinX = ebh->GetMinX();
	double dMaxX = ebh->GetMaxX();
	double dMinY = ebh->GetMinY();
	double dMaxY = ebh->GetMaxY();
	long lOffset = -1L;

	//
	// check that the site is within the grid extent
	//
	if ((dX < dMinX) || (dX > dMaxX) || (dY < dMinY) || (dY > dMaxY))
	{
		return(lOffset);
	}

	else
	{
		// convert into zero based row and column indices
		int nX = int(floor((dX - dMinX) / dSize));
		int nY = int(floor((dMaxY - dY) / dSize));
		lOffset = (nY * nCols) + nX;
		return(lOffset);
	}
}


//
// Returns number of rows in text file less one for the expected header
//
int GetNumRows(char *InPath)
{
	char RowData[BUFFMEG];
	int nCount = 0;
	FILE *fp = fopen(InPath, "r+t");
	fgets(RowData, BUFFMEG, fp);  // get header
	while (1)
	{
		if (NULL == fgets(RowData, BUFFMEG, fp))
			break;
		++nCount;
	}
	fclose(fp);
	return(nCount);
}


//
// Returns number of rows in text file less one for the expected header
//
int GetNumCols(char *InPath)
{
	char RowData[BUFFMEG];
	FILE *fp = fopen(InPath, "r+t");
	fgets(RowData, BUFFMEG, fp);  // get header
	fclose(fp);
	int nCount = 0;
	char *p = strtok(RowData, ",\n");    // get the pointid
	if (p)
		++nCount;
	else
		return(nCount);
	while (1)
	{
		if (NULL == strtok(NULL, ",\n"))
			break;
		++nCount;
	}
	return(nCount);
}


//
// Line for Line copy of the PCT input table
//
bool DirectCopyPCTFile(char *InPath)
{
	char outPath[BUFF1024];
	char RowData[BUFFMEG];

	gmpath gmp;
	sprintf(outPath, "%s\\%s_ENV.csv", gmp.GetDirectoryPath(InPath), gmp.GetName(InPath));
	FILE *fpIn = fopen(InPath, "r+t");
	FILE *fpOut = fopen(outPath, "w+t");
	while (1)
	{
		if (NULL == fgets(RowData, BUFFMEG, fpIn))
			break;
		fputs(RowData, fpOut);
	}
	fclose(fpIn);
	fclose(fpOut);
	return(true);
}



//
// Line for Line copy to create ID,SC,X,Y,PCT table
//
bool CreateShapeTable(char *InPath, int NumPCTs)
{
	char outPath[BUFF1024];
	char RowData[BUFFMEG];

	gmpath gmp;
	sprintf(outPath, "%s\\%s_SHP.csv", gmp.GetDirectoryPath(InPath), gmp.GetName(InPath));
	FILE *fpIn = fopen(InPath, "r+t");
	fgets(RowData, BUFFMEG, fpIn);

	FILE *fpOut = fopen(outPath, "w+t");
	fprintf(fpOut, "ID,SC,X,Y,PCT\n");
	
	char *seps = ",\n";
	while (1)
	{
		if (NULL == fgets(RowData, BUFFMEG, fpIn))
			break;

		// ID
		char *p = strtok(RowData, seps);
		int nID = atoi(p);

		// SC
		p = strtok(NULL, seps);
		int nSC = atoi(p);

		// X
		p = strtok(NULL, seps);
		double dX = atof(p);

		// Y
		p = strtok(NULL, seps);
		double dY = atof(p);

		// optional
		p = strtok(NULL, seps);

		// PCT probabilities
		for (int i = 0; i < NumPCTs; i++)
		{
			p = strtok(NULL, seps);
		}

		// Max Probability
		p = strtok(NULL, seps);

		// best pct
		p = strtok(NULL, seps);
		int nBestPCT = atoi(p);

		// write output row
		fprintf(fpOut, "%d,%d,%lf,%lf,%d\n", nID, nSC, dX, dY, nBestPCT);

	}
	fclose(fpIn);
	fclose(fpOut);
	return(true);
}