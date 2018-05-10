//
// GDMTimeSeriesDissimilarity.cpp
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
//#include "ann.h" 
#include "clsANNTree.h" 


float **CreateGDMTransformMatrix( char *PointFile, int *pRows, int *pCols );

float **CreateDemandSiteDataFromCSVTableV5( char *PointFile );
ANNpointArray CreateDemandSiteDataFromCSVTableANN( char *PointFile );
int GetDemandSiteColumnCount(char *PointFile);
int GetDemandSiteCount(char *PointFile);

float CalculateBCSimilarity(float **ppData, int x, int y, int nCols);
float CalculateBCDissimilarity(float **ppData, int x, int y, int nCols);


bool DoTimeSeriesDissimilarityGrid(char *pStartParams,
	                               char *pEndParams,
                                   char *OutPath,
								   FPTR fptr)
{
	//Message("At entry to DoTimeSeriesDissimilarityGrid", "INFO");
	//Message(pStartParams, "pStartParams");
	//Message(pEndParams, "pEndParams");
	//Message(OutPath, "OutPath");
	//return(false);

	//
	// setup the grid stacks
	//
	BinaryFileStack *bfs0 = new BinaryFileStack(pStartParams, BFS_TranGrids);
	if (NULL == bfs0)
	{
		Message("Cannot construct BinaryFileStack bfs0", "ERROR");
		return(false);
	}
	BinaryFileStack *bfs1 = new BinaryFileStack(pEndParams, BFS_TranGrids);
	if (NULL == bfs0)
	{
		Message("Cannot construct BinaryFileStack bfs1", "ERROR");
		return(false);
	}


	//
	// create the output grid
	//
	gmpath gmPath;
	BinaryFileClass *bfc_Out = new BinaryFileClass(gmPath.ChangeExtension(OutPath, ".flt"), BFC_CreateOrTrunc);	

	float *pOutRow = new float [bfs0->GetNumCols()];
	int nCurrent = 0;
	fptr("Creating dissimilarity Grid...", nCurrent);
	for ( int i=0; i<bfs0->GetNumRows(); i++ )
	{
		//
		// update progress
		//
		if ( i * 100 / bfs0->GetNumRows() > nCurrent)
		{
			nCurrent = i * 100 / bfs0->GetNumRows();
			fptr("Creating dissimilarity Grid...", nCurrent);
		}

		//
		// get the dissimilarities
		//
		bfs0->GetNextRow();
		bfs1->GetNextRow();
		
		for (int j=0; j<bfs0->GetNumCols(); j++ )
		{
			if (bfs0->GetRowData()[0][j] != bfs0->GetNoDataVal())
			{
				bfs0->GetCellAt(j);
				bfs1->GetCellAt(j);
				pOutRow[j] = (float)bfs0->GetGDMDissimilarity(bfs1->GetCellData());
			}
			else
			{
				pOutRow[j] = bfs0->GetNoDataVal();
			}
		}

		//
		// write the output row
		//
		bfc_Out->WriteFloat(pOutRow, bfs0->GetNumCols());
	}
	
	//
	// Clean up
	//
	bfc_Out->Close();
	bfs0->CopyHeaderTo(gmPath.ChangeExtension(OutPath, ".hdr"));
	if (pOutRow) delete[] pOutRow;
	if (bfc_Out) delete bfc_Out;
	if (bfs0) delete bfs0;
	if (bfs1) delete bfs1;
	nCurrent = 0;
	fptr("Status: Ready", nCurrent);
	return(true);
}




bool DoTimeSeriesDissimilarityGridAsBlocks(char *pStartParams,
	                                       char *pEndParams,
                                           char *OutPath,
										   int nBlockSize,
								           FPTR fptr)
{
	//Message("At entry to DoTimeSeriesDissimilarityGridAsBlocks", "INFO");
	//Message(pStartParams, "pStartParams");
	//Message(pEndParams, "pEndParams");
	//Message(OutPath, "OutPath");
	//Message(nBlockSize, "nBlockSize");
	//return(false);

	//
	// setup the grid stacks
	//
	BinaryFileBlockStack *bfs0 = new BinaryFileBlockStack(pStartParams, BFS_TranGrids, nBlockSize);
	if (NULL == bfs0)
	{
		Message("Cannot construct BinaryFileStack bfs0", "ERROR");
		return(false);
	}
	BinaryFileBlockStack *bfs1 = new BinaryFileBlockStack(pEndParams, BFS_TranGrids, nBlockSize);
	if (NULL == bfs1)
	{
		Message("Cannot construct BinaryFileStack bfs1", "ERROR");
		return(false);
	}


	//
	// create the output grid
	//
	gmpath gmPath;
	BinaryFileClass *bfc_Out = new BinaryFileClass(gmPath.ChangeExtension(OutPath, ".flt"), BFC_CreateOrTrunc);	
	float *pOutRow = new float [bfs0->GetNumCols()];

	int nThisRow = 0;
	int nBlocks = bfs0->GetNumBlocks();
	int nRemainder = bfs0->GetBlockRemainder();
	int nRows = bfs0->GetNumRows();
	int nCols = bfs0->GetNumCols();

	//
	// do the main blocks...
	//
	int nCurrent = 0;
	fptr("Creating dissimilarity Grid...", nCurrent);
	for ( int i=0; i<nBlocks; i++ )
	{
		//
		// update progress
		//
		if ( i * 100 / (nBlocks-1) > nCurrent)
		{
			nCurrent = i * 100 / (nBlocks-1);
			fptr("Creating dissimilarity Grid...", nCurrent);
		}

		bfs0->GetNextRow();
		bfs1->GetNextRow();

		for ( int j=0; j<nBlockSize; j++ )
		{
			for (int k=0; k<nCols; k++ )
			{
				if (bfs0->GetRowData()[0][(j*nCols)+k] != bfs0->GetNoDataVal())
				{
					bfs0->GetCellAt(j, k);
					bfs1->GetCellAt(j, k);
					pOutRow[k] = (float)bfs0->GetGDMDissimilarity(bfs1->GetCellData());
					//pOutRow[k] = 1.0F;
				}
				else
				{
					pOutRow[k] = bfs0->GetNoDataVal();
				}
			}

			//
			// write the output row
			//
			bfc_Out->WriteFloat(pOutRow, nCols);
		}
	}

	//
	// finally, do the remainder...
	//
	bfs0->GetNextRow();
	bfs1->GetNextRow();
	fptr("Creating dissimilarity Grid...", 90);
	for ( int j=0; j<nRemainder; j++ )
	{
		for (int k=0; k<nCols; k++ )
		{
			if (bfs0->GetRowData()[0][(j*nCols)+k] != bfs0->GetNoDataVal())
			{
				bfs0->GetCellAt(j, k);
				bfs1->GetCellAt(j, k);
				pOutRow[k] = (float)bfs0->GetGDMDissimilarity(bfs1->GetCellData());
				//pOutRow[k] = 1.0F;
			}
			else
			{
				pOutRow[k] = bfs0->GetNoDataVal();
			}
		}

		//
		// write the output row
		//
		bfc_Out->WriteFloat(pOutRow, nCols);
	}
	
	//
	// Clean up
	//
	fptr("Creating dissimilarity Grid...", 95);
	bfc_Out->Close();
	bfs0->CopyHeaderTo(gmPath.ChangeExtension(OutPath, ".hdr"));
	fptr("Creating dissimilarity Grid...", 100);
	if (pOutRow) delete[] pOutRow;
	if (bfc_Out) delete bfc_Out;
	if (bfs0) delete bfs0;
	if (bfs1) delete bfs1;
	nCurrent = 0;
	fptr("Status: Ready", nCurrent);
	return(true);
}


//
// Implimentation of Dan Faith's Climate Change Impact Approach 
// based on the sum of minimum distances between demand points and ALL valid gridcells
//
bool DoTimeSeriesImpactTable(char *pDemandPointPath,
	                         char *pStartParams,
	                         char *pEndParams,
                             char *OutPath,
							 FPTR fptr)
{
	Message("At entry to DoTimeSeriesImpactTable", "INFO");

	//
	// Fill in the demand point data frame
	//
	int nDemandSites = GetDemandSiteCount(pDemandPointPath);
	//Message(nDemandSites, "info");

	float **ppDemandData = CreateDemandSiteDataFromCSVTableV5(pDemandPointPath);
	if (NULL == ppDemandData)
	{
		Message("Cannot construct ppDemandData", "ERROR");
		return(false);
	}


	//
	// setup the grid stacks
	//
	BinaryFileStack *bfs0 = new BinaryFileStack(pStartParams, BFS_TranGrids);
	if (NULL == bfs0)
	{
		Message("Cannot construct BinaryFileStack bfs0", "ERROR");
		return(false);
	}
	BinaryFileStack *bfs1 = new BinaryFileStack(pEndParams, BFS_TranGrids);
	if (NULL == bfs0)
	{
		Message("Cannot construct BinaryFileStack bfs1", "ERROR");
		return(false);
	}


	//
	// create the demand class min-value vectors
	//	
	float *pMinRow0 = new float [nDemandSites];
	float *pMinRow1 = new float [nDemandSites];
	for ( int i=0; i<nDemandSites; i++ )
	{
		pMinRow0[i] = pMinRow1[i] = FLT_MAX; 
	}


	//
	// Run the analysis
	//
	int nCurrent = 0;
	fptr("Creating CC Impact Table...", nCurrent);
	for ( int i=0; i<bfs0->GetNumRows(); i++ )
	{
		//
		// update progress
		//
		if ( i * 100 / bfs0->GetNumRows() > nCurrent)
		{
			nCurrent = i * 100 / bfs0->GetNumRows();
			char banner[256];
			sprintf(banner, "Creating CC Impact Table...Row %d of %d", i, bfs0->GetNumRows());
			fptr(banner, nCurrent);
		}

		// get the dissimilarities
		bfs0->GetNextRow();
		bfs1->GetNextRow();
		for (int j=0; j<bfs0->GetNumCols(); j++ )
		{
			if ((bfs0->GetRowData()[0][j] != bfs0->GetNoDataVal()) && (bfs1->GetRowData()[0][j] != bfs1->GetNoDataVal()))
			{
				bfs0->GetCellAt(j);
				bfs1->GetCellAt(j);

				/////////////////////////////////////////////////////////////////////////////////////
				//
				// This is the MINIMUM dissimilarity version (as per Dan Faith's original calc)
				//
				//
				// work thru the demand sites and get the minimum ED (expressed as a Manhattan distance) for each epoch
				for ( int k=0; k<nDemandSites; k++ )
				{
					float f0 = (float)bfs0->GetManhattanDist(ppDemandData[k]);
					float f1 = (float)bfs1->GetManhattanDist(ppDemandData[k]);   

					if ( f0 < pMinRow0[k] ) pMinRow0[k] = f0;
					if ( f1 < pMinRow1[k] ) pMinRow1[k] = f1;
				}
			}
		}
	}
	

	//
	// write the output table and calculate the proportion of species loss/gain between epochs
	//
	float fSum_0 = 0.0F;
	float fSum_1 = 0.0F;
	gmpath gmPath;
	FILE *fp = fopen(gmPath.ChangeExtension(OutPath, ".csv"), "w+t");	
	fprintf(fp, "Demand_Point,Start_dissim,EndDissim,Difference\n");
	for ( int i=0; i<nDemandSites; i++ )
	{
		fprintf(fp, "%d,%f,%f,%f\n", 
			    i+1, 
				1.0-exp(-pMinRow0[i]), 
			    1.0-exp(-pMinRow1[i]), 
				(1.0-exp(-pMinRow1[i])) - (1.0-exp(-pMinRow0[i])));

		fSum_0 += (float)1.0-exp(-pMinRow0[i]);
		fSum_1 += (float)1.0-exp(-pMinRow1[i]);
	}
	fprintf(fp, "Totals:,%f,%f,%f\n", fSum_0, fSum_1, fSum_1 - fSum_0);
	fprintf(fp, "   ,   ,Percent Lost:,%f\n", (fSum_1 - fSum_0) * 100 / nDemandSites);
	if (fp) fclose(fp);

	//
	// Clean up
	//
	if (pMinRow0) delete[] pMinRow0;
	if (pMinRow1) delete[] pMinRow1;
	for (int i=0; i<nDemandSites; i++ )
	{
		if (ppDemandData[i]) delete[] ppDemandData[i];
	}
	if (ppDemandData) delete[] ppDemandData;
	if (bfs0) delete bfs0;
	if (bfs1) delete bfs1;
	nCurrent = 0;
	fptr("Status: Ready", nCurrent);
	return(true);
}


//
// Create a GDM Transform Similarity Table from a K-Means cluster centroid Table
//
bool CreateImpactSimilarityTable(char *pDemandPointPath, char *OutPath, FPTR fptr)
{
	//Message("At entry to CreateImpactSimilarityTable", "INFO");

	//
	// get the Transform Data
	//
	int nRows = 0;
	int nCols = 0;
	float **ppData = CreateGDMTransformMatrix(pDemandPointPath, &nRows, &nCols);
	//Message(nRows, "nRows");
	//Message(nCols, "nCols");


	//
	// create the similarity matrix
	//
	float **ppSimTable = new float * [nRows]; 
	for (int i=0; i<nRows; i++)
	{
		ppSimTable[i] = new float [nRows];
		for (int j=0; j<nRows; j++)
		{
			if (i==j)
				ppSimTable[i][j] = 1.0F;
			else
				ppSimTable[i][j] = 0.0F;
		}
	}


	//
	// calculate bray-curtis similarities
	//
	for ( int i=1; i<nRows; i++)
	{
		for ( int j=0; j<i; j++ )
		{
			ppSimTable[i][j] = ppSimTable[j][i] = CalculateBCSimilarity(ppData, i, j, nCols);
		}
	}


	//
	// Write the output table file
	//
	FILE *fp = fopen(OutPath, "w+t");
	// write header row
	for ( int i=0; i<nRows; i++ ) fprintf(fp, "%d,", i);
	fprintf(fp, "%d\n", nRows);

	// write data
	for ( int i=0; i<nRows; i++ )
	{
		fprintf(fp, "%d,", i+1); // write row index

		for ( int j=0; j<nRows; j++ )
		{
			fprintf(fp, "%f", ppSimTable[i][j]);
			if ( j<nRows-1)
				fprintf(fp, ",");
			else
				fprintf(fp, "\n");
		}
	}
	fclose(fp);
	

	//
	// cleanup
	//
	if (ppData)
	{
		for ( int i=0; i<nRows; i++ )
		{
			if (ppData[i]) delete[] ppData[i];
			if (ppSimTable[i]) delete[] ppSimTable[i];
		}
		delete[] ppData;
		if (ppSimTable) delete[] ppSimTable;
	}

	return(true);
}


//
// Create a GDM Transform Dissimilarity Table from a K-Means cluster centroid Table
//
bool CreateImpactDissimilarityTable(char *pDemandPointPath, char *OutPath, FPTR fptr)
{
	//Message("At entry to CreateImpactDissimilarityTable", "INFO");
	
	//
	// get the Transform Data
	//
	int nRows = 0;
	int nCols = 0;
	float **ppData = CreateGDMTransformMatrix(pDemandPointPath, &nRows, &nCols);
	//Message(nRows, "nRows");
	//Message(nCols, "nCols");

	//
	// create the dissimilarity matrix
	//
	float **ppDissimTable = new float * [nRows]; 
	for (int i=0; i<nRows; i++)
	{
		ppDissimTable[i] = new float [nRows];
		for (int j=0; j<nRows; j++)
		{
			if (i==j)
				ppDissimTable[i][j] = 0.0F;
			else
				ppDissimTable[i][j] = 1.0F;
		}
	}


	//
	// calculate bray-curtis dissimilarities
	//
	for ( int i=1; i<nRows; i++)
	{
		for ( int j=0; j<i; j++ )
		{
			ppDissimTable[i][j] = ppDissimTable[j][i] = CalculateBCDissimilarity(ppData, i, j, nCols);
		}
	}


	//
	// Write the output table file
	//
	FILE *fp = fopen(OutPath, "w+t");
	
	// write header row
	for ( int i=0; i<nRows; i++ ) fprintf(fp, "%d,", i);
	fprintf(fp, "%d\n", nRows);

	// write data
	for ( int i=0; i<nRows; i++ )
	{
		fprintf(fp, "%d,", i+1); // write row index

		for ( int j=0; j<nRows; j++ )
		{
			fprintf(fp, "%f", ppDissimTable[i][j]);
			if ( j<nRows-1)
				fprintf(fp, ",");
			else
				fprintf(fp, "\n");
		}
	}
	fclose(fp);
	

	//
	// cleanup
	//
	if (ppData)
	{
		for ( int i=0; i<nRows; i++ )
		{
			if (ppData[i]) delete[] ppData[i];
			if (ppDissimTable[i]) delete[] ppDissimTable[i];
		}
		delete[] ppData;
		if (ppDissimTable) delete[] ppDissimTable;
	}

	return(true);
}


//
// Create and populate the GDM Transform Data matrix
//
float **CreateGDMTransformMatrix( char *PointFile, int *pRows, int *pCols )
{
	//
	// leading cols are ID,X,Y so there MUST be 3 more columns than there are predictors
	//
	int nCols = GetDemandSiteColumnCount(PointFile);
	nCols -= 3; // so we only have the predictor count
	*pCols = nCols;
	//Message(nCols, "nCols");

	int nRows = GetDemandSiteCount(PointFile);
	*pRows = nRows;
	//Message(nRows, "nRows");


	//
	// allocate the ANN array
	//
	float **ppDemandData = new float *[nRows];
	for (int i=0; i<nRows; i++ )
	{
		ppDemandData[i] = new float [nCols];		
	}
	
	if ( NULL == ppDemandData )
	{
		Message("Could not allocate GDM Transform Data", "CreateGDMTransformMatrix" );
		return(NULL);
	}


	//
	// open the demand point file
	//
	FILE *fp = fopen( PointFile , "r+t" );

	// get the header
	char buff[BUFFLEN];
	fgets( buff, BUFFLEN, fp );

	int nThisRow = 0;
	while(1)
	{
		if ( NULL == fgets( buff, BUFFLEN, fp ) )
			break;

		char *p = strtok(buff, ",\n"); // get the ID
		p = strtok(NULL, ",\n");       // get X
		p = strtok(NULL, ",\n");       // get Y

		for ( int j=0; j<nCols; j++ ) // get the predictor values
		{
			p = strtok(NULL, ",\n");
			ppDemandData[nThisRow][j] = (float)atof(p);
		}

		++nThisRow;
	}

	// close file
	fclose( fp );
	return(ppDemandData);
}


//
// return Bray-Curtis Similarity for a pair of vectors
//
float CalculateBCSimilarity(float **ppData, int x, int y, int nCols)
{
	float fVal = 0.0F;
	for ( int i=0; i<nCols; i++ )
	{
		fVal += fabs(ppData[x][i] - ppData[y][i]);
	}
	return((float)exp(-fVal));
}


//
// return Bray-Curtis Dissimilarity for a pair of vectors
//
float CalculateBCDissimilarity(float **ppData, int x, int y, int nCols)
{
	float fVal = 0.0F;
	for ( int i=0; i<nCols; i++ )
	{
		fVal += fabs(ppData[x][i] - ppData[y][i]);
	}
	return((float)(1.0F - exp(-fVal)));
}


////
//// Implimentation of Dan Faith's Climate Change Impact Approach
////
//bool DoTimeSeriesImpactGrids(char *pDemandPointPath,
//	                         char *pStartParams,
//	                         char *pEndParams,
//                             char *OutPathStart,
//							 char *OutPathEnd,
//							 FPTR fptr)
//{
//	//Message("At entry to DoTimeSeriesImpactGrids", "INFO");
//	//Message(pDemandPointPath, "pDemandPointPath");
//	//Message(pStartParams, "pStartParams");
//	//Message(pEndParams, "pEndParams");
//	//Message(OutPathStart, "OutPathStart");
//	//Message(OutPathEnd, "OutPathEnd");
//	//return(false);
//
//	//
//	// Fill in the demand point data frame
//	//
//	int nDemandSites = GetDemandSiteCount(pDemandPointPath);
//	//Message(nDemandSites, "info");
//
//	float **ppDemandData = CreateDemandSiteDataFromCSVTableV5(pDemandPointPath);
//	if (NULL == ppDemandData)
//	{
//		Message("Cannot construct ppDemandData", "ERROR");
//		return(false);
//	}
//
//
//	//
//	// setup the grid stacks
//	//
//	BinaryFileStack *bfs0 = new BinaryFileStack(pStartParams, BFS_TranGrids);
//	if (NULL == bfs0)
//	{
//		Message("Cannot construct BinaryFileStack bfs0", "ERROR");
//		return(false);
//	}
//	BinaryFileStack *bfs1 = new BinaryFileStack(pEndParams, BFS_TranGrids);
//	if (NULL == bfs0)
//	{
//		Message("Cannot construct BinaryFileStack bfs1", "ERROR");
//		return(false);
//	}
//
//
//	//
//	// create the output grids
//	//
//	gmpath gmPath;
//	//BinaryFileClass *bfc_OutStart = new BinaryFileClass(gmPath.ChangeExtension(OutPathStart, ".flt"), BFC_CreateOrTrunc);	
//	//BinaryFileClass *bfc_OutEnd = new BinaryFileClass(gmPath.ChangeExtension(OutPathEnd, ".flt"), BFC_CreateOrTrunc);	
//
//
//	FILE *fp0 = fopen(gmPath.ChangeExtension(OutPathStart, ".csv"), "w+t");	
//	FILE *fp1 = fopen(gmPath.ChangeExtension(OutPathEnd, ".csv"), "w+t");
//	fprintf(fp0, "Point,Dissimilarity\n");
//	fprintf(fp1, "Point,Dissimilarity\n");
//
//
//	//float *pOutRow0 = new float [bfs0->GetNumCols()];
//	//float *pOutRow1 = new float [bfs1->GetNumCols()];
//	int nCurrent = 0;
//	fptr("Creating CC Impact Grid...", nCurrent);
//	for ( int i=0; i<bfs0->GetNumRows(); i++ )
//	{
//		//
//		// update progress
//		//
//		if ( i * 100 / bfs0->GetNumRows() > nCurrent)
//		{
//			nCurrent = i * 100 / bfs0->GetNumRows();
//			char banner[256];
//			sprintf(banner, "Creating CC Impact Grid...Row %d of %d", i, bfs0->GetNumRows());
//			fptr(banner, nCurrent);
//		}
//
//		//
//		// get the dissimilarities
//		//
//		//bfs0->GetNextRow();
//		//bfs1->GetNextRow();
//		
//		if ( i == 6368)
//		{
//			bfs0->GetRowAt(i);
//			bfs1->GetRowAt(i);
//
//			bfs0->GetCellAt(9908);
//			bfs1->GetCellAt(9908);
//
//			for ( int k=0; k<nDemandSites; k++ )
//			{
//				float f0 = (float)bfs0->GetGDMDissimilarity(ppDemandData[k]);
//				float f1 = (float)bfs1->GetGDMDissimilarity(ppDemandData[k]);
//				fprintf(fp0, "%d,%f\n", k+1, f0);
//				fprintf(fp1, "%d,%f\n", k+1, f1);
//			}
//
//
//			//for (int j=0; j<bfs0->GetNumCols(); j++ )
//			//{
//			//	if ((bfs0->GetRowData()[0][j] != bfs0->GetNoDataVal()) && (bfs1->GetRowData()[0][j] != bfs1->GetNoDataVal()))
//			//	{
//			//		bfs0->GetCellAt(j);
//			//		bfs1->GetCellAt(j);
//
//			//		/////////////////////////////////////////////////////////////////////////////////////
//			//		//
//			//		// This is the MINIMUM dissimilarity version (as per Dan Faith's original calc)
//			//		//
//			//		//
//			//		// work thru the demand sites and get the minimum dissimilarity for each epoch
//			//		float fMin0 = 1.0F;
//			//		float fMin1 = 1.0F;
//			//		for ( int k=0; k<nDemandSites; k++ )
//			//		{
//			//			float f0 = (float)bfs0->GetGDMDissimilarity(ppDemandData[k]);
//			//			float f1 = (float)bfs1->GetGDMDissimilarity(ppDemandData[k]);
//
//			//			if ( f0 < fMin0 ) fMin0 = f0;
//			//			if ( f1 < fMin1 ) fMin1 = f1;
//			//		}
//			//		pOutRow0[j] = fMin0;
//			//		pOutRow1[j] = fMin1;
//			//	}
//			//	else
//			//	{
//			//		pOutRow0[j] = bfs0->GetNoDataVal();
//			//		pOutRow1[j] = bfs1->GetNoDataVal();
//			//	}
//			//}
//		}
//
//		//
//		// write the output row
//		//
//		//bfc_OutStart->WriteFloat(pOutRow0, bfs0->GetNumCols());
//		//bfc_OutEnd->WriteFloat(pOutRow1, bfs1->GetNumCols());
//	}
//	
//
//	if (fp0) fclose(fp0);
//	if (fp1) fclose(fp1);
//
//
//	//
//	// Clean up
//	//
//	for (int i=0; i<nDemandSites; i++ )
//	{
//		if (ppDemandData[i]) delete[] ppDemandData[i];
//	}
//	if (ppDemandData) delete[] ppDemandData;
//	//bfc_OutStart->Close();
//	//bfc_OutEnd->Close();
//	bfs0->CopyHeaderTo(gmPath.ChangeExtension(OutPathStart, ".hdr"));
//	bfs0->CopyHeaderTo(gmPath.ChangeExtension(OutPathEnd, ".hdr"));
//	//if (pOutRow0) delete[] pOutRow0;
//	//if (pOutRow1) delete[] pOutRow1;
//	//if (bfc_OutStart) delete bfc_OutStart;
//	//if (bfc_OutEnd) delete bfc_OutEnd;
//	if (bfs0) delete bfs0;
//	if (bfs1) delete bfs1;
//	nCurrent = 0;
//	fptr("Status: Ready", nCurrent);
//	return(true);
//}
//
//
////
//// Implimentation of Dan Faith's Climate Change Impact Approach
////
//bool DoTimeSeriesImpactGridsANN(char *pDemandPointPath,
//	                            char *pStartParams,
//	                            char *pEndParams,
//                                char *OutPathStart,
//							    char *OutPathEnd,
//							    FPTR fptr)
//{
//	Message("At entry to DoTimeSeriesImpactGridsANN", "INFO");
//
//	//
//	// Fill in the demand point data frame
//	//
//	int nDemandSites = GetDemandSiteCount(pDemandPointPath);
//	//Message(nDemandSites, "info");
//
//	int nCols = GetDemandSiteColumnCount(pDemandPointPath);
//	nCols -= 3; // so we only have the predictor count
//
//	ANNpointArray ppDemandData = CreateDemandSiteDataFromCSVTableANN(pDemandPointPath);
//	if (NULL == ppDemandData)
//	{
//		Message("Cannot construct ppDemandData", "ERROR");
//		return(false);
//	}
//
//	//
//	// now build the kd tree search structure with 
//	// the data points, number of data points and dimension of search space
//	//
//	ANNkd_tree *theTree = new ANNkd_tree ( ppDemandData, nDemandSites, nCols ); 
//
//
//	//
//	// find the nearest neighbour classes
//	//
//	// 
//	// setup some structures for the queries to follow
//	ANNidxArray nn_idx = new ANNidx[1];
//	ANNdistArray dist = new ANNdist[1];
//	ANNpoint query_pt0 = annAllocPt( nCols );
//	ANNpoint query_pt1 = annAllocPt( nCols );
//
//
//	//
//	// setup the grid stacks
//	//
//	BinaryFileStack *bfs0 = new BinaryFileStack(pStartParams, BFS_TranGrids);
//	if (NULL == bfs0)
//	{
//		Message("Cannot construct BinaryFileStack bfs0", "ERROR");
//		return(false);
//	}
//	BinaryFileStack *bfs1 = new BinaryFileStack(pEndParams, BFS_TranGrids);
//	if (NULL == bfs0)
//	{
//		Message("Cannot construct BinaryFileStack bfs1", "ERROR");
//		return(false);
//	}
//
//
//	//
//	// create the output grids
//	//
//	gmpath gmPath;
//	BinaryFileClass *bfc_OutStart = new BinaryFileClass(gmPath.ChangeExtension(OutPathStart, ".flt"), BFC_CreateOrTrunc);	
//	BinaryFileClass *bfc_OutEnd = new BinaryFileClass(gmPath.ChangeExtension(OutPathEnd, ".flt"), BFC_CreateOrTrunc);	
//
//	float *pOutRow0 = new float [bfs0->GetNumCols()];
//	float *pOutRow1 = new float [bfs1->GetNumCols()];
//	int nCurrent = 0;
//	fptr("Creating CC Impact Grid...", nCurrent);
//	for ( int i=0; i<bfs0->GetNumRows(); i++ )
//	{
//		//
//		// update progress
//		//
//		if ( i * 100 / bfs0->GetNumRows() > nCurrent)
//		{
//			nCurrent = i * 100 / bfs0->GetNumRows();
//			char banner[256];
//			sprintf(banner, "Creating CC Impact Grid...Row %d of %d", i, bfs0->GetNumRows());
//			fptr(banner, nCurrent);
//		}
//
//		//
//		// get the dissimilarities
//		//
//		bfs0->GetNextRow();
//		bfs1->GetNextRow();
//		
//		for (int j=0; j<bfs0->GetNumCols(); j++ )
//		{
//			if ((bfs0->GetRowData()[0][j] != bfs0->GetNoDataVal()) && (bfs1->GetRowData()[0][j] != bfs1->GetNoDataVal()))
//			{
//				bfs0->GetCellAt(j);
//				bfs1->GetCellAt(j);
//
//				for ( int k=0; k<nCols; k++ )
//				{
//					query_pt0[k] = bfs0->GetCellData()[k];
//					query_pt1[k] = bfs1->GetCellData()[k];
//				}
//
//				//
//				// find the nearest neighbour class in the ANN tree
//				// to get the minimum dissimilarity for each epoch
//				//
//				theTree->annkSearch( query_pt0, 1, nn_idx, dist, 0 );	
//				pOutRow0[j] = (float)(1.0F - exp(-dist[nn_idx[0]]));   // apply the GDM link function
//
//				theTree->annkSearch( query_pt1, 1, nn_idx, dist, 0 );	
//				pOutRow1[j] = (float)(1.0F - exp(-dist[nn_idx[0]]));   // apply the GDM link function
//			}
//			else
//			{
//				pOutRow0[j] = bfs0->GetNoDataVal();
//				pOutRow1[j] = bfs1->GetNoDataVal();
//			}
//		}
//
//		//
//		// write the output row
//		//
//		bfc_OutStart->WriteFloat(pOutRow0, bfs0->GetNumCols());
//		bfc_OutEnd->WriteFloat(pOutRow1, bfs1->GetNumCols());
//	}
//	
//	//
//	// Clean up
//	//
//	if (theTree) delete theTree;
//	bfc_OutStart->Close();
//	bfc_OutEnd->Close();
//	bfs0->CopyHeaderTo(gmPath.ChangeExtension(OutPathStart, ".hdr"));
//	bfs0->CopyHeaderTo(gmPath.ChangeExtension(OutPathEnd, ".hdr"));
//	if (pOutRow0) delete[] pOutRow0;
//	if (pOutRow1) delete[] pOutRow1;
//	if (bfc_OutStart) delete bfc_OutStart;
//	if (bfc_OutEnd) delete bfc_OutEnd;
//	if (bfs0) delete bfs0;
//	if (bfs1) delete bfs1;
//	nCurrent = 0;
//	fptr("Status: Ready", nCurrent);
//	return(true);
//}


//
// Create and populate the demand site ANN data matrix
//


ANNpointArray CreateDemandSiteDataFromCSVTableANN( char *PointFile )
{
	//
	// leading cols are ID,X,Y so there MUST be 3 more columns than there are predictors
	//
	int nCols = GetDemandSiteColumnCount(PointFile);
	nCols -= 3; // so we only have the predictor count
	//Message(nCols, "nCols");

	int nRows = GetDemandSiteCount(PointFile);
	//Message(nRows, "nRows");

	//
	// allocate the ANN array
	//
	ANNpointArray ppDemandData = annAllocPts( nRows, nCols );
	
	//
	// open the demand point file
	//
	FILE *fp = fopen( PointFile , "r+t" );

	// get the header
	char buff[BUFFLEN];
	fgets( buff, BUFFLEN, fp );

	int nThisRow = 0;
	while(1)
	{
		if ( NULL == fgets( buff, BUFFLEN, fp ) )
			break;

		char *p = strtok(buff, ",\n"); // get the ID
		p = strtok(NULL, ",\n");       // get X
		p = strtok(NULL, ",\n");       // get Y

		for ( int j=0; j<nCols; j++ ) // get the predictor values
		{
			p = strtok(NULL, ",\n");
			ppDemandData[nThisRow][j] = atof(p);
		}

		++nThisRow;
	}

	// close file
	fclose( fp );
	return(ppDemandData);
}



//
// Create and populate the demand site data matrix
//
float **CreateDemandSiteDataFromCSVTableV5( char *PointFile )
{
	//
	// leading cols are ID,X,Y so there MUST be 3 more columns than there are predictors
	//
	int nCols = GetDemandSiteColumnCount(PointFile);
	nCols -= 3; // so we only have the predictor count
	//Message(nCols, "nCols");

	int nRows = GetDemandSiteCount(PointFile);
	//Message(nRows, "nRows");


	//
	// allocate the ANN array
	//
	float **ppDemandData = new float *[nRows];
	for (int i=0; i<nRows; i++ )
	{
		ppDemandData[i] = new float [nCols];		
	}
	
	if ( NULL == ppDemandData )
	{
		Message("Could not allocate demandData", "CreateDemandSiteDataFromCSVTableV5" );
		return(NULL);
	}


	//
	// open the demand point file
	//
	FILE *fp = fopen( PointFile , "r+t" );

	// get the header
	char buff[BUFFLEN];
	fgets( buff, BUFFLEN, fp );

	int nThisRow = 0;
	while(1)
	{
		if ( NULL == fgets( buff, BUFFLEN, fp ) )
			break;

		char *p = strtok(buff, ",\n"); // get the ID
		p = strtok(NULL, ",\n");       // get X
		p = strtok(NULL, ",\n");       // get Y

		for ( int j=0; j<nCols; j++ ) // get the predictor values
		{
			p = strtok(NULL, ",\n");
			ppDemandData[nThisRow][j] = (float)atof(p);
		}

		++nThisRow;
	}

	// close file
	fclose( fp );
	return(ppDemandData);
}


//
// Get the number of columns in the demand site file 
//
int GetDemandSiteColumnCount(char *PointFile)
{
	// open the site file
	FILE *fp = fopen( PointFile , "r+t" );

	// get the header
	char buff[BUFFLEN];
	fgets( buff, BUFFLEN, fp );

	// determine the number of columns in the header
	int nColumns = 0;

	char *p = strtok(buff, ",\n");
	if ( p ) 
		++nColumns;
	
	do 
	{
		p = strtok(NULL, ",\n");
		if (p) 
			++nColumns;
	} while(p);

	// close file
	fclose( fp );

	return(nColumns);
}


//
// Get the number of records in the demand site file
//
int GetDemandSiteCount(char *PointFile)
{
	// open the site file
	FILE *fp = fopen( PointFile , "r+t" );

	// get the header
	char buff[BUFFLEN];
	fgets( buff, BUFFLEN, fp );

	// determine the number of records in the site file
	int nItems = 0;
	while(1)
	{
		if ( NULL == fgets( buff, BUFFLEN, fp ) )
			break;

		++nItems;
	}

	// close file
	fclose( fp );

	return(nItems);
}



//
// A little test driver to try the ANN class
//
bool DoANNTest(char *pDemandPointPath, int nToSkip)
{
	Message("At entry to DoANNTest", "INFO");
	Message(pDemandPointPath, "pDemandPointPath");
	//Message(nToSkip, "nToSkip");
	ANNTree *myTree = new ANNTree(pDemandPointPath, nToSkip);
	if (myTree) delete myTree;
	Message("At End", "INFO");
	return(true);
}


//
// Classify a GDM stack with a K-Means Centroids table (Version 1.0 - No banding)
//
bool ClassifyGDMViaKMeans(char *pDemandPointPath, int nToSkip, char *pGDMParams, char *OutPath, FPTR fptr)
{
	//Message("At entry to ClassifyGDMViaKMeans", "INFO");

	//
	// setup the ANN tree
	//
	ANNTree *myTree = new ANNTree(pDemandPointPath, nToSkip);
	ANNpoint query_pt = myTree->GetQueryPt();
	int nANNCols = myTree->GetNumCols();


	//
	// setup the grid stack
	//
	BinaryFileStack *bfs = new BinaryFileStack(pGDMParams, BFS_TranGrids);
	if (NULL == bfs)
	{
		Message("Cannot construct BinaryFileStack bfs", "ERROR");
		if (myTree) delete myTree;
		return(false);
	}
	int nRows = bfs->GetNumRows();
	int nCols = bfs->GetNumCols();
	float NoData = bfs->GetNoDataVal();


	//
	// create the output grid and row vector
	//
	gmpath gmPath;
	BinaryFileClass *bfc_Out = new BinaryFileClass(gmPath.ChangeExtension(OutPath, ".flt"), BFC_CreateOrTrunc);	
	float *pOutRow = new float [nCols];


	//
	// do the classification
	//
	int nCurrent = 0;
	fptr("Classify GDM Via K-Means...", nCurrent);	
	for ( int i=0; i<nRows; i++ )
	{
		//
		// update progress
		//
		if ( i * 100 / nRows > nCurrent)
		{
			nCurrent = i * 100 / nRows;
			char banner[256];
			sprintf(banner, "Classify GDM Via K-Means...Row %d of %d", i, nRows);
			fptr(banner, nCurrent);
		}

		// get the grid row data
		bfs->GetNextRow();
		for (int j=0; j<nCols; j++ )
		{
			if (bfs->GetRowData()[0][j] != NoData)
			{
				bfs->GetCellAt(j);
				float *pCellData = bfs->GetCellData();

				for ( int k=0; k<nANNCols; k++ )
				{
					query_pt[k] = pCellData[k];
				}

				// search for the nearest neighbour class
				myTree->Search();

				// assign nearest neighbour class
				pOutRow[j] = (float)myTree->GetIndices()[0] + 1;  // class indices are one based
			}
			else
			{
				pOutRow[j] = NoData;
			}
		}

		//
		// write the output row
		//		
		bfc_Out->WriteFloat(pOutRow, nCols);
	}

	// finish the output grid write
	bfc_Out->Close();
	bfs->CopyHeaderTo(gmPath.ChangeExtension(OutPath, ".hdr"));	

	//
	// create a metadata file for this classification
	//
	FILE *fpMeta = fopen( gmPath.ChangeExtension(OutPath, ".mta"), "w+t");
	fprintf(fpMeta, "##\n## Metadata for Full GDM Classification written to \n## %s\n##\n", OutPath);
	fprintf(fpMeta, "Input ANN Lookup Table: %s\n", pDemandPointPath );
	fprintf(fpMeta, "ANN Records To Search: %d\n", myTree->GetNumRows());
	fprintf(fpMeta, "GDM Predictors: %d\n", myTree->GetNumCols());
	fprintf(fpMeta, "Input GDM Parameter file: %s\n", pGDMParams );
	fprintf(fpMeta, "Output Grid Path: %s\n", OutPath );
	fclose(fpMeta);

	//
	// cleanup
	//
	if (myTree) delete myTree;
	if (bfs) delete bfs;
	if (bfc_Out) delete bfc_Out;
	nCurrent = 0;
	fptr("Status: Ready", nCurrent);
	return(true);
}

//
// Classify a GDM stack with a K-Means Centroids table (Version 2.0 - Use banding)
//
bool ClassifyGDMViaKMeansBanded(char *pDemandPointPath, int nToSkip, 
								char *pGDMParams, 
								char *OutPath, int nLower, int nUpper,
								FPTR fptr)
{
	//Message("At entry to ClassifyGDMViaKMeansBanded", "INFO");
	int nCurrent = 0;
	fptr("Classify GDM Via K-Means...", nCurrent);	

	//
	// setup the ANN tree
	//
	ANNTree *myTree = new ANNTree(pDemandPointPath, nToSkip);
	ANNpoint query_pt = myTree->GetQueryPt();
	int nANNCols = myTree->GetNumCols();
	nCurrent = 20;
	fptr("Classify GDM Via K-Means...", nCurrent);	


	//
	// setup the grid stack
	//
	BinaryFileStack *bfs = new BinaryFileStack(pGDMParams, BFS_TranGrids);
	if (NULL == bfs)
	{
		Message("Cannot construct BinaryFileStack bfs", "ERROR");
		if (myTree) delete myTree;
		return(false);
	}
	int nRows = bfs->GetNumRows();
	int nCols = bfs->GetNumCols();
	float NoData = bfs->GetNoDataVal();
	nCurrent = 50;
	fptr("Classify GDM Via K-Means...", nCurrent);	


	//
	// create the output grid and row vector
	//
	gmpath gmPath;
	BinaryFileClass *bfc_Out = new BinaryFileClass(gmPath.ChangeExtension(OutPath, ".flt"), BFC_CreateOrTrunc);	
	float *pOutRow = new float [nCols];

	//
	// create a noData Row for the out-of-band outputs
	//
	float *pNoDataRow = new float [nCols];
	for ( int i=0; i<nCols; i++ ) pNoDataRow[i] = NoData;
	nCurrent = 80;
	fptr("Classify GDM Via K-Means...", nCurrent);	

	//
	// do the classification
	//
	int nRange = nUpper - nLower;	
	nCurrent = 0;
	fptr("Classify GDM Via K-Means...", nCurrent);	
	for ( int i=0; i<nRows; i++ )
	{
		if ((i >= nLower) && (i <= nUpper))
		{
			//
			// update progress
			//
			if ( i * 100 / nRange > nCurrent)
			{
				nCurrent = (i-nLower) * 100 / nRange;
				char banner[256];
				sprintf(banner, "Classify GDM Via K-Means...Banded Row %d of %d", i, nUpper);
				fptr(banner, nCurrent);
			}

			// get the grid row data
			bfs->GetRowAt(i);


			for (int j=0; j<nCols; j++ )
			{
				if (bfs->GetRowData()[0][j] != NoData)
				{
					bfs->GetCellAt(j);
					float *pCellData = bfs->GetCellData();

					for ( int k=0; k<nANNCols; k++ )
					{
						query_pt[k] = pCellData[k];
					}

					// search for the nearest neighbour class
					myTree->Search();

					// assign nearest neighbour class
					pOutRow[j] = (float)myTree->GetIndices()[0] + 1;  // class indices are one based
				}
				else
				{
					pOutRow[j] = NoData;
				}
			}

			//
			// write the output row
			//
			bfc_Out->WriteFloat(pOutRow, nCols);
		}
		else
		{
			//
			// write the output row
			//
			bfc_Out->WriteFloat(pNoDataRow, nCols);
		}
	}

	// finish the output grid write
	bfc_Out->Close();
	bfs->CopyHeaderTo(gmPath.ChangeExtension(OutPath, ".hdr"));	

	//
	// create a metadata file for this classification
	//
	FILE *fpMeta = fopen( gmPath.ChangeExtension(OutPath, ".mta"), "w+t");
	fprintf(fpMeta, "##\n## Metadata for Banded GDM Classification written to \n## %s\n##\n", OutPath);
	fprintf(fpMeta, "Input ANN Lookup Table: %s\n", pDemandPointPath );
	fprintf(fpMeta, "ANN Records To Search: %d\n", myTree->GetNumRows());
	fprintf(fpMeta, "GDM Predictors: %d\n", myTree->GetNumCols());
	fprintf(fpMeta, "Input GDM Parameter file: %s\n", pGDMParams );
	fprintf(fpMeta, "Output Grid Path: %s\n", OutPath );
	fprintf(fpMeta, "Lower Band Index: %d\n", nLower);
	fprintf(fpMeta, "Upper Band Index: %d\n", nUpper);
	fclose(fpMeta);

	//
	// cleanup
	//
	if (pOutRow) delete[] pOutRow;
	if (pNoDataRow) delete[] pNoDataRow;
	if (myTree) delete myTree;
	if (bfs) delete bfs;
	if (bfc_Out) delete bfc_Out;
	nCurrent = 0;
	fptr("Status: Ready", nCurrent);
	return(true);
}




////////////////////////////////////// Test processing time for version of TS Impact code ///////////////////////////////////////
//
//
// Implimentation of Dan Faith's Climate Change Impact Approach 
// based on the sum of minimum distances between demand points and ALL valid gridcells
//
bool DoTimeSeriesImpactTable_Version_1(char *pDemandPointPath, char *pStartParams, FPTR fptr)
{
	Message("At entry to DoTimeSeriesImpactTable_Version_1", "INFO");

	//
	// Fill in the demand point data frame
	//
	int nDemandSites = GetDemandSiteCount(pDemandPointPath);
	//Message(nDemandSites, "info");

	float **ppDemandData = CreateDemandSiteDataFromCSVTableV5(pDemandPointPath);
	if (NULL == ppDemandData)
	{
		Message("Cannot construct ppDemandData", "ERROR");
		return(false);
	}


	//
	// setup the grid stacks
	//
	BinaryFileStack *bfs0 = new BinaryFileStack(pStartParams, BFS_TranGrids);
	if (NULL == bfs0)
	{
		Message("Cannot construct BinaryFileStack bfs0", "ERROR");
		return(false);
	}


	//
	// create the demand class min-value vectors
	//	
	float *pMinRow0 = new float [nDemandSites];
	for ( int i=0; i<nDemandSites; i++ )
	{
		pMinRow0[i] = 100.0F;
	}


	//
	// Run the analysis
	//
	int nCurrent = 0;
	fptr("Creating CC Impact Table Version 1...", nCurrent);
	for ( int i=0; i<bfs0->GetNumRows(); i++ )
	{
		//
		// update progress
		//
		if ( i * 100 / bfs0->GetNumRows() > nCurrent)
		{
			nCurrent = i * 100 / bfs0->GetNumRows();
			char banner[256];
			sprintf(banner, "Creating CC Impact Table Version 1...Row %d of %d", i, bfs0->GetNumRows());
			fptr(banner, nCurrent);
		}

		// get the dissimilarities
		bfs0->GetNextRow();
		for (int j=0; j<bfs0->GetNumCols(); j++ )
		{
			if (bfs0->GetRowData()[0][j] != bfs0->GetNoDataVal())
			{
				bfs0->GetCellAt(j);

				/////////////////////////////////////////////////////////////////////////////////////
				//
				// This is the MINIMUM dissimilarity version (as per Dan Faith's original calc)
				//
				//
				// work thru the demand sites and get the minimum dissimilarity for each epoch
				for ( int k=0; k<nDemandSites; k++ )
				{
					float f0 = (float)bfs0->GetGDMDissimilarity(ppDemandData[k]);
					if ( f0 < pMinRow0[k] ) pMinRow0[k] = f0;
				}
			}
		}
	}
	
	//
	// Clean up
	//
	if (pMinRow0) delete[] pMinRow0;
	for (int i=0; i<nDemandSites; i++ )
	{
		if (ppDemandData[i]) delete[] ppDemandData[i];
	}
	if (ppDemandData) delete[] ppDemandData;
	if (bfs0) delete bfs0;
	nCurrent = 0;
	fptr("Status: Ready", nCurrent);
	return(true);
}


//
// Implimentation of Dan Faith's Climate Change Impact Approach 
// based on the sum of minimum distances between demand points and ALL valid gridcells
//
bool DoTimeSeriesImpactTable_Version_2(char *pDemandPointPath, char *pStartParams, FPTR fptr)
{
	Message("At entry to DoTimeSeriesImpactTable_Version_2", "INFO");

	//
	// Fill in the demand point data frame
	//
	int nDemandSites = GetDemandSiteCount(pDemandPointPath);
	//Message(nDemandSites, "info");

	float **ppDemandData = CreateDemandSiteDataFromCSVTableV5(pDemandPointPath);
	if (NULL == ppDemandData)
	{
		Message("Cannot construct ppDemandData", "ERROR");
		return(false);
	}


	//
	// setup the grid stacks
	//
	BinaryFileStack *bfs0 = new BinaryFileStack(pStartParams, BFS_TranGrids);
	if (NULL == bfs0)
	{
		Message("Cannot construct BinaryFileStack bfs0", "ERROR");
		return(false);
	}


	//
	// create the demand class min-value vectors
	//	
	float *pMinRow0 = new float [nDemandSites];
	for ( int i=0; i<nDemandSites; i++ )
	{
		pMinRow0[i] = 100.0F;
	}


	//
	// Run the analysis
	//
	int nCurrent = 0;
	fptr("Creating CC Impact Table Version 2...", nCurrent);
	for ( int i=0; i<bfs0->GetNumRows(); i++ )
	{
		//
		// update progress
		//
		if ( i * 100 / bfs0->GetNumRows() > nCurrent)
		{
			nCurrent = i * 100 / bfs0->GetNumRows();
			char banner[256];
			sprintf(banner, "Creating CC Impact Table Version 2...Row %d of %d", i, bfs0->GetNumRows());
			fptr(banner, nCurrent);
		}

		// get the dissimilarities
		bfs0->GetNextRow();
		for (int j=0; j<bfs0->GetNumCols(); j++ )
		{
			if (bfs0->GetRowData()[0][j] != bfs0->GetNoDataVal())
			{
				bfs0->GetCellAt(j);

				/////////////////////////////////////////////////////////////////////////////////////
				//
				// This is the MINIMUM dissimilarity version (as per Dan Faith's original calc)
				//
				//
				// work thru the demand sites and get the minimum dissimilarity for each epoch
				for ( int k=0; k<nDemandSites; k++ )
				{
					float f0 = (float)bfs0->GetManhattanDist(ppDemandData[k]);
					if ( f0 < pMinRow0[k] ) pMinRow0[k] = f0;
				}
			}
		}
	}
	
	//
	// Clean up
	//
	if (pMinRow0) delete[] pMinRow0;
	for (int i=0; i<nDemandSites; i++ )
	{
		if (ppDemandData[i]) delete[] ppDemandData[i];
	}
	if (ppDemandData) delete[] ppDemandData;
	if (bfs0) delete bfs0;
	nCurrent = 0;
	fptr("Status: Ready", nCurrent);
	return(true);
}


