//
// GetDemandPoints.cpp
//

#include "stdafx.h"
#include "GetDemandPoints.h"
#include "Message.h"
#include "GDMBufferSizeDefs.h"
#include "GDMClassification.h"
#include "ParamsW16.h"
#include "clsDoPath.h"
#include "ann.h" 
#include "clsBinaryFileClass.h"
#include "clsEsriBinaryHeader.h"
#include "ranlib.h"
#include "cluster.h"

double CalcED(double *data_0, double *data_1, int nCols);
bool AlreadySelected(int *pExtractedIndices, int NumDemandPointsExtracted, int nn_idx);


bool GetDemandPoints_01( char *pParams, char *OutPath, int NumSamples, int NumDemandPoints, FPTR fptr)
{
	//Message(pParams, "pParams");
	//Message(OutPath, "OutPath");
	//Message(NumDemandPoints, "NumDemandPoints");
	//Message(NumSamples, "NumSamples");
	//
	// Get a path to the domain (use the first transform grid)
	//
	char DomainPath[FILEPATH_BUFFSIZE];
	GetProfileString( "TRANSPREDS", "PredTran1", DomainPath, pParams );
	//Message(DomainPath, "DomainPath");

	
	//
	// Extract a set of samples points for the nearest neighbour classification and write to
	// a comma delimited text file in the Working directory called ClassificationSample.csv
	//
	gmpath GmPath;
	char lpSamplePath[FILEPATH_BUFFSIZE];
	sprintf(lpSamplePath, "%s", GmPath.ChangeExtension(OutPath, ".txt"));
	if (!CreateSamplePointMesh( pParams, lpSamplePath, DomainPath, NumSamples, false, fptr ))
	{
		Message("Cannot CreateSamplePointFile()", "GetDemandPoints_01");
		return(false);
	}
	//Message(lpSamplePath, "lpSamplePath");


	//
	// Extract transform grid data at the sample point locations
	// into a matrix to pass into the classification routines
	//
	int nRows, nCols;
	double **ppData = ExtractDataFromTransformGrids( pParams, DomainPath, lpSamplePath, &nRows, &nCols, false );
	if ( NULL == ppData )
	{
		Message( "ExtractDataFromTransformGrids failure", "GetDemandPoints_01" );
		return( false );
	}


	//
	// Write output
	//
	//FILE *fp = fopen(OutPath, "w+t");
	////Message(OutPath, "OutPath");
	//for ( int i=0; i<nRows; i++ )
	//{
	//	for ( int j=0; j<nCols; j++ )
	//	{
	//		fprintf(fp, "%lf", ppData[i][j]);

	//		if (j < nCols-1)
	//			fprintf(fp, ",");
	//		else
	//			fprintf(fp, "\n");
	//	}
	//}
	//fclose(fp);


	//
	// Get Transformed sample points into an ANN tree
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
	ANNidxArray nn_idx = new ANNidx[nRows];
	ANNdistArray dist = new ANNdist[nRows];
	ANNpoint query_pt = annAllocPt(nCols);


	//
	// build the kd tree search structure with the data points, number of data points and dimension of search space
	//
	ANNkd_tree *theTree = new ANNkd_tree ( data_pnts, nRows, nCols ); 


	//
	// Initialise the random number generator
	//
	initran();
	int nRand = int(genunf(0.0F, (float)nRows));
	//Message(nRand);
	

	//
	// Initialise for the demand point extraction by using the midpoint record
	//
	int *pExtractedIndices = new int [NumDemandPoints];
	double *pExtractedDistances = new double [NumDemandPoints];
	//pExtractedIndices[0] = nRows / 2;
	pExtractedIndices[0] = nRand;
	int NumDemandPointsExtracted = 1;


	//
	// Now find the rest...
	//
	int nCurrent = 0;
	fptr("Extracting Demand Points...", nCurrent);	
	while (NumDemandPointsExtracted < NumDemandPoints)
	{
		if ( NumDemandPointsExtracted * 100 / NumDemandPoints > nCurrent )
		{
			nCurrent = NumDemandPointsExtracted * 100 / NumDemandPoints;
			fptr("Extracting Demand Points...", nCurrent);	
		}

		// initialise from the current demand point
		for ( int i=0; i<nCols; i++ )
		{
			query_pt[i] = data_pnts[NumDemandPointsExtracted-1][i];
		}
		
		// sort the tree by manhattan distance
		theTree->annkSearch(query_pt, nRows, nn_idx, dist, 0);		

		//
		// locate the cutpoint in the results that partitions the 
		// sorted data into nearer/further ED from current demand point
		//
		//double CurrentCutpoint = 1.5;   // EDs above this will be 80% or more dissimilar
		//double CurrentCutpoint = 2.7;   // EDs above this will be 90% or more dissimilar
		double CurrentCutpoint = 3.0;     // EDs above this will be 95% or more dissimilar
		for ( int j=nRows-1; j>0; j-- )
		{
			//
			// This will occur if the nearest neighbours get closer than the cutpoint
			//
			if (dist[j-1] < CurrentCutpoint)
			{
				//
				// This will occur if the nearest neighbours get closer than the cutpoint
				//
				// choose a random item from this selection
				nRand = int(genunf(float(j-1), (float)nRows));

				pExtractedIndices[NumDemandPointsExtracted] = nn_idx[nRand];
				pExtractedDistances[NumDemandPointsExtracted] = dist[nRand];
				++NumDemandPointsExtracted;
				break;
			} // if (dist[j-1] < CurrentCutpoint)
		} // for ( int j=nRows-1; j>0; j-- )

	} // while (NumDemandPointsExtracted < NumDemandPoints)

	
	//
	// Extract the X and Y coords
	//
	char SiteBuff[FILEPATH_BUFFSIZE];
	double *dX = new double [nRows];
	double *dY = new double [nRows];
	FILE *fpSites = fopen(lpSamplePath, "r+t");
	fgets(SiteBuff, FILEPATH_BUFFSIZE, fpSites); // get header
	for ( int i=0; i<nRows; i++ )
	{
		fgets(SiteBuff, FILEPATH_BUFFSIZE, fpSites);

		char *p = strtok( SiteBuff, ",\n"); // get site ID

		p = strtok( NULL, ",\n");           // get X
		dX[i] = atof(p);

		p = strtok( NULL, ",\n");           // get Y
		dY[i] = atof(p);
	}
	fclose(fpSites);


	//
	// Dump out the Demand Point Records
	//
	FILE *fpOut = fopen(OutPath, "w+t");
	fprintf(fpOut, "_ID,Index,Dist,X,Y,");
	for (int i=0; i<nCols; i++ )
	{
		fprintf(fpOut, "Pred_%d", i+1);
		if (i<nCols-1)
			fprintf(fpOut, ",");
		else
			fprintf(fpOut, "\n");
	}
	for ( int i=0; i<NumDemandPoints; i++ )
	{
		int thisIndex = pExtractedIndices[i];
		fprintf(fpOut, "%d,%d,%lf,%lf,%lf,", i+1, thisIndex, pExtractedDistances[i], dX[thisIndex], dY[thisIndex]);
		for (int j=0; j<nCols; j++ )
		{
			fprintf(fpOut, "%lf", data_pnts[pExtractedIndices[i]][j]);
			if (j<nCols-1)
				fprintf(fpOut, ",");
			else
				fprintf(fpOut, "\n");
		}
	}
	fclose(fpOut);


	//
	// dump out a distance table
	//
	gmpath gmPath;
	char pDistPath[BUFFLEN];
	sprintf(pDistPath, "%s\\%s_dist.csv", gmPath.GetDirectoryPath(OutPath), gmPath.GetName(OutPath));
	//Message(pDistPath, "Distance Matrix");
	FILE *fpDist = fopen(pDistPath, "w+t");
	fprintf(fpDist, "X,Y,Average_Dist,Count < 5%% Sim.,");
	for ( int i=0; i<NumDemandPoints; i++ )
	{
		fprintf(fpDist, "P%d", i+1);
		if (i<NumDemandPoints-1)
			fprintf(fpOut, ",");
		else
			fprintf(fpOut, "\n");
	}
	for ( int i=0; i<NumDemandPoints; i++ )
	{
		double dAverage = 0.0;
		int nCount = 0;
		for ( int j=0; j<NumDemandPoints; j++ )
		{
			if (i != j)
			{
				double dED = CalcED(data_pnts[pExtractedIndices[i]], data_pnts[pExtractedIndices[j]], nCols);
				dAverage += dED;

				if (exp(-dED) <= 0.05) ++nCount;
			}
		}
		dAverage /= 100.0;
		dAverage = exp(-dAverage);

		fprintf(fpDist, "%lf,%lf,%lf,%d,", dX[pExtractedIndices[i]], dY[pExtractedIndices[i]], dAverage, nCount);
		for ( int j=0; j<NumDemandPoints; j++ )
		{
			double dED = CalcED(data_pnts[pExtractedIndices[i]], data_pnts[pExtractedIndices[j]], nCols);
			fprintf(fpOut, "%lf", exp(-dED));
						
			if (j<NumDemandPoints-1)
				fprintf(fpOut, ",");
			else
				fprintf(fpOut, "\n");
		}
	}
	fclose(fpDist);

	
	//
	// Cleanup
	//
	if (dX) delete[] dX;
	if (dY) delete[] dY;
	if (pExtractedDistances) delete[] pExtractedDistances;
	if (pExtractedIndices) delete[] pExtractedIndices;
	if (theTree) delete theTree;
	for ( int i=0; i<nRows; i++ ) if (ppData[i]) delete[] ppData[i];
	if (ppData) delete[] ppData;
	return(true);
}


double CalcED(double *data_0, double *data_1, int nCols)
{
	double dCalc = 0.0;
	for (int i=0; i<nCols; i++ )
	{
		dCalc += fabs(data_0[i] - data_1[i]);
	}
	return(dCalc);
}


bool AlreadySelected(int *pExtractedIndices, int NumDemandPointsExtracted, int nn_idx)
{
	for ( int i=0; i<NumDemandPointsExtracted; i++ )
	{
		if (pExtractedIndices[i] == nn_idx)
			return(true);
	}
	return(false);
}





//#region Old_Code


	////
	//// Extract the filepaths for the transform grids in the parameter file
	////
	//bool fHaveEuclidean = GetProfileInt( "GDMODEL", "UseEuclidean", pParams ) == 1 ? true : false;
	//int nGrids = 0;
	//if ( fHaveEuclidean ) nGrids += 2;

	//int nEnvGrids = GetProfileInt( "PREDICTORS", "NumPredictors", pParams );
	//char lpKey[64];
	//for ( int i=1; i<=nEnvGrids; i++ )
	//{
	//	//
	//	// we ONLY transform grids that have NON-ZERO coefficients
	//	//
	//	char lpPredPath[256];
	//	sprintf( lpKey, "PredTran%d", i );
	//	GetProfileString( "TRANSPREDS", lpKey, lpPredPath, pParams );
	//	if ( strlen( lpPredPath ) > 0 )      //if ( fPredictorHasNonZeroCoeffs( lpParams, i ) )
	//	{
	//		++nGrids;
	//	}
	//}
	//
	//
	////
	//// create an array for the file handles to the transform grids
	////
	//gmpath gmPath;
	//char lpKey[64];
	//char pGridPath[BUFFLEN];
	//BinaryFileClass **bfcTranGrids = new BinaryFileClass * [nGrids];
	//for ( int i=0; i<nGrids; i++ )
	//{
	//	if (fHaveEuclidean)
	//	{
	//		if (i==0)
	//		{
	//			GetProfileString( "TRANSPREDS", "EuclXTran", pGridPath, pParams );
	//			bfcTranGrids[i] = new BinaryFileClass(gmPath.ChangeExtension(pGridPath, ".flt"));
	//		}

	//		else if (i==1)
	//		{
	//			GetProfileString( "TRANSPREDS", "EuclYTran", pGridPath, pParams );
	//			bfcTranGrids[i] = new BinaryFileClass(gmPath.ChangeExtension(pGridPath, ".flt"));
	//		}

	//		else
	//		{
	//			sprintf( lpKey, "PredTran%d", i+1 );
	//			GetProfileString( "TRANSPREDS", lpKey, pGridPath, pParams );
	//			if ( strlen( pGridPath ) > 0 )
	//				bfcTranGrids[i] = new BinaryFileClass(gmPath.ChangeExtension(pGridPath, ".flt"));
	//		}
	//	}

	//	else
	//	{
	//		sprintf( lpKey, "PredTran%d", i+1 );
	//		GetProfileString( "TRANSPREDS", lpKey, pGridPath, pParams );
	//		if ( strlen( pGridPath ) > 0 )
	//			bfcTranGrids[i] = new BinaryFileClass(gmPath.ChangeExtension(pGridPath, ".flt"));
	//	}
	//}

	//
	////
	//// Open domain grid
	////
	//char lpDomainPath[BUFFLEN];
	//strcpy(lpDomainPath,DomainPath);
	//sprintf(lpDomainPath, "%s", gmPath.ChangeExtension(lpDomainPath, ".flt"));
	//BinaryFileClass *bfcDomain = new BinaryFileClass(lpDomainPath);

	////
	//// get some grid metrics from the domain grid
	////	
	//EsriBinaryHeader *DomainHeader = new EsriBinaryHeader(gmPath.ChangeExtension(lpDomainPath, ".hdr"));
	//int nDomainGridRows = DomainHeader->GetNumRows();
	//int nDomainGridCols = DomainHeader->GetNumCols();
	//float fDomainNoData = DomainHeader->GetNoDataValue();
	//double dDomainResX = DomainHeader->GetCellSize();
	//double dDomainMinX = DomainHeader->GetXllCorner();
	//double dDomainMinY = DomainHeader->GetYllCorner();
	//double dDomainMaxX = dDomainMinX + (dDomainResX * nDomainGridCols);
	//double dDomainMaxY = dDomainMinY + (dDomainResX * nDomainGridRows);

	//float **ppInputRows = new float * [nGrids];
	//for ( int i=0; i<nGrids; i++ ) ppInputRows[i] = new float [nDomainGridCols];
	//int *pDomainRow = new int [nDomainGridCols];


	////
	//// Extract the demand points
	////



	////
	//// Cleanup
	////
	//if (pDomainRow) delete[] pDomainRow;
	//for (int i=0; i<nGrids; i++)
	//{
	//	if (ppInputRows[i]) delete[] ppInputRows[i];
	//	bfcTranGrids[i]->Close();
	//	if (bfcTranGrids[i]) delete[] bfcTranGrids[i];
	//}
	//if (bfcTranGrids) delete[] bfcTranGrids;
	//if (ppInputRows) delete ppInputRows;


#undef DO_THIS
#ifdef DO_THIS

		for ( int j=nRows-1; j>0; j-- )
		{
			//
			// This will occur if the nearest neighbours get closer than the cutpoint
			//
			if (dist[j-1] < CurrentCutpoint)
			{
				//
				// locate the record that maximises the average ED between the already 
				// selected demand points and those that remain unselected beyond the cutpoint
				//
				int nBestIndex = -1;
				double dBestDist = 0.0;
								
				for ( int k=j; k<nRows; k++ )
				{
					// make sure that we don't reselect the same sites
					if (AlreadySelected(pExtractedIndices, NumDemandPointsExtracted, nn_idx[k]))
						continue;

					double dTmpDist = 0.0;
					
					for (int nCrntExt=0; nCrntExt<NumDemandPointsExtracted; nCrntExt++ )
					{
						for ( int q=0; q<nCols; q++ )
						{
							dTmpDist += fabs(data_pnts[pExtractedIndices[nCrntExt]][q] - data_pnts[nn_idx[k]][q]);
						}											
					}
					
					if (dTmpDist > dBestDist)
					{
						dBestDist = dTmpDist;
						nBestIndex = k;
					}
				}

				if (nBestIndex == -1)
				{
					Message("Cannot locate new demand point", "ERROR");
				}
				else
				{
					pExtractedIndices[NumDemandPointsExtracted] = nn_idx[nBestIndex];
					pExtractedDistances[NumDemandPointsExtracted] = dist[nBestIndex];
					++NumDemandPointsExtracted;
				}

				break;

			} // if (dist[j-1] < CurrentCutpoint)
		} // for ( int j=nRows-1; j>0; j-- )
		//
		// locate the furthest unselected site to the current selected demand point
		//
		for ( int j=nRows-1; j>0; j-- )
		{
			//
			// This will occur if the nearest neighbours get closer than the cutpoint
			//
			if (dist[j-1] < CurrentCutpoint)
			{
				// choose a random item from this selection
				nRand = int(genunf(float(j-1), (float)nRows));

				/*pExtractedIndices[NumDemandPointsExtracted] = nn_idx[j];
				pExtractedDistances[NumDemandPointsExtracted] = dist[j];*/
				
				pExtractedIndices[NumDemandPointsExtracted] = nn_idx[nRand];
				pExtractedDistances[NumDemandPointsExtracted] = dist[nRand];
				++NumDemandPointsExtracted;
				break;
			}
		}

		for ( int j=nRows-1; j>0; j-- )
		{
			int nCurrent = nn_idx[j];
			bool found = false;
			for ( int k=0; k<NumDemandPointsExtracted; k++ )
			{
				if (pExtractedIndices[k] == nCurrent)
				{
					found = true;
					break;
				}
			}

			if (!found) // use it
			{
				pExtractedIndices[NumDemandPointsExtracted] = nn_idx[j];
				pExtractedDistances[NumDemandPointsExtracted] = dist[j];
				++NumDemandPointsExtracted;
				break;
			}
		}
#endif // DO_THIS

//#endregion