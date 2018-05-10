//
// GDMDensity.cpp
//
#include "stdafx.h"
#include "GdmBinLib.h"
#include "GDMDensity.h"
#include "GDMBufferSizeDefs.h"
#include "clsDoPath.h"
#include "ParamsW16.h"
#include "Message.h"
#include "clsEsriBinaryHeader.h"
#include "clsBinaryFileClass.h"
#include "GDMClassification.h"

////////////////////////////////////////////////////////////////
//
// Using an ANN tree is good if wwe are applying cutpoints
// if not, then iterate thru all the samples points
//
//#define USE_ANN_TREE
#undef USE_ANN_TREE
////////////////////////////////////////////////////////////////

//
// Create GDM Density Grid (called externally)
//
bool DoGDMDensitySim(char *Params, 
					 char *DomainPath,
					 char *OutName,
					 int NumSamples,
 				     FPTR fptr)
{
	int nValidCells = GetValidDataCellCount(DomainPath, fptr);
	if ( nValidCells < 1 )
	{
		Message( "No Valid Non-Zero cells to use", "DoGDMDensitySim" );
		return(false);
	}
	double SampleMultiple = (double)nValidCells / (double)NumSamples;
	//Message(SampleMultiple, "SampleMultiple");

	char myBuff[BUFFLEN];
	int NumPredGrids = GetPrivateProfileInt("PREDICTORS", "NumPredictors", 0, Params);
	if ( NumPredGrids < 1 )
	{
		Message( "No Predictor Grids to use???", "DoGDMDensitySim" );
		return(false);
	}

	int numTranGrids = 0;
	char testBuff[64];
	for ( int i=1; i<=NumPredGrids; i++ )
	{
		sprintf( testBuff, "PredTran%d", i );
		GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, Params);
		if ( 0 != strcmp("0",myBuff) )
		{
			++numTranGrids;
		}
	}

	if ( numTranGrids < 1 )
	{
		Message( "No Transform Grids to use???", "DoGDMDensitySim" );
		return(false);
	}


	//
	// allocate strings for transformed input files
	//
	char **ppTranGridPaths = new char * [numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		ppTranGridPaths[i] = new char [BUFFLEN];
	}


	//
	// Collect the Tran grid paths
	//
	int nThis = 0;
	for ( int i=1; i<=NumPredGrids; i++ )
	{
		sprintf( testBuff, "PredTran%d", i );
		GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, Params);
		if ( 0 != strcmp("0",myBuff) )
		{
			strcpy (ppTranGridPaths[nThis++], myBuff );
		}
	}


	//
	// Generate the sample data points
	//
	char pWorkspace[BUFFLEN];
	GetProfileString( "GDMODEL", "WorkspacePath", pWorkspace, Params );
	char pSamplePath[BUFFLEN];
	sprintf( pSamplePath, "%s\\%s\\%s.csv", pWorkspace, OutName, OutName );
	bool GotSamplePoints = CreateSamplePointMeshForDensity( Params, pSamplePath, DomainPath, NumSamples, false, fptr );
	if (!GotSamplePoints)
	{
		Message("Cannot CreateSamplePointFile", "DoGDMDensitySim");
		return(false);
	}


	// 
	// determine the number of rows ( the number of records in the points file )
	//
	int nRecords = 0;
	char buff[BUFFLEN];
	FILE *fpPoints = fopen( pSamplePath, "r+t" );
	// get the header
	fgets( buff, BUFFLEN, fpPoints );
	
	while( 1 )
	{
		if ( NULL == fgets( buff, BUFFLEN, fpPoints ) )
			break;

		++nRecords;
	}

	double **ppPoints = new double * [nRecords];
	for ( int i=0; i<nRecords; i++ )
	{
		ppPoints[i] = new double [2];
	}

	rewind(fpPoints);
	// get the header
	fgets( buff, BUFFLEN, fpPoints );
	char seps[] = ",\n";
	int nCurrent = 0;
	fptr("Creating GDM Density Grid (A)", nCurrent);	
	for ( int i=0; i<nRecords; i++ )
	{
		if ( i * 100 / nRecords > nCurrent )
		{
			nCurrent = i * 100 / nRecords;
			fptr("Creating GDM Density Grid (A)", nCurrent);	
		}
		
		if ( NULL == fgets( buff, BUFFLEN, fpPoints ) )
			break;

		char *p = strtok( buff, seps ); // get ID
		
		p = strtok( NULL, seps );       // get easting
		ppPoints[i][0] = atof(p);

		p = strtok( NULL, seps );       // get northing
		ppPoints[i][1] = atof(p);
	}
	fclose(fpPoints);	


	double **data_pnts = new double * [nRecords];
	for ( int i=0; i<nRecords; i++ )
	{
		data_pnts[i] = new double [numTranGrids];
	}


	// now build the kd tree search structure
	ANNkd_tree *theTree = new ANNkd_tree ( data_pnts,		// the data points
										   nRecords,		// number of data points
										   numTranGrids );	// dimension of search space

	ANNpoint query_pt = annAllocPt( numTranGrids );
	ANNidxArray nn_idx = new ANNidx[nRecords];
	ANNdistArray dist = new ANNdist[nRecords];


	nCurrent = 0;
	fptr("Creating GDM Density Grid (B)", nCurrent);	
	for ( int i=0; i<numTranGrids; i++ )
	{
		if ( i * 100 / numTranGrids > nCurrent )
		{
			nCurrent = i * 100 / numTranGrids;
			fptr("Creating GDM Density Grid (B)", nCurrent);	
		}

		DoESRIBinaryDataColumn( data_pnts, i, ppTranGridPaths[i], ppPoints, nRecords );
	}


	//
	// setup a vector for the condition of each sample point in the condition grid
	//
	float *fCondition = new float [nRecords];
	GetESRIBinaryConditionVals( DomainPath, fCondition, ppPoints, nRecords );


	//
	// write the contents of the sample training data table
	//
	nCurrent = 0;
	fptr("Writing Sample Data Table...", nCurrent);	
	FILE *fpData = fopen( pSamplePath, "w+t" );
	fprintf(fpData, "Condition,X,Y,");
	for (int i=1; i<=numTranGrids; i++ )
	{
		if ( i < numTranGrids)
			fprintf(fpData, "Pred%d,", i);
		else
			fprintf(fpData, "Pred%d\n", i);
	}
	for ( int i=0; i<nRecords; i++)
	{
		if ( i * 100 / nRecords > nCurrent )
		{
			nCurrent = i * 100 / nRecords;
			fptr("Writing Sample Data Table...", nCurrent);	
		}

		fprintf(fpData, "%f,%lf,%lf,", fCondition[i], ppPoints[i][0], ppPoints[i][1]);
		for ( int j=0; j<numTranGrids; j++ )
		{
			fprintf(fpData, "%lf", data_pnts[i][j]);
			if ( j<numTranGrids-1 )
				fprintf(fpData, ",");
			else
				fprintf(fpData, "\n");
		}
	}
	fclose(fpData);


	//
	// open the header file (.hdr)
	//
	gmpath gmPath;
	char tmpBuf[500];
	strcpy(tmpBuf, gmPath.ChangeExtension(DomainPath, ".hdr"));
	EsriBinaryHeader *header = new EsriBinaryHeader(tmpBuf);
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	float fNoData = header->GetNoDataValue();
	double dMinX = header->GetMinX();
	double dMaxX = header->GetMaxX();
	double dMinY = header->GetMinY();
	double dMaxY = header->GetMaxY();
	double dCellSize = header->GetCellSize();

	BinaryFileClass *bfcDomain = new BinaryFileClass(gmPath.ChangeExtension(DomainPath, ".flt"));
	if ( !bfcDomain->IsValid() )
	{
		Message("Cannot open Domain Grid", "DoGDMDensitySim");
		return(false);
	}
	
	//
	// create output files
	//
	char pOutPath[BUFFLEN];
	sprintf( pOutPath, "%s\\%s\\%s.flt", pWorkspace, OutName, OutName );

	header->CopyTo(gmPath.ChangeExtension(pOutPath, ".hdr"));
	delete header;
	
	BinaryFileClass *bfcOut = new BinaryFileClass(gmPath.ChangeExtension(pOutPath, ".flt"), BFC_CreateOrTrunc);
	if ( !bfcOut->IsValid() )
	{
		Message("Cannot create OutGrid", "DoGDMDensitySim");
		return(false);
	}

	float *domainRow = new float [nCols];
	float *pOutRow = new float [nCols];

	//
	// setup for the transform grids
	//
	float **ppRows = new float * [numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		ppRows[i] = new float [nCols];
	}
	
	BinaryFileClass **bfcTranGrids = new BinaryFileClass * [numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		bfcTranGrids[i] = new BinaryFileClass(gmPath.ChangeExtension(ppTranGridPaths[i], ".flt"));
		if ( !bfcTranGrids[i]->IsValid() )
		{
			Message("Cannot open Transform Grid", "DoGDMDensitySim");
			return(false);
		}
	}


	//
	//  The domain can be used in a number of ways....
	//
	//  Only cells in the domain that have data are used for the density calculations.
	//	Domain values of 0.0 will have their density calculated but will NOT have any sample points placed in them.
	//  Domain values > 0.0 will have sample points placed in them.
	//  Domain values > 0.0 and <= 1.0 will also be used to weight the density calculations. 
	//  Domain values < 0.0 and > 1.0 will be ignored and set to No-Data
	//
	nCurrent = 0;
	fptr("Creating GDM Density Grid", nCurrent);	
	for ( int i=0; i<nRows; i++ )
	{
		if ( i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			sprintf(myBuff, "Creating GDM Density Grid <%s> Row %d of %d", OutName, i, nRows );
			fptr(myBuff, nCurrent);	
		}

		//
		// Initialise the output row to no-Data
		//
		for ( int j=0; j< nCols; j++ )
		{
			pOutRow[j] = fNoData;
		}

		//
		// Get the domain data row
		//
		bfcDomain->ReadFloat(domainRow, nCols);

		//
		// Get the transform grid data rows
		//
		for ( int j=0; j<numTranGrids; j++ )
		{
			bfcTranGrids[j]->ReadFloat(ppRows[j], nCols);
		}	

		//
		// Do the density calculations for the current row
		//		
		for ( int j=0; j<nCols; j++ )
		{
			if (( domainRow[j] >= 0.0F ) && ( domainRow[j] <= 1.0F ))
			{
				//
				// Do the nearest neighbour search...
				//
#ifdef USE_ANN_TREE
				for ( int k=0; k<numTranGrids; k++ )
				{
					query_pt[k] = ppRows[k][j];
				}

				theTree->annkSearch( query_pt, nRecords, nn_idx, dist, 0 );
#endif // USE_ANN_TREE

				// sum the EDs
				double dVal = 0.0;
				for ( int r=0; r<nRecords; r++ )
				{
#ifdef USE_ANN_TREE
					double dSim = exp(-dist[r]);
					double dCutpoint = 0.5;
					if ( dSim < dCutpoint ) 
					{		
						break;
					}

#else
					double fDist = CalculateED(ppRows, j, data_pnts[r], numTranGrids);
					double dSim = exp(-fDist);
#endif // USE_ANN_TREE

					// now sum after squaring the similarity
					dVal += dSim * dSim;
				}
					
				//
				// set density value after scaling sum of densities 
				// by the proportion of sample sites to all grid sites
				//
				//pOutRow[j] = (float)(dVal / dCount * SampleMultiple);

				//
				// set density value applying the sum of the 
				// similarities between this cell and the sample sites
				//
				pOutRow[j] = (float)(dVal); 
				
				//
				// set density value applying the sum of the 
				// similarities between this cell and the sample sites weighted by the domain grid
				//
				//pOutRow[j] = (float)(dVal) * domainRow[j]; 

			}  // if (( domainRow[j] >= 0.0F ) && ( domainRow[j] <= 1.0F ))
		} // for ( j=0; j<nCols; j++ )
		bfcOut->WriteFloat(pOutRow, nCols);
	}

	bfcDomain->Close();
	for ( int i=0; i<numTranGrids; i++ )
	{
		bfcTranGrids[i]->Close();
	}
	bfcOut->Close();
	if ( domainRow ) delete[] domainRow;

	//
	// Clean up
	//
	if ( fCondition ) delete[] fCondition;

	for ( int i=0; i<nRecords; i++ )
	{
		if (ppPoints[i]) delete[] ppPoints[i];
		if (data_pnts[i]) delete[] data_pnts[i];
	}
	if ( ppPoints ) delete[] ppPoints;
	if (data_pnts) delete[] data_pnts;
	
	if (pOutRow) delete[] pOutRow;
	for ( int i=0; i<numTranGrids; i++ )
	{
		if ( ppRows[i] ) delete[] ppRows[i];
	}
	if ( ppRows ) delete[] ppRows;
	for ( int i=0; i<numTranGrids; i++ )
	{
		if ( ppTranGridPaths[i] ) delete[] ppTranGridPaths[i];
	}
	if ( ppTranGridPaths ) delete[] ppTranGridPaths;

	if (theTree) delete theTree;
	if (nn_idx) delete[] nn_idx;
	if (dist) delete[] dist;

	nCurrent = 0;
	fptr("Creating GDM Density Grid", nCurrent);	
	return(true);
}



//
// Create standard GDM Density Grid (called externally)
//
bool DoStandardGDMDensity(char *pParams, 
						  char *lpDomainPath,
						  char *lpOutPath,
						  char *lpSamplePath,
						  int NumSamples,
						  bool UseNearest,
						  float fCutPoint,
						  bool DoSOS,
						  char *TimeString,
						  FPTR fptr)
{
	int nValidCells = GetValidDataCellCount(lpDomainPath, fptr);
	if ( nValidCells < 1 )
	{
		Message( "No Valid Non-Zero cells to use", "DoStandardGDMDensity" );
		return(false);
	}
	
	double SampleProp = (double)NumSamples / (double)nValidCells;

	char myBuff[BUFFLEN];
	int NumPredGrids = GetPrivateProfileInt("PREDICTORS", "NumPredictors", 0, pParams);
	if ( NumPredGrids < 1 )
	{
		Message( "No Predictor Grids to use???", "DoStandardGDMDensity" );
		return(false);
	}


	int numTranGrids = 0;
	char testBuff[64];
	for ( int i=1; i<=NumPredGrids; i++ )
	{
		sprintf( testBuff, "PredTran%d", i );
		GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, pParams);
		if ( 0 == strcmp("0",myBuff) )
		{
			//Message("does not exist", testBuff );
		}
		else
		{
			//Message("TransPath", myBuff );
			++numTranGrids;
		}
	}


	if ( numTranGrids < 1 )
	{
		Message( "No Transform Grids to use???", "DoStandardGDMDensity" );
		return(false);
	}


	//
	// allocate strings for transformed input files
	//
	char **ppTranGridPaths = new char * [numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		ppTranGridPaths[i] = new char [BUFFLEN];
	}


	//
	// Collect the Tran grid paths
	//
	int nThis = 0;
	for ( int i=1; i<=NumPredGrids; i++ )
	{
		sprintf( testBuff, "PredTran%d", i );
		GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, pParams);
		if ( 0 != strcmp("0",myBuff) )
		{
			strcpy (ppTranGridPaths[nThis++], myBuff );
		}
	}

	
	//
	// Generate the sample data points
	//
	bool GotSamplePoints = CreateSamplePointMeshForDensity( pParams, lpSamplePath, lpDomainPath, NumSamples, false, fptr );

	if (!GotSamplePoints)
	{
		Message("Cannot CreateSamplePointFile()", "DoStandardGDMDensity");
		return(false);
	}


	// 
	// determine the number of rows ( the number of records in the points file )
	//
	int nRecords = 0;

	// create the classification sample point file
	char pWorkspace[BUFFLEN];
	GetProfileString( "GDMODEL", "WorkspacePath", pWorkspace, pParams );
	char pPointsFile[BUFFLEN];
	sprintf( pPointsFile, "%s\\%s", pWorkspace, lpSamplePath );

	char buff[BUFFLEN];
	FILE *fpPoints = fopen( pPointsFile, "r+t" );
	// get the header
	fgets( buff, BUFFLEN, fpPoints );
	
	while( 1 )
	{
		if ( NULL == fgets( buff, BUFFLEN, fpPoints ) )
			break;

		++nRecords;
	}

	double **ppPoints = new double * [nRecords];
	for ( int i=0; i<nRecords; i++ )
	{
		ppPoints[i] = new double [2];
	}

	
	rewind(fpPoints);
	// get the header
	fgets( buff, BUFFLEN, fpPoints );
	char seps[] = ",\n";
	int nCurrent = 0;
	fptr("Creating GDM Density Grid (A)", nCurrent);	
	for ( int i=0; i<nRecords; i++ )
	{
		if ( i * 100 / nRecords > nCurrent )
		{
			nCurrent = i * 100 / nRecords;
			fptr("Creating GDM Density Grid (A)", nCurrent);	
		}
		
		if ( NULL == fgets( buff, BUFFLEN, fpPoints ) )
			break;

		char *p = strtok( buff, seps ); // get ID
		
		p = strtok( NULL, seps ); // get easting
		ppPoints[i][0] = atof(p);

		p = strtok( NULL, seps ); // get northing
		ppPoints[i][1] = atof(p);
	}
	fclose(fpPoints);	

	//sprintf( qqq, "nRecords: %d   numTranGrids: %d", nRecords, numTranGrids );
	//DEBUG( qqq, "INFO" );

	ANNpointArray data_pnts = annAllocPts( nRecords, numTranGrids );

	nCurrent = 0;
	fptr("Creating GDM Density Grid (B)", nCurrent);	
	for ( int i=0; i<numTranGrids; i++ )
	{
		if ( i * 100 / numTranGrids > nCurrent )
		{
			nCurrent = i * 100 / numTranGrids;
			fptr("Creating GDM Density Grid (B)", nCurrent);	
		}

		DoESRIBinaryANNColumn( data_pnts, i, ppTranGridPaths[i], ppPoints, nRecords );
	}


	// setup a vector for the condition of each sample point in the condition grid
	float *fCondition = new float [nRecords];
	GetESRIBinaryConditionVals( lpDomainPath, fCondition, ppPoints, nRecords );
	//char pConditionFile[BUFFLEN];
	//sprintf( pConditionFile, "%s\\%s", pWorkspace, "ConditionDebug.csv" );
	//FILE *fpCond = fopen( pConditionFile, "w+t" );
	//fprintf( fpCond, "_ID,X,Y,Cond\n" );
	//for ( i=0; i<nRecords; i++ )
	//{
	//	fprintf( fpCond, "%d,%lf,%lf,%f\n", i+1, ppPoints[i][0], ppPoints[i][1], fCondition[i] );
	//}
	//fclose(fpCond);


	//
	// open the header file (.hdr)
	//
	gmpath gmPath;
	char tmpBuf[500];
	strcpy(tmpBuf, gmPath.ChangeExtension(lpDomainPath, ".hdr"));
	EsriBinaryHeader *header = new EsriBinaryHeader(tmpBuf);
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	float fNoData = header->GetNoDataValue();
	double dMinX = header->GetMinX();
	double dMaxX = header->GetMaxX();
	double dMinY = header->GetMinY();
	double dMaxY = header->GetMaxY();
	double dCellSize = header->GetCellSize();

	BinaryFileClass *bfcDomain = new BinaryFileClass(gmPath.ChangeExtension(lpDomainPath, ".flt"));
	if ( !bfcDomain->IsValid() )
	{
		Message("Cannot open Domain Grid", "DoStandardGDMDensity");
		return(false);
	}
		
	//
	// create output files
	//
	header->CopyTo(gmPath.ChangeExtension(lpOutPath, ".hdr"));
	delete header;
	
	
	BinaryFileClass *bfcOut = new BinaryFileClass(gmPath.ChangeExtension(lpOutPath, ".flt"), BFC_CreateOrTrunc);
	if ( !bfcOut->IsValid() )
	{
		Message("Cannot create OutGrid", "DoStandardGDMDensity");
		return(false);
	}
	float *domainRow = new float [nCols];


	// now build the kd tree search structure
	ANNkd_tree *theTree = new ANNkd_tree ( data_pnts,		// the data points
										   nRecords,		// number of data points
										   numTranGrids );	// dimension of search space

	ANNpoint query_pt = annAllocPt( numTranGrids );
	ANNidxArray nn_idx = new ANNidx[nRecords];
	ANNdistArray dist = new ANNdist[nRecords];
	float *pOutRow = new float [nCols];


	//
	// setup for the transform grids
	//
	float **ppRows = new float * [numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		ppRows[i] = new float [nCols];
	}
	
	BinaryFileClass **bfcTranGrids = new BinaryFileClass * [numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		bfcTranGrids[i] = new BinaryFileClass(gmPath.ChangeExtension(ppTranGridPaths[i], ".flt"));
		if ( !bfcTranGrids[i]->IsValid() )
		{
			Message("Cannot open Transform Grid", "DoVanillaGDMDensity");
			return(false);
		}
	}


	//
	//  The domain can be used in a number of ways....
	//
	//  Only cells in the domain that have data are used for the density calculations.
	//	Domain values of 0.0 will have their density calculated but will NOT have any sample points placed in them.
	//  Domain values > 0.0 will have sample points placed in them.
	//  Domain values > 0.0 and <= 1.0 will also be used to weight the density calculations. 
	//  Domain values < 0.0 and > 1.0 will be ignored and set to No-Data
	//
	float f1Min = 1000000.0F;
	float f1Max = -1.0F;

	/*char pDebugFile[BUFFLEN];
	sprintf( pDebugFile, "%s\\%s", pWorkspace, gmPath.ChangeExtension(lpSamplePath,".txt"));
	FILE *fp = fopen(pDebugFile, "w+t");
	fprintf(fp, "Min,Max,Cutpoint\n" );*/

	nCurrent = 0;
	fptr("Creating GDM Density Grid", nCurrent);	
	for ( int i=0; i<nRows; i++ )
	{
		if ( i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			fptr("Creating GDM Density Grid", nCurrent);	
		}


		//
		// Initialise the output row to no-Data
		//
		for ( int j=0; j< nCols; j++ )
		{
			pOutRow[j] = fNoData;
		}


		//
		// Get the domain data row
		//
		bfcDomain->ReadFloat(domainRow, nCols);


		//
		// Get the transform grid data rows
		//
		for ( int j=0; j<numTranGrids; j++ )
		{
			bfcTranGrids[j]->ReadFloat(ppRows[j], nCols);
		}	

		
		//
		// Do the density calculations for the current row
		//		
		for ( int j=0; j<nCols; j++ )
		{
			if (( domainRow[j] >= 0.0F ) && ( domainRow[j] <= 1.0F ))
			{
				for ( int k=0; k<numTranGrids; k++ )
				{
					query_pt[k] = ppRows[k][j];
				}

				theTree->annkSearch( query_pt, nRecords, nn_idx, dist, 0 );

				double f1Num = 0.0;

				//
				// if we are using the set of NEAREST similarities
				//
				if ( UseNearest )
				{
					int nc = 0;
					int r;
					for ( r=0; r<nRecords; r++ )
					{
						//double fDist = dist[r] * dist[r];
						double fDist = dist[r];

						// only use sites with similarities > 0.9
						//if ( fDist > 2.3 )    // similarities > 0.1
						//if ( fDist > 0.7 )    // similarities > 0.5
						//if ( fDist > 0.1 )      // similarities > 0.9
						//{		
						//	break;
						//}
	
						//
						// only consider sites with user defined Cutpoint similarity
						//
						//if ( fDist > fCutPoint )				
						//	break;

						double fInc = exp(-fDist) * (double)fCondition[nn_idx[r]];

						//if ( DoSOS )
						//{
						//	// square the distance for this version
						//	fInc *= fInc;
						//}

						// increment the running total
						f1Num += fInc;	
					}

					/*if (r<nRecords)
						nc = r;
					else
						nc = nRecords;

					fprintf(fp, "%lf,%lf,%d\n", dist[0], dist[nc], nc);*/
				}

				//
				// else we are using the set of FURTHEST similarities
				//
				else
				{
					Message("Not running this code at the moment for FLOAT domains for FURTHEST similarities", "DoStandardGDMDensity" );
					return(false);

					for ( int r=nRecords-1; r>=0; r-- )
					{
						float fDist = (float)dist[r];
	
						//
						// only consider sites with user defined Cutpoint similarity
						//
						if ( fDist < fCutPoint )				
							break;

						float fInc = (float)(exp(-fDist)) * fCondition[r];

						if ( DoSOS )
						{
							// square the distance for this version
							fInc *= fInc;
						}

						// increment the running total
						f1Num += fInc;	
					}
				}

					
				// set density value and weight by domain, then divide by the number of non-zero cells
				pOutRow[j] = (float)(f1Num / SampleProp);

				// update min and max vals
				if ( f1Min > pOutRow[j] ) f1Min = pOutRow[j];
				if ( f1Max < pOutRow[j] ) f1Max = pOutRow[j];

			}  // if ( domainRow[j] > -1.0F )
		} // for ( j=0; j<nCols; j++ )

		bfcOut->WriteFloat(pOutRow, nCols);
	}

	//fclose(fp);
	
	bfcDomain->Close();
	for ( int i=0; i<numTranGrids; i++ )
	{
		bfcTranGrids[i]->Close();
	}
	bfcOut->Close();
	if ( domainRow ) delete[] domainRow;


	//
	// write the metadata for the density calculation
	//
	char myFilePath[BUFFLEN];
	sprintf(myFilePath, "%s_meta.txt", lpOutPath);
	FILE *fpMeta = fopen(myFilePath, "w+t");
	fprintf(fpMeta, "/////////////////////////////////////////////////////////////////////////////////////////////////\n");
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "Metadata for GDM Density Grid Calculation: %s\n", TimeString);
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "/////////////////////////////////////////////////////////////////////////////////////////////////\n");
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "Classification Domain: %s\n", lpDomainPath);
	fprintf(fpMeta, "Output Path: %s\n", lpOutPath);
	fprintf(fpMeta, "Sample Path: %s\n", lpSamplePath);
	fprintf(fpMeta, "Number of Samples: %d\n", NumSamples);
	//fprintf(fpMeta, "Use Nearest: %d\n", UseNearest);
	//fprintf(fpMeta, "Cutpoint: %f\n", fCutPoint);
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "/////////////////////////////////////////////////////////////////////////////////////////////////\n");


	//
	// Clean up
	//
	if ( fCondition ) delete[] fCondition;

	for ( int i=0; i<nRecords; i++ )
	{
		if (ppPoints[i]) delete[] ppPoints[i];
	}
	if ( ppPoints ) delete[] ppPoints;
	
	if (pOutRow) delete[] pOutRow;
	for ( int i=0; i<numTranGrids; i++ )
	{
		if ( ppRows[i] ) delete[] ppRows[i];
	}
	if ( ppRows ) delete[] ppRows;
	for ( int i=0; i<numTranGrids; i++ )
	{
		if ( ppTranGridPaths[i] ) delete[] ppTranGridPaths[i];
	}
	if ( ppTranGridPaths ) delete[] ppTranGridPaths;

	nCurrent = 0;
	fptr("Creating GDM Density Grid", nCurrent);	
	return(true);
}




//
//
// Count the Valid data cells in a Binary Export Grid 
//
int GetValidDataCellCount(char *GridPath, FPTR fptr)
{
	gmpath gmPath;
	int nCurrent = 0;
	fptr("Counting Valid Data Cells in Domain Grid...", nCurrent);	


	//
	// open the header file (.hdr)
	//
	char tmpBuf[500];
	strcpy(tmpBuf, gmPath.ChangeExtension(GridPath, ".hdr"));
	EsriBinaryHeader *header = new EsriBinaryHeader(tmpBuf);
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	float fNoData = header->GetNoDataValue();
	delete header;


	//
	// open the data file
	//
	strcpy(tmpBuf, gmPath.ChangeExtension(GridPath, ".flt"));
	BinaryFileClass *bfc = new BinaryFileClass(tmpBuf, BFC_ReadOnly);
	if (!bfc->IsValid())
	{
		Message("Cannot create domain grid", "GetValidDataCellCount");
		return(false);
	}


	//
	// create a grid row buffer
	//
	float *pData = new float [nCols];


	//
	// Count Valid data cells
	//
	int nValid = 0;
	fptr("Counting Valid Data Cells...", nCurrent);	
	for ( int y=0; y<nRows; y++ )
	{
		if ( y * 100 / nRows > nCurrent )
		{
			nCurrent = y * 100 / nRows;
			fptr("Counting Valid Data Cells...", nCurrent);	
		}

		// read a row of data 
		bfc->ReadFloat(pData, nCols);

		// find the Valid data cells
		for ( int x=0; x<nCols; x++ )
		{
			if ( pData[x] != fNoData )
			{
				++nValid;
			}
		} // for ( int x=0; x<nCols; x++ )

	} // for ( int y=0; y<nRows; y++ )
	nCurrent = 0;
	fptr("Counting Valid Data Cells...", nCurrent);	

	if ( pData ) delete[] pData;
	return(nValid);
}



//
// Populate ANN columns from transform grid
//
void DoESRIBinaryANNColumn( ANNpointArray data_pnts, int Index, char *GridPath, double **ppPoints, int nRecords )
{
	//
	// open the header file (.hdr)
	//
	gmpath gmPath;
	char tmpBuf[500];
	strcpy(tmpBuf, gmPath.ChangeExtension(GridPath, ".hdr"));
	EsriBinaryHeader *header = new EsriBinaryHeader(tmpBuf);
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	float fNoData = header->GetNoDataValue();
	double dMinX = header->GetMinX();
	double dMaxX = header->GetMaxX();
	double dMinY = header->GetMinY();
	double dMaxY = header->GetMaxY();
	double dCellSize = header->GetCellSize();
	delete header;
	

	strcpy(tmpBuf, gmPath.ChangeExtension(GridPath, ".flt"));
	BinaryFileClass *bfc = new BinaryFileClass(tmpBuf, BFC_ReadOnly);
	if (!bfc->IsValid())
	{
		Message("Cannot open grid", "DoESRIBinaryANNColumn");
		return;
	}	
	

	//
	// populate the column with data from the transform grid
	//
	float fVal;
	for ( int i=0; i<nRecords; i++ )
	{
		if ( ( ppPoints[i][0] > dMinX ) && ( ppPoints[i][0] < dMaxX ) && 
			 ( ppPoints[i][1] > dMinY ) && ( ppPoints[i][1] < dMaxY ) )
		{
			// calculate the linear memory block offset from the X and Y values
			int nX = int( floor( ( ppPoints[i][0] - dMinX ) / dCellSize ) );
			int nY = int( floor( ( dMaxY - ppPoints[i][1] ) / dCellSize ) );

			bfc->SeekTo(((nY * nCols) + nX) * sizeof(float));
			bfc->ReadFloat(&fVal,1);
			if ( fVal != fNoData )
			{
				data_pnts[i][Index] = fVal;
			}
		}
	}
}



//
// Get condition data from a GDM transform grid
//
void GetESRIBinaryConditionVals( char *lpDomainPath, float *fCondition, double **ppPoints, int nRecords )
{
	//
	// open the header file (.hdr)
	//
	gmpath gmPath;
	char tmpBuf[500];
	strcpy(tmpBuf, gmPath.ChangeExtension(lpDomainPath, ".hdr"));
	EsriBinaryHeader *header = new EsriBinaryHeader(tmpBuf);
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	float fNoData = header->GetNoDataValue();
	double dMinX = header->GetMinX();
	double dMaxX = header->GetMaxX();
	double dMinY = header->GetMinY();
	double dMaxY = header->GetMaxY();
	double dCellSize = header->GetCellSize();
	delete header;
	

	strcpy(tmpBuf, gmPath.ChangeExtension(lpDomainPath, ".flt"));
	BinaryFileClass *bfc = new BinaryFileClass(tmpBuf, BFC_ReadOnly);
	if (!bfc->IsValid())
	{
		Message("Cannot open grid", "GetESRIBinaryConditionVals");
		return;
	}	
	

	//
	// populate the column with data from the transform grid
	//
	float fVal;
	for ( int i=0; i<nRecords; i++ )
	{
		if ( ( ppPoints[i][0] > dMinX ) && ( ppPoints[i][0] < dMaxX ) && 
			 ( ppPoints[i][1] > dMinY ) && ( ppPoints[i][1] < dMaxY ) )
		{
			// calculate the linear memory block offset from the X and Y values
			int nX = int( floor( ( ppPoints[i][0] - dMinX ) / dCellSize ) );
			int nY = int( floor( ( dMaxY - ppPoints[i][1] ) / dCellSize ) );

			bfc->SeekTo(((nY * nCols) + nX) * sizeof(float));
			bfc->ReadFloat(&fVal,1);
			if ( fVal != fNoData )
			{
				fCondition[i] = fVal;
			}
		}
	}
}



//
// Creates a sample point file from Floating Point Domain Grid for ther Density Calculations
// 
bool CreateSamplePointMeshForDensity( char *lpParams, 
	                                  char *lpSamplePointPath, 
	                                  char *lpDomainPath, 
									  int nSamples, bool DoBatch, FPTR fptr )
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
			Message( "Cannot open Domain file for READ", "CreateSamplePointMeshForDensity" );
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
			if ( (ppDomain[y][x] != fNoData) && (ppDomain[y][x] > 0.0F) )
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
		char pWorkspace[BUFFLEN];
		GetProfileString( "GDMODEL", "WorkspacePath", pWorkspace, lpParams );
		char pPointsFile[BUFFLEN];
		sprintf( pPointsFile, "%s\\%s", pWorkspace, lpSamplePointPath );
		FILE *fp = fopen( pPointsFile, "w+t" );

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
				if ( (ppDomain[y][x] != fNoData) && (ppDomain[y][x] > 0.0F) )
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
					if ( ppDomain[ nYIndex ][ nXIndex ] > 0.0F ) 
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
// Populate DATA columns from transform grid
//
void DoESRIBinaryDataColumn( double **data_pnts, int Index, char *GridPath, double **ppPoints, int nRecords )
{
	//
	// open the header file (.hdr)
	//
	gmpath gmPath;
	char tmpBuf[500];
	strcpy(tmpBuf, gmPath.ChangeExtension(GridPath, ".hdr"));
	EsriBinaryHeader *header = new EsriBinaryHeader(tmpBuf);
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	float fNoData = header->GetNoDataValue();
	double dMinX = header->GetMinX();
	double dMaxX = header->GetMaxX();
	double dMinY = header->GetMinY();
	double dMaxY = header->GetMaxY();
	double dCellSize = header->GetCellSize();
	delete header;
	

	strcpy(tmpBuf, gmPath.ChangeExtension(GridPath, ".flt"));
	BinaryFileClass *bfc = new BinaryFileClass(tmpBuf, BFC_ReadOnly);
	if (!bfc->IsValid())
	{
		Message("Cannot open grid", "DoESRIBinaryDataColumn");
		return;
	}	
	

	//
	// populate the column with data from the transform grid
	//
	float fVal;
	for ( int i=0; i<nRecords; i++ )
	{
		if ( ( ppPoints[i][0] > dMinX ) && ( ppPoints[i][0] < dMaxX ) && 
			 ( ppPoints[i][1] > dMinY ) && ( ppPoints[i][1] < dMaxY ) )
		{
			// calculate the linear memory block offset from the X and Y values
			int nX = int( floor( ( ppPoints[i][0] - dMinX ) / dCellSize ) );
			int nY = int( floor( ( dMaxY - ppPoints[i][1] ) / dCellSize ) );

			bfc->SeekTo(((nY * nCols) + nX) * sizeof(float));
			bfc->ReadFloat(&fVal,1);
			if ( fVal != fNoData )
			{
				data_pnts[i][Index] = fVal;
			}
		}
	}
}



//
// get the environmental distance (ED)
//
double CalculateED(float **ppGridData, int ColIndex, double *pSample, int nTranGrids)
{
	double dResult = 0.0;
	for ( int i=0; i<nTranGrids; i++ )
	{
		dResult += fabs((double)ppGridData[i][ColIndex] - pSample[i]);
	}
	return(dResult);
}

