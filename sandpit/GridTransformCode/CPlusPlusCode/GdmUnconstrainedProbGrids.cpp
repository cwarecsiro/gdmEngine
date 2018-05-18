//
// GdmUnconstrainedProbGrids.cpp
//
#include "stdafx.h"
#include "GdmBinLib.h"
#include "GdmUnconstrainedProbGrids.h"
#include "GdmBinLib.h"
#include "GDMBufferSizeDefs.h"
#include "clsDoPath.h"
#include "ParamsW16.h"
#include "Message.h"
#include "clsEsriBinaryHeader.h"
#include "clsBinaryFileClass.h"



//
//  Creates a time series set of unconstrained probability grids using an external JRTraining File
//
bool DoTimeSeriesProbGridsExternal(char *pStartParams,
	                               char *pEndParams,
	                               bool UserDefinedParams,
                                   int JParam,
                                   double RParam,
                                   char *TrainingDataPath,
                                   char *OutDir,
								   int nStartYear,
                                   int nEndYear,
                                   int nIncrement,
                                   bool IgnoreEmptyClasses,
								   bool DoNormalisation,
								   char *TimeString,
                                   FPTR fptr)
{
	Message("DoTimeSeriesProbGridsExternal", "INFO");
	for ( int i=0; i<100; i++)
	{
		fptr("Testing 1,2,3...", i);
		for ( int j=0; j<0xFFFFFFFF; j++ );
	}
	return(false);
}



//
//  Creates a time series set of unconstrained probability grids using an 
//  external JRTraining File but applies a mask.
//
bool DoMaskedTimeSeriesProbGridsExternal(char *pStartParams,
	                                     char *pEndParams,
	                                     bool UserDefinedParams,
                                         int JParam,
                                         double RParam,
                                         char *TrainingDataPath,
                                         char *MaskPath,
                                         char *OutDir,
										 int nStartYear,
                                         int nEndYear,
                                         int nIncrement,
                                         bool IgnoreEmptyClasses,
										 bool DoNormalisation,
										 char *TimeString,
                                         FPTR fptr)
{
	//
	// get maximum class from training data
	//
	int MaxClass = GetMaxTrainingClass(TrainingDataPath, fptr);
	//Message(MaxClass, "MaxClass");

	//
	// setup domain file paths
	//
	char pDomainPath[BUFFLEN];
	if (false == GetTransformGridAsDomain( pStartParams, pDomainPath ))
	{
		Message("Cannot extract transform domain path", "DoMaskedTimeSeriesProbGridsExternal");
		return(false);
	}

	gmpath gmPath;
	//
	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(pDomainPath, ".hdr"));
	//
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	double ResX = header->GetCellSize();
	float fNoData = header->GetNoDataValue();
	double dMinX = header->GetXllCorner();
	double dMinY = header->GetYllCorner();
	double dMaxX = dMinX + (ResX * nCols);
	double dMaxY = dMinY + (ResX * nRows);


	//
	// setup Mask file paths
	//
	EsriBinaryHeader *maskHeader = new EsriBinaryHeader(gmPath.ChangeExtension(MaskPath, ".hdr"));
	//
	int nMaskRows = maskHeader->GetNumRows();
	int nMaskCols = maskHeader->GetNumCols();
	double MaskResX = maskHeader->GetCellSize();
	float fMaskNoData = maskHeader->GetNoDataValue();
	double dMaskMinX = maskHeader->GetXllCorner();
	double dMaskMinY = maskHeader->GetYllCorner();
	double dMaskMaxX = dMaskMinX + (MaskResX * nMaskCols);
	double dMaskMaxY = dMaskMinY + (MaskResX * nMaskRows);


	// sanity check to make sure that the mask is a subset of the domain
	if ((dMaskMinX < (dMinX-MaskResX)) || (dMaskMaxX > (dMaxX+MaskResX)) || (dMaskMinY < (dMinY-MaskResX)) || (dMaskMaxY > (dMaxY+MaskResX)))
	{
		Message("The Mask Grid is NOT a subgrid of the Domain Grid", "DoMaskedTimeSeriesProbGridsExternal");
		char qqq[64];
		sprintf(qqq, "dMinX: %lf  dMaskMinX: %lf", dMinX, dMaskMinX);
		Message(qqq, "MinX");

		sprintf(qqq, "dMaxX: %lf  dMaskMaxX: %lf", dMaxX, dMaskMaxX);
		Message(qqq, "MaxX");

		sprintf(qqq, "dMinY: %lf  dMaskMinY: %lf", dMinY, dMaskMinY);
		Message(qqq, "MinY");

		sprintf(qqq, "dMaxY: %lf  dMaskMaxY: %lf", dMaxY, dMaskMaxY);
		Message(qqq, "MaxY");
		if (header) delete header;
		if (maskHeader) delete maskHeader;
		return(false);
	}


	//
	// get the bounding rectangle indices of the mask grid relative to the domain grid
	//
	int nMinRow = int((dMaxY - dMaskMaxY) / ResX);
	int nMaxRow = int((dMaxY - dMaskMinY) / ResX);
	int nMinCol = int((dMaskMinX - dMinX) / ResX);
	int nMaxCol = int((dMaskMaxX - dMinX) / ResX);

	//
	// Now setup the ANN tree from the training data
	//
	int nCurrent = 0;
	fptr("Creating the ANN tree from the training data...", nCurrent);	
	char *scratch = new char [TABLE_ROW_BUFFSIZE];
	char seps[] = ",\n";

	FILE *fpEnv = fopen( TrainingDataPath, "r+t" );
	// assume that there is a HEADER ROW and the first column is the class,
	// the second column is the easting and the third column is the northing,
	// all the following columns correspond to the input grid file list.
	// assume that all the grid env values are floating point so the grids are also.

	// determine the number of columns
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ); // get the header	

	int nEnvCols = 0;
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv );
	char *p = strtok( scratch, seps ); // this will be the class column
	p = strtok( NULL, seps );          // this will be the easting column
	p = strtok( NULL, seps );          // this will be the northing column
	while( 1 ) 
	{
		if ( NULL == strtok( NULL, seps ) )
			break;

		else
			++nEnvCols;
	}

	// determine the number of rows
	int nEnvRows = 0;
	rewind( fpEnv );
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ); // get the header
	while( 1 )
	{
		if ( NULL == fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ) )
			break;

		else
			++nEnvRows;
	}
	

	///////////////////////////////////////////////////////////////////////////////////////
	// setup the approximate nearest neighbour data structs for the queries to follow
	ANNpoint query_pt = annAllocPt( nEnvCols );
	// this is for the whole data set
	ANNpointArray data_pnts = annAllocPts( nEnvRows, nEnvCols );

	int *pClasses = new int [ nEnvRows ];
	double *pTrainX = new double [ nEnvRows ];
	double *pTrainY = new double [ nEnvRows ];

	
	// now read the training table data into the point array
	rewind( fpEnv );
	// read the header
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv );


	// extract the class and the training data fields
	int nClasses = 0;
	int nThisRow = 0;
	while( 1 )
	{
		if ( NULL == fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ) )
			break;

		// extract the class field 
		char *t = strtok( scratch, "," );
		pClasses[ nThisRow ] = atoi( t );

		t = strtok( NULL, "," );	// skip past the easting field....
		pTrainX[ nThisRow ] = atof( t );

		t = strtok( NULL, "," );	// skip past the northing field....
		pTrainY[ nThisRow ] = atof( t );

		// now extract the rest of the environmental data
		for ( int d=0; d<nEnvCols; d++ )
		{
			t = strtok( NULL, ",\n" );

			data_pnts[nThisRow][d] = atof( t );

		}

		++nThisRow;
			
	}
	fclose(fpEnv);


	// now build the kd tree search structure
	ANNkd_tree *theTree = new ANNkd_tree ( data_pnts,		// the data points
										   nEnvRows,		// number of data points
										   nEnvCols );		// dimension of search space

	//
	// Setup for the Unconstrained Predictions
	//
	double *explu = new double [JParam];
	// setup some structures for the queries to follow
	ANNidxArray nn_idx = new ANNidx[JParam];
	ANNdistArray dist = new ANNdist[JParam];

	// setup transform grid handles
	char myBuff[BUFFLEN];
	int NumPredGrids = GetPrivateProfileInt("PREDICTORS", "NumPredictors", 0, pStartParams);
	if ( NumPredGrids < 1 )
	{
		Message( "No Predictor Grids to use???", "DoMaskedTimeSeriesProbGridsExternal" );
		return(false);
	}


	// open the domain and mask grids
	BinaryFileClass *bfc_Domain = new BinaryFileClass(gmPath.ChangeExtension(pDomainPath, ".flt"), BFC_ReadOnly);	
	//
	BinaryFileClass *bfc_Mask = new BinaryFileClass(gmPath.ChangeExtension(MaskPath, ".flt"), BFC_ReadOnly);	

	//
	// create nClasses floating point output grids
	//
	float *pRowData = new float[ nMaskCols];
	for ( int i=0; i<nMaskCols; i++ )
	{
		pRowData[i] = fMaskNoData;   // init to no-data
	}


	// setup transform grid handles
	int numTranGrids = 0;
	char testBuff[64];
	for ( int i=1; i<=NumPredGrids; i++ )
	{
		sprintf( testBuff, "PredTran%d", i );
		GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, pStartParams);
		if ( 0 == strcmp("0",myBuff) )
		{
			//printf( "%s does not exist\n", testBuff );
		}
		else
		{
			//printf( "TransPath: %s\n", myBuff );
			++numTranGrids;
		}
	}

	// check if there are geographic distance transform grids
	if (1 == GetPrivateProfileInt("GDMODEL", "UseEuclidean", 0, pStartParams))
	{
		GetPrivateProfileString("TRANSPREDS", "EuclXTran", "0", myBuff, BUFFLEN, pStartParams);
		if (0 == strcmp("0",myBuff))
		{
			Message( "EuclXTran does not exist", "DoMaskedTimeSeriesProbGridsExternal" );
			return(false);
		}
		else
		{
			++numTranGrids;
		}

		GetPrivateProfileString("TRANSPREDS", "EuclYTran", "0", myBuff, BUFFLEN, pStartParams);
		if (0 == strcmp("0",myBuff))
		{
			Message( "EuclYTran does not exist", "DoMaskedTimeSeriesProbGridsExternal" );
			return(false);
		}
		else
		{
			++numTranGrids;
		}
	}

	if ( numTranGrids < 1 )
	{
		Message( "No Transform Grids to use???", "DoMaskedTimeSeriesProbGridsExternal" );
		return(false);
	}


	//
	// allocate strings for the start step transformed input files
	//
	char **ppStartTranGridPaths = new char * [numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		ppStartTranGridPaths[i] = new char [BUFFLEN];
	}

	//
	// allocate strings for the end step transformed input files
	//
	char **ppEndTranGridPaths = new char * [numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		ppEndTranGridPaths[i] = new char [BUFFLEN];
	}


	//
	// Collect the Tran grid paths
	//
	if (1 == GetPrivateProfileInt("GDMODEL", "UseEuclidean", 0, pStartParams))
	{
		int nThis = 0;
		GetPrivateProfileString("TRANSPREDS", "EuclXTran", "0", myBuff, BUFFLEN, pStartParams);
		if ( 0 != strcmp("0",myBuff) )
		{
			strcpy (ppStartTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
		}
		GetPrivateProfileString("TRANSPREDS", "EuclYTran", "0", myBuff, BUFFLEN, pStartParams);
		if ( 0 != strcmp("0",myBuff) )
		{
			strcpy (ppStartTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
		}
		for ( int i=1; i<=NumPredGrids; i++ )
		{
			sprintf( testBuff, "PredTran%d", i );
			GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, pStartParams);
			if ( 0 != strcmp("0",myBuff) )
			{
				strcpy (ppStartTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
			}
		}
	}
	else
	{
		int nThis = 0;
		for ( int i=1; i<=NumPredGrids; i++ )
		{
			sprintf( testBuff, "PredTran%d", i );
			GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, pStartParams);
			if ( 0 != strcmp("0",myBuff) )
			{
				strcpy (ppStartTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
			}
		}
	}


	if (1 == GetPrivateProfileInt("GDMODEL", "UseEuclidean", 0, pEndParams))
	{
		int nThis = 0;
		GetPrivateProfileString("TRANSPREDS", "EuclXTran", "0", myBuff, BUFFLEN, pEndParams);
		if ( 0 != strcmp("0",myBuff) )
		{
			strcpy (ppEndTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
		}
		GetPrivateProfileString("TRANSPREDS", "EuclYTran", "0", myBuff, BUFFLEN, pEndParams);
		if ( 0 != strcmp("0",myBuff) )
		{
			strcpy (ppEndTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
		}
		for ( int i=1; i<=NumPredGrids; i++ )
		{
			sprintf( testBuff, "PredTran%d", i );
			GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, pEndParams);
			if ( 0 != strcmp("0",myBuff) )
			{
				strcpy (ppEndTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
			}
		}
	}
	else
	{
		int nThis = 0;
		for ( int i=1; i<=NumPredGrids; i++ )
		{
			sprintf( testBuff, "PredTran%d", i );
			GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, pEndParams);
			if ( 0 != strcmp("0",myBuff) )
			{
				strcpy (ppEndTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
			}
		}
	}

	
	//
	// open the start step transform grids
	//
	BinaryFileClass **bfcStartTran = new BinaryFileClass * [numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		bfcStartTran[i] = new BinaryFileClass(ppStartTranGridPaths[i]);
		if (!bfcStartTran[i]->IsValid())
		{
			Message( "Cannot open Transform Grid", "DoMaskedTimeSeriesProbGridsExternal" );
			return(false);
		}
	}

	//
	// open the end step transform grids
	//
	BinaryFileClass **bfcEndTran = new BinaryFileClass * [numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		bfcEndTran[i] = new BinaryFileClass(ppEndTranGridPaths[i]);
		if (!bfcEndTran[i]->IsValid())
		{
			Message( "Cannot open Transform Grid", "DoMaskedTimeSeriesProbGridsExternal" );
			return(false);
		}
	}


	//
	// allocate a set of row vectors for the start and end transform rows
	//
	float **ppStartRowVectors = new float *[numTranGrids];
	float **ppEndRowVectors = new float *[numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		ppStartRowVectors[i] = new float [nCols];
		ppEndRowVectors[i] = new float [nCols];
	}


	//
	// Create the No-Data vector
	//
	float *pNoDataVec = new float [nMaskCols];
	for ( int x=0; x<nMaskCols; x++ )
	{
		pNoDataVec[x] = fMaskNoData;
	}


	// allocate a vector to the masked output grid rows( one for each class )
	float **ppOutVec = new float * [MaxClass];
	for ( int i=0; i<MaxClass; i++ )
	{
		ppOutVec[i] = new float [nMaskCols];
	}
	

	//
	// for each year in the time series, create the initialised probability grids
	//
	for ( int year=nStartYear; year<=nEndYear; year+=nIncrement )
	{
		int nCurrent = 0;
		fptr("Creating initialised probability grids...", nCurrent);	
		int nThis = 0;
		char myFilePath[BUFFLEN];
		for ( int i=1; i<=MaxClass; i++ )
		{
			if ( ++nThis * 100 / MaxClass > nCurrent )
			{
				nCurrent = nThis * 100 / MaxClass;
				if ( nCurrent > 100 ) nCurrent = 100;
				fptr("Creating initialised probability grids...", nCurrent);	
			}

			//
			// ignore any classes that are not represented in the training data table
			//
			//if (IgnoreEmptyClasses && EmptyClass(i, emptyClasses, NumEmptyClasses))
			//{
			//	continue;
			//}

			//
			// Create the .FLT grid
			//
			// create a filename for the grid
			sprintf( myFilePath, "%s\\T%d\\PGrids\\P%02d%03d.flt", OutDir, year, year%100, i );
		
			// create the grid
			BinaryFileClass *bfc_Tmp = new BinaryFileClass(myFilePath, BFC_CreateOrTrunc);

			// initialise to no-data
			for ( int j=0; j<nMaskRows; j++ )
			{
				bfc_Tmp->WriteFloat( pRowData, nMaskCols );
			}

			// close the grid
			bfc_Tmp->Close();
			if (bfc_Tmp) delete bfc_Tmp;

		
			//
			// Create the .HDR grid
			//
			// create a filename for the grid
			sprintf( myFilePath, "%s\\T%d\\PGrids\\P%02d%03d.hdr", OutDir, year, year%100, i );
			maskHeader->CopyTo(myFilePath);

			// sanity check for existence of the header file
			if ( ( _access( myFilePath, 0 ) ) == -1 ) 
			{
				char myErrorString[BUFFLEN];
				sprintf(myErrorString, "Cannot find %s", myFilePath );
				Message(myErrorString, "ERROR");
				return(false);
			}
		} // for ( i=1; i<=nMaxClass; i++ )


		//////////////////////////////////////////////////////////////////////////////////////////////////
		// Finally do the unconstrained predictions
		//
		//int nThisYear = (year - nStartYear) / nIncrement;
		int nThisYear = year - nStartYear;
		BinaryFileClass **bfcOut = new BinaryFileClass * [MaxClass];
		for ( int i=1; i<=MaxClass; i++ )
		{
			//
			// ignore any classes that are not represented in the training data table
			//
			/*if (IgnoreEmptyClasses && EmptyClass(i, emptyClasses, NumEmptyClasses))
			{
				continue;
			}*/

			// open the out grid
			sprintf( myFilePath, "%s\\T%d\\PGrids\\P%02d%03d.flt", OutDir, year, year%100, i );
			bfcOut[i-1] = new BinaryFileClass(myFilePath, BFC_ReadWrite);
			if (!bfcOut[i-1]->IsValid())
			{
				Message("Cannot open outgrid", "DoMaskedTimeSeriesProbGridsExternal");
				return(false);
			}
		} // for ( i=1; i<=nMaxClass; i++ )


		nCurrent = 0;
		char banner[64];
		sprintf(banner, "Calculating probability grids for %d...", year);
		fptr(banner, nCurrent);	
		for ( int y=0; y<nRows; y++ )
		{
			// update status
			if ( y * 100 / nRows > nCurrent )
			{
				nCurrent = y * 100 / nRows;
				fptr(banner, nCurrent);	
			}


			// skip any rows that are out of the Mask extent
			if (( y < nMinRow ) || (y >= nMaxRow) ) continue;


			// read a row from the mask grid
			bfc_Mask->SeekTo((y-nMinRow) * nMaskCols * sizeof(float));
			bfc_Mask->ReadFloat(pRowData, nMaskCols);


			// read in the row data from the transform grids
			for ( int i=0; i<numTranGrids; i++ )
			{
				bfcStartTran[i]->SeekTo(y * nCols * sizeof(float));
				bfcStartTran[i]->ReadFloat(ppStartRowVectors[i], nCols);

				bfcEndTran[i]->SeekTo(y * nCols * sizeof(float));
				bfcEndTran[i]->ReadFloat(ppEndRowVectors[i], nCols);
			}


			// initialise the output rows to all No_Data
			for ( int nThisClass=0; nThisClass<MaxClass; nThisClass++ )
			{
				// initialise the output rows' probabilities to No-Data
				memcpy( ppOutVec[nThisClass], pNoDataVec, nMaskCols * sizeof(float) ); 
			}


			for ( int x=0; x<nCols; x++ )
			{
				// only if we have data and the cell is inside the mask grid
				if ((x < nMinCol) || (x >= nMaxCol)) continue;

				//if ( ppRowVectors[0][x] != fNoData )
				if (pRowData[x-nMinCol] != fMaskNoData)
				{
					/////////////////////////////////////////////////////
					/////////////////////// STEP 1 //////////////////////
					// get the nJParam nearest neighbours to this site
					for ( int c=0; c<nEnvCols; c++ ) 
					{
						//
						// This is where we apply the incremental shift in the transform data
						// to accomodate this time step in the climate change model
						//
						//float fGap = (ppEndRowVectors[c][x] - ppStartRowVectors[c][x]) / (float)nIncrement;
						float fGap = (ppEndRowVectors[c][x] - ppStartRowVectors[c][x]) / (float)(nEndYear - nStartYear);
						query_pt[ c ] = ppStartRowVectors[c][x] + (fGap * nThisYear);
					}

					theTree->annkSearch( query_pt, JParam, nn_idx, dist, 0 );

					/////////////////////////////////////////////////////
					/////////////////////// STEP 2 //////////////////////
			
					// CALCULATE SOME PROBABILITIES HERE FOR THIS POINT
					/////////////////////////////////////////////////////
					// section A
					double sumtheta = 0.0;

					for ( int m=0; m<JParam; m++ )
					{
						//dist[m] = sqrt(dist[m]);		// taking the square root of the ANN distance 
						sumtheta += dist[m];
					}

					double theta = sumtheta / (double)JParam * RParam;
					// square above value and multiply by 2 ready for section B;
					theta = 2 * theta * theta;	
					/////////////////////////////////////////////////////

					/////////////////////////////////////////////////////
					// section B
					double sumJValues = 0.0;
				
					int j;
					for ( j=0; j<JParam; j++ )
					{
						if ( theta > 0.0 )
						{
							explu[j] = exp((-(dist[j] * dist[j]))/theta);
							sumJValues += explu[j];
						}
					}

					/////////////////////////////////////////////////////
					// now apply the conditional probabilities to the row vectors
					for ( int nThisClass=0; nThisClass<MaxClass; nThisClass++ )
					{
						// initialise the output rows' probabilities to zero
						ppOutVec[nThisClass][x-nMinCol] = 0L; 
					}

					for ( j=0; j<JParam; j++ )
					{					
						// copy the probabilities over to the out vector
						if ( sumJValues > 0.0 )
							ppOutVec[(pClasses[nn_idx[j]])-1][x-nMinCol] += (float)(explu[j] / sumJValues);
					}

				} // if ( ppRowVectors[0][x] != fNoData )

			} // for ( int x=0; x<nCols; x++ )


			//
			// adjust probabilities so that they add to 1.0 for each valid data cell
			// (or 0.0 if there are no non-zero probabilities)
			//
			if (DoNormalisation)
			{
				for ( int x=0; x<nMaskCols; x++ )
				{
					if (ppOutVec[0][x] != fNoData) // a valid data cell
					{
						// sum probabilities
						float fProbSum = 0.0F;
						for ( int i=0; i<MaxClass; i++ )
						{
							fProbSum += ppOutVec[i][x];
						}

						// normalise
						if (fProbSum > 0.0F)
						{
							for ( int i=0; i<MaxClass; i++ )
							{
								ppOutVec[i][x] /= fProbSum;
							}
						}

					} // if (ppOutVec[0][x] != fNoData)
				} // for ( int x=0; x<nCols; x++ )
			} // if (DoNormalisation)


			// now write each row in the sum array to the output grid
			for ( int i=1; i<=MaxClass; i++ )
			{
				//
				// ignore any classes that are not represented in the training data table
				//
				//if (IgnoreEmptyClasses && EmptyClass(i, emptyClasses, NumEmptyClasses))
				//{
				//	continue;
				//}

				//
				// update the output grids
				//
				bfcOut[i-1]->SeekTo((y-nMinRow) * nMaskCols * sizeof(float));
				bfcOut[i-1]->WriteFloat(ppOutVec[i-1], nMaskCols);
			} // for ( i=1; i<=nMaxClass; i++ )

		} // for ( int y=0; y<nRows; y++ )

		for ( int i=1; i<=MaxClass; i++ )
		{
			//
			// ignore any classes that are not represented in the training data table
			//
			//if (IgnoreEmptyClasses && EmptyClass(i, emptyClasses, NumEmptyClasses))
			//{
			//	continue;
			//}

			// close the grid
			bfcOut[i-1]->Close();
			if (bfcOut[i-1]) delete bfcOut[i-1];
		} // for ( i=0; i<nMaxClass; i++ )
		if (bfcOut) delete[] bfcOut;
	} // for ( int i=nStartYear; i<=nEndYear; i+=nIncrement )

	if ( header) delete header;
	if ( maskHeader) delete maskHeader;


	//
	// Cleanup
	//
	for ( int i=0; i<MaxClass; i++ )
	{
		if (ppOutVec[i]) delete[] ppOutVec[i];
	}
	if (ppOutVec) delete[] ppOutVec;
	if (pNoDataVec) delete[] pNoDataVec;
	for ( int i=0; i<numTranGrids; i++ )
	{
		if ( ppStartRowVectors[i] ) delete[] ppStartRowVectors[i];
		if ( ppEndRowVectors[i] ) delete[] ppEndRowVectors[i];
	}
	if ( ppStartRowVectors ) delete[] ppStartRowVectors;
	if ( ppEndRowVectors ) delete[] ppEndRowVectors;
	for ( int i=0; i<numTranGrids; i++ )
	{
		bfcStartTran[i]->Close();
		if (bfcStartTran[i]) delete bfcStartTran[i];
		bfcEndTran[i]->Close();
		if (bfcEndTran[i]) delete bfcEndTran[i];
	}
	if (bfcEndTran) delete[] bfcEndTran;
	for ( int i=0; i<numTranGrids; i++ )
	{
		if (ppStartTranGridPaths[i]) delete[] ppStartTranGridPaths[i];
		if (ppEndTranGridPaths[i]) delete[] ppEndTranGridPaths[i];
	}
	if (ppStartTranGridPaths) delete[] ppStartTranGridPaths;
	if (ppEndTranGridPaths) delete[] ppEndTranGridPaths;
	if (pClasses) delete[] pClasses;
	if (pTrainX) delete[] pTrainX;
	if (pTrainY) delete[] pTrainY;
	if ( scratch ) delete[] scratch;
	if ( explu ) delete[] explu;
	if (theTree) delete theTree;
	return(true);
}


//
// Get maximum class from training data
//
int GetMaxTrainingClass(char *TrainingDataPath, FPTR fptr)
{
	int nMax = 0;

	FILE *fp = fopen(TrainingDataPath, "r+t");
	if (NULL == fp)
	{
		Message("Unable to open Table Path", "GetMaxTrainingClass");
		return(nMax);
	}

	int nCurrent = 0;
	fptr("GetMaxTrainingClass",nCurrent);
	char *pBuff = new char [TABLE_ROW_BUFFSIZE];
	char seps[] = ",\n";
	fgets(pBuff, TABLE_ROW_BUFFSIZE, fp); // get header
	fptr("Get Maximum Training Class", 50);
	while(1)
	{
		if (NULL == fgets(pBuff, TABLE_ROW_BUFFSIZE, fp))
			break;

		char *p = strtok(pBuff, ",\n");
		int nClass = atoi(p);

		// determine if maximum class
		if (nClass > nMax) nMax = nClass;
	}
	fptr("Get Maximum Training Class", 0);
	fclose(fp);
	if (pBuff) delete[] pBuff;
	return(nMax);
}



//
//  Creates a time series set of GDM Unsupervised Classification grids using an 
//  external JRTraining File but applies a mask.
//
bool DoMaskedTimeSeriesGDMClassificationExternal(char *pStartParams,
	                                             char *pEndParams,
                                                 char *TrainingDataPath,
                                                 char *MaskPath,
                                                 char *OutDir,
												 char *OutName,
										         int nStartYear,
                                                 int nEndYear,
                                                 int nIncrement,
										         char *TimeString,
                                                 FPTR fptr)
{
	//
	// setup domain file paths
	//
	char pDomainPath[BUFFLEN];
	if (false == GetTransformGridAsDomain( pStartParams, pDomainPath ))
	{
		Message("Cannot extract transform domain path", "DoMaskedTimeSeriesGDMClassificationExternal");
		return(false);
	}

	gmpath gmPath;
	//
	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(pDomainPath, ".hdr"));
	//
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	double ResX = header->GetCellSize();
	float fNoData = header->GetNoDataValue();
	double dMinX = header->GetXllCorner();
	double dMinY = header->GetYllCorner();
	double dMaxX = dMinX + (ResX * nCols);
	double dMaxY = dMinY + (ResX * nRows);


	//
	// setup Mask file paths
	//
	EsriBinaryHeader *maskHeader = new EsriBinaryHeader(gmPath.ChangeExtension(MaskPath, ".hdr"));
	//
	int nMaskRows = maskHeader->GetNumRows();
	int nMaskCols = maskHeader->GetNumCols();
	double MaskResX = maskHeader->GetCellSize();
	float fMaskNoData = maskHeader->GetNoDataValue();
	double dMaskMinX = maskHeader->GetXllCorner();
	double dMaskMinY = maskHeader->GetYllCorner();
	double dMaskMaxX = dMaskMinX + (MaskResX * nMaskCols);
	double dMaskMaxY = dMaskMinY + (MaskResX * nMaskRows);


	// sanity check to make sure that the mask is a subset of the domain
	if ((dMaskMinX < (dMinX-MaskResX)) || (dMaskMaxX > (dMaxX+MaskResX)) || (dMaskMinY < (dMinY-MaskResX)) || (dMaskMaxY > (dMaxY+MaskResX)))
	{
		Message("The Mask Grid is NOT a subgrid of the Domain Grid", "DoMaskedTimeSeriesGDMClassificationExternal");
		char qqq[64];
		sprintf(qqq, "dMinX: %lf  dMaskMinX: %lf", dMinX, dMaskMinX);
		Message(qqq, "MinX");

		sprintf(qqq, "dMaxX: %lf  dMaskMaxX: %lf", dMaxX, dMaskMaxX);
		Message(qqq, "MaxX");

		sprintf(qqq, "dMinY: %lf  dMaskMinY: %lf", dMinY, dMaskMinY);
		Message(qqq, "MinY");

		sprintf(qqq, "dMaxY: %lf  dMaskMaxY: %lf", dMaxY, dMaskMaxY);
		Message(qqq, "MaxY");
		if (header) delete header;
		if (maskHeader) delete maskHeader;
		return(false);
	}


	//
	// get the bounding rectangle indices of the mask grid relative to the domain grid
	//
	int nMinRow = int((dMaxY - dMaskMaxY) / ResX);
	int nMaxRow = int((dMaxY - dMaskMinY) / ResX);
	int nMinCol = int((dMaskMinX - dMinX) / ResX);
	int nMaxCol = int((dMaskMaxX - dMinX) / ResX);

	//
	// Now setup the ANN tree from the training data
	//
	int nCurrent = 0;
	fptr("Creating the ANN tree from the training data...", nCurrent);	
	char *scratch = new char [TABLE_ROW_BUFFSIZE];
	char seps[] = ",\n";

	//FILE *fpdebug = fopen(gmPath.ChangeExtension(TrainingDataPath, ".CSW"), "w+t");

	FILE *fpEnv = fopen( TrainingDataPath, "r+t" );
	// assume that there is a HEADER ROW and the first column is the class,
	// the second column is the easting and the third column is the northing,
	// all the following columns correspond to the input grid file list.
	// assume that all the grid env values are floating point so the grids are also.

	// determine the number of columns
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ); // get the header	

	int nEnvCols = 0;
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv );
	char *p = strtok( scratch, seps ); // this will be the class column
	p = strtok( NULL, seps );          // this will be the easting column
	p = strtok( NULL, seps );          // this will be the northing column
	while( 1 ) 
	{
		if ( NULL == strtok( NULL, seps ) )
			break;

		else
			++nEnvCols;
	}
		

	// determine the number of rows
	int nEnvRows = 0;
	rewind( fpEnv );
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ); // get the header
	while( 1 )
	{
		if ( NULL == fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ) )
			break;

		else
			++nEnvRows;
	}
	

	///////////////////////////////////////////////////////////////////////////////////////
	// setup the approximate nearest neighbour data structs for the queries to follow
	ANNpoint query_pt = annAllocPt( nEnvCols );
	// this is for the whole data set
	ANNpointArray data_pnts = annAllocPts( nEnvRows, nEnvCols );

	int *pClasses = new int [ nEnvRows ];
	double *pTrainX = new double [ nEnvRows ];
	double *pTrainY = new double [ nEnvRows ];

	
	// now read the training table data into the point array
	rewind( fpEnv );
	// read the header
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv );


	// extract the class and the training data fields
	int nClasses = 0;
	int nThisRow = 0;
	while( 1 )
	{
		if ( NULL == fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ) )
			break;

		// extract the class field 
		char *t = strtok( scratch, "," );
		pClasses[ nThisRow ] = atoi( t );

		//fprintf(fpdebug, "%d,", pClasses[ nThisRow ]);

		t = strtok( NULL, "," );	// skip past the easting field....
		pTrainX[ nThisRow ] = atof( t );

		//fprintf(fpdebug, "%lf,", pTrainX[ nThisRow ]);

		t = strtok( NULL, "," );	// skip past the northing field....
		pTrainY[ nThisRow ] = atof( t );

		//fprintf(fpdebug, "%lf,", pTrainY[ nThisRow ]);

		// now extract the rest of the environmental data
		for ( int d=0; d<nEnvCols; d++ )
		{
			t = strtok( NULL, ",\n" );

			data_pnts[nThisRow][d] = atof( t );

			//fprintf(fpdebug, "%lf,", data_pnts[nThisRow][d]);
			//if (d < nEnvCols-1)
			//	fprintf(fpdebug, ",");
			//else
			//	fprintf(fpdebug, "\n");
		}

		++nThisRow;
			
	}
	fclose(fpEnv);

	//if (fpdebug) 
	//	fclose(fpdebug);
	//return(false);

	// now build the kd tree search structure
	ANNkd_tree *theTree = new ANNkd_tree ( data_pnts,		// the data points
										   nEnvRows,		// number of data points
										   nEnvCols );		// dimension of search space

	//
	// Setup for the Unconstrained Predictions
	//
	double *explu = new double [nEnvRows];
	// setup some structures for the queries to follow
	ANNidxArray nn_idx = new ANNidx[nEnvRows];
	ANNdistArray dist = new ANNdist[nEnvRows];

	// setup transform grid handles
	char myBuff[BUFFLEN];
	int NumPredGrids = GetPrivateProfileInt("PREDICTORS", "NumPredictors", 0, pStartParams);
	if ( NumPredGrids < 1 )
	{
		Message( "No Predictor Grids to use???", "DoMaskedTimeSeriesGDMClassificationExternal" );
		return(false);
	}


	// open the domain and mask grids
	BinaryFileClass *bfc_Domain = new BinaryFileClass(gmPath.ChangeExtension(pDomainPath, ".flt"), BFC_ReadOnly);	
	//
	BinaryFileClass *bfc_Mask = new BinaryFileClass(gmPath.ChangeExtension(MaskPath, ".flt"), BFC_ReadOnly);	

	//
	// create nClasses floating point output grids
	//
	float *pRowData = new float[ nMaskCols];
	for ( int i=0; i<nMaskCols; i++ )
	{
		pRowData[i] = fMaskNoData;   // init to no-data
	}


	// setup transform grid handles
	int numTranGrids = 0;
	char testBuff[64];
	for ( int i=1; i<=NumPredGrids; i++ )
	{
		sprintf( testBuff, "PredTran%d", i );
		GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, pStartParams);
		if ( 0 == strcmp("0",myBuff) )
		{
			//printf( "%s does not exist\n", testBuff );
		}
		else
		{
			//printf( "TransPath: %s\n", myBuff );
			++numTranGrids;
		}
	}

	// check if there are geographic distance transform grids
	if (1 == GetPrivateProfileInt("GDMODEL", "UseEuclidean", 0, pStartParams))
	{
		GetPrivateProfileString("TRANSPREDS", "EuclXTran", "0", myBuff, BUFFLEN, pStartParams);
		if (0 == strcmp("0",myBuff))
		{
			Message( "EuclXTran does not exist", "DoMaskedTimeSeriesGDMClassificationExternal" );
			return(false);
		}
		else
		{
			++numTranGrids;
		}

		GetPrivateProfileString("TRANSPREDS", "EuclYTran", "0", myBuff, BUFFLEN, pStartParams);
		if (0 == strcmp("0",myBuff))
		{
			Message( "EuclYTran does not exist", "DoMaskedTimeSeriesGDMClassificationExternal" );
			return(false);
		}
		else
		{
			++numTranGrids;
		}
	}

	if ( numTranGrids < 1 )
	{
		Message( "No Transform Grids to use???", "DoMaskedTimeSeriesGDMClassificationExternal" );
		return(false);
	}


	//
	// allocate strings for the start step transformed input files
	//
	char **ppStartTranGridPaths = new char * [numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		ppStartTranGridPaths[i] = new char [BUFFLEN];
	}

	//
	// allocate strings for the end step transformed input files
	//
	char **ppEndTranGridPaths = new char * [numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		ppEndTranGridPaths[i] = new char [BUFFLEN];
	}


	//
	// Collect the Tran grid paths
	//
	if (1 == GetPrivateProfileInt("GDMODEL", "UseEuclidean", 0, pStartParams))
	{
		int nThis = 0;
		GetPrivateProfileString("TRANSPREDS", "EuclXTran", "0", myBuff, BUFFLEN, pStartParams);
		if ( 0 != strcmp("0",myBuff) )
		{
			strcpy (ppStartTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
		}
		GetPrivateProfileString("TRANSPREDS", "EuclYTran", "0", myBuff, BUFFLEN, pStartParams);
		if ( 0 != strcmp("0",myBuff) )
		{
			strcpy (ppStartTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
		}
		for ( int i=1; i<=NumPredGrids; i++ )
		{
			sprintf( testBuff, "PredTran%d", i );
			GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, pStartParams);
			if ( 0 != strcmp("0",myBuff) )
			{
				strcpy (ppStartTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
			}
		}
	}
	else
	{
		int nThis = 0;
		for ( int i=1; i<=NumPredGrids; i++ )
		{
			sprintf( testBuff, "PredTran%d", i );
			GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, pStartParams);
			if ( 0 != strcmp("0",myBuff) )
			{
				strcpy (ppStartTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
			}
		}
	}


	if (1 == GetPrivateProfileInt("GDMODEL", "UseEuclidean", 0, pEndParams))
	{
		int nThis = 0;
		GetPrivateProfileString("TRANSPREDS", "EuclXTran", "0", myBuff, BUFFLEN, pEndParams);
		if ( 0 != strcmp("0",myBuff) )
		{
			strcpy (ppEndTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
		}
		GetPrivateProfileString("TRANSPREDS", "EuclYTran", "0", myBuff, BUFFLEN, pEndParams);
		if ( 0 != strcmp("0",myBuff) )
		{
			strcpy (ppEndTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
		}
		for ( int i=1; i<=NumPredGrids; i++ )
		{
			sprintf( testBuff, "PredTran%d", i );
			GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, pEndParams);
			if ( 0 != strcmp("0",myBuff) )
			{
				strcpy (ppEndTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
			}
		}
	}
	else
	{
		int nThis = 0;
		for ( int i=1; i<=NumPredGrids; i++ )
		{
			sprintf( testBuff, "PredTran%d", i );
			GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, pEndParams);
			if ( 0 != strcmp("0",myBuff) )
			{
				strcpy (ppEndTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
			}
		}
	}

	
	//
	// open the start step transform grids
	//
	BinaryFileClass **bfcStartTran = new BinaryFileClass * [numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		bfcStartTran[i] = new BinaryFileClass(ppStartTranGridPaths[i]);
		if (!bfcStartTran[i]->IsValid())
		{
			Message( "Cannot open Transform Grid", "DoMaskedTimeSeriesGDMClassificationExternal" );
			return(false);
		}
	}

	//
	// open the end step transform grids
	//
	BinaryFileClass **bfcEndTran = new BinaryFileClass * [numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		bfcEndTran[i] = new BinaryFileClass(ppEndTranGridPaths[i]);
		if (!bfcEndTran[i]->IsValid())
		{
			Message( "Cannot open Transform Grid", "DoMaskedTimeSeriesGDMClassificationExternal" );
			return(false);
		}
	}


	//
	// allocate a set of row vectors for the start and end transform rows
	//
	float **ppStartRowVectors = new float *[numTranGrids];
	float **ppEndRowVectors = new float *[numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		ppStartRowVectors[i] = new float [nCols];
		ppEndRowVectors[i] = new float [nCols];
	}


	//
	// Create the No-Data vector
	//
	float *pNoDataVec = new float [nMaskCols];
	for ( int x=0; x<nMaskCols; x++ )
	{
		pNoDataVec[x] = fMaskNoData;
	}


	// allocate a vector to the masked output grid rows( one for each class )
	float *pOutVec = new float [nMaskCols];


	//
	// for each year in the time series, create the initialised probability grids
	//
	for ( int year=nStartYear; year<=nEndYear; year+=nIncrement )
	{
		int nCurrent = 0;
		fptr("Creating GDM Classification...", nCurrent);	
		int nThis = 0;
		char OutFileHDR[BUFFLEN];

		//
		// Create the .HDR grid
		//
		// create a filename for the grid
		//sprintf( OutFileHDR, "%s\\T%d\\CGrids\\%s%02d.hdr", OutDir, year, OutName, year%100 );
		sprintf( OutFileHDR, "%s\\T%d\\CGrids\\%s%02d.hdr", OutDir, year, OutName, year );
		maskHeader->CopyTo(OutFileHDR);

		// sanity check for existence of the header file
		if ( ( _access( OutFileHDR, 0 ) ) == -1 ) 
		{
			char myErrorString[BUFFLEN];
			sprintf(myErrorString, "Cannot find %s", OutFileHDR );
			Message(myErrorString, "ERROR");
			return(false);
		}

				
		//
		// Create the .FLT grid
		//
		// create a filename for the grid
		char OutFileFLT[BUFFLEN];
		//sprintf( OutFileFLT, "%s\\T%d\\CGrids\\%s%02d.flt", OutDir, year, OutName, year%100 );
		sprintf( OutFileFLT, "%s\\T%d\\CGrids\\%s%02d.flt", OutDir, year, OutName, year );
		
		// create the grid
		BinaryFileClass *bfcOut = new BinaryFileClass(OutFileFLT, BFC_CreateOrTrunc);
		if (!bfcOut->IsValid())
		{
			Message("Cannot open outgrid", "DoMaskedTimeSeriesGDMClassificationExternal");
			return(false);
		}


		//int nThisYear = (year - nStartYear) / nIncrement;
		int nThisYear = year - nStartYear;
		nCurrent = 0;
		char banner[64];
		sprintf(banner, "Creating GDM Classification for %d...", year);
		fptr(banner, nCurrent);	
		for ( int y=0; y<nRows; y++ )
		{
			// update status
			if ( y * 100 / nRows > nCurrent )
			{
				nCurrent = y * 100 / nRows;
				sprintf(banner, "Creating GDM Classification for %d (Row %d of %d)", year, y, nRows);
				fptr(banner, nCurrent);	
			}


			// skip any rows that are out of the Mask extent
			if (( y < nMinRow ) || (y >= nMaxRow) ) continue;


			// read a row from the mask grid
			bfc_Mask->SeekTo((y-nMinRow) * nMaskCols * sizeof(float));
			bfc_Mask->ReadFloat(pRowData, nMaskCols);


			// read in the row data from the transform grids
			for ( int i=0; i<numTranGrids; i++ )
			{
				bfcStartTran[i]->SeekTo(y * nCols * sizeof(float));
				bfcStartTran[i]->ReadFloat(ppStartRowVectors[i], nCols);

				bfcEndTran[i]->SeekTo(y * nCols * sizeof(float));
				bfcEndTran[i]->ReadFloat(ppEndRowVectors[i], nCols);
			}


			// initialise the output row to all No_Data
			memcpy( pOutVec, pNoDataVec, nMaskCols * sizeof(float) ); 

			for ( int x=0; x<nCols; x++ )
			{
				// only if we have data and the cell is inside the mask grid
				if ((x < nMinCol) || (x >= nMaxCol)) continue;

				if (pRowData[x-nMinCol] != fMaskNoData)
				{
					// get the nJParam nearest neighbours to this site
					for ( int c=0; c<nEnvCols; c++ ) 
					{
						//
						// This is where we apply the incremental shift in the transform data
						// to accomodate this time step in the climate change model
						//
						float fGap = (ppEndRowVectors[c][x] - ppStartRowVectors[c][x]) / (float)(nEndYear - nStartYear);
						query_pt[ c ] = ppStartRowVectors[c][x] + (fGap * nThisYear);
					}

					// find nearest neighbour
					theTree->annkSearch( query_pt, 1, nn_idx, dist, 0 );

					// assign nearest neighbour
					pOutVec[x-nMinCol] = (float)pClasses[nn_idx[0]];

				} // if ( ppRowVectors[0][x] != fNoData )
			} // for ( int x=0; x<nCols; x++ )
						
			// write output row
			//bfcOut->SeekTo((y-nMinRow) * nMaskCols * sizeof(float));
			bfcOut->WriteFloat(pOutVec, nMaskCols);

		} // for ( int y=0; y<nRows; y++ )

		// close the grid
		bfcOut->Close();
		if (bfcOut) delete bfcOut;
		
	} // for ( int i=nStartYear; i<=nEndYear; i+=nIncrement )

	if ( header) delete header;
	if ( maskHeader) delete maskHeader;


	//
	// Cleanup
	//
	if (pOutVec) delete[] pOutVec;
	if (pNoDataVec) delete[] pNoDataVec;
	for ( int i=0; i<numTranGrids; i++ )
	{
		if ( ppStartRowVectors[i] ) delete[] ppStartRowVectors[i];
		if ( ppEndRowVectors[i] ) delete[] ppEndRowVectors[i];
	}
	if ( ppStartRowVectors ) delete[] ppStartRowVectors;
	if ( ppEndRowVectors ) delete[] ppEndRowVectors;
	for ( int i=0; i<numTranGrids; i++ )
	{
		bfcStartTran[i]->Close();
		if (bfcStartTran[i]) delete bfcStartTran[i];
		bfcEndTran[i]->Close();
		if (bfcEndTran[i]) delete bfcEndTran[i];
	}
	if (bfcEndTran) delete[] bfcEndTran;
	for ( int i=0; i<numTranGrids; i++ )
	{
		if (ppStartTranGridPaths[i]) delete[] ppStartTranGridPaths[i];
		if (ppEndTranGridPaths[i]) delete[] ppEndTranGridPaths[i];
	}
	if (ppStartTranGridPaths) delete[] ppStartTranGridPaths;
	if (ppEndTranGridPaths) delete[] ppEndTranGridPaths;
	if (pClasses) delete[] pClasses;
	if (pTrainX) delete[] pTrainX;
	if (pTrainY) delete[] pTrainY;
	if ( scratch ) delete[] scratch;
	if ( explu ) delete[] explu;
	if (theTree) delete theTree;
	return(true);
}



//
//  Creates a set of unconstrained probability table for the Input Sites using an external JRTraining File 
//
bool DoUnconstrainedBinaryProbTableExternal(char *pParams,
                                            int JParam,
                                            double RParam,
											int MaxClassIndex,
                                            char *TrainingDataPath,
											char *SiteDataPath,
                                            char *OutTablePath,
											char *TimeString,
                                            FPTR fptr)
{
	gmpath gmPath;

	//
	// Now setup the ANN tree from the training data
	//
	int nCurrent = 0;
	fptr("Creating the ANN tree from the training data...", nCurrent);	
	char *scratch = new char [TABLE_ROW_BUFFSIZE];
	char seps[] = ",\n";

	FILE *fpEnv = fopen( TrainingDataPath, "r+t" );
	// assume that there is a HEADER ROW and the first column is the class,
	// the second column is the easting and the third column is the northing,
	// all the following columns correspond to the input grid file list.
	// assume that all the grid env values are floating point so the grids are also.

	// determine the number of columns
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ); // get the header	

	int nEnvCols = 0;
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv );
	char *p = strtok( scratch, seps ); // this will be the class column
	p = strtok( NULL, seps );          // this will be the easting column
	p = strtok( NULL, seps );          // this will be the northing column
	while( 1 ) 
	{
		if ( NULL == strtok( NULL, seps ) )
			break;

		else
			++nEnvCols;
	}

	// determine the number of rows
	int nEnvRows = 0;
	rewind( fpEnv );
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ); // get the header
	while( 1 )
	{
		if ( NULL == fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ) )
			break;

		else
			++nEnvRows;
	}
	

	///////////////////////////////////////////////////////////////////////////////////////
	// setup the approximate nearest neighbour data structs for the queries to follow
	ANNpoint query_pt = annAllocPt( nEnvCols );
	// this is for the whole data set
	ANNpointArray data_pnts = annAllocPts( nEnvRows, nEnvCols );

	int *pClasses = new int [ nEnvRows ];
	double *pTrainX = new double [ nEnvRows ];
	double *pTrainY = new double [ nEnvRows ];

	
	// now read the training table data into the point array
	rewind( fpEnv );
	// read the header
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv );


	// extract the class and the training data fields
	int nClasses = 0;
	int nThisRow = 0;
	while( 1 )
	{
		if ( NULL == fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ) )
			break;

		// extract the class field 
		char *t = strtok( scratch, "," );
		pClasses[ nThisRow ] = atoi( t );

		t = strtok( NULL, "," );	// skip past the easting field....
		pTrainX[ nThisRow ] = atof( t );

		t = strtok( NULL, "," );	// skip past the northing field....
		pTrainY[ nThisRow ] = atof( t );

		// now extract the rest of the environmental data
		for ( int d=0; d<nEnvCols; d++ )
		{
			t = strtok( NULL, ",\n" );

			data_pnts[nThisRow][d] = atof( t );

		}

		++nThisRow;
			
	}
	fclose(fpEnv);


	// now build the kd tree search structure
	ANNkd_tree *theTree = new ANNkd_tree ( data_pnts,		// the data points
										   nEnvRows,		// number of data points
										   nEnvCols );		// dimension of search space

	//
	// Setup for the Unconstrained Predictions
	//
	double *explu = new double [JParam];
	// setup some structures for the queries to follow
	ANNidxArray nn_idx = new ANNidx[JParam];
	ANNdistArray dist = new ANNdist[JParam];

	// setup transform grid handles
	char myBuff[BUFFLEN];
	int NumPredGrids = GetPrivateProfileInt("PREDICTORS", "NumPredictors", 0, pParams);
	if ( NumPredGrids < 1 )
	{
		Message( "No Predictor Grids to use???", "DoUnconstrainedBinaryProbTableExternal" );
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
			//printf( "%s does not exist\n", testBuff );
		}
		else
		{
			//printf( "TransPath: %s\n", myBuff );
			++numTranGrids;
		}
	}

	// check if there are geographic distance transform grids
	if (1 == GetPrivateProfileInt("GDMODEL", "UseEuclidean", 0, pParams))
	{
		GetPrivateProfileString("TRANSPREDS", "EuclXTran", "0", myBuff, BUFFLEN, pParams);
		if (0 == strcmp("0",myBuff))
		{
			Message( "EuclXTran does not exist", "ERROR" );
			return(false);
		}
		else
		{
			++numTranGrids;
		}

		GetPrivateProfileString("TRANSPREDS", "EuclYTran", "0", myBuff, BUFFLEN, pParams);
		if (0 == strcmp("0",myBuff))
		{
			Message( "EuclYTran does not exist", "ERROR" );
			return(false);
		}
		else
		{
			++numTranGrids;
		}
	}

	if ( numTranGrids < 1 )
	{
		Message( "No Transform Grids to use???", "DoUnconstrainedBinaryProbTableExternal" );
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
	if (1 == GetPrivateProfileInt("GDMODEL", "UseEuclidean", 0, pParams))
	{
		nThis = 0;
		GetPrivateProfileString("TRANSPREDS", "EuclXTran", "0", myBuff, BUFFLEN, pParams);
		if ( 0 != strcmp("0",myBuff) )
		{
			strcpy (ppTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
		}
		GetPrivateProfileString("TRANSPREDS", "EuclYTran", "0", myBuff, BUFFLEN, pParams);
		if ( 0 != strcmp("0",myBuff) )
		{
			strcpy (ppTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
		}
		for ( int i=1; i<=NumPredGrids; i++ )
		{
			sprintf( testBuff, "PredTran%d", i );
			GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, pParams);
			if ( 0 != strcmp("0",myBuff) )
			{
				strcpy (ppTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
			}
		}
	}
	else
	{
		nThis = 0;
		for ( int i=1; i<=NumPredGrids; i++ )
		{
			sprintf( testBuff, "PredTran%d", i );
			GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, pParams);
			if ( 0 != strcmp("0",myBuff) )
			{
				strcpy (ppTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
			}
		}
	}


	//
	// open the transform grids
	//
	BinaryFileClass **bfcTran = new BinaryFileClass * [numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		bfcTran[i] = new BinaryFileClass(ppTranGridPaths[i]);
		if (!bfcTran[i]->IsValid())
		{
			Message( "Cannot open Transform Grid", "DoUnconstrainedBinaryProbTableExternal" );
			return(false);
		}
	}


	//
	// Use the first transform grid to get the grid metrics from...
	//
	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(ppTranGridPaths[0], ".hdr"));
	//
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	double ResX = header->GetCellSize();
	float fNoData = header->GetNoDataValue();
	double dMinX = header->GetXllCorner();
	double dMinY = header->GetYllCorner();
	double dMaxX = dMinX + (ResX * nCols);
	double dMaxY = dMinY + (ResX * nRows);

	
	//
	// allocate a set of row transform grid row vectors
	//
	float *pSiteTransformVector = new float [numTranGrids];


	//
	// allocate the output probability vector
	//
	float *pOutVec = new float [MaxClassIndex];
	for ( int i=0; i<MaxClassIndex; i++ )
	{
		pOutVec[i] = 0.0F;
	}


	//
	// Create the output table file
	//
	FILE *fpOut = fopen(OutTablePath, "w+t");
	fprintf(fpOut, "X,Y,");
	for ( int i=1; i<=MaxClassIndex; i++ )
	{
		fprintf(fpOut, "Prob_%d", i);
		if (i < MaxClassIndex)
			fprintf(fpOut, ",");
		else
			fprintf(fpOut, "\n");
	}
	

	//
	// Open the Site Table and count the rows
	//
	char sitebuff[BUFFLEN];
	FILE *fpSites = fopen(SiteDataPath, "r+t");
	fgets(sitebuff,BUFFLEN,fpSites); // get header
	int nSites = 0;
	while(1)
	{
		if ( NULL == fgets(sitebuff,BUFFLEN,fpSites) )
			break;

		++nSites;
	}
	rewind(fpSites);
	fgets(sitebuff,BUFFLEN,fpSites); // get header
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Finally do the unconstrained predictions
	//
	long MaxOffset = nRows * nCols * sizeof(float);
	double dX,dY;
	//char qqq[64];
	nCurrent = 0;
	fptr("Calculating probability grids...", nCurrent);	
	for ( int y=0; y<nSites; y++ )
	{
		// update status
		if ( y * 100 / nSites > nCurrent )
		{
			nCurrent = y * 100 / nSites;
			fptr("Calculating probability grids...", nCurrent);	
		}

		// initialise the output row to 0.0F
		for ( int nThisClass=0; nThisClass<MaxClassIndex; nThisClass++ )
		{
			pOutVec[nThisClass] = 0.0F; 
		}

		// get the sites for this row
		fgets(sitebuff,BUFFLEN,fpSites);
		char *p = strtok(sitebuff, ",\n");
		dX = atof(p);
		p = strtok(NULL, ",\n");
		dY = atof(p);
		
		// calculate the offset for this site
		long lOffset = (long((dMaxY - dY) / ResX) * nCols) + long((dX - dMinX)/ResX);

		// if the site is valid...
		if (( lOffset >= 0L ) && (lOffset < MaxOffset))
		{
			// test for data validity using the first transform grid
			bfcTran[0]->SeekTo(lOffset * sizeof(float));
			float fTest;
			bfcTran[0]->ReadFloat(&fTest, 1);
			if (fTest == fNoData)
			{
				continue;
			}

			// read in the row data from the transform grids
			for ( int i=0; i<numTranGrids; i++ )
			{
				bfcTran[i]->SeekTo(lOffset * sizeof(float));
				bfcTran[i]->ReadFloat(&pSiteTransformVector[i], 1);
			}
						
			/////////////////////////////////////////////////////
			/////////////////////// STEP 1 //////////////////////
			// get the nJParam nearest neighbours to this site
			for ( int c=0; c<nEnvCols; c++ ) 
			{
				query_pt[ c ] = pSiteTransformVector[c];
			}
			theTree->annkSearch( query_pt, JParam, nn_idx, dist, 0 );

			/////////////////////////////////////////////////////
			/////////////////////// STEP 2 //////////////////////
			//
			// CALCULATE SOME PROBABILITIES HERE FOR THIS POINT
			/////////////////////////////////////////////////////
			// section A
			double sumtheta = 0.0;

			for ( int m=0; m<JParam; m++ )
			{
				sumtheta += dist[m];
			}

			double theta = sumtheta / (double)JParam * RParam;
			// square above value and multiply by 2 ready for section B;
			theta = 2 * theta * theta;	
			/////////////////////////////////////////////////////

			/////////////////////////////////////////////////////
			// section B
			double sumJValues = 0.0;
				
			int j;
			for ( j=0; j<JParam; j++ )
			{
				if ( theta > 0.0 )
				{
					explu[j] = exp((-(dist[j] * dist[j]))/theta);
					sumJValues += explu[j];
				}
			}

			/////////////////////////////////////////////////////
			for ( j=0; j<JParam; j++ )
			{					
				// copy the probabilities over to the out vector
				if ( sumJValues > 0.0 )
					pOutVec[(pClasses[nn_idx[j]])-1] += (float)(explu[j] / sumJValues);
			}
		} // if (( lOffset >= 0L ) && (lOffset < MaxOffset))

		// now write each site prediction to the output table file
		fprintf(fpOut, "%lf,%lf,", dX, dY );
		for ( int i=1; i<=MaxClassIndex; i++ )
		{
			fprintf(fpOut, "%f", pOutVec[i-1]);
			if (i < MaxClassIndex)
				fprintf(fpOut, ",");
			else
				fprintf(fpOut, "\n");
		} // for ( i=1; i<=nMaxClass; i++ )

	} // for ( int y=0; y<nRows; y++ )


	//
	// clean up
	//	
	fclose(fpOut);
	fclose(fpSites);

	if (pClasses) delete[] pClasses;
	if (pTrainX) delete[] pTrainX;
	if (pTrainY) delete[] pTrainY;
	if ( scratch ) delete[] scratch;
	if ( explu ) delete[] explu;
	if (theTree) delete theTree;

	for ( int i=0; i<numTranGrids; i++ )
	{
		bfcTran[i]->Close();
		if (bfcTran[i]) delete bfcTran[i];
	}
	if (bfcTran) delete[] bfcTran;
	for ( int i=0; i<numTranGrids; i++ )
	{
		if (ppTranGridPaths[i]) delete[] ppTranGridPaths[i];
	}
	if (ppTranGridPaths) delete[] ppTranGridPaths;
	if ( pOutVec ) delete[] pOutVec;
	return(true);
}



//
//  Creates a set of unconstrained probability grids using an external JRTraining File
//
bool DoUnconstrainedBinaryProbGridsExternal(char *pParams,
	                                        bool UserDefinedParams,
                                            int JParam,
                                            double RParam,
											int MaxClassIndex,
                                            char *TrainingDataPath,
                                            char *OutDir,
											bool IgnoreEmptyClasses,
											bool DoNormalisation,
											char *TimeString,
                                            FPTR fptr)
{
	//
	// open the training data table and extract the maximum class index
	//
	int NumEmptyClasses = 0;
	int *emptyClasses = GetEmptyClasses(&NumEmptyClasses, TrainingDataPath, MaxClassIndex);

	//
	// setup domain file paths
	//	
	char pDomainPath[BUFFLEN];
	if (false == GetTransformGridAsDomain( pParams, pDomainPath ))
	{
		Message("Cannot extract transform domain path", "DoUnconstrainedBinaryProbGridsExternal");
		if (emptyClasses) delete[] emptyClasses;
		return(false);
	}

	gmpath gmPath;
	//
	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(pDomainPath, ".hdr"));	

	//
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	double ResX = header->GetCellSize();
	float fNoData = header->GetNoDataValue();
	double dMinX = header->GetXllCorner();
	double dMinY = header->GetYllCorner();
	double dMaxX = dMinX + (ResX * nCols);
	double dMaxY = dMinY + (ResX * nRows);

	//
	// create nClasses floating point output grids
	//
	float *pRowData = new float[ nCols];
	for ( int i=0; i<nCols; i++ )
	{
		pRowData[i] = fNoData;   // init to no-data
	}

	int nCurrent = 0;
	fptr("Creating initialised probability grids...", nCurrent);	
	int nThis = 0;
	char myFilePath[BUFFLEN];
	for ( int i=1; i<=MaxClassIndex; i++ )
	{
		if ( ++nThis * 100 / MaxClassIndex > nCurrent )
		{
			nCurrent = nThis * 100 / MaxClassIndex;
			if ( nCurrent > 100 ) nCurrent = 100;
			fptr("Creating initialised probability grids...", nCurrent);	
		}

		//
		// ignore any classes that are not represented in the training data table
		//
		if (IgnoreEmptyClasses && EmptyClass(i, emptyClasses, NumEmptyClasses))
		{
			continue;
		}

		//
		// Create the .GRI grid
		//
		// create a filename for the grid
		sprintf( myFilePath, "%s\\uncon%03d.flt", OutDir, i );
		
		// create the grid
		BinaryFileClass *bfc_Tmp = new BinaryFileClass(myFilePath, BFC_CreateOrTrunc);

		// initialise to no-data
		for ( int j=0; j<nRows; j++ )
		{
			bfc_Tmp->WriteFloat( pRowData, nCols );
		}

		// close the grid
		bfc_Tmp->Close();
		if (bfc_Tmp) delete bfc_Tmp;

		
		//
		// Create the .GRD grid
		//
		// create a filename for the grid
		sprintf( myFilePath, "%s\\uncon%03d.hdr", OutDir, i );
		header->CopyTo(myFilePath);

		// sanity check for existence of the header file
		if ( ( _access( myFilePath, 0 ) ) == -1 ) 
		{
			char myErrorString[BUFFLEN];
			sprintf(myErrorString, "Cannot find %s", myFilePath );
			Message(myErrorString, "ERROR");
			return(false);
		}
	} // for ( i=1; i<=nMaxClass; i++ )

	if ( header) delete header;
	if ( pRowData ) delete[] pRowData;

	//
	// Now setup the ANN tree from the training data
	//
	nCurrent = 0;
	fptr("Creating the ANN tree from the training data...", nCurrent);	
	char *scratch = new char [TABLE_ROW_BUFFSIZE];
	char seps[] = ",\n";

	FILE *fpEnv = fopen( TrainingDataPath, "r+t" );
	// assume that there is a HEADER ROW and the first column is the class,
	// the second column is the easting and the third column is the northing,
	// all the following columns correspond to the input grid file list.
	// assume that all the grid env values are floating point so the grids are also.

	// determine the number of columns
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ); // get the header	

	int nEnvCols = 0;
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv );
	char *p = strtok( scratch, seps ); // this will be the class column
	p = strtok( NULL, seps );          // this will be the easting column
	p = strtok( NULL, seps );          // this will be the northing column
	while( 1 ) 
	{
		if ( NULL == strtok( NULL, seps ) )
			break;

		else
			++nEnvCols;
	}

	// determine the number of rows
	int nEnvRows = 0;
	rewind( fpEnv );
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ); // get the header
	while( 1 )
	{
		if ( NULL == fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ) )
			break;

		else
			++nEnvRows;
	}
	

	///////////////////////////////////////////////////////////////////////////////////////
	// setup the approximate nearest neighbour data structs for the queries to follow
	ANNpoint query_pt = annAllocPt( nEnvCols );
	// this is for the whole data set
	ANNpointArray data_pnts = annAllocPts( nEnvRows, nEnvCols );

	int *pClasses = new int [ nEnvRows ];
	double *pTrainX = new double [ nEnvRows ];
	double *pTrainY = new double [ nEnvRows ];

	
	// now read the training table data into the point array
	rewind( fpEnv );
	// read the header
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv );


	// extract the class and the training data fields
	int nClasses = 0;
	int nThisRow = 0;
	while( 1 )
	{
		if ( NULL == fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ) )
			break;

		// extract the class field 
		char *t = strtok( scratch, "," );
		pClasses[ nThisRow ] = atoi( t );

		t = strtok( NULL, "," );	// skip past the easting field....
		pTrainX[ nThisRow ] = atof( t );

		t = strtok( NULL, "," );	// skip past the northing field....
		pTrainY[ nThisRow ] = atof( t );

		// now extract the rest of the environmental data
		for ( int d=0; d<nEnvCols; d++ )
		{
			t = strtok( NULL, ",\n" );

			data_pnts[nThisRow][d] = atof( t );

		}

		++nThisRow;
			
	}
	fclose(fpEnv);


	// now build the kd tree search structure
	ANNkd_tree *theTree = new ANNkd_tree ( data_pnts,		// the data points
										   nEnvRows,		// number of data points
										   nEnvCols );		// dimension of search space

	//
	// Setup for the Unconstrained Predictions
	//
	double *explu = new double [JParam];
	// setup some structures for the queries to follow
	ANNidxArray nn_idx = new ANNidx[JParam];
	ANNdistArray dist = new ANNdist[JParam];

	// setup transform grid handles
	char myBuff[BUFFLEN];
	int NumPredGrids = GetPrivateProfileInt("PREDICTORS", "NumPredictors", 0, pParams);
	if ( NumPredGrids < 1 )
	{
		Message( "No Predictor Grids to use???", "DoMaskedUnconstrainedBinaryProbGridsExternal" );
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
			//printf( "%s does not exist\n", testBuff );
		}
		else
		{
			//printf( "TransPath: %s\n", myBuff );
			++numTranGrids;
		}
	}

	// check if there are geographic distance transform grids
	if (1 == GetPrivateProfileInt("GDMODEL", "UseEuclidean", 0, pParams))
	{
		GetPrivateProfileString("TRANSPREDS", "EuclXTran", "0", myBuff, BUFFLEN, pParams);
		if (0 == strcmp("0",myBuff))
		{
			Message( "EuclXTran does not exist", "ERROR" );
			return(false);
		}
		else
		{
			++numTranGrids;
		}

		GetPrivateProfileString("TRANSPREDS", "EuclYTran", "0", myBuff, BUFFLEN, pParams);
		if (0 == strcmp("0",myBuff))
		{
			Message( "EuclYTran does not exist", "ERROR" );
			return(false);
		}
		else
		{
			++numTranGrids;
		}
	}

	if ( numTranGrids < 1 )
	{
		Message( "No Transform Grids to use???", "DoMaskedUnconstrainedBinaryProbGridsExternal" );
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
	if (1 == GetPrivateProfileInt("GDMODEL", "UseEuclidean", 0, pParams))
	{
		nThis = 0;
		GetPrivateProfileString("TRANSPREDS", "EuclXTran", "0", myBuff, BUFFLEN, pParams);
		if ( 0 != strcmp("0",myBuff) )
		{
			strcpy (ppTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
		}
		GetPrivateProfileString("TRANSPREDS", "EuclYTran", "0", myBuff, BUFFLEN, pParams);
		if ( 0 != strcmp("0",myBuff) )
		{
			strcpy (ppTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
		}
		for ( int i=1; i<=NumPredGrids; i++ )
		{
			sprintf( testBuff, "PredTran%d", i );
			GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, pParams);
			if ( 0 != strcmp("0",myBuff) )
			{
				strcpy (ppTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
			}
		}
	}
	else
	{
		nThis = 0;
		for ( int i=1; i<=NumPredGrids; i++ )
		{
			sprintf( testBuff, "PredTran%d", i );
			GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, pParams);
			if ( 0 != strcmp("0",myBuff) )
			{
				strcpy (ppTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
			}
		}
	}


	//
	// open the transform grids
	//
	BinaryFileClass **bfcTran = new BinaryFileClass * [numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		bfcTran[i] = new BinaryFileClass(ppTranGridPaths[i]);
		if (!bfcTran[i]->IsValid())
		{
			Message( "Cannot open Transform Grid", "DoMaskedUnconstrainedBinaryProbGridsExternal" );
			return(false);
		}
	}


	//
	// allocate a set of row vectors
	//
	float **ppRowVectors = new float *[numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		ppRowVectors[i] = new float [nCols];
	}


	// allocate a vector to the masked output grid rows( one for each class )
	float **ppOutVec = new float * [MaxClassIndex];
	for ( int i=0; i<MaxClassIndex; i++ )
	{
		ppOutVec[i] = new float [nCols];
		for ( int q=0; q<nCols; q++ ) ppOutVec[i][q] = fNoData;
	}


	float *pNoDataVec = new float [nCols];
	for ( int x=0; x<nCols; x++ )
	{
		pNoDataVec[x] = fNoData;
	}


	float *pMaxProbs = new float [MaxClassIndex];
	BinaryFileClass **bfcOut = new BinaryFileClass * [MaxClassIndex];
	for ( int i=1; i<=MaxClassIndex; i++ )
	{
		//
		// Init probs to zero
		//
		pMaxProbs[i-1] = 0.0;
		
		//
		// ignore any classes that are not represented in the training data table
		//
		if (IgnoreEmptyClasses && EmptyClass(i, emptyClasses, NumEmptyClasses))
		{
			continue;
		}

		// open the out grid
		sprintf( myFilePath, "%s\\uncon%03d.flt", OutDir, i );
		bfcOut[i-1] = new BinaryFileClass(myFilePath, BFC_ReadWrite);
		if (!bfcOut[i-1]->IsValid())
		{
			Message("Cannot open outgrid", "DoMaskedUnconstrainedBinaryProbGridsExternal");
			return(false);
		}
	} // for ( i=1; i<=nMaxClass; i++ )


	//
	// Create a metadata file detailing the creation parameters and the output grid metrics
	//
	sprintf(myFilePath, "%s\\metadata.txt", OutDir);
	FILE *fpMeta = fopen(myFilePath, "w+t");
	fprintf(fpMeta, "/////////////////////////////////////////////////////////////////////////////////////////////////\n");
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "Metadata for unconstrained probability grid creation: %s\n", TimeString);
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "/////////////////////////////////////////////////////////////////////////////////////////////////\n");
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "InputTable: %s\n", TrainingDataPath);
	fprintf(fpMeta, "Output Directory: %s\n", OutDir);
	fprintf(fpMeta, "Maximum Class Index: %d\n", MaxClassIndex);
	fprintf(fpMeta, "J Param: %d       R Param: %lf\n", JParam, RParam);
	if (IgnoreEmptyClasses)
		fprintf(fpMeta, "Empty Classes are excluded\n");
	else
		fprintf(fpMeta, "Empty Classes are included\n");
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "/////////////////////////////////////////////////////////////////////////////////////////////////\n");
	

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Finally do the unconstrained predictions
	//
	//char qqq[64];
	nCurrent = 0;
	fptr("Calculating probability grids...", nCurrent);	
	for ( int y=0; y<nRows; y++ )
	{
		// update status
		if ( y * 100 / nRows > nCurrent )
		{
			nCurrent = y * 100 / nRows;
			//sprintf(qqq, "Calculating probability grids: Row %d of %d", y, nRows);
			//fptr(qqq, nCurrent);	
			fptr("Calculating probability grids...", nCurrent);	
		}


		// read in the row data from the transform grids
		for ( int i=0; i<numTranGrids; i++ )
		{
			bfcTran[i]->SeekTo(y * nCols * sizeof(float));
			bfcTran[i]->ReadFloat(ppRowVectors[i], nCols);
		}


		// initialise the output rows to all No_Data
		for ( int nThisClass=0; nThisClass<MaxClassIndex; nThisClass++ )
		{
			// initialise the output rows' probabilities to No-Data
			memcpy( ppOutVec[nThisClass], pNoDataVec, nCols * sizeof(float) ); 
		}


		for ( int x=0; x<nCols; x++ )
		{
			if ( ppRowVectors[0][x] != fNoData )
			{
				/////////////////////////////////////////////////////
				/////////////////////// STEP 1 //////////////////////
				// get the nJParam nearest neighbours to this site
				for ( int c=0; c<nEnvCols; c++ ) 
				{
					query_pt[ c ] = ppRowVectors[c][x];
				}

				theTree->annkSearch( query_pt, JParam, nn_idx, dist, 0 );

				/////////////////////////////////////////////////////
				/////////////////////// STEP 2 //////////////////////
			
				// CALCULATE SOME PROBABILITIES HERE FOR THIS POINT
				/////////////////////////////////////////////////////
				// section A
				double sumtheta = 0.0;

				for ( int m=0; m<JParam; m++ )
				{
					//dist[m] = sqrt(dist[m]);		// taking the square root of the ANN distance 
					sumtheta += dist[m];
				}

				double theta = sumtheta / (double)JParam * RParam;
				// square above value and multiply by 2 ready for section B;
				theta = 2 * theta * theta;	
				/////////////////////////////////////////////////////

				/////////////////////////////////////////////////////
				// section B
				double sumJValues = 0.0;
				
				int j;
				for ( j=0; j<JParam; j++ )
				{
					if ( theta > 0.0 )
					{
						explu[j] = exp((-(dist[j] * dist[j]))/theta);
						sumJValues += explu[j];
					}
				}

				/////////////////////////////////////////////////////
				// now apply the conditional probabilities to the row vectors
				for ( int nThisClass=0; nThisClass<MaxClassIndex; nThisClass++ )
				{
					// initialise the output rows' probabilities to zero
					ppOutVec[nThisClass][x] = 0L; 
				}

				for ( j=0; j<JParam; j++ )
				{					
					// copy the probabilities over to the out vector
					if ( sumJValues > 0.0 )
						ppOutVec[(pClasses[nn_idx[j]])-1][x] += (float)(explu[j] / sumJValues);
				}

			} // if ( ppRowVectors[0][x] != fNoData )

		} // for ( int x=0; x<nCols; x++ )


		//
		// adjust probabilities so that they add to 1.0 for each valid data cell
		// (or 0.0 if there are no non-zero probabilities)
		//
		if (DoNormalisation)
		{
			for ( int x=0; x<nCols; x++ )
			{
				if (ppOutVec[0][x] != fNoData) // a valid data cell
				{
					// sum probabilities
					float fProbSum = 0.0F;
					for ( int i=0; i<MaxClassIndex; i++ )
					{
						fProbSum += ppOutVec[i][x];
					}

					// normalise
					if (fProbSum > 0.0F)
					{
						for ( int i=0; i<MaxClassIndex; i++ )
						{
							ppOutVec[i][x] /= fProbSum;
						}
					}

				} // if (ppOutVec[0][x] != fNoData)
			} // for ( int x=0; x<nCols; x++ )
		} // if (DoNormalisation)


		// now write each row in the sum array to the output grid
		for ( int i=1; i<=MaxClassIndex; i++ )
		{
			//
			// ignore any classes that are not represented in the training data table
			//
			if (IgnoreEmptyClasses && EmptyClass(i, emptyClasses, NumEmptyClasses))
			{
				continue;
			}

			//
			// update the output grids
			//
			bfcOut[i-1]->SeekTo(y * nCols * sizeof(float));
			bfcOut[i-1]->WriteFloat(ppOutVec[i-1], nCols);
		} // for ( i=1; i<=nMaxClass; i++ )


		//
		// update the min and max probibilities
		//
		for ( int i=0; i<MaxClassIndex; i++ )
		{
			for (int j=0; j<nCols; j++)
			{
				if ( ppOutVec[i][j] != fNoData )
				{
					if (pMaxProbs[i] < ppOutVec[i][j]) pMaxProbs[i] = ppOutVec[i][j];
				}
			}
		}

	} // for ( int y=0; y<nRows; y++ )


	//
	// Finish writing the metadata 
	//
	for ( int i=0; i<MaxClassIndex; i++ )
	{
		fprintf(fpMeta, "Class: %d   MaxProb: %f\n", i+1, pMaxProbs[i]);
	}
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "/////////////////////////////////////////////////////////////////////////////////////////////////\n");
	if (fpMeta) fclose(fpMeta);


	//
	// clean up
	//	
	if (pMaxProbs) delete[] pMaxProbs;

	for ( int i=1; i<=MaxClassIndex; i++ )
	{
		//
		// ignore any classes that are not represented in the training data table
		//
		if (IgnoreEmptyClasses && EmptyClass(i, emptyClasses, NumEmptyClasses))
		{
			continue;
		}

		// close the grid
		bfcOut[i-1]->Close();
		if (bfcOut[i-1]) delete bfcOut[i-1];
	} // for ( i=0; i<0nMaxClass; i++ )

	if (bfcOut) delete[] bfcOut;
	if (emptyClasses) delete[] emptyClasses;
	if (pClasses) delete[] pClasses;
	if (pTrainX) delete[] pTrainX;
	if (pTrainY) delete[] pTrainY;
	if ( scratch ) delete[] scratch;
	if ( explu ) delete[] explu;
	if ( pNoDataVec ) delete[] pNoDataVec;
	if (theTree) delete theTree;

	for ( int i=0; i<numTranGrids; i++ )
	{
		bfcTran[i]->Close();
		if (bfcTran[i]) delete bfcTran[i];
	}
	if (bfcTran) delete[] bfcTran;
	for ( int i=0; i<numTranGrids; i++ )
	{
		if (ppTranGridPaths[i]) delete[] ppTranGridPaths[i];
	}
	if (ppTranGridPaths) delete[] ppTranGridPaths;
	for ( int i=0; i<numTranGrids; i++ )
	{
		if ( ppRowVectors[i] ) delete[] ppRowVectors[i];
	}
	if ( ppRowVectors ) delete[] ppRowVectors;
	for ( int i=0; i<MaxClassIndex; i++ )
	{
		if (ppOutVec[i]) delete[] ppOutVec[i];
	}
	if ( ppOutVec ) delete[] ppOutVec;
	return(true);
}



//
//  Creates a set of unconstrained probability grids using an external JRTraining File but applies a mask
//
bool DoMaskedUnconstrainedBinaryProbGridsExternal(char *pParams,
	                                              bool UserDefinedParams,
                                                  int JParam,
                                                  double RParam,
												  int MaxClassIndex,
                                                  char *TrainingDataPath,
                                                  char *MaskPath,
                                                  char *OutDir,
                                                  bool IgnoreEmptyClasses,
												  bool DoNormalisation,
												  char *TimeString,
                                                  FPTR fptr)
{
	//
	// open the training data table and extract the maximum class index
	//
	int NumEmptyClasses = 0;
	int *emptyClasses = GetEmptyClasses(&NumEmptyClasses, TrainingDataPath, MaxClassIndex);

	//
	// setup domain file paths
	//
	char pDomainPath[BUFFLEN];
	if (false == GetTransformGridAsDomain( pParams, pDomainPath ))
	{
		Message("Cannot extract transform domain path", "DoMaskedUnconstrainedBinaryProbGridsExternal");
		if (emptyClasses) delete[] emptyClasses;
		return(false);
	}

	gmpath gmPath;
	//
	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(pDomainPath, ".hdr"));
	//
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	double ResX = header->GetCellSize();
	float fNoData = header->GetNoDataValue();
	double dMinX = header->GetXllCorner();
	double dMinY = header->GetYllCorner();
	double dMaxX = dMinX + (ResX * nCols);
	double dMaxY = dMinY + (ResX * nRows);


	//
	// setup Mask file paths
	//
	EsriBinaryHeader *maskHeader = new EsriBinaryHeader(gmPath.ChangeExtension(MaskPath, ".hdr"));
	//
	int nMaskRows = maskHeader->GetNumRows();
	int nMaskCols = maskHeader->GetNumCols();
	double MaskResX = maskHeader->GetCellSize();
	float fMaskNoData = maskHeader->GetNoDataValue();
	double dMaskMinX = maskHeader->GetXllCorner();
	double dMaskMinY = maskHeader->GetYllCorner();
	double dMaskMaxX = dMaskMinX + (MaskResX * nMaskCols);
	double dMaskMaxY = dMaskMinY + (MaskResX * nMaskRows);


	// sanity check to make sure that the mask is a subset of the domain
	if ((dMaskMinX < (dMinX-MaskResX)) || (dMaskMaxX > (dMaxX+MaskResX)) || (dMaskMinY < (dMinY-MaskResX)) || (dMaskMaxY > (dMaxY+MaskResX)))
	{
		Message("The Mask Grid is NOT a subgrid of the Domain Grid", "DoMaskedUnconstrainedBinaryProbGridsExternal");
		char qqq[64];
		sprintf(qqq, "dMinX: %lf  dMaskMinX: %lf", dMinX, dMaskMinX);
		Message(qqq, "MinX");

		sprintf(qqq, "dMaxX: %lf  dMaskMaxX: %lf", dMaxX, dMaskMaxX);
		Message(qqq, "MaxX");

		sprintf(qqq, "dMinY: %lf  dMaskMinY: %lf", dMinY, dMaskMinY);
		Message(qqq, "MinY");

		sprintf(qqq, "dMaxY: %lf  dMaskMaxY: %lf", dMaxY, dMaskMaxY);
		Message(qqq, "MaxY");
		if (header) delete header;
		if (maskHeader) delete maskHeader;
		return(false);
	}


	// open the domain and mask grids
	BinaryFileClass *bfc_Domain = new BinaryFileClass(gmPath.ChangeExtension(pDomainPath, ".flt"), BFC_ReadOnly);	
	//
	BinaryFileClass *bfc_Mask = new BinaryFileClass(gmPath.ChangeExtension(MaskPath, ".flt"), BFC_ReadOnly);	

	//
	// create nClasses floating point output grids
	//
	float *pRowData = new float[ nMaskCols];
	for ( int i=0; i<nMaskCols; i++ )
	{
		pRowData[i] = fMaskNoData;   // init to no-data
	}

	int nCurrent = 0;
	fptr("Creating initialised probability grids...", nCurrent);	
	int nThis = 0;
	char myFilePath[BUFFLEN];
	for ( int i=1; i<=MaxClassIndex; i++ )
	{
		if ( ++nThis * 100 / MaxClassIndex > nCurrent )
		{
			nCurrent = nThis * 100 / MaxClassIndex;
			if ( nCurrent > 100 ) nCurrent = 100;
			fptr("Creating initialised probability grids...", nCurrent);	
		}

		//
		// ignore any classes that are not represented in the training data table
		//
		if (IgnoreEmptyClasses && EmptyClass(i, emptyClasses, NumEmptyClasses))
		{
			continue;
		}

		//
		// Create the .GRI grid
		//
		// create a filename for the grid
		sprintf( myFilePath, "%s\\uncon%03d.flt", OutDir, i );
		
		// create the grid
		BinaryFileClass *bfc_Tmp = new BinaryFileClass(myFilePath, BFC_CreateOrTrunc);

		// initialise to no-data
		for ( int j=0; j<nMaskRows; j++ )
		{
			bfc_Tmp->WriteFloat( pRowData, nMaskCols );
		}

		// close the grid
		bfc_Tmp->Close();
		if (bfc_Tmp) delete bfc_Tmp;

		
		//
		// Create the .GRD grid
		//
		// create a filename for the grid
		sprintf( myFilePath, "%s\\uncon%03d.hdr", OutDir, i );
		maskHeader->CopyTo(myFilePath);

		// sanity check for existence of the header file
		if ( ( _access( myFilePath, 0 ) ) == -1 ) 
		{
			char myErrorString[BUFFLEN];
			sprintf(myErrorString, "Cannot find %s", myFilePath );
			Message(myErrorString, "ERROR");
			return(false);
		}
	} // for ( i=1; i<=nMaxClass; i++ )

	if ( header) delete header;
	if ( maskHeader) delete maskHeader;


	//
	// get the bounding rectangle indices of the mask grid relative to the domain grid
	//
	int nMinRow = int((dMaxY - dMaskMaxY) / ResX);
	int nMaxRow = int((dMaxY - dMaskMinY) / ResX);
	int nMinCol = int((dMaskMinX - dMinX) / ResX);
	int nMaxCol = int((dMaskMaxX - dMinX) / ResX);

	//
	// Now setup the ANN tree from the training data
	//
	nCurrent = 0;
	fptr("Creating the ANN tree from the training data...", nCurrent);	
	char *scratch = new char [TABLE_ROW_BUFFSIZE];
	char seps[] = ",\n";

	FILE *fpEnv = fopen( TrainingDataPath, "r+t" );
	// assume that there is a HEADER ROW and the first column is the class,
	// the second column is the easting and the third column is the northing,
	// all the following columns correspond to the input grid file list.
	// assume that all the grid env values are floating point so the grids are also.

	// determine the number of columns
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ); // get the header	

	int nEnvCols = 0;
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv );
	char *p = strtok( scratch, seps ); // this will be the class column
	p = strtok( NULL, seps );          // this will be the easting column
	p = strtok( NULL, seps );          // this will be the northing column
	while( 1 ) 
	{
		if ( NULL == strtok( NULL, seps ) )
			break;

		else
			++nEnvCols;
	}

	// determine the number of rows
	int nEnvRows = 0;
	rewind( fpEnv );
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ); // get the header
	while( 1 )
	{
		if ( NULL == fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ) )
			break;

		else
			++nEnvRows;
	}
	

	///////////////////////////////////////////////////////////////////////////////////////
	// setup the approximate nearest neighbour data structs for the queries to follow
	ANNpoint query_pt = annAllocPt( nEnvCols );
	// this is for the whole data set
	ANNpointArray data_pnts = annAllocPts( nEnvRows, nEnvCols );

	int *pClasses = new int [ nEnvRows ];
	double *pTrainX = new double [ nEnvRows ];
	double *pTrainY = new double [ nEnvRows ];

	
	// now read the training table data into the point array
	rewind( fpEnv );
	// read the header
	fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv );


	// extract the class and the training data fields
	int nClasses = 0;
	int nThisRow = 0;
	while( 1 )
	{
		if ( NULL == fgets( scratch, TABLE_ROW_BUFFSIZE, fpEnv ) )
			break;

		// extract the class field 
		char *t = strtok( scratch, "," );
		pClasses[ nThisRow ] = atoi( t );

		t = strtok( NULL, "," );	// skip past the easting field....
		pTrainX[ nThisRow ] = atof( t );

		t = strtok( NULL, "," );	// skip past the northing field....
		pTrainY[ nThisRow ] = atof( t );

		// now extract the rest of the environmental data
		for ( int d=0; d<nEnvCols; d++ )
		{
			t = strtok( NULL, ",\n" );

			data_pnts[nThisRow][d] = atof( t );

		}

		++nThisRow;
			
	}
	fclose(fpEnv);


	// now build the kd tree search structure
	ANNkd_tree *theTree = new ANNkd_tree ( data_pnts,		// the data points
										   nEnvRows,		// number of data points
										   nEnvCols );		// dimension of search space

	//
	// Setup for the Unconstrained Predictions
	//
	double *explu = new double [JParam];
	// setup some structures for the queries to follow
	ANNidxArray nn_idx = new ANNidx[JParam];
	ANNdistArray dist = new ANNdist[JParam];

	// setup transform grid handles
	char myBuff[BUFFLEN];
	int NumPredGrids = GetPrivateProfileInt("PREDICTORS", "NumPredictors", 0, pParams);
	if ( NumPredGrids < 1 )
	{
		Message( "No Predictor Grids to use???", "DoMaskedUnconstrainedBinaryProbGridsExternal" );
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
			//printf( "%s does not exist\n", testBuff );
		}
		else
		{
			//printf( "TransPath: %s\n", myBuff );
			++numTranGrids;
		}
	}

	// check if there are geographic distance transform grids
	if (1 == GetPrivateProfileInt("GDMODEL", "UseEuclidean", 0, pParams))
	{
		GetPrivateProfileString("TRANSPREDS", "EuclXTran", "0", myBuff, BUFFLEN, pParams);
		if (0 == strcmp("0",myBuff))
		{
			Message( "EuclXTran does not exist", "ERROR" );
			return(false);
		}
		else
		{
			++numTranGrids;
		}

		GetPrivateProfileString("TRANSPREDS", "EuclYTran", "0", myBuff, BUFFLEN, pParams);
		if (0 == strcmp("0",myBuff))
		{
			Message( "EuclYTran does not exist", "ERROR" );
			return(false);
		}
		else
		{
			++numTranGrids;
		}
	}

	if ( numTranGrids < 1 )
	{
		Message( "No Transform Grids to use???", "DoMaskedUnconstrainedBinaryProbGridsExternal" );
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
	if (1 == GetPrivateProfileInt("GDMODEL", "UseEuclidean", 0, pParams))
	{
		nThis = 0;
		GetPrivateProfileString("TRANSPREDS", "EuclXTran", "0", myBuff, BUFFLEN, pParams);
		if ( 0 != strcmp("0",myBuff) )
		{
			strcpy (ppTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
		}
		GetPrivateProfileString("TRANSPREDS", "EuclYTran", "0", myBuff, BUFFLEN, pParams);
		if ( 0 != strcmp("0",myBuff) )
		{
			strcpy (ppTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
		}
		for ( int i=1; i<=NumPredGrids; i++ )
		{
			sprintf( testBuff, "PredTran%d", i );
			GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, pParams);
			if ( 0 != strcmp("0",myBuff) )
			{
				strcpy (ppTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
			}
		}
	}
	else
	{
		nThis = 0;
		for ( int i=1; i<=NumPredGrids; i++ )
		{
			sprintf( testBuff, "PredTran%d", i );
			GetPrivateProfileString("TRANSPREDS", testBuff, "0", myBuff, BUFFLEN, pParams);
			if ( 0 != strcmp("0",myBuff) )
			{
				strcpy (ppTranGridPaths[nThis++], gmPath.ChangeExtension(myBuff, ".flt") );
			}
		}
	}



	//
	// open the transform grids
	//
	BinaryFileClass **bfcTran = new BinaryFileClass * [numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		bfcTran[i] = new BinaryFileClass(ppTranGridPaths[i]);
		if (!bfcTran[i]->IsValid())
		{
			Message( "Cannot open Transform Grid", "DoMaskedUnconstrainedBinaryProbGridsExternal" );
			return(false);
		}
	}



	//
	// allocate a set of row vectors
	//
	float **ppRowVectors = new float *[numTranGrids];
	for ( int i=0; i<numTranGrids; i++ )
	{
		ppRowVectors[i] = new float [nCols];
	}



	// allocate a vector to the masked output grid rows( one for each class )
	float **ppOutVec = new float * [MaxClassIndex];
	for ( int i=0; i<MaxClassIndex; i++ )
	{
		ppOutVec[i] = new float [nMaskCols];
		for ( int q=0; q<nMaskCols; q++ ) ppOutVec[i][q] = fMaskNoData;
	}



	float *pNoDataVec = new float [nMaskCols];
	for ( int x=0; x<nMaskCols; x++ )
	{
		pNoDataVec[x] = fMaskNoData;
	}


	float *pMaxProbs = new float [MaxClassIndex];
	BinaryFileClass **bfcOut = new BinaryFileClass * [MaxClassIndex];
	for ( int i=1; i<=MaxClassIndex; i++ )
	{
		//
		// Init probs to zero
		//
		pMaxProbs[i-1] = 0.0;
		
		//
		// ignore any classes that are not represented in the training data table
		//
		if (IgnoreEmptyClasses && EmptyClass(i, emptyClasses, NumEmptyClasses))
		{
			continue;
		}

		// open the out grid
		sprintf( myFilePath, "%s\\uncon%03d.flt", OutDir, i );
		bfcOut[i-1] = new BinaryFileClass(myFilePath, BFC_ReadWrite);
		if (!bfcOut[i-1]->IsValid())
		{
			Message("Cannot open outgrid", "DoMaskedUnconstrainedBinaryProbGridsExternal");
			return(false);
		}
	} // for ( i=1; i<=nMaxClass; i++ )


	//
	// Create a metadata file detailing the creation parameters and the output grid metrics
	//
	sprintf(myFilePath, "%s\\metadata.txt", OutDir);
	FILE *fpMeta = fopen(myFilePath, "w+t");
	fprintf(fpMeta, "/////////////////////////////////////////////////////////////////////////////////////////////////\n");
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "Metadata for unconstrained probability grid creation: %s\n", TimeString);
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "/////////////////////////////////////////////////////////////////////////////////////////////////\n");
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "InputTable: %s\n", TrainingDataPath);
	fprintf(fpMeta, "Mask: %s\n", MaskPath);
	fprintf(fpMeta, "Output Directory: %s\n", OutDir);
	fprintf(fpMeta, "Maximum Class Index: %d\n", MaxClassIndex);
	fprintf(fpMeta, "J Param: %d       R Param: %lf\n", JParam, RParam);
	if (IgnoreEmptyClasses)
		fprintf(fpMeta, "Empty Classes are excluded\n");
	else
		fprintf(fpMeta, "Empty Classes are included\n");
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "/////////////////////////////////////////////////////////////////////////////////////////////////\n");
	

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Finally do the unconstrained predictions
	//
	nCurrent = 0;
	fptr("Calculating probability grids...", nCurrent);	
	for ( int y=0; y<nRows; y++ )
	{
		// update status
		if ( y * 100 / nRows > nCurrent )
		{
			nCurrent = y * 100 / nRows;
			fptr("Calculating probability grids...", nCurrent);	
		}


		// skip any rows that are out of the Mask extent
		if (( y < nMinRow ) || (y >= nMaxRow) ) continue;


		// read a row from the mask grid
		bfc_Mask->SeekTo((y-nMinRow) * nMaskCols * sizeof(float));
		bfc_Mask->ReadFloat(pRowData, nMaskCols);


		// read in the row data from the transform grids
		for ( int i=0; i<numTranGrids; i++ )
		{
			bfcTran[i]->SeekTo(y * nCols * sizeof(float));
			bfcTran[i]->ReadFloat(ppRowVectors[i], nCols);
		}


		// initialise the output rows to all No_Data
		for ( int nThisClass=0; nThisClass<MaxClassIndex; nThisClass++ )
		{
			// initialise the output rows' probabilities to No-Data
			memcpy( ppOutVec[nThisClass], pNoDataVec, nMaskCols * sizeof(float) ); 
		}


		for ( int x=0; x<nCols; x++ )
		{
			// only if we have data and the cell is inside the mask grid
			if ((x < nMinCol) || (x >= nMaxCol)) continue;

			//if ( ppRowVectors[0][x] != fNoData )
			if (pRowData[x-nMinCol] != fMaskNoData)
			{
				/////////////////////////////////////////////////////
				/////////////////////// STEP 1 //////////////////////
				// get the nJParam nearest neighbours to this site
				for ( int c=0; c<nEnvCols; c++ ) 
				{
					query_pt[ c ] = ppRowVectors[c][x];
				}

				theTree->annkSearch( query_pt, JParam, nn_idx, dist, 0 );

				/////////////////////////////////////////////////////
				/////////////////////// STEP 2 //////////////////////
			
				// CALCULATE SOME PROBABILITIES HERE FOR THIS POINT
				/////////////////////////////////////////////////////
				// section A
				double sumtheta = 0.0;

				for ( int m=0; m<JParam; m++ )
				{
					//dist[m] = sqrt(dist[m]);		// taking the square root of the ANN distance 
					sumtheta += dist[m];
				}

				double theta = sumtheta / (double)JParam * RParam;
				// square above value and multiply by 2 ready for section B;
				theta = 2 * theta * theta;	
				/////////////////////////////////////////////////////

				/////////////////////////////////////////////////////
				// section B
				double sumJValues = 0.0;
				
				int j;
				for ( j=0; j<JParam; j++ )
				{
					if ( theta > 0.0 )
					{
						explu[j] = exp((-(dist[j] * dist[j]))/theta);
						sumJValues += explu[j];
					}
				}

				/////////////////////////////////////////////////////
				// now apply the conditional probabilities to the row vectors
				for ( int nThisClass=0; nThisClass<MaxClassIndex; nThisClass++ )
				{
					// initialise the output rows' probabilities to zero
					ppOutVec[nThisClass][x-nMinCol] = 0L; 
				}

				for ( j=0; j<JParam; j++ )
				{					
					// copy the probabilities over to the out vector
					if ( sumJValues > 0.0 )
						ppOutVec[(pClasses[nn_idx[j]])-1][x-nMinCol] += (float)(explu[j] / sumJValues);
				}

			} // if ( ppRowVectors[0][x] != fNoData )

		} // for ( int x=0; x<nCols; x++ )


		//
		// adjust probabilities so that they add to 1.0 for each valid data cell
		// (or 0.0 if there are no non-zero probabilities)
		//
		if (DoNormalisation)
		{
			for ( int x=0; x<nMaskCols; x++ )
			{
				if (ppOutVec[0][x] != fNoData) // a valid data cell
				{
					// sum probabilities
					float fProbSum = 0.0F;
					for ( int i=0; i<MaxClassIndex; i++ )
					{
						fProbSum += ppOutVec[i][x];
					}

					// normalise
					if (fProbSum > 0.0F)
					{
						for ( int i=0; i<MaxClassIndex; i++ )
						{
							ppOutVec[i][x] /= fProbSum;
						}
					}

				} // if (ppOutVec[0][x] != fNoData)
			} // for ( int x=0; x<nCols; x++ )
		} // if (DoNormalisation)


		// now write each row in the sum array to the output grid
		for ( int i=1; i<=MaxClassIndex; i++ )
		{
			//
			// ignore any classes that are not represented in the training data table
			//
			if (IgnoreEmptyClasses && EmptyClass(i, emptyClasses, NumEmptyClasses))
			{
				continue;
			}

			//
			// update the output grids
			//
			bfcOut[i-1]->SeekTo((y-nMinRow) * nMaskCols * sizeof(float));
			bfcOut[i-1]->WriteFloat(ppOutVec[i-1], nMaskCols);
		} // for ( i=1; i<=nMaxClass; i++ )


		//
		// update the min and max probibilities
		//
		for ( int i=0; i<MaxClassIndex; i++ )
		{
			for (int j=0; j<nMaskCols; j++)
			{
				if ( ppOutVec[i][j] != fNoData )
				{
					if (pMaxProbs[i] < ppOutVec[i][j]) pMaxProbs[i] = ppOutVec[i][j];
				}
			}
		}

	} // for ( int y=0; y<nRows; y++ )


	//
	// Finish writing the metadata 
	//
	for ( int i=0; i<MaxClassIndex; i++ )
	{
		fprintf(fpMeta, "Class: %d   MaxProb: %f\n", i+1, pMaxProbs[i]);
	}
	fprintf(fpMeta, "//\n");
	fprintf(fpMeta, "/////////////////////////////////////////////////////////////////////////////////////////////////\n");
	if (fpMeta) fclose(fpMeta);


	//
	// clean up
	//	
	bfc_Domain->Close();
	if (bfc_Domain) delete bfc_Domain;
	bfc_Mask->Close();
	if (bfc_Mask) delete bfc_Mask;	
	if (pMaxProbs) delete[] pMaxProbs;
	if ( pRowData ) delete[] pRowData;

	for ( int i=1; i<=MaxClassIndex; i++ )
	{
		//
		// ignore any classes that are not represented in the training data table
		//
		if (IgnoreEmptyClasses && EmptyClass(i, emptyClasses, NumEmptyClasses))
		{
			continue;
		}

		// close the grid
		bfcOut[i-1]->Close();
		if (bfcOut[i-1]) delete bfcOut[i-1];
	} // for ( i=0; i<0nMaxClass; i++ )

	if (bfcOut) delete[] bfcOut;
	if (emptyClasses) delete[] emptyClasses;
	if (pClasses) delete[] pClasses;
	if (pTrainX) delete[] pTrainX;
	if (pTrainY) delete[] pTrainY;
	if ( scratch ) delete[] scratch;
	if ( explu ) delete[] explu;
	if ( pNoDataVec ) delete[] pNoDataVec;
	if (theTree) delete theTree;

	for ( int i=0; i<numTranGrids; i++ )
	{
		bfcTran[i]->Close();
		if (bfcTran[i]) delete bfcTran[i];
	}
	if (bfcTran) delete[] bfcTran;
	for ( int i=0; i<numTranGrids; i++ )
	{
		if (ppTranGridPaths[i]) delete[] ppTranGridPaths[i];
	}
	if (ppTranGridPaths) delete[] ppTranGridPaths;
	for ( int i=0; i<numTranGrids; i++ )
	{
		if ( ppRowVectors[i] ) delete[] ppRowVectors[i];
	}
	if ( ppRowVectors ) delete[] ppRowVectors;
	for ( int i=0; i<MaxClassIndex; i++ )
	{
		if (ppOutVec[i]) delete[] ppOutVec[i];
	}
	if ( ppOutVec ) delete[] ppOutVec;
	return(true);
}



//
// Get the maximum one-based class index from the training data table
//
bool GetMaxClassIndexFromTrainingData(int *pMaxClassIndex, char *JRInputPath, FPTR fptr)
{
	//
	// Get the training data table metrics
	//
	int nRows = 0;
	int nCols = 0;
	if (!GetANNTableMetrics(JRInputPath, &nRows, &nCols))
	{
		return(false);
	}


	//
	// Sanity check
	//
	if (nRows < 3)
	{
		Message("There MUST be at least 2 records in the training data table", "GetBestJRParams");
		return(false);
	}

	//
	// Get the maximum class index
	//
	int *pClasses = GetClassesFromANN(JRInputPath, nRows, pMaxClassIndex, fptr);
	if (pClasses) delete[] pClasses;
	return(true);
}


//
// Wrapper function called from the Kernel Regression Form in GDM
//
bool GetBestJRValsViaParamFile(char *ParamFilePath, FPTR fptr)
{
	int nJParam = 0;
	double dRParam = 0.0;
	int MaxClass = 0;
	char lpTrainingDataPath[BUFF1024];
	GetPrivateProfileString("INPUTS", "TrainingData", "", lpTrainingDataPath, BUFF1024, ParamFilePath);
	bool retVal = GetBestJRParams(&nJParam, &dRParam, &MaxClass, lpTrainingDataPath, NULL, fptr);
	if (retVal == true)
	{
		// update parameter file
		SetProfileInt("INPUTS", "JParameter", nJParam, ParamFilePath);
		SetProfileDouble("INPUTS", "RParameter", dRParam, ParamFilePath);
		SetProfileInt("INPUTS", "MaxClass", MaxClass, ParamFilePath);
	}
	return(retVal);
}


//
// Master Test Function for SOS calculations
//
bool GetBestJRParams(int *pJParam, double *pRParam, int *pMaxClassIndex, char *JRInputPath, char *JROutputPath, FPTR fptr)
{
	//
	// Get the training data table metrics
	//
	int nRows = 0;
	int nCols = 0;
	if (!GetANNTableMetrics(JRInputPath, &nRows, &nCols))
	{
		return(false);
	}

	//
	// Sanity check
	//
	if (nRows < 3)
	{
		Message("There MUST be at least 2 records in the training data table", "GetBestJRParams");
		return(false);
	}

	//
	// Allocate ANN support vectors
	//
	ANNidxArray nn_idx = new ANNidx[nRows];
	ANNdistArray dist = new ANNdist[nRows];	
	ANNpoint query_pnt = annAllocPt(nCols);
	ANNpointArray data_pnts = GetDataAsANNTable(JRInputPath, nRows, nCols, fptr);
	double *explu = new double [nRows];
	double *pOutVec = new double [nRows];
	int MaxClassIndex = 0;
	int *pClasses = GetClassesFromANN(JRInputPath, nRows, &MaxClassIndex, fptr);

	// build the kd tree search structure
	ANNkd_tree *theTree = new ANNkd_tree ( data_pnts,		// the data points
										   nRows,		    // number of data points
								           nCols,		    // dimension of search space
										   20,				// bucket size
										   ANN_KD_SUGGEST );// splitting method

	//
	// Create output file
	//
	FILE *fpOut = NULL;
	if (NULL != JROutputPath)
	{
		fpOut = fopen(JROutputPath, "w+t");
		if (NULL == fpOut)
		{
			Message("Cannot create open JROutputPath", "ERROR");
			if (nn_idx) delete[] nn_idx;
			if (dist) delete[] dist;
			if (explu) delete[] explu;
			if (pClasses) delete[] pClasses;
			if (pOutVec) delete[] pOutVec;
			annDeallocPt(query_pnt);
			annDeallocPts(data_pnts);
			if (theTree) delete theTree;
			return(false);
		}
		fprintf(fpOut, "J_Param,R_Param,SOS\n");
	}

	//
	// Adjust for Max J Parameter
	//
	int MinJVal = 30;
	int MaxJVal = 100;
	if (nRows < MaxJVal) 
		MaxJVal = nRows-1;

	//
	// Run the repeated sum of squares calculations
	//
	char buff[64];
	//double dMinSOS = 1000000000.0;
	double dMaxSOS = 0.0;
	int nBestJ = -1;
	double dBestR = -1.0;
	for (double dR = 0.1; dR < 1.0; dR += 0.1)
	{
		sprintf(buff, "Extracting Optimal Kernel Parameters for R=%0.1lf", dR);
		fptr(buff, 0);

		double dLocalMax = 0.0;
		int nLocalJ = -1;
		int nCurrent = 0;
		for (int nJ = MinJVal; nJ <= MaxJVal; nJ++)
		{
			if (nCurrent < nJ * 100 / MaxJVal)
			{
				nCurrent = nJ * 100 / MaxJVal;
				fptr(buff, nCurrent);
			}

			double dSOS = GetSumOfSquaresANN(
				theTree,
				nRows, nCols,
				nJ, dR,
				nn_idx,
				dist,
				explu,
				query_pnt,
				data_pnts,
				MaxClassIndex,
				pClasses,
				pOutVec,
				fptr);

			if (NULL != JROutputPath)
			{
				if (fpOut)
					fprintf(fpOut, "%d,%lf,%lf\n", nJ, dR, dSOS);
			}

			// Update
			if (dSOS > dMaxSOS)
			{
				dMaxSOS = dSOS;
				nBestJ = nJ;
				dBestR = dR;
			}
			if (dSOS > dLocalMax) 
			{
				dLocalMax = dSOS;
				nLocalJ = nJ;
			}
		} // for (int nJ = MinJVal; nJ <= MaxJVal; nJ++)
		fptr(buff, 100);

		//
		// we may have reached the inflexion point for this R 
		// so we can terminate as we have found our optimal J and R...
		//
		if ((nLocalJ > MinJVal) && (nLocalJ < MaxJVal))
		{
			if (fpOut) fclose(fpOut);
			//
			// set J and R params and the max class index
			//
			*pJParam = nLocalJ;
			*pRParam = dR;
			*pMaxClassIndex = MaxClassIndex;

			//
			// Clean up
			//
			fptr("Ready: ", 0);
			if (nn_idx) delete[] nn_idx;
			if (dist) delete[] dist;
			if (explu) delete[] explu;
			if (pClasses) delete[] pClasses;
			if (pOutVec) delete[] pOutVec;
			annDeallocPt(query_pnt);
			annDeallocPts(data_pnts);
			if (theTree) delete theTree;
			return(true);
		} // if ((nLocalJ > MinJVal) && (nLocalJ < MaxJVal)) terminating condition is met
	} // for (double dR = 0.1; dR < 1.0; dR += 0.1)

	//
	// set J and R params and the max class index
	//
	*pJParam = nBestJ;
	*pRParam = dBestR;
	*pMaxClassIndex = MaxClassIndex;

	//
	// Clean up
	//
	if (NULL != JROutputPath)
	{
		if (fpOut) fclose(fpOut);
	}
	fptr("Ready: ", 0);
	if (nn_idx) delete[] nn_idx;
	if (dist) delete[] dist;			
	if (explu) delete[] explu;
	if (pClasses) delete[] pClasses;
	if (pOutVec) delete[] pOutVec;
	annDeallocPt(query_pnt);
	annDeallocPts(data_pnts);
	if (theTree) delete theTree;
	return(true);
}



//
// Calculate the Sum Of Squares for an ANN dataset given 
// a particular number of nearest neighbours (nJ) to find
//
double GetSumOfSquaresANN( ANNkd_tree *theTree, 
	                       int nDataRows, int nDataCols, 
						   int nJ, double dR, 
						   ANNidxArray nn_idx,
						   ANNdistArray dist,
						   double *explu,
						   ANNpoint query_pnt,
						   ANNpointArray data_pnts,
						   int NumClasses, int *pClasses, double *pOutVec,
						   FPTR fptr)
{
	double dSOS = 0.0;
	for (int i=0; i<nDataRows; i++ )
	{
		//
		// Initialise the query vector with the current data row
		//
		for (int j=0; j<nDataCols; j++)
		{
			query_pnt[j] = data_pnts[i][j];
		}

		//
		// sort ANN tree
		//
		theTree->annkSearch(query_pnt, nJ, nn_idx, dist, 0.0);

		//
		// Calculate Probabilities
		//
		double dSumTheta = 0.0;
		for (int j=0; j<nJ; j++)
		{
			dSumTheta += dist[j];
		}
		double dTmpTheta = dR * dSumTheta / nJ;

		//
		// Sum Probabilities
		//
		double theta = 2.0 * dTmpTheta * dTmpTheta;
		double dSumJValues = 0.0;
		for (int j=0; j<nJ; j++)
		{
			if (theta > 0.0)
			{
				explu[j] = exp((-(dist[j] * dist[j]))/theta);
				dSumJValues += explu[j];
			}
			else
			{
				explu[j] = 0.0;
			}
		}

		//
		// Initialise Probabilities To Zero
		//
		for (int j=0; j<NumClasses; j++)
		{
			pOutVec[j] = 0.0;
		}

		//
		// Apply Conditional Probabilities To Each Class
		//
		for (int j=0; j<nJ; j++)
		{
			if (dSumJValues > 0.0)
			{
				pOutVec[(pClasses[nn_idx[j]])-1] += explu[j] / dSumJValues;
			}
		}

		//
		// Finally, Accumulate The Sum Of Squares
		//
		for (int j = 0; j < NumClasses; j++)
		{
			// if this class matches the training class
			double dProb;
			if ((pClasses[i]-1) == j)
			{
				dProb = (pOutVec[j] - 1.0) *(pOutVec[j] - 1.0);
			}

			else
			{
				dProb = pOutVec[j] *pOutVec[j];
			}
			dSOS += dProb;
		}
	} // for (int i=0; i<nDataRows; i++ )

	//
	// Divide by nJ to get the average Sum Of Squares
	//
	dSOS /= (double)nJ;
	return(dSOS);
}




//
// Determine the number of data rows from a comma delimited text file with header
// and determine the number of predictor columns after skipping Class,X and Y.
//
bool GetANNTableMetrics(char *TablePath, int *pRows, int *pCols)
{
	FILE *fp = fopen(TablePath, "r+t");
	if (NULL == fp)
	{
		Message("Unable to open Table Path", "GetANNTableMetrics");
		return(false);
	}

	//
	// extract the number of cols from the header
	//
	char *pBuff = new char [TABLE_ROW_BUFFSIZE];
	char seps[] = ",\n";
	int nCols = 0;
	fgets(pBuff, TABLE_ROW_BUFFSIZE, fp); // get header
	char *p = strtok(pBuff, seps);        // get class
	p = strtok(NULL, seps);               // get X
	p = strtok(NULL, seps);               // get Y	
	while(1)
	{
		if (NULL == strtok(NULL, seps))
			break;
		else 
			++nCols;
	}

	//
	// now determine the number of data rows
	//
	int nRows = 0;
	while(1)
	{
		if (NULL == fgets(pBuff, TABLE_ROW_BUFFSIZE, fp))
			break;
		else
			++nRows;
	}

	if(pBuff) delete[] pBuff;
	if (fp) fclose(fp);
	*pRows = nRows;
	*pCols = nCols;
	return(true);
}



//
// Extract the predictor columns from a comma delimited training data table 
// after skipping the header and after skipping the first three columns (class,X,Y)
//
ANNpointArray GetDataAsANNTable(char *TablePath, int nRows, int nCols, FPTR fptr)
{
	FILE *fp = fopen(TablePath, "r+t");
	if (NULL == fp)
	{
		Message("Unable to open Table Path", "GetDataAsANNTable");
		return(NULL);
	}

	ANNpointArray pData = annAllocPts(nRows,nCols);
	if (NULL == pData)
	{
		Message("Unable to allocate ANNpointArray", "GetDataAsANNTable");
		return(NULL);
	}

	int nCurrent = 0;
	fptr("GetDataAsANNTable",nCurrent);
	char *pBuff = new char [TABLE_ROW_BUFFSIZE];
	char seps[] = ",\n";
	fgets(pBuff, TABLE_ROW_BUFFSIZE, fp); // get header
	for (int i=0; i<nRows; i++)
	{
		if (i * 100 / nRows > nCurrent)
		{
			nCurrent = i * 100 / nRows;
			fptr("GetDataAsANNTable",nCurrent);
		}

		fgets(pBuff, TABLE_ROW_BUFFSIZE, fp); // read a row
		char *p = strtok(pBuff, seps);        // skip class field
		p = strtok(NULL, seps);               // skip X field
		p = strtok(NULL, seps);               // skip Y field

		for ( int j=0; j<nCols; j++)
		{
			p = strtok(NULL, seps);      // current predictor field
			pData[i][j] = atof(p);
		}
	}
	if(pBuff) delete[] pBuff;
	if (fp) fclose(fp);
	return(pData);
}



//
// Extract the class field data [0] from a comma delimited training data table 
// after skipping the header and also determining the Maximum Class Index
//
int *GetClassesFromANN(char *TablePath, int nRows, int *MaxClassIndex, FPTR fptr)
{
	FILE *fp = fopen(TablePath, "r+t");
	if (NULL == fp)
	{
		Message("Unable to open Table Path", "GetDataAsANNTable");
		return(NULL);
	}

	int nCurrent = 0;
	fptr("GetClassesFromANN",nCurrent);

	int *pClasses = new int [nRows];
	char *pBuff = new char [TABLE_ROW_BUFFSIZE];
	char seps[] = ",\n";
	fgets(pBuff, TABLE_ROW_BUFFSIZE, fp); // get header

	int MaxClass = 0;
	for (int i=0; i<nRows; i++)
	{
		if (i * 100 / nRows > nCurrent)
		{
			nCurrent = i * 100 / nRows;
			fptr("GetClassesFromANN",nCurrent);
		}

		fgets(pBuff, TABLE_ROW_BUFFSIZE, fp);

		char *p = strtok(pBuff, seps);   // select class field
		pClasses[i] = atoi(p);

		// update the Max Class Index
		if (pClasses[i] > MaxClass) MaxClass = pClasses[i];
	}
	*MaxClassIndex = MaxClass;
	if(pBuff) delete[] pBuff;
	if (fp) fclose(fp);
	return(pClasses);
}


//
// Extract a vector of class indices representing classes to skip in the probability grid creation
//
int *GetEmptyClasses(int *pNumEmptyClasses, char *pTrainingFile, int nMaxClass)
{
	//
	// allocate a presence/absence vector
	//
	bool *pPA = new bool [nMaxClass];
	for ( int i=0; i<nMaxClass; i++ ) pPA[i] = false;


	//
	// determine Presence/absence
	//
	FILE *fp = fopen(pTrainingFile, "r+t");
	char *rowbuff = new char [TABLE_ROW_BUFFSIZE];
	fgets( rowbuff, TABLE_ROW_BUFFSIZE, fp ); // get header
	while(1)
	{
		if ( NULL == fgets( rowbuff, TABLE_ROW_BUFFSIZE, fp ) )
			break;

		char *p = strtok(rowbuff, ",");
		int nVal = atoi(p);

		// convert one-based class index to zero based array index...
		pPA[nVal-1] = true;
	}
	if ( fp ) fclose(fp);
	if ( rowbuff ) delete[] rowbuff;


	//
	// count the number of empty classes
	//
	int tmpNumEmptyClasses = 0;
	int i;
	for ( i=0; i<nMaxClass; i++ )
	{
		if (!pPA[i]) ++tmpNumEmptyClasses;
	}

	//
	// there may be NO empty classes so just return NULL
	//
	if ( 0 == tmpNumEmptyClasses )
	{
		*pNumEmptyClasses = 0;
		if (pPA) delete[] pPA;
		return(NULL);
	}

	else
	{
		//
		// allocate a vector to hold the empty class indices
		//
		int *EmptyClassVector = new int [tmpNumEmptyClasses];

		//
		// set the empty class indices
		//
		int nThis = 0;
		for ( int i=0; i<nMaxClass; i++ )
		{
			if (!pPA[i])
			{
				EmptyClassVector[nThis] = i+1;  // one-based class indices
				++nThis;
			}
		}
		*pNumEmptyClasses = tmpNumEmptyClasses;
		if (pPA) delete[] pPA;
		return(EmptyClassVector);
	}
}



//
// Returns true if veg class in empty in the training data table
//
bool EmptyClass(int Class, int *pClassVector, int nSize)
{
	// we may have an empty vector (no empty classes) so just return false
	if ( nSize < 1 )
		return(false);

	for ( int i=0; i<nSize; i++ )
	{
		if (Class == pClassVector[i])
			return(true);
	}
	return(false);
}


//
// Returns true if we can extract a valid transform grid path from a GDM parameter file
//
bool GetTransformGridAsDomain( char *pParams, char *pDomainPath )
{
	//
	// check if we have eucidean distance transforms
	//
	if ( 1 == GetProfileInt( "GDMODEL", "UseEuclidean", pParams ) )
	{
		GetPrivateProfileString("TRANSPREDS", "EuclXTran", "NOVALUE", pDomainPath, BUFFLEN, pParams);
		if (0 != strcmp("NOVALUE", pDomainPath))
		{
			return(true);
		}

		GetPrivateProfileString("TRANSPREDS", "EuclYTran", "NOVALUE", pDomainPath, BUFFLEN, pParams);
		if (0 != strcmp("NOVALUE", pDomainPath))
		{
			return(true);
		}
	}

	//
	// try the environmental predictors
	// 
	int nPreds = GetProfileInt( "PREDICTORS", "NumPredictors", pParams );
	char lpKey[64];
	for ( int i=1; i<= nPreds; i++ )
	{
		sprintf(lpKey, "PredTran%d", i);

		GetPrivateProfileString("TRANSPREDS", lpKey, "NOVALUE", pDomainPath, BUFFLEN, pParams);

		if (0 != strcmp("NOVALUE", pDomainPath))
		{
			return(true);
		}
	}

	//
	// if we get this far then we have a problem
	//
	return(false);
}

