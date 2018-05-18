//
// GdmTransform.cpp
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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
#define GDM_Extrapolation_Default10			0
#define GDM_Extrapolation_WholeGradient		1
#define GDM_Extrapolation_MostConservative	2
#define GDM_Extrapolation_None  			3
//
// The above extrapolation modes are set in the parameter file as [] ExtrapolationMode=0,1,2,3
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// create GDM Parameter File from R model data and user defined predictor file list and perform GDM transformations on the In Grid List
//
bool CreateTransformsFromR(char **ppInputGridTable,
	                       char **ppOutputGridsWorkDir,
	                       int *Extrap_Method,
	                       bool *DoGeo, 
	                       int *NumPredictors,
	                       int *pSplines,
	                       double *pQuantiles, 
	                       double *pCoeffs)
{
	int nPreds = *NumPredictors;
	int nDoGeo = *DoGeo;

	char buff0[BUFF1024];
	std::sprintf(buff0, "%s\\GDM_Params.txt", *ppOutputGridsWorkDir);
	FILE *fpOut = fopen(buff0, "w+t");
	std::fprintf(fpOut, "###\n### Parameter file create from R to transform grids via the GDM dll\n###\n\n");

	std::fprintf(fpOut, "[GDMODEL]\n");
	std::fprintf(fpOut, "UseEuclidean=%d\n", nDoGeo);
	std::fprintf(fpOut, "ExtrapolationMethod=%d\n", *Extrap_Method);
	std::fprintf(fpOut, "WorkspacePath=%s\n", *ppOutputGridsWorkDir);
	std::fprintf(fpOut, "\n\n");
		
	std::fprintf(fpOut, "[PREDICTORS]\n");
	std::fprintf(fpOut, "NumPredictors=%d\n", nPreds);
	std::fprintf(fpOut, "\n\n");
		
	int nGeoOffset = 0;
	if (nDoGeo)
	{
		nGeoOffset = pSplines[0] - 1;

		std::fprintf(fpOut, "EuclSpl=%d\n", pSplines[0]);  // number of splines

		for (int j = 0; j < pSplines[0]; j++)
		{
			std::fprintf(fpOut, "EuclSplVal%d=%lf\n", j + 1, pQuantiles[j]);  // quantile
		}

		for (int j = 0; j < pSplines[0]; j++)
		{
			std::fprintf(fpOut, "EuclCoef%d=%lf\n", j + 1, pCoeffs[j]);  // coefficient
		}

		std::fprintf(fpOut, "\n\n");
	}
	
	char buff1[BUFF1024];
	FILE *fpIn = fopen(*ppInputGridTable, "r+t");
	for (int i = 1; i <= nPreds; i++)
	{
		std::fgets(buff1, BUFF1024, fpIn);
		std::fprintf(fpOut, "EnvGrid%d=%s", i, buff1);
	}
	std::fprintf(fpOut, "\n\n");
	fclose(fpIn);

	for (int i = 1; i <= nPreds; i++)
	{
		std::fprintf(fpOut, "PredType%d=%d\n", i, 0);  // PredType = GRID
	}
	std::fprintf(fpOut, "\n\n");

	for (int i = nGeoOffset; i < nPreds; i++)
	{
		std::fprintf(fpOut, "PredSpl%d=%d\n", i+1, pSplines[i]);  // number of splines
	}
	std::fprintf(fpOut, "\n\n");

	for (int i = 0; i < nPreds; i++)
	{
		for (int j = 0; j < pSplines[i]; j++)
		{
			std::fprintf(fpOut, "PredSplVal%d.%d=%lf\n", i+1, j+1, pQuantiles[((nGeoOffset+i)*pSplines[i]) + j]);  // quantile
		}
		for (int j = 0; j < pSplines[i]; j++)
		{
			std::fprintf(fpOut, "PredCoef%d.%d=%lf\n", i+1, j+1, pCoeffs[((nGeoOffset + i)*pSplines[i]) + j]);  // coefficient
		}
		std::fprintf(fpOut, "\n\n");
	}
	std::fprintf(fpOut, "\n\n");

	std::fprintf(fpOut, "[TRANSPREDS]\n\n");
	if (nDoGeo)
	{
		std::fprintf(fpOut, "%s\\gdmXtran\n", *ppOutputGridsWorkDir);
		std::fprintf(fpOut, "%s\\gdmYtran\n", *ppOutputGridsWorkDir);
	}

	gmpath gmp;
	char buff2[BUFF1024];
	fpIn = fopen(*ppInputGridTable, "r+t");
	for (int i = 0; i < nPreds; i++)
	{
		std::fgets(buff2, BUFF1024, fpIn);
		double CoeffSum = 0.0;

		for (int j = 0; j < pSplines[i]; j++)
		{
			CoeffSum += pCoeffs[(i*pSplines[i]) + j];
		}

		if (CoeffSum != 0.0)  // we have a valid GDM Transform File
		{
			std::fprintf(fpOut, "PredTran%d=%s\\%s", i + 1, *ppOutputGridsWorkDir, gmp.GetName(buff2));
		}
	}
	fclose(fpIn);
	fclose(fpOut);

	
	//
	// Do the Transform
	//
	std::sprintf(buff0, "%s\\GDM_Params.txt", *ppOutputGridsWorkDir);	
	//Message(buff0, "buff0");
	TransformPredictors(buff0, false, NULL);
	
	return(true);
}
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//#define DO_PROGRESS 
#undef DO_PROGRESS   // turn all the FPTR stuff OFF to run from R


//
// Extract min and max values from a GDM transform grid
//
bool GetTransformMinMax( char *sPath, double *pMin, double *pMax, FPTR fptr)
{
	gmpath gmPath;

	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(sPath, ".hdr"));
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	float fNoData = header->GetNoDataValue();

	BinaryFileClass *data = new BinaryFileClass(gmPath.ChangeExtension(sPath, ".flt"));
	if (!data->IsValid())
	{
		Message( "Cannot open Transform Grid", "GetTransformMinMax" );
		return(false);
	}

	float *pRow = new float [nCols];

	float dTmpMin = 1000000.0F;
	float dTmpMax = -1000000.0F;
	int nCurrent = 0;
#ifdef	DO_PROGRESS
	if (fptr) fptr(" ", nCurrent);
#endif
	for ( int i=0; i<nRows; i++ )
	{
#ifdef	DO_PROGRESS
		if (i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			if (fptr) fptr(" ", nCurrent);
		}
#endif

		data->ReadFloat(pRow, nCols);

		for ( int j=0; j<nCols; j++ )
		{
			if ( pRow[j] != fNoData )
			{
				if (pRow[j] < dTmpMin) dTmpMin = pRow[j];
				if (pRow[j] > dTmpMax) dTmpMax = pRow[j];
			}
		}
	}

	*pMin = (double)dTmpMin;
	*pMax = (double)dTmpMax;
	if (pRow) delete[] pRow;
	if (header) delete header;
	if (data) delete data;
	return(true);
}


bool TransformPredictors( char *pParams, bool DoBatch, FPTR fptr)
{
	//if (!DoBatch) Message(pParams, "At entry to TransformPredictors");

	return(TransposeGridsFromParamFile(pParams, fptr));
}


//
// Transform geographic distance to a GDM Transformed grid if DoEuclidean is defined 
// in the model and transform binary grid predictors used in tyhe GDM model.
//
bool TransposeGridsFromParamFile( char *pParams, FPTR fptr )
{
	if (fptr) fptr("Transforming Predictors...", 0);

	//
	// do the Euclidean predictors first, if there are any.....
	//	
	if ( 1 == GetProfileInt( "GDMODEL", "UseEuclidean", pParams ) )
	{
		//
		// use the domain grid to create Easting and Northing grids
		// transposed by the euclidean coeficients in the param file
		//		
		if ( false == CreateEastingGrid( pParams, fptr ) )
		{
			Message( "Cannot CreateEastingGrid", "TransposeGridsFromParamFile()" );
			return false;
		}

		if ( false == CreateNorthingGrid( pParams, fptr ) )
		{
			Message( "Cannot CreateNorthingGrid", "TransposeGridsFromParamFile()" );
			return false;
		}		
	}	

	//
	// if there aren't any euclidean grids, then remove 
	// any residue references from the parameter file
	//
	else
	{
		SetProfileString( "TRANSPREDS", "EuclXTran", NULL, pParams );			
		SetProfileString( "TRANSPREDS", "EuclYTran", NULL, pParams );			

		SetProfileString( "ZEROPREDS", "EuclXTran", NULL, pParams );			
		SetProfileString( "ZEROPREDS", "EuclYTran", NULL, pParams );			
	}

	//
	// using the parameters file, loop thru any predictors that
	// have NO-ZERO coefficients and create a GDM transposed
	// Binary grid using the Domain grid as a mask.
	//
	if ( false == CreateEnvPredictorGrids( pParams, fptr ) )
	{
		Message( "Cannot CreateEnvPredictorGrids", "TransposeGridsFromParamFile()" );
		return false;
	}


	//
	// Loop thru the _ERR grids that are produced for each non-zero-coefficient-sum transform
	// and create a grid that is the sum of the absolute values contained in all the transformed _ERR grids
	//
	if (false == CreateSummedAbsERRGrid( pParams, fptr))
	{
		Message( "Cannot CreateSummedAbsERRGrid", "TransposeGridsFromParamFile()" );
		return false;
	}


	//fptr("Done Transforming Predictors...", 0);
	return( true );
}



//
// Loop thru the _ERR grids that are produced for each non-zero-coefficient-sum transform
// and create a grid that is the sum of the absolute values contained in all the transformed _ERR grids
//
bool CreateSummedAbsERRGrid( char *pParams, FPTR fptr)
{
	//
	// transformed grids are written to the Workspace
	//
	char pWorkspacePath[BUFFLEN];
	GetProfileString( "GDMODEL", "WorkspacePath", pWorkspacePath, pParams );

	//
	// create the output paths
	//
	char pSumGridFLT[BUFFLEN];
	//char pSumGridHDR[BUFFLEN];
	sprintf( pSumGridFLT, "%s\\ABS_ERR_SUM.flt", pWorkspacePath);								
	//sprintf( pSumGridHDR, "%s\\ABS_ERR_SUM.hdr", pWorkspacePath);								

	//
	// Count the _ERR Grids
	//
	char pGridName[BUFFLEN];
	char pGridPath[BUFFLEN];
	char lpEnvGrid[BUFFLEN];
	char OutlierFLT[BUFFLEN];
	char myName[64];
	char pKey[64];
	int nErrGridCount = 0;

	int nPreds = GetProfileInt( "PREDICTORS", "NumPredictors", pParams );
	for ( int i=1; i<=nPreds; i++ )
	{
		// get the environmental grid path
		sprintf( pKey, "EnvGrid%d", i );
		GetProfileString( "PREDICTORS", pKey, lpEnvGrid, pParams );			
		
		//
		// we ONLY transform grids that have NON-ZERO coefficients
		//
		if ( bPredictorHasNonZeroCoeffs( pParams, i ) )
		{
			switch(GetPredictorTypeAt(pParams, i))
			{
				case PRED_TYPE_GRID:
					// construct path to this transform grid set
					_splitpath(lpEnvGrid,NULL,NULL,myName,NULL);
					sprintf( pGridName, "%sTran", myName );
					sprintf( pGridPath, "%s\\%s", pWorkspacePath, pGridName );								
					sprintf(OutlierFLT, "%s_ERR.flt", pGridPath);
					gmpath gmPath;
					if (gmPath.FileExists(OutlierFLT)) ++nErrGridCount;
			}
		}		
	}
	//Message(nErrGridCount, "nErrGridCount");
	if (nErrGridCount > 0)
	{
		//Message(pSumGridFLT, "pSumGridFLT");
		//Message(pSumGridHDR, "pSumGridHDR");
		bool HaveFirstTran = false;
		for ( int i=1; i<=nPreds; i++ )
		{
			// get the environmental grid path
			sprintf( pKey, "EnvGrid%d", i );
			GetProfileString( "PREDICTORS", pKey, lpEnvGrid, pParams );			
		
			//
			// we ONLY transform grids that have NON-ZERO coefficients
			//
			if ( bPredictorHasNonZeroCoeffs( pParams, i ) )
			{
				switch(GetPredictorTypeAt(pParams, i))
				{
					case PRED_TYPE_GRID:
						// construct path to this transform grid set
						_splitpath(lpEnvGrid,NULL,NULL,myName,NULL);
						sprintf( pGridName, "%sTran", myName );
						sprintf( pGridPath, "%s\\%s", pWorkspacePath, pGridName );								
						sprintf(OutlierFLT, "%s_ERR.flt", pGridPath);
						gmpath gmPath;
						if (gmPath.FileExists(OutlierFLT))
						{
							if ( HaveFirstTran )
							{
								// sum the absolute values of _ERR values
								// create the initialised output grid
								if (false == SumAbsoluteGridValues(pSumGridFLT, OutlierFLT, fptr))
								{
									Message( "Cannot SumAbsoluteGridValues", "CreateSummedAbsERRGrid()" );
									return( false );
								}
							}
							else
							{
								// create the initialised output grid
								if (false == InitBinaryGridFromGrid(pSumGridFLT, OutlierFLT, fptr))
								{
									Message( "Cannot InitBinaryGridFromGrid", "CreateSummedAbsERRGrid()" );
									return( false );
								}
								
								// now sum the absolute values of _ERR values
								// create the initialised output grid
								if (false == SumAbsoluteGridValues(pSumGridFLT, OutlierFLT, fptr))
								{
									Message( "Cannot SumAbsoluteGridValues", "CreateSummedAbsERRGrid()" );
									return( false );
								}

								// reset flag
								HaveFirstTran = true;
							}
						} // if (gmPath.FileExists(OutlierFLT))
				} // switch(GetPredictorTypeAt(pParams, i))
			} // if ( bPredictorHasNonZeroCoeffs( pParams, i ) )		
		} // for ( int i=1; i<=nPreds; i++ )
	} // if (nErrGridCount > 0) 
	return(true);
}


//
// add the valid Outlier values as absolute values to the values in Sum Grid
//
bool SumAbsoluteGridValues(char *pSumGridFLT, char *OutlierFLT, FPTR fptr)
{
	gmpath gmPath;
	//
	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(pSumGridFLT, ".hdr"));
	//
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	float fNoData = header->GetNoDataValue();
	//
	BinaryFileClass *bfc_In = new BinaryFileClass(gmPath.ChangeExtension(OutlierFLT, ".flt"), BFC_ReadOnly);	
	//
	BinaryFileClass *bfc_Out = new BinaryFileClass(gmPath.ChangeExtension(pSumGridFLT, ".flt"), BFC_ReadWrite);	
	bfc_Out->SeekToStart();
	//
	float *pInData = new float [nCols];
	float *pOutData = new float [nCols];

	int nCurrent = 0;
#ifdef	DO_PROGRESS
	if (fptr) fptr("Summing ERR Grids...", 0);
#endif
	for ( int i=0; i<nRows; i++ )
	{
#ifdef	DO_PROGRESS
		if ( i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			if (fptr) fptr("Summing ERR Grids...", nCurrent);
		}
#endif

		bfc_In->ReadFloat(pInData, nCols);
		bfc_Out->ReadFloat(pOutData, nCols);

		for ( int j=0; j<nCols; j++ )
		{
			if (pOutData[j] != fNoData)
			{
				if (pInData[j] > 0.0F)
				{
					pOutData[j] += pInData[j];
				}
				else
				{
					pOutData[j] -= pInData[j];   // double minus acts as absolute value
				}
			}
		}
		bfc_Out->SeekBackward(nCols * sizeof(float));
		bfc_Out->WriteFloat(pOutData, nCols);
	}
	bfc_In->Close();
	bfc_Out->Close();
	if (pInData) delete[] pInData;
	if (pOutData) delete[] pOutData;
	return(true);
}



//
// Create a binary grid by applying the geographic transform coordinates to the easting coodinates
//
bool CreateEastingGrid( char *pParams, FPTR fptr )
{
	//
	// Create the full grid path with the workspace path
	//
	char pWorkspacePath[BUFFLEN];
	GetProfileString( "GDMODEL", "WorkspacePath", pWorkspacePath, pParams );

	char pEastingPath[BUFFLEN];
	sprintf( pEastingPath, "%s\\gdmXtran", pWorkspacePath );

	//
	// Initialilise XTran grid to the same metrics as the domain grid
	// No-Data in the domain is No-Data in XTran and 
	// data in Domain is set to 0.0 in XTran
	//
	char pDomainPath[BUFFLEN];
	if (!GetGridAsDomain(pParams, pDomainPath))
	{
		Message( "Cannot Extract Domain Grid Path", "CreateEastingGrid()" );
		return( false );
	}

	if (false == InitBinaryGridFromDomain(pEastingPath, pDomainPath, fptr))
	{
		Message( "Cannot InitBinaryGridFromDomain", "CreateEastingGrid()" );
		return( false );
	}

	//
	// Apply the spline coefficients from the GDM model
	//
	if (false == ApplyEastWestGradient( pEastingPath, pParams, fptr ))
	{
		Message( "Cannot ApplyEastWestGradient", "CreateEastingGrid()" );
		return( false );
	}

	//
	// update the parameter file
	//
	SetProfileString( "TRANSPREDS", "EuclXTran", pEastingPath, pParams );
	return(true);
}



//
// Create a binary grid by applying the geographic transform coordinates to the northing coodinates
//
bool CreateNorthingGrid( char *pParams, FPTR fptr )
{
	//
	// Create the full grid path with the workspace path
	//
	char pWorkspacePath[BUFFLEN];
	GetProfileString( "GDMODEL", "WorkspacePath", pWorkspacePath, pParams );

	char pNorthingPath[BUFFLEN];
	sprintf( pNorthingPath, "%s\\gdmYtran", pWorkspacePath );

	//
	// Initialilise XTran grid to the same metrics as the domain grid
	// No-Data in the domain is No-Data in XTran and 
	// data in Domain is set to 0.0 in XTran
	//
	char pDomainPath[BUFFLEN];
	if (!GetGridAsDomain(pParams, pDomainPath))
	{
		Message( "Cannot Extract Domain Grid Path", "CreateEastingGrid()" );
		return( false );
	}

	if (false == InitBinaryGridFromDomain(pNorthingPath, pDomainPath, fptr))
	{
		Message( "Cannot InitBinaryGridFromDomain", "CreateEastingGrid()" );
		return( false );
	}

	//
	// Apply the spline coefficients from the GDM model
	//
	if ( false == ApplyNorthSouthGradient( pNorthingPath, pParams, fptr ) )
	{
		Message( "Cannot ApplyNorthSouthGradient", "CreateNorthingGrid()" );
		return( false );
	}

	//
	// update the parameter file
	//
	SetProfileString( "TRANSPREDS", "EuclYTran", pNorthingPath, pParams );
	return(true);
}



//
// Create the transformed predictors grids from grid, covariance grid or categorical grid lookup
//
bool CreateEnvPredictorGrids( char *pParams, FPTR fptr )
{
	//Message("At Entry: CreateEnvPredictorGrids()", "INFO");
	//
	// transformed grids are written to the Workspace
	//
	char pWorkspacePath[BUFFLEN];
	GetProfileString( "GDMODEL", "WorkspacePath", pWorkspacePath, pParams );

	int nPreds = GetProfileInt( "PREDICTORS", "NumPredictors", pParams );
	for ( int i=1; i<=nPreds; i++ )
	{
		char pGridName[BUFFLEN];
		char pGridPath[BUFFLEN];
		char myName[64];

		// get the environmental grid path
		char lpEnvGrid[BUFFLEN];
		char pKey[64];
		sprintf( pKey, "EnvGrid%d", i );
		GetProfileString( "PREDICTORS", pKey, lpEnvGrid, pParams );			
		
		//Message(lpEnvGrid, "environmental grid path");

		//
		// we ONLY transform grids that have NON-ZERO coefficients
		//
		if ( bPredictorHasNonZeroCoeffs( pParams, i ) )
		{
			//Message("bPredictorHasNonZeroCoeffs", "INFO");

			switch(GetPredictorTypeAt(pParams, i))
			{
				case PRED_TYPE_GRID:
					// construct path to this transform grid set
					_splitpath(lpEnvGrid,NULL,NULL,myName,NULL);
					sprintf( pGridName, "%sTran", myName );
					sprintf( pGridPath, "%s\\%s", pWorkspacePath, pGridName );								

					//
					// just copy the environment grid to initialise
					//
					if (false == InitBinaryGridFromGrid(pGridPath, lpEnvGrid, fptr))
					{
						Message( "Cannot InitBinaryGridFromGrid", "CreateEnvPredictorGrids()" );
						return( false );
					}


					//Message(pGridPath, "pGridPath");


					//
					// Apply the spline coefficients from the GDM model
					//
					if ( false == ApplySplineGradient( pGridPath, pParams, i, fptr ) )
					{
						Message( "Cannot ApplySplineGradient", "CreateEnvPredictorGrids()" );
						return( false );
					}

					//
					// update the parameter file
					//
					sprintf( pKey, "PredTran%d", i );
					SetProfileString( "TRANSPREDS", pKey, pGridPath, pParams );
					break;

				case PRED_TYPE_COVARIANT:
				case PRED_TYPE_TABLE:
					//
					// skip any creation of transforms...
					//
					//// construct path to this transform grid set
					//_splitpath(lpEnvGrid,NULL,NULL,myName,NULL);
					//sprintf( pGridName, "%sTran", myName );
					//sprintf( pGridPath, "%s\\%s", pWorkspacePath, pGridName );			

					////
					//// just copy the environment grid to initialise
					////
					//if (false == InitBinaryGridFromGrid(pGridPath, lpEnvGrid, fptr))
					//{
					//	Message( "Cannot InitBinaryGridFromGrid", "CreateEnvPredictorGrids()" );
					//	return( false );
					//}

					////
					//// Apply the spline coefficients from the GDM model
					////
					//if ( false == ApplySplineGradient( pGridPath, pParams, i, fptr ) )
					//{
					//	Message( "Cannot ApplySplineGradient", "CreateEnvPredictorGrids()" );
					//	return( false );
					//}

					////
					//// update the parameter file
					////
					//sprintf( pKey, "PredTran%d", i );
					//SetProfileString( "TRANSPREDS", pKey, pGridPath, pParams );
					break;


				case PRED_TYPE_CATEGORICAL:
					char lpCatGrid[BUFFLEN];
					char lpLookup[BUFFLEN];
					char *p = strtok(lpEnvGrid, "+");
					strcpy(lpCatGrid, p);
					//Message(lpCatGrid, "lpCatGrid");

					p = strtok(NULL, "\n");
					strcpy(lpLookup, p);
					// Message(lpLookup, "lpLookup");


					//
					// copy the environment grid to initialise
					//
					_splitpath(lpCatGrid,NULL,NULL,myName,NULL);
					sprintf( pGridName, "%sTran", myName );
					sprintf( pGridPath, "%s\\%s", pWorkspacePath, pGridName );			
					if (false == InitBinaryGridFromGrid(pGridPath, lpCatGrid, fptr))
					{
						Message( "Cannot InitBinaryGridFromGrid", "CreateEnvPredictorGrids()" );
						return( false );
					}

										
					//
					// Create the transformed lookup table
					//
					char lpTranLookup[BUFFLEN];
					_splitpath(lpLookup,NULL,NULL,myName,NULL);
					sprintf( pGridName, "%sTran.csv", myName );
					sprintf( lpTranLookup, "%s\\%s", pWorkspacePath, pGridName );			
					if (false == TransformLookupTable(lpTranLookup, lpLookup, pParams, i))
					{
						Message( "Cannot TransformLookupTable", "CreateEnvPredictorGrids()" );
						return( false );
					}
					

					//
					// update the parameter file
					//
					sprintf( pKey, "PredTran%d", i );
					sprintf(pGridPath, "%s+%s", pGridPath, lpTranLookup);
					SetProfileString( "TRANSPREDS", pKey, pGridPath, pParams );
					break;
			} // switch(GetPredictorTypeAt(pParams, i))
		} // if ( bPredictorHasNonZeroCoeffs( pParams, i ) )

		else  // if the predictor has coefficient sum of zero
		{
			// if it was a grid predictor
			if (PRED_TYPE_GRID == GetPredictorTypeAt(pParams, i))
			{
				// and it was regressed in the model (ie didn't contribute anything)
				if (1 == GetPredictorInUseAt(pParams, i))  
				{
					// then create a zero transform grid from the domain for Tom Harward
					_splitpath(lpEnvGrid,NULL,NULL,myName,NULL);
					sprintf( pGridName, "%sZero", myName );
					sprintf( pGridPath, "%s\\%s", pWorkspacePath, pGridName );								

					if (!InitBinaryGridFromDomain( pGridPath, lpEnvGrid, fptr ))
					{
						Message("Cannot create a zero coefficient grid for Tom Harward", "Oops?!?");
					}

					//
					// update the parameter file
					//
					sprintf( pKey, "ZeroTran%d", i );
					SetProfileString( "ZEROPREDS", pKey, pGridPath, pParams );
				}
			}
		}
	} // for ( int i=1; i<=nPreds; i++ )

	return(true);
}



//
// Transform lookup table with GDM coefficients
//
bool TransformLookupTable(char *To, char *From, char *pParams, int nIndex)
{
	//
	// setup the spline coefficients
	//
	char pKey[64];
	sprintf( pKey, "PredSpl%d", nIndex );
	int nSplines = GetProfileInt( "PREDICTORS", pKey, pParams );
	double *pCoeffs = new double [ nSplines ];
	for ( int i=0; i<nSplines; i++ )
	{
		sprintf( pKey, "PredCoef%d.%d", nIndex, i+1 );
		pCoeffs[i] = GetProfileDouble( "PREDICTORS", pKey, pParams );
		//Message(pCoeffs[i], pKey);
	}


	//
	// setup the spline quantiles
	//
	double *pQuants = new double [ nSplines ];
	for ( int i=0; i<nSplines; i++ )
	{
		sprintf( pKey, "PredSplVal%d.%d", nIndex, i+1 );
		pQuants[i] = GetProfileDouble( "PREDICTORS", pKey, pParams );
		//Message(pQuants[i], pKey);
	}


	//
	// Setup a matrix from the lookup file
	//
	int nRows = 0;
	char *buff = new char [TABLE_ROW_BUFFSIZE];
	FILE *fp = fopen(From, "r+t");
	fgets(buff, TABLE_ROW_BUFFSIZE, fp); // get header
	while(1)
	{
		if (NULL == fgets(buff, TABLE_ROW_BUFFSIZE, fp))
			break;

		++nRows;
	}
	double **ppData = new double * [nRows];
	for ( int i=0; i<nRows; i++ )
	{
		ppData[i] = new double [nRows];
	}
	rewind(fp);
	fgets(buff, TABLE_ROW_BUFFSIZE, fp); // get header
	for ( int i=0; i<nRows; i++ )
	{
		fgets(buff, TABLE_ROW_BUFFSIZE, fp);
		char *p = strtok(buff, ",\n");

		for ( int j=0; j<nRows; j++ )
		{
			p = strtok(NULL, ",\n");
			ppData[i][j] = atof(p);
		}
	}
	fclose(fp);


	// 
	// Transform the matrix
	//
	for ( int i=0; i<nRows; i++ )
	{
		for ( int j=0; j<nRows; j++ )
		{
			if ( i != j )
			{
				double dCalc = 0.0;
				double dVal = ppData[i][j];

				for ( int s=0; s<nSplines; s++ )
				{
					// skip calcs for zero coefficient
					if ( pCoeffs[s] == 0.0 ) continue;

					// otherwise got something to do here
					else if ( s == 0 )
					{
						double d0 = pCoeffs[s] * DoTranSplineCalc( dVal, pQuants[0], pQuants[0], pQuants[1] );
						dCalc += d0;
					}

					else if ( s == (nSplines-1) )
					{
						double d2 = pCoeffs[s] * DoTranSplineCalc( dVal, pQuants[nSplines-2], pQuants[nSplines-1], pQuants[nSplines-1] );
						dCalc += d2;
					}

					else
					{
						double d1 = pCoeffs[s] * DoTranSplineCalc( dVal, pQuants[s-1], pQuants[s], pQuants[s+1] );
						dCalc += d1;
					}
				}
				ppData[i][j] = dCalc;
			}
		}
	}


	//
	// Write the output lookup table
	//
	fp = fopen(To, "w+t");
	// write header
	fprintf(fp, "0,");
	for (int i=0; i<nRows; i++ )
	{
		fprintf(fp, "Tran_Class_%d", i+1);
		if (i < nRows-1)
			fprintf(fp, ",");
		else
			fprintf(fp, "\n");
	}
	// write data
	for (int i=0; i<nRows; i++ )
	{
		fprintf(fp, "Tran_Class_%d,", i+1);
		for ( int j=0; j<nRows; j++ )
		{
			fprintf(fp, "%lf", ppData[i][j]);
			if (j < nRows-1)
				fprintf(fp, ",");
			else
				fprintf(fp, "\n");
		}
	}
	fclose(fp);
	if (buff) delete[] buff;
	for ( int i=0; i<nRows; i++ ) if (ppData[i] ) delete[] ppData[i];
	if (ppData) delete[] ppData;
	return(true);
}


//
// Return true if a predictor has ONE or MORE non-zero coefficients
//
bool bPredictorHasNonZeroCoeffs( char *pParams, int index )
{
	char pKey[64];
	sprintf( pKey, "PredSpl%d", index );
	int nSplines = GetProfileInt( "PREDICTORS", pKey, pParams );

	double dCoeffSum = 0.0;
	for ( int i=1; i<=nSplines; i++ )
	{
		sprintf( pKey, "PredCoef%d.%d", index, i );
		dCoeffSum += GetProfileDouble( "PREDICTORS", pKey, pParams );
	}
	return( ( dCoeffSum > 0.0 ) ? true : false );
}


//
// Extract the first environmental grid path in the parameter file to use as a domain path
//
bool GetGridAsDomain(char *pParams, char *pDomainPath)
{
	GetProfileString("PREDICTORS", "EnvGrid1", pDomainPath, pParams);
	return(pDomainPath != NULL);
}



//
// Create a Binary grid mapped to a Domain grid initialised to 0.0F or No-Data
//
bool InitBinaryGridFromDomain( char *pPath, char *pDomain, FPTR fptr )
{
	gmpath gmPath;
	//
	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(pDomain, ".hdr"));
	//
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	float fNoData = header->GetNoDataValue();
	//
	BinaryFileClass *bfc_Domain = new BinaryFileClass(gmPath.ChangeExtension(pDomain, ".flt"), BFC_ReadOnly);	
	//
	BinaryFileClass *bfc_Out = new BinaryFileClass(gmPath.ChangeExtension(pPath, ".flt"), BFC_CreateOrTrunc);	
	//
	float *pData = new float [nCols];
	int nCurrent = 0;
	char Banner[128];
	sprintf(Banner, "Initialising Binary Grid <%s>", gmPath.GetName(pPath));
#ifdef	DO_PROGRESS
	if (fptr) fptr(Banner, nCurrent);
#endif
	for ( int i=0; i<nRows; i++ )
	{
#ifdef	DO_PROGRESS
		if ( i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			if (fptr) fptr(Banner, nCurrent);
		}
#endif

		// read a row
		bfc_Domain->ReadFloat( pData, nCols );

		for ( int j=0; j<nCols; j++ )
		{
			if ( pData[j] != fNoData ) pData[j] = 0.0F;
		}

		// write a row
		bfc_Out->WriteFloat( pData, nCols );
	}
	bfc_Domain->Close();
	bfc_Out->Close();
	header->CopyTo(gmPath.ChangeExtension(pPath, ".hdr"));

	//
	// clean up
	//
	if (pData) delete[] pData;
	if (header) delete header;
	return( true );
}


//
// Create a Binary grid copied from another grid
//
bool InitBinaryGridFromGrid( char *pPath, char *pGrid, FPTR fptr )
{
	gmpath gmPath;
	//
	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(pGrid, ".hdr"));
	//
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	float fNoData = header->GetNoDataValue();
	//
	BinaryFileClass *bfc_Domain = new BinaryFileClass(gmPath.ChangeExtension(pGrid, ".flt"), BFC_ReadOnly);	
	//
	BinaryFileClass *bfc_Out = new BinaryFileClass(gmPath.ChangeExtension(pPath, ".flt"), BFC_CreateOrTrunc);	
	//
	float *pData = new float [nCols];
	int nCurrent = 0;
	char Banner[128];
	sprintf(Banner, "Initialising Binary Grid <%s>", gmPath.GetName(pPath));
#ifdef	DO_PROGRESS
	if (fptr) fptr(Banner, nCurrent);
#endif
	for ( int i=0; i<nRows; i++ )
	{
#ifdef	DO_PROGRESS
		if ( i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			if (fptr) fptr(Banner, nCurrent);
		}
#endif

		// read a row
		bfc_Domain->ReadFloat( pData, nCols );

		// write a row
		bfc_Out->WriteFloat( pData, nCols );
	}
	bfc_Domain->Close();
	bfc_Out->Close();
	header->CopyTo(gmPath.ChangeExtension(pPath, ".hdr"));

	//
	// clean up
	//
	if (pData) delete[] pData;
	if (header) delete header;
	return( true );
}



//
// Apply the GDM transformation to the East-West Geographic predictor
//
bool ApplyEastWestGradient( char *pPath, char *pParams, FPTR fptr )
{
	gmpath gmPath;
	//
	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(pPath, ".hdr"));
	//
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	double ResX = header->GetCellSize();
	float fNoData = header->GetNoDataValue();
	//
	BinaryFileClass *bfc_Out = new BinaryFileClass(gmPath.ChangeExtension(pPath, ".flt"), BFC_ReadWrite);	
	//
	float *pData = new float [nCols];
	float fCellSize = (float)ResX;

	
	//
	// Setup a linear scalar to approximate the general slope of 
	// the model curve to apply to the distance from the left edge. 
	// Distances exceeding this are simply extrapolated by the linear scalar.
	//
	double dTranX = 0.0;
	double dTranY = 0.0;

	int nSplines = GetProfileInt( "PREDICTORS", "EuclSpl", pParams );
	char pKey[64];
	for ( int i=0; i<nSplines; i++ )
	{
		sprintf( pKey, "EuclCoef%d", i+1 );
		dTranY += GetProfileDouble( "PREDICTORS", pKey, pParams );

		// get the maximum quantile value here
		if ( i == nSplines-1 )
		{
			sprintf( pKey, "EuclSplVal%d", i+1 );
			dTranX = GetProfileDouble( "PREDICTORS", pKey, pParams );
		}
	}

	float fScalar = (float)( dTranY / dTranX );

	//
	// Cycle thru the data, cell by cell, calculating the distance
	// after applying the linear scalar to the distances
	// and resetting values back into data block.
	// 
	//
	char Banner[128];
	sprintf(Banner, "Transforming Binary Grid <%s>", gmPath.GetName(pPath));
	int nCurrent = 0;
	for ( int i=0; i<nRows; i++ )
	{
		if ( i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
#ifdef	DO_PROGRESS			
			if (fptr) fptr(Banner, nCurrent);
#endif
		}

		// read a row of data
		bfc_Out->ReadFloat(pData, nCols);	

		for ( int j=0; j<nCols; j++ )
		{
			// Data cells will currently have the value of 0.0
			if (pData[j] == 0.0)
			{
				// calculate data value as distance from east edge
				pData[j] = fScalar * j * fCellSize;
			}
		}

		// write a row of data
		bfc_Out->SeekBackward(nCols * sizeof(float));
		bfc_Out->WriteFloat(pData, nCols);	
	}
	bfc_Out->Close();

	//
	// clean up
	//
	if (pData) delete[] pData;
	if (header) delete header;
	return( true );
}


//
// Apply the GDM transformation to the East-West Geographic predictor
//
bool ApplyNorthSouthGradient( char *pPath, char *pParams, FPTR fptr )
{
	gmpath gmPath;
	//
	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(pPath, ".hdr"));
	//
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	double ResX = header->GetCellSize();
	float fNoData = header->GetNoDataValue();
	//
	BinaryFileClass *bfc_Out = new BinaryFileClass(gmPath.ChangeExtension(pPath, ".flt"), BFC_ReadWrite);	
	//
	float *pData = new float [nCols];
	float fCellSize = (float)ResX;

	//
	// Setup a linear scalar to approximate the general slope of 
	// the model curve to apply to the distance from the left edge. 
	// Distances exceeding this are simply extrapolated by the linear scalar.
	//
	double dTranX = 0.0;
	double dTranY = 0.0;

	int nSplines = GetProfileInt( "PREDICTORS", "EuclSpl", pParams );
	char lpKey[64];
	for ( int i=0; i<nSplines; i++ )
	{
		sprintf( lpKey, "EuclCoef%d", i+1 );
		dTranY += GetProfileDouble( "PREDICTORS", lpKey, pParams );

		// get the maximum quantile value here
		if ( i == nSplines-1 )
		{
			sprintf( lpKey, "EuclSplVal%d", i+1 );
			dTranX = GetProfileDouble( "PREDICTORS", lpKey, pParams );
		}
	}

	float fScalar = (float)( dTranY / dTranX );

	//
	// Cycle thru the data, cell by cell, calculating the distance
	// after applying the linear scalar to the distances
	// and resetting values back into data block.
	// 
	//
	char Banner[128];
	sprintf(Banner, "Transforming Binary Grid <%s>", gmPath.GetName(pPath));
	int nCurrent = 0;
	for ( int i=0; i<nRows; i++ )
	{
#ifdef	DO_PROGRESS
		if ( i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			if (fptr) fptr(Banner, nCurrent);
		}
#endif

		// read a row
		bfc_Out->ReadFloat( pData, nCols);	

		for ( int j=0; j<nCols; j++ )
		{
			// Data cells will currently have the value of 0.0
			if ( pData[j] == 0.0 )
			{
				// calculate data value as distance from Top edge
				pData[j] = fScalar * i * fCellSize;
			}
		}

		// write a row
		bfc_Out->SeekBackward(nCols * sizeof(float));
		bfc_Out->WriteFloat(pData, nCols);	
	}
	bfc_Out->Close();

	//
	// clean up
	//
	if (pData) delete[] pData;
	if (header) delete header;
	return( true );
}


//
// Apply the GDM transformation to a grid predictor
//
bool ApplySplineGradient( char *pPath, char *pParams, int nIndex, FPTR fptr )
{
	//Message(GetProfileInt("GDMODEL", "ExtrapolationMethod", pParams), "Method");
	
	//
	// setup the spline coefficients
	//
	char pKey[64];
	sprintf( pKey, "PredSpl%d", nIndex );
	int nSplines = GetProfileInt( "PREDICTORS", pKey, pParams );
	double *pCoeffs = new double [ nSplines ];
	for ( int i=0; i<nSplines; i++ )
	{
		sprintf( pKey, "PredCoef%d.%d", nIndex, i+1 );
		pCoeffs[i] = GetProfileDouble( "PREDICTORS", pKey, pParams );
		//Message(pCoeffs[i], pKey);
	}


	//
	// setup the spline quantiles
	//
	double *pQuants = new double [ nSplines ];
	for ( int i=0; i<nSplines; i++ )
	{
		sprintf( pKey, "PredSplVal%d.%d", nIndex, i+1 );
		pQuants[i] = GetProfileDouble( "PREDICTORS", pKey, pParams );
	}

	double dMin = pQuants[0];
	double dMid = pQuants[nSplines/2];
	double dMax = pQuants[nSplines-1];

	//
	// setup for possible outliers less than pQuants[0](dMin)
	// or greater than pQuants[nSplines-1](dMax)
	//
	double dFullSlopeGap = (dMax - dMin);
	double dSlopeGap = ( dMax - dMin ) / 10.0;
	double dVal10 = dMin + dSlopeGap;
	double dVal90 = dMax - dSlopeGap;
	
	double dCalc10 = 0.0;
	double dCalc90 = 0.0;
	double dCalcMax = 0.0;
	double dCalcMin = 0.0;

	for ( int s=0; s<nSplines; s++ )
	{
		if ( s == 0 )
		{
			dCalcMin += pCoeffs[s] * DoTranSplineCalc( dMin, pQuants[0], pQuants[0], pQuants[1] );
			dCalc10  += pCoeffs[s] * DoTranSplineCalc( dVal10, pQuants[0], pQuants[0], pQuants[1] );
			dCalc90  += pCoeffs[s] * DoTranSplineCalc( dVal90, pQuants[0], pQuants[0], pQuants[1] );
			dCalcMax += pCoeffs[s] * DoTranSplineCalc( dMax, pQuants[0], pQuants[0], pQuants[1] );
		}

		else if ( s == (nSplines-1) )
		{
			dCalcMin += pCoeffs[s] * DoTranSplineCalc( dMin, pQuants[nSplines-2], pQuants[nSplines-2], pQuants[nSplines-1] );
			dCalc10  += pCoeffs[s] * DoTranSplineCalc( dVal10, pQuants[nSplines-2], pQuants[nSplines-2], pQuants[nSplines-1] );
			dCalc90  += pCoeffs[s] * DoTranSplineCalc( dVal90, pQuants[nSplines-2], pQuants[nSplines-1], pQuants[nSplines-1] );
			dCalcMax += pCoeffs[s] * DoTranSplineCalc( dMax, pQuants[nSplines-2], pQuants[nSplines-1], pQuants[nSplines-1] );
		}

		else
		{
			dCalcMin += pCoeffs[s] * DoTranSplineCalc( dMin, pQuants[s-1], pQuants[s], pQuants[s+1] );
			dCalc10  += pCoeffs[s] * DoTranSplineCalc( dVal10, pQuants[s-1], pQuants[s], pQuants[s+1] );
			dCalc90  += pCoeffs[s] * DoTranSplineCalc( dVal90, pQuants[s-1], pQuants[s], pQuants[s+1] );
			dCalcMax += pCoeffs[s] * DoTranSplineCalc( dMax, pQuants[s-1], pQuants[s], pQuants[s+1] );
		}
	}

	double dFullSlope = (dCalcMax - dCalcMin) / dFullSlopeGap;
	double dMinSlope = ( dCalc10 - dCalcMin ) / dSlopeGap;
	double dMaxSlope = ( dCalcMax - dCalc90 ) / dSlopeGap;

	gmpath gmPath;
	//
	EsriBinaryHeader *header = new EsriBinaryHeader(gmPath.ChangeExtension(pPath, ".hdr"));
	//
	int nRows = header->GetNumRows();
	int nCols = header->GetNumCols();
	float fNoData = header->GetNoDataValue();
	//
	BinaryFileClass *bfc_Out = new BinaryFileClass(gmPath.ChangeExtension(pPath, ".flt"), BFC_ReadWrite);	
	//
	float *pData = new float [nCols];

	//float myMin = 10000000.0F;
	//float myMax = -10000000.0F;

	// outlier value grid path
	char OutlierFLT[512];
	sprintf(OutlierFLT, "%s_ERR.flt", pPath);
	BinaryFileClass *bfc_Err = new BinaryFileClass(OutlierFLT, BFC_CreateOrTrunc);	
	//BinaryFileClass *bfc_Err = new BinaryFileClass(GetOutlierFilename(pPath, ".flt"), BFC_CreateOrTrunc);	
	float *pErrData = new float [nCols];

	char headerr [512];
	//strcpy(headerr, GetOutlierFilename(pPath, ".hdr"));
	sprintf(headerr, "%s_ERR.hdr", pPath);
	header->CopyTo(headerr);

	char Banner[128];
	sprintf(Banner, "Transforming Binary Grid <%s>", gmPath.GetName(pPath));
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Cycle thru the data, cell by cell, calculating spline 
	// transforms and resetting values back into data block.
	// 
	//
	int nCurrent = 0;
	for ( int i=0; i<nRows; i++ )
	{
#ifdef	DO_PROGRESS
		if ( i * 100 / nRows > nCurrent )
		{
			nCurrent = i * 100 / nRows;
			if (fptr) fptr(Banner, nCurrent);
		}
#endif

		// read in the row data
		bfc_Out->ReadFloat(pData, nCols);	

		//
		// The extrapolation modes are set in the parameter file as [] ExtrapolationMode=0,1,2,3
		//
		// GDM_Extrapolation_Default10			0
		// GDM_Extrapolation_WholeGradient		1
		// GDM_Extrapolation_MostConservative	2
		// GDM_Extrapolation_None  			    3
		//
		int XtrapMethod = GetProfileInt("GDMODEL", "ExtrapolationMethod", pParams);
		for ( int j=0; j<nCols; j++ )
		{
			// Data cells will currently have the value of 0.0
			if ( pData[j] != fNoData )
			{
				// Use the existing environmental value as a data point
				double dVal = (double)pData[j];
				double dCalc = 0.0;				

				// calculate outliers less than minimum quantile
				if ( dVal < dMin )
				{
					switch (XtrapMethod)
					{
					case 0:
						dCalc = -(dMin - dVal) * dMinSlope;
						pErrData[j] = (float)(dCalc - dCalcMin);
						break;

					case 1:
						dCalc = -(dMin - dVal) * dFullSlope;
						pErrData[j] = (float)(dCalc - dCalcMin);
						break;

					case 2:
						if (dMinSlope < dFullSlope)
							dCalc = -(dMin - dVal) * dMinSlope;
						else
							dCalc = -(dMin - dVal) * dFullSlope;
						pErrData[j] = (float)(dCalc - dCalcMin);
						break;

					case 3:
						dCalc = dMin;
						pErrData[j] = (float)(dCalc);
						break;

					default: // GDM_Extrapolation_Default10
						dCalc = -(dMin - dVal) * dMinSlope;
						pErrData[j] = (float)(dCalc - dCalcMin);
						break;
					} // switch (XtrapMethod)
				} // if ( dVal < dMin )

				
				// calculate outliers greater than maximum quantile
				else if (dVal > dMax)
				{
					switch (XtrapMethod)
					{
					case 0:
						dCalc = -(dVal - dMax) * dMaxSlope;
						pErrData[j] = (float)(dCalc - dCalcMax);
						break;

					case 1:
						dCalc = -(dVal - dMax) * dFullSlope;
						pErrData[j] = (float)(dCalc - dCalcMax);
						break;

					case 2:
						if (dMinSlope < dFullSlope)
							dCalc = -(dVal - dMax) * dMaxSlope;
						else
							dCalc = -(dVal - dMax) * dFullSlope;
						pErrData[j] = (float)(dCalc - dCalcMax);
						break;

					case 3:
						dCalc = dMax;
						pErrData[j] = (float)(dCalc);
						break;

					default: // GDM_Extrapolation_Default10
						dCalc = -(dVal - dMax) * dMaxSlope;
						pErrData[j] = (float)(dCalc - dCalcMax);
						break;
					} // switch (XtrapMethod)
				} // if ( dVal > dMax )

				else
				{
					for ( int s=0; s<nSplines; s++ )
					{
						// skip calcs for zero coefficient
						if ( pCoeffs[s] == 0.0 ) continue;

						// otherwise got something to do here
						else if ( s == 0 )
						{
							double d0 = pCoeffs[s] * DoTranSplineCalc( dVal, pQuants[0], pQuants[0], pQuants[1] );
							dCalc += d0;
						}

						else if ( s == (nSplines-1) )
						{
							double d2 = pCoeffs[s] * DoTranSplineCalc( dVal, pQuants[nSplines-2], pQuants[nSplines-1], pQuants[nSplines-1] );
							dCalc += d2;
						}

						else
						{
							double d1 = pCoeffs[s] * DoTranSplineCalc( dVal, pQuants[s-1], pQuants[s], pQuants[s+1] );
							dCalc += d1;
						}
					}
					pErrData[j] = 0.0F;
				}

				//if ( myMin > (float)dCalc ) myMin = (float)dCalc;
				//if ( myMax < (float)dCalc ) myMax = (float)dCalc;

				// set the summed spline values in the data block
				pData[j] = (float)dCalc;
			} // if ( pData[j] != fNoData )

			else
			{
				pErrData[j] = fNoData;
			}

		} // for ( int j=0; j<nCols; j++ )

		// write the row
		bfc_Out->SeekBackward(nCols * sizeof(float));
		bfc_Out->WriteFloat(pData, nCols);
		bfc_Err->WriteFloat(pErrData, nCols);
	} // for ( i=0; i<nRows; i++ )	
	
	//
	// clean up
	//
	bfc_Out->Close();
	bfc_Err->Close();
	if (pQuants) delete[] pQuants;
	if (pCoeffs) delete[] pCoeffs;
	if (pData) delete[] pData;
	if (header) delete header;
	return( true );
}


//
// Create a filename from the transform filename for the outlying prediction grid
//
char *GetOutlierFilename(char *filename, char *ext)
{	
	char outDir[BUFFLEN];
	char outGrid[BUFFLEN];
	char outName[BUFFLEN];
	gmpath gmPath;
	strcpy(outDir, gmPath.GetDirectoryPath(filename));
	strcpy(outGrid, gmPath.GetName(filename));
	sprintf(outName, "%s\\%s_ERR%s", outDir, outGrid, ext);
	return(outName);
}


//
// Calculate the I-Spline value for dVal given quantiles q1, q2, q3
//
double DoTranSplineCalc( double dVal, double q1, double q2, double q3 )
{
	if ( dVal <= q1 )
	{
		return( 0.0 );
	}

	else if ( dVal >= q3 ) 
	{
		return( 1.0 );
	}

	else if ( ( q1 < dVal ) && ( dVal < q2 ) )
	{
		return( ( ( ( dVal - q1 ) * ( dVal - q1 ) ) / ( ( q2 - q1 ) * ( q3 - q1 ) ) ) );
	}

	else
	{
		return( ( 1.0 - ( ( ( q3 - dVal ) * ( q3 - dVal ) ) / ( ( q3 - q2 ) * ( q3 - q1 ) ) ) ) );
	}
}


#define Grids_OK   0
#define Diff_Size  1
#define Diff_X_Min 2
#define Diff_X_Max 3
#define Diff_Y_Min 4
#define Diff_Y_Max 5

//
// Test that the two grids have the same extent and CellSize
//
int TestExtentAndCellsize(char *pRef, char *pTest)
{
	EsriBinaryHeader *h0 = new EsriBinaryHeader(pRef);
	EsriBinaryHeader *h1 = new EsriBinaryHeader(pTest);

	double dSize0 = h0->GetCellSize();
	double dSize1 = h1->GetCellSize();

	double dXMin0 = h0->GetMinX();
	double dXMin1 = h1->GetMinX();

	double dXMax0 = h0->GetMaxX();
	double dXMax1 = h1->GetMaxX();

	double dYMin0 = h0->GetMinY();
	double dYMin1 = h1->GetMinY();

	double dYMax0 = h0->GetMaxY();
	double dYMax1 = h1->GetMaxY();

	if (dSize0 != dSize1)
	{
		return(Diff_Size);
	}

	else if (dXMin0 != dXMin1)
	{
		return(Diff_X_Min);
	}

	else if (dXMax0 != dXMax1)
	{
		return(Diff_X_Max);
	}

	else if (dYMin0 != dYMin1)
	{
		return(Diff_Y_Min);
	}

	else if (dYMax0 != dYMax1)
	{
		return(Diff_Y_Max);
		return(false);
	}
	
	else 
		return(Grids_OK);
}
