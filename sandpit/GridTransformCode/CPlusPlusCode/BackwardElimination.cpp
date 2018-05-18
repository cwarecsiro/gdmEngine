//
// BackwardElimination.cpp
//
#include "stdafx.h"
#include "GdmBinLib.h"
#include "myCallback.h"
#include "Message.h"
#include "BackwardElimination.h"
#include "GdmGuiLib.h"
#include "GDMBufferSizeDefs.h"
#include "ParamsW16.h"
#include "clsDoPath.h"
#include "GdmTabLib.h"
#include "NNLS_Double.h"


//
// BATCH MODE ONLY
// Run a GDM progressively dropping each of the predictors
// from the analysis to determine its the contribution to the fit.
//
bool PerformBatchedBackwardEliminationGDM(char *pParams, FPTR fptr, UPTR uptr)
{
	//
	// Do the full GDM and retain the binary matrix file
	// for sub-setting to use in the backward elimination
	//
	double *pResponse = NULL;
	double *pWeights = NULL;
	int NumRows, NumCols;
	if ( !DoBaseGDM(pParams, true, false, &pResponse, &pWeights, &NumRows, &NumCols, fptr) )
	{
		Message("Cannot DoBaseGDM in PerformBatchedBackwardEliminationGDM", "ERROR");
		return(false);
	}
	//Message(NumRows, "NumRows");
	//Message(NumCols, "NumCols");

	if ( pResponse == NULL )
	{
		Message("Got a NULL Response Vector", "ERROR");
		return(false);
	}

	if ( pWeights == NULL )
	{
		Message("Got a NULL Weights Vector", "ERROR");
		return(false);
	}


	//
	// Update the parameter file to start a new run with this result being the benchmark
	//		and the dropped predictor runs below being subsets of this run
	//
	int NumRuns = GetProfileInt( "Elimination_Runs", "NumRuns", pParams );
	if ( NumRuns < 1 )
		NumRuns = 1;
	else
		++NumRuns;
	//Message(NumRuns, "NumRuns");

	// setup the paragraph name for this run
	AppendRow(pParams);
	SetProfileInt( "Elimination_Runs", "NumRuns", NumRuns, pParams );
	AppendRow(pParams);
	
	char myRunName [64];
	sprintf( myRunName, "Elimination_Run_%d", NumRuns );
	SetProfileDouble( myRunName, "NULLDeviance", GetNullDeviance(pParams), pParams );
	SetProfileDouble( myRunName, "GDMDeviance", GetGDMDeviance(pParams), pParams );
	SetProfileDouble( myRunName, "ExplainedDeviance", GetDevianceExplained(pParams), pParams );
	/*SetProfileInt( myRunName, "Columns", nCols, pParams );
	double dCoeffSum = 0.0;
	for ( int i=0; i<nCols; i++ )
	{
		char myKey[64];
		sprintf( myKey, "Coeff%d", i );
		SetProfileDouble( myRunName, myKey, pFullCoeffs[i], pParams );
		dCoeffSum += pFullCoeffs[i];
	}*/
	SetProfileDouble( myRunName, "ElimRunSumOfCoeffs", GetSumOfCoefficients(pParams), pParams );
	AppendRow(pParams);

		
	// update the backward elimination .NET interface
	double dExplainedDeviance = GetDevianceExplained(pParams);
	if (uptr)
		uptr(dExplainedDeviance, dExplainedDeviance, -1, -1, -1);


	// setup path to full GDM binary matrix file
	char *lpTmpFile = GetWorkspacePath(pParams);
	strcat(lpTmpFile, "\\gdmtmp.bin");


	// setup path to partial GDM binary matrix file for use in backward elimination
	char *lpElimFile = GetWorkspacePath(pParams);
	strcat(lpElimFile, "\\tmpPredElim.bin");


	////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Having done the full GDM with all user selected predictors,
	// progressively drop each predictor and update the 
	// amount of deviance explained by this dropped predictor.
	//
	// If the predictor is set to 'Locked' (ie Dropped Predictor value == -1)
	//        then always include it in the model and don't drop it
	//
	//
	double *pCopyCol = new double [NumRows];
	int nPreds = GetNumPredictors(pParams);
	//char myKey [64];
	int hIn, hOut;

	//
	// do the euclidean predictor first
	//
	if ( (UseEuclidean(pParams)) && (0 == GetEuclideanPredDroppedState(pParams)))
	{
		// change colour to processing (orange)
		if (uptr)
			uptr(-5, -5, -5, -5, 0);

		
		// open the binary file for the FULL GDM 
		hIn = _open( lpTmpFile, _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE );
		if (hIn < 0)
		{
			Message("Cannot create open binary file for in PerformBatchedBackwardEliminationGDM()", "ERROR");
			if (lpTmpFile) delete[] lpTmpFile;
			if (lpElimFile) delete[] lpElimFile;
			if (pResponse) delete[] pResponse;
			if (pWeights) delete[] pWeights;
			if (pCopyCol) delete[] pCopyCol;
			return(false);
		}


		// create new binary matrix using intercept column followed by predictor columns
		hOut = _open( lpElimFile, _O_BINARY | _O_CREAT | _O_TRUNC | _O_RDWR, S_IREAD | S_IWRITE );
		if (hOut < 0)
		{
			Message("Cannot create binary file for in PerformBatchedBackwardEliminationGDM()", "ERROR");
			_close(hIn);
			if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );
			if (lpTmpFile) delete[] lpTmpFile;
			if (lpElimFile) delete[] lpElimFile;
			if (pResponse) delete[] pResponse;
			if (pWeights) delete[] pWeights;
			if (pCopyCol) delete[] pCopyCol;
			return(false);
		}

		// copy the intercept
		_lseek(hIn, 0L, SEEK_SET);
		_read(hIn, pCopyCol, NumRows * sizeof(double));
		_write(hOut, pCopyCol, NumRows * sizeof(double));
		int NumTmpCols = 1;
				
		// skip past the euclidean distance columns to drop the euclidean distance predictor
		int nSplines2Drop = GetEuclideanSplines(pParams);
		_lseek(hIn, (long)(nSplines2Drop * NumRows * sizeof(double)), SEEK_CUR);

		// now write the rest of the selected predictor data
		int Cols2Copy = NumCols - nSplines2Drop - 1;
		for ( int i=0; i<Cols2Copy; i++ )
		{
			_read(hIn, pCopyCol, NumRows * sizeof(double));
			_write(hOut, pCopyCol, NumRows * sizeof(double));
			++NumTmpCols;
		}

		_close(hIn);
		_close(hOut);

		if (!DoIncrementalGDM(pParams, true, true, 
					          lpElimFile, myRunName, 
							  pResponse, pWeights, 
							  NumRows, NumTmpCols, 0,
							  fptr, uptr))
		{
			Message("Cannot DoIncrementalGDM in PerformBatchedBackwardEliminationGDM()", "ERROR");
			if (lpTmpFile) delete[] lpTmpFile;
			if (lpElimFile) delete[] lpElimFile;
			if (pResponse) delete[] pResponse;
			if (pWeights) delete[] pWeights;
			if (pCopyCol) delete[] pCopyCol;
			return(false);
		}
	} // if ( (UseEuclidean(pParams)) && (0 == GetEuclideanPredDroppedState(pParams)))


	//
	// now do the non-dropped predictors
	//
	for ( int ThisPred=1; ThisPred<=nPreds; ThisPred++ )
	{
		if ( -1 == GetPredictorDroppedAt(pParams, ThisPred) )
		{
			continue;
		}


		// if this predictor is a candidate to be dropped...
		if ( 0 < GetPredictorInUseAt(pParams, ThisPred) ) 
		{
			// change colour to processing (orange)
			if (uptr)
				uptr(-5, -5, -5, -5, ThisPred);

			// open the binary file for the FULL GDM 
			hIn = _open( lpTmpFile, _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE );
			if (hIn < 0)
			{
				Message("Cannot create open binary file for in PerformBatchedBackwardEliminationGDM()", "ERROR");
				if (lpTmpFile) delete[] lpTmpFile;
				if (lpElimFile) delete[] lpElimFile;
				if (pResponse) delete[] pResponse;
				if (pWeights) delete[] pWeights;
				if (pCopyCol) delete[] pCopyCol;
				return(false);
			}


			// create new binary matrix using intercept column followed by predictor columns
			hOut = _open( lpElimFile, _O_BINARY | _O_CREAT | _O_TRUNC | _O_RDWR, S_IREAD | S_IWRITE );
			if (hOut < 0)
			{
				Message("Cannot create binary file for in PerformBatchedBackwardEliminationGDM()", "ERROR");
				_close(hIn);
				if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );
				if (lpTmpFile) delete[] lpTmpFile;
				if (lpElimFile) delete[] lpElimFile;
				if (pResponse) delete[] pResponse;
				if (pWeights) delete[] pWeights;
				if (pCopyCol) delete[] pCopyCol;
				return(false);
			}

			// copy the intercept
			_lseek(hIn, 0L, SEEK_SET);
			_read(hIn, pCopyCol, NumRows * sizeof(double));
			_write(hOut, pCopyCol, NumRows * sizeof(double));

			int NumTmpCols = 1;

			// if we have geographic distance as a predictor
			if ( (UseEuclidean(pParams)) && (0 == GetEuclideanPredDroppedState(pParams)))
			{
				int nTmpSplines = GetEuclideanSplines(pParams);
				for (int i=0; i<nTmpSplines; i++ )
				{
					_read(hIn, pCopyCol, NumRows * sizeof(double));
					_write(hOut, pCopyCol, NumRows * sizeof(double));
				}
				NumTmpCols += nTmpSplines; // increment the columns count in the regression matrix
			}

			// include all the other 'non-dropped' predictors
			for ( int j=1; j<=nPreds; j++ )
			{
				int nTmpSplines = GetPredictorSplinesAt(pParams, j);

				// ignore if the predictor is NOT 'in-use'
				if ( 0 == GetPredictorInUseAt(pParams, j) )
				{
					continue;
				}

				// don't use the dropped predictor
				if ( ThisPred == j ) 
				{
					_lseek(hIn, (long)(nTmpSplines * NumRows * sizeof(double)), SEEK_CUR);
					continue;
				}

				// else copy the relevent columns to the temporary binary matrix
				for ( int k=0; k<nTmpSplines; k++ )
				{
					_read(hIn, pCopyCol, NumRows * sizeof(double));
					_write(hOut, pCopyCol, NumRows * sizeof(double));
				}
				NumTmpCols += nTmpSplines; // increment the columns count in the regression matrix
			} // for ( int j=1; j<=nPreds; j++ )

			_close(hIn);
			_close(hOut);

			if (!DoIncrementalGDM(pParams, true, true, 
							      lpElimFile, myRunName, 
								  pResponse, pWeights, 
							      NumRows, NumTmpCols, ThisPred,
							      fptr, uptr))
		    {
				Message("Cannot DoIncrementalGDM in PerformBatchedBackwardEliminationGDM()", "ERROR");
				if (lpTmpFile) delete[] lpTmpFile;
				if (lpElimFile) delete[] lpElimFile;
				if (pResponse) delete[] pResponse;
				if (pWeights) delete[] pWeights;
				if (pCopyCol) delete[] pCopyCol;
				return(false);
			} // if (!DoIncrementalGDM(...
		} // if ( 0 < GetPredictorInUseAt(pParams, ThisPred) )
	} // for ( int ThisPred=1; ThisPred<=nPreds; ThisPred++ )


	//
	// Cleanup
	//
	if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );
	if ( ( _access( lpElimFile, 0 ) ) != -1 ) remove( lpElimFile );
	if (lpTmpFile) delete[] lpTmpFile;
	if (lpElimFile) delete[] lpElimFile;
	if (pResponse) delete[] pResponse;
	if (pWeights) delete[] pWeights;
	if (pCopyCol) delete[] pCopyCol;
	return(true);
}



//
// Run a GDM progressively dropping each of the predictors
// from the analysis to determine its the contribution to the fit.
//
bool PerformBackwardEliminationGDM(char *pParams, bool BatchRun, FPTR fptr, UPTR uptr)
{
	//
	// Do the full GDM and retain the binary matrix file
	// for sub-setting to use in the backward elimination
	//
	double *pResponse = NULL;
	double *pWeights = NULL;
	int NumRows, NumCols;
	if ( !DoBaseGDM(pParams, BatchRun, false, &pResponse, &pWeights, &NumRows, &NumCols, fptr) )
	{
		if (!BatchRun)
		{
			Message("Cannot DoBaseGDM in PerformBackwardEliminationGDM", "ERROR");
		}
		return(false);
	}
	//Message(NumRows, "NumRows");
	//Message(NumCols, "NumCols");

	if ( pResponse == NULL )
	{
		Message("Got a NULL Response Vector", "ERROR");
		return(false);
	}

	if ( pWeights == NULL )
	{
		Message("Got a NULL Weights Vector", "ERROR");
		return(false);
	}


	//
	// Update the parameter file to start a new run with this result being the benchmark
	//		and the dropped predictor runs below being subsets of this run
	//
	int NumRuns = GetProfileInt( "Elimination_Runs", "NumRuns", pParams );
	if ( NumRuns < 1 )
		NumRuns = 1;
	else
		++NumRuns;
	//Message(NumRuns, "NumRuns");

	// setup the paragraph name for this run
	AppendRow(pParams);
	SetProfileInt( "Elimination_Runs", "NumRuns", NumRuns, pParams );
	AppendRow(pParams);
	
	char myRunName [64];
	sprintf( myRunName, "Elimination_Run_%d", NumRuns );
	SetProfileDouble( myRunName, "NULLDeviance", GetNullDeviance(pParams), pParams );
	SetProfileDouble( myRunName, "GDMDeviance", GetGDMDeviance(pParams), pParams );
	SetProfileDouble( myRunName, "ExplainedDeviance", GetDevianceExplained(pParams), pParams );
	/*SetProfileInt( myRunName, "Columns", nCols, pParams );
	double dCoeffSum = 0.0;
	for ( int i=0; i<nCols; i++ )
	{
		char myKey[64];
		sprintf( myKey, "Coeff%d", i );
		SetProfileDouble( myRunName, myKey, pFullCoeffs[i], pParams );
		dCoeffSum += pFullCoeffs[i];
	}*/
	SetProfileDouble( myRunName, "ElimRunSumOfCoeffs", GetSumOfCoefficients(pParams), pParams );
	AppendRow(pParams);

		
	// update the backward elimination .NET interface
	double dExplainedDeviance = GetDevianceExplained(pParams);
	if (uptr)
		uptr(dExplainedDeviance, dExplainedDeviance, -1, -1, -1);


	// setup path to full GDM binary matrix file
	char *lpTmpFile = GetWorkspacePath(pParams);
	strcat(lpTmpFile, "\\gdmtmp.bin");


	// setup path to partial GDM binary matrix file for use in backward elimination
	char *lpElimFile = GetWorkspacePath(pParams);
	strcat(lpElimFile, "\\tmpPredElim.bin");


	////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Having done the full GDM with all user selected predictors,
	// progressively drop each predictor and update the 
	// amount of deviance explained by this dropped predictor.
	//
	// If the predictor is set to 'Locked' (ie Dropped Predictor value == -1)
	//        then always include it in the model and don't drop it
	//
	//
	double *pCopyCol = new double [NumRows];
	int nPreds = GetNumPredictors(pParams);
	//char myKey [64];
	int hIn, hOut;

	//
	// do the euclidean predictor first
	//
	if ( (UseEuclidean(pParams)) && (0 == GetEuclideanPredDroppedState(pParams)))
	{
		// change colour to processing (orange)
		if (uptr)
			uptr(-5, -5, -5, -5, 0);

		
		// open the binary file for the FULL GDM 
		hIn = _open( lpTmpFile, _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE );
		if (hIn < 0)
		{
			if (!BatchRun)
				Message("Cannot create open binary file for in PerformBackwardEliminationGDM()", "ERROR");
			if (lpTmpFile) delete[] lpTmpFile;
			if (lpElimFile) delete[] lpElimFile;
			if (pResponse) delete[] pResponse;
			if (pWeights) delete[] pWeights;
			if (pCopyCol) delete[] pCopyCol;
			return(false);
		}


		// create new binary matrix using intercept column followed by predictor columns
		hOut = _open( lpElimFile, _O_BINARY | _O_CREAT | _O_TRUNC | _O_RDWR, S_IREAD | S_IWRITE );
		if (hOut < 0)
		{
			if (!BatchRun)
				Message("Cannot create binary file for in PerformBackwardEliminationGDM()", "ERROR");
			_close(hIn);
			if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );
			if (lpTmpFile) delete[] lpTmpFile;
			if (lpElimFile) delete[] lpElimFile;
			if (pResponse) delete[] pResponse;
			if (pWeights) delete[] pWeights;
			if (pCopyCol) delete[] pCopyCol;
			return(false);
		}

		// copy the intercept
		_lseek(hIn, 0L, SEEK_SET);
		_read(hIn, pCopyCol, NumRows * sizeof(double));
		_write(hOut, pCopyCol, NumRows * sizeof(double));
		int NumTmpCols = 1;
				
		// skip past the euclidean distance columns to drop the euclidean distance predictor
		int nSplines2Drop = GetEuclideanSplines(pParams);
		_lseek(hIn, (long)(nSplines2Drop * NumRows * sizeof(double)), SEEK_CUR);

		// now write the rest of the selected predictor data
		int Cols2Copy = NumCols - nSplines2Drop - 1;
		for ( int i=0; i<Cols2Copy; i++ )
		{
			_read(hIn, pCopyCol, NumRows * sizeof(double));
			_write(hOut, pCopyCol, NumRows * sizeof(double));
			++NumTmpCols;
		}

		_close(hIn);
		_close(hOut);

		if (!DoIncrementalGDM(pParams, BatchRun, true, 
					          lpElimFile, myRunName, 
							  pResponse, pWeights, 
							  NumRows, NumTmpCols, 0,
							  fptr, uptr))
		{
			if (!BatchRun)
				Message("Cannot DoIncrementalGDM in PerformBackwardEliminationGDM()", "ERROR");

			if (lpTmpFile) delete[] lpTmpFile;
			if (lpElimFile) delete[] lpElimFile;
			if (pResponse) delete[] pResponse;
			if (pWeights) delete[] pWeights;
			if (pCopyCol) delete[] pCopyCol;
			return(false);
		}
	} // if ( (UseEuclidean(pParams)) && (0 == GetEuclideanPredDroppedState(pParams)))


	//
	// now do the non-dropped predictors
	//
	for ( int ThisPred=1; ThisPred<=nPreds; ThisPred++ )
	{
		if ( -1 == GetPredictorDroppedAt(pParams, ThisPred) )
		{
			continue;
		}


		// if this predictor is a candidate to be dropped...
		if ( 0 < GetPredictorInUseAt(pParams, ThisPred) ) 
		{
			// change colour to processing (orange)
			if (uptr)
				uptr(-5, -5, -5, -5, ThisPred);

			// open the binary file for the FULL GDM 
			hIn = _open( lpTmpFile, _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE );
			if (hIn < 0)
			{
				if (!BatchRun)
					Message("Cannot create open binary file for in PerformBackwardEliminationGDM()", "ERROR");
				if (lpTmpFile) delete[] lpTmpFile;
				if (lpElimFile) delete[] lpElimFile;
				if (pResponse) delete[] pResponse;
				if (pWeights) delete[] pWeights;
				if (pCopyCol) delete[] pCopyCol;
				return(false);
			}


			// create new binary matrix using intercept column followed by predictor columns
			hOut = _open( lpElimFile, _O_BINARY | _O_CREAT | _O_TRUNC | _O_RDWR, S_IREAD | S_IWRITE );
			if (hOut < 0)
			{
				if (!BatchRun)
					Message("Cannot create binary file for in PerformBackwardEliminationGDM()", "ERROR");
				_close(hIn);
				if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );
				if (lpTmpFile) delete[] lpTmpFile;
				if (lpElimFile) delete[] lpElimFile;
				if (pResponse) delete[] pResponse;
				if (pWeights) delete[] pWeights;
				if (pCopyCol) delete[] pCopyCol;
				return(false);
			}

			// copy the intercept
			_lseek(hIn, 0L, SEEK_SET);
			_read(hIn, pCopyCol, NumRows * sizeof(double));
			_write(hOut, pCopyCol, NumRows * sizeof(double));

			int NumTmpCols = 1;

			// if we have geographic distance as a predictor
			if ( (UseEuclidean(pParams)) && (0 == GetEuclideanPredDroppedState(pParams)))
			{
				int nTmpSplines = GetEuclideanSplines(pParams);
				for (int i=0; i<nTmpSplines; i++ )
				{
					_read(hIn, pCopyCol, NumRows * sizeof(double));
					_write(hOut, pCopyCol, NumRows * sizeof(double));
				}
				NumTmpCols += nTmpSplines; // increment the columns count in the regression matrix
			}

			// include all the other 'non-dropped' predictors
			for ( int j=1; j<=nPreds; j++ )
			{
				int nTmpSplines = GetPredictorSplinesAt(pParams, j);

				// ignore if the predictor is NOT 'in-use'
				if ( 0 == GetPredictorInUseAt(pParams, j) )
				{
					continue;
				}

				// don't use the dropped predictor
				if ( ThisPred == j ) 
				{
					_lseek(hIn, (long)(nTmpSplines * NumRows * sizeof(double)), SEEK_CUR);
					continue;
				}

				// else copy the relevent columns to the temporary binary matrix
				for ( int k=0; k<nTmpSplines; k++ )
				{
					_read(hIn, pCopyCol, NumRows * sizeof(double));
					_write(hOut, pCopyCol, NumRows * sizeof(double));
				}
				NumTmpCols += nTmpSplines; // increment the columns count in the regression matrix
			} // for ( int j=1; j<=nPreds; j++ )

			_close(hIn);
			_close(hOut);

			if (!DoIncrementalGDM(pParams, BatchRun, true, 
							      lpElimFile, myRunName, 
								  pResponse, pWeights, 
							      NumRows, NumTmpCols, ThisPred,
							      fptr, uptr))
		    {
				if (!BatchRun)
					Message("Cannot DoIncrementalGDM in PerformBackwardEliminationGDM()", "ERROR");

				if (lpTmpFile) delete[] lpTmpFile;
				if (lpElimFile) delete[] lpElimFile;
				if (pResponse) delete[] pResponse;
				if (pWeights) delete[] pWeights;
				if (pCopyCol) delete[] pCopyCol;
				return(false);
			} // if (!DoIncrementalGDM(...
		} // if ( 0 < GetPredictorInUseAt(pParams, ThisPred) )
	} // for ( int ThisPred=1; ThisPred<=nPreds; ThisPred++ )


	//
	// Cleanup
	//
	if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );
	if ( ( _access( lpElimFile, 0 ) ) != -1 ) remove( lpElimFile );
	if (lpTmpFile) delete[] lpTmpFile;
	if (lpElimFile) delete[] lpElimFile;
	if (pResponse) delete[] pResponse;
	if (pWeights) delete[] pWeights;
	if (pCopyCol) delete[] pCopyCol;
	return(true);
}




//
// Run the core GDM for a Backward Elimination
//
bool DoBaseGDM(char *pParams, bool BatchRun, bool DeleteBinaryFile, 
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
		if (!BatchRun)
			Message("DataTablePath does NOT exist in GDM_FitFromParamFile", "ERROR");
		if (DataTablePath) delete[] DataTablePath;
		return(false);
	}

	
	//
	// Get a memory image of the composite data file
	//
	GDM_INT nTableRows = 0;
	GDM_INT nTableCols = 0;
	double *pTableData = GetCompositeTableIntoMemory(DataTablePath, &nTableRows, &nTableCols, BatchRun, fptr);
	if (NULL == pTableData)
	{
		if (!BatchRun)
			Message("Cannot GetCompositeTableIntoMemory() in GDM_FitFromParamFile", "ERROR");
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
							   &nMatrixCols, BatchRun, fptr))
	{
		if (!BatchRun)
			Message("Cannot CreatePredictorBinary() in GDM_FitFromParamFile", "ERROR");
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		if (pTableData) delete[] pTableData;
		return(false);
	}
	//PrintBinary2CSV(lpTmpFile, pParams, theSplineCounts, nTableRows, nMatrixCols, fptr);

	//
	// Extract the response column and the weights column before freeing the table data memory block
	//
	*pResponse = new double [nTableRows];
	memcpy(*pResponse, &pTableData[COL_RESPONSE * nTableRows], nTableRows * sizeof(double));
	if (NULL == *pResponse)
	{
		if (!BatchRun)
			Message("Cannot extract Response column in GDM_FitFromParamFile", "ERROR");
		if (pTableData) delete[] pTableData;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}

	*pWeights = new double [nTableRows];
	memcpy(*pWeights, &pTableData[COL_WEIGHTS * nTableRows], nTableRows * sizeof(double));
	if (NULL == *pWeights)
	{
		if (!BatchRun)
			Message("Cannot extract Weights column in GDM_FitFromParamFile", "ERROR");
		if (*pResponse) delete[] *pResponse;
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
	double *pPredictorMatrix = new double [nTableRows * nMatrixCols];
	if (NULL == pPredictorMatrix)
	{
		if (!BatchRun)
			Message("Cannot allocate pPredictorMatrix in GDM_FitFromParamFile", "ERROR");
		if (*pResponse) delete[] *pResponse;
		if (*pWeights) delete[] *pWeights;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}


	//
	// Populate the Predictor Matrix from the binary image created in CreatePredictorBinary()
	//
	int h = _open( lpTmpFile, _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE );
	if (h < 0)
	{
		if (!BatchRun)
			Message("Cannot open binary file for READ in GDM_FitFromParamFile", "ERROR");
		if (pPredictorMatrix) delete[] pPredictorMatrix;
		if (*pResponse) delete[] *pResponse;
		if (*pWeights) delete[] *pWeights;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}
	//_read(h, pPredictorMatrix, nTableRows * nMatrixCols * sizeof( double ));
	double *pTmp = pPredictorMatrix;
	for ( int i=0; i<nMatrixCols; i++ )
	{
		_read(h, pTmp, (unsigned)nTableRows * sizeof( double ));
		pTmp += nTableRows;
	}
	_close(h);


	//
	// Do the matrix regression
	//
	fptr("Performing GDM regression...", 0);
	double dGDMDeviance;
	double *pCoefficients = WeightedNNLSRegression( lpTmpFile, 
		                                            pPredictorMatrix, 
											        nTableRows, 
											        nMatrixCols, 
											        *pResponse, 
											        &dGDMDeviance, 
											        *pWeights,
													fptr);
	if ( NULL == pCoefficients )
	{
		if (!BatchRun)
			Message("pCoefficients are NULL", "ERROR in GDM_FitFromParamFile");
		if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );
		if (pPredictorMatrix) delete[] pPredictorMatrix;
		if (*pResponse) delete[] *pResponse;
		if (*pWeights) delete[] *pWeights;
		if (theSplineCounts) delete[] theSplineCounts; 
		if (theQuantiles) delete[] theQuantiles;
		return(false);
	}


	//
	// remove the temporary matrix file if function parameters define removal
	//
	if (DeleteBinaryFile)
	{
		if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );
	}
	

	//
	// create a NULL model and return the deviance
	//
	double dNullDeviance = GetWeightedNULLDeviance( nTableRows, *pResponse, *pWeights );


	//
	// calculate deviance explained as a percentage
	//
	double dDevianceExplained = ( 1.0 - ( dGDMDeviance / dNullDeviance ) ) * 100;


	//
	// write relevent outputs back to the parameter file
	//	
	fptr("Updating parameter file...", 0);
	//PrintGDMResults(pParams, dDevianceExplained, nMatrixCols, pCoefficients, false);
	UpdateParameterFile(pParams, dNullDeviance, dGDMDeviance, dDevianceExplained,
						pCoefficients[0], &pCoefficients[1], (int)nMatrixCols-1, theSplineCounts, (int)nTableRows, fptr);

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
// Do an incremental GDM for the Backward Elimination Process
//
bool DoIncrementalGDM(char *pParams, bool BatchRun, bool DeleteBinaryFile, 
	                  char *lpTmpFile, char *myRunName,
	                  double *pResponse, double *pWeights, 
			          int NumRows, int NumCols, int Index,
			          FPTR fptr, UPTR uptr)
{
	double *pPredictorMatrix = new double [NumRows * NumCols];
	if (NULL == pPredictorMatrix)
	{
		if (!BatchRun)
			Message("Cannot allocate pPredictorMatrix in DoIncrementalGDM()", "ERROR");
		return(false);
	}


	//
	// Populate the Predictor Matrix from the binary image created in CreatePredictorBinary()
	//
	int h = _open( lpTmpFile, _O_BINARY | _O_RDONLY, S_IREAD | S_IWRITE );
	if (h < 0)
	{
		if (!BatchRun)
			Message("Cannot open binary file for READ in DoIncrementalGDM()", "ERROR");
		if (pPredictorMatrix) delete[] pPredictorMatrix;
		return(false);
	}
	//_read(h, pPredictorMatrix, NumRows * NumCols * sizeof( double ));
	double *pTmp = pPredictorMatrix;
	for ( int i=0; i<NumCols; i++ )
	{
		_read(h, pTmp, NumRows * sizeof( double ));
		pTmp += NumRows;
	}
	_close(h);


	//
	// Do the matrix regression
	//
	fptr("Performing GDM regression...", 0);
	double dGDMDeviance;
	double *pCoefficients = WeightedNNLSRegression( lpTmpFile, 
		                                            pPredictorMatrix, 
											        NumRows, 
											        NumCols, 
											        pResponse, 
											        &dGDMDeviance, 
											        pWeights,
													fptr);
	if ( NULL == pCoefficients )
	{
		if (!BatchRun)
			Message("pCoefficients are NULL", "ERROR in GDM_FitFromParamFile");
		if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );
		if (pPredictorMatrix) delete[] pPredictorMatrix;
		return(false);
	}


	//
	// remove the temporary matrix file if function parameters define removal
	//
	if (DeleteBinaryFile)
	{
		if ( ( _access( lpTmpFile, 0 ) ) != -1 ) remove( lpTmpFile );
	}


	//
	// create a NULL model and return the deviance
	//
	double dNullDeviance = GetWeightedNULLDeviance( NumRows, pResponse, pWeights );


	//
	// calculate deviance explained as a percentage
	//
	double dDevianceExplained = ( 1.0 - ( dGDMDeviance / dNullDeviance ) ) * 100;


	//
	// get the sum of coefficients for the dropped predictor
	// 
	double dSum = 0.0;
	for ( int i=1; i<NumCols; i++ )
	{
		dSum += pCoefficients[i];
	}

	
	//
	// update parameter file for this GDM predictor's elimination run
	//
	char myKey[64];
	char pName[64];
	sprintf( myKey, "DroppedPred_%d", Index );
	if (Index == 0)
		strcpy(pName, "Geographic Distance");
	else
		GetPredictorNameAt(pParams, pName, Index);	
	SetProfileString( myRunName, myKey, pName, pParams );
	//
	sprintf( myKey, "Intercept_%d", Index );
	SetProfileDouble( myRunName, myKey, pCoefficients[0], pParams );
	//
	sprintf( myKey, "NULLDeviance_%d", Index );
	SetProfileDouble( myRunName, myKey, dNullDeviance, pParams );
	//
	sprintf( myKey, "GDMDeviance_%d", Index );
	SetProfileDouble( myRunName, myKey, dGDMDeviance, pParams );
	//
	sprintf( myKey, "DevianceExplained_%d", Index );
	SetProfileDouble( myRunName, myKey, dDevianceExplained, pParams );
	//
	double dEliminated = GetProfileDouble("GDMODEL", "DevExplained", pParams) - dDevianceExplained;
	sprintf( myKey, "PredictorExplainedDeviance_%d", Index );
	SetProfileDouble( myRunName, myKey, dEliminated, pParams );
	//
	sprintf( myKey, "SumOfCoeffs_%d", Index);
	SetProfileDouble( myRunName, myKey, dSum, pParams );
	//
	sprintf( myKey, "NumPreds_%d", Index);
	SetProfileInt( myRunName, myKey, GetNumPredictors(pParams) - GetNumDroppedPreds(pParams), pParams );
	//
	int nPreds = GetNumPredictors(pParams);
	for ( int qq=0; qq<=nPreds; qq++ )
	{
		sprintf( myKey, "TestPred_%d.%d", Index, qq);

		if ( qq == Index )
			SetProfileInt( myRunName, myKey, 0, pParams );
		else
		{
			if (qq == 0)
			{
				if (0 == GetEuclideanPredDroppedState(pParams))
					SetProfileInt( myRunName, myKey, 1, pParams );
				else
					SetProfileInt( myRunName, myKey, 0, pParams );
			}
			else
			{
				SetProfileInt( myRunName, myKey, GetPredictorInUseAt(pParams, qq), pParams );
			}
		}
	}
	AppendRow(pParams);

	//
	// update the list box with the difference in deviance explained
	//
	if (uptr)
		uptr(dDevianceExplained, dEliminated, GetSumOfCoefficients(pParams), dSum, Index);

	// change colour to done (green)
	if (uptr)
		uptr(-10, -10, -10, -10, Index);

	//
	// Cleanup
	//
	if (pPredictorMatrix) delete[] pPredictorMatrix;
	return(true);
}



