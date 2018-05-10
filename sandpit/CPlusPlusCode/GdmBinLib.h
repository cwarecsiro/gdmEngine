//
// GdmBinLib.h
//
#ifndef	__GDMBINLIB_H__
#define	__GDMBINLIB_H__

#include "stdafx.h"
#include "myCallback.h"

extern "C" 
{

_declspec(dllexport) void Test( char *ParcelPath, 
	                            char *DataPath, 
								char *InputPath, 
								char *OutputPath, 
								FPTR fptr );

_declspec(dllexport) bool RunGi2Csv( char *ParcelPath, 
	                            char *DataPath, 
								char *InputPath, 
								char *OutputPath, 
								bool InitialStep,
								FPTR fptr );



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These are specifically called from the R-Stats package for the PCT modelling project
//
_declspec(dllexport) bool FilterPCTsFromR(char **ppPCTInputTable, char **ppWorkspace, int *pNumPCTs, int *ppPCTs);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These are specifically called from the R-Stats package
//
_declspec(dllexport) void SaveGDMParams(char **WDPath, char **paramFileName, char **tablePath, int *NumPreds, bool UseGeoDist);

_declspec(dllexport) bool ExtractAndUpdateQuantilesSigTest(char **ppParams);

_declspec(dllexport) bool DoSigTestGDM(char **ppParams, int *pIterations);




 _declspec(dllexport) void GDM_FitFromTable( char **wspath, double *pData, 
	                                         int *pDoGeo, int *pPreds, 
										     int *pRows, int *pCols, 
										     int *pSplines, double *pQuantiles,
										     double *pGDMDev, double *pNullDev, double *pExpDev, 
										     double *pIntercept, double *pCoeffs,
										     double *pY, double *pX, double *pE );



 _declspec(dllexport) void GDM_PredictFromTable(double *pData, 
	                                           int *pDoGeo, int *pPreds, int *pRows, 
				                               double *pQuantiles,  int *pSplines, double *pCoeffs,
				                               double *pX);



_declspec(dllexport)void GetPredictorPlotData( double *pPredData, int *pLength,
					                           double *pCoefficients,
					                           double *pQuantiles,
					                           int *pSplines );

//
// These are specifically called from the .NET interface
//
_declspec(dllexport) void TestDLLBitSize();

_declspec(dllexport) bool GDM_FitFromResiduals( char *pParams, bool BatchRun, FPTR fptr );

_declspec(dllexport) bool GDM_FitFromParamFile( char *pParams, bool BatchRun, FPTR fptr );

_declspec(dllexport) bool GDM_PredictFromParamFile( char *pModelParams, char *pInTablePath, char *pOutPath, FPTR fptr );

_declspec(dllexport) bool ExtractAndUpdateQuantiles( char *pParams, bool BatchRun, FPTR fptr);

_declspec(dllexport) bool ExtractAndUpdateAllQuantiles( char *pParams, bool BatchRun, FPTR fptr);

_declspec(dllexport) bool GetGridMinMax(char *GridPath, double *pMin, double *pMax, FPTR fptr);

_declspec(dllexport) bool FilterTableThruGrids(char *pParams, FPTR fptr);

_declspec(dllexport) bool CreateCompositeFromComposite(char *InputPath, char *OutputPath, FPTR fptr);

_declspec(dllexport) bool CreateMaskedCompositeFromComposite(char *InputPath, 
	                                                         char *MaskPath, 
										                     char *OutputPath, 
										                     bool DoAggregation,
										                     FPTR fptr);

_declspec(dllexport) bool CreateCompositeFromSitePlusSpecies( char *InputPath, 
	                                                          char *OutputPath, 
															  bool HaveAbundance, 
															  int WeightType,
															  FPTR fptr);

_declspec(dllexport) bool CreateCompositeFromMaskedSitePlusSpecies( char *InputPath, 
	                                                                char *MaskPath,
	                                                                char *OutputPath, 
											                        bool DoAggregate,
										                            bool HaveAbundance, 
										                            int WeightType,
										                            FPTR fptr);

_declspec(dllexport) bool CreateCompositeFromSiteBySpecies( char *InputPath, 
	                                                        char *OutputPath, 
															int Cutpoint, 
															int WeightType,
															FPTR fptr);

_declspec(dllexport) bool CreateCompositeFromSiteBySpeciesAbundance( char *InputPath, 
	                                                                 char *OutputPath, 
																	 int WeightType,
																	 FPTR fptr);

_declspec(dllexport) bool MaskSiteBySpeciesTable(char *InputPath, 
	                                             char *OutputPath, 
							                     char *MaskPath, 
							                     bool DoAggregation, 
							                     FPTR fptr);

_declspec(dllexport) bool SxSDataHasFloatingPoint( char *InputPath );

_declspec(dllexport) bool ExtractPredictors( char *pParams, bool DoBatch, FPTR fptr);

_declspec(dllexport) bool PerformBatchedBackwardEliminationGDM(char *pParams, FPTR fptr, UPTR uptr);

_declspec(dllexport) bool PerformBackwardEliminationGDM(char *pParams, bool BatchRun, FPTR fptr, UPTR uptr);

_declspec(dllexport) int GDM_CountRowsInTable(char *pPath, FPTR fptr);

_declspec(dllexport) void GDM_GetGeographicMinMax(char *pPath, double *pMin, double *pMax, FPTR fptr);

_declspec(dllexport) void GDM_GetPredictorMinMax(char *pPath, int nIndex, int nPreds, double *pMin, double *pMax, FPTR fptr);

_declspec(dllexport) bool TransformPredictors( char *pParams, bool DoBatch, FPTR fptr);

_declspec(dllexport) bool CreateTransformsFromR(char **ppInputGridTable,
	                                            char **ppOutputGridsWorkDir,
	                                            int *Extrap_Method,
	                                            bool *DoGeo, 
	                                            int *NumPredictors, 
	                                            int *pSplines,
	                                            double *pQuantiles, 
	                                            double *pCoeffs);

_declspec(dllexport) bool GetTransformMinMax( char *sPath, double *pMin, double *pMax, FPTR fptr);

_declspec(dllexport) bool DoGDMClassification(char *pParams, 
						                      char *lpDomainPath,
						                      char *lpOutName,
						                      char *lpOutPath,
   					                          int NumClasses, 
											  int NumSamples,
											  bool DoBatch,
											  char *TimeString,
						                      FPTR fptr );

_declspec(dllexport) bool DomainIsValidSubGridForClassification(char *SubDomainPath, char *TransformPath);


_declspec(dllexport) bool CreateClassTrainingTable(char *pParams, 
	                                               char *lpDomainPath, 
	                                               char *lpPredDomainPath, 
												   char *lpOutPath, 
												   int NumSamples, 
												   FPTR fptr );


_declspec(dllexport) bool ClassifyFromTable(char *pParams, 
                                            char *TablePath, 
                                            int  nSamples, 
                                            int  nPreds, 
                                            char *lpDomainPath, 
                                            char *lpOutName, 
                                            FPTR fptr);


_declspec(dllexport) bool ClassifyFromCompositeTableData( char *ParamPath,
	                                                      int nSamples, 
	                                                      int nRows,
                                                          int nPreds, 
	                                                      char *InputTablePath, 
                                                          char *OutputTablePath, 
                                                          FPTR fptr);



_declspec(dllexport) bool CreateSpatiallyBalancedSamplesFromGDMClassification(char *pParams, int NumSamples, char *pOutPath, FPTR fptr);



_declspec(dllexport) bool CreateJRTrainingFile(char *pParams, char *JRInputXY, char *JROutputPath, FPTR fptr);



_declspec(dllexport) bool DoUnconstrainedBinaryProbTableExternal(char *pParams,
                                                                 int JParam,
                                                                 double RParam,
											                     int MaxClassIndex,
                                                                 char *TrainingDataPath,
											                     char *SiteDataPath,
                                                                 char *OutTablePath,
											                     char *TimeString,
                                                                 FPTR fptr);


_declspec(dllexport) bool DoUnconstrainedBinaryProbGridsExternal(char *pParams,
	                                                             bool UserDefinedParams,
                                                                 int JParam,
                                                                 double RParam,
																 int MaxClassIndex,
                                                                 char *TrainingDataPath,
                                                                 char *OutDir,
                                                                 bool IgnoreEmptyClasses,
																 bool DoNormalisation,
																 char *TimeString,
                                                                 FPTR fptr);


_declspec(dllexport) bool DoMaskedUnconstrainedBinaryProbGridsExternal(char *pParams,
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
                                                                       FPTR fptr);


//
// Get the maximum one-based class index from the training data table
//
_declspec(dllexport) bool GetMaxClassIndexFromTrainingData(int *pMaxClassIndex, char *JRInputPath, FPTR fptr);


//
// Master Test Function for SOS calculations
//
_declspec(dllexport) bool GetBestJRValsViaParamFile(char *ParamFilePath, FPTR fptr);

_declspec(dllexport) bool GetBestJRParams(int *pJParam, double *pRParam, int *pMaxClassIndex,
	                                      char *JRInputPath, char *JROutputPath, FPTR fptr);



//
// Create standard GDM Density Grid
//
_declspec(dllexport) bool DoStandardGDMDensity(char *pParams, 
						                       char *lpDomainPath,
						                       char *lpOutPath,
						                       char *lpSamplePath,
						                       int NumSamples,
						                       bool UseNearest,
						                       float fCutPoint,
						                       bool DoSOS,
											   char *TimeStamp,
						                       FPTR fptr);


//
// Create GDM Density Grid (called externally)
//
_declspec(dllexport) bool DoGDMDensitySim(char *Params, 
					                      char *DomainPath,
					                      char *OutName,
					                      int NumSamples,
 				                          FPTR fptr);


//
// Test that data cells in the domain and test grid are ALL coindident
//
_declspec(dllexport) bool TestEsriGrids(char *pDomainPath, char *TestGridPath, char *pLogFilePath, 
				                        int nRows, int nCols, 
				                        float fDomainNoData, float fGridNoData, 
										double dCellSize, double dXMin, double dYMax,
				                        FPTR fptr);


//
// Test that the two grids have the same extent and CellSize
//
_declspec(dllexport) int TestExtentAndCellsize(char *pRef, char *pTest);


//
// Determine the maximum size block that can be allocated
//
_declspec(dllexport) double GetMaxMemBlock(FPTR fptr);


//
// Filter composite table thru binary export grids called from .NET batch dialog
//
_declspec(dllexport) bool FilterCompositeTables(char *pParams, bool DoBatch, FPTR fptr);


//
//  Creates a time series set of unconstrained probability grids using an external JRTraining File
//
_declspec(dllexport) bool DoTimeSeriesProbGridsExternal(char *pStartParams,
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
                                                        FPTR fptr);


//
//  Creates a time series set of unconstrained probability grids using an 
//  external JRTraining File but applies a mask.
//
_declspec(dllexport) bool DoMaskedTimeSeriesProbGridsExternal(char *pStartParams,
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
                                                              FPTR fptr);


//
//  Creates a time series set of GDM Unsupervised Classification grids using an 
//  external JRTraining File but applies a mask.
//
_declspec(dllexport) bool DoMaskedTimeSeriesGDMClassificationExternal(char *pStartParams,
	                                                                  char *pEndParams,
                                                                      char *TrainingDataPath,
                                                                      char *MaskPath,
                                                                      char *OutDir,
												                      char *OutName,
										                              int nStartYear,
                                                                      int nEndYear,
                                                                      int nIncrement,
										                              char *TimeString,
                                                                      FPTR fptr);


//
//  Creates a dissimilarity grid from two time-step gdm models (both GDMs MUST have the SAME transform predictors)
//
_declspec(dllexport) bool DoTimeSeriesDissimilarityGrid(char *pStartParams,
	                                                    char *pEndParams,
                                                        char *OutPath,
										                FPTR fptr);

_declspec(dllexport) bool DoTimeSeriesDissimilarityGridAsBlocks(char *pStartParams,
	                                                            char *pEndParams,
                                                                char *OutPath,
														        int nBlockSize,
										                        FPTR fptr);

_declspec(dllexport) bool DoTimeSeriesImpactTable(char *pDemandPointPath,
	                                              char *pStartParams,
	                                              char *pEndParams,
                                                  char *OutPath,
							                      FPTR fptr);


_declspec(dllexport) bool CreateImpactSimilarityTable(char *pDemandPointPath, char *OutPath, FPTR fptr);


_declspec(dllexport) bool CreateImpactDissimilarityTable(char *pDemandPointPath, char *OutPath, FPTR fptr);


_declspec(dllexport) bool ClassifyGDMViaKMeans(char *pDemandPointPath, int nToSkip, char *pGDMParams, char *OutPath, FPTR fptr);

_declspec(dllexport) bool ClassifyGDMViaKMeansBanded(char *pDemandPointPath, 
	                                                 int nToSkip, 
													 char *pGDMParams, 
													 char *OutPath, 
													 int nLower,
													 int nUpper,
													 FPTR fptr);


_declspec(dllexport) bool DoTimeSeriesImpactGrids(char *pDemandPointPath,
	                                              char *pStartParams,
	                                              char *pEndParams,
                                                  char *OutPathStart,
							                      char *OutPathEnd,
							                      FPTR fptr);

_declspec(dllexport) bool DoTimeSeriesImpactGridsANN(char *pDemandPointPath,
	                                                 char *pStartParams,
	                                                 char *pEndParams,
                                                     char *OutPathStart,
							                         char *OutPathEnd,
							                         FPTR fptr);

_declspec(dllexport) bool DoANNTest(char *pDemandPointPath, int nToSkip);

_declspec(dllexport) bool DoTimeSeriesImpactTable_Version_1(char *pDemandPointPath, char *pStartParams, FPTR fptr);

_declspec(dllexport) bool DoTimeSeriesImpactTable_Version_2(char *pDemandPointPath, char *pStartParams, FPTR fptr);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This group is for use in the BFT Composite JPEG creator...
//
_declspec(dllexport) bool RGBGridsOK(char *RGrid, char *GGrid, char *BGrid, FPTR fptr);

_declspec(dllexport) bool GetGridMetrics(char *pPath, double *pMin, double *pMax, FPTR fptr);

_declspec(dllexport) bool CreateCompositeFile(char *RedPath, double dRedMin, double dRedMax,
                                              char *GreenPath, double dGreenMin, double dGreenMax,
                                              char *BluePath, double dBlueMin, double dBlueMax, 
                                              char *CompPath, FPTR fptr);

_declspec(dllexport) bool GetRowColCount(char *pPath, int *pRows, int *pCols);
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


_declspec(dllexport) bool CreateReport(char *ModelPath, 
	                                   char *IbraPath, 
									   char *LLSPath,
									   char *OutPath, 
									   char *IbraLookupPath,
									   int LLSIndex, 
									   bool CountCells, 
									   FPTR fptr);


_declspec(dllexport) bool CreateRVCSimilarityTable(char *pParams, FPTR fptr);

_declspec(dllexport) bool DoRVClassification(char *pParams, FPTR fptr);

_declspec(dllexport) bool DoANNClassificationViaParamFile(char *ParamFilePath, FPTR fptr);

_declspec(dllexport) bool CreateRasterCentroids(double dCellSize, double dLeftEdge, double dTopEdge, int nNumRows, int nNumCols, char *pOutPath, FPTR fptr);

}
#endif // __GDMBINLIB_H__


