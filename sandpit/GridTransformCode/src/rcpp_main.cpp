#include <Rcpp.h>
using namespace Rcpp;
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


// [[Rcpp::export]] 
int rcpp_CreateTransformsFromR(std::string ppInputGridTable,
                               std::string ppOutputGridsWorkDir,
                               int Extrap_Method,
                               bool DoGeo, 
                               int NumPredictors,
                               int* pSplines,
                               double* pQuantiles, 
                               double* pCoeffs) 
  {
  int nPreds = NumPredictors;
  int nDoGeo = DoGeo;
  //int nPreds = *NumPredictors;
  //int nDoGeo = *DoGeo;
  
  char buff0[BUFF1024];
  std::sprintf(buff0, "%s\\GDM_Params.txt", ppOutputGridsWorkDir);
  //std::sprintf(buff0, "%s\\GDM_Params.txt", *ppOutputGridsWorkDir);
  FILE *fpOut = fopen(buff0, "w+t");
  std::fprintf(fpOut, "###\n### Parameter file create from R to transform grids via the GDM dll\n###\n\n");
  
  std::fprintf(fpOut, "[GDMODEL]\n");
  std::fprintf(fpOut, "UseEuclidean=%d\n", nDoGeo);
  std::fprintf(fpOut, "ExtrapolationMethod=%d\n", Extrap_Method);
  std::fprintf(fpOut, "WorkspacePath=%s\n", ppOutputGridsWorkDir);
  //std::fprintf(fpOut, "ExtrapolationMethod=%d\n", *Extrap_Method);
  //std::fprintf(fpOut, "WorkspacePath=%s\n", *ppOutputGridsWorkDir);
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
  FILE *fpIn = fopen(ppInputGridTable, "r+t");
  //FILE *fpIn = fopen(*ppInputGridTable, "r+t");
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
    std::fprintf(fpOut, "%s\\gdmXtran\n", ppOutputGridsWorkDir);
    std::fprintf(fpOut, "%s\\gdmYtran\n", ppOutputGridsWorkDir);
    //std::fprintf(fpOut, "%s\\gdmXtran\n", *ppOutputGridsWorkDir);
    //std::fprintf(fpOut, "%s\\gdmYtran\n", *ppOutputGridsWorkDir);
  }
  
  gmpath gmp;
  char buff2[BUFF1024];
  fpIn = fopen(ppInputGridTable, "r+t");
  //fpIn = fopen(*ppInputGridTable, "r+t");
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
      std::fprintf(fpOut, "PredTran%d=%s\\%s", i + 1, ppOutputGridsWorkDir, gmp.GetName(buff2));
      //std::fprintf(fpOut, "PredTran%d=%s\\%s", i + 1, *ppOutputGridsWorkDir, gmp.GetName(buff2));
    }
  }
  fclose(fpIn);
  fclose(fpOut);
  
  
  //
  // Do the Transform
  //
  std::sprintf(buff0, "%s\\GDM_Params.txt", ppOutputGridsWorkDir);
  //std::sprintf(buff0, "%s\\GDM_Params.txt", *ppOutputGridsWorkDir);	
  //Message(buff0, "buff0");
  TransformPredictors(buff0, false, NULL);
  
  return(1);
  } // end rcpp_CreateTransformsFromR


