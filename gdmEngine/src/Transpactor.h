//---------------------------------------------------------------------------

#ifndef CompyH
#define CompyH

/*****************************************************************************
 * TRANSPACTOR.H
 * Object to manage compacting of grids, and transformation
 *
 * Written October 2017 by Tom Harwood.
 *
 ***************************************************************************/
class TTranspactor
 {
 public:
   // constructor
   TTranspactor(void);
   // constructor2
   TTranspactor(const int arg_i_extrapolation);
   // destructor
   ~TTranspactor(void);
//   int TransCompactBinaryGrids(string arg_s_res_path);
   int TransBinaryGridsToFloat(void);//(string arg_s_res_path);       // note names of transgrids are in s_binary_out
   void Refresh(void);
//   int LoadGlobalGDMParameters(string arg_s_path);
   void SetNLayers(const int arg_n_layers){n_layers=arg_n_layers;}
   void SetNRows(const int arg_n_rows){n_rows=arg_n_rows;}
   void SetNCols(const int arg_n_cols){n_cols=arg_n_cols;}
   void SetExtrapolation(const int arg_i_extrapolation){i_extrapolation=arg_i_extrapolation;}
   void SetGeo(const int arg_i_geo){i_geo=arg_i_geo;}
   
   void SetSplines( const int arg_i_layer,           // layer to set splines for (default just use 0)
                    const int arg_n_splines,          // number of splines
                    const double *arg_d_coeffs,            //ptr to array of coeffs[arg_n_splines]
                    const double *arg_d_quants);            // ptr to array of quants[arg_n_splines]
   void SetUpExtrapolation(const int arg_i_layer); // set up linear parameters
   double GetExtraSplineValue(const double arg_d_value,     // value to transform
                              const int arg_i_layer);       // layer to use 0 by default
   double GetExtrapolation(const double arg_d_value,     // value to transform
                           const int arg_i_layer);       // layer to use 0 by default
   double GetSplineValue(const double arg_d_value,     // value to transform
                         const int arg_i_layer);       // layer to use 0 by default
   double GetSplineValueGeo(const double arg_d_value);
   // Backward compatible version
   double CalculateSplineValue(const double arg_d_position,     // value to transform
                               const int arg_n_splines,          // number of splines
                               const double *arg_d_pCoeffs,            //ptr to array of coeffs[arg_n_splines]
                               const double *arg_d_pQuants);            // ptr to array of quants[arg_n_splines]
   //Backward compatible version
   double DoSplineCalc(const double arg_d_val,        // value to transform
                       const double arg_d_q1,     // three spline quantiles
                       const double arg_d_q2,
                       const double arg_d_q3);
//	void SetZero(const int arg_index){b_iszero[arg_index]=1;}
//	int IsZero(const int arg_index){return b_iszero[arg_index];}
//	void SetBare(const int arg_index){b_isbare[arg_index]=1;}
//	void SetBareScale(const double arg_d_scale){d_bare_scale=arg_d_scale;}
//	double GetBareValue(const double arg_d_value)     // value to transform
//	  {return (arg_d_value*d_bare_scale);}
//	void SetUpGlobFilePaths(const string arg_s_clim_path,
//                         const string arg_s_land_path,
//                         const string arg_s_subs_path,
//                         const string arg_s_resid_path);
//	void SetUpGlobOutFLTFilePaths(const string arg_s_res_path);
	void SetModelParameters(int arg_geo,
                         int arg_n_layers,
                         int* arg_n_splines,                         
                         float* arg_knots,
                         float* arg_coefficients); 
	void SetGridPaths(int arg_geo,
                   int arg_n_layers,
                   std::string arg_InFilePaths[],
                   std::string arg_OutFilePaths[]);
	
protected:

 private:
   int n_layers;
   long int n_rows, n_cols;
   int i_geo;
   string s_binary_in[50]; // full pathnames for all input grids (IN ALPHABETICAL ORDER of grid name)
   string s_binary_out[5]; // full pathnames for output transgrids as FLT
//   int b_isbare[30]; // flag special transform for binary BARE ground
//   int b_iszero[30]; // is this a zero grid for this biome? Check before doing anything
//   double d_bare_scale; // simple scaling if bare
   int n_splines[50];
   int i_extrapolation; //0 don't extrapolate 1 do capped minimum extrapolation 2 10% extrapolation
   int i_transtype[50]; //Type of transform for each layer: 0 set all to zero, 1 GDM transform, 2 Binary scaling
   int b_active[50]; //flag for active spline layers in the list.. e.g. constant grids, repeat grids
   double d_quant[50][10];  //Assume nobody uses >5 splines for overfitting reasons
   double d_coeff[50][10];  //Assume nobody uses >5 splines for overfitting reasons
   double d_m[50][2];  // slope parameter y=mx+c     [layer][bot or top]
   double d_c[50][2];  // constant parameter y=mx+c   [layer][bot or top]
   int i_geo_splines; // number of splines for geo
   double d_geo_quant[10]; // quantiles for geo
   double d_geo_coeff[10]; // coefficients for geo
   string s_geo_binary_out[2]; // full pathnames for geo output transgrids as FLT
   }; // end class TTranspactor
//---------------------------------------------------------------------------
#endif
