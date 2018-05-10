//
// NNLSDoubleCore.h
//
#ifndef __NNLSDOUBLECORE_H__
#define __NNLSDOUBLECORE_H__

#include "stdafx.h"

#if defined _M_X64

int nnls_c( double* a, 
            GDM_INT* mda, 
			GDM_INT* m, 
			GDM_INT* n, 
			double* b, 
			double* x, 
			double* rnorm, 
			double* w, 
			double* zz, 
			GDM_INT* index, 
			int* mode );


int nnls_Weighted( double* a, 
				   GDM_INT* mda, 
				   GDM_INT* m, 
				   GDM_INT* n, 
				   double* b, 
				   double *pWghts,
	               double* x, 
				   double* rnorm, 
				   double* w, 
				   double* zz, 
				   GDM_INT* index, 
	               int* mode );


#elif defined _WIN32

int nnls_c( double* a, 
            int* mda, 
			int* m, 
			int* n, 
			double* b, 
			double* x, 
			double* rnorm, 
			double* w, 
			double* zz, 
			int* index, 
			int* mode );


int nnls_Weighted( double* a, 
				   int* mda, 
				   int* m, 
				   int* n, 
				   double* b, 
				   double *pWghts,
	               double* x, 
				   double* rnorm, 
				   double* w, 
				   double* zz, 
				   int* index, 
	               int* mode );

#endif

#endif // __NNLSDOUBLECORE_H__