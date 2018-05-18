//
// pcoMatrix.h
//
#ifndef MATRIX

#define MATRIX
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS	0
#endif
#ifndef EXIT_FAILURE
#define EXIT_FAILURE	1
#endif

#ifndef SCALAR
	#define SCALAR double
#endif
#ifndef	EPS
	#define	EPS		1E-6
#endif
#ifndef TINY
	#define TINY	1E-20
#endif
#ifndef MAX_ITER
	#define MAX_ITER 100
#endif
#define forall(i)	for(i = 0; i < n; i++)
#define rotate(a, i, j, k, l) {		 \
	double x = a[i][j], y = a[k][l]; \
	a[i][j] = x * c - y * s;		 \
	a[k][l] = x * s + y * c;		}

typedef SCALAR *vector, **matrix;
typedef int *ivector, **imatrix;


vector newvec(int n);

matrix newmat(int nrow, int ncol);

vector new_vector(int n);

matrix new_matrix(int nrow, int ncol);

void init_vector(vector v, int size);

void init_matrix(matrix a, int ncol);

void free_vector(vector v);

void free_matrix(matrix a);

double innerproduct(int n, vector u, vector v);

void vec_print(vector v, int n, int perline, char *format);

void vec_fprint(FILE *fp, vector v, int n, int perline, char *format);

int vec_fscan(FILE *fp, vector v, int n);

void mat_print(matrix a, int ncol, int perline, char *format);

void mat_fprint(FILE *fp, matrix a, int ncol, int perline, char *format);

void mat_fprint_with_header(FILE *fp, matrix a, int ncol, int perline, char *format);

int mat_fscan(FILE *fp, matrix a, int ncol);

double gjmatinv(int n, matrix a);

double house(int n, vector x);

void tridiagonalize(int n, matrix a, vector d, vector e);

int eigen(int n, matrix a, vector d, vector e);

int jacobi(int n, matrix a, vector e, matrix w);

imatrix newimat(int nrow, int ncol);

imatrix new_imatrix(int nrow, int ncol);

void free_imatrix(imatrix a);

ivector newivec(int n);

ivector new_ivector(int n);

void free_ivector(ivector v);

#endif

