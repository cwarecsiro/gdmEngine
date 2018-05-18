//
// pcoMatrix.cpp
//
#include "stdafx.h"
#include "pcoMatrix.h"
#include "Message.h"

vector newvec(int n)
{
	return (vector)malloc(sizeof(SCALAR) * n);
}

matrix newmat(int nrow, int ncol)
{
	int i;
	matrix a;
	
	a = (matrix)malloc((nrow+1) * sizeof(void *));
	if(a == NULL) return NULL;
	for(i = 0; i < nrow; i++) {
		a[i] = (vector)malloc(sizeof(SCALAR) * ncol);
		if(a[i] == NULL) {
			while(--i >= 0) free(a[i]);
			free(a); return NULL;
		}
	}
	a[nrow] = NULL;
	return a;
}

vector new_vector(int n)
{
	vector v;
	
	v = newvec(n);
	if(v == NULL) {
		fprintf(stderr, "insufficent memory available\n");
		exit(EXIT_FAILURE);
	}
	return v;
}

matrix new_matrix(int nrow, int ncol)
{
	matrix a;
	
	a = newmat(nrow, ncol);
	if(a == NULL) {
		fprintf(stderr, "insufficient memory available\n");
		exit(EXIT_FAILURE);
	}
	return a;
}

void init_vector(vector v, int size)
{
	int i;
	
	for(i = 0; i < size; i++)
	{
		v[i] = 0;
	}
}

void init_matrix(matrix a, int ncol)
{
	int i, j;
	
	for(i = 0; a[i] != NULL; i++)
	{
		for(j = 0; j < ncol; j++)
		{
			a[i][j] = 0;
		}
	}
}


void free_vector(vector v)
{
	free(v);
}

void free_matrix(matrix a)
{
	matrix b;
	
	b = a;
	while(*b != NULL) free(*b++);
	free(a);
}

double innerproduct(int n, vector u, vector v)
{
	int i, n5;
	double s;
	
	s = 0;
	n5 = n % 5;
	for(i = 0; i < n5; i++) s += u[i]*v[i];
	for(i = n5; i < n; i += 5)
		s += u[i]*v[i] + u[i+1]*v[i+1] + u[i+2]*v[i+2]
		               + u[i+3]*v[i+3] + u[i+4]*v[i+4];
	return s;
}

void vec_print(vector v, int n, int perline, char *format)
{
	int j, k;
	
	k = 0;
	for(j = 0; j < n; j++) {
		printf(format, v[j]);
		if(++k >= perline) {
			k = 0;
			printf("\n");
		}
	}
	if(k != 0) printf("\n");
}

void vec_fprint(FILE *fp, vector v, int n, int perline, char *format)
{
	int j, k;
	
	k = 0;
	for(j = 0; j < n; j++) {
		fprintf(fp, format, v[j]);
		if(++k >= perline) {
			k = 0;
			fprintf(fp, "\n");
		}
	}
	if(k != 0) fprintf(fp, "\n");
}

int vec_fscan(FILE *fp, vector v, int n)
{
	int i;
	double temp;
	
	for(i = 0; i < n; i++) {
		fscanf(fp, "%lf", &temp);
		if(ferror(fp) || feof(fp))
			return EXIT_FAILURE;
		v[i] = (SCALAR)temp;
	}
	return EXIT_SUCCESS;
}

void mat_print(matrix a, int ncol, int perline, char *format)
{
	int i;
	
	for(i = 0; a[i] != NULL; i++) {
		vec_print(a[i], ncol, perline, format);
		if(ncol > perline) printf("\n");
	}
}

void mat_fprint(FILE *fp, matrix a, int ncol, int perline, char *format)
{
	int i;
	
	for(i = 0; a[i] != NULL; i++) {
		vec_fprint(fp, a[i], ncol, perline, format);
		if(ncol > perline) fprintf(fp, "\n");
	}
}

void mat_fprint_with_header(FILE *fp, matrix a, int ncol, int perline, char *format)
{
	int i;
	
	for(i = 0; a[i] != NULL; i++) {}
		fprintf(fp, "%d %d\n\n", i, ncol);

	for(i = 0; a[i] != NULL; i++) {
		vec_fprint(fp, a[i], ncol, perline, format);
		if(ncol > perline) fprintf(fp, "\n");
	}
}

int mat_fscan(FILE *fp, matrix a, int ncol)
{
	int i;
	int	row_size, col_size;
	
	fscanf(fp, "%d %d", &row_size, &col_size);
	
	for(i = 0; a[i] != NULL; i++) {}
	
	if((i != row_size) || (ncol != col_size))
		return EXIT_FAILURE;
	
	for(i = 0; a[i] != NULL; i++)
		if(vec_fscan(fp, a[i], ncol))
			return EXIT_FAILURE;
	
	return EXIT_SUCCESS;
}

double gjmatinv(int n, matrix a)
{
	int i, j, k;
	double t, u, det;
	
	det = 1;
	for(k = 0; k < n; k++) {
		t = a[k][k];
		det *= t;
		for(i = 0; i < n; i++) a[k][i] /= t;
		a[k][k] = 1 / t;
		for(j = 0; j < n; j++)
			if(j != k) {
				u = a[j][k];
				for(i = 0; i < n; i++)
					if(i != k) a[j][i] -= a[k][i] * u;
					else       a[j][i] = -u / t;
			}
	}
	return det;
}

double house(int n, vector x)
{
	int i;
	double s, t;
	
	s = sqrt(innerproduct(n, x, x));
	
	if(s != 0) {
		if(x[0] < 0) s = -s;
		x[0] += s;
		t = 1 / sqrt(x[0] * s);
		for(i = 0; i < n; i++)
		{
			x[i] *= t;
		}
	}
	return -s;
}


void tridiagonalize(int n, matrix a, vector d, vector e)
{
	int i, j, k;
	double s, t, p, q;
	vector v, w;
	
	for(k = 0; k < n - 2; k++)
	{
		v = a[k];
		d[k] = v[k];
		e[k] = house(n - k - 1, &v[k + 1]);
		if(e[k] == 0) continue;
		
		for(i = k + 1; i < n; i++)
		{
			s = 0;
			for(j = k + 1; j < i; j++)
			{
				s += a[j][i] * v[j];
			}
			for(j = i; j< n; j++)
			{
				s += a[i][j] * v[j];
			}
			d[i] = s;
		}
		
		t = innerproduct(n-k-1, &v[k + 1], &d[k + 1]) / 2;
		
		for(i = n - 1; i > k; i--)
		{
			p = v[i];
			q = d[i] - t * p;
			d[i] = q;
			for(j = i; j < n; j++)
			{
				a[i][j] -= p * d[j] + q * v[j];
			}
		}
	}
	if(n >= 2) {
		d[n - 2] = a[n - 2][n - 2];
		e[n - 2] = a[n - 2][n - 1];
	}
	if(n >= 1) {
		d[n - 1] = a[n - 1][n - 1];
	}
	for(k = n - 1; k >= 0; k--)
	{
		v = a[k];
		if(k < n - 2) {
			for(i = k + 1; i < n; i++)
			{
				w = a[i];
				t = innerproduct(n-k-1, &v[k + 1], &w[k + 1]);
				for(j = k + 1; j < n; j++)
				{
					w[j] -= t * v[j];
				}
			}
		}
		for(i = 0; i < n; i++)
		{
			v[i] = 0;
		}
		v[k] = 1;
	}
}

int eigen(int n, matrix a, vector d, vector e)
{
	int i, j, k, h, iter;
	double c, s, t, w, x, y;
	vector v;

	tridiagonalize(n, a, d, &e[1]);

	e[0] = 0;
	for(h = n - 1; h > 0; h--) {
		j = h;
		while(fabs(e[j]) > EPS * (fabs(d[j - 1]) + fabs(d[j])))
			j--;
		if(j == h) continue;
		iter = 0;
		do {
			if(++iter > MAX_ITER) return EXIT_FAILURE;
			w = (d[h - 1] - d[h]) / 2;
			t = e[h] * e[h];
			s = sqrt(w * w + t);
			if(w < 0)
				s = -s;
			x = d[j] - d[h] + t / (w + s);
			y = e[j + 1];
			for(k = j; k < h; k++)
			{
				if(fabs(x) >= fabs(y)) {
					t = -y / x;
					c = 1 / sqrt(t * t + 1);
					s = t * c;
				}
				else {
					t = -x / y;
					s = 1 / sqrt(t * t + 1);
					c = t * s;
				}
				w = d[k] - d[k + 1];
				t = (w * s + 2 * c * e[k + 1]) * s;
				d[k] -= t;
				d[k + 1] += t;
				if(k > j)
					e[k] = c * e[k] - s * y;
				e[k + 1] += s * (c * w - 2 * s * e[k + 1]);
				
				for(i = 0; i < n; i++)
				{
					x = a[k][i];
					y = a[k + 1][i];
					a[k][i] = c * x - s * y;
					a[k + 1][i] = s * x + c * y;
				}
				if(k < h - 1) {
					x = e[k + 1];
					y = -s * e[k + 2];
					e[k + 2] *= c;
				}
			}
		}while(fabs(e[h]) > EPS * (fabs(d[h - 1]) + fabs(d[h])));
	
	}


	for(k = 0; k < n - 1; k++)
	{
		h = k;
		t = d[h];
		
		for(i = k + 1; i < n; i++)
		{
			if(d[i] > t) {
				h = i;
				t = d[h];
			}
		}
		d[h] = d[k];
		d[k] = t;
		v = a[h];
		a[h] = a[k];
		a[k] = v;
	}

	return EXIT_SUCCESS;
}

int jacobi(int n, matrix a, vector e, matrix w)
{
	int i, j, k, iter;
	double t, c, s, tolerance, offdiag;
	vector v;
	
	s = offdiag = 0;
	forall(j) {
		forall(k) w[j][k] = 0;
		w[j][j] = 1; s += a[j][j] * a[j][j];
		for(k = j + 1; k < n; k++)
			offdiag += a[j][k] * a[j][k];
	}
	tolerance = EPS * EPS * (s / 2 + offdiag);
	for(iter = 1; iter <= MAX_ITER; iter++) {
		offdiag = 0;
		for(j = 0; j < n - 1; j++)
			for(k = j + 1; k < n; k++)
				offdiag += a[j][k] * a[j][k];
		if(offdiag < tolerance) break;
		for(j = 0; j < n - 1; j++) {
			for(k = j + 1; k < n; k++) {
				if(fabs(a[j][k]) < TINY) continue;
				t = (a[k][k] - a[j][j]) / (2 * a[j][k]);
				if(t >= 0)  t = 1 / (t + sqrt(t * t + 1));
				else		t = 1 / (t - sqrt(t * t + 1));
				c = 1 / sqrt(t * t + 1);
				s = t * c;	t *= a[j][k];
				a[j][j] -= t; 	a[k][k] += t;	a[j][k] = 0;
				for(i = 0; i < j; i++)		rotate(a, i, j, i, k)
				for(i = j + 1; i < k; i++)	rotate(a, j, i, i, k)
				for(i = k + 1; i < n; i++)	rotate(a, j, i, k, i)
				forall(i)					rotate(w, j, i, k, i)
			}
		}
	}

	if(iter > MAX_ITER) return EXIT_FAILURE;

	forall(i) {
		e[i] = a[i][i];
		forall(j)
			a[i][j] = w[i][j];
	}
	
	for(i = 0; i < n - 1; i++) {
		k = i; t = e[i];
		for(j = i + 1; j < n; j++)
			if(e[j] > t) { k = j; t = e[k]; }
		e[k] = e[i]; e[i] = t;
		v = a[k]; a[k] = a[i]; a[i] = v;
	}
	return EXIT_SUCCESS;
}


imatrix newimat(int nrow, int ncol)
{
	int i;
	imatrix a;
	
	a = (imatrix)malloc((nrow+1) * sizeof(void *));
	if(a == NULL) return NULL;
	for(i = 0; i < nrow; i++) {
		a[i] = (ivector)malloc(sizeof(int) * ncol);
		if(a[i] == NULL) {
			while(--i >= 0) free(a[i]);
			free(a); return NULL;
		}
	}
	a[nrow] = NULL;
	return a;
}

imatrix new_imatrix(int nrow, int ncol)
{
	imatrix a;
	
	a = newimat(nrow, ncol);
	if(a == NULL) {
		/* fprintf(stderr, "insufficient memory available\n"); */
		exit(EXIT_FAILURE);
	}
	return a;
}

void free_imatrix(imatrix a)
{
	imatrix b;
	
	b = a;
	while(*b != NULL) free(*b++);
	free(a);
}

ivector newivec(int n)
{
	return (ivector)malloc(sizeof(int) * n);
}

ivector new_ivector(int n)
{
	ivector v;
	
	v = newivec(n);
	if(v == NULL) {
		exit(EXIT_FAILURE);
	}
	return v;
}

void free_ivector(ivector v)
{
	free(v);
}
