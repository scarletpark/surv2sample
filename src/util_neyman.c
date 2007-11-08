#include <R.h>
#include <Rmath.h>
#include "my_cholesky.h"
#include "util_neyman.h"

/* QUADRATIC SCORE TEST STATISTIC */

void quadrstat(double *score, double *sigma, int *d, double *stat, double *work, double *choltol)
{
	/* computes the quadratic statistic score^T * sigma^{-1} * score by cholesky factorisation */
	/* score and upper triangle of sigma left unchanged, work destroyed, diag of sigma destroyed */
	int cholflag,j;
	
	for (j=0; j<*d; j++)
		work[j] = score[j];
	cholflag = my_cholesky2(sigma,d,choltol);
	/* sigma may have the lower triangle undefined, my_cholesky2 copies upper to lower */
	
	my_chsolve2(sigma,d,work); /* work is overwritten, score left unchanged */
	*stat = 0.;
	for (j=0; j<*d; j++)
		*stat += work[j]*score[j];
}

/* BASES OF SMOOTH FUNCTIONS used in Neyman's tests */

void cosine_basis(double *x, int n, double *p, int deg, int degzero,
	int dummy_shifted, int dummy_normalized)
{
	/* dummy variables in order to have the same formal parameters as legendre; because of pointers to these fncs */
	int i,j,k;
	int ncolp = ((degzero==0)?(deg):(deg+1)); /* number of columns of p */
	
	for (i=0; i<n; i++) {
		if (degzero==0) {
			for (j=0; j<ncolp; j++) { /* j is the column in p */
				p[i+j*n] = M_SQRT2 * cos(x[i]*(j+1.)*M_PI);
			}
		} else {
			p[i+0*n] = 1.;
			for (j=1; j<ncolp; j++) { /* j is the column in p */
				p[i+j*n] = M_SQRT2 * cos(x[i]*j*M_PI);
			}
		}
	}
}

void cosine1_basis(double *x, int n, double *p, int deg, int degzero,
	int dummy_shifted, int dummy_normalized)
{
	/* cos((j-1/2)*pi*u) */
	
	/* dummy variables in order to have the same formal parameters as legendre; because of pointers to these fncs */
	int i,j,k;
	int ncolp = ((degzero==0)?(deg):(deg+1)); /* number of columns of p */
	
	for (i=0; i<n; i++) {
		if (degzero==0) {
			for (j=0; j<ncolp; j++) { /* j is the column in p */
				p[i+j*n] = M_SQRT2 * cos(x[i]*(j+.5)*M_PI);
			}
		} else {
			p[i+0*n] = 1.;
			for (j=1; j<ncolp; j++) { /* j is the column in p */
				p[i+j*n] = M_SQRT2 * cos(x[i]*(j-.5)*M_PI);
			}
		}
	}
}

void legendre_basis(double *x, int n, double *p, int deg, int degzero,
	int shifted, int normalized)
{
	/*
	x		points where polynomials are to be evaluated, vector of length n
	n		length of x
	p		values of polynomials, matrix (n by deg or n by deg+1 depending on degzero); output, must be allocated
	deg		polynomials of degree up to deg will be computed
	degzero		should p contain also the polynomial of degree 0? 0=no, 1=yes
	shifted		should we compute shifted polynomials on [0,1]? 0=no (then x is in [-1,1]), 1=yes (x in [0,1])
	normalized	should polynomials be normalized? 0=no, 1=yes
	*/
	
	int i,j,k;
	int ncolp = ((degzero==0)?(deg):(deg+1)); /* number of columns of p */
	double u;
	double p_2,p_1; /* value of P_{j-1} and P_{j-1} */
	double *norms; /* norms of polynomials */
	
	if (normalized != 0) {
		norms = (double *) R_alloc(ncolp,sizeof(double));
		if (shifted != 0) {
			if (degzero !=0 ) {
				for (j=0; j<ncolp; j++) norms[j] = sqrt(1./(2.*j+1.));
			} else {
				for (j=0; j<ncolp; j++) norms[j] = sqrt(1./(2.*(j+1.)+1.));
			}
		} else {
			if (degzero !=0 ) {
				for (j=0; j<ncolp; j++) norms[j] = sqrt(2./(2.*j+1.));
			} else {
				for (j=0; j<ncolp; j++) norms[j] = sqrt(2./(2.*(j+1.)+1.));
			}
		}
	}
	
	for (i=0; i<n; i++) {
		u = ( (shifted == 0) ? x[i] : 2.*x[i]-1. );
		j = 0;
		if (degzero != 0) {
			p[i+j*n] = 1.;
			++j;
		}
		p_1 = 1.;
		p_2 = 0.;
		k = 1; /* k is the degree of the polynomial being just computed */
		for ( ; j<ncolp; j++) { /* j is the column in p */
			p[i+j*n] = (2.*k-1.)/k*u*p_1 - (k-1.)/k*p_2;
			p_2 = p_1;
			p_1 = p[i+j*n];
			++k;
		}
		
		if (normalized !=0) {
			for (j=0; j<ncolp; j++) p[i+j*n] /= norms[j];
		}
		
	}
}
