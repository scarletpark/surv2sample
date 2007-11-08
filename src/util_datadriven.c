#include <R.h>
#include <Rmath.h>
#include "util_neyman.h"
#include "util_datadriven.h"

/* DATA-DRIVEN STATISTICS with nested subsets */

void datadriven_stats_nested(double *score_d, double *sigma_d, int *d, int *d0,
	double *stats, double *stats_penal, int *selected,
	double *stat_d, double *stat_d0, double *stat_selected,
	double *sigmawork, double *scorework, double *penalty, double *choltol)
{
	/* computes score statistics for nested subsets with dimension d0 to d, */
	/* penalised statistics, their maximum and maximiser */
	/* penalty = log(n) for Schwarz's BIC */
	int i,j,k;
	
	*selected = 0;
	for (k=*d0; k<=*d; k++) {
		/* copy the k by k submatrix of sigma_d to sigmawork */
		for (i=0; i<k; i++)
			for (j=i; j<k; j++)
				sigmawork[i+j*k] = sigma_d[i+j**d];
		quadrstat(score_d, sigmawork, &k, stats+k-*d0, scorework, choltol);
		stats_penal[k-*d0] = stats[k-*d0] - (double) k**penalty;
		if (stats_penal[k-*d0] > stats_penal[*selected]) {
			*selected = k-*d0;
		}
	}
	*stat_d = stats[*d-*d0];
	*stat_d0 = stats[0];
	*stat_selected = stats[*selected];
}


/* TWO-TERM APPROXIMATIONS */
/* of the distribution of the data-driven statistic */
/* with nested subsets with minimum dimension 1 */
/* see Kraus (2007, Lifetime Data Anal., eq. (12)) */

/* chi-square approximations */

void h_approx_nested(double *stat_bic, int *n, double *pval_bic_h)
{
	double logn,logn2;
	double *p_logn,*p_logn2;
	logn = log((double) *n);
	p_logn = &logn;
	logn2 = 2.*logn;
	p_logn2 = &logn2;
	
	if (*stat_bic <= logn) {
		*pval_bic_h = 1. - h1_approx_nested(stat_bic,n);
	} else {
		if (*stat_bic >= logn2) {
			*pval_bic_h = 1. - h2_approx_nested(stat_bic,n);
		} else {
			*pval_bic_h = 1. - ( (logn2-*stat_bic)/(logn)*h1_approx_nested(p_logn,n) + (*stat_bic-logn)/(logn)*h2_approx_nested(p_logn2,n) );
		}
	}
}

double h1_approx_nested(double *x, int *n)
{
	return (2.*pnorm(sqrt(*x),0.,1.,1,0)-1.)*(2.*pnorm(sqrt(log((double) *n)),0.,1.,1,0)-1.);
}
double h2_approx_nested(double *x, int *n)
{
	return (2.*pnorm(sqrt(*x),0.,1.,1,0)-1.)*(2.*pnorm(sqrt(log((double) *n)),0.,1.,1,0)-1.) + 2.*pnorm(sqrt(log((double) *n)),0.,1.,0,0);
}

/* approximations with F-distributions instead of chi-square (currently not used) */

void h_approx_nested_f(double *stat_bic, int *n, double *pval_bic_h_f)
{
	double logn,logn2;
	double *p_logn,*p_logn2;
	logn = log((double) *n);
	p_logn = &logn;
	logn2 = 2.*logn;
	p_logn2 = &logn2;
	
	if (*stat_bic <= logn) {
		*pval_bic_h_f = 1. - h1_approx_nested_f(stat_bic,n);
	} else {
		if (*stat_bic >= logn2) {
			*pval_bic_h_f = 1. - h2_approx_nested_f(stat_bic,n);
		} else {
			*pval_bic_h_f = 1. - ( (logn2-*stat_bic)/(logn)*h1_approx_nested_f(p_logn,n) + (*stat_bic-logn)/(logn)*h2_approx_nested_f(p_logn2,n) );
		}
	}
}

double h1_approx_nested_f(double *x, int *n)
{
	return (2.*pt(sqrt(*x),*n-1.,1,0)-1.)*(2.*pt(sqrt(log((double) *n)),*n-1.,1,0)-1.);
}
double h2_approx_nested_f(double *x, int *n)
{
	return (2.*pt(sqrt(*x),*n-1.,1,0)-1.)*(2.*pt(sqrt(log((double) *n)),*n-1.,1,0)-1.) + 2.*pt(sqrt(log((double) *n)),*n-1.,0,0);
}

