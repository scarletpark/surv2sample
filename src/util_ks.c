#include <R.h>
#include <Rmath.h>
#include "util_ks.h"

double ks_stat_incr(double *dx, int *n)
{
	/* computes the KS statistic from increments of the test process */
	/* i.e., max(abs(cumsum(dx))) */
	double y = 0., x = 0.;
	int i;
	for (i=0; i<*n; i++) {
		x += dx[i];
		if (fabs(x)>y) {
			y = fabs(x);
		}
	}
	return y;
}

double ks_stat_cum(double *x, int *n)
{
	/* computes the KS statistic from a cumulative process */
	/* x is a (cumulative) test process, i.e., returns max|x| */
	double y = 0.;
	int i;
	for (i=0; i<*n-1; i++) {
		if (fabs(x[i])>y) {
			y = fabs(x[i]);
		}
	}
	return y;
}

double psupw(double x, int kmax)
{
	/* computes c.d.f. of sup_{[0,1]}|W| (B is a Brownian motion) (Billingsley, 1968) */
	/* kmax = where to cut the series for the c.d.f. */
	unsigned int k;
	double temp = -1. + 2.*pnorm(x,0.,1.,1,0);
	
	for (k=1; k<=kmax; k++) {
		temp += ((k & 1) ? (-2.) : 2.) * (pnorm((2.*k+1.)*x,0.,1.,1,0)-pnorm((2.*k-1.)*x,0.,1.,1,0)); /* plus if k even, minus if odd */
	}
	return temp;
}

double psupb(double x, double a, int kmax)
{
	/* computes c.d.f. of sup_{[0,a]}|B| (B is a Brownian bridge) (Hall & Wellner, 1980) */
	/* kmax = where to cut the series for the c.d.f. */
	unsigned int k;
	double temp;
	
	if (a<1.) {
		temp = -1. + 2.*pnorm(x*pow(a*(1.-a),-.5),0.,1.,1,0);
		for (k=1; k<=kmax; k++) {
			temp += ((k & 1) ? (-2.) : 2.) * exp(-2.*k*k*x*x) * ( pnorm(x*sqrt((1.-a)/a)*(2.*k+1./(1.-a)),0.,1.,1,0) - pnorm(x*sqrt((1.-a)/a)*(2.*k-1./(1.-a)),0.,1.,1,0) );
		}
	} else {
		temp = 1.;
		for (k=1; k<=kmax; k++) {
			temp += ((k & 1) ? (-2.) : 2.) * exp(-2.*k*k*x*x);
		}
	}
	
	return temp;
}

