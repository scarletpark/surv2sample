#include <R.h>
#include <Rmath.h>
#include "cif_misc.h"

void twosample_incidence_f01(int *event, int *y1, int *y2, double *s1, double *s2,
	double *f01, int *n)
{
	int i;
	f01[0] = ((double) (event[0]==1)) / ((double) y1[0]/1. + (double) y2[0]/1.);
	for (i=1; i<*n; i++)
		f01[i] = f01[i-1] + ((double) (event[i]==1)) / ((double) y1[i]/s1[i-1] + (double) y2[i]/s2[i-1]);
}

double incr2sym(double *a, int n, int i, int j)
{
	/* computes the increment of a symmetric function of 2 variables (e.g., covariance fnc) */
	/* a is n by n, returns a[i,j]-a[i-1,j]-a[i,j-1]+a[i-1,j-1] */
	/* assumes a defined for j<=i (for j>i a[i,j] may be undefined, not accessed) */
	/* (if j>i they're swapped) */
	/* handles edges as if a[i,j]=0 for i<0 or j<0 */

	if (j>i) {
		int ii = i;
		i = j;
		j = ii;
	}
	if (j==0) {
		if (i==0) {
			return a[0+0*n]; /* a[0,0] */
		} else {
			return a[i+0*n]-a[(i-1)+0*n]; /* a[i,0]-a[i-1,0] */
		}
	} else { /* a[i,j]-a[i-1,j]-a[i,j-1]+a[i-1,j-1] */
		if (i==j) { /* diagonal */
			return a[i+j*n]-2.*a[i+(i-1)*n]+a[(i-1)+(i-1)*n]; /* symmetry */
		} else {
			return a[i+j*n]-a[(i-1)+j*n]-a[i+(j-1)*n]+a[(i-1)+(j-1)*n];
		}
	}
}

double incr1sym(double *a, int n, int j)
{
	/* computes the increment a[n-1,j]-a[n-1,j-1] */
	/* a is n by n, for j>i a[i,j] is may be undefined, it's symmetric */
	/* handles a[n-1,-1] as 0 */
	if (j==0) {
		return a[(n-1)+0*n];
	} else {
		return a[(n-1)+j*n]-a[(n-1)+(j-1)*n];
	}
}

void fj1_covariance(double *rho, int *n, int grp, int *event, int *group,
	double *fj1, double *fj2, int *yj)
{
	/* computes the covariance function for the limit of \hat F_{j1} - F_{j1} for the jth group (j=grp) */
	/* see formula (3) of Lin (1997), beware of misprints */
	/* rho[i,j] computed only for j<=i (symmetric) */
	int i,j;
	double temp1,temp2;
	double a,b,c,d,e; /* values of 5 integrals in rho at time min(t_i,t_j) which is t_j here (j<=i) */
	a = b = c = d = e = 0.;
	for (j=0; j<*n; j++) {
		temp1 = ((yj[j]>0) ? 1./(yj[j]*yj[j]) : 0.) * ((group[j]==grp)*(event[j]==1));
		temp2 = ((yj[j]>0) ? 1./(yj[j]*yj[j]) : 0.) * ((group[j]==grp)*(event[j]==2));
		a += (1.-fj2[j])*(1.-fj2[j])*temp1;
		b += fj1[j]*fj1[j]*temp2;
		c += temp1+temp2;
		d += (1.-fj2[j])*temp1;
		e += fj1[j]*temp2;
		for (i=j; i<*n; i++) {
			rho[i+j**n] = a + b + fj1[i]*fj1[j]*c - (fj1[i]+fj1[j])*(d+e);
		}
	}
}
