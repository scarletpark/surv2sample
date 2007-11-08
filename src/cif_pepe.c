#include <R.h>
#include <Rmath.h>
#include "cif_misc.h"
#include "cif_ks.h"
#include "cif_pepe.h"

/* TWO-SAMPLE TEST FOR EQUALITY OF CUMULATIVE INCIDENCE FUNCTIONS */
/* Pepe (1991, JASA), Bajorunaite & Klein (2007, CSDA) */

void twosample_incidence_pepe(double *time, int *event, int *group, int *n, double *tau,
	int *nsim, double *f11, double *f12, double *f21, double *f22, double *stat,
	double *pval_asympt, double *pval_sim)
{
	int i,j;
	int *y1,*y2;
	y1 = (int *) R_alloc(2**n,sizeof(int));
	y2 = y1 + *n;
	double temp,var,stat_sim;
	double *ks_process,*s1,*s2,*f01,*g;
	ks_process = (double *) R_alloc(5**n,sizeof(double));
	s1 = ks_process + *n;
	s2 = s1 + *n;
	f01 = s2 + *n;
	g = f01 + *n;
	
	/* compute the test process f21-f11 */
	twosample_incidence_ks_process(event, group, n, y1, y2, f11, f12, f21, f22,
		s1, s2, ks_process);
	/* compute the Pepe test statistic int of f21(t)-f11(t) dt from 0 to tau */
	*stat = int_stat_pepe(ks_process, time, n, tau);
	/* compute the null estimate of f01 */
	twosample_incidence_f01(event, y1, y2, s1, s2, f01, n);

	twosample_incidence_pepe_variance(time,event,group,n,tau,1,f01,f12,y1,&temp);
	var = temp;
	twosample_incidence_pepe_variance(time,event,group,n,tau,2,f01,f22,y2,&temp);
	var += temp;
	
	/*
	twosample_incidence_pepe_variance(time,event,group,n,tau,1,f11,f12,y1,&temp);
	var = temp;
	twosample_incidence_pepe_variance(time,event,group,n,tau,2,f21,f22,y2,&temp);
	var += temp;
	*/
	
	temp = fabs(*stat);
	*stat /= sqrt(var);
	*pval_asympt = 2.*pnorm(fabs(*stat),0.,1.,0,0);
	
	/* LWY simulation */
	if (*nsim > 0) {
		GetRNGstate();
		*pval_sim = 0.;
		for (j=0; j<*nsim; j++) {
			for (i=0; i<*n; i++)
				g[i] = norm_rand();
			twosample_incidence_lwy_process(event, group, n, y1, y2, f01, f12, f01, f22,
				g, ks_process);
			stat_sim = fabs(int_stat_pepe(ks_process, time, n, tau));
			*pval_sim += (double) (stat_sim > temp);
		}
		*pval_sim /= *nsim;
		PutRNGstate();
	}
}

double twosample_incidence_pepe_variance(double *time, int *event, int *group, int *n,
	double *tau, int grp, double *fj1, double *fj2, int *yj, double *var)
{
	/* martingale-based variance estimator of Pepe's statistic (Bajorunaite & Klein, 2007, CSDA) */
	int i;
	double intfj1 = 0.;
	double t1 = *tau;
	double temp,temp1;
	
	*var = 0.;
	for (i=*n-1; i>=0; i--) {
		intfj1 += fj1[i-1]*(t1-time[i]);
		t1 = time[i];
		temp = (yj[i]>0 ? 1./yj[i] : 0.);
		temp1 = (*tau-time[i])*(1.-fj2[i])*temp - temp*intfj1;
		*var += temp1*temp1*((group[i]==grp)*(event[i]==1));
		temp1 = (*tau-time[i])*fj1[i]*temp - temp*intfj1;
		*var += temp1*temp1*((group[i]==grp)*(event[i]==2));
	}
}

double int_stat_pepe(double *f, double *time, int *n, double *tau)
{
	/* Pepe's (1991) integral statistic */
	/* computes \int_0^\tau f(t)dt where f(t) is piecewise constant with jumps at time */
	/* must be tau >= time[n-1] */
	int i;
	double y = 0.;
	
	for (i=0; i<*n-1; i++)
		y += f[i]*(time[i+1]-time[i]);
	y += f[*n-1]*(*tau-time[*n-1]);
	return y;
}
