#include <R.h>
#include <Rmath.h>
#include "util_ks.h"
#include "cif_misc.h"
#include "cif_ks.h"

/* KOLMOGOROV--SMIRNOV TEST FOR COMPARING CUMULATIVE INCIDENCE FUNCTIONS IN TWO SAMPLES */
/* (Lin, 1997, Statistics in  Medicine) */
/* the test compares cif's for cause 1, i.e., H0: f11=f21 */

void twosample_incidence_ks(int *event, int *group, int *n, int *nsim,
	double *f11, double *f12, double *f21, double *f22,	double *test_process,
	double *stat, double *pval_sim, double *test_process_plot_sim, int *nsim_plot)
{
	GetRNGstate();
	
	int i,j,n1,n2;
	double temp;
	int *y1,*y2;
	y1 = (int *) R_alloc(2**n,sizeof(int));
	y2 = y1 + *n;

	double *s1,*s2,*test_process_sim,*g,*f01;
	s1 = (double *) R_alloc(5**n,sizeof(double));
	s2 = s1 + *n;
	test_process_sim = s2 + *n;
	g = test_process_sim + *n;
	f01 = g + *n;

	double stat_sim;
	
	/* OBSERVED test statistic */
	twosample_incidence_ks_process(event, group, n, y1, y2, f11, f12, f21, f22,
		s1, s2, test_process);
	*stat = ks_stat_cum(test_process, n);
	
	n1 = y1[0];
	n2 = y2[0];
	
	/* null (pooled sample) estimator of the cause 1 cif */
	twosample_incidence_f01(event, y1, y2, s1, s2, f01, n);
	
	/* LWY SIMULATION */
	/* simulated processes are computed with the pooled sample estimator f01 */
	/* (simulations with pooled sample f01 lead to conservative test, */
	/* recommended by Bajorunatite & Klein (2007, CSDA); */
	/* individual f11, f21 give anticonserv. approximation) */
	if (*nsim>0) {
		*pval_sim = 0.;
		/* *pval_sim_indiv = 0.; */
		for (j=0; j<*nsim_plot; j++) { /* always must be *nsim>=*nsim_plot */
			for (i=0; i<*n; i++) {
				g[i] = norm_rand();
			}
			/* compute the resampled test process with pooled sample f01 */
			twosample_incidence_lwy_process(event, group, n, y1, y2, f01, f12, f01, f22,
				g, test_process_sim);
			stat_sim = ks_stat_cum(test_process_sim, n);
			*pval_sim += (double) (stat_sim > *stat);
			/* save the simulated process for plotting */
			for (i=0; i<*n; i++)
				test_process_plot_sim[i+j**n] = test_process_sim[i];
			/* resampling with individual f11,f21; not used */
			/*
			twosample_incidence_lwy_process(event, group, n, y1, y2, f11, f12, f21, f22,
				g, test_process_sim);
			stat_sim = ks_stat_cum(test_process_sim, n, i1, i2);
			*pval_sim_indiv += (double) (stat_sim > *stat);
			*/
		}
		for (j=*nsim_plot; j<*nsim; j++) {
			for (i=0; i<*n; i++) {
				g[i] = norm_rand();
			}
			/* compute the resampled test process with pooled sample f01 */
			twosample_incidence_lwy_process(event, group, n, y1, y2, f01, f12, f01, f22,
				g, test_process_sim);
			stat_sim = ks_stat_cum(test_process_sim, n);
			*pval_sim += (double) (stat_sim > *stat);
			/* resampling with individual f11,f21; not used */
			/*
			twosample_incidence_lwy_process(event, group, n, y1, y2, f11, f12, f21, f22,
				g, test_process_sim);
			stat_sim = ks_stat_cum(test_process_sim, n, i1, i2);
			*pval_sim_indiv += (double) (stat_sim > *stat);
			*/
		}
		*pval_sim /= *nsim;
		/* *pval_sim_indiv /= *nsim; */
	}
	
	PutRNGstate();
}

void twosample_incidence_ks_process(int *event, int *group, int *n, int *y1, int *y2,
	double *f11, double *f12, double *f21, double *f22, double *s1, double *s2,
	double *test_process)
{
	int i;
	
	/* compute number at risk (y1,y2) and increments of nelson--aalen (s1,s2) */
	if (group[*n-1]==1) {
		y1[*n-1] = 1;
		y2[*n-1] = 0;
		s1[*n-1] = (double) (event[*n-1]>0);
		s2[*n-1] = 0.;
	} else {
		y1[*n-1] = 0;
		y2[*n-1] = 1;
		s1[*n-1] = 0.;
		s2[*n-1] = (double) (event[*n-1]>0);
	}
	for (i=*n-2; i>=0; i--) {
		if (group[i]==1) {
			y1[i] = y1[i+1] + 1;
			y2[i] = y2[i+1];
			s1[i] = (double) (event[i]>0)/y1[i];
			s2[i] = 0.;
		} else {
			y1[i] = y1[i+1];
			y2[i] = y2[i+1] + 1;
			s1[i] = 0.;
			s2[i] = (double) (event[i]>0)/y2[i];
		}
	}
	
	/* compute exp of minus nelson--aalen (s1,s2), cif's for event type 1 (f11,f21) */
	/* (for i=0 and i=*n-1 computed separately) */

	f11[0] = 1. * (double) ((group[0]==1)*(event[0]==1)) * (y1[0]>0 ? 1./y1[0] : 0.);
	f12[0] = 1. * (double) ((group[0]==1)*(event[0]==2)) * (y1[0]>0 ? 1./y1[0] : 0.);
	f21[0] = 1. * (double) ((group[0]==2)*(event[0]==1)) * (y2[0]>0 ? 1./y2[0] : 0.);
	f22[0] = 1. * (double) ((group[0]==2)*(event[0]==2)) * (y2[0]>0 ? 1./y2[0] : 0.);
	test_process[0] = f21[0] - f11[0];
	for (i=1; i<*n; i++) {
		s1[i] += s1[i-1];
		s1[i-1] = exp(-s1[i-1]);
		f11[i] = f11[i-1] + s1[i-1] * (double) ((group[i]==1)*(event[i]==1)) * (y1[i]>0 ? 1./y1[i] : 0.);
		f12[i] = f12[i-1] + s1[i-1] * (double) ((group[i]==1)*(event[i]==2)) * (y1[i]>0 ? 1./y1[i] : 0.);
		s2[i] += s2[i-1];
		s2[i-1] = exp(-s2[i-1]);
		f21[i] = f21[i-1] + s2[i-1] * (double) ((group[i]==2)*(event[i]==1)) * (y2[i]>0 ? 1./y2[i] : 0.);
		f22[i] = f22[i-1] + s2[i-1] * (double) ((group[i]==2)*(event[i]==2)) * (y2[i]>0 ? 1./y2[i] : 0.);
 		test_process[i] = f21[i] - f11[i];
	}
	s1[*n-1] = exp(-s1[*n-1]);
	s2[*n-1] = exp(-s2[*n-1]);
}

void twosample_incidence_lwy_process(int *event, int *group, int *n, int *y1, int *y2,
	double *f11, double *f12, double *f21, double *f22, double *g, double *test_process_sim)
{
	/* computes the Lin--Wei--Ying simulated path of the test process under H0 */
	/* see the martingale representation in Lin (1997), beware of misprints there */
	/* n independent N(0,1) variables are on input in g */
	int i;
	double temp1, temp2;
	double sum1a = 0., sum1b = 0., sum2a = 0., sum2b = 0.;
	
	for (i=0; i<*n; i++) {
		temp1 = (y1[i]>0 ? g[i]/y1[i]*(group[i]==1) : 0.);
		sum1a += (1.-f12[i])*temp1*(event[i]==1) + f11[i]*temp1*(event[i]==2);
		sum1b += temp1*(event[i]>0);
		temp2 = (y2[i]>0 ? g[i]/y2[i]*(group[i]==2) : 0.);
		sum2a += (1.-f22[i])*temp2*(event[i]==1) + f21[i]*temp2*(event[i]==2);
		sum2b += temp2*(event[i]>0);
		test_process_sim[i] = (sum2a - f21[i]*sum2b)
			- (sum1a - f11[i]*sum1b);
	}
}

