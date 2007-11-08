#include <R.h>
#include <Rmath.h>
#include "util_ks.h"
#include "surv_kscmad.h"


/* two-sample weighted score (logrank) process tests (Kolmogorov--Smirnov, Cramer--von Mises,
   Anderson--Darling) based on both the original weighted logrank process (asymptotically
   Brownian motion) and on a Hall--Wellner transformed transformed process (asympt. Brownian
   bridge)
*/
/* p-values can be obtained from asymptotics (for KS only), LWY simulations, permutations,
   bootstrap
*/


void twosample_ks_cm_ad(double *time, int *event, int *group, int *n, double *rho, double *gamma,
	int *nsim, int *nperm, int *nboot, double *du, int *nsim_plot, int *nperm_plot, int *nboot_plot,
	double *du_sim_plot, double *du_perm_plot, double *du_boot_plot,
	double *stat_ks_w, double *pval_ks_w_sim, double *pval_ks_w_perm, double *pval_ks_w_boot, double *pval_ks_w_asympt,
	double *stat_ks_b, double *pval_ks_b_sim, double *pval_ks_b_perm, double *pval_ks_b_boot, double *pval_ks_b_asympt,
	double *stat_cm_w, double *pval_cm_w_sim, double *pval_cm_w_perm, double *pval_cm_w_boot,
	double *stat_cm_b, double *pval_cm_b_sim, double *pval_cm_b_perm, double *pval_cm_b_boot,
	double *stat_ad_w, double *pval_ad_w_sim, double *pval_ad_w_perm, double *pval_ad_w_boot,
	double *stat_ad_b, double *pval_ad_b_sim, double *pval_ad_b_perm, double *pval_ad_b_boot)
{
	
	int i,j,k,r;
	int y1,y2,y1bak,y2bak,first_event,first_event_g;
	double temp1,temp2,temp3,temp4,temp5,temp6;
	
	/* various work arrays of length n */
	double *dsigma; /* increments of variance */
	double *ks_b; /* weight for KS for Brownian bridge */
	double *cm_w,*cm_b,*ad_w,*ad_b; /* integrators for CM and AD for Brownian motion (.w) and bridge (.b) */
	double *g; /* normal random numbers for LWY simulation */
	double *du_temp; /* increments of simulated, permutation or bootstrap test process */
	int *eventg; /* bootstrap sample */
	int *freqg; /* bootstrap sample frequencies */
	
	GetRNGstate();
	
	/* allocate everything at once */
	dsigma = (double *) R_alloc(7*(*n),sizeof(double)); /* 6*n, not 7*n, because cm.w is dsigma */
	ks_b = dsigma + *n;
	cm_w = dsigma;
	cm_b = dsigma + *n*2;
	ad_w = dsigma + *n*3;
	ad_b = dsigma + *n*4;
	g = dsigma + *n*5;
	du_temp = dsigma + *n*6;
	
	/* initialise at-risk values */
	y1 = 0;
	y2 = 0;
	for (i=0; i<*n; i++) {
		if (group[i] == 1) {
			++y1;
		} else {
			++y2;
		}
	}
	y1bak = y1;
	y2bak = y2;
	/* find the first event (needed in AD test to avoid division by zero) */
	for (first_event=0; (event[first_event]==0)&&(first_event<*n) ; first_event++)
		;
	
	/* compute the OBSERVED process and statistics (g is 1) */
	ks_cm_ad_process(event, group, n, &y1, &y2, rho, gamma, &first_event, du, dsigma, ks_b, cm_w, cm_b, ad_w, ad_b);
	y1 = y1bak;
	y2 = y2bak;
	for (i=0; i<*n; i++)
		g[i] = 1.;
	ks_cm_ad_stat(n,du,g,ks_b,cm_w,cm_b,ad_w,ad_b,stat_ks_w,stat_ks_b,stat_cm_w,stat_cm_b,stat_ad_w,stat_ad_b);
	
	/* now do the true LWY SIMULATIONS (g is random) */
	/* LWY must be done before permutations and bootstrap (permutations and bootstrap destroy group) */
	*pval_ks_w_sim = *pval_ks_b_sim = *pval_cm_w_sim = *pval_cm_b_sim = *pval_ad_w_sim = *pval_ad_b_sim = 0.;
	if (*nsim>0) {
		for (j=0; j<*nsim; j++) {
			for (i=0; i<*n; i++) {
				g[i] = norm_rand();
			}
			/* store this LWY-simulated process for plotting */
			if (j<*nsim_plot) {
				for (i=0; i<*n; i++)
					du_sim_plot[i+j**n] = g[i]*du[i];
			}
			ks_cm_ad_stat(n,du,g,ks_b,cm_w,cm_b,ad_w,ad_b,&temp1,&temp2,&temp3,&temp4,&temp5,&temp6);
			*pval_ks_w_sim += (temp1>*stat_ks_w);
			*pval_ks_b_sim += (temp2>*stat_ks_b);
			*pval_cm_w_sim += (temp3>*stat_cm_w);
			*pval_cm_b_sim += (temp4>*stat_cm_b);
			*pval_ad_w_sim += (temp5>*stat_ad_w);
			*pval_ad_b_sim += (temp6>*stat_ad_b);
		}
		*pval_ks_w_sim /= (double) *nsim;
		*pval_ks_b_sim /= (double) *nsim;
		*pval_cm_w_sim /= (double) *nsim;
		*pval_cm_b_sim /= (double) *nsim;
		*pval_ad_w_sim /= (double) *nsim;
		*pval_ad_b_sim /= (double) *nsim;
	}
	
	/* BOOTSTRAP */
	*pval_ks_w_boot = *pval_ks_b_boot = *pval_cm_w_boot = *pval_cm_b_boot = *pval_ad_w_boot = *pval_ad_b_boot = 0.;
	for (i=0; i<*n; i++)
		g[i] = 1.;
	if (*nboot>0) {
		/* allocate arrays for the bootstrap sample */
		eventg = (int *) R_alloc(2**n,sizeof(int));
		freqg = eventg + *n;
		for (i=0; i<*nboot; i++) {
			/* compute frequencies of the original data in the bootstrap sample (we want to avoid sorting the bootstrap sample) */
			for (j=0; j<*n; j++)
				freqg[j] = 0;
			for (j=0; j<*n; j++) { /* n times randomly select an index and increase its frequency */
				r = (unif_rand()**n);
				freqg[r] += 1;
			}
			/* now take corresponding number of observations from the data */
			r = 0;
			for (j=0; j<*n; j++) {
				for (k=0; k<freqg[j]; k++) { /* put the j-th observation freqg[j]-times to the bootstrap sample */
					eventg[r] = event[j];
					++r;
				}
			}
			for (first_event_g=0; (event[first_event_g]==0)&&(first_event_g<*n) ; first_event_g++)
				;
			for (j=0; j<*n; j++) { /* each observation is randomly assigned to one of the groups; group from input is reused */
				if (unif_rand() < ((double) y1/(y1+y2))) {
					group[j] = 1;
					--y1;
				} else {
					group[j] = 2;
					--y2;
				}
			}
			y1 = y1bak;
			y2 = y2bak;
			ks_cm_ad_process(eventg, group, n, &y1, &y2, rho, gamma, &first_event_g, du_temp, dsigma, ks_b, cm_w, cm_b, ad_w, ad_b);
			/* store this bootstrapped process for plotting */
			if (i<*nboot_plot) {
				for (j=0; j<*n; j++)
					du_boot_plot[j+i**n] = du_temp[j];
			}
			y1 = y1bak;
			y2 = y2bak;
			ks_cm_ad_stat(n,du_temp,g,ks_b,cm_w,cm_b,ad_w,ad_b,&temp1,&temp2,&temp3,&temp4,&temp5,&temp6);
			*pval_ks_w_boot += (temp1>*stat_ks_w);
			*pval_ks_b_boot += (temp2>*stat_ks_b);
			*pval_cm_w_boot += (temp3>*stat_cm_w);
			*pval_cm_b_boot += (temp4>*stat_cm_b);
			*pval_ad_w_boot += (temp5>*stat_ad_w);
			*pval_ad_b_boot += (temp6>*stat_ad_b);
		
		}
		*pval_ks_w_boot /= (double) *nboot;
		*pval_ks_b_boot /= (double) *nboot;
		*pval_cm_w_boot /= (double) *nboot;
		*pval_cm_b_boot /= (double) *nboot;
		*pval_ad_w_boot /= (double) *nboot;
		*pval_ad_b_boot /= (double) *nboot;
	}
	
	/* PERMUTATION tests (g is 1) */
	*pval_ks_w_perm = *pval_ks_b_perm = *pval_cm_w_perm = *pval_cm_b_perm = *pval_ad_w_perm = *pval_ad_b_perm = 0.;
	for (i=0; i<*n; i++)
		g[i] = 1.;
	if (*nperm>0) {
		for (i=0; i<*nperm; i++) {
			for (j=0; j<*n; j++) { /* each observation is randomly assigned to one of the groups; group from input is reused (destroyed) */
				if (unif_rand() < ((double) y1/(y1+y2))) {
					group[j] = 1;
					--y1;
				} else {
					group[j] = 2;
					--y2;
				}
			}
			y1 = y1bak;
			y2 = y2bak;
			ks_cm_ad_process(event, group, n, &y1, &y2, rho, gamma, &first_event, du_temp, dsigma, ks_b, cm_w, cm_b, ad_w, ad_b);
			/* store this permutation process for plotting */
			if (i<*nperm_plot) {
				for (j=0; j<*n; j++)
					du_perm_plot[j+i**n] = du_temp[j];
			}
			y1 = y1bak;
			y2 = y2bak;
			ks_cm_ad_stat(n,du_temp,g,ks_b,cm_w,cm_b,ad_w,ad_b,&temp1,&temp2,&temp3,&temp4,&temp5,&temp6);
			*pval_ks_w_perm += (temp1>*stat_ks_w);
			*pval_ks_b_perm += (temp2>*stat_ks_b);
			*pval_cm_w_perm += (temp3>*stat_cm_w);
			*pval_cm_b_perm += (temp4>*stat_cm_b);
			*pval_ad_w_perm += (temp5>*stat_ad_w);
			*pval_ad_b_perm += (temp6>*stat_ad_b);
		
		}
		*pval_ks_w_perm /= (double) *nperm;
		*pval_ks_b_perm /= (double) *nperm;
		*pval_cm_w_perm /= (double) *nperm;
		*pval_cm_b_perm /= (double) *nperm;
		*pval_ad_w_perm /= (double) *nperm;
		*pval_ad_b_perm /= (double) *nperm;
	}
	
	/* ASYMPTOTIC p-values of KS */
	*pval_ks_w_asympt = 1. - psupw(*stat_ks_w,50);
	*pval_ks_b_asympt = 1. - psupb(*stat_ks_b,.5,50);
	
	PutRNGstate();
	
}

void ks_cm_ad_stat(int *n, double *du, double *g, double *ks_b, double *cm_w, double *cm_b,
	double *ad_w, double *ad_b, double *stat_ks_w, double *stat_ks_b,
	double *stat_cm_w, double *stat_cm_b, double *stat_ad_w, double *stat_ad_b)
{
	double u = 0., u2;
	int i;
	
	*stat_ks_w = *stat_ks_b = *stat_cm_w = *stat_cm_b = *stat_ad_w = *stat_ad_b = 0.;
	for (i=0; i<*n; i++) {
		u += g[i]*du[i];
		if (fabs(u)>*stat_ks_w)
			*stat_ks_w = fabs(u);
		if (fabs(u*ks_b[i])>*stat_ks_b)
			*stat_ks_b = fabs(u*ks_b[i]);
		u2 = u*u;
		*stat_cm_w += u2*cm_w[i];
		*stat_cm_b += u2*cm_b[i];
		*stat_ad_w += u2*ad_w[i];
		*stat_ad_b += u2*ad_b[i];
	}
}

void ks_cm_ad_process(int *event, int *group, int *n, int *y1, int *y2, double *rho, double *gamma,
	int *first_event, double *du, double *dsigma, double *ks_b, double *cm_w, double *cm_b,
	double *ad_w, double *ad_b)
{
	double temp1,temp2,temp3;
	int i;
	
	/* compute increments du of the logrank process and dsigma of the variance function */
	temp1 = 0.; /* this will be current value of left-continuous nelson--aalen (computed on the fly, not in advance) */
	temp2 = 0.; /* this will contain sigma(tau) */
	/* instead of left-continuous kaplan--meier we use exp of minus left-cont. nelson--aalen */
	for (i=0; i<*n; i++) {
		if (event[i] == 1) { /* event, hence u and sigma have increments */
			if (group[i] == 1) {
				du[i] = - exp(-*rho*temp1)*pow(1.-exp(-temp1),*gamma) / (*y1+*y2);
				dsigma[i] = du[i]*du[i]**y1**y2;
				du[i] *= *y2;
			} else {
				du[i] = exp(-*rho*temp1)*pow(1.-exp(-temp1),*gamma) / (*y1+*y2);
				dsigma[i] = du[i]*du[i]**y1**y2;
				du[i] *= *y1;
			}
			temp2 += dsigma[i];
			temp1 += 1./(*n-i); /* this is used in the next step, thus left continuous */
		} else { /* no event, no increment */
			du[i] = 0.;
			dsigma[i] = 0.;
		}
		if (group[i] == 1) {
			--*y1;
		} else {
			--*y2;
		}
	}
	
	/* standardise by variance at the end (stored in temp2) */
	/* also prepare various processes for computation of ks, cm, ad from brownian motion and bridge (HW) */
	temp1 = sqrt(temp2); /* sqrt of sigma(tau) */
	temp3 = 0.; /* temp3 will be the current value of cumulative sigma */
	for (i=0; i<*n; i++) {
		du[i] /= temp1;
		dsigma[i] /= temp2;
		temp3 += dsigma[i];
		ks_b[i] = 1./(1.+temp3);
		cm_b[i] = dsigma[i]*pow(ks_b[i],4);
		if (i>=*first_event) {
			ad_w[i] = dsigma[i]/temp3;
			ad_b[i] = ad_w[i]*ks_b[i]*ks_b[i];
		} else {
			ad_w[i] = 0.;
			ad_b[i] = 0.;
		}
	}
	
}

