#include <R.h>
#include <Rmath.h>
#include "proprate_estim.h"
#include "util_ks.h"
#include "proprate_ks.h"

/* KOLMOGOROV--SMIRNOV TEST of fit of the two-sample proportional rate model */

void twosample_proprate_ks(int *event, int *group, int *n, int *model, double *beta,
	int *maxiter, double *eps, int *flag, int *nsim, double *u1_process,
	double *stat_ks, double *pval_ks_sim, double *test_process_plot_sim, int *nsim_plot)
{
	/* u1_process (observed path) and test_process_plot_sim (simulated) are increments */
	int i,j;
	int y1,y2,y1work,y2work;
	double stat_ks_sim,stat_ks_sim_2;
	double theta,temp;
	
	double (*p_q)(),(*p_qdot)();
	if (*model==0) { /* prop. hazards */
		p_q = q_prophaz;
		p_qdot = qdot_prophaz;
	} else { /* prop. odds */
		p_q = q_propodds;
		p_qdot = qdot_propodds;
	}

	double *q1,*q2,*q1dot,*q2dot;
	q1 = (double *) R_alloc(4**n,sizeof(double));
	q2 = q1+*n;
	q1dot = q2+*n;
	q2dot = q1dot+*n;
	
	double *u1,*d11_process,*d11,*sigma11,*u1_process_sim;
	u1 = (double *) R_alloc(1+*n+1+1+*n,sizeof(double));
	d11_process = u1+1;
	d11 = d11_process+*n;
	sigma11 = d11+1;
	u1_process_sim = sigma11+1;
	
	double *logliks = (double *) R_alloc(2,sizeof(double));
	
	y1 = 0;
	y2 = 0;
	for (i=0; i<*n; i++) {
		if (group[i] == 1) {
			++y1;
		} else {
			++y2;
		}
	}
	
	twosample_proprate_estim(event, group, n, beta, maxiter, eps, flag, p_q, p_qdot, q1, q2, q1dot, q2dot,
		u1_process, u1, d11_process, d11, sigma11, logliks);
	*stat_ks = ks_stat_incr(u1_process, n); /* u1_process contains increments */
	
	/* LWY SIMULATION */
	/* must always be *nsim>=*nsim_plot */
	if (*nsim>0) {
		GetRNGstate();
		*pval_ks_sim = 0.;
		theta = exp(*beta);
		y1work = y1;
		y2work = y2;
		/* run simulations */
		for (j=0; j<*nsim_plot; j++) { /* first *nsim_plot simulated paths stored for plotting */
			y1work = y1;
			y2work = y2;
			twosample_proprate_ks_sim(event, group, n, &y1work, &y2work, &theta, q1, q2,
					q1dot, q2dot, d11_process, d11, u1_process_sim);
			stat_ks_sim = ks_stat_incr(u1_process_sim, n);
			*pval_ks_sim += (stat_ks_sim > *stat_ks);
			y1work = y1;
			y2work = y2;
			/* save the simulated process for plotting */
			for (i=0; i<*n; i++)
				test_process_plot_sim[i+j**n] = u1_process_sim[i];
		}
		for (j=*nsim_plot; j<*nsim; j++) {
			y1work = y1;
			y2work = y2;
			twosample_proprate_ks_sim(event, group, n, &y1work, &y2work, &theta, q1, q2,
					q1dot, q2dot, d11_process, d11, u1_process_sim);
			stat_ks_sim = ks_stat_incr(u1_process_sim, n);
			*pval_ks_sim += (stat_ks_sim > *stat_ks);
			y1work = y1;
			y2work = y2;
		}
		*pval_ks_sim /= *nsim;
		PutRNGstate();
	}
}

void twosample_proprate_ks_sim(int *event, int *group, int *n, int *y1, int *y2, double *theta,
	double *q1, double *q2, double *q1dot, double *q2dot, double *d11_process, double *d11,
	double *u1_process_sim)
{
	/* Lin--Wei--Ying style simulation of the test process */
	/* (martingale increments replaced by N(0,1) variables) */
	int i;
	double g, u1_sim = 0., na1 = 0., na2 = 0., temp;
	
	u1_sim = 0.;
	for (i=0; i<*n; i++) {
		if (event[i] == 1) {
			g = norm_rand();
			temp = 1./(q1[i]**y1+*theta*q2[i]**y2);
			if (group[i]==1) {
				na1 += g/(*y1); /* nelson--aalen */
				u1_process_sim[i] = - *theta*q2[i]**y2*temp*g;
			} else {
				na2 += g/(*y2); /* nelson--aalen */
				u1_process_sim[i] = q1[i]**y1*temp*g;
			}
			u1_process_sim[i] += (na1*q1dot[i]*q2[i] - na2*q1[i]*q2dot[i]) * *theta**y1**y2*temp*temp;
		} else {
			u1_process_sim[i] = 0.;
		}
		u1_sim += u1_process_sim[i];
		if (group[i]==1) {
			--*y1;
		} else {
			--*y2;
		}
	}
	
	for (i=0; i<*n; i++) {
		u1_process_sim[i] -= d11_process[i]/(*d11)*u1_sim;
	}
}

