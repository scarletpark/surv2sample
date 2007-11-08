#include <R.h>
#include <Rmath.h>
#include "proprate_estim.h"
#include "util_neyman.h"
#include "util_datadriven.h"
#include "proprate_neyman.h"

/* NEYMAN'S TEST of fit of the two-sample proportional rate transformation model */

void twosample_proprate_neyman(double *time, int *event, int *group, int *n, int *model,
	double *beta, int *maxiter, double *eps, int *flag, int *d, int *d0, int *timetransf,
	int *basis, double *penalty, double *choltol, double *score, double *sigma,
	double *stat_d, double *pval_d, double *stat_d0, double *pval_d0,
	double *stat_nested, double *stats_nested, double *stats_nested_penal, int *s_nested,
	double *pval_nested_asympt, double *pval_nested_h)
{
	int i,j,k;
	int y1,y2,y1work,y2work;
	double temp1;
	double theta;
	void (*p_basis)();
	
	double (*p_q)(),(*p_qdot)();
	if (*model==0) { /* Cox */
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
	
	double *u1_process,*u1,*d11_process,*d11,*sigma11,*u1_process_sim,*d21,*sigma21;
	u1_process = (double *) R_alloc(*n+1+*n+1+1+*n+*d+*d,sizeof(double));
	u1 = u1_process+*n;
	d11_process = u1+1;
	d11 = d11_process+*n;
	sigma11 = d11+1;
	u1_process_sim = sigma11+1;
	d21 = u1_process_sim+*n;
	sigma21 = d21+*d;
	
	double *logliks = (double *) R_alloc(2,sizeof(double));
	
	double *sigmawork;
	sigmawork = (double *) R_alloc(*d**d,sizeof(double));
	double *workd2;
	workd2 = (double *) R_alloc(*d,sizeof(double));
	
	y1 = 0;
	y2 = 0;
	for (i=0; i<*n; i++) {
		if (group[i] == 1) {
			++y1;
		} else {
			++y2;
		}
	}
	
	double *a, *psi;
	a = (double *) R_alloc(*n+*n**d,sizeof(double));
	psi = a + *n;
	
	/* compute time-transformation */
	if (*timetransf == 2) { /* nelson--aalen (left-continuous) */
		a[0] = 0.;
		for (i=1; i<*n; i++)
			a[i] = a[i-1] + (double) event[i-1]/(*n-i+1.); /* Efron's handling of ties */
		for (i=0; i<*n; i++)
			a[i] /= a[*n-1];
	} else {
		if (*timetransf == 1) { /* distribution function (left-continuous), that is one minus exp of minus left-cont. nelson--aalen */
			a[0] = 0.;
			for (i=1; i<*n; i++) {
				a[i] = a[i-1] + (double) event[i-1]/(*n-i+1.);
				a[i-1] = 1. - exp(-a[i-1]);
			}
			a[*n-1] = 1. - exp(-a[*n-1]);
			for (i=0; i<*n; i++)
				a[i] /= a[*n-1];
		} else {
			if (*timetransf == 3) { /* identity */
				for (i=0; i<*n; i++) a[i] = time[i]/time[*n-1];
			} else {
				if (*timetransf == 4) { /* sigma (left-cont.); this is impossible for add. haz.; sigma21 != d21 (LS estimation is not ML) (unlike for Cox) */
				}
			}
		}
	}
	
	if (*basis == 1) {
		p_basis = legendre_basis;
	} else {
		if (*basis == 2) {
			p_basis = cosine_basis;
		} else {
			if (*basis == 3) {
				p_basis = cosine1_basis;
			} else {
				if (*basis == 4) {
/* 					p_basis = indicator_basis; */
				}
			}
		}
	}
	
	/* compute the basis at all times */
	(*p_basis)(a, *n, psi, *d, 0, 1, 0);
	
	/* estimate the model */
	twosample_proprate_estim(event, group, n, beta, maxiter, eps, flag, p_q, p_qdot, q1, q2, q1dot, q2dot,
		u1_process, u1, d11_process, d11, sigma11, logliks);
	
	/* compute the score vector and its variance */
	twosample_proprate_neyman_score(score, sigma, event, group, n, beta, q1, q2, q1dot, q2dot,
		psi, d, u1_process, d11_process, d11, sigma11, d21, sigma21);
	
	/* store diag of sigma */
	/* (not necessary, datadriven_stats_nested doesn't destroy sigma) */
	for (j=0; j<*d; j++)
		a[j] = sigma[j+j**d];

	/* fixed and data-driven test */
	datadriven_stats_nested(score, sigma, d, d0, stats_nested, stats_nested_penal,
		s_nested, stat_d, stat_d0, stat_nested, sigmawork, workd2, penalty, choltol);
	*pval_d = pchisq(*stat_d, (double) *d, 0, 0);
	*pval_d0 = pchisq(*stat_d0, (double) *d0, 0, 0);
	*pval_nested_asympt = pchisq(*stat_nested, (double) *d0, 0, 0);
	if (*d0==1)
		h_approx_nested(stat_nested, n, pval_nested_h); /* two-term approx. only for d0=1 */
	
	for (j=0; j<*d; j++) {
		sigma[j+j**d] = a[j]; /* reconstruct diag of sigma */
		for (k=0; k<j; k++) {
			sigma[j+k**d] = sigma[k+j**d]; /* copy upper triangle to lower */
		}
	}
}

void twosample_proprate_neyman_score(double *score, double *sigma, int *event, int *group, int *n,
	double *beta, double *q1, double *q2, double *q1dot, double *q2dot, double *psi, int *d,
	double *u1_process, double *d11_process, double *d11, double *sigma11, double *d21,
	double *sigma21)
{
	/* computes increments of the score process (u1_process) and the score vector (u1) at beta (U1 in the paper) */
	/* the same for its derivative (D11 in the paper) */
	/* loglik is the log of the modified partial likelihood */
	/* sigma11 is the cumulative sigma11 at the end (time tau) */
	/* q1, q2 are q(u)=lambda0(Lambda0^{-1}(u)) evaluated at Nelson--Aalen estimators for samples 1, 2;
	/* q(u)=1 for prop. haz., q(u)=exp(-u) for prop. odds */
	int i,j,k;
	int y1work = 0, y2work = 0;
	double theta = exp(*beta), temp, h11 = 0., h12 = 0., r11 = 0., r12 = 0., *r21, *r22, *h21, *h22;
	
	r21 = (double *) R_alloc(4**d,sizeof(double));
	r22 = r21 + *d;
	h21 = r22 + *d;
	h22 = h21 + *d;
	
	for (j=0; j<*d; j++) {
		score[j] = d21[j] = sigma21[j] = r21[j] = r22[j] = h21[j] = h22[j] = 0.;
		for (k=0; k<*d; k++) {
			sigma[j+k**d] = 0.;
		}
	}
	y1work = y2work = 0;
	*sigma11=0.;
	for (i=*n-1; i>=0; i--) {
		/* loop over times; from last to first because r1,r2 in h1,h2 are
		   cumulative from t to tau; notation different from the paper */
		if(group[i]==1) {
			++y1work;
		} else {
			++y2work;
		}
		if (event[i]==1) {
			temp = 1./(q1[i]*y1work+theta*q2[i]*y2work);
			r11 += y1work*q1dot[i]*theta*y2work*q2[i]*temp*temp;
			r12 += y1work*q1[i]*theta*y2work*q2dot[i]*temp*temp;
			h11 = theta*q2[i]*y2work*temp - (y1work>0 ? 1./y1work : 0.) * r11;
			h12 = q1[i]*y1work*temp - (y2work>0 ? 1./y2work : 0.) * r12;
			*sigma11 += (h11*h11*y1work*q1[i]+h12*h12*theta*y2work*q2[i])*temp;
			for (j=0; j<*d; j++) {
				score[j] += psi[i+j**n]*u1_process[i];
				d21[j] += psi[i+j**n]*d11_process[i];
				r21[j] += psi[i+j**n]*y1work*q1dot[i]*theta*y2work*q2[i]*temp*temp;
				r22[j] += psi[i+j**n]*y1work*q1[i]*theta*y2work*q2dot[i]*temp*temp;
				h21[j] = psi[i+j**n]*theta*q2[i]*y2work*temp - (y1work>0 ? 1./y1work : 0.) * r21[j];
				h22[j] = psi[i+j**n]*q1[i]*y1work*temp - (y2work>0 ? 1./y2work : 0.) * r22[j];
				sigma21[j] += (h21[j]*h11*y1work*q1[i]+h22[j]*h12*theta*y2work*q2[i])*temp;
				for (k=0; k<=j; k++) { /* sigma22 */
					sigma[k+j**d] += (h21[k]*h21[j]*y1work*q1[i]+h22[k]*h22[j]*theta*y2work*q2[i])*temp;
				}
			}
		}
	}

	for (j=0; j<*d; j++) {
		for (k=j; k<*d; k++) {
			sigma[j+k**d] += - d21[j]/(*d11)*sigma21[k] - sigma21[j]/(*d11)*d21[k] + d21[j]/(*d11)*(*sigma11)/(*d11)*d21[k];
		}
	}
}

