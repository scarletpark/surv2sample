#include <R.h>
#include <Rmath.h>
#include "util_neyman.h"
#include "util_datadriven.h"
#include "cif_misc.h"
#include "cif_neyman.h"

/* NEYMAN'S SMOOTH TEST OF EQUALITY OF CUMULATIVE INCIDENCE FUNCTIONS IN TWO SAMPLES */
/* (comparing CIF's for failure cause 1) */

/* also used for logrank's test: if logrank==1, d=1, G^r weighted subdistribution logrank-type test */

void twosample_incidence_neyman(double *time, int *event, int *group, int *n,
	int *d, int *d0, int *basis, double *penalty, int *logrank, double *logrank_r,
	double *choltol, double *score,	double *sigma, double *f11, double *f12, double *f21,
	double *f22, double *stat_d, double *pval_d, double *stat_d0, double *pval_d0,
	double *stat_nested, double *stats_nested, double *stats_nested_penal,
	int *s_nested, double *pval_nested_asympt, double *pval_nested_h)
{
	/* penalty = log(n) for bic */
	
	int i,j,k,n1,n2;
	double temp,temp1,temp2;
	void (*p_basis)(); /* pointer to the function evaluating the basis (legendre, cos,...) */

	int *y1,*y2; /* vectors of length n */
	double *s1,*s2; /* vectors of length n */
	double *h1,*h2; /* vectors of length d */
	double *q1,*q2,*psi; /* n by d matrices */
	double *rho; /* n by n matrix */
	double *workd1,*workd2; /* length d */
	double *f01; /* length n */
	double *sigmawork; /* d by d */
	/* allocate these arrays */
	y1 = (int *) R_alloc(2**n,sizeof(int));
	y2 = y1 + *n;
	s1 = (double *) R_alloc(2**n+2**d+3**n**d+*n**n+2**d+*n+*d**d,sizeof(double));
	s2 = s1 + *n;
	h1 = s2 + *n;
	h2 = h1 + *d;
	q1 = h2 + *d;
	q2 = q1 + *n**d;
	psi = q2 + *n**d;
	rho = psi + *n**d;
	workd1 = rho + *n**n;
	workd2 = workd1 + *d;
	f01 = workd2 + *d;
	sigmawork = f01 + *n;
	
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
	
	/* computing the test statistic */
	/* compute the basis psi */
	twosample_incidence_psi(time, event, group, n, y1, y2, s1, s2, f11, psi, d,
		p_basis, logrank, logrank_r);
	
	/* individual samples estimator of sigma */
	/* compute the score vector */
	twosample_incidence_score(score, sigma, d, /* pooled = */ 0, time, event, group,
		n, y1, y2, s1, s2, f11, f12, f21, f22, f01, psi, h1, h2, q1, q2, rho, workd1);
	
	/* fixed and data-driven test */
	datadriven_stats_nested(score, sigma, d, d0, stats_nested, stats_nested_penal,
		s_nested, stat_d, stat_d0, stat_nested, sigmawork, workd2, penalty, choltol);
	*pval_d = pchisq(*stat_d, (double) *d, 0, 0);
	*pval_d0 = pchisq(*stat_d0, (double) *d0, 0, 0);
	*pval_nested_asympt = pchisq(*stat_nested, (double) *d0, 0, 0);
	if (*d0==1)
		h_approx_nested(stat_nested, n, pval_nested_h); /* two-term approx. only for d0=1 */
}

void twosample_incidence_psi(double *time, int *event, int *group, int *n,
	int *y1, int *y2, double *a1, double *a2, double *condf01,
	double *psi, int *d, void (*p_basis)(), int *logrank, double *logrank_r)
{
	/* computes the psi basis for Neyman's smooth test */
	/* the time-transformation is based on the null estimate of f01 */
	/* if logrank==1, logrank's (1988) G^r logrank-type test */
	
	int i;
	
	/* compute y1,y2 (# at risk), and increments of a1,a2 (nelson--aalen) */
	if (group[*n-1]==1) {
		y1[*n-1] = 1;
		y2[*n-1] = 0;
		a1[*n-1] = (double) (event[*n-1]>0);
		a2[*n-1] = 0.;
	} else {
		y1[*n-1] = 0;
		y2[*n-1] = 1;
		a1[*n-1] = 0.;
		a2[*n-1] = (double) (event[*n-1]>0);
	}
	for (i=*n-2; i>=0; i--) {
		if (group[i]==1) {
			y1[i] = y1[i+1] + 1;
			y2[i] = y2[i+1];
			a1[i] = (double) (event[i]>0)/y1[i];
			a2[i] = 0.;
		} else {
			y1[i] = y1[i+1];
			y2[i] = y2[i+1] + 1;
			a1[i] = 0.;
			a2[i] = (double) (event[i]>0)/y2[i];
		}
	}
	
	/* time-transformation = condf01 = standardised null estimate of cum. inc. fnc. for type 1 */
	condf01[0] = ((double) (event[0]==1)) / (y1[0]+y2[0]);
	for (i=1; i<*n; i++) {
		condf01[i] = condf01[i-1] + ((double) (event[i]==1)) / ((double) y1[i]*exp(a1[i-1]) + (double) y2[i]*exp(a2[i-1]));
		a1[i] += a1[i-1];
		a2[i] += a2[i-1];
	}
	
	if (*logrank==1) { /* logrank's test */
		for (i=0; i<*n; i++)
			psi[i] = pow(1.-condf01[i],*logrank_r);
	} else { /* Neyman's test */
		for (i=0; i<*n; i++)
			condf01[i] /= condf01[*n-1];
	 	(*p_basis)(condf01, *n, psi, *d-1, 1, 1, 0);
	 }
}

void twosample_incidence_score(double *score, double *sigma, int *d, int pooled,
	double *time, int *event, int *group, int *n, int *y1, int *y2,
	double *s1, double *s2, double *f11, double *f12, double *f21, double *f22, double *f01,
	double *psi, double *h1, double *h2, double *q1, double *q2, double *rho, double *work)
{
	/* computes the score vector for the Neyman smooth test */
	/* fij = cum. inc. fnc. for group i, failure type j */
	/* work of size d is used in the variance */
	
	int i,j;
	double temp,temp1,temp2, temp0;
	double r1,r2, r0;
	
	for (i=0; i<*d; i++) {
		h1[i] = h2[i] = score[i] = 0.;
		for (j=0; j<*d; j++)
			sigma[i+j**d] = 0.;
	}
	
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
	
	/* compute exp of minus nelson--aalen (s1,s2), cum. incidence for event type 1 (f11,f21) */
	/* and score, h1, h2, q1, q2 */
	/* (for i=0 and i=*n-1 computed separately) */

	s1[1] += s1[0];
	s1[0] = exp(-s1[0]);
	f11[0] = 1. * (double) ((group[0]==1)*(event[0]==1)) * (y1[0]>0 ? 1./y1[0] : 0.);
	f12[0] = 1. * (double) ((group[0]==1)*(event[0]==2)) * (y1[0]>0 ? 1./y1[0] : 0.);
	r1 = (1.-0.)/1.*y1[0];
	s2[1] += s2[0];
	s2[0] = exp(-s2[0]);
	f21[0] = 1. * (double) ((group[0]==2)*(event[0]==1)) * (y2[0]>0 ? 1./y2[0] : 0.);
	f22[0] = 1. * (double) ((group[0]==2)*(event[0]==2)) * (y2[0]>0 ? 1./y2[0] : 0.);
	r2 = (1.-0.)/1.*y2[0];
	f01[0] = (double) (event[0]==1) / (y1[0]+y2[0]);
	/* compute score etc. */
	i = 0;
	temp = r1*r2 * ((y1[i]+y2[i])>0 ? 1./(r1+r2) : 0.);
	temp1 = temp * (double)((group[i]==1)*(event[i]==1)) * (y1[i]>0 ? 1./r1 : 0.);
	temp2 = temp * (double)((group[i]==2)*(event[i]==1)) * (y2[i]>0 ? 1./r2 : 0.);
	temp0 = temp * (double)(event[i]==1) / (r1+r2);
	for (j=0; j<*d; j++) {
		score[j] += psi[i+j**n]*(temp2-temp1);
		if (pooled==0) { /* individual */
			h1[j] += psi[i+j**n]*temp1/(1.-f11[i]);
			h2[j] += psi[i+j**n]*temp2/(1.-f21[i]);
			q1[i+j**n] = psi[i+j**n]*temp/(1.-f11[i]) - h1[j];
			q2[i+j**n] = psi[i+j**n]*temp/(1.-f21[i]) - h2[j];
		} else { /* pooled */
			h1[j] += psi[i+j**n]*temp0/(1.-f01[i]);
			h2[j] += psi[i+j**n]*temp0/(1.-f01[i]);
			q1[i+j**n] = psi[i+j**n]*temp/(1.-f01[i]) - h1[j];
			q2[i+j**n] = psi[i+j**n]*temp/(1.-f01[i]) - h2[j];
		}
	}
	for (i=1; i<*n-1; i++) {
		s1[i+1] += s1[i];
		s1[i] = exp(-s1[i]);
		f11[i] = f11[i-1] + s1[i-1] * (double) ((group[i]==1)*(event[i]==1)) * (y1[i]>0 ? 1./y1[i] : 0.);
		f12[i] = f12[i-1] + s1[i-1] * (double) ((group[i]==1)*(event[i]==2)) * (y1[i]>0 ? 1./y1[i] : 0.);
		r1 = (1.-f11[i-1])/s1[i-1]*y1[i];
		s2[i+1] += s2[i];
		s2[i] = exp(-s2[i]);
		f21[i] = f21[i-1] + s2[i-1] * (double) ((group[i]==2)*(event[i]==1)) * (y2[i]>0 ? 1./y2[i] : 0.);
		f22[i] = f22[i-1] + s2[i-1] * (double) ((group[i]==2)*(event[i]==2)) * (y2[i]>0 ? 1./y2[i] : 0.);
		r2 = (1.-f21[i-1])/s2[i-1]*y2[i];
		f01[i] = f01[i-1] + (double) (event[i]==1) / (y1[i]/s1[i-1]+y2[i]/s2[i-1]);
		/* compute score increment etc. */
		temp = r1*r2 * ((y1[i]+y2[i])>0 ? 1./(r1+r2) : 0.);
		temp1 = temp * (double)((group[i]==1)*(event[i]==1)) * (y1[i]>0 ? 1./r1 : 0.);
		temp2 = temp * (double)((group[i]==2)*(event[i]==1)) * (y2[i]>0 ? 1./r2 : 0.);
		temp0 = temp * (double)(event[i]==1) / (r1+r2);
		for (j=0; j<*d; j++) {
			score[j] += psi[i+j**n]*(temp2-temp1);
			if (pooled==0) { /* individual */
				h1[j] += psi[i+j**n]*temp1/(1.-f11[i]);
				h2[j] += psi[i+j**n]*temp2/(1.-f21[i]);
				q1[i+j**n] = psi[i+j**n]*temp/(1.-f11[i]) - h1[j];
				q2[i+j**n] = psi[i+j**n]*temp/(1.-f21[i]) - h2[j];
			} else { /* pooled */
				h1[j] += psi[i+j**n]*temp0/(1.-f01[i]);
				h2[j] += psi[i+j**n]*temp0/(1.-f01[i]);
				q1[i+j**n] = psi[i+j**n]*temp/(1.-f01[i]) - h1[j];
				q2[i+j**n] = psi[i+j**n]*temp/(1.-f01[i]) - h2[j];
			}
		}
	}
	s1[*n-1] = exp(-s1[*n-1]);
	f11[*n-1] = f11[*n-2] + s1[*n-2] * (double) ((group[*n-1]==1)*(event[*n-1]==1)) * (y1[*n-1]>0 ? 1./y1[*n-1] : 0.);
	f12[*n-1] = f12[*n-2] + s1[*n-2] * (double) ((group[*n-1]==1)*(event[*n-1]==2)) * (y1[*n-1]>0 ? 1./y1[*n-1] : 0.);
	r1 = (1.-f11[*n-2])/s1[*n-2]*y1[*n-1];
	s2[*n-1] = exp(-s2[*n-1]);
	f21[*n-1] = f21[*n-2] + s2[*n-2] * (double) ((group[*n-1]==2)*(event[*n-1]==1)) * (y2[*n-1]>0 ? 1./y2[*n-1] : 0.);;
	f22[*n-1] = f22[*n-2] + s2[*n-2] * (double) ((group[*n-1]==2)*(event[*n-1]==2)) * (y2[*n-1]>0 ? 1./y2[*n-1] : 0.);;
	r2 = (1.-f21[*n-2])/s2[*n-2]*y2[*n-1];
	f01[*n-1] = f01[*n-2] + (double) (event[*n-1]==1) / (y1[*n-1]/s1[*n-2]+y2[*n-1]/s2[*n-2]);
	/* compute score increment etc. */
	i = *n-1;
	temp = r1*r2 * ((y1[i]+y2[i])>0 ? 1./(r1+r2) : 0.);
	temp1 = temp * (double)((group[i]==1)*(event[i]==1)) * (y1[i]>0 ? 1./r1 : 0.);
	temp2 = temp * (double)((group[i]==2)*(event[i]==1)) * (y2[i]>0 ? 1./r2 : 0.);
	temp0 = temp * (double)(event[i]==1) / (r1+r2);
	for (j=0; j<*d; j++) {
		score[j] += psi[i+j**n]*(temp2-temp1);
		if (pooled==0) { /* individual */
			h1[j] += psi[i+j**n]*temp1/(1.-f11[i]);
			h2[j] += psi[i+j**n]*temp2/(1.-f21[i]);
			q1[i+j**n] = psi[i+j**n]*temp/(1.-f11[i]) - h1[j];
			q2[i+j**n] = psi[i+j**n]*temp/(1.-f21[i]) - h2[j];
		} else { /* pooled */
			h1[j] += psi[i+j**n]*temp0/(1.-f01[i]);
			h2[j] += psi[i+j**n]*temp0/(1.-f01[i]);
			q1[i+j**n] = psi[i+j**n]*temp/(1.-f01[i]) - h1[j];
			q2[i+j**n] = psi[i+j**n]*temp/(1.-f01[i]) - h2[j];
		}
	}
	
	/* compute the variance matrix of the score (only the upper triangle) */
	/* compute covariance function rho_1 of f11 */
	if (pooled==0) { /* individual */
		/* compute the first sample's contribution to sigma */
		fj1_covariance(rho,n,1,event,group,f11,f12,y1);
		twosample_incidence_score_variance(q1,h1,rho,n,d,sigma,work);
		/* compute the second sample's contribution to sigma */
		fj1_covariance(rho,n,2,event,group,f21,f22,y2);
		twosample_incidence_score_variance(q2,h2,rho,n,d,sigma,work);
	} else { /* pooled */
		/* compute the first sample's contribution to sigma */
		fj1_covariance(rho,n,1,event,group,f01,f12,y1);
		twosample_incidence_score_variance(q1,h1,rho,n,d,sigma,work);
		/* compute the second sample's contribution to sigma */
	 	fj1_covariance(rho,n,2,event,group,f01,f22,y2);
		twosample_incidence_score_variance(q2,h2,rho,n,d,sigma,work);
	}
}

void twosample_incidence_score_variance(double *q, double *h, double *rho, int *n, int *d,
	double *sigma, double *intqrho)
{
	/* computes variance matrix of the score vector in Neyman's test */
	
	int i,j,k,kk;
	double temp;
	
	/* sigma is not initialised here; just add values to its current content */
	for (k=0; k<*d; k++)
		intqrho[k] = 0.;
	
	for (i=0; i<*n; i++) {
		temp = incr1sym(rho,*n,i);
		for (k=0; k<*d; k++)
			intqrho[k] += q[i+k**n]*temp;
		for (j=0; j<*n; j++) {
			temp = incr2sym(rho,*n,i,j);
			for (k=0; k<*d; k++) 
				for (kk=k; kk<*d; kk++)
					sigma[k+kk**d] += q[i+k**n]*q[j+kk**n]*temp;
		}
	}
	
	for (k=0; k<*d; k++) 
		for (kk=k; kk<*d; kk++) {
			sigma[k+kk**d] += + h[k]*intqrho[kk] + intqrho[k]*h[kk] + rho[(*n-1)+(*n-1)**n]*h[k]*h[kk];
		}
}

