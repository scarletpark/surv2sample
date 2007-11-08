#include <R.h>
#include <Rmath.h>
#include "proprate_estim.h"
#include "proprate_gs_weights.h"
#include "proprate_gs.h"

/* GILL--SCHUMACHER TYPE TEST for the two-sample proportional rate model */
/* originally for prop. haz. proposed by Gill & Schumacher (1987, Biometrika) */

void twosample_proprate_gs(int *event, int *group, int *n, int *model,
	int *weight1, int *weight2, double *beta1, double *beta2, double *stat, double *pval)
{
	/*
	weight1,2	possible values: 0 = logrank, 1 = Prentice--Wilcoxon (G^1), -1 = Gehan
	Gill & Schumacher recommend combination 0 and 1
	*/
	int i,j;
	int *y1,*y2; /* numbers at risk at each time */
	int n1,n2;
	double *k1,*k2; /* weights */
	double *da1,*da2; /* nelson--aalen estimators (or their weird-bootstrap resamples) */
	double *q1,*q2,*q1dot,*q2dot;
	double *b11,*b12,*b21,*b22;
	double rho11, rho12, rho21, rho22,var,v11,v12,v22,temp,temp2;
	
	y1 = (int *) R_alloc(2**n,sizeof(int));
	y2 = y1 + *n;
	k1 = (double *) R_alloc(12**n,sizeof(double));
	k2 = k1 + *n;
	da1 = k2 + *n;
	da2 = da1 + *n;
	q1 = da2 + *n;
	q2 = q1 + *n;
	q1dot = q2 + *n;
	q2dot = q1dot + *n;
	b11 = q2dot + *n;
	b12 = b11 + *n;
	b21 = b12 + *n;
	b22 = b21 + *n;
	
	double (*p_q)(),(*p_qdot)();
	if (*model==0) { /* prop. hazards */
		p_q = q_prophaz;
		p_qdot = qdot_prophaz;
	} else { /* prop. odds */
		p_q = q_propodds;
		p_qdot = qdot_propodds;
	}

	void (*p_weight1)(), (*p_weight2)();
	/* weight k1 */
	if (*weight1 == 0) { /* logrank */
		p_weight1 = logrank_weight_q;
	} else {
		if (*weight1 == -1) { /* Gehan */
			p_weight1 = gehan_weight_q;
		} else { /* Prentice--Wilcoxon (G^1) */
			p_weight1 = g1_weight_q;
		}
	}
	/* weight k2 */
	if (*weight2 == 0) { /* logrank */
		p_weight2 = logrank_weight_q;
	} else {
		if (*weight2 == -1) { /* Gehan */
			p_weight2 = gehan_weight_q;
		} else { /* Prentice--Wilcoxon (G^1) */
			p_weight2 = g1_weight_q;
		}
	}
	
	/* compute auxilliary processes */
	twosample_proprate_gs_aux(event, group, k1, k2, da1, da2, q1, q2, q1dot, q2dot,
		b11, b12, b21, b22, y1, y2, n, p_q, p_qdot, p_weight1, p_weight2);
	n1 = y1[0];
	n2 = y2[0];

	/* compute the test statistic */
	twosample_proprate_gs_stat(k1, k2, da1, da2, q1, q2, n, &rho11, &rho12, &rho21, &rho22, stat);
	/* and its variance */
	twosample_proprate_gs_var(k1, k2, da1, da2, q1, q2, q1dot, q2dot, b11, b12, b21, b22,
		y1, y2, n, &rho11, &rho12, &rho21, &rho22, &var);
	
	temp = *stat;
	*stat = *stat/sqrt(var);
	*pval = 2.*pnorm(fabs(*stat),0.,1.,0,0);
	*beta1 = rho12/rho11;
	*beta2 = rho22/rho21;
}

void twosample_proprate_gs_aux(int *event, int *group, double *k1, double *k2,
	double *da1, double *da2, double *q1, double *q2, double *q1dot, double *q2dot,
	double *b11, double *b12, double *b21, double *b22, int *y1, int *y2, int *n,
	double (*p_q)(), double (*p_qdot)(), void (*p_weight1)(), void (*p_weight2)())
{
	/* computes auxilliary processes for the GS test: y1,y2,k1,k2,da1,da2,q1,q2 */
	int i;
	double temp1,temp2;
	
	/* numbers at risk */
	y1[*n-1] = (group[*n-1]==1);
	y2[*n-1] = (group[*n-1]==2);
	for (i=*n-2; i>=0; i--) {
		y1[i] = y1[i+1] + (group[i]==1);
		y2[i] = y2[i+1] + (group[i]==2);
	}
	
	/* increments of nelson--aalen estimators da1, da2, and values q(a1), q(a2) */
	temp1 = temp2 = 0.;
	for (i=0; i<*n; i++) {
		if ((event[i]==1)&&(group[i]==1)&&(y1[i]>0)) {
			da1[i] = 1./y1[i];
		} else {
			da1[i] = 0.;
		}
		if ((event[i]==1)&&(group[i]==2)&&(y2[i]>0)) {
			da2[i] = 1./y2[i];
		} else {
			da2[i] = 0.;
		}
		q1[i] = (*p_q)(temp1);
		q2[i] = (*p_q)(temp2);
		q1dot[i] = (*p_qdot)(temp1);
		q2dot[i] = (*p_qdot)(temp2);
		temp1 += da1[i];
		temp2 += da2[i];
	}
	
	/* weights k1, k2 */
	/* weights are unnormalised, q1 is also used as a dummy argument in place of time and is not altered */
	p_weight1(y1, y2, event, q1, q1, q2, n, 0, k1);
	p_weight2(y1, y2, event, q1, q1, q2, n, 0, k2);
	
	i = 0;
	b11[i] = -k1[i]*q1dot[i]*da1[i]/(q1[i]*q1[i]);
	b12[i] = -k1[i]*q2dot[i]*da2[i]/(q2[i]*q2[i]);
	b21[i] = -k2[i]*q1dot[i]*da1[i]/(q1[i]*q1[i]);
	b22[i] = -k2[i]*q2dot[i]*da2[i]/(q2[i]*q2[i]);
	for (i=1; i<*n; i++) {
		b11[i] = b11[i-1] + (-k1[i]*q1dot[i]*da1[i]/(q1[i]*q1[i]));
		b12[i] = b12[i-1] + (-k1[i]*q2dot[i]*da2[i]/(q2[i]*q2[i]));
		b21[i] = b21[i-1] + (-k2[i]*q1dot[i]*da1[i]/(q1[i]*q1[i]));
		b22[i] = b22[i-1] + (-k2[i]*q2dot[i]*da2[i]/(q2[i]*q2[i]));
	}
}

void twosample_proprate_gs_stat(double *k1, double *k2, double *da1, double *da2, double *q1, double *q2,
	int *n, double *rho11, double *rho12, double *rho21, double *rho22, double *stat)
{
	int i;
	double temp1 = 0., tempq1, temp2 = 0., tempq2;
	
	double a1 = 0., a1prev = 0., a2 = 0., a2prev = 0.;
	
	*rho11 = *rho12 = *rho21 = *rho22 = 0.;
	for (i=0; i<*n; i++) {
		*rho11 += k1[i]*da1[i]/q1[i];
		*rho12 += k1[i]*da2[i]/q2[i];
		*rho21 += k2[i]*da1[i]/q1[i];
		*rho22 += k2[i]*da2[i]/q2[i];
	}
	*stat = *rho22**rho11-*rho12**rho21;
}

void twosample_proprate_gs_var(double *k1, double *k2, double *da1, double *da2,
	double *q1, double *q2, double *q1dot, double *q2dot,
	double *b11, double *b12, double *b21, double *b22, int *y1, int *y2, int *n,
	double *rho11, double *rho12, double *rho21, double *rho22, double *var)
{
	int i;
	double theta_mean, temp, temp1, temp2;
	
	/* analog of Gill & Schumacher's variance estimator; may be and sometimes is negative */
	/* DANGEROUS, don't use */
	/*
	double v11, v12, v22, temp1, temp2;
	double temp111, temp112, temp221, temp222, temp121, temp122;
	v11 = v12 = v22 = 0.;
	for (i=0; i<*n; i++) {
		temp1 = ( (y2[i]>0) ? q2[i]*da1[i]/(q1[i]*y2[i]) : 0. );
		temp2 = ( (y1[i]>0) ? q1[i]*da2[i]/(q2[i]*y1[i]) : 0. );
		temp111 = (k1[i]/q1[i]-b11[i]+b11[*n-1])*(k1[i]/q1[i]-b11[i]+b11[*n-1]);
		temp112 = (k1[i]/q2[i]-b12[i]+b12[*n-1])*(k1[i]/q2[i]-b12[i]+b12[*n-1]);
		temp221 = (k2[i]/q1[i]-b21[i]+b21[*n-1])*(k2[i]/q1[i]-b21[i]+b21[*n-1]);
		temp222 = (k2[i]/q2[i]-b22[i]+b22[*n-1])*(k2[i]/q2[i]-b22[i]+b22[*n-1]);
		temp121 = (k1[i]/q1[i]-b11[i]+b11[*n-1])*(k2[i]/q1[i]-b21[i]+b21[*n-1]);
		temp122 = (k1[i]/q2[i]-b12[i]+b12[*n-1])*(k2[i]/q2[i]-b22[i]+b22[*n-1]);
		v11 += temp111*temp2 + temp112*temp1;
		v12 += temp121*temp2 + temp122*temp1;
		v22 += temp221*temp2 + temp222*temp1;
	}
	*var = *rho21**rho22*v11 - (*rho12**rho21+*rho11**rho22)*v12 + *rho11**rho12*v22;
	*/
	
	/* this variance estimator is always positive (sum of squares) */
	/* this estimator is different from that of Gill & Schumacher (1987), */
	/* which may be and sometimes is negative */
	*var = 0.;
	theta_mean = .5*(*rho12/(*rho11)+*rho22/(*rho21));
	temp = *rho21/(*rho11);
	for (i=0; i<*n; i++) {
		temp1 = temp*(k1[i]/q1[i]-b11[i]+b11[*n-1]) - (k2[i]/q1[i]-b21[i]+b21[*n-1]); /* rho21/rho11*r11(t) - r21(t) */
		temp2 = (k2[i]/q2[i]-b22[i]+b22[*n-1]) - temp*(k1[i]/q2[i]-b12[i]+b12[*n-1]); /* r22(t) - rho21/rho11*r12(t) */
		*var += ( (y1[i]>0) ? temp1*temp1*theta_mean*q1[i]*da2[i]/(q2[i]*y1[i]) : 0. )
			+ ( (y2[i]>0) ? temp2*temp2*theta_mean*q2[i]*da1[i]/(q1[i]*y2[i]) : 0. );
	}
	*var *= *rho11**rho11;
}

