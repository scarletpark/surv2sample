#include <R.h>
#include <Rmath.h>
#include "proprate_estim.h"

/* ESTIMATION PROCEDURE for two-sample transformation models */

void twosample_proprate_estim_r_call(int *event, int *group, int *n, double *beta, int *maxiter,
	double *eps, int *flag, int *model, double *q1, double *q2,
	double *q1dot, double *q2dot, double *u1_process, double *u1, double *d11_process,
	double *d11, double *sigma11, double *logliks)
{
	/* interface for calling twosample_proprate_estim from R */
	/* because of pointers p_q, p_qdot to functions */
	/* TODO: user defined q, qdot */
	
	double (*p_q)(),(*p_qdot)();
	if (*model==0) { /* prop. hazards */
		p_q = q_prophaz;
		p_qdot = qdot_prophaz;
	} else { /* prop. odds */
		p_q = q_propodds;
		p_qdot = qdot_propodds;
	}
	
	twosample_proprate_estim(event, group, n, beta, maxiter, eps, flag, p_q, p_qdot, q1, q2, q1dot, q2dot,
		u1_process, u1, d11_process, d11, sigma11, logliks);
}

void twosample_proprate_estim(int *event, int *group, int *n, double *beta, int *maxiter,
	double *eps, int *flag, double (*p_q)(), double (*p_qdot)(), double *q1, double *q2,
	double *q1dot, double *q2dot, double *u1_process, double *u1, double *d11_process,
	double *d11, double *sigma11, double *logliks)
{
	/*
	simplified partial likelihood estimation procedure
	score equation solved by the Newton--Raphson method with step-halving
	(similar to agfit3 in agfit3.c in package survival)
	beta on input is the init. value for iteration, on output the estimate
	logliks: [0] at beta init, [1] at the estimate
	*/
	int i,j;
	int y1,y2,y1work,y2work;
	int iter,halving;
	double temp1,temp2;
	
	double loglik,logliknew,betanew;
	
	y1 = 0;
	y2 = 0;
	for (i=0; i<*n; i++) {
		if (group[i] == 1) {
			++y1;
		} else {
			++y2;
		}
	}
	y1work = y1;
	y2work = y2;
	temp1 = temp2 = 0.; /* this will now be left-cont. nelson--aalens for the two samples */
	for (i=0; i<*n; i++) {
		q1[i] = (*p_q)(temp1);
		q2[i] = (*p_q)(temp2);
		q1dot[i] = (*p_qdot)(temp1);
		q2dot[i] = (*p_qdot)(temp2);
		if (group[i]==1) {
			temp1 += (double) (event[i])/y1work;
			--y1work;
		} else {
			temp2 += (double) (event[i])/y2work;
			--y2work;
		}
	}

	betanew = *beta; /* initial value */
	for (iter=0; iter<*maxiter; iter++) {
		twosample_proprate_score(event, group, n, &betanew, q1, q2, q1dot, q2dot, u1_process, u1, d11_process, d11, sigma11, &logliknew);
		if (iter==0) { /* in iteration 0 we don't test for convergence */
			loglik = logliknew;
			logliks[0] = loglik;
			betanew += 1./(*d11)**u1;
		} else { /* not in iteration 0, we test for convergence */
			if (fabs(1-(loglik/logliknew))<=*eps ) { /* all done, converged */
				loglik = logliknew;
				logliks[1] = loglik;
				*beta = betanew;
				if (halving==1) *flag = 1000; /* didn't converge after all */
				*maxiter = iter;
				return;
			}
			if (iter==*maxiter) break;  /* didn't converge after maxiter; skip the step halving and etc */
			if (logliknew < loglik) { /* it is not converging; betanew doesn't improve loglik */
				halving = 1;
				betanew = (betanew+*beta)/2.; /* half of old increment (go to the middle between old and new) */
			} else { /* OK, accept betanew as beta and propose next betanew (newton step) */
				halving = 0;
				loglik = logliknew;
				*beta = betanew;
				betanew += 1./(*d11)**u1;
			}
		}
	}
	
	/* this is done only if iter is maxiter (didn't converge) */
	loglik = logliknew;
	*beta = betanew;
	*flag = 1000;
}

void twosample_proprate_score(int *event, int *group, int *n, double *beta, double *q1, double *q2,
	double *q1dot, double *q2dot, double *u1_process, double *u1, double *d11_process, double *d11,
	double *sigma11, double *loglik)
{
	/* computes increments of the score process (u1_process) and the score vector (u1) at beta (U1 in the paper) */
	/* the same for its derivative (D11 in the paper) */
	/* loglik is the log of the modified partial likelihood */
	/* sigma11 is the cumulative sigma11 at the end (time tau) */
	/* q1, q2 are q(u)=lambda0(Lambda0^{-1}(u)) evaluated at Nelson--Aalen estimators for samples 1, 2;
	/* q(u)=1 for prop. haz., q(u)=exp(-u) for prop. odds */
	int i,j;
	int y1work = 0, y2work = 0;
	double theta = exp(*beta), temp, h1 = 0., h2 = 0., r1 = 0., r2 = 0.;
	
	*u1 = *d11 = *sigma11 = *loglik = 0.;
	y1work = y2work = 0;
	for (i=*n-1; i>=0; i--) { /* loop over times; from last to first because r1,r2 in h1,h2 are cumulative from t to tau */
		if(group[i]==1) {
			++y1work;
		} else {
			++y2work;
		}
		if (event[i]==1) {
			temp = 1./(q1[i]*y1work+theta*q2[i]*y2work);
			u1_process[i] = temp*((group[i]==2)*q1[i]*y1work - (group[i]==1)*theta*q2[i]*y2work);
			*u1 += u1_process[i];
			d11_process[i] = q1[i]*y1work*theta*q2[i]*y2work*temp*temp;
			*d11 += d11_process[i];
			*loglik += *beta*(group[i]==2) + log(temp);
			r1 += -y1work*q1dot[i]*theta*y2work*q2[i]*temp*temp;
			r2 += y1work*q1[i]*theta*y2work*q2dot[i]*temp*temp;
			h1 = -theta*q2[i]*y2work*temp - (y1work>0 ? 1./y1work : 0.) * r1;
			h2 = q1[i]*y1work*temp - (y2work>0 ? 1./y2work : 0.) * r2;
			*sigma11 += (h1*h1*y1work*q1[i]+h2*h2*theta*y2work*q2[i])*temp;
		} else {
			u1_process[i] = 0.;
			d11_process[i] = 0.;
		}
	}
}

/* DEFINITION of q_\omega and its derivative for MAIN MODELS */

/* proportional hazards */

double q_prophaz(double t)
{
	return 1.;
}

double qdot_prophaz(double t)
{
	return 0.;
}

/* proportional odds */

double q_propodds(double t)
{
	return exp(-t);
}

double qdot_propodds(double t)
{
	return -exp(-t);
}


