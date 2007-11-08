#include <R.h>
#include <Rmath.h>
#include "proprate_gs_weights.h"

/* WEIGHTS FOR GILL--SCHUMACHER TESTS for proportional rate transformation models */

/* for proportional rates weights needn't be normalised */
/* (time and normalisation would be necessary only for additive models, not implemented) */

/* logrank weight */
void logrank_weight_q(int *y1, int *y2, int *event, double *time, double *q1, double *q2,
	int *n, int normalised, double *weight)
{
	int i;
	double integral,t0,dt;
	
	if (normalised==0) {
		for (i=0; i<*n; i++)
			weight[i] = (double) q1[i]*y1[i]*q2[i]*y2[i]/(q1[i]*y1[i]+q2[i]*y2[i]);
	} else {
		integral = 0.;
		t0 = 0.;
		for (i=0; i<*n; i++) {
			dt = time[i]-t0;
			t0 = time[i];
			weight[i] = (double) q1[i]*y1[i]*q2[i]*y2[i]/(q1[i]*y1[i]+q2[i]*y2[i]);
			integral += weight[i]*dt;
		}
		for (i=0; i<*n; i++)
			weight[i] /= integral;
	}
}

/* Gehan weight */
void gehan_weight_q(int *y1, int *y2, int *event, double *time, double *q1, double *q2,
	int *n, int normalised, double *weight)
{
	int i;
	double integral,t0,dt;
	
	if (normalised==0) {
		for (i=0; i<*n; i++)
			weight[i] = (double) q1[i]*y1[i]*q2[i]*y2[i];
	} else {
		integral = 0.;
		t0 = 0.;
		for (i=0; i<*n; i++) {
			dt = time[i]-t0;
			t0 = time[i];
			weight[i] = (double) q1[i]*y1[i]*q2[i]*y2[i];
			integral += weight[i]*dt;
		}
		for (i=0; i<*n; i++)
			weight[i] /= integral;
	}
}

/* Prentice--Wilcoxon G^1 weight */
void g1_weight_q(int *y1, int *y2, int *event, double *time, double *q1, double *q2,
	int *n, int normalised, double *weight)
{
	int i;
	double integral,t0,dt,nelsonaalen;
	
	if (normalised==0) {
		nelsonaalen = 0.;
		for (i=0; i<*n; i++) {
			weight[i] =  (double) q1[i]*y1[i]*q2[i]*y2[i]/(q1[i]*y1[i]+q2[i]*y2[i])*exp(-nelsonaalen);
			nelsonaalen += (double) event[i]/(y1[i]+y2[i]);
		}
	} else {
		nelsonaalen = 0.;
		integral = 0.;
		t0 = 0.;
		for (i=0; i<*n; i++) {
			dt = time[i]-t0;
			t0 = time[i];
			weight[i] = (double) q1[i]*y1[i]*q2[i]*y2[i]/(q1[i]*y1[i]+q2[i]*y2[i])*exp(-nelsonaalen);
			nelsonaalen += (double) event[i]/(y1[i]+y2[i]);
			integral += weight[i]*dt;
		}
		for (i=0; i<*n; i++)
			weight[i] /= integral;
	}
}

