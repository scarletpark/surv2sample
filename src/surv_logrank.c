#include <R.h>
#include <Rmath.h>
#include "my_cholesky.h"
#include "surv_logrank.h"

/* TWO-SAMPLE G^{\rho,\gamma} WEIGHTED LOGRANK TESTS AND THEIR COMBINATIONS (sum, max) */

void twosample_grhogamma(double *time, int *event, int *group, int *n, double *rho, double *gamma,
	double *b, int *m, int *nsim_asympt, int *nperm, int *nboot, double *choltol,
	double *stats, double *sigma, double *pvals, double *pvals_perm, double *pvals_boot,
	double *stat_sum, double *pval_sum, double *pval_sum_perm, double *pval_sum_boot,
	double *stat_max, double *pval_max, double *pval_max_perm, double *pval_max_boot)
{
	/* computes G^{\rho,\gamma} tests rho and gamma (rho, gamma are vectors of length m, */
	/* hence m statistics are computed at once) */
	/* further computes max and b-weighted sum of them (sum of b[j]*fabs(stats[j]) */
	/* p-values approximated by asymptotics, permutations, or bootstrap */
	
	int i,j,k,r;
	int y1,y2,y1bak,y2bak;
	double temp1,temp2,temp3;
	double *incr; /* work array to store increments of weighted logranks and other things */
	double *stats_work;
	double *sigma_work;
	int *eventg; /* bootstrap sample */
	int *freqg; /* bootstrap sample frequencies */
	
	GetRNGstate();
	
	incr = (double *) R_alloc(*m+*m+*m**m,sizeof(double));
	stats_work = incr+*m;
	sigma_work = stats_work+*m;
	
	
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
	
	/* OBSERVED test statistics */
	grhogamma_stat(event, group, n, rho, gamma, b, m, &y1, &y2, incr, stats, stat_sum, stat_max, sigma);
	y1 = y1bak;
	y2 = y2bak;
	
	/* BOOTSTRAP */
	/* i.e., sampling from the data with replacement */
	for (j=0; j<*m; j++)
		pvals_boot[j] = 0.;
	*pval_sum_boot = *pval_max_boot = 0.;
	if (*nboot>0) {
		/* allocate arrays for the bootstrap sample */
		eventg = (int *) R_alloc(2**n,sizeof(int));
		freqg = eventg + *n;
		for (i=0; i<*nboot; i++) {
			/* compute frequencies of the original data in the bootstrap sample (we want to avoid sorting the bootstrap sample) */
			/* frequencies will be temporarily stored in eventg */
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
			/* compute the statistics from the PERMUTED data */
			grhogamma_stat(eventg, group, n, rho, gamma, b, m, &y1, &y2, incr, stats_work, &temp2, &temp3, sigma_work);
			y1 = y1bak;
			y2 = y2bak;
			for (j=0; j<*m; j++)
				/* pvals_boot[j] += (stats_work[j] > stats[j]); */
				pvals_boot[j] += (fabs(stats_work[j]) > fabs(stats[j]));
			*pval_sum_boot += (temp2 > *stat_sum);
			*pval_max_boot += (temp3 > *stat_max);
		}
		for (j=0; j<*m; j++)
			pvals_boot[j] /= *nboot;
		*pval_sum_boot /= *nboot;
		*pval_max_boot /= *nboot;
	}
	
	/* PERMUTATION test */
	/* i.e., sampling from the data without replacement */
	if (*nperm>0) {
		for (j=0; j<*m; j++)
			pvals_perm[j] = 0.;
		*pval_sum_perm = *pval_max_perm = 0.;
		for (i=0; i<*nperm; i++) {
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
			/* compute the statistics from the PERMUTED data */
			grhogamma_stat(event, group, n, rho, gamma, b, m, &y1, &y2, incr, stats_work, &temp2, &temp3, sigma_work);
			y1 = y1bak;
			y2 = y2bak;
			for (j=0; j<*m; j++)
				/* pvals_perm[j] += (stats_work[j] > stats[j]); */
				pvals_perm[j] += (fabs(stats_work[j]) > fabs(stats[j]));
			*pval_sum_perm += (temp2 > *stat_sum);
			*pval_max_perm += (temp3 > *stat_max);
		}
		for (j=0; j<*m; j++)
			pvals_perm[j] /= *nperm;
		*pval_sum_perm /= *nperm;
		*pval_max_perm /= *nperm;
	}
	
	/* standardise observed variance matrix (make cor from cov) */
	for (j=0; j<*m; j++) {
		incr[j] = sqrt(sigma[j+j**m]); /* use incr as storage for sqrt of diagonal of sigma */
		sigma[j+j**m] = 1.;
		for (i=0; i<j; i++)
			sigma[i+j**m] /= (incr[i]*incr[j]); /* make cov from cor */
	}
	
	/* ASYMPTOTIC P-VALUES */
	
	/* asymptotic p-values of individual weighted logrank tests */
	/* standardise individual test statistics, compute combined stats, correlation matrix from covariance matrix etc. */
	for (j=0; j<*m; j++) {
		pvals[j] = 2.*pnorm(fabs(stats[j]),0.,1.,0,0);
	}

	/* asymptotic p-values of sum and max */
	/* now approximate p-values of sum and max stat by simulation from the limiting distribution */
	/* cholesky of sigma is stored in the lower triangle including diagonal, the upper triangle
	is unchanged; diag of the original sigma is 1 (cor matrix), hence sigma can be reconstructed */
	my_cholesky2(sigma,m,choltol);
	/* cholesky of sigma is F*D*F^T, D on diag, F in lower triangle
	now change it into E*E^T, i.e., E = F*D^(1/2) (spare multiplications in simulation loop) */
	for (j=0; j<*m; j++) {
		sigma[j+j**m] = sqrt(sigma[j+j**m]);
		for (i=j+1; i<*m; i++)
			sigma[i+j**m] *= sigma[j+j**m];
	}
	*pval_sum = *pval_max = 0.;
	if (*nsim_asympt>0) {
		for (i=0; i<*nsim_asympt; i++) {
			temp1 = 0.; /* sum */
			temp2 = 0.; /* max */
			for (j=0; j<*m; j++) {
				incr[j] = norm_rand(); /* use incr as work array for simulated normal vector with cov = sigma */
			}
			for (j=0; j<*m; j++) {
				temp3 = 0.; /* j-th line of E times incr */
				for (k=0; k<=j; k++)
					temp3 += sigma[j+k**m]*incr[k];
				temp1 += b[j]*fabs(temp3);
				if (fabs(temp3) > temp2)
					temp2 = fabs(temp3);
			}
			*pval_sum += (*stat_sum < temp1);
			*pval_max += (*stat_max < temp2);
		}
		*pval_sum /= *nsim_asympt;
		*pval_max /= *nsim_asympt;
	}
		
	/* finally reconstruct sigma */
	for (j=0; j<*m; j++) {
		sigma[j+j**m] = 1.;
		for (k=0; k<j; k++)
			sigma[j+k**m] = sigma[k+j**m];
	}
	
	PutRNGstate();
	
}

void grhogamma_stat(int *event, int *group, int *n, double *rho, double *gamma,
	double *b, int *m, int *y1, int *y2, double *incr, double *stats, double *stat_sum,
	double *stat_max, double *sigma)
{
	/* computes test statistics */
	/* stats contains G(rho,gamma) statistics (signed, not abs) */
	/* stat_sum is weighted sum of abs(stats), max their maximum */
	/* sigma is cov matrix of stats */
	
	double temp1;
	int i,j,k;
	
	for (j=0; j<*m; j++) {
		stats[j] = 0.;
		for (k=0; k<*m; k++)
			sigma[j+k**m] = 0.;
	}
	
	/* compute the weighted logrank test statistic(s) and their covariance */
	temp1 = 0.; /* this will be current value of left-continuous nelson--aalen (computed on the fly, not in advance) */
	/* instead of left-continuous kaplan--meier we use exp of minus left-cont. nelson--aalen */
	for (i=0; i<*n; i++) {
		if (event[i] == 1) { /* event, hence u and sigma have nonzero increments */
			if (group[i] == 1) {
				for (j=*m-1; j>=0; j--) { /* we'll compute upper triangle of sigma, that's why reversed order */
					incr[j] = - exp(-rho[j]*temp1)*pow(1.-exp(-temp1),gamma[j]) / (*y1+*y2);
					stats[j] += incr[j]**y2;
					for (k=j; k<*m; k++)
						sigma[j+k**m] += incr[j]*incr[k]**y1**y2;
				}
			} else {
				for (j=*m-1; j>=0; j--) { /* we'll compute upper triangle of sigma, that's why reversed order */
					incr[j] = + exp(-rho[j]*temp1)*pow(1.-exp(-temp1),gamma[j]) / (*y1+*y2);
					stats[j] += incr[j]**y1;
					for (k=j; k<*m; k++)
						sigma[j+k**m] += incr[j]*incr[k]**y1**y2;
				}
			}
			temp1 += 1./(*n-i); /* nelson--aalen for the next step, thus left continuous */
		}
		if (group[i] == 1) {
			--*y1;
		} else {
			--*y2;
		}
	}
	
	/* compute test statistics: standardised individual and combined test statistics */
	*stat_sum = 0.;
	*stat_max = 0.;
	for (j=0; j<*m; j++) {
		stats[j] = stats[j]/sqrt(sigma[j+j**m]);
		*stat_sum += b[j]*fabs(stats[j]); /* sum combination statistic */
		if (fabs(stats[j]) > *stat_max)
			*stat_max = fabs(stats[j]); /* max combination statistic */
	}

}
