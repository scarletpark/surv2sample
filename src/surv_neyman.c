#include <R.h>
#include <Rmath.h>
#include "util_neyman.h"
#include "util_datadriven.h"
#include "surv_neyman.h"

/* Neyman's smooth tests of euqality of survival distributions in two samples */

void twosample_neyman(double *time, int *event, int *group, int *n, int *d, int *d0,
	int *timetransf, int *basis, int *datadriven, double *penalty, int *emptyset,
	double *choltol, int *nsim_asympt, int *nperm, int *nboot, double *score, double *sigma,
	double *stat_d, double *pval_d, double *pval_d_perm, double *pval_d_boot,
	double *stat_nested, double *pval_nested_asympt, double *pval_nested_h,
	double *pval_nested_perm, double *pval_nested_boot,
	double *stats_nested, double *stats_penal_nested, int *s_nested,
	double *stat_all, double *pval_all_asympt, double *pval_all_perm, double *pval_all_boot,
	double *stats_all, double *stats_penal_all, int *s_all, int *all_subsets)
{
	/*
	time		vector of sorted times (ascending), length n
	event		event indicator (0=censored, 1=death), length n
	group		group 1 or 2, length n
	n			number of observations
	d			number of basis functions for neyman's test
	d0			minimum number of basis functions for nested subsets (i.e., 1,...,d0 always included)
	timetransf	time transformation in basis functions (1 = distrib fnc, 2 = cumul hazard, 3 = identity)
	basis		basis of smooth functions (1 = legendre, 2 = cosine, 3 = indicators?)
	datadriven	should a data-driven test be performed? (0=no, 1=nested subsets, 2=allsubsets, 3=both)
	choltol		cholesky tolerance
	nsim		number of simulations
	score		score vector, length d, must be allocated; output
	sigma		variance matrix of score, d by d, must be allocated; output
	stat_d		chi-square test statistic score^T * sigma^{-1} * score (test with fixed dimension d); output
	pval_d		p-value; output
	stat_nested	nested-subsets BIC data-driven test statistic; output
	pval_nested_asympt	p-value based on asymptotic chisq_{max(d0,1)}, works only theoretically, for small to moderate samples not; output
	pval_nested_h	p-val based on the two-term h-approximation, works correctly; output
	stats_nested	score statistics for nested subsets, length d; output
	stats_penal_nested	penalised score statistics, length d; output
	s_nested	selected subset; output
	stat_all	all-subsets BIC data-driven test statistic; output
	pval_all_asympt		p-value based on asymptotic max of chisq_1 (for d0=0, simulated) or chisq_d0 (d0>0), works correctly; output
	pval_all_maxchisq1_indep	[not used] only for d0=0, p-val based on asymptotic max(chisq_1) as if they were indep., works for timetransf sigma
	 				(orthogonality), sometimes works correctly even for other but don't rely on this; output
	pval_all_h	[not used] p-value based on the two-term approximation (with simulation); output
	stats_all	score statistics for nested subsets, length d; output
	stats_penal_all	penalised score statistics, length d; output
	s_all		selected subset; output
	all_subsets	all 2^(d-d0) (for d0>0) or 2^d-1 (d0=0) subsets: each row consists of d 0/1 indicators); output
	*/
	
	int y1,y2,y1bak,y2bak;
	int i,j,k,r,kk;
	double temp1,temp2,temp3,temp4;
	double logn = log((double) *n);
	double *a; /* transformed times for basis functions */
	double *psi; /* vector of values of basis fncs at the i-th (transformed) time, and later work arr */
	void (*p_basis)(); /* pointer to the function evaluating the basis (legendre_basis, cosine_basis or ...) */
	double *sigmawork; /* work array for sigma, allocated only if data-driven test is to be done */
	int *subsets;  /* work array giving a subset of {1,...,d} and its size, for all-subsets BIC, dimension (2^d-1) x (d+1) */
	double *scorework; /* work array for subsets of score, for all-subsets BIC */
	double *scoreg; /* LWY-simulated or permutation or bootstrap score vector */
	double *sigmag; /* corresponding var matrix */
	double *stats_work; /* work array for score statistics in LWY simulations */
	double *stats_penal_work; /* work array for penalised score statistics in LWY simulations */
	double *timeg; /* times in the bootstrap sample */
	int *eventg; /* event indicators in the bootstrap sample */
	int *freqg; /* frequencies of the observed data in the bootstrap sample */
	double *sigmad01; /* var matrix for (d0+1)-dimensional statistics, used in the two-term approx. with all subsets */
	
	GetRNGstate();
	
	a = (double *) R_alloc(*n+*d+*d+*d**d,sizeof(double));
	psi = a+*n;
	scoreg = psi+*d;
	sigmag = scoreg+*d;
	
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
	
	/* allocate some work arrays */
	if (*datadriven != 0) /* any data-driven */
		sigmawork = (double *) R_alloc(*d*(*d),sizeof(double));
	if ( *datadriven == 1 ) { /* nested */
		/* work arrays of length kk=d-max(d0,1)+1 for simulation, permutations,... */
		kk = *d - (*d0>1 ? *d0 : 1) + 1;
		stats_work = (double *) R_alloc(2*kk,sizeof(double));
		stats_penal_work = stats_work + kk;
	}
	if ( (*datadriven == 2) || (*datadriven == 3)) { /* all subsets or both */
		if (*d0 == 0) {
			kk = (int) pow(2,*d)-1; /* save number of all subsets */
		} else {
			kk = (int) pow(2,*d-*d0);  /* save number of all subsets */
		}
		scorework = (double *) R_alloc(*d+2*kk,sizeof(double));
		stats_work = scorework + *d;
		stats_penal_work = stats_work + kk;
		subsets = (int *) R_alloc(kk*(*d+1),sizeof(int)); /* i-th row contains indexes of elements of i-th subset and its size */
		if (*d0 == 0) {
			for (i=1; i<=kk; i++) { /* that is i goes from 1 to 2^d-1 */
				k = 0; /* will be # of elements of the i-th subset now */
				for (j=0; j<*d; j++) {
					/* if j is in the i-th subset, add j to subset and increase the counter kk */
					if ( ((((unsigned) i) >> j) & 1) == 1 ) { /* tricky: j-th bit of i from the right */
						subsets[(i-1)+k*kk] = j;
						++k;
					}
					all_subsets[(i-1)+j*kk] = ((((unsigned) i) >> j) & 1); /* return all subsets */
				}
				subsets[(i-1)+(*d)*kk] = k; /* the last entry on the i-th row is the number of elements of the i-th subset */
			}
		} else {
			for (i=0; i<=kk-1; i++) { /* that is i goes from 0 to 2^(d-d0)-1 */
				for (j=0; j<*d0; j++) { /* include first d0 functions */
					subsets[i+j*kk] = j;
					all_subsets[i+j*kk] = 1;
				}
				k = *d0; /* will be # of elements of the i-th subset now */
				for (j=0; j<*d-*d0; j++) {
					/* if j is in the i-th subset, add j to subset and increase the counter kk */
					if ( ((((unsigned) i) >> j) & 1) == 1 ) { /* tricky: j-th bit of i from the right */
						subsets[i+k*kk] = *d0+j;
						++k;
					}
					all_subsets[i+(*d0+j)*kk] = ((((unsigned) i) >> j) & 1); /* return all subsets */
				}
				subsets[i+(*d)*kk] = k; /* the last entry on the i-th row is the number of elements of the i-th subset */
			}
		}
	}
	
	
	/* compute the score vector from the OBSERVED data and save it in score */
	neyman_score(score, sigma, d, time, event, group, n, &y1, &y2, p_basis, timetransf, a, psi);
	y1 = y1bak;
	y2 = y2bak;
	/* compute the score statistics (fixed-dim. or data-driven) from the OBSERVED score */
	neyman_stat(score,scorework,sigma,sigmawork,d,d0,stat_d,stat_nested,stats_nested,stats_penal_nested,s_nested,
		stat_all,stats_all,stats_penal_all,s_all,subsets,&kk,penalty,psi,choltol,datadriven,emptyset);
	
	
	/* BOOTSTRAP */
	/* i.e., sampling from the data with replacement */
	*pval_d_boot = *pval_nested_boot = *pval_all_boot = 0.;
	if (*nboot>0) {
		/* allocate arrays for the bootstrap sample */
		timeg = (double *) R_alloc(*n,sizeof(double));
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
					timeg[r] = time[j];
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
			/* compute the score vector from the BOOTSTRAP sample and save it in scoreg */
			neyman_score(scoreg, sigmag, d, timeg, eventg, group, n, &y1, &y2, p_basis, timetransf, a, psi);
			y1 = y1bak;
			y2 = y2bak;
			/* compute the score statistics (fixed-dim. or data-driven) from the BOOTSTRAP score */
			neyman_stat(scoreg,scorework,sigmag,sigmawork,d,d0,&temp1,&temp2,stats_work,stats_penal_work,&j,
				    &temp3,stats_work,stats_penal_work,&j,subsets,&kk,penalty,psi,choltol,datadriven,emptyset);
			*pval_d_boot += (temp1 > *stat_d);
			*pval_nested_boot += (temp2 > *stat_nested);
			*pval_all_boot += (temp3 > *stat_all);
		}
		*pval_d_boot /= *nboot;
		*pval_nested_boot /= *nboot;
		*pval_all_boot /= *nboot;
	}
	
	/* PERMUTATION test */
	/* i.e., sampling from the data without replacement */
	*pval_d_perm = *pval_nested_perm = *pval_all_perm = 0.;
	if (*nperm>0) {
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
			/* compute the score vector from the PERMUTED data and save it in scoreg */
			neyman_score(scoreg, sigmag, d, time, event, group, n, &y1, &y2, p_basis, timetransf, a, psi);
			y1 = y1bak;
			y2 = y2bak;
			/* compute the score statistics (fixed-dim. or data-driven) from the PERMUTATION score */
			neyman_stat(scoreg,scorework,sigmag,sigmawork,d,d0,&temp1,&temp2,stats_work,stats_penal_work,&j,
				    &temp3,stats_work,stats_penal_work,&j,subsets,&kk,penalty,psi,choltol,datadriven,emptyset);
			*pval_d_perm += (temp1 > *stat_d);
			*pval_nested_perm += (temp2 > *stat_nested);
			*pval_all_perm += (temp3 > *stat_all);
		}
		*pval_d_perm /= *nperm;
		*pval_nested_perm /= *nperm;
		*pval_all_perm /= *nperm;
	}	
	
	
	/* ASYMPTOTIC p-values */
	
	/* asymptotic p-values of fixed-dimensional test */
	*pval_d = pchisq(*stat_d,*d,0,0); /* lower_tail = 0, give_log = 0 */
	
	/* asymptotic p-values of nested-subsets data-driven test */
	if ( (*datadriven == 1) || (*datadriven == 3)) {
		kk = (*d0>=1 ? *d0 : 1);
		*pval_nested_asympt = pchisq(*stat_nested,kk,0,0);
		if ((*d0<=1) && (*emptyset==0)) /* two-term approx. if nested subsets with min. dim. 1 and BIC */
			h_approx_nested(stat_nested,n,pval_nested_h);
/* 		h_approx_nested_k0(stat_nested,n,&kk,pval_nested_h); */
/*		h_approx_nested_k0_2(stat_nested,n,&kk,pval_nested_h); */
	}
	
	/* asymptotic p-values of all-subsets data-driven test */
	if ( (*datadriven == 2) || (*datadriven == 3)) {
		/* now compute the same by simulations (simulate the distribution of max of N(0,sigma)^2 */
		/* sigmawork still (from the last call of neyman) contains cholesky of sigma which is F*D*F^T, D on diag, F in lower triangle
		now change it into E*E^T, i.e., E = F*D^(1/2) (this will spare multiplications in simulation loop) */
		for (j=0; j<*d; j++) {
			sigmawork[j+j**d] = sqrt(sigmawork[j+j**d]);
			for (i=j+1; i<*d; i++)
				sigmawork[i+j**d] *= sigmawork[j+j**d];
		}
		if (*d0 == 0) { /* d0=0, hence the asympt. distr. is max of chisq1 which is simualted */
			/* p-value of max of indep. chisq1 (independence holds only for time transformation "sigma"); not used */
			/* 
			*pval_all_maxchisq1_indep = 1. - pow(pchisq(*stat_all,1.,1,0),*d);
			*/
			*pval_all_asympt = 0.;
			/* *pval_all_h = 0.; */
			kk = 2; /* dimension of 2-dim stats needed in quadrstat (dimension is passed as a pointer) */
			if (*nsim_asympt>0) {
				for (i=0; i<*nsim_asympt; i++) {
					for (j=0; j<*d; j++) {
						scorework[j] = norm_rand();
					}
					temp1 = 0.; /* temp1 will be the max of scorework (1-dimensional score statistics) */
					temp2 = 0.; /* temp2 will be the max of two-dimensional statistics */
					for (j=*d-1; j>=0; j--) {
						scorework[j] *= sigmawork[j+j**d];
						for (k=0; k<j; k++)
							scorework[j] += sigmawork[j+k**d]*scorework[k];
						/* j-th component of scorework is ready; scorework is E*scorework, i.e., it is N(0,sigma) */
						if (scorework[j]*scorework[j] > temp1)
							temp1 = scorework[j]*scorework[j]; /* finding max of scorework, i.e., max of one-dim stats */
						for (k=j+1; k<*d; k++) { /* loop over 2x2 submatrices */
							/* submatrix is stored in first 4 elements of a and contains j-th and k-th row and col of sigma */
							a[0] = a[3] = 1.; /* diag of sigma, hence of a, is 1 */
							a[1] = a[2] = sigma[j+k**d]; /* upper triangle would be enough in cholesky but for clarity... */
							quadrstat(scorework, a, &kk, &temp3, psi, choltol); /* scorework left unchanged, psi is the work array */
							if (temp3 > temp2)
								temp2 = temp3;
						} /**  tady bylo zapomenuto spocitat max (temp2) **/
					}
					*pval_all_asympt += (*stat_all < temp1);
					/* *pval_all_h += ( (temp1>*stat_all)&&((temp2-temp1)<logn) ) + ( (temp2>*stat_all)&&((temp2-temp1)>=logn) ); */
				}
				*pval_all_asympt /= *nsim_asympt;
				/* *pval_all_h /= *nsim_asympt; */
			}
		} else { /* d0>0, hence asymptotic distribution is chisq_d0 */
			*pval_all_asympt = pchisq(*stat_all,*d0,0,0);
		}
	}
	
	/* reconstruct diag of sigma (cholesky destroyed it) and copy the upper triangle to lower */
	for (j=0; j<*d; j++) {
		sigma[j+j**d] = 1.;
		for (k=0; k<j; k++)
			sigma[j+k**d] = sigma[k+j**d];
	}
	
	PutRNGstate();
	
}

void neyman_score(double *score, double *sigma, int *d, double *time, int *event, int *group,
	int *n, int *y1, int *y2, void (*p_basis)(), int *timetransf, double *a, double *psi)
{
	/* computes the score vector and sigma (may be applied to the observed data or to permuted data [only group differs]) */
	
	/* y1 and y2 are destroyed */
	
	int i,j,k,y1bak,y2bak;
	double temp1,temp2;
	
	/* compute the predictable (left-continuous) time transformation for the basis functions and store it in a */
	/* for permutation tests: with timetransf=1,2,3, a could be computed only once outside this function
	(a doesn't depend on the permutation)
	but for bootstrap: everything must be done here */
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
				if (*timetransf == 4) { /* sigma (left-cont.) */
					y1bak = *y1;
					y2bak = *y2;
					a[0] = 0.;
					for (i=1; i<*n; i++) { /* left-cont., thus group[i-1], event[i-1], and y1,y2 decremented afterwards */
						a[i] = a[i-1] + (double) event[i-1]**y1**y2/((*y1+*y2)*(*y1+*y2));
						if (group[i-1] == 1) {
							--*y1;
						} else {
							--*y2;
						}
					}
					for (i=0; i<*n; i++)
						a[i] /= a[*n-1];
					*y1 = y1bak;
					*y2 = y2bak;
				}
			}
		}
	}
	
	for (i=0; i<*d; i++) {
		score[i] = 0.;
		for (j=0; j<*d; j++) sigma[i+j**d] = 0.;
	}
	
	for (i=0; i<*n; i++) {
		if (event[i] == 1) {
			(*p_basis)(a+i, 1, psi, *d-1, 1, 1, 0);
			if (group[i] == 1) {
				temp1 = - (double) *y2/(*y1+*y2);
				temp2 = - (double) *y1/(*y1+*y2);
			} else {
				temp1 = + (double) *y1/(*y1+*y2);
				temp2 = + (double) *y2/(*y1+*y2);
			}
			temp2 *= temp1;
			for (j=0; j<*d; j++) {
				score[j] += psi[j]*temp1;
				for (k=j; k<*d; k++) {
					sigma[j+k**d] += psi[j]*psi[k]*temp2; /* only the upper triangle; cholesky doesn't need the lower one */
				}
			}
		}
		if (group[i] == 1) {
			--*y1;
		} else {
			--*y2;
		}
	}
	
	/* cholesky of sigma is stored in the lower triangle including diagonal, the upper triangle
	is unchanged;
	advantage of standardisation: diag of sigma is 1 (cor matrix), hence it can be reconstructed
	after cholesky, no need of work array for sigma */
	
	/* standardise everything, including increments of the score process; this facilitates computation of score stats (in cholesky) */
	for (j=0; j<*d; j++) {
		psi[j] = sqrt(sigma[j+j**d]); /* using psi as temporary storage for sqrt of diagonal of sigma */
		score[j] /= psi[j];
		sigma[j+j**d] = 1.;
		for (i=0; i<j; i++) {
			sigma[i+j**d] /= (psi[i]*psi[j]); /* make cor from cov */
		}
	}
}

void neyman_stat(double *score, double *scorework, double *sigma, double *sigmawork,
	int *d, int *d0, double *stat_d,
	double *stat_nested, double *stats_nested, double *stats_penal_nested, int *s_nested,
	double *stat_all, double *stats_all, double *stats_penal_all, int *s_all,
	int *subsets, int *nsubsets, double *penalty, double *work, double *choltol,
	int *datadriven, int *emptyset)
{
	/*
	emptyset	should the empty set be included (1, for AIC), or excluded (0, for BIC)
	penalty		penalty in the selection rule (BIC: log(n), AIC: 2)
	*/
	
	/* computes smooth test statistics (fixed-dim and data-driven) from the score u */
	/* u is either the observed score or the LWY-simulated score or the permutation score or... */
	int i,j,k,kk;
	double temp1;
	
	/* fixed dimension test */
	/* compute the test statistic score^T * sigma^{-1} * score */
	quadrstat(score, sigma, d, stat_d, work, choltol);
	
	/* reconstruct diag of sigma (cholesky in quadrstat destroyed it) (data-driven tests will need it) */
	for (j=0; j<*d; j++) {
		sigma[j+j**d] = 1.;
	}
	
	/* nested subsets, or both nested and all */
	if ((*datadriven == 1) || (*datadriven == 3)) {
		*s_nested = 0;
		if (*d0 > 0) {
			kk = *d0;
		} else {
			kk = 1;
		}
		temp1 = - *penalty**d;
		for (k=kk ; k<=*d; k++) { /* k=1,...,d, not 0,...,d-1, because of quadrstat which requires k>=1; in fact k=d0,...,d */
			for (i=0; i<k; i++) {
				for (j=i; j<k; j++) { /* sigma is only in the upper triangle but cholesky doesn't need lower on input */
					sigmawork[i+j*k] = sigma[i+j**d];
				}
			}
			quadrstat(score, sigmawork, &k, stats_nested+k-kk, work, choltol); /* use psi as the work array, sigmawork is changed */
			stats_penal_nested[k-kk] = stats_nested[k-kk] - *penalty*k;
			if (stats_penal_nested[k-kk] > temp1) {
				*s_nested = k-kk; /* that is (k-kk+1)-th of possible subsets was selected, i.e., phi_1,...,phi_k */
				temp1 = stats_penal_nested[k-kk];
			}
		}
		if ((*emptyset > 0) && (*d0 == 0)) { /* empty set included (AIC) */
			if (0 > temp1) {
				*s_nested = -1;
				*stat_nested = 0.;
			} else {
				*stat_nested = stats_nested[*s_nested];
			}
		} else { /* empty set excluded (BIC) */
			*stat_nested = stats_nested[*s_nested];
		}
	}
		
	/* all subsets, or both nested and all */
	if ( (*datadriven == 2) || (*datadriven == 3)) {
		*s_all = 0;
		temp1 = - *penalty**d;
		for (i=0; i<*nsubsets; i++) {
			kk = subsets[i+(*d)**nsubsets]; /* # of elements of the i-th subset */
			/* sigmawork and scorework will be the submatrix and subvector of sigma and score corresponding to subsets[i,] */
			for (j=0; j<kk; j++) {
				for (k=j; k<kk; k++) { /* the upper tringle is enough, lower not needed by cholesky on input */
					sigmawork[j+k*kk] = sigma[subsets[i+j**nsubsets]+subsets[i+k**nsubsets]**d];
				}
				scorework[j] = score[subsets[i+j**nsubsets]];
			}
			quadrstat(scorework, sigmawork, &kk, stats_all+i, work, choltol); /* use psi as work array, sigmawork is changed */
			stats_penal_all[i] = stats_all[i] - kk**penalty;
			if (stats_penal_all[i] > temp1) {
				*s_all = i;
				temp1 = stats_penal_all[i];
			}
		}
		if (*emptyset == 0) { /* empty set excluded (BIC) */
			*stat_all = stats_all[*s_all];
		} else { /* empty set included (AIC) */
			if (0 > temp1) {
				*s_all = -1;
				*stat_all = 0.;
			} else {
				*stat_all = stats_all[*s_all];
			}
		}
	}
}
