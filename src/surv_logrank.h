
void twosample_grhogamma(double *time, int *event, int *group, int *n, double *rho, double *gamma,
	double *b, int *m, int *nsim_asympt, int *nperm, int *nboot, double *choltol,
	double *stats, double *sigma, double *pvals, double *pvals_perm, double *pvals_boot,
	double *stat_sum, double *pval_sum, double *pval_sum_perm, double *pval_sum_boot,
	double *stat_max, double *pval_max, double *pval_max_perm, double *pval_max_boot);

void grhogamma_stat(int *event, int *group, int *n, double *rho, double *gamma,
	double *b, int *m, int *y1, int *y2, double *incr, double *stats, double *stat_sum,
	double *stat_max, double *sigma);

