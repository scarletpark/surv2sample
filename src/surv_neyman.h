
void twosample_neyman(double *time, int *event, int *group, int *n, int *d, int *d0,
	int *timetransf, int *basis, int *datadriven, double *penalty, int *emptyset,
	double *choltol, int *nsim_asympt, int *nperm, int *nboot, double *score, double *sigma,
	double *stat_d, double *pval_d, double *pval_d_perm, double *pval_d_boot,
	double *stat_nested, double *pval_nested_asympt, double *pval_nested_h,
	double *pval_nested_perm, double *pval_nested_boot,
	double *stats_nested, double *stats_penal_nested, int *s_nested,
	double *stat_all, double *pval_all_asympt, double *pval_all_perm, double *pval_all_boot,
	double *stats_all, double *stats_penal_all, int *s_all, int *all_subsets);

void neyman_score(double *score, double *sigma, int *d, double *time, int *event, int *group,
	int *n, int *y1, int *y2, void (*p_basis)(), int *timetransf, double *a, double *psi);

void neyman_stat(double *score, double *scorework, double *sigma, double *sigmawork,
	int *d, int *d0, double *stat_d,
	double *stat_nested, double *stats_nested, double *stats_penal_nested, int *s_nested,
	double *stat_all, double *stats_all, double *stats_penal_all, int *s_all,
	int *subsets, int *nsubsets, double *penalty, double *work, double *choltol,
	int *datadriven, int *emptyset);

