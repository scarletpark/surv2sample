
void twosample_incidence_pepe(double *time, int *event, int *group, int *n, double *tau,
	int *nsim, double *f11, double *f12, double *f21, double *f22, double *stat,
	double *pval_asympt, double *pval_sim);

double twosample_incidence_pepe_variance(double *time, int *event, int *group, int *n,
	double *tau, int grp, double *fj1, double *fj2, int *yj, double *var);

double int_stat_pepe(double *f, double *time, int *n, double *tau);
