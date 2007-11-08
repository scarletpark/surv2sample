
void twosample_proprate_neyman(double *time, int *event, int *group, int *n, int *model,
	double *beta, int *maxiter, double *eps, int *flag, int *d, int *d0, int *timetransf,
	int *basis, double *penalty, double *choltol, double *score, double *sigma,
	double *stat_d, double *pval_d, double *stat_d0, double *pval_d0,
	double *stat_nested, double *stats_nested, double *stats_nested_penal, int *s_nested,
	double *pval_nested_asympt, double *pval_nested_h);

void twosample_proprate_neyman_score(double *score, double *sigma, int *event, int *group, int *n,
	double *beta, double *q1, double *q2, double *q1dot, double *q2dot, double *psi, int *d,
	double *u1_process, double *d11_process, double *d11, double *sigma11, double *d21,
	double *sigma21);

