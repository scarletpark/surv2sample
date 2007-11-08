
void twosample_incidence_neyman(double *time, int *event, int *group, int *n,
	int *d, int *d0, int *basis, double *penalty, int *logrank, double *logrank_r,
	double *choltol, double *score,	double *sigma, double *f11, double *f12, double *f21,
	double *f22, double *stat_d, double *pval_d, double *stat_d0, double *pval_d0,
	double *stat_nested, double *stats_nested, double *stats_nested_penal,
	int *s_nested, double *pval_nested_asympt, double *pval_nested_h);

void twosample_incidence_psi(double *time, int *event, int *group, int *n,
	int *y1, int *y2, double *a1, double *a2, double *condf01,
	double *psi, int *d, void (*p_basis)(), int *logrank, double *logrank_r);

void twosample_incidence_score(double *score, double *sigma, int *d, int pooled,
	double *time, int *event, int *group, int *n, int *y1, int *y2,
	double *s1, double *s2, double *f11, double *f12, double *f21, double *f22, double *f01,
	double *psi, double *h1, double *h2, double *q1, double *q2, double *rho, double *work);

void twosample_incidence_score_variance(double *q, double *h, double *rho, int *n, int *d,
	double *sigma, double *intqrho);
