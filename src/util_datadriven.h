
void datadriven_stats_nested(double *score_d, double *sigma_d, int *d, int *d0,
	double *stats, double *stats_penal, int *selected,
	double *stat_d, double *stat_d0, double *stat_selected,
	double *sigmawork, double *scorework, double *penalty, double *choltol);

void h_approx_nested(double *stat_bic, int *n, double *pval_bic_h);
double h1_approx_nested(double *x, int *n);
double h2_approx_nested(double *x, int *n);

void h_approx_nested_f(double *stat_bic, int *n, double *pval_bic_h_f);
double h1_approx_nested_f(double *x, int *n);
double h2_approx_nested_f(double *x, int *n);

