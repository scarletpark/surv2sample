
void twosample_incidence_ks(int *event, int *group, int *n, int *nsim,
	double *f11, double *f12, double *f21, double *f22,	double *test_process,
	double *stat, double *pval_sim, double *test_process_plot_sim, int *nsim_plot);

void twosample_incidence_ks_process(int *event, int *group, int *n, int *y1, int *y2,
	double *f11, double *f12, double *f21, double *f22, double *s1, double *s2,
	double *test_process);

void twosample_incidence_lwy_process(int *event, int *group, int *n, int *y1, int *y2,
	double *f11, double *f12, double *f21, double *f22, double *g, double *test_process_sim);
