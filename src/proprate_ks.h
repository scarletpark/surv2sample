
void twosample_proprate_ks(int *event, int *group, int *n, int *model, double *beta,
	int *maxiter, double *eps, int *flag, int *nsim, double *u1_process,
	double *stat_ks, double *pval_ks_sim, double *test_process_plot_sim, int *nsim_plot);

void twosample_proprate_ks_sim(int *event, int *group, int *n, int *y1, int *y2, double *theta,
	double *q1, double *q2, double *q1dot, double *q2dot, double *d11_process, double *d11,
	double *u1_process_sim);
