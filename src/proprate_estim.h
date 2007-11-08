
void twosample_proprate_estim_r_call(int *event, int *group, int *n, double *beta, int *maxiter,
	double *eps, int *flag, int *model, double *q1, double *q2,
	double *q1dot, double *q2dot, double *u1_process, double *u1, double *d11_process,
	double *d11, double *sigma11, double *logliks);

void twosample_proprate_estim(int *event, int *group, int *n, double *beta, int *maxiter,
	double *eps, int *flag, double (*p_q)(), double (*p_qdot)(), double *q1, double *q2,
	double *q1dot, double *q2dot, double *u1_process, double *u1, double *d11_process,
	double *d11, double *sigma11, double *logliks);

void twosample_proprate_score(int *event, int *group, int *n, double *beta, double *q1, double *q2,
	double *q1dot, double *q2dot, double *u1_process, double *u1, double *d11_process, double *d11,
	double *sigma11, double *loglik);

double q_prophaz(double t);

double qdot_prophaz(double t);

double q_propodds(double t);

double qdot_propodds(double t);
