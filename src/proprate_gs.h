
void twosample_proprate_gs(int *event, int *group, int *n, int *model,
	int *weight1, int *weight2, double *beta1, double *beta2, double *stat, double *pval);

void twosample_proprate_gs_aux(int *event, int *group, double *k1, double *k2,
	double *da1, double *da2, double *q1, double *q2, double *q1dot, double *q2dot,
	double *b11, double *b12, double *b21, double *b22, int *y1, int *y2, int *n,
	double (*p_q)(), double (*p_qdot)(), void (*p_weight1)(), void (*p_weight2)());

void twosample_proprate_gs_stat(double *k1, double *k2, double *da1, double *da2, double *q1, double *q2,
	int *n, double *rho11, double *rho12, double *rho21, double *rho22, double *stat);

void twosample_proprate_gs_var(double *k1, double *k2, double *da1, double *da2,
	double *q1, double *q2, double *q1dot, double *q2dot,
	double *b11, double *b12, double *b21, double *b22, int *y1, int *y2, int *n,
	double *rho11, double *rho12, double *rho21, double *rho22, double *var);

