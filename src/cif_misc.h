
void twosample_incidence_f01(int *event, int *y1, int *y2, double *s1, double *s2,
	double *f01, int *n);

double incr2sym(double *a, int n, int i, int j);

double incr1sym(double *a, int n, int j);

void fj1_covariance(double *rho, int *n, int grp, int *event, int *group,
	double *fj1, double *fj2, int *yj);
