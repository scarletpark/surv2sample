
void quadrstat(double *score, double *sigma, int *d, double *stat, double *work, double *choltol);

void cosine_basis(double *x, int n, double *p, int deg, int degzero,
	int dummy_shifted, int dummy_normalized);
void cosine1_basis(double *x, int n, double *p, int deg, int degzero,
	int dummy_shifted, int dummy_normalized);
void legendre_basis(double *x, int n, double *p, int deg, int degzero,
	int shifted, int normalized);
