
void twosample_ks_cm_ad(double *time, int *event, int *group, int *n, double *rho, double *gamma,
	int *nsim, int *nperm, int *nboot, double *du, int *nsim_plot, int *nperm_plot, int *nboot_plot,
	double *du_sim_plot, double *du_perm_plot, double *du_boot_plot,
	double *stat_ks_w, double *pval_ks_w_sim, double *pval_ks_w_perm, double *pval_ks_w_boot, double *pval_ks_w_asympt,
	double *stat_ks_b, double *pval_ks_b_sim, double *pval_ks_b_perm, double *pval_ks_b_boot, double *pval_ks_b_asympt,
	double *stat_cm_w, double *pval_cm_w_sim, double *pval_cm_w_perm, double *pval_cm_w_boot,
	double *stat_cm_b, double *pval_cm_b_sim, double *pval_cm_b_perm, double *pval_cm_b_boot,
	double *stat_ad_w, double *pval_ad_w_sim, double *pval_ad_w_perm, double *pval_ad_w_boot,
	double *stat_ad_b, double *pval_ad_b_sim, double *pval_ad_b_perm, double *pval_ad_b_boot);

void ks_cm_ad_stat(int *n, double *du, double *g, double *ks_b, double *cm_w, double *cm_b,
	double *ad_w, double *ad_b, double *stat_ks_w, double *stat_ks_b,
	double *stat_cm_w, double *stat_cm_b, double *stat_ad_w, double *stat_ad_b);

void ks_cm_ad_process(int *event, int *group, int *n, int *y1, int *y2, double *rho, double *gamma,
	int *first_event, double *du, double *dsigma, double *ks_b, double *cm_w, double *cm_b,
	double *ad_w, double *ad_b);

