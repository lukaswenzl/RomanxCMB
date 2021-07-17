int growth_de (int n, double *a, double *d, double *f);

double w (double a);
double w_int(double z, void *param);
//double DarkEnergy_a( double a );

double mg_mu_t(double a);
double mg_mu_t_propto_de(double a);
double mg_mu_t_f_of_R(double a);

double Xde_int (double a ,void * params);
double Xde (double a);

int func (double t, const double y[], double f[], void *params);
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params);

//TODO change arg if needed
int get_growthfactor(int n, double *a,double om, double ov, double w, double w2, double *d, double *f, int model, double mu0, int f_of_R_n, double f_of_R_fR, double k_in);
//int get_growthfactor(double a,double om, double w, double w2, double *gf,double mu);

double get_f_of_R_scale_dependence(double a);
double get_f_of_R_mass_of_a(double a);