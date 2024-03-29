//Adapted from existing code in CosmoSIS by Agnes Ferte to include mu
//Further adapted (and optimized) by Chen Heinrich (2021)
//

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include "growthfactor.h"

//Code to calculate the scale-indep. and dependent 
//  - linear growth factor D, 
//  - linear growth rate, f = dlnD/dlna for modified gravity models.

//Currently supporting:
//   Model 1: (Sigma,Mu) parametrisation
//     with a time dependance mu = 1 + mu(t) where mu(t) = mu0 * omega_DE(a)
//   Model 2: f(R) model parametrized by n and f_R

double D;
double linear_f;
double w0,wa;
double omega_m,omega_lambda;

int mg_model; 
double mg_mu0;
int mg_f_of_R_n;
double mg_f_of_R_fR; 

double k;  

#define LIMIT_SIZE 1000

int get_growthfactor(int n, double *a,double om, double ov, double w, double w2, double *d, double *f, int model, double mu0, int f_of_R_n, double f_of_R_fR, double k_in)
{
	w0 = w;
	wa = w2;
	omega_m = om;
	omega_lambda = 1.0 - omega_m; 

	mg_model = model; 
	mg_mu0 = mu0;
	mg_f_of_R_n = f_of_R_n;
	mg_f_of_R_fR = f_of_R_fR;
	k = k_in;

    int tmp;
    tmp = growth_de(n, a, d, f);

	return 0;
}

//careful: this is the mu(t) in 1+mu(t) (mu(t) = 0 = GR)
double mg_mu_t(double a)
{
	if (mg_model == 0){
		return mg_mu0;
	}
	else if (mg_model == 1){
		return mg_mu_t_propto_de(a);
	}
	else if (mg_model == 2){
		return mg_mu_t_f_of_R(a);
	}
}

double mg_mu_t_propto_de(double a)
{
	return mg_mu0/(omega_lambda + (1.0 - omega_lambda)*pow(a,-3.0));
}

double mg_mu_t_f_of_R(double a)
{
    double scale_dependence;
    scale_dependence = get_f_of_R_scale_dependence(a); 
    return (1.0 + 4.0/3.0*scale_dependence)/(1.0 + scale_dependence) - 1.0;
}

double w (double a)
{
return w0 + (1.0-a)*wa;
}

double w_int(double z, void *param)
{
	return (1. + w(1./(z+1)) )/( 1. + z);
}

double Xde_int (double a,void *params )
{
    if (a == 0.) a = 1e-3;
    double Xde_int = w(a)/a;
    return Xde_int;
}

double Xde (double a)
{
    gsl_integration_workspace * worksp  = gsl_integration_workspace_alloc (LIMIT_SIZE);
    gsl_function F;
    F.function = &Xde_int;
    F.params =0;
    double integral,error;
    if (a == 0.) 
		a = 1e-3;
    gsl_integration_qags (&F, a, 1,  1.0e-20, 1.0e-10, LIMIT_SIZE,worksp, &integral, &error); 
    //adaptive integration to get an estimate of the integral of Xde_int over a to 1. 
    gsl_integration_workspace_free (worksp);

    return omega_m/(1.-omega_m)*exp(-3.*integral);
}

int func (double t, const double y[], double f[], void *params)
{
	f[0] = y[1];
	f[1] = -( 3.5 - 1.5 * w(t)/( 1 + Xde(t) ) )*y[1]/t - 1.5 *( 1 - w(t) - Xde(t)*mg_mu_t(t) )/( 1 + Xde(t))*(y[0]/t/t);
	return GSL_SUCCESS;
}

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
	//double mu = *(double *)params; 
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2); 
	//matrix_view = temporary object, can be used to work on subset of matrix elements ----> here, dfdy_mat
	gsl_matrix * m = &dfdy_mat.matrix; 
	gsl_matrix_set (m, 0, 0, 0.0); //dy1/dx1
	gsl_matrix_set (m, 0, 1, 1.0);
	gsl_matrix_set (m, 1, 0, - 1.5 *( 1 - w(t) - Xde(t)*mg_mu_t(t) )/( 1 + Xde(t))*(1./t/t));
	gsl_matrix_set (m, 1, 1, -( 3.5 - 1.5 * w(t)/( 1 + Xde(t) ) )/t);
	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	return GSL_SUCCESS;
}
//gsl_matrix_set(matrix *m, size i, size j, double x) : set the value of (i,j) of matrix m to x. 
//---> here, m is df/dy, a 2 by 2 matrix.      


int growth_de(int n, double *a, double *d, double *f)
{
	const gsl_odeiv_step_type * T 
		= gsl_odeiv_step_rk4;
     
	gsl_odeiv_step * s 
         = gsl_odeiv_step_alloc (T, 2); //stepping function which advance solution from t to t+h
	gsl_odeiv_control * c 
         = gsl_odeiv_control_y_new (1e-6, 0.0); //control the error with respect to the solution
	gsl_odeiv_evolve * e 
         = gsl_odeiv_evolve_alloc (2); //evolve the solution (system of 2 dimension)
     
	double mu = 10;
	gsl_odeiv_system sys = {func, jac, 2, &mu}; //solve dy_i/dt = f_i(y_1, ..., y_n). 
	//func = store the elements of f(t,y,params) in the array dy/dt
	//jac = store the elements of f(t,y,params) in the array df/dt
	//2 is the dimension of the system
	//&mu is a pointer to the parameters of the system.
     
	double t = 1.e-3;
    double h = 1e-6;
	double y[2] = {1.,0.}; 
     
    int i;

    for (i = 0; i < n; i++)
    {
        double ti = a[i];
        while (t < ti)
        {
            gsl_odeiv_evolve_apply (e, c, s, 
                                    &sys, 
                                    &t, ti, &h,
                                    y);
        }
        d[i] = y[0]*a[i]; //D growth factor (= G*a)
	    f[i] = y[1]*a[i]*a[i]/(y[0]*a[i]) +1.; // f = d ln D/ d ln a
    }

	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);

    return 0;

}

//f(R) model

//TODO suboptimal implementation for now where for different 
// k values, the m(a) is calculated all over again for the same a.
double get_f_of_R_scale_dependence(double a)
{
    double mass; 
    double scale_dependence; 

    mass = get_f_of_R_mass_of_a(a);
    scale_dependence = pow(k/(a*mass), 2);

    return scale_dependence;
}

double get_f_of_R_mass_of_a(double a)
{
    double term1, term2, term3, term4; 
    double c_in_km_per_sec = 2.99792458e5;
    double H0_over_c_in_hoverMpc = 100./c_in_km_per_sec;

    term1 = pow(omega_m + 4.*omega_lambda, -0.5*(mg_f_of_R_n+1));
    term2 = pow((mg_f_of_R_n+1)*mg_f_of_R_fR, -0.5);
    term3 = pow(omega_m * pow(a, -3) + 4.*omega_lambda, 0.5 * (mg_f_of_R_n+2));
    term4 = H0_over_c_in_hoverMpc;

    return term1 * term2 * term3 * term4; 
}