//Adapted from existing code in CosmoSIS by Agnes Ferte to include mu
//Further adapted (and optimized) by Chen Heinrich (2021)
//TODO need to check if this works for varying w(a)

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include "growthfactor.h"

//Code to calculate the linear growth factor D, 
//and linear growth rate, f = dlnD/dlna for modified gravity.
//Currently supporting (Sigma,Mu) parametrisation
//with a time dependance mu = 1 + mu(t) where
// mu(t) = mu0 * omega_DE(a)

double D;
double linear_f;
double w0,wa;
double omega_m,omega_lambda;
double mg_mu;

#define LIMIT_SIZE 1000

int get_growthfactor(int n, double *a,double om, double ov, double w, double w2, double *d, double *f, double mu0)
{
	w0 = w;
	wa = w2;
	omega_m = om;
	omega_lambda = 1.0 - omega_m; 
	mg_mu = mu0;

    int tmp;
    tmp = growth_de(n, a, d, f);

	return 0;
}

//careful: this is the mu(t) in 1+mu(t) (mu(t) = 0 = GR)
double mg_mu_t(double a)
{
return mg_mu/(omega_lambda + (1.0 - omega_lambda)*pow(a,-3.0));
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
