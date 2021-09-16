#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cosmosis/datablock/c_datablock.h"
#include "cosmosis/datablock/section_names.h"
#include "growthfactor.h"

const char * cosmo = COSMOLOGICAL_PARAMETERS_SECTION;
const char * like = LIKELIHOODS_SECTION;
const char * growthparameters = GROWTH_PARAMETERS_SECTION;
const char * mg = "modified_gravity_parameters";

typedef struct growth_config {
	double zmin;
	double zmax;
	double dz;
	int nz_lin;
	double zmax_log;
	int nz_log;
} growth_config;

void reverse(double * x, int n)
{
	for (int i=0; i<n/2; i++){
		double tmp = x[i];
		x[i] = x[n-1-i];
		x[n-1-i] = tmp;
	}
}

growth_config * setup(c_datablock * options)
{
	int status = 0;
	growth_config * config = malloc(sizeof(growth_config));
	status |= c_datablock_get_double_default(options, OPTION_SECTION, "zmin", 0.0, &(config->zmin));
	status |= c_datablock_get_double_default(options, OPTION_SECTION, "zmax", 3.0, &(config->zmax));
	status |= c_datablock_get_double_default(options, OPTION_SECTION, "dz", 0.01, &(config->dz));
	status |= c_datablock_get_double_default(options, OPTION_SECTION, "zmax_log", 1100.0, &(config->zmax_log));
	status |= c_datablock_get_int_default(options, OPTION_SECTION, "nz_log", 0, &(config->nz_log));

	config->nz_lin = (int)((config->zmax-config->zmin)/config->dz)+1;

	printf("Will calculate f(z) and d(z) in %d bins (%lf:%lf:%lf)\n", config->nz_lin, config->zmin, config->zmax, config->dz);
	if (config->nz_log > 0){
	  printf("Will extend growth calculation using %d logarithmically spaced z-bins to zmax = %lf \n", config->nz_log, config->zmax_log);
	}

	if (config->zmax_log <= config->zmax){
		fprintf(stderr, "zmax_log parameter must be more than zmax in growth function module\n");
		exit(1);
	}
	
	// status |= c_datablock_get_double(options, OPTION_SECTION, "redshift", config);
	if (status){
		fprintf(stderr, "Please specify the redshift in the growth function module.\n");
		exit(status);
	}
	return config;

} 
    
int execute(c_datablock * block, growth_config * config)
{

	int i,status=0;
	double w,wa,omega_m,omega_v;

	int nz_lin = config->nz_lin;
	int nz_log = config->nz_log;
	int nz = nz_lin + nz_log;
	double zmin_log = config->zmax + config->dz;

	int mg_model; 
	double mu0;
	int f_of_R_n; 
	double f_of_R_fR;

	//TODO this is hard coded for now. Might want to make a parameter
	//If changed for other modules, need to change here by hand
	double kmin=1e-5; 
	double kmax=10.0;
	double k_large_scale = kmin;
	
	//read cosmological params from datablock
	status |= c_datablock_get_double_default(block, cosmo, "w", -1.0, &w);
	status |= c_datablock_get_double_default(block, cosmo, "wa", 0.0, &wa);
	status |= c_datablock_get_double(block, cosmo, "omega_m", &omega_m);
	status |= c_datablock_get_double_default(block, cosmo, "omega_lambda", 1-omega_m, &omega_v);

	status |= c_datablock_get_int(block, mg, "model", &mg_model);
	status |= c_datablock_get_double(block, mg, "mu0", &mu0);
	status |= c_datablock_get_int(block, mg, "f_of_R_n", &f_of_R_n);
	status |= c_datablock_get_double(block, mg, "f_of_R_fR", &f_of_R_fR);

	if (status){
		fprintf(stderr, "Could not get required parameters for growth function (%d)\n", status);
		return status;
	}

	//allocate memory for single D, f and arrays as function of z
	double *a = malloc(nz*sizeof(double));
	double *dz = malloc(nz*sizeof(double));
	double *fz = malloc(nz*sizeof(double));
	double *z = malloc(nz*sizeof(double));


	// First specify the z bins in increasing z, for my own sanity
	for (i=0;i<nz_lin;i++){
		z[i] = config->zmin + i*config->dz;
	}	

	for (i=0;i<nz_log;i++){
		z[nz_lin + i] = zmin_log*exp(i*(log(config->zmax_log) - log(zmin_log))/(nz_log-1));
	}


	for (i=0;i<nz;i++){
		a[i] = 1.0/(1+z[i]);
	}

	// Now reverse them to increasing a, decreasing z.
	// Note that we don't need to reverse z as we don't use it
	// in the function below
	reverse(a,nz);

	// Compute D and f
	status = get_growthfactor(nz, a, omega_m, omega_v, w, wa, dz, fz, mg_model, mu0, f_of_R_n, f_of_R_fR, k_large_scale);
	
	// Now reverse everything back to increasing z
	// Note that we do not unreverse z as we never reversed it in the first place.
	reverse(a,nz);
	reverse(dz,nz);
	reverse(fz,nz);

	status |= c_datablock_put_double_array_1d(block,growthparameters, "d_z", dz, nz);
	status |= c_datablock_put_double_array_1d(block,growthparameters, "f_z", fz, nz);
	status |= c_datablock_put_double_array_1d(block,growthparameters, "z", z, nz);
	status |= c_datablock_put_double_array_1d(block,growthparameters, "a", a, nz);

	reverse(a,nz);

	double k_value;
	double dk;
	int nk_steps=200;

	double *k = malloc(nk_steps*sizeof(double));

	int c = nz, r = nk_steps, j;
	
	double **d_k_z;
	double **f_k_z;
	double *d_ptr; 
	double *f_ptr; 

	int len;
	len = sizeof(double *) * r + sizeof(double) * c * r;
	d_k_z = (double**) malloc(len);
	f_k_z = (double**) malloc(len);

	d_ptr = (double *)(d_k_z + r);
	f_ptr = (double *)(f_k_z + r);

	// for loop to point rows pointer to appropriate location in 2D array
    for(i = 0; i < r; i++){
        d_k_z[i] = (d_ptr + c * i);
		f_k_z[i] = (f_ptr + c * i);
	}

	dk = (log(kmax) - log(kmin))/(nk_steps-1);
	for (i = 0; i < nk_steps; i++)
    {
		k[i] = exp(log(kmin) + i*dk);
	}

	for (i = 0; i < nk_steps; i++)
    {
		k_value = k[i];

		status = get_growthfactor(nz, a, omega_m, omega_v, w, wa, dz, fz, mg_model, mu0, f_of_R_n, f_of_R_fR, k_value);
		reverse(dz,nz);
		reverse(fz,nz);

		for (j = 0; j < c; j++){
			memcpy(&d_k_z[i][j], &dz[j], sizeof(dz[j]));
			memcpy(&f_k_z[i][j], &fz[j], sizeof(fz[j]));
		}

	}
		
	// Now reverse everything back to increasing z
	// Note that we do not unreverse z as we never reversed it in the first place.
	reverse(a,nz);

	status |= c_datablock_put_double_grid(block,growthparameters, "k_for_d_k_z", nk_steps, k, "z_for_d_k_z", nz, z, "d_k_z", d_k_z);
	status |= c_datablock_put_double_grid(block,growthparameters, "k_for_f_k_z", nk_steps, k, "z_for_f_k_z", nz, z, "f_k_z", f_k_z);
	
	free(d_k_z);
	free(f_k_z);
	free(k);

	free(fz);
	free(dz);
	free(z);
	free(a);
	
return status;
}


int cleanup(growth_config * config)
{
	free(config);
	return 0;
}
