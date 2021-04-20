from numpy import log, pi, interp, linalg, sqrt, log, exp, sin, cos
from scipy import integrate
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
import numpy as np
import scipy.integrate as integrate

cosmo = section_names.cosmological_parameters
cmb = section_names.cmb_cl
matter_powspec = section_names.matter_power_lin


def setup(options):
    config = {}
    config["output_section"] = options.get_string(option_section, "output_section", "growth_parameters")

    config["gamma_parametrization"] = options.get_bool(option_section, "gamma_parametrization", False)


    config["zmin"] = options.get_double(option_section, "zmin")
    config["zmax"] = options.get_double(option_section, "zmax")
    config["dz"] = options.get_double(option_section, "dz")
    config["zmax_log"] = options.get_double(option_section, "zmax_log")
    config["nz_log"] = options.get_int(option_section, "nz_log")

    sampling = np.arange(config["zmin"], config["zmax"]+config["dz"], config["dz"])
    sampling = np.append(sampling, np.logspace(np.log10(config["zmax"]+config["dz"]),np.log10(config["zmax_log"]), config["nz_log"] ))
    config["sampling"] = sampling
    print(sampling)
    return config


def execute(block, config):
    gamma_parametrization = config["gamma_parametrization"]

    z_output = config["sampling"]
    a_output = 1/(1+z_output)

    # Get parameters from sampler
    omega_m0 = block[cosmo, 'omega_m']
    omega_lambda0 = block[cosmo, 'omega_lambda']
    w0 = block[cosmo, 'w']
    wa = block[cosmo, 'wa']

    #internal sampling
    a = np.linspace(0.,1., 100)

    if(gamma_parametrization):
        gamma0 = block[cosmo, 'gamma0']
        gammaa = block[cosmo, 'gammaa'] #set to 0 if only want to sample over gamma0
        #https://arxiv.org/pdf/1606.00180.pdf  1.8.5 (page 108)
        exponent_factor = gamma0 + (1.- a)*gammaa
    else:
        w0 = block[cosmo, 'w']
        wa = block[cosmo, 'wa']
        #Linder2008: 0.55 + 0.05* (w(z=1)) 
        #and for w<-1: 0.55 + 0.02* (w(z=1)) 
        #z = 1 -> a = 1/2
        weff = w0 + (1.- 1./2.)*wa
        if (weff >= -1): 
            exponent_factor = 0.55 + 0.05 * (1+weff)
        else:
            exponent_factor = 0.55 + 0.05 * (1+weff)

    #integrating growth rate to get growth function
    #$$D(a)/a = exp(\int_0^a da/a [f(a) -1])$$

    a[0] = 0.0001 #to avoid dividing by zero, will handle case a = 0 by hand further down

    #friedman, for each component time scaling is a^(-3 (1+w(a)))
    Hsquare = omega_m0 * a**(-3) + omega_lambda0 * a**(-3 * (1+   (w0 + (1. -a)*wa)))
    omega_m = omega_m0 * a**(-3) /(Hsquare)

    f = omega_m**exponent_factor

    integrand = (f-1)/(a)

    #limit a=0
    f[0] = 1
    integrand[0] = 0
    a[0] = 0

    D = (integrate.cumtrapz(integrand, a, initial=0.))
    D = np.exp(D)
    D = D*a



    # Save back to block
    output_section = config["output_section"]
    block[output_section, 'a'] = a_output
    block[output_section, 'z'] = z_output
    block[output_section, 'f_z'] = np.interp(a_output,a,f)
    block[output_section, 'd_z'] = np.interp(a_output,a,D)

    # signal that everything went fine
    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
