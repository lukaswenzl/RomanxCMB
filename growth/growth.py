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

    z = config["sampling"]
    a = 1/(1+z)

    # Get parameters from sampler
    omega_m0 = block[cosmo, 'omega_m']
    omega_lambda0 = block[cosmo, 'omega_lambda']
    
    if(gamma_parametrization):
        gamma0 = block[cosmo, 'gamma0']
        gammaa = block[cosmo, 'gammaa'] #set to 0 if only want to sample over gamma0
        #https://arxiv.org/pdf/1606.00180.pdf  1.8.5 (page 108)
        exponent_factor = gamma0 + (1.- a)*gammaa
    else:
        w0 = block[cosmo, 'w']
        wa = block[cosmo, 'wa']
        #Linder2008: 0.55 + 0.05* (w(z=1))
        #z = 1 -> a = 1/2
        exponent_factor = 0.55 + 0.05 * (w0 + (1.- 1./2.)*wa)
        print("consider case where this gets positive? or am i avoiding this here already")

    

    Hsquare = omega_m0 * a**(-3) + omega_lambda0
    omega_m = omega_m0 * a**(-3) /(Hsquare)

    f = omega_m**exponent_factor

    loga = np.log(a)
    print("i might need finer sampling here and then downsample!")
    print("THis code still needs to be tested, I am not sure if it is accurate (maybe z is reversed!)")

    #D = np.exp(integrate.cumtrapz(f[::-1], loga[::-1], initial=0.))
    #D = D[::-1]
    D = np.exp(integrate.cumtrapz(f, loga, initial=0.))
    D = 5.*(omega_m0)/2. * D #normalization

    # Save back to block
    output_section = config["output_section"]
    block[output_section, 'a'] = a 
    block[output_section, 'z'] = z
    block[output_section, 'f_z'] = f
    block[output_section, 'd_z'] = D

    print("might need to do more sampling?!, but code should work now -> test")
    # signal that everything went fine
    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
