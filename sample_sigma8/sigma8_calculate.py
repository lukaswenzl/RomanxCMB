from numpy import log, pi, interp, linalg, sqrt, log, exp, sin, cos
from scipy import integrate
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section

cosmo = section_names.cosmological_parameters
cmb = section_names.cmb_cl
matter_powspec = section_names.matter_power_lin


def setup(options):
    return 0


def execute(block, config):

    # Get parameters from sampler and CAMB output

    A_s = block[cosmo, 'A_s']
    P_k = block[matter_powspec, 'P_k']

    zmin = block[matter_powspec, 'z'].min()
    if zmin != 0.0:
        raise ValueError(
            "You need to set zmin=0 in CAMB to use the sigma8_rescale module.")
    
    #need to calculat the sigma8 value ourselfs
    k_h = block[matter_powspec, 'k_h']
    kmin = k_h.min()
    kmax = k_h.max()
    def W(k, R=8.):
        #fourier transformation of theta(||k|| < R)
        x= k *R
        return 3.*(sin(x)-x*cos(x))/x**3
    def integrand (logk):
        k = exp(logk)
        return k*k*k* W(k)**2 *interp(k, k_h, P_k[0, :]) 
        #note: extra factor of k because dlogk = k dk
        #we also assume index 0 is z=0
    result= integrate.quad(integrand, log(kmin), log(kmax), epsrel=3e-4, limit=100)[0]
    sigma8_calculated = sqrt(result * 1./2./pi/pi)          
        
    # print("hello check")
    # print(sigma8_calculated)
    # print(block[cosmo, 'sigma8_input'])

    block[cosmo, 'sigma_8_aftermodgrav'] = sigma8_calculated

    # signal that everything went fine
    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
