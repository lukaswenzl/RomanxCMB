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
    sigma8_input = block[cosmo, 'sigma8_input']
    if block.has_value(cosmo, 'sigma_8'):
        sigma8_before = block[cosmo, 'sigma_8']
        boltzmann = True
    else:
        boltzmann = False

    A_s = block[cosmo, 'A_s']
    if(boltzmann):
        TT = block[cmb, 'TT']
        EE = block[cmb, 'EE']
        BB = block[cmb, 'BB']
        TE = block[cmb, 'TE']
    P_k = block[matter_powspec, 'P_k']

    zmin = block[matter_powspec, 'z'].min()
    if zmin != 0.0:
        raise ValueError(
            "You need to set zmin=0 in CAMB to use the sigma8_rescale module.")
    
    if(not boltzmann):
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
        sigma8_before = sqrt(result * 1./2./pi/pi)          
        #print(sigma8_before)


    # Calculate rescale factor
    r = (sigma8_input**2) / (sigma8_before**2)

    # Rescale CMB Cl and matter power outputs
    A_s *= r
    if(boltzmann):
        TT *= r
        EE *= r
        BB *= r
        TE *= r
    P_k *= r

    # Save back to block
    block[cosmo, 'A_s'] = A_s
    if(boltzmann):
        block[cmb, 'TT'] = TT
        block[cmb, 'EE'] = EE
        block[cmb, 'BB'] = BB
        block[cmb, 'TE'] = TE
    block[matter_powspec, 'P_k'] = P_k

    block[cosmo, 'sigma_8'] = sigma8_input

    # signal that everything went fine
    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
