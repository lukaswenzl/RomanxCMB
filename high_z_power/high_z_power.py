"""
Use the growth factor to extrapolate the P(k) we have to higher redshifts.

"""
import numpy as np
from cosmosis.datablock import option_section, names
import scipy.interpolate
from cosmosis.runtime.utils import Timer

def setup(options):
    sections = options.get_string(option_section, "sections", names.matter_power_nl).split()
    zmax = options.get_double(option_section, "zmax", 1150.0)
    nz = options.get_int(option_section, "nz", 50)
    verbose = options.get_bool(option_section, "verbose", False)

    return {"sections":sections, "zmax":zmax, "nz":nz, "verbose":verbose}


def extrapolate(block, section, D, zmax, nz):
    #Load the current matter power spectrum in this section
    z_current, k, P = block.get_grid(section, 'z', 'k_h', 'p_k')

    #Decide the new points at which to sample in redshift, between the existing
    #maximum and the new one
    zmax_current = z_current.max()
    z_sample = np.logspace(np.log10(zmax_current*1.1), np.log10(zmax), nz)

    #Make space for the 
    nz_current = len(z_current)
    nz_total = nz_current + nz
    nk = len(k)
    P_new = np.zeros((nz_total, nk))

    # Fill in the new output values with current data
    z_new = np.concatenate((z_current, z_sample))
    P_new[:nz_current,:] = P

    #Loop through the new points using the ratio of the squared growth function
    #to get the remaining values
    P_zmax_current = P[-1, :]
    d_zmax_current = D(zmax_current)
    for i,zi in enumerate(z_sample):
        P_new[nz_current+i, :] = P_zmax_current * (D(zi)/d_zmax_current)**2

    #Update the grid
    block.replace_grid(section, 'z', z_new, 'k_h', k, 'p_k', P_new)


class LogZInterpolator(object):
    "An interpolator for log(D) vs log(1+z)"
    def __init__(self, z, d):
        self.log1z = np.log(1+z)
        self.logd = np.log(d)
        self.interp = scipy.interpolate.interp1d(self.log1z, self.logd, kind='linear')
    def __call__(self, z):
        log1z = np.log(1+z)
        try:
            logd = self.interp(log1z)
        except:
            print ("z = "+str( z))
            print ("Log(1+z) = ", str(log1z))
            print ("Upper limit Log(1+z) = ", str(self.log1z[:-10]))
            raise
        return np.exp(logd)



def execute(block, config):
    #Extract options
    sections = config['sections']
    zmax = config['zmax']
    nz = config['nz']
    verbose = config['verbose']

    #Get the growth function interpolator
    z = block[names.growth_parameters, 'z']
    d_z = block[names.growth_parameters, 'd_z']
    D = LogZInterpolator(z, d_z)

    #For all desired sections run the extrapolation
    for section in sections:
        if verbose:
            print ("Using growth function to extrapolate {} to z={} in {} steps".format(section, zmax, nz))
        extrapolate(block, section, D, zmax, nz)

    return 0