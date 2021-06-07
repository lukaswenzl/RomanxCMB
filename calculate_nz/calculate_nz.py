from __future__ import print_function
from builtins import range
import numpy as np
from cosmosis.datablock import option_section, names as section_names
from scipy import integrate


def setup(options):
    # only one parameter - filepath
    filename = options[option_section, "filepath"]

    input_section = options.get_string(option_section, "input_section", "")

    output_section = options.get_string(
        option_section, "output_section", default=section_names.wl_number_density)


    nbin = options.get_int(option_section, "nbin")
    zmin = options.get_double(option_section, "z_min", 0.0)
    zmax = options.get_double(option_section, "z_max", 4.0)
    nsamples = options.get_int(option_section, "nsamples", 100)
    sampling = np.linspace(0., zmax, nsamples)
    if(nsamples > 400):
        print("WARNING: the requested sampling accuracy is quite high. The function has only been tested up to about a sampling of 100 per redshift range of 1")

    full_distribution = Distribution(filename)
    use_only_one_sigma = options.get_bool(option_section, "use_only_one_sigma", False)

    # print(nsamples)
    # print(sampling)


    print("Loaded file %s" % (filename))
    return (full_distribution, nbin, sampling, input_section, zmin, output_section,use_only_one_sigma)


def execute(block, config):
    (full_distribution, nbin, sampling, input_section, zmin, output_section,use_only_one_sigma) = config
    bias = [block[input_section, "bias_%d" % i] for i in range(1, nbin+1)]
    if(use_only_one_sigma):
        overal_sigma_value = block[input_section, "sigma"]
        sigma = [overal_sigma_value for i in range(1, nbin+1)]
    else: 
        sigma = [block[input_section, "sigma_%d" % i] for i in range(1, nbin+1)]
    for i in range(nbin): 
        if(sigma[i] < 0.0008):
            sigma[i] = 0.0008 #smaller sigmas lead to numerical inaccuracy. If you want to lower this you need to test the integral convergence

    n_of_z, z = calc_bins(full_distribution.devide_into_equal_bins(n=nbin, z_min=zmin), full_distribution, sampling=sampling, sigma = sigma, Delta_i = bias)
    for col in n_of_z:
        norm = np.trapz(col, z)
        col /= norm

    block[output_section, 'nz'] = len(sampling)
    block[output_section, 'nbin'] = nbin
    block[output_section, 'z'] = sampling

    for (bin, bin_n_of_z) in enumerate(n_of_z):
        name = "bin_%d" % (bin + 1)
        block[output_section, name] = bin_n_of_z

    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0




##classes and functions
class Distribution():
    #redshift info
    low = np.array([])
    mean = np.array([])
    high = np.array([])
    #bin value
    value = np.array([])

    z_min = 0
    z_max = 4

    def __init__(self,filename):
        data =np.loadtxt(filename)
        self.low = data[:,0]
        self.mean = data[:,1]
        self.high = data[:,2]
        self.value = data[:,3]

    #interpolate result
    def interpolate(self, z):
        return np.interp(z, self.mean, self.value)
    def __call__(self, z):
        if(z< self.z_min or z > self.z_max): #seems to be never called since z in bin integral by definition >= zmin
            return 0
        return self.interpolate(z)


    def devide_into_equal_bins(self, n=10, z_min=0, z_max=4):
        #cumulative distribution (possibly not the most accurate...)
        self.z_min = z_min 
        self.z_max = z_max
        mask = np.logical_and(self.mean>=z_min,  self.mean<=z_max)
        cumsum = np.cumsum(self.value[mask])
        step = cumsum[-1]/n
        bin_cuts = np.arange(0, cumsum[-1]+step, step)
        return np.interp(bin_cuts, cumsum, self.high[mask])#high is working better than mean (makes sense)
        # low_cuts = np.interp(bin_cuts, cumsum, self.low)
        # high_cuts = np.interp(bin_cuts, cumsum, self.high)
        # return low_cuts, high_cuts


def calc_bins(equal_bins, distribution, sampling=np.linspace(0,4, 100), sigma = 0.01, Delta_i = 0.0): #need to make Delta_i an array!
    n_bins = len(equal_bins)-1
    if (not isinstance(Delta_i, list)):
        Delta_i = np.zeros(n_bins) + Delta_i
    if (not isinstance(sigma, list)):
        sigma = np.zeros(n_bins) + sigma
    bins = np.zeros((n_bins, len(sampling)))
    for i in range(n_bins):
        z_min = equal_bins[i]
        z_max = equal_bins[i+1]
        p = P(sampling[0], sigma[i], Delta_i[i])
        #for zph_index in range(len(sampling)):
        #    if(sampling[zph_index] > z_min - sigma*5 and sampling[zph_index] < z_max + sigma*5):
         #       p.set_z(sampling[zph_index])
          #      bins[i,zph_index ]  = integrate.quad(lambda z: distribution(z)*p(z), z_min, z_max, epsabs=1e-4, epsrel=1e-4)[0]

        loc_range = np.arange(z_min, z_max, 1./400)
        loc_dist = np.array([distribution(z) for z in loc_range])
        loc_mask = np.logical_and(sampling > z_min - sigma[i]*5 - np.abs(Delta_i[i]), sampling < z_max + sigma[i]*5 + np.abs(Delta_i[i]))
        for zph_index in range(len(sampling)):
            if(loc_mask[zph_index]):
                p.set_z(sampling[zph_index])
                #bins[i,zph_index ]  = integrate.quad(lambda z: distribution(z)*p(z), z_min, z_max, epsabs=1e-4, epsrel=1e-4)[0]
                bins[i,zph_index ]  = integrate.trapz(loc_dist*p(loc_range) , loc_range)

    return bins, sampling
class P():
    #Gaussian distribution weighted by redshift
    zph = 0.
    sigma = 0.
    Delta_i = 0.
    count = 0
    def __init__(self,zph, sigma, Delta_i):
        self.zph = zph 
        self.sigma = sigma 
        self.Delta_i = Delta_i
    
    def set_z(self,zph):
        self.zph = zph
    def __call__(self,z):
        self.count += 1
        prefactor = 1/(np.sqrt(2*np.pi)*self.sigma*(1+z))
        expfactor = -1/(2*(self.sigma* (1+z))**2)
        return prefactor * np.exp(expfactor * (z - self.zph - self.Delta_i)**2)
