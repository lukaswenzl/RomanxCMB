from cosmosis.datablock import option_section, names
import numpy as np
import pdb
from scipy.interpolate import CubicSpline, interp1d
from scipy.ndimage import gaussian_filter1d
from cosmosis.datablock import option_section, names



def get_cov(x, scale, amplitude, regularizer=1e-10):
    n = len(x)
    cov = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            cov[i,j] = amplitude * np.exp(-0.5*(x[i]-x[j])**2 / scale**2)

    return cov + regularizer*np.identity(n)

def get_smoothness_prior(block, fz_s, zmax, scale, amplitude):
    zz_test = np.linspace(0,zmax,200)
    fz = fz_s(zz_test)
    smooth = gaussian_filter1d(fz, sigma=scale/(zz_test[1]-zz_test[0]))
    cov = get_cov(zz_test, scale, amplitude, regularizer=1e-10)
    icov = np.linalg.inv(cov)
    
    diff = fz - smooth
    chi2 = np.dot(diff, np.dot(icov, diff))
    lnlike = -0.5*chi2
    
    block[names.likelihoods, "rescale_Pk_fz_smoothness_like"] = lnlike
    block[names.data_vector, "rescale_Pk_fz_smoothness_theory"] = diff # all happens as if data=0 and theory=diff
    block[names.data_vector, "rescale_Pk_fz_smoothness_inverse_covariance"] = icov
    block[names.data_vector, "rescale_Pk_fz_smoothness_chi2"] = chi2

    return lnlike

def setup(options):
    mode = options.get(option_section, "mode")
    assert mode == 1

    smoothness_prior_scale = options.get_double(option_section, "smoothness_prior_scale", -1.0)
    smoothness_prior_amplitude = options.get_double(option_section, "smoothness_prior_amplitude", -1.0)

    config = {'smoothness_prior_scale':smoothness_prior_scale, 'smoothness_prior_amplitude':smoothness_prior_amplitude, 'mode':mode}

    if mode==0 or mode=='pk_grid':
        print("[rescale_Pk_fz] : mode 0 ")
    
    if mode==1 or mode=='interp1d':
        zmax = options.get(option_section, "zmax")
        nknots = options.get(option_section, "nknots")
        
        print("[rescale_Pk_fz] : mode 1 ")
        print(" - Using {} knots from 0 to {}".format( nknots, zmax))

        config['zmax'] = zmax     
        config['nknots'] = nknots     

    if mode==2 or mode=='pca':
        filepath = options.get(option_section, "filepath")
        neigen = options.get(option_section, "neigen")
        zmax = options.get(option_section, "zmax")

        v = np.loadtxt(filepath)
        zz = np.linspace(0., zmax, v.shape[0])
        ez = []
        for i in range(neigen):
            ez.append(CubicSpline(zz,v[:,i]/np.sqrt(np.trapz(v[:,i]**2,zz))))

        print("[rescale_Pk_fz] : mode 2 ")
        print(" - Using splines from file", filepath)
        print(" - Using number of eigen modes", neigen)
        
        config['neigen'] = neigen     
        config['ez'] = ez  
    
    return config

def execute(block, config):
    mode = config['mode']
    p_k = block["matter_power_lin", "p_k"]
    z = block["matter_power_lin", "z"]
    
    if mode == 0:
        # for i in range(p_k.shape[0]):
        #     p_k[i,:] *= block["rescale_Pk_fz", "alpha_{}".format(i)]**2
        block["rescale_Pk_fz", "fz"] = np.array([block["rescale_Pk_fz", "alpha_{}".format(i)] for i in range(p_k.shape[0])])

    if mode == 1:
        # _, zmax, nknots = config
        zmax = config["zmax"]
        nknots = config["nknots"]
        z_s = np.linspace(0., zmax, nknots)
        s_z = np.zeros(nknots)
        for j in range(nknots):
            s_z[j] = block["rescale_pk_fz", "alpha_{}".format(j)]
        
        # fz = np.ones_like(z) + interp1d(z_s, s_z, kind='linear', bounds_error=False, fill_value=0.0)(z)
        #fz_s_ = CubicSpline(z_s, s_z, bc_type='clamped', extrapolate=False)
        # fz_s_ = CubicSpline(z_s, s_z, bc_type='not-a-knot', extrapolate=False)
        # fz_s = lambda z : 1. if z>zmax else fz_s_(z)
        #fz_s = np.vectorize(fz_s)

        fz_s = extraCublicSpline(zmax, z_s, s_z, bc_type=((2,0.),(1.,0.)), extrapolate=False)

        block["rescale_Pk_fz", "fz"] =  fz_s(z)

        if config["smoothness_prior_amplitude"] > 0.:
            block[names.likelihoods, "rescale_Pk_fz_smoothness_like"] = get_smoothness_prior(block, fz_s, zmax, config["smoothness_prior_scale"], config["smoothness_prior_amplitude"])
        else:
            block[names.likelihoods, "rescale_Pk_fz_smoothness_like"] = 0.
            
        # for i in range(p_k.shape[0]):
        #     p_k[i,:] *= fz[i]**2
        d_z = block["growth_parameters", "d_z"]
        f_z = block["growth_parameters", "f_z"]
        z_growth = block["growth_parameters", "z"]
        # idx_zmax = np.where(z_growth < zmax)[0]

        # #note: different f_z here: growth factor
        # block["growth_parameters", "d_z"][idx_zmax] *= fz_s(z[idx_zmax]) / fz_s(0.)
        # block["growth_parameters", "f_z"][idx_zmax] -= (1.+z[idx_zmax]) / fz_s(z[idx_zmax]) * fz_s.derivative(1)(z[idx_zmax])
        z = block["growth_parameters", "z"]
        block["growth_parameters", "d_z"] *= fz_s(z) / fz_s(0.)
        block["growth_parameters", "f_z"] -= (1.+z) / fz_s(z) * fz_s.derivative(z)


    if mode == 2:
        # _, neigen, ez = config
        neigen = config["neigen"]
        ez = config["ez"]
        fz = np.ones_like(z)
        for j in range(neigen):
            fz += block["rescale_Pk_fz", "beta_{}".format(j)] * ez[j](z)
        block["rescale_Pk_fz", "fz"] = fz
        
        # for i in range(p_k.shape[0]):
        #     p_k[i,:] *= fz[i]**2

    for i in range(p_k.shape[0]):
        p_k[i,:] *= block["rescale_Pk_fz", "fz"][i]**2
        # also multiply growth! by f(z) not squared
        # growth = growth *fz(z)/ fz(0)
        # how to change f? 
        # make a quick plot: 

    block["matter_power_lin", "p_k"] = p_k

    return 0

def cleanup(config):
    pass


class extraCublicSpline:
    def __init__(self, zmax, *args, **kwargs):
        self.cubicspline = CubicSpline(*args, **kwargs)
        self.zmax = zmax
        self.first_derivative = self.cubicspline.derivative(1)
        self.extrap = self.cubicspline(zmax)

    def __call__(self, z):
        return np.where(z<self.zmax, self.cubicspline(z), np.ones_like(z)*self.extrap)

    def derivative(self, z):
        return np.where(z<self.zmax, self.first_derivative(z), np.zeros_like(z))