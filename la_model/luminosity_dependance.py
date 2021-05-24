#Author Lukas Wenzl
#Script to calculate the luminosity functions used for the krause_eifler_blazek Intrinsic Alignment model
import numpy as np
from scipy import integrate
from scipy import special
#goal is to calculate phi_all and phi_red
#input parameters


#utility functions
def Magnitude_to_Luminosity(M):
    #convert absolute Magnitudes to Luminosity
    L0 =  3.0128e28 #W
    return L0 * 10.**(-0.4*M)

def resample(F, z_old, z_new, extrapolate_linearly_at_highz=True):
    if extrapolate_linearly_at_highz:
        F_new = np.interp(z_new, z_old, F)

        #extrapolate at linearly at high redshift from the highest 3 datapoints
        #print(z_old[-3:])
        fit_high = np.poly1d(np.polyfit(z_old[-3:], F[-3:], 1))
        mask_high = np.argwhere(z_new>np.max(z_old))
        F_new[mask_high] = fit_high(z_new[mask_high])

        return F_new
    else:
        return np.interp(z_new, z_old, F)

class func_container():
    def __init__(self, values, x):
        self.values_interp = values
        self.x_interp = x
        #self.calls = 0
    def __call__(self, logx):
        #to handle scales use log spacing
        x = np.exp(logx)
        #self.calls += 1
        return np.interp(x, self.x_interp, self.values_interp)

#integrate from lower cutoff to infinity
# def luminosity_integral(lower_lim, phi, L0=None,beta=None):
#     #phi depends on z and L, in order [L,z]
#     L = phi["L"]
#     z = phi["z"]
#     if(L0 != None and beta != None):#TEST THIS: todo
#         tmp = (L/L0)**beta 
#         _, L_over_L0_mesh = np.meshgrid(z, tmp)
#         integrand = phi["phi"]*L_over_L0_mesh
#     else:
#         integrand = phi["phi"]
#     integral = np.zeros(len(z))
#     print(lower_lim)
#     for i in range(len(z)):
#         f = func_container(integrand[:,i], L)
#         print(i)
#         integral[i]= integrate.quad(f, np.log(lower_lim[i]), np.log(L[-1]), epsrel=1e-6)[0] ##be careful with the upper limit. Only works if phi is sampled to high enough L.
#         #relative tolerance is tested with gamma function comparison, see test(), 
#     result = {"result":integral, "z":phi["z"]}
#     return result



#preloading files for faster cosmosis
def setup_luminosity_dependance(band = "r", pathtok_e_corr="k_e_correction/"):
    
    if(band == "r"):
        k_corr = np.loadtxt(pathtok_e_corr+"kcorr.dat", dtype= {'names': ('z', 'band', 'E', 'E2', 'Sa','Sc'), 'formats': (np.float, '|S2', np.float, np.float, np.float, np.float)})
        charar = np.chararray((len(k_corr)))
        charar[:] = 'r'
        k_corr = k_corr[charar == k_corr["band"]]
        
        e_corr = np.loadtxt(pathtok_e_corr+"ecorr.dat", dtype= {'names': ('z', 'band', 'E', 'E2', 'Sa','Sc'), 'formats': (np.float, '|S2', np.float, np.float, np.float, np.float)})
        charar = np.chararray((len(e_corr)))
        charar[:] = 'r'
        e_corr = e_corr[charar == e_corr["band"]]
    else: 
        print("Specified band not implemented")
        return None
    
    return k_corr, e_corr

class Luminosity_function():
    def __init__(self, phi_star_0, M_star, alpha, P, Q, z):
        self.phi_star_0 = phi_star_0
        self.M_star = M_star
        self.alpha = alpha
        self.P = P
        self.Q = Q
        self.phi_star = self.phi_star_0 * 10.**(0.4*self.P*z)
        self.L_star = Magnitude_to_Luminosity(self.M_star - self.Q*(z-0.1))     
        


def calculate_L_lim(m_lim, D_L,k_corr,e_corr, h, galaxy_type="Sa"):
    # $M_{\lim }\left(z, m_{\lim }\right)=m_{\lim }-\left(5 \log _{10} \frac{D_{\mathrm{L}}(z)}{\mathrm{Mpc} / \mathrm{h}}+25+k(z)\right)$
    # needs the limiting magnitude of the survey and the luminosity distance function for cosmology D_L(z) [Mpc/h]
    
    z = D_L["z"]
    k = resample(k_corr[galaxy_type],k_corr["z"], z)+resample(e_corr[galaxy_type],e_corr["z"], z)# choice in galaxy type for correction makes little difference
    distances = np.array(D_L["D_L"])/h
    distances = np.where(distances > 1., distances, 1.) #avoid taking log of 0
    M_lim = m_lim - ( 5.*np.log10(distances)+25.+k) # important factor of h! camb distance is in Mpc so need to divide by h here
    L_lim_values = Magnitude_to_Luminosity(M_lim)
    if(z[0]==0):
        L_lim_values[0] = 1e30
    L_lim = {"L_lim":L_lim_values,"z":z}
    #import pdb; pdb.set_trace()
    return L_lim


def gamma_incomplete(a, x):
    if(a>0):
        #gammaincc is the regularized upper incomplete gamma function 
        gamma_inc = special.gammaincc(a, x)  * special.gamma(a)
    else:
        #the implementation in scipy doesn't handle negative a
        #use recursion relation:
        #gamma (a+1, x) = a * gamma(a,x) + x**a * exp(-x)
        # -> gamma(a, x) = (gamma(a+1, x) - x**a * np.exp(-x))/a
        gamma_inc = (gamma_incomplete(a+1, x) - x**a * np.exp(-x))/a
    return gamma_inc

def luminosity_dependence(z, phi_red, L_lim_red, M0, beta):
    L0 = Magnitude_to_Luminosity(M0 - phi_red.Q * z)

    #C.3 in Joachimi 2011
    gam1 = gamma_incomplete(phi_red.alpha+beta+1,L_lim_red["L_lim"]/phi_red.L_star)
    gam2 = gamma_incomplete(phi_red.alpha+1,L_lim_red["L_lim"]/phi_red.L_star)

    LoverL0_beta = (phi_red.L_star / L0)**beta * gam1/gam2

    return LoverL0_beta

def red_fraction(z, L_lim_red, L_lim_all, phi_red, phi_all):

    gam_red = gamma_incomplete(phi_red.alpha+1,L_lim_red["L_lim"]/phi_red.L_star)
    gam_all = gamma_incomplete(phi_all.alpha+1,L_lim_all["L_lim"]/phi_all.L_star)

    #import pdb; pdb.set_trace()
    f_red = phi_red.phi_star / phi_all.phi_star * gam_red / gam_all
    return f_red

