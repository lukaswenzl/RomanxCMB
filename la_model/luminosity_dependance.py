#Author Lukas Wenzl
#Script to calculate the luminosity functions used for the krause_eifler_blazek Intrinsic Alignment model
import numpy as np
from scipy import integrate
#goal is to calculate phi_all and phi_red
#input parameters




#utility functions
def Magnitude_to_Luminosity(M):
    #convert absolute Magnitudes to Luminosity
    L0 =  3.0128e28 #W
    return L0 * 10**(-0.4*M)

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
def luminosity_integral(lower_lim, phi, L0=None,beta=None):
    #phi depends on z and L, in order [L,z]
    L = phi["L"]
    z = phi["z"]
    if(L0 != None and beta != None):#TEST THIS: todo
        tmp = (L/L0)**beta 
        _, L_over_L0_mesh = np.meshgrid(z, tmp)
        integrand = phi["phi"]*L_over_L0_mesh
    else:
        integrand = phi["phi"]
    integral = np.zeros(len(z))
    print(lower_lim)
    for i in range(len(z)):
        f = func_container(integrand[:,i], L)
        print(i)
        integral[i]= integrate.quad(f, np.log(lower_lim[i]), np.log(L[-1]), epsrel=1e-6)[0] ##be careful with the upper limit. Only works if phi is sampled to high enough L.
        #relative tolerance is tested with gamma function comparison, see test(), 
    result = {"result":integral, "z":phi["z"]}
    return result



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

def calculate_luminosity_function(z,phi_star_0,M_star, alpha,P,Q,h, smallest_logL=31):
    #$\phi(L, z)=\phi^{*}(z)\left(\frac{L}{L^{*}(z)}\right)^{\alpha} \exp \left(-\frac{L}{L^{*}(z)}\right)$
    #with $\phi^{*}(z)=\phi_{0}^{*} 10^{0.4 P z}$

    phi_star = phi_star_0* 10**(0.4*P*z)  #(h/Mpc)^3 !!!
  
    L = np.logspace(smallest_logL,40,100)#W #31 to 40; 100 samples is enough for to reach 1e-6 accuracy for the integral
    #on the high end the exponential takes over and the values are much smaller than machine precision for 10^40W
    #on the low end consider the faintest dwarfs currently know are 10^5 L_sun https://arxiv.org/abs/1901.05465, this is of order 10^31 W
    L_star = Magnitude_to_Luminosity(M_star -5*np.log10(h)- Q*(z-0.1)) #typically 10^36- 10^37W (order of magnitude of Milkyway luminosity) ##WRONG:np.log(h)?!? Why not +
    _, L_mesh = np.meshgrid(z, L)
    phi_star_mesh, _ = np.meshgrid(phi_star, L)
    L_star_mesh, _ = np.meshgrid(L_star, L)

    phi_values = phi_star_mesh *(L_mesh/L_star_mesh)**alpha * np.exp(-L_mesh/L_star_mesh)
    phi = {"phi":phi_values,  "z":z, "L":L}
    return phi


def calculate_L_lim(m_lim, D_L,k_corr,e_corr, h, galaxy_type="Sa"):
    # $M_{\lim }\left(z, m_{\lim }\right)=m_{\lim }-\left(5 \log _{10} \frac{D_{\mathrm{L}}(z)}{\mathrm{Mpc} / \mathrm{h}}+25+k(z)\right)$
    # needs the limiting magnitude of the survey and the luminosity distance function for cosmology D_L(z) [Mpc/h]
    
    z = D_L["z"]
    print("right now we are using Sa galaxy k and e correction. This might be a pretty bad choice, only for code testing purposes!")
    k = resample(k_corr[galaxy_type],k_corr["z"], z)+resample(e_corr[galaxy_type],e_corr["z"], z)# still unclear which column to use!!!!!!!!!!!!!!!!!!!!!!!!
    M_lim = m_lim - ( 5.*np.log10(np.array(D_L["D_L"])/h)+25.+k)#factor h!! is this correct? TODO
    L_lim_values = Magnitude_to_Luminosity(M_lim)
    if(z[0]==0):
        L_lim_values[0] = 1e30
    L_lim = {"L_lim":L_lim_values,"z":z}
    #import pdb; pdb.set_trace()
    return L_lim
    
def red_fraction(L_lim_red, L_lim_all, phi_red, phi_all):
    #$f_{\mathrm{red}}\ (m_{\lim} , z  ) = \frac{\int_{L\ (m_{\lim} , z)}^\infty d L    \phi_{\mathrm{red}}(L, z)}{\int_{L\ (m_{\lim} , z)}^\infty  dL \phi_{\mathrm{all}}(L, z)  }$
    #only for testing. Not used in final pipeline.

    red = luminosity_integral(L_lim_red["L_lim"], phi_red) 
    all = luminosity_integral(L_lim_all["L_lim"], phi_all) 
    f_red = {"f_red": red["result"]/all["result"], "z":red["z"]}
    return f_red,red,all

def integrate_luminosity_dependence(A_L0, phi_red, phi_all,L_lim_red, L_lim_all,M0,beta):
    #integrate luminosity dependence and multiply with red fraction (saves on calculation since one integral cancels)
    #red = luminosity_integral(lower_lim, phi_red) 
    L0 = Magnitude_to_Luminosity(M0)
    lum_scaling = luminosity_integral(L_lim_red["L_lim"], phi_red,L0,beta) 
    all = luminosity_integral(L_lim_all["L_lim"], phi_all,L0,beta) 

    # #resample in redshift
    # print(len(all["z"]))
    # print(len(A_L0["z"]))
    # all["result"] = resample(all["result"], all["z"], A_L0["z"], extrapolate_linearly_at_highz=False)
    # all["z"] = A_L0["z"]
    # lum_scaling["result"] = resample(lum_scaling["result"], lum_scaling["z"], A_L0["z"], extrapolate_linearly_at_highz=False)
    # lum_scaling["z"] = A_L0["z"]

    A_mlim_values = A_L0["A_L0"]* lum_scaling["result"]/all["result"]
    A_mlim = {"A_mlim":A_mlim_values, "z":A_L0["z"]} #might have to resample in z at some point here!
    return A_mlim


def test():
    h = 0.68 
    k_corr, e_corr = setup_luminosity_dependance(band = "r", pathtok_e_corr="k_e_correction/")
    #check integral, when setting lower limit to 0 we should get \gamma(alpha+1)= alpha (when normalizing the unit constants to 1)
    z = np.linspace(0,3,20)
    M_star = -20.70
    Q = 0.7
    alpha = 1
    phi = calculate_luminosity_function(z,phi_star_0=1,M_star=M_star, alpha=alpha,P=0,Q=Q,h=h, smallest_logL=1)
    lower_limit = np.zeros_like(phi["z"])
    result = luminosity_integral(lower_limit, phi, L0=None,beta=None)
    print("mean error from \gamma(2) = 1 test: "+str(np.mean(result["result"]-1)))

    alpha = 0.4
    goal = 2.21815954375768822
    phi = calculate_luminosity_function(z,phi_star_0=1,M_star=M_star, alpha=alpha,P=0,Q=Q,h=h, smallest_logL=1)
    lower_limit = np.zeros_like(phi["z"])
    result = luminosity_integral(lower_limit, phi, L0=None,beta=None)
    print("mean error from \gamma(0.4) = 2.21815954375768822 test: "+str(np.mean(result["result"]-goal)))

    #from a fiducial cosmology run with CAMB we prepared a test luminosity distance function
    D_L_values = [    0.        ,   672.57723759,  1455.48030101,  2331.48753462,
        3285.68795113,  4305.61805689,  5381.09254979,  6503.78506304,
        7667.02475848,  8865.42007404, 10094.600975  , 11350.98890521,
       12631.60526837, 13933.98750759, 15256.06520632, 16596.09059777,
       17952.55788669, 19324.16005975, 20709.76860546, 22108.39459514,
       23519.17381332, 24941.32636342, 26374.15655182, 27817.04384621,
       29269.42850506, 30730.81106704, 32200.72303016, 33678.74138781,
       35164.4798896 , 36657.58331913]
    z = np.linspace(0,4,30)
    if(z[0]==0):
        z = z[1:]
        D_L_values = D_L_values[1:]
    D_L = {"D_L":D_L_values,"z":z}
    #testrun calculating A and red fraction for GAMA case in Krause et al 2016
    #all galaxies case, GAMA case
    phi_star_0 = 9.4e-3 #(h/Mpc)^3 
    M_star = -20.70  #  +5*np.log(0.7)#weird h factor but I think I do not have to add that here.
    alpha = -1.23
    P=1.8
    Q = 0.7
    phi_all = calculate_luminosity_function(D_L["z"],phi_star_0,M_star, alpha,P,Q, h)###
    phi_star_0 = 1.1e-2 #(h/Mpc)^3 
    M_star = -20.34 #+5*np.log(0.7)
    alpha = -0.57
    P=-1.2 
    Q = 1.8
    
    phi_red = calculate_luminosity_function(D_L["z"],phi_star_0,M_star, alpha,P,Q, h)###

    m_lim = 27.5 
    print("testing with m_lim = 27.5 in r band for Roman HLS survey as in Krause et al 2016")
    L_lim_all = calculate_L_lim(m_lim, D_L, k_corr, h, e_corr,galaxy_type="Sa")
    L_lim_red = calculate_L_lim(m_lim, D_L, k_corr, h, e_corr,galaxy_type="E")

    print("here are the limiting Luminosities, they should be between 1e30 and 1e40 and increasing")
    print(L_lim_all["L_lim"])
    f_red,red,all = red_fraction(L_lim_red,L_lim_all , phi_red, phi_all)

    #testing for lower limiting mag: 
    m_lim = 24.5
    L_lim = calculate_L_lim(m_lim, D_L, k_corr, h, e_corr,galaxy_type="E")#WRONG
    L_lim_all = calculate_L_lim(m_lim, D_L, k_corr, h, e_corr,galaxy_type="Sa")
    L_lim_red = calculate_L_lim(m_lim, D_L, k_corr, h, e_corr,galaxy_type="E")
    f_red_gamma_24_5,_,_ = red_fraction(L_lim_red,L_lim_all , phi_red, phi_all)

    #deep2 survey
    phi_star_0 = 9.4e-3 #(h/Mpc)^3 
    M_star = -20.70  #  +5*np.log(0.7)#weird h factor but I think I do not have to add that here.
    alpha = -1.23
    P=-0.3
    Q=1.23
    phi_all = calculate_luminosity_function(D_L["z"],phi_star_0,M_star, alpha,P,Q, h)
    phi_star_0 = 1.1e-2 #(h/Mpc)^3 
    M_star = -20.34 #+5*np.log(0.7)
    alpha = -0.57
    P=-1.15
    Q=1.20
    phi_red = calculate_luminosity_function(D_L["z"],phi_star_0,M_star, alpha,P,Q ,h)
    m_lim = 27.5
    L_lim_all = calculate_L_lim(m_lim, D_L, k_corr, e_corr, h, galaxy_type="Sa")##THIS HAS A HUGE IMPACT??!
    L_lim_red = calculate_L_lim(m_lim, D_L, k_corr, e_corr, h, galaxy_type="E")##THIS HAS A HUGE IMPACT??!
    #L_lim["L_lim"] = L_lim["L_lim"]*5#WRONG
    f_red_deep2_27_5,_,_ = red_fraction(L_lim_red,L_lim_all , phi_red, phi_all)
    m_lim = 24.5
    L_lim_all = calculate_L_lim(m_lim, D_L, k_corr, e_corr, h, galaxy_type="Sa")##THIS HAS A HUGE IMPACT??!
    L_lim_red = calculate_L_lim(m_lim, D_L, k_corr, e_corr, h, galaxy_type="E")##THIS HAS A HUGE IMPACT??!
    #L_lim["L_lim"] = L_lim["L_lim"]*5#WRONG
    f_red_deep2_24_5,_,_ = red_fraction(L_lim_red,L_lim_all , phi_red, phi_all)


    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.figure(figsize=(7,4))
    plt.plot(f_red["z"],f_red["f_red"],"--", color="red", label="GAMA for Roman")
    plt.plot(f_red_gamma_24_5["z"],f_red_gamma_24_5["f_red"], color="red", label="GAMA for m_lim 24.5")
    plt.plot(f_red_deep2_27_5["z"],f_red_deep2_27_5["f_red"], "--", color="black", label="DEEP2 for Roman")
    plt.plot(f_red_deep2_24_5["z"],f_red_deep2_24_5["f_red"], color="black",label="DEEP2 for m_lim 24.5")

    plt.xlabel("z")
    plt.ylabel("fraction of red galaxies")
    plt.ylim(0.,0.3)
    plt.plot( [0.2],[0.055],"o", label="Krause et al 2016 peak")
    plt.legend()
    plt.xlim(0,3.5)
    plt.savefig("test_fraction_of_red_galaxies.png")

    plt.figure()
    plt.plot(all["z"],all["result"], "--", label="all")
    plt.plot(red["z"],red["result"], "--", label="red")
    plt.xlabel("z")
    plt.legend()
    plt.savefig("test_number_of_galaxies_vs_redshift.png")

    #Luminosity function looks accurate
    # plt.figure()
    # plt.plot(D_L["z"], D_L["D_L"])
    # plt.yscale("log")
    # plt.savefig("test.png")


    import pdb; pdb.set_trace()

if __name__ == '__main__':
    test()
    import pdb; pdb.set_trace()