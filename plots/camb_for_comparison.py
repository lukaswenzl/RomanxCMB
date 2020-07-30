import camb
from camb import model

import numpy as np


#calculate accurate Cl phiphi with camb
pars = camb.CAMBparams()
#h0 0.6727, 
#omb 0.0491685
#omc 0.3156-0.0491685-omnu = 0.3156-0.0491685- (0.0006155/0.6727/0.6727) =  0.265071355, omch2 = 0.265071355 * 0.6727*0.6727 = 0.1199514918
#omnuh2 0.0006155
#ombh2 0.0491685*0.6727*0.6727 =0.02224998972
#omch2 (0.3156-0.0491685)*0.6727*0.6727=0.1205669918
#neutrinos: 2 massless and 1 massive
#mnu = omnuh2 * 93.14 = 0.05732767 (no other factor since only one is assumed to be massive)

#for testing: fnu = 0.04
#omnuh2 = 0.005712679261
#mnu = 0.005712679261* 93.14 = 0.5320789464

#pars.set_cosmology(H0=69, ombh2=0.022852799999999996, omch2=0.1191472, num_massive_neutrinos=3, mnu=0.06, nnu=3.046, omk=0, tau=0.0697186, YHe=0.245341)#, num_massive_neutrinos=3, mnu=.06, nnu=0.046, YHe=0.245341)
#wrong interpretation of omega m: pars.set_cosmology(H0=67.27, ombh2=0.02224998972, omch2=0.1205669918, num_massive_neutrinos=1, mnu=0.05732767, nnu=3.046, omk=0, tau=0.08, YHe=0.245341)#, num_massive_neutrinos=3, mnu=.06, nnu=0.046, YHe=0.245341)
pars.set_cosmology(H0=67.27, ombh2=0.02224998972, omch2=0.1199514918, num_massive_neutrinos=1, mnu=0.05732767, nnu=3.046, omk=0, tau=0.08, YHe=0.245341)#, num_massive_neutrinos=3, mnu=.06, nnu=0.046, YHe=0.245341)

#pars.set_cosmology(H0=67.27, ombh2=0.02224998972, omch2=0.1205669918, num_massive_neutrinos=1, mnu=0.5320789464, nnu=3.046, omk=0, tau=0.08, YHe=0.245341)

pars.InitPower.set_params(As=2.1e-9, ns=0.9645, nrun=0., pivot_scalar=0.05)
pars.set_matter_power(kmax=10.0)
pars.set_for_lmax(10000, lens_potential_accuracy=4)
pars.NonLinear = model.NonLinear_both
pars.NonLinearModel.set_params(halofit_version='takahashi')
results = camb.get_results(pars)
cl = results.get_lens_potential_cls(lmax=10000, raw_cl=True)
ell = np.arange(len(cl[:,0]))
#plt.plot(cl[:,0]*(ell*(ell+1.)/2.)**2, label='pycamb (lens_potential_accuracy=1)')

np.savetxt("cmblensing_accuracy/camb_cl_phiphi.txt",cl[:,0])
#note this is actually Clphiphi, usually plotted with (L(L+1))**2 /2./pi factor that still needs to be added
np.savetxt("cmblensing_accuracy/ell.txt",ell)

#results.get_nonlinear_matter_power_spectrum()
pars.set_matter_power(redshifts=[1100., 10., 4.0, 0.5, 0.], kmax=10.0)
kh_nonlin, z, pk_takahashi = results.get_nonlinear_matter_power_spectrum(params=pars)
np.savetxt("matterpower_accuracy/k_h.txt",kh_nonlin)
np.savetxt("matterpower_accuracy/z.txt",z)
np.savetxt("matterpower_accuracy/p_k.txt",pk_takahashi)
np.savetxt("matterpower_accuracy/sigma8.txt",[results.get_sigma8()[-1]])


print(results.get_sigma8()[-1])


#####test of limber, is more consistent with limber approx code in cosmosis 
# nz = 100 #number of steps to use for the radial/redshift integration
# kmax=10  #kmax to use
# #First set up parameters as usual
# pars = camb.CAMBparams()
# pars.set_cosmology(H0=67.27, ombh2=0.02224998972, omch2=0.1199514918, num_massive_neutrinos=1, mnu=0.05732767, nnu=3.046, omk=0, tau=0.08, YHe=0.245341)#, num_massive_neutrinos=3, mnu=.06, nnu=0.046, YHe=0.245341)

# pars.InitPower.set_params(As=2.1e-9, ns=0.9645, nrun=0., pivot_scalar=0.05)

# #For Limber result, want integration over \chi (comoving radial distance), from 0 to chi_*.
# #so get background results to find chistar, set up a range in chi, and calculate corresponding redshifts
# results= camb.get_background(pars)
# chistar = results.conformal_time(0)- results.tau_maxvis
# chis = np.linspace(0,chistar,nz)
# zs=results.redshift_at_comoving_radial_distance(chis)
# #Calculate array of delta_chi, and drop first and last points where things go singular
# dchis = (chis[2:]-chis[:-2])/2
# chis = chis[1:-1]
# zs = zs[1:-1]

# #Get the matter power spectrum interpolation object (based on RectBivariateSpline). 
# #Here for lensing we want the power spectrum of the Weyl potential.
# PK = camb.get_matter_power_interpolator(pars, nonlinear=True, 
#     hubble_units=False, k_hunit=False, kmax=kmax,
#     var1=model.Transfer_Weyl,var2=model.Transfer_Weyl, zmax=zs[-1])

# #Get lensing window function (flat universe)
# win = ((chistar-chis)/(chis**2*chistar))**2
# #Do integral over chi
# ls = np.arange(2,3000+1, dtype=np.float64)
# cl_kappa=np.zeros(ls.shape)
# w = np.ones(chis.shape) #this is just used to set to zero k values out of range of interpolation
# for i, l in enumerate(ls):
#     k=(l+0.5)/chis
#     w[:]=1
#     w[k<1e-4]=0
#     w[k>=kmax]=0
#     cl_kappa[i] = np.dot(dchis, w*PK.P(zs, k, grid=False)*win/k**4)
# cl_kappa*= (ls*(ls+1))**2
# cl_phiphi = 4. * cl_kappa / ((ls*(ls+1))**2)

# np.savetxt("cmblensing_accuracy/camb_cl_phiphi.txt",cl_phiphi)
# #note this is actually Clphiphi, usually plotted with (L(L+1))**2 /2./pi factor that still needs to be added
# np.savetxt("cmblensing_accuracy/ell.txt",ls)
