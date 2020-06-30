import camb
from camb import model

import numpy as np


#calculate accurate Cl phiphi with camb
pars = camb.CAMBparams()
#h0 0.6727, STILL NEED TO MATCH THOSE
#omb 0.0491685
#omc 0.3156-0.0491685
#omnuh2 0.0006155
#ombh2 0.0491685*0.6727*0.6727=0.02224998972
#omch2 (0.3156-0.0491685)*0.6727*0.6727=0.1205669918
#mnu UNCLEAR HOW TO CALCULATE THAT FROM omnuh2 (will leave standard value for now)

#pars.set_cosmology(H0=69, ombh2=0.022852799999999996, omch2=0.1191472, num_massive_neutrinos=3, mnu=0.06, nnu=3.046, omk=0, tau=0.0697186, YHe=0.245341)#, num_massive_neutrinos=3, mnu=.06, nnu=0.046, YHe=0.245341)
pars.set_cosmology(H0=67.27, ombh2=0.02224998972, omch2=0.1205669918, num_massive_neutrinos=3, mnu=0.06, nnu=3.046, omk=0, tau=0.08, YHe=0.245341)#, num_massive_neutrinos=3, mnu=.06, nnu=0.046, YHe=0.245341)

pars.InitPower.set_params(As=2.1e-9, ns=0.9645)
pars.set_for_lmax(10000, lens_potential_accuracy=1)
pars.NonLinear = model.NonLinear_both
pars.NonLinearModel.set_params(halofit_version='takahashi')
results = camb.get_results(pars)
cl = results.get_lens_potential_cls(lmax=10000, raw_cl=True)
ell = np.arange(len(cl[:,0]))
#plt.plot(cl[:,0]*(ell*(ell+1.)/2.)**2, label='pycamb (lens_potential_accuracy=1)')

np.savetxt("cmblensing_accuracy/camb_cl_phiphi.txt",cl[:,0])
#note this is actually Clphiphi, usually plotted with (L(L+1))**2 /2./pi factor that still needs to be added
np.savetxt("cmblensing_accuracy/ell.txt",ell)

