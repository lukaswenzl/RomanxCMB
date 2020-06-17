import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

output = "cmblensing_accuracy/"
# path_to_datablocks = "../../../6x2pt_WFIRST_SO_plot_naive_cmbkappa/"
#path_to_datablocks = "../../../6x2pt_WFIRST_SO_plot_better_cmbkappa_nosigma8rescale/"
path_to_datablocks = "../../../6x2pt_WFIRST_SO_plot/"




def plot(folder_prefix, name, label, factor = 0, absolute=False):
    Cl = np.loadtxt("{}/{}.txt".format(path_to_datablocks+folder_prefix,name))
    ell = np.loadtxt("{}/ell.txt".format(path_to_datablocks+folder_prefix))

    plt.xlabel(r"$\ell$")

    if (factor == 1):
        factor_density = ell*(ell+1.)/2./np.pi
    elif (factor == 2):
        factor_density = 4. /2./np.pi
    elif (factor == 3):
        factor_density = ell*(ell+1.)
    else:
        factor_density = 1 # CMB already saved with factor included

    if(not absolute):
        plt.ylabel(r"$[L(L+1)]^2C_L^{\phi}/2\pi$")
        plt.plot(ell,factor_density*(Cl),label=r'{}'.format(label))
    else:
        plt.ylabel(r"$[L(L+1)]^2C_L^{\phi}/2\pi$")
        plt.plot(ell,factor_density*np.abs(Cl),label=r'{}'.format(label))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(3,10000)
    plt.legend()


#cmb
#plot("cmb_cl", "tt", "tt", factor=False)


#KK CMB lensing
plt.figure()
plot("cmb_cl", "pp", "camb", factor=3)
plot("cmbkappa_cl", "bin_1_1", "limber $(4 C_l^{\kappa \kappa} /2\pi)$", factor = 2)

# plt.savefig(output+"cmblensing_plot_naive_cmbkappa.png")
# plt.savefig(output+"cmblensing_plot.png")
#plt.savefig(output+"cmblensing_plot_both_nonlinear.png")
plt.savefig(output+"cmblensing_plot_test.png")



