import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

output = "example/"
path_to_datablocks = "../../../6x2pt_WFIRST_SO_plot_naive_cmbkappa/"

def plot(folder_prefix, name, label, factor = True, absolute=False):
    Cl = np.loadtxt("{}/{}.txt".format(path_to_datablocks+folder_prefix,name))
    ell = np.loadtxt("{}/ell.txt".format(path_to_datablocks+folder_prefix))

    plt.xlabel(r"$\ell$")

    if (factor == True):
        factor_density = ell*(ell+1.)/2./np.pi
    else:
        factor_density = 1 # CMB already saved with factor included

    if(not absolute):
        plt.ylabel(r"$\ell (\ell+1) C_l^{XY} / 2 \pi $")
        plt.plot(ell,factor_density*(Cl),label=r'${}$'.format(label))
    else:
        plt.ylabel(r"$\ell (\ell+1) |C_l^{XY}| / 2 \pi $")
        plt.plot(ell,factor_density*np.abs(Cl),label=r'${}$'.format(label))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(3,10000)
    plt.legend()


#cmb
plot("cmb_cl", "tt", "tt", factor=False)
plot("cmb_cl", "ee", "ee", factor=False)
plot("cmb_cl", "bb", "bb", factor=False)
plot("cmb_cl", "te", "te", factor=False)
plt.title("cmb")
plt.savefig(output+"cmb.png")


#KK CMB lensing
plt.figure()
plot("cmbkappa_cl", "bin_1_1", "\kappa \kappa")

plt.savefig(output+"cmblensing.png")


#galaxy
plt.figure()
plot("galaxy_cl", "bin_3_2", "g3 g2")
plot("galaxy_cl", "bin_4_3", "g4 g3")

plt.savefig(output+"galaxy_galaxy.png")


#shear
plt.figure()
plot("shear_cl", "bin_3_2", "\gamma 3 \gamma 2")
plot("shear_cl", "bin_4_3", "\gamma 4 \gamma 3")

plt.savefig(output+"shear.png")

#intrinsic
plt.figure()
plot("shear_cl_ii", "bin_3_2", "i 3 i 2")
plot("shear_cl_ii", "bin_4_3", "i 4 i 3")

plt.savefig(output+"intrinsic.png")

#galaxy shear
plt.figure()
plot("galaxy_shear_cl", "bin_3_2", "g 3 \gamma 2")
plot("galaxy_shear_cl", "bin_4_3", "g 4 \gamma 3")

plt.savefig(output+"galaxy_shear.png")

#shear intrinsic
plt.figure()
plot("shear_cl_gi", "bin_3_2", "\gamma 3 i 2", absolute=True)
plot("shear_cl_gi", "bin_4_3", "\gamma 4 i 3", absolute=True)
plt.savefig(output+"shear_intrinsic.png")

#galaxy intrinsic
plt.figure()
plot("galaxy_intrinsic_cl", "bin_3_2", "g 3 i 2", absolute=True)
plot("galaxy_intrinsic_cl", "bin_4_3", "g 4 i 3", absolute=True)
plt.savefig(output+"galaxy_intrinsic.png")


#galaxy cmb lensing
plt.figure()
plot("galaxy_cmbkappa_cl", "bin_3_1", "g 3 \kappa", absolute=True)
plot("galaxy_cmbkappa_cl", "bin_4_1", "g 4 \kappa", absolute=True)
plt.yscale('log')
plt.savefig(output+"galaxy_cmblensing.png")

#shear cmb lensing
plt.figure()
plot("shear_cmbkappa_cl", "bin_3_1", "\gamma 3 \kappa", absolute=True)
plot("shear_cmbkappa_cl", "bin_4_1", "\gamma 4 \kappa", absolute=True)
plt.yscale('log')
plt.savefig(output+"shear_cmblensing.png")


#intrinsic cmb lensing
plt.figure()
plot("intrinsic_cmbkappa_cl", "bin_3_1", "i 3 \kappa", absolute=True)
plot("intrinsic_cmbkappa_cl", "bin_4_1", "i 4 \kappa", absolute=True)
plt.savefig(output+"intrinsic_cmblensing.png")
