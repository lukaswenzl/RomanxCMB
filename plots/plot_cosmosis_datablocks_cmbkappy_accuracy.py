import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

output = "cmblensing_accuracy/"
path_to_datablocks = "../../../6x2pt_Roman_SO_camb_plot/"
path_to_datablocks_ehu = "../../../6x2pt_Roman_SO_plot/"


path_to_datablocks_pipeline_camb = "../../../6x2pt_Roman_SO_camb/"
path_to_datablocks_pipeline_ehu = "../../../6x2pt_Roman_SO/"
path_to_datablocks_pipeline_EH99 = "../../../6x2pt_Roman_SO_EH99/"






def plot(path_to_datablocks, folder_prefix, name, label, factor = 0, absolute=False, ):
    Cl = np.loadtxt("{}/{}.txt".format(path_to_datablocks+folder_prefix,name))
    ell = np.loadtxt("{}/ell.txt".format(path_to_datablocks+folder_prefix))

    plt.xlabel(r"$\ell$")

    if (factor == 1):
        factor_density = ell*(ell+1.)/2./np.pi
    elif (factor == 2):
        factor_density = 4. /2./np.pi
    elif (factor == 3):
        #factor_density = ell*(ell+1.)
        factor_density = (ell*(ell+1.))**2 /2./np.pi

    else:
        factor_density = 1 # CMB already saved with factor included

    if(not absolute):
        plt.ylabel(r"$[L(L+1)]^2C_L^{\phi \phi}/2\pi$")
        plt.plot(ell,factor_density*(Cl),label=r'{}'.format(label))
    else:
        plt.ylabel(r"$[L(L+1)]^2C_L^{\phi \phi}/2\pi$")
        plt.plot(ell,factor_density*np.abs(Cl),label=r'{}'.format(label), alpha=0.5)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(3,10000)
    plt.legend()

    return ell,factor_density*(Cl)


#cmb
#plot("cmb_cl", "tt", "tt", factor=False)


#KK CMB lensing
plt.figure()
ellcamb, Cl_camb = plot("",                output, "camb_cl_phiphi",  "camb", factor=3)
ellcosmosis, Cl_cosmosis = plot(path_to_datablocks,"cmbkappa_cl", "bin_1_1", r"camb limber $(4 C_l^{\kappa \kappa} /2\pi)$", factor = 2)
ellcosmosis_ehu, Cl_cosmosis_ehu = plot(path_to_datablocks_ehu,"cmbkappa_cl", "bin_1_1", r"ehu limber $(4 C_l^{\kappa \kappa} /2\pi)$", factor = 2)

ellcosmosis_pipeline, Cl_cosmosis_pipeline = plot(path_to_datablocks_pipeline_camb,"cmbkappa_cl", "bin_1_1", r"camb limber $(4 C_l^{\kappa \kappa} /2\pi)$", factor = 2)
ellcosmosis_pipeline_ehu, Cl_cosmosis_pipeline_ehu = plot(path_to_datablocks_pipeline_ehu,"cmbkappa_cl", "bin_1_1", r"ehu limber $(4 C_l^{\kappa \kappa} /2\pi)$", factor = 2)
ellcosmosis_pipeline_EH99, Cl_cosmosis_pipeline_EH99 = plot(path_to_datablocks_pipeline_EH99,"cmbkappa_cl", "bin_1_1", r"ehu limber $(4 C_l^{\kappa \kappa} /2\pi)$", factor = 2)

# plt.savefig(output+"cmblensing_plot_naive_cmbkappa.png")
# plt.savefig(output+"cmblensing_plot.png")
#plt.savefig(output+"cmblensing_plot_both_nonlinear.png")
plt.savefig(output+"cmblensing_plotRAW.png")

plt.figure()

plt.plot(ellcosmosis, np.interp(ellcosmosis, ellcamb, Cl_camb)/Cl_cosmosis , label="camb/ cosmosis camb")
plt.plot(ellcosmosis_ehu, Cl_cosmosis_ehu/np.interp(ellcosmosis_ehu, ellcosmosis, Cl_cosmosis) , label="cosmosis_ehu/ cosmosis camb")


plt.plot(ellcosmosis_pipeline_ehu, Cl_cosmosis_pipeline_ehu/np.interp(ellcosmosis_pipeline_ehu, ellcosmosis_pipeline, Cl_cosmosis_pipeline) , label="cosmosis_ehu/ cosmosis camb (pipeline)")
plt.plot(ellcosmosis_pipeline_EH99, Cl_cosmosis_pipeline_EH99/np.interp(ellcosmosis_pipeline_EH99, ellcosmosis_pipeline, Cl_cosmosis_pipeline) , label="cosmosis EH99/ cosmosis camb (pipeline)")


#plt.plot(ellcosmosis, Cl_cosmosis/np.interp(ellcosmosis, ellcamb, Cl_camb) , label="cosmosis/camb")
#plt.plot(ellcosmosis_ehu, Cl_cosmosis_ehu/np.interp(ellcosmosis_ehu, ellcamb, Cl_camb) , label="cosmosis_ehu/camb")

plt.ylabel(r"ratio of $C_l$")
plt.xlabel(r"$\ell$")
plt.legend()
plt.savefig(output+"cmbleninsg_plot_ratioRAW.png")



