import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

output = "matterpower_accuracy/"
path_to_datablocks = "../../../6x2pt_Roman_SO_camb/" #_plot
path_to_datablocks_ehu = "../../../6x2pt_Roman_SO/" #_plot
path_to_datablocks_EH99 = "../../../6x2pt_Roman_SO_EH99/"




def plot(path_to_datablocks, folder_prefix, name, label, redshift=0.):
    Pk = np.loadtxt("{}/P_k.txt".format(path_to_datablocks+folder_prefix,name))
    k = np.loadtxt("{}/k_h.txt".format(path_to_datablocks+folder_prefix))
    z = np.loadtxt("{}/z.txt".format(path_to_datablocks+folder_prefix))

    Pk_resamplez = interp1d(z, Pk, axis=0)



    plt.xlabel(r"$k/h$ Mpc")


    plt.ylabel(r"$P_k (z= {})$".format(redshift))

    plt.plot(k,(Pk_resamplez(redshift)),"--",label=r'{}'.format(label))

    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(5e-6,1000)
    plt.legend()

    return Pk_resamplez, k


#cmb
#plot("cmb_cl", "tt", "tt", factor=False)


#matter power spectrum, redshift 0
plt.figure()
plt.title("z = 0")
Pk_camb, kcamb = plot("",                output, "",  "camb")
Pk_cosmosis_camb, k_cosmosis_camb = plot(path_to_datablocks,"matter_power_nl", "", r"cosmosis camb")
Pk_cosmosis_ehu, k_cosmosis_ehu = plot(path_to_datablocks_ehu,"matter_power_nl", "", r"cosmosis ehu")
Pk_cosmosis_EH99, k_cosmosis_EH99 = plot(path_to_datablocks_EH99,"matter_power_nl", "", r"cosmosis EH99")





plt.savefig(output+"matterpower_plot_z0RAW.pdf")

plt.figure()
plt.title("z = 0")
plt.plot(kcamb, Pk_camb(0.)/np.interp(kcamb,k_cosmosis_camb, Pk_cosmosis_camb(0.)) , label="camb/ cosmosis camb")
plt.plot(k_cosmosis_camb, np.interp(k_cosmosis_camb,k_cosmosis_ehu, Pk_cosmosis_ehu(0.))/Pk_cosmosis_camb(0.) , label="cosmosis_ehu/cosmosis camb")
plt.plot(k_cosmosis_camb, np.interp(k_cosmosis_camb,k_cosmosis_EH99, Pk_cosmosis_EH99(0.))/Pk_cosmosis_camb(0.) ,"--", label="cosmosis_EH99/cosmosis camb")


plt.xscale("log")
plt.ylabel(r"ratio of $P_k$")
plt.xlabel(r"$k/h$ Mpc")
plt.legend()
plt.savefig(output+"matterpower_ratio_plot_z0RAW.pdf")


#matter power spectrum, redshift 0.5
plt.figure()
plt.title("z = 0.5")
Pk_camb, kcamb = plot("",                output, "",  "camb", redshift=0.5)
Pk_cosmosis_camb, k_cosmosis_camb = plot(path_to_datablocks,"matter_power_nl", "", r"cosmosis camb", redshift=0.5)
Pk_cosmosis_ehu, k_cosmosis_ehu = plot(path_to_datablocks_ehu,"matter_power_nl", "", r"cosmosis ehu", redshift=0.5)
Pk_cosmosis_EH99, k_cosmosis_EH99 = plot(path_to_datablocks_EH99,"matter_power_nl", "", r"cosmosis EH99", redshift=0.5)



plt.savefig(output+"matterpower_plot_z0_5RAW.pdf")

plt.figure()
plt.title("z = 0.5")
plt.plot(kcamb, Pk_camb(0.5)/np.interp(kcamb,k_cosmosis_camb, Pk_cosmosis_camb(0.5)) , label="camb/ cosmosis camb")
plt.plot(k_cosmosis_camb, np.interp(k_cosmosis_camb,k_cosmosis_ehu, Pk_cosmosis_ehu(0.5))/Pk_cosmosis_camb(0.5) , label="cosmosis_ehu/cosmosis camb")
plt.plot(k_cosmosis_camb, np.interp(k_cosmosis_camb,k_cosmosis_EH99, Pk_cosmosis_EH99(0.5))/Pk_cosmosis_camb(0.5) ,"--", label="cosmosis_EH99/cosmosis camb")


plt.xscale("log")
plt.ylabel(r"ratio of $P_k$")
plt.xlabel(r"$k/h$ Mpc")
plt.legend()
plt.savefig(output+"matterpower_ratio_plot_z0_5RAW.pdf")


#matter power spectrum, redshift 4
plt.figure()
plt.title("z = 4")
Pk_camb, kcamb = plot("",                output, "",  "camb", redshift=4.)
Pk_cosmosis_camb, k_cosmosis_camb = plot(path_to_datablocks,"matter_power_nl", "", r"cosmosis camb", redshift=4.)
Pk_cosmosis_ehu, k_cosmosis_ehu = plot(path_to_datablocks_ehu,"matter_power_nl", "", r"cosmosis ehu", redshift=4.)
Pk_cosmosis_EH99, k_cosmosis_EH99 = plot(path_to_datablocks_EH99,"matter_power_nl", "", r"cosmosis EH99", redshift=4.)



plt.savefig(output+"matterpower_plot_z4RAW.pdf")

plt.figure()
plt.title("z = 4")
plt.plot(kcamb, Pk_camb(4.)/np.interp(kcamb,k_cosmosis_camb, Pk_cosmosis_camb(4.)) , label="camb/ cosmosis camb")
plt.plot(k_cosmosis_camb, np.interp(k_cosmosis_camb,k_cosmosis_ehu, Pk_cosmosis_ehu(4.))/Pk_cosmosis_camb(4.) , label="cosmosis_ehu/cosmosis camb")
plt.plot(k_cosmosis_camb, np.interp(k_cosmosis_camb,k_cosmosis_EH99, Pk_cosmosis_EH99(4.))/Pk_cosmosis_camb(4.) ,"--", label="cosmosis_EH99/cosmosis camb")


plt.xscale("log")
plt.ylabel(r"ratio of $P_k$")
plt.xlabel(r"$k/h$ Mpc")
plt.legend()
plt.savefig(output+"matterpower_ratio_plot_z4RAW.pdf")

#matter power spectrum, redshift 10
plt.figure()
plt.title("z = 10")
Pk_camb, kcamb = plot("",                output, "",  "camb", redshift=10.)
Pk_cosmosis_camb, k_cosmosis_camb = plot(path_to_datablocks,"matter_power_nl", "", r"cosmosis camb", redshift=10.)
Pk_cosmosis_ehu, k_cosmosis_ehu = plot(path_to_datablocks_ehu,"matter_power_nl", "", r"cosmosis ehu", redshift=10.)
Pk_cosmosis_EH99, k_cosmosis_EH99 = plot(path_to_datablocks_EH99,"matter_power_nl", "", r"cosmosis EH99", redshift=10.)



plt.savefig(output+"matterpower_plot_z10RAW.pdf")

plt.figure()
plt.title("z = 10")
plt.plot(kcamb, Pk_camb(10.)/np.interp(kcamb,k_cosmosis_camb, Pk_cosmosis_camb(10.)) , label="camb/ cosmosis camb")
plt.plot(k_cosmosis_camb, np.interp(k_cosmosis_camb,k_cosmosis_ehu, Pk_cosmosis_ehu(10.))/Pk_cosmosis_camb(10.) , label="cosmosis_ehu/cosmosis camb")
plt.plot(k_cosmosis_camb, np.interp(k_cosmosis_camb,k_cosmosis_EH99, Pk_cosmosis_EH99(10.))/Pk_cosmosis_camb(10.) ,"--", label="cosmosis_EH99/cosmosis camb")

plt.xscale("log")
plt.ylabel(r"ratio of $P_k$")
plt.xlabel(r"$k/h$ Mpc")
plt.legend()
plt.savefig(output+"matterpower_ratio_plot_z10RAW.pdf")

#matter power spectrum, redshift 1100
plt.figure()
plt.title("z = 1100")
Pk_camb, kcamb = plot("",                output, "",  "camb", redshift=1100.)
Pk_cosmosis_camb, k_cosmosis_camb = plot(path_to_datablocks,"matter_power_nl", "", r"cosmosis camb", redshift=1100.)
Pk_cosmosis_ehu, k_cosmosis_ehu = plot(path_to_datablocks_ehu,"matter_power_nl", "", r"cosmosis ehu", redshift=1100.)
Pk_cosmosis_EH99, k_cosmosis_EH99 = plot(path_to_datablocks_EH99,"matter_power_nl", "", r"cosmosis EH99", redshift=1100.)


plt.savefig(output+"matterpower_plot_z1100RAW.pdf")

plt.figure()
plt.title("z = 1100")
plt.plot(kcamb, Pk_camb(1100.)/np.interp(kcamb,k_cosmosis_camb, Pk_cosmosis_camb(1100.)) , label="camb/ cosmosis camb")
plt.plot(k_cosmosis_camb, np.interp(k_cosmosis_camb,k_cosmosis_ehu, Pk_cosmosis_ehu(1100.))/Pk_cosmosis_camb(1100.) , label="cosmosis_ehu/cosmosis camb")
plt.plot(k_cosmosis_camb, np.interp(k_cosmosis_camb,k_cosmosis_EH99, Pk_cosmosis_EH99(1100.))/Pk_cosmosis_camb(1100.) ,"--", label="cosmosis_EH99/cosmosis camb")


plt.xscale("log")
plt.ylabel(r"ratio of $P_k$")
plt.xlabel(r"$k/h$ Mpc")
plt.legend()
plt.savefig(output+"matterpower_ratio_plot_z1100RAW.pdf")


#demo 8 matter power spectrum, redshift 
plt.figure()
plt.title("z = 0")
Pk_cosmosis_camb, k_cosmosis_camb = plot("../../../output/demo8/","matter_power_lin", "", r"demo8 lin")
Pk_cosmosis_ehu, k_cosmosis_ehu = plot("../../../output/demo8/","matter_power_no_bao", "", r"demo8 ehu no bao")


plt.savefig(output+"demo8_linear_matterpower_plot_z0RAW.pdf")

plt.figure()
plt.title("z = 0")
plt.plot(k_cosmosis_ehu, 1/(np.interp(k_cosmosis_ehu,k_cosmosis_camb, Pk_cosmosis_camb(0.))/Pk_cosmosis_ehu(0.)) , label="ehu/camb (cosmosis)")
#plt.plot(kcamb, np.interp(kcamb,k_cosmosis_ehu, Pk_cosmosis_ehu(0.))/Pk_camb(0.) , label="demo9 ehu no bao")

plt.xscale("log")
plt.ylabel(r"ratio of $P_k$")
plt.xlabel(r"$k/h$ Mpc")
plt.legend()
plt.savefig(output+"demo8_linear_matterpower_ratio_plot_z0RAW.pdf")

#linear matter power spectrum, redshift 0
plt.figure()
plt.title("z = 0")
Pk_camb, kcamb = plot("",                output, "",  "camb")
Pk_cosmosis_camb, k_cosmosis_camb = plot(path_to_datablocks,"matter_power_lin", "", r"cosmosis camb")
Pk_cosmosis_ehu, k_cosmosis_ehu = plot(path_to_datablocks_ehu,"matter_power_lin", "", r"cosmosis ehu")


plt.savefig(output+"linear_matterpower_plot_z0RAW.pdf")

plt.figure()
plt.title("z = 0")
plt.plot(k_cosmosis_ehu, 1/(np.interp(k_cosmosis_ehu,k_cosmosis_camb, Pk_cosmosis_camb(0.))/Pk_cosmosis_ehu(0.) ), label="ehu/camb (cosmosis)")
#plt.plot(kcamb, np.interp(kcamb,k_cosmosis_ehu, Pk_cosmosis_ehu(0.))/Pk_camb(0.) , label="demo9 ehu no bao")

plt.xscale("log")
plt.ylabel(r"ratio of $P_k$")
plt.xlabel(r"$k/h$ Mpc")
plt.legend()
plt.savefig(output+"linear_matterpower_ratio_plot_z0RAW.pdf")
