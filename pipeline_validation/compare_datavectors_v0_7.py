import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd
import astropy.io.fits as fits
#from astropy import units as u
#from astropy import constants as const


plt.style.use(['science'])#,'no-latex'])

#load cosmolike interface tools
import sys
sys.path.append('../2pt_modified/')

import cosmolike_metadata
import twopoint

def loop_comosis_datavector(two_point_data, n_ell=20):
    cl = []
    ell = []
    bin1 = []
    bin2 = []
    names = []
    ell_idx = []
    for spectrum in two_point_data.spectra:
        cl.extend(spectrum.value)
        ell.extend(spectrum.angle)
        for b in spectrum.bin_pairs:
            for i in range(n_ell): #ell index starting from 0
                #galaxy_cl110 
                bin1.append(b[0])
                bin2.append(b[1])
                names.append(spectrum.name)
                ell_idx.append(i)
    return {"cl":cl, "ell":ell, "bin1":bin1, "bin2":bin2, "names":names}


filename = "../6x2pt_Roman_SO_v0_7_483f98a.fits"
two_point_data = twopoint.TwoPointFile.from_fits(filename)
cosmosis_ehu = loop_comosis_datavector(two_point_data)

filename = "../6x2pt_Roman_SO_camb_v0_7_483f98a.fits"
two_point_data = twopoint.TwoPointFile.from_fits(filename)
cosmosis_camb = loop_comosis_datavector(two_point_data)

cosmolike = cosmolike_metadata.rearange_cosmolike_datavec(filename, two_point_data.spectra, "../cosmolike_data/cov_indices_apr9.txt")



x = np.arange(len(cosmosis_ehu["cl"]))
fig, axs = plt.subplots(2,1, figsize=(9,3))

#plt.figure(figsize=(16,4))
fig.suptitle("datavector comparison")
axs[0].set_yscale("log")
axs[0].plot(x, cosmosis_ehu["cl"], label="cosmosis EH99")
axs[0].plot(x, cosmosis_camb["cl"], "--",label="cosmosis camb")
axs[0].plot(x, cosmolike["cl"], ":", label="cosmolike")


axs[0].text(300, 5e-12, "shear")
axs[0].axvline(x=1100)
axs[0].text(1105, 5e-12, "galaxy x shear")
axs[0].axvline(x=1740)
axs[0].text(1745, 5e-7, "shear x \ncmb $\kappa$")
axs[0].axvline(x=1940)
axs[0].text(1945, 5e-12, "galaxy")
axs[0].axvline(x=2140)
axs[0].text(2145, 5e-12, "galaxy x \ncmb $\kappa$")
axs[0].axvline(x=2340)
axs[0].text(2340, 5e-12, "$\kappa$")

axs[0].legend(ncol=2)
axs[0].set_ylabel("$C_l$ (times l factors)")
#axs[0].set_xlabel("Index of datavector")
axs[0].set_xlim(0,2370)

#plt.savefig("v0_5_compare_datavector_logspaceRAW.png")

##########################
#plt.figure(figsize=(16,4))
#plt.title("Ratio of datavectors")

axs[1].plot(x, cosmosis_ehu["cl"]/cosmolike["cl"], label="cosmosis EH99/cosmolike")
#axs[1].plot(x, cosmosis_camb["cl"]/cosmolike["cl"],"--", label="cosmosis camb/cosmolike")
axs[1].plot(x, np.array(cosmosis_camb["cl"])/cosmosis_ehu["cl"],label="cosmosis camb/cosmosis EH99", linewidth=2, color="grey", alpha=0.4)




axs[1].set_ylabel("$C_l$ ratio")
axs[1].set_xlabel("Index of datavector")
axs[1].legend(loc=(0.2,0.52))
axs[1].set_xlim(0,2370)#making it a little wider so the kappa fits

fig.savefig("v0_7_compare_datavectorRAW.pdf")


# # #compare high and low resolution sampling of ehu. I want to make sure that by reducing the sampling the Cl don't change dramatically
# filename = "../6x2pt_Roman_SO_251nz_v0_6.fits"
# two_point_data = twopoint.TwoPointFile.from_fits(filename)
# cosmosis_ehu_highsampling = loop_comosis_datavector(two_point_data)

# plt.figure(figsize=(9,2))

# #plt.figure(figsize=(16,4))
# plt.title("ehu sampling in z")
# plt.plot(x, np.array(cosmosis_ehu["cl"])/cosmosis_ehu_highsampling["cl"], label="cosmosis ehu nz=101 / nz=301")

# plt.legend()
# plt.savefig("v0_6_compare_z_sampling_for_ehuRAW.pdf")