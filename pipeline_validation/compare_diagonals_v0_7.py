import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd
import astropy.io.fits as fits
#from astropy import units as u
#from astropy import constants as const

plt.style.use(['science','no-latex'])

print("not final")

def load_cov(filename, name):
    hdul = fits.open(filename)
    print(hdul.info())
    cov = hdul[1].data
    cov_log = np.log10(cov)
    data = {"cov":cov, "cov_log":cov_log, "name":name}
    return data


g_cosmolike = load_cov("../6x2pt_Roman_SO_g_cosmolike.fits", "g cosmolike")
g1 = load_cov("../6x2pt_Roman_SO_gaussian.fits", "g cosmosis EH99 type1 (sum l bins)")
g2 = load_cov("../6x2pt_Roman_SO_gaussian_type2.fits", "g cosmosis EH99 type2 (multiply delta l)")



x = np.arange(len(np.diagonal(g_cosmolike["cov"])))
fig, axs = plt.subplots(2,1, figsize=(9,4))
#axs[0].set_yscale("log")
fig.suptitle("Diag of Covariance, only gaussian part, fsky larger for cmb autocorr")
axs[0].plot(x, np.diagonal(g_cosmolike["cov_log"]), "--",label=g_cosmolike["name"])
axs[0].plot(x, np.diagonal(g1["cov_log"]), "--",label=g1["name"])
axs[0].plot(x, np.diagonal(g2["cov_log"]), "--",label=g2["name"])


axs[0].text(0, -17, "shear ->")
axs[0].axvline(x=1100)
axs[0].text(1100, -11, "galaxy x shear ->")
axs[0].axvline(x=1740)
axs[0].text(1740, -11, "shear x cmb")
axs[0].axvline(x=1940)
axs[0].text(1940, -20, "galaxy")
axs[0].axvline(x=2140)
axs[0].text(2140, -11, "galaxy x cmb")
axs[0].axvline(x=2340)
axs[0].text(2310, -22, "cmb")

axs[0].legend()
axs[0].set_ylabel("Diagonal log Cov")
#axs[0].set_xlabel("Index of covariance")
axs[0].set_xlim(0,2370)

##plt.savefig("v0_4_compare_diag_logspaceRAW.png")
##########################
# plt.figure(figsize=(20,5))
# plt.title("ratio")
axs[1].plot(x,np.diagonal(g1["cov"])/np.diagonal(g_cosmolike["cov"]), "--", label=g1["name"]+" / "+g_cosmolike["name"] );
axs[1].plot(x,np.diagonal(g2["cov"])/np.diagonal(g_cosmolike["cov"]), label=g2["name"]+" / "+g_cosmolike["name"]);


# plt.text(0, 1.2, "shear ->")
# plt.axvline(x=1100)
# plt.text(1100, 1.2, "galaxy x shear ->")
# plt.axvline(x=1740)
# plt.text(1740, 1.2, "shear x cmb")
# plt.axvline(x=1940)
# plt.text(1940, 1.2, "galaxy")
# plt.axvline(x=2140)
# plt.text(2140, 1.2, "galaxy x cmb")
# plt.axvline(x=2340)
# plt.text(2340, 1.2, "cmb")
axs[1].set_ylabel("Diagonal Cov ratio")
axs[1].set_xlabel("Index of covariance")
axs[1].legend(ncol=2)
axs[1].set_xlim(0,2370)

plt.savefig("v0_7_compare_diag_ratioRAW.pdf")