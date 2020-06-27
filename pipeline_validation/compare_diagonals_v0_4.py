import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd
import astropy.io.fits as fits
#from astropy import units as u
#from astropy import constants as const

def load_cov(filename, name):
    hdul = fits.open(filename)
    print(hdul.info())
    cov = hdul[1].data
    cov_log = np.log10(cov)
    data = {"cov":cov, "cov_log":cov_log, "name":name}
    return data


g_cosmolike = load_cov("../6x2pt_WFIRST_SO_g_cosmolike_v0_4_02887fd.fits", "g cosmolike")
g1 = load_cov("../6x2pt_WFIRST_SO_gaussian_v0_4_02887fd.fits", "g cosmosis type1 (sum l bins)")
g2 = load_cov("../6x2pt_WFIRST_SO_gaussian_type2_v0_4_02887fd.fits", "g cosmosis type2 (multiply delta l)")






x = np.arange(len(np.diagonal(g_cosmolike["cov"])))
plt.figure(figsize=(20,5))
plt.title("comparison")
plt.plot(x, np.diagonal(g_cosmolike["cov_log"]), "--",label=g_cosmolike["name"])
plt.plot(x, np.diagonal(g1["cov_log"]), "--",label=g1["name"])
plt.plot(x, np.diagonal(g2["cov_log"]), "--",label=g2["name"])


plt.text(0, -15, "shear ->")
plt.axvline(x=1100)
plt.text(1100, -11, "galaxy x shear ->")
plt.axvline(x=1740)
plt.text(1740, -11, "shear x cmb")
plt.axvline(x=1940)
plt.text(1940, -20, "galaxy")
plt.axvline(x=2140)
plt.text(2140, -11, "galaxy x cmb")
plt.axvline(x=2340)
plt.text(2340, -11, "cmb")

plt.legend()
plt.ylabel("Diagonal log Cov")
plt.xlabel("Index of covariance")

##plt.savefig("v0_4_compare_diag_logspaceRAW.png")
##########################
plt.figure(figsize=(20,5))
plt.title("ratio")
plt.plot(x,np.diagonal(g1["cov"])/np.diagonal(g_cosmolike["cov"]), "--", label=g1["name"]+" / "+g_cosmolike["name"] );
plt.plot(x,np.diagonal(g2["cov"])/np.diagonal(g_cosmolike["cov"]), label=g2["name"]+" / "+g_cosmolike["name"]);


plt.text(0, 1.2, "shear ->")
plt.axvline(x=1100)
plt.text(1100, 1.2, "galaxy x shear ->")
plt.axvline(x=1740)
plt.text(1740, 1.2, "shear x cmb")
plt.axvline(x=1940)
plt.text(1940, 1.2, "galaxy")
plt.axvline(x=2140)
plt.text(2140, 1.2, "galaxy x cmb")
plt.axvline(x=2340)
plt.text(2340, 1.2, "cmb")
plt.ylabel("Diagonal Cov ratio")
plt.xlabel("Index of covariance")
plt.legend()

##plt.savefig("v0_4_compare_diag_ratioRAW.png")