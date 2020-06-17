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


g_cosmolike = load_cov("../6x2pt_WFIRST_SO_g_cosmolike_v0_3_1996900.fits", "g cosmolike")
g1 = load_cov("../6x2pt_WFIRST_SO_gaussian_v0_3_a98b754.fits", "g cosmosis type1")
g2 = load_cov("../6x2pt_WFIRST_SO_gaussian_type2test.fits", "g cosmosis type2")
g3 = load_cov("../6x2pt_WFIRST_SO_gaussian_type3test.fits", "g cosmosis type3")








x = np.arange(len(np.diagonal(g_cosmolike["cov"])))
plt.figure(figsize=(20,5))
plt.title("comparison")
plt.plot(x, np.diagonal(g_cosmolike["cov_log"]), "--",label=g_cosmolike["name"])
plt.plot(x, np.diagonal(g1["cov_log"]), "--",label=g1["name"])
plt.plot(x, np.diagonal(g2["cov_log"]), "--",label=g2["name"])
plt.plot(x, np.diagonal(g3["cov_log"]), "--",label=g3["name"])


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

plt.savefig("tmp1.png")
##########################
plt.figure(figsize=(20,5))
plt.title("ratio")
plt.plot(x,np.diagonal(g1["cov"])/np.diagonal(g_cosmolike["cov"]), "--", label=g1["name"]+" / "+g_cosmolike["name"] );
plt.plot(x,np.diagonal(g2["cov"])/np.diagonal(g_cosmolike["cov"]), label=g2["name"]+" / "+g_cosmolike["name"]);
#plt.plot(x,np.diagonal(g2["cov"])/np.diagonal(g1["cov"]), label=g2["name"]+" / "+g1["name"]);

plt.plot(x,np.diagonal(g3["cov"])/np.diagonal(g_cosmolike["cov"]), label=g3["name"]+" / "+g_cosmolike["name"]);


plt.text(0, 1.3, "shear ->")
plt.axvline(x=1100)
plt.text(1100, 1.3, "galaxy x shear ->")
plt.axvline(x=1740)
plt.text(1740, 1.3, "shear x cmb")
plt.axvline(x=1940)
plt.text(1940, 1.3, "galaxy")
plt.axvline(x=2140)
plt.text(2140, 1.3, "galaxy x cmb")
plt.axvline(x=2340)
plt.text(2340, 1.3, "cmb")
plt.ylabel("Diagonal Cov ratio")
plt.xlabel("Index of covariance")
plt.legend()

plt.savefig("tmp2.png")