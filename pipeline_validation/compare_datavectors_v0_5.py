import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd
import astropy.io.fits as fits
#from astropy import units as u
#from astropy import constants as const

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


filename = "../6x2pt_Roman_SO.fits"
two_point_data = twopoint.TwoPointFile.from_fits(filename)


cosmosis = loop_comosis_datavector(two_point_data)
cosmolike = cosmolike_metadata.rearange_cosmolike_datavec(filename, two_point_data.spectra, "../cosmolike_data/cov_indices_apr9.txt")


print("not final")



x = np.arange(len(cosmosis["cl"]))
plt.figure(figsize=(20,5))
plt.title("datavector comparison")
plt.yscale("log")
plt.plot(x, cosmosis["cl"], label="cosmosis")
plt.plot(x, cosmolike["cl"], label="cosmolike")

plt.text(0, 1e-5, "shear ->")
plt.axvline(x=1100)
plt.text(1100, 1e-5, "galaxy x shear ->")
plt.axvline(x=1740)
plt.text(1740, 1e-5, "shear x cmb")
plt.axvline(x=1940)
plt.text(1940, 1e-5, "galaxy")
plt.axvline(x=2140)
plt.text(2140, 1e-5, "galaxy x cmb")
plt.axvline(x=2340)
plt.text(2340, 1e-5, "cmb")

plt.legend()
plt.ylabel("C_l (times l factors)")
plt.xlabel("Index of datavector")


plt.savefig("v0_5_compare_datavector_logspaceRAW.png")

##########################
plt.figure(figsize=(20,5))
plt.title("Ratio of datavectors")

plt.plot(x, cosmosis["cl"]/cosmolike["cl"], label="cosmosis/cosmolike")

plt.text(0, 1.1, "shear ->")
plt.axvline(x=1100)
plt.text(1100, 1.1, "galaxy x shear ->")
plt.axvline(x=1740)
plt.text(1740, 1.1, "shear x cmb")
plt.axvline(x=1940)
plt.text(1940, 1.1, "galaxy")
plt.axvline(x=2140)
plt.text(2140, 1.1, "galaxy x cmb")
plt.axvline(x=2340)
plt.text(2340, 1.1, "cmb")
plt.ylabel("C_l ratio")
plt.xlabel("Index of datavector")
plt.legend()

plt.savefig("v0_5_compare_datavector_ratioRAW.png")