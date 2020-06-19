#plot the datavector cl's and errors
#assumes l 30 to 3000 with 20 log spaced bins
import numpy as np

import matplotlib
matplotlib.use('Agg')


import matplotlib.pyplot as plt
import astropy.io.fits as fits
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


#cosmosis covariance, assuming fsky 0.05 throughout
filename = "6x2pt_Roman_SO_gaussian.fits"
outputfolder = ""

def make_error_boxes(ax, xdata, ydata, xerror, yerror, facecolor='r',
                     edgecolor='None', alpha=0.5, label="Measurements"):

    # Create list for all the error patches
    errorboxes = []

    # Loop over data points; create box from errors at each point
    for x, y, xe, ye in zip(xdata, ydata, xerror.T, yerror.T):
        rect = Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum())
        errorboxes.append(rect)

    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha,
                         edgecolor=edgecolor, label=label)

    # Add collection to axes
    ax.add_collection(pc)

    # Plot errorbars
    #artists = ax.errorbar(xdata, ydata, xerr=xerror, yerr=yerror,
    #                      fmt='None', ecolor='k')

    #return artists

fig, ax = plt.subplots(1)

hdul = fits.open(filename)
print(hdul.info())
cov = hdul["COVMAT"].data
cov_hdr = hdul["COVMAT"].header

starts = {}
for i in range(6):
    starts[cov_hdr["NAME_"+str(i)]] = cov_hdr["STRT_"+str(i)]


#########################################################################################################################
##cmb lensing
#########################################################################################################################


name = "cmbkappa_cl"
hdr = hdul[name].header
length = hdr["NAXIS2"]
variance = np.diagonal(cov)[starts[name]:starts[name]+length]
Cl = hdul[name].data

print("approximate ranges (NOT EXACT)")
bins = np.logspace(np.log10(30),np.log10(3000),21)
ranges = bins[1:]-bins[:-1]


#plt.errorbar(Cl["ANG"], Cl["VALUE"], yerr=np.sqrt(variance), label="measurements")
make_error_boxes(ax, Cl["ANG"], Cl["VALUE"], np.array([ranges/2, ranges/2]), np.array([np.sqrt(variance),np.sqrt(variance)]), facecolor="blue")
plt.plot(Cl["ANG"], Cl["VALUE"], label=r"fiducial $C_l^{\kappa \kappa}$", color="orange")


plt.plot(Cl["ANG"], np.sqrt(ranges*variance), label="Unbinned error")
plt.plot(Cl["ANG"], np.sqrt(variance),".", label="Binned error", color="blue", alpha=0.4)

#load external noise curve for SO
cov_SO = np.loadtxt("SO_noise_curves/lensing_v3_0_0/Apr17_mv_nlkk_deproj0_SENS1_fsky_16000_iterOn.csv")
plt.plot(Cl["ANG"], np.interp(Cl["ANG"], cov_SO[:,0], cov_SO[:,1]), "--", color="grey", label="SO noise curve $f_{sky}=0.4$")

#going from Cl error back to original ( With this you can reconstruct the noise curve from the error, was only used to check calculation)
#fsky = 0.4
#noisecurve = np.sqrt(variance * fsky* (2*Cl["ANG"]+1)*ranges/2) - Cl["VALUE"]
#plt.plot(Cl["ANG"], noisecurve, label="TEST")
#plt.errorbar(Cl["ANG"],Cl["VALUE"], yerr=noisecurve/ranges, label="TESTnoise")



plt.ylabel(r"$C_l^{\kappa \kappa}$")
plt.xlabel("l")
plt.yscale("log")
plt.xscale("log")
plt.legend()

plt.savefig(outputfolder+"Cl_cmblensing_with_errorRAW.png")

#calculate SNR
cov_cutout = cov[starts[name]:starts[name]+length, starts[name]:starts[name]+length]
inv = np.linalg.inv(cov_cutout) 
SNR = np.einsum('i,ij,j', Cl["VALUE"], inv, Cl["VALUE"])
print("SNR for "+name)
print(SNR)



#########################################################################################################################
##weak lensing shear
#########################################################################################################################

fig, ax = plt.subplots(1)
name = "shear_cl"
hdr = hdul[name].header
length = hdr["NAXIS2"]
variances = np.diagonal(cov)[starts[name]:starts[name]+length]
Cls = hdul[name].data

print("approximate ranges (NOT EXACT)")
bins = np.logspace(np.log10(30),np.log10(3000),21)
ranges = bins[1:]-bins[:-1]

#select bin combinations:
cmap =  [ "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
idx = 0
for (i,j) in [(1,1), (4,4),(4,10), (10,10)]:
    mask = np.logical_and(Cls["BIN1"] == i, Cls["BIN2"]==j)
    Cl = Cls[mask]
    variance = variances[mask]
    

    make_error_boxes(ax, Cl["ANG"], Cl["VALUE"], np.array([ranges/2, ranges/2]), np.array([np.sqrt(variance),np.sqrt(variance)]), facecolor="blue")
    plt.plot(Cl["ANG"], Cl["VALUE"], label="$C_l^{\gamma"+str(i)+" \gamma "+str(j)+"}$", color = cmap[idx])


    # plt.plot(Cl["ANG"], np.sqrt(ranges*variance), label="Unbinned error")
    plt.plot(Cl["ANG"], np.sqrt(variance),".", label="Binned error", alpha=0.4, color = cmap[idx])
    idx = idx +1

plt.ylabel(r"$C_l^{\gamma \gamma}$")
plt.xlabel("l")
plt.yscale("log")
plt.xscale("log")
plt.legend()

plt.savefig(outputfolder+"Cl_shear_with_errorRAW.png")

#calculate SNR
cov_cutout = cov[starts[name]:starts[name]+length, starts[name]:starts[name]+length]
inv = np.linalg.inv(cov_cutout) 
SNR = np.einsum('i,ij,j', Cls["VALUE"], inv, Cls["VALUE"])
print("SNR for "+name)
print(SNR)