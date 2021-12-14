from __future__ import print_function
from builtins import str
from ctypes import c_ushort
from cosmosis.runtime.config import Inifile
from cosmosis.runtime.pipeline import LikelihoodPipeline
import numpy as np
import matplotlib
matplotlib.use('agg')

import astropy.io.fits as fits
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

import matplotlib.pyplot as plt

plt.style.use(['modules/RomanxCMB/plots/paper.mplstyle'])

import matplotlib as mpl
#mpl.rcParams['ytick.labelsize'] = 7
#mpl.rcParams['xtick.labelsize'] = 7

all_data = []
all_data_label = []


from cycler import cycler
cmap_wBLACK = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
plt.rc('axes', prop_cycle=(cycler('color', cmap_wBLACK) ))

fig, ax = plt.subplots(2, 6, sharex="col", sharey="row", figsize=(8.5,4), gridspec_kw={'height_ratios': [1, 4]})

ax[0,0].set_title("$\delta_{\textrm{g 5}}  \delta_{\textrm{g 5}}$")
ax[0,0].set_ylabel("$|C_l|$")
#ax[0,0].set_xlabel("l")

ax[0,1].set_title("$\delta_{\textrm{ g 3}} \gamma_7$")
ax[0,2].set_title("$\gamma_5 \gamma_5$")
ax[0,3].set_title("$\delta_{\textrm{ g 5}} \kappa$")
ax[0,4].set_title("$\gamma_5 \kappa$")
ax[0,5].set_title("$\kappa \kappa$")


ax[1,0].set_ylabel("$C_\ell$ ratio")
ax[1,0].set_xlabel("Multipole $\ell$")
ax[1,1].set_xlabel("Multipole $\ell$")
ax[1,2].set_xlabel("Multipole $\ell$")
ax[1,3].set_xlabel("Multipole $\ell$")
ax[1,4].set_xlabel("Multipole $\ell$")
ax[1,5].set_xlabel("Multipole $\ell$")


ax[1,0].set_xscale("log")
ax[1,1].set_xscale("log")
ax[1,2].set_xscale("log")
ax[1,3].set_xscale("log")
ax[1,4].set_xscale("log")
ax[1,5].set_xscale("log")


#share single yaxis
#ax[1,0].set_ylim([0.925, 1.075])
ax[1,0].set_ylim([0.88, 1.12])

#ax[1,0].set_yticks([0.94,0.96, 0.98, 1.0, 1.02, 1.04, 1.06])
#ax[1,0].set_yticklabels(["-6%","-4%", "-2%", "  ", "2%", "4%", "6%"])

ax[1,0].set_yticks([0.9,0.95, 1.0, 1.05, 1.10])
ax[1,0].set_yticklabels(["-10%","-5%", "0%", "5%", "10%"])



ax[1,0].set_xlim([30, 3000])
ax[1,1].set_xlim([30, 3000])
ax[1,2].set_xlim([30, 3000])
ax[1,3].set_xlim([30, 3000])
ax[1,4].set_xlim([30, 3000])
ax[1,5].set_xlim([30, 3000])


# ax[1,0].get_shared_y_axes().join(ax[1,0], ax[1,1])
# ax[1,0].get_shared_y_axes().join(ax[1,0], ax[1,2])
# ax[1,0].get_shared_y_axes().join(ax[1,0], ax[1,3])
# ax[1,0].get_shared_y_axes().join(ax[1,0], ax[1,4])
# ax[1,0].get_shared_y_axes().join(ax[1,0], ax[1,5])

# ax[1,1].set_yticks([0.94,0.96, 0.98, 1.0, 1.02, 1.04, 1.06])
# ax[1,2].set_yticks([0.94,0.96, 0.98, 1.0, 1.02, 1.04, 1.06])
# ax[1,3].set_yticks([0.94,0.96, 0.98, 1.0, 1.02, 1.04, 1.06])
# ax[1,4].set_yticks([0.94,0.96, 0.98, 1.0, 1.02, 1.04, 1.06])
# ax[1,5].set_yticks([0.94,0.96, 0.98, 1.0, 1.02, 1.04, 1.06])
# ax[1,1].set_yticklabels([])
# ax[1,2].set_yticklabels([])
# ax[1,3].set_yticklabels([])
# ax[1,4].set_yticklabels([])
# ax[1,5].set_yticklabels([])



###########show errors for Cls
#loading covariance from file
filename = "modules/RomanxCMB/6x2pt_Roman_SO_v1_2_bf26108.fits"
hdul = fits.open(filename)
print(hdul.info())
cov = hdul["COVMAT"].data
cov_hdr = hdul["COVMAT"].header

starts = {}
for i in range(6):
    starts[cov_hdr["NAME_"+str(i)]] = cov_hdr["STRT_"+str(i)]

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
def get_scales( x_min, x_max, nbins, logspaced=True, integer_lims=False, two_thirds_midpoint=False):
    """
    Get scales
    """
    if logspaced:
        log_lims = np.linspace( np.log(x_min), np.log(x_max), nbins+1 )
        lims = np.exp(log_lims)
        if two_thirds_midpoint:
            xmin = lims[:-1]
            xmax = lims[1:]
            mids = (2./3.) * (xmax**3 - xmin**3) / (xmax**2 - xmin**2)
        else:
            xmin = log_lims[:-1]
            xmax = log_lims[1:]
            log_mids = 0.5*(xmin + xmax)
            mids = np.exp(log_mids)
    else:
        lims = np.linspace(x_min, x_max, nbins+1)
        xmin = lims[:-1]
        xmax = lims[1:]
        if two_thirds_midpoint:
            mids = (2./3.) * (xmax**3 - xmin**3) / (xmax**2 - xmin**2)
        else:
            mids = 0.5 * (xmin + xmax)

    return lims, mids

def get_ell_scales( ell_min, ell_max, nbins, logspaced=True, two_thirds_midpoint=False, ell_as_int=True):
    #ells should be integers. So get scales using get_scales and then 
    #convert
    lims, mids = get_scales( ell_min, ell_max, nbins, logspaced=logspaced, two_thirds_midpoint=two_thirds_midpoint)
    if (ell_as_int):
        ell_lims = (np.floor(lims)).astype(int)
        ell_mids = np.exp( 0.5 * ( np.log(ell_lims[1:]) + np.log(ell_lims[:-1]) ) )
    else:
        ell_lims = lims 
        ell_mids = mids
    return ell_lims, ell_mids


bins, _ = get_ell_scales( 30, 3000, 20,
            logspaced=True, two_thirds_midpoint=False, ell_as_int=True )


for (name,idx, i,j) in [("galaxy_cl", 0, 5,5), ("galaxy_shear_cl", 1, 3,7), ("shear_cl", 2, 5,5), ("galaxy_cmbkappa_cl", 3, 5,1), ("shear_cmbkappa_cl", 4, 5,1), ("cmbkappa_cl", 5, 1,1)]:
    #name = "galaxy_cl"
    hdr = hdul[name].header
    length = hdr["NAXIS2"]
    variances = np.diagonal(cov)[starts[name]:starts[name]+length]
    Cls = hdul[name].data

    #print("approximate ranges (NOT EXACT)")
    #bins = np.logspace(np.log10(30),np.log10(3000),21)
    ranges = bins[1:]-bins[:-1]
    mask = np.logical_and(Cls["BIN1"] == i, Cls["BIN2"]==j)
    Cl = Cls[mask]
    variance = variances[mask]
    

    make_error_boxes(ax[0,idx], Cl["ANG"], Cl["VALUE"], np.array([ranges/2, ranges/2]), np.array([np.sqrt(variance),np.sqrt(variance)]), facecolor="grey")

    make_error_boxes(ax[1,idx], Cl["ANG"], Cl["VALUE"]/Cl["VALUE"], np.array([ranges/2, ranges/2]), np.array([np.sqrt(variance)/np.abs(Cl["VALUE"]),np.sqrt(variance)/np.abs(Cl["VALUE"])]), facecolor="grey", alpha=0.2)
    #plt.plot(Cl["ANG"], Cl["VALUE"], label="$C_l^{\gamma"+str(i)+" \gamma "+str(j)+"}$", color = cmap[idx])


    # plt.plot(Cl["ANG"], np.sqrt(ranges*variance), label="Unbinned error")
    #plt.plot(Cl["ANG"], np.sqrt(variance),".", label="Binned error", alpha=0.4, color = cmap[idx])
    #idx = idx +1

############
######

output_folder = "modules/RomanxCMB/plots/modified_gravity/"

# The easiest way to start a pipeline it from a parameter file.
ini = Inifile("modules/RomanxCMB/params_modgrav.ini")

# You can modify things in the ini file object after loading.
# In this case we will switch off some verbose output
#
ini.set("pipeline", "values",  "%(ROMANxCMB_SRC_DIR)s/plots/modified_gravity/values_modgrav_musigma.ini")
ini.set("DEFAULT", "2PT_FILE",  "%(ROMANxCMB_SRC_DIR)s/6x2pt_Roman_SO_v1_2_bf26108.fits")
#6x2pt_Roman_SO_fR6_v1_2_809a767.fits


# Make the pipeline itself
pipeline = LikelihoodPipeline(ini)

# pipeline.set_varied("modified_gravity_parameters", "mu0", -3., 3.)
# pipeline.set_varied("modified_gravity_parameters", "sigma0", -3., 3.)

# pipeline.set_varied("modified_gravity_parameters", "model", 0, 3)
#pipeline.set_varied("modified_gravity_parameters", "f_of_R_fR", 0., 0.001)



#pipeline.set_fixed("cosmological_parameters", "h0", 0.72)

# You can also override these properties if useful
pipeline.quiet = True
pipeline.debug = False
pipeline.timing = False




params_names = pipeline.varied_params
index = params_names.index("modified_gravity_parameters--mu0")
print("the index is")
print(index)
for mu0 in [ 0.,-0.2, 0.2]:

    # In this method of running the pipeline we
    # pass it a value for each of the parameters 
    # we have told it to vary.
    # We could check what these are by looking at
    #pipeline.varied_params

    vector = pipeline.start_vector()
    vector[index] = mu0

    data = pipeline.run_parameters(vector)
    loc_dict = {}
    for key in data.keys():
        loc_dict[key] = data[key]
    all_data.append(loc_dict)

    # data is a DataBlock - can get things out of it as in any
    # cosmosis module:
    ell = data['galaxy_cl', 'ell']
    cl1  = data['galaxy_cl', 'bin_5_5'] #* ell * (ell+1) / 2. /np.pi
    cl2  = data['galaxy_shear_cl', 'bin_3_7'] #* ell * (ell+1) / 2. /np.pi
    cl3  = data['shear_cl', 'bin_5_5'] #* ell * (ell+1) / 2. /np.pi
    cl4  = data['galaxy_cmbkappa_cl', 'bin_5_1'] #* ell * (ell+1) / 2. /np.pi
    cl5  = data['shear_cmbkappa_cl', 'bin_5_1'] #* ell * (ell+1) / 2. /np.pi
    cl6  = data['cmbkappa_cl', 'bin_1_1'] #* ell * (ell+1) / 2. /np.pi

    # Make a plot for this value
    if(mu0 == 0.):
        all_data_label.append ("fiducial")
        cl1_fiducial = cl1
        cl2_fiducial = cl2
        cl3_fiducial = cl3
        cl4_fiducial = cl4
        cl5_fiducial = cl5
        cl6_fiducial = cl6


        ax[0,0].loglog(ell, np.abs(cl1), label="fiducial")
        ax[1,0].plot(ell, np.abs(cl1/cl1_fiducial), label="fiducial")
        ax[0,1].loglog(ell, np.abs(cl2), label="fiducial ")
        ax[1,1].plot(ell, np.abs(cl2/cl2_fiducial), label="fiducial")
        ax[0,2].loglog(ell, np.abs(cl3), label="fiducial")
        ax[1,2].plot(ell, np.abs(cl3/cl3_fiducial))
        ax[0,3].loglog(ell, np.abs(cl4), label="fiducial")
        ax[1,3].plot(ell, np.abs(cl4/cl4_fiducial))
        ax[0,4].loglog(ell, np.abs(cl5), label="fiducial")
        ax[1,4].plot(ell, np.abs(cl5/cl5_fiducial))
        ax[0,5].loglog(ell, np.abs(cl6), label="fiducial")
        ax[1,5].plot(ell, np.abs(cl6/cl6_fiducial))
    else:
        all_data_label.append("$\mu_0$="+str(mu0))
        ax[0,0].loglog(ell, np.abs(cl1), label="$\mu_0$="+str(mu0))
        ax[1,0].plot(ell, np.abs(cl1/cl1_fiducial), label="$\mu_0$="+str(mu0))
        ax[0,1].loglog(ell, np.abs(cl2), label="$\mu_0$="+str(mu0))
        ax[1,1].plot(ell, np.abs(cl2/cl2_fiducial), label="$\mu_0$="+str(mu0))
        ax[0,2].loglog(ell, np.abs(cl3), label="$\mu_0$="+str(mu0))
        ax[1,2].plot(ell, np.abs(cl3/cl3_fiducial))
        ax[0,3].loglog(ell, np.abs(cl4), label="$\mu_0$="+str(mu0))
        ax[1,3].plot(ell, np.abs(cl4/cl4_fiducial))
        ax[0,4].loglog(ell, np.abs(cl5), label="$\mu_0$="+str(mu0))
        ax[1,4].plot(ell, np.abs(cl5/cl5_fiducial))
        ax[0,5].loglog(ell, np.abs(cl6), label="$\mu_0$="+str(mu0))
        ax[1,5].plot(ell, np.abs(cl6/cl6_fiducial))


    print("Done ", mu0)

index = params_names.index("modified_gravity_parameters--sigma0")
print("the index is")
print(index)
for sigma0 in [ -0.05, 0.05]: #approximately the 1 sigma constraint

    # In this method of running the pipeline we
    # pass it a value for each of the parameters 
    # we have told it to vary.
    # We could check what these are by looking at
    #pipeline.varied_params

    vector = pipeline.start_vector()
    vector[index] = sigma0

    data = pipeline.run_parameters(vector)
    loc_dict = {}
    for key in data.keys():
        loc_dict[key] = data[key]
    all_data.append(loc_dict)
    all_data_label.append("$\Sigma_0$="+str(sigma0))

    # data is a DataBlock - can get things out of it as in any
    # cosmosis module:
    ell = data['galaxy_cl', 'ell']
    cl1  = data['galaxy_cl', 'bin_5_5'] #* ell * (ell+1) / 2. /np.pi
    cl2  = data['galaxy_shear_cl', 'bin_3_7'] #* ell * (ell+1) / 2. /np.pi
    cl3  = data['shear_cl', 'bin_5_5'] #* ell * (ell+1) / 2. /np.pi
    cl4  = data['galaxy_cmbkappa_cl', 'bin_5_1'] #* ell * (ell+1) / 2. /np.pi
    cl5  = data['shear_cmbkappa_cl', 'bin_5_1'] #* ell * (ell+1) / 2. /np.pi
    cl6  = data['cmbkappa_cl', 'bin_1_1'] #* ell * (ell+1) / 2. /np.pi

    # Make a plot for this value
    ax[0,0].loglog(ell, np.abs(cl1), "--", label="$\Sigma_0$="+str(sigma0))
    ax[1,0].plot(ell, np.abs(cl1/cl1_fiducial), "--", label="$\Sigma_0$="+str(sigma0))
    ax[0,1].loglog(ell, np.abs(cl2), "--", label="$\Sigma_0$="+str(sigma0))
    ax[1,1].plot(ell, np.abs(cl2/cl2_fiducial), "--", label="$\Sigma_0$="+str(sigma0))
    ax[0,2].loglog(ell, np.abs(cl3), "--", label="$\Sigma_0$="+str(sigma0))
    ax[1,2].plot(ell, np.abs(cl3/cl3_fiducial), "--")
    ax[0,3].loglog(ell, np.abs(cl4), "--", label="$\Sigma_0$="+str(sigma0))
    ax[1,3].plot(ell, np.abs(cl4/cl4_fiducial), "--")
    ax[0,4].loglog(ell, np.abs(cl5), "--", label="$\Sigma_0$="+str(sigma0))
    ax[1,4].plot(ell, np.abs(cl5/cl5_fiducial), "--")
    ax[0,5].loglog(ell, np.abs(cl6), "--", label="$\Sigma_0$="+str(sigma0))
    ax[1,5].plot(ell, np.abs(cl6/cl6_fiducial), "--")

    print("Done ", sigma0)


######## f(R)



ini = Inifile("modules/RomanxCMB/params_modgrav.ini")

# You can modify things in the ini file object after loading.
# In this case we will switch off some verbose output
#
ini.set("pipeline", "values",  "%(ROMANxCMB_SRC_DIR)s/plots/modified_gravity/values_modgrav_fR.ini")
ini.set("DEFAULT", "2PT_FILE",  "%(ROMANxCMB_SRC_DIR)s/6x2pt_Roman_SO_v1_2_bf26108.fits")
ini.set("fR_sampling", "sample_logspace",  "F")


#6x2pt_Roman_SO_fR6_v1_2_809a767.fits

# Make the pipeline itself
pipeline = LikelihoodPipeline(ini)
pipeline.quiet = True
pipeline.debug = False
pipeline.timing = False

params_names = pipeline.varied_params
index = params_names.index("modified_gravity_parameters--f_of_r_fr")
print("the index is")
print(index)
for fR in [ 1.e-6]: #approximately the 1 sigma constraint

    # In this method of running the pipeline we
    # pass it a value for each of the parameters 
    # we have told it to vary.
    # We could check what these are by looking at
    #pipeline.varied_params

    vector = pipeline.start_vector()
    vector[index] = fR

    data = pipeline.run_parameters(vector)
    loc_dict = {}
    for key in data.keys():
        loc_dict[key] = data[key]
    all_data.append(loc_dict)

    # data is a DataBlock - can get things out of it as in any
    # cosmosis module:
    ell = data['galaxy_cl', 'ell']
    cl1  = data['galaxy_cl', 'bin_5_5'] #* ell * (ell+1) / 2. /np.pi
    cl2  = data['galaxy_shear_cl', 'bin_3_7'] #* ell * (ell+1) / 2. /np.pi
    cl3  = data['shear_cl', 'bin_5_5'] #* ell * (ell+1) / 2. /np.pi
    cl4  = data['galaxy_cmbkappa_cl', 'bin_5_1'] #* ell * (ell+1) / 2. /np.pi
    cl5  = data['shear_cmbkappa_cl', 'bin_5_1'] #* ell * (ell+1) / 2. /np.pi
    cl6  = data['cmbkappa_cl', 'bin_1_1'] #* ell * (ell+1) / 2. /np.pi

    # Make a plot for this value
    if(fR == 1e-6):
        fR_name = 6 
        all_data_label.append("fR6")
    else:
        print("error! need to define name")
    ax[0,0].loglog(ell, np.abs(cl1), ":", label="fR"+str(fR_name))
    ax[1,0].plot(ell, np.abs(cl1/cl1_fiducial), ":", label="fR"+str(fR_name))
    ax[0,1].loglog(ell, np.abs(cl2), ":", label="fR"+str(fR_name))
    ax[1,1].plot(ell, np.abs(cl2/cl2_fiducial), ":", label="fR"+str(fR_name))
    ax[0,2].loglog(ell, np.abs(cl3), ":", label="fR"+str(fR_name))
    ax[1,2].plot(ell, np.abs(cl3/cl3_fiducial), ":")
    ax[0,3].loglog(ell, np.abs(cl4), ":", label="fR"+str(fR_name))
    ax[1,3].plot(ell, np.abs(cl4/cl4_fiducial), ":")
    ax[0,4].loglog(ell, np.abs(cl5), ":", label="fR"+str(fR_name))
    ax[1,4].plot(ell, np.abs(cl5/cl5_fiducial), ":")
    ax[0,5].loglog(ell, np.abs(cl6), ":", label="fR"+str(fR_name))
    ax[1,5].plot(ell, np.abs(cl6/cl6_fiducial), ":")

    print("Done ", fR)


#scale cuts
#g5g5 cutoff is 1702.1
ax[1,0].axvline(1702.1, color="black", linestyle="--")
#ax[1,0].fill_between([1702.1, 3000], [2, 2], [0.5, 0.5], facecolor=None, hatch='//', alpha=0.99)

#g3gamma7 cutoff is 1253.5
ax[1,1].axvline(1253.5, color="black", linestyle="--")
#ax[1,1].axvspan(1253.5, 3000, facecolor=None, hatch="X", linewidth=2, edgecolor="black", alpha = 0.5)
#g5kappa cutoff is 1702.1
ax[1,3].axvline(1702.1, color="black", linestyle="--")
#ax[1,3].axvspan(1702.1, 3000, facecolor=None, hatch="X", edgecolor="black", alpha = 0.01)


# Save our plot.
#ax[1,0].legend(fontsize="xx-small")
ax[1,1].legend(fontsize="xx-small")

#ax[1,5].set_yticks([0.95, 1.0, 1.05])
#plt.tight_layout()
plt.subplots_adjust(hspace=.0, wspace=0.)
plt.savefig(output_folder+"Clall_modgrav_dependenceRAW.pdf")


import pickle as pkl
pkl.dump(all_data,  open( output_folder+"modified_gravity_plot_data.pkl", "wb" ) )
pkl.dump(all_data_label, open( output_folder+"modified_gravity_plot_data_labels.pkl", "wb"))
np.savetxt(output_folder+"modified_gravity_plot_data_labels.npy", all_data_label,  fmt='%s')
print(all_data_label)