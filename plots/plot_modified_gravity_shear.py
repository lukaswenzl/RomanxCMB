from __future__ import print_function
from builtins import str
from cosmosis.runtime.config import Inifile
from cosmosis.runtime.pipeline import LikelihoodPipeline
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt


fig, ax = plt.subplots(2, 3)
ax[0,0].set_title("$\gamma_1 \gamma_1$")
ax[0,0].set_ylabel("$C_l \cdot l(l+1) (2\pi)^{-1}$")
ax[0,0].set_xlabel("l")

ax[0,1].set_title("$\gamma_5 \gamma_5$")
ax[0,2].set_title("$\gamma_{10} \gamma_{10}$")

ax[1,0].set_xscale("log")
ax[1,0].set_ylabel("Cl ratio")
ax[1,0].set_xlabel("l")

ax[1,1].set_xscale("log")
ax[1,2].set_xscale("log")





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

# Let's look through different values of omega_m
# and get a Galaxy Galaxy-Lensing spectrum for each of them


params_names = pipeline.varied_params
# print(params_names)
# params_values = pipeline.start_vector()
# print(params_values)
# params_dict = {}
# for k,v in zip(params_names,params_values):
#     print(k, v)
#     params_dict[str(k)] = v

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

    # data is a DataBlock - can get things out of it as in any
    # cosmosis module:
    ell = data['shear_cl', 'ell']
    cl11  = data['shear_cl', 'bin_1_1'] * ell * (ell+1) / 2. /np.pi
    cl55  = data['shear_cl', 'bin_5_5'] * ell * (ell+1) / 2. /np.pi
    cl1010  = data['shear_cl', 'bin_10_10'] * ell * (ell+1) / 2. /np.pi



    # Make a plot for this value
    if(mu0 == 0.):
        cl11_fiducial = cl11
        cl55_fiducial = cl55
        cl1010_fiducial = cl1010

        ax[0,0].loglog(ell, np.abs(cl11), label="fiducial")
        ax[1,0].plot(ell, np.abs(cl11/cl11_fiducial))
        ax[0,1].loglog(ell, np.abs(cl55), label="fiducial ")
        ax[1,1].plot(ell, np.abs(cl55/cl55_fiducial))
        ax[0,2].loglog(ell, np.abs(cl1010), label="fiducial")
        ax[1,2].plot(ell, np.abs(cl1010/cl1010_fiducial))
    else:
        ax[0,0].loglog(ell, np.abs(cl11), label="$\mu_0$="+str(mu0))
        ax[1,0].plot(ell, np.abs(cl11/cl11_fiducial))
        ax[0,1].loglog(ell, np.abs(cl55), label="mu0="+str(mu0))
        ax[1,1].plot(ell, np.abs(cl55/cl55_fiducial))
        ax[0,2].loglog(ell, np.abs(cl1010), label="mu0="+str(mu0))
        ax[1,2].plot(ell, np.abs(cl1010/cl1010_fiducial))



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

    # data is a DataBlock - can get things out of it as in any
    # cosmosis module:
    ell = data['shear_cl', 'ell']
    cl11  = data['shear_cl', 'bin_1_1'] * ell * (ell+1) / 2. /np.pi
    cl55  = data['shear_cl', 'bin_5_5'] * ell * (ell+1) / 2. /np.pi
    cl1010  = data['shear_cl', 'bin_10_10'] * ell * (ell+1) / 2. /np.pi

    # Make a plot for this value
    ax[0,0].loglog(ell, np.abs(cl11),"--", label="$\Sigma_0$="+str(sigma0))
    ax[1,0].plot(ell, np.abs(cl11/cl11_fiducial),"--")
    ax[0,1].loglog(ell, np.abs(cl55),"--", label="$\Sigma_0$="+str(sigma0))
    ax[1,1].plot(ell, np.abs(cl55/cl55_fiducial),"--")
    ax[0,2].loglog(ell, np.abs(cl1010),"--", label="$\Sigma_0$="+str(sigma0))
    ax[1,2].plot(ell, np.abs(cl1010/cl1010_fiducial),"--")

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

    # data is a DataBlock - can get things out of it as in any
    # cosmosis module:
    ell = data['shear_cl', 'ell']
    cl11  = data['shear_cl', 'bin_1_1'] * ell * (ell+1) / 2. /np.pi
    cl55  = data['shear_cl', 'bin_5_5'] * ell * (ell+1) / 2. /np.pi
    cl1010  = data['shear_cl', 'bin_10_10'] * ell * (ell+1) / 2. /np.pi

    # Make a plot for this value
    ax[0,0].loglog(ell, np.abs(cl11),":", label="fR="+str(fR))
    ax[1,0].plot(ell, np.abs(cl11/cl11_fiducial),":")
    ax[0,1].loglog(ell, np.abs(cl55),":", label="fR="+str(fR))
    ax[1,1].plot(ell, np.abs(cl55/cl55_fiducial),":")
    ax[0,2].loglog(ell, np.abs(cl1010),":", label="fR="+str(fR))
    ax[1,2].plot(ell, np.abs(cl1010/cl1010_fiducial),":")

    print("Done ", fR)



# Save our plot.
ax[0,0].legend()
plt.savefig(output_folder+"Clgammagamma_modgrav_dependenceRAW.pdf")
