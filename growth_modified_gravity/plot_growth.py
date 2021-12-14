from __future__ import print_function
from builtins import str
from cosmosis.runtime.config import Inifile
from cosmosis.runtime.pipeline import LikelihoodPipeline
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

output_folder = "modules/RomanxCMB/growth_modified_gravity/"

# The easiest way to start a pipeline it from a parameter file.
ini = Inifile("modules/RomanxCMB/params_modgrav.ini")

# You can modify things in the ini file object after loading.
# In this case we will switch off some verbose output
#ini.set("pipeline", "values",  "%(ROMANxCMB_SRC_DIR)s/values.ini")


# Make the pipeline itself
pipeline = LikelihoodPipeline(ini)

pipeline.set_varied("modified_gravity_parameters", "mu0", -3., 3.)
pipeline.set_varied("modified_gravity_parameters", "sigma0", -3., 3.)


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
for mu0 in [-1, 0., 1]:

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
    ell = data['galaxy_shear_cl', 'ell']
    cl  = data['galaxy_shear_cl', 'bin_1_1']

    # Make a plot for this value
    plt.loglog(ell, np.abs(cl), label="mu0="+str(mu0))
    print("Done ", mu0)

index = params_names.index("modified_gravity_parameters--sigma0")
print("the index is")
print(index)
for sigma0 in [ -1.,0., 1.]:

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
    ell = data['galaxy_shear_cl', 'ell']
    cl  = data['galaxy_shear_cl', 'bin_1_1']

    # Make a plot for this value
    plt.loglog(ell, np.abs(cl),"--", label="sigma0="+str(sigma0))
    print("Done ", sigma0)

# Save our plot.
plt.legend()
plt.savefig(output_folder+"Clg1g1_modgrav_dependenceRAW.pdf")
