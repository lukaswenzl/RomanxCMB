from __future__ import print_function
from builtins import str
from cosmosis.runtime.config import Inifile
from cosmosis.runtime.pipeline import LikelihoodPipeline
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

output_folder = "modules/RomanxCMB/plots/"

# The easiest way to start a pipeline it from a parameter file.
ini = Inifile("modules/RomanxCMB/params.ini")

# You can modify things in the ini file object after loading.
# In this case we will switch off some verbose output
#ini.set("pipeline", "values",  "%(ROMANxCMB_SRC_DIR)s/values.ini")


# Make the pipeline itself
pipeline = LikelihoodPipeline(ini)

# You can modify which parameters you vary now
pipeline.set_varied("cosmological_parameters", "omega_m", 0.2, 0.4)
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

index = params_names.index("cosmological_parameters--omega_m")
print("the index is")
print(index)
for omega_m in [0.2, 0.25, 0.3, 0.35, 0.4]:

    # In this method of running the pipeline we
    # pass it a value for each of the parameters 
    # we have told it to vary.
    # We could check what these are by looking at
    #pipeline.varied_params

    vector = pipeline.start_vector()
    vector[index] = omega_m

    data = pipeline.run_parameters(vector)

    # data is a DataBlock - can get things out of it as in any
    # cosmosis module:
    ell = data['galaxy_shear_cl', 'ell']
    cl  = data['galaxy_shear_cl', 'bin_1_1']

    # Make a plot for this value
    plt.loglog(ell, np.abs(cl), label=str(omega_m))
    print("Done ", omega_m)

# Save our plot.
plt.legend()
plt.savefig(output_folder+"test_scriptingRAW.pdf")
