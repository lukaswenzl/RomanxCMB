from __future__ import print_function
from builtins import str
from cosmosis.runtime.config import Inifile
from cosmosis.runtime.pipeline import LikelihoodPipeline
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

output_folder = "modules/RomanxCMB/sample_sigma8_of_z/"

# The easiest way to start a pipeline it from a parameter file.
ini = Inifile("modules/RomanxCMB/params_sigma8ofz.ini")

# You can modify things in the ini file object after loading.
# In this case we will switch off some verbose output
#ini.set("pipeline", "values",  "%(ROMANxCMB_SRC_DIR)s/values.ini")


# Make the pipeline itself
pipeline = LikelihoodPipeline(ini)


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

idx = [params_names.index("rescale_pk_fz--alpha_{}".format(i)) for i in range(5)]

print("the indices for alpha i are")
print(idx)

# for shift in [-0.1, -0.05, 0, 0.05, 0.1]:

#     # In this method of running the pipeline we
#     # pass it a value for each of the parameters 
#     # we have told it to vary.
#     # We could check what these are by looking at
#     #pipeline.varied_params

#     vector = pipeline.start_vector()
#     vector[idx[2]] = vector[idx[2]] + shift

#     data = pipeline.run_parameters(vector)

#     # data is a DataBlock - can get things out of it as in any
#     # cosmosis module:
#     p_k = data["matter_power_lin", "p_k"]
#     z = data["matter_power_lin", "z"]

#     # Make a plot for this value
#     plt.loglog(z, p_k[:,10], label=str(shift))
#     print("Done ", shift)

# # Save our plot.
# plt.xlabel("z")
# plt.ylabel("Pk lin (large scale k)")
# plt.legend()
# plt.savefig(output_folder+"power_spectrum_sigma8ofz_RAW.pdf")



for shift in [-0.1, -0.05, 0, 0.05, 0.1]:

    # In this method of running the pipeline we
    # pass it a value for each of the parameters 
    # we have told it to vary.
    # We could check what these are by looking at
    #pipeline.varied_params

    vector = pipeline.start_vector()
    vector[idx[2]] = vector[idx[2]] + shift

    data = pipeline.run_parameters(vector)

    # data is a DataBlock - can get things out of it as in any
    # cosmosis module:
    p_k = data["matter_power_nl", "p_k"]
    z = data["matter_power_nl", "z"]
    k = data["matter_power_nl", "k_h"]

    # Make a plot for this value
    plt.loglog(z, p_k[:,10], label=str(shift)+", k/h={0:.2f}".format(k[10]))
    plt.loglog(z, p_k[:,190], "--", label=str(shift)+", k/h={0:.2f}".format(k[190]))

    print("Done ", shift)

# Save our plot.
plt.xlabel("z")
plt.ylabel("Pk nl (large scale k)")
plt.legend()
plt.savefig(output_folder+"power_spectrum_nl_sigma8ofz_RAW.pdf")
