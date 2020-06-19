# Import lots of bits from across cosmosis.
from cosmosis.samplers.fisher.fisher_sampler import FisherSampler
from cosmosis.runtime.pipeline import LikelihoodPipeline
from cosmosis.runtime.config import Inifile
from cosmosis.output.in_memory_output import InMemoryOutput
from cosmosis.output.text_output import TextColumnOutput
import numpy as np

from cosmosis.postprocessing.postprocess import FisherProcessor

#plots
from cosmosis.postprocessing import plots
from cosmosis.postprocessing import lazy_pylab as pylab


output_folder = "TESTscript_gaussian_vary_tracer_combinations/"
import os
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

#written for v0.1
#run with mpirun -n 14 python modules/WFIRSTxCMB/script_gaussian_vary_tracer_combinations.py
from mpi4py import MPI #without this line mpi won't work but if mpi is not installed there will be an error: can also run without it


runs = []
#--------------------------------------------------------------------
# #RUN 1

output = "6x2pt.txt"
# Create the configuration based on a file.
override = {}
override[("output", "filename")] = output_folder + "6x2pt_WFIRST_SO_gaussian-chain.txt"
#override[("pipeline", "values")] = "modules/WFIRSTxCMB/values_few_variables.ini" 
#override[("pipeline", "priors")] = ""


ini = Inifile("modules/WFIRSTxCMB/params_gaussian_covariance.ini", override = override)

# Build the likelihood pipeline based on the config
pipeline = LikelihoodPipeline(ini)

#set_parameter_file(pipeline, cosmos="LCDM", m="vary", b="fixed")#todo

# change parameters to fixed to make work
# pipeline.set_fixed("shear_calibration_parameters", "m1", 0.0)
# pipeline.set_fixed("shear_calibration_parameters", "m2", 0.0)
# pipeline.set_fixed("shear_calibration_parameters", "m3", 0.0)
# pipeline.set_varied("cosmological_parameters", "omega_m", 0.2, 0.4)
#pipeline.set_fixed("cosmological_parameters", "omega_m", 0.3)

# tell cosmosis to keep the output in memory instead
# # of writing to a text file
output = InMemoryOutput()
#output = TextColumnOutput(output_folder + output)

# Create the sampler.  In this case for speed
# we will just use a Fisher Matrix
sampler = FisherSampler(ini, pipeline, output)#ini)
sampler.config()

while not sampler.is_converged(): 
    sampler.execute()
# import pdb; pdb.set_trace()
# The output object now contains the data
# that would have been output to file
#fisher_matrix = np.array(output.rows)

# You could now do stuff with this matrix if you wanted.
#print(fisher_matrix)

runs.append({"label":"test", "output":output, "ini":ini})


#--------------------------------------------------------------------
# #RUN 2
# output = "gg.txt"

# override = {}
# override[("output", "filename")] = output_folder + "6x2pt_WFIRST_SO_gaussian-chain.txt"
# override[("pipeline", "values")] = "modules/WFIRSTxCMB/values_few_variables.ini" 
# override[("pipeline", "priors")] = ""
# ini2 = Inifile("modules/WFIRSTxCMB/params_gaussian_covariance.ini", override={("DEFAULT", "2pt_data_sets"):"test"})
# pipeline = LikelihoodPipeline(ini)
# output = TextColumnOutput(output_folder + output)
# sampler = FisherSampler(ini, pipeline, output)
# sampler.config()
# while not sampler.is_converged(): 
#     sampler.execute()


import os
os.system('ls -l')


outputs = {}
for i in range(len(runs)):#??
    #plots
    #processor = FisherProcessor(runs[i]["ini"], runs[i]["label"], i)#, **vars(args))
    processor = FisherProcessor(runs[i]["output"], runs[i]["label"], i)#, **vars(args))
    import pdb; pdb.set_trace()
    #Inherit any plots from the previous postprocessor
    #so we can make plots with multiple datasets on
    processor.outputs.update(outputs)
    processor.run()

# if args.tweaks:
# 		tweaks = Tweaks.instances_from_file(args.tweaks)
# 		for tweak in tweaks:
# 			processor.apply_tweaks(tweak)
#Save all the image files and close the text files
processor.finalize()