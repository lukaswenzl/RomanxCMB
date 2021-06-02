"""
Apply a cutoff to the sum of w0 and wa. 

"""
import numpy as np
from cosmosis.datablock import option_section, names
from cosmosis.runtime.utils import Timer

cosmo = names.cosmological_parameters

def setup(options):
    w_infinite_cutoff = options.get_double(option_section, "w_infinite_cutoff", -0.33334)
    #new_name = options.get_double(option_section, "new_name", "matter_power_nl")

    return {"w_infinite_cutoff":w_infinite_cutoff}




def execute(block, config):

    w = block[cosmo, 'w']
    wpluswa = block[cosmo, 'wpluswa']

    if(wpluswa > config["w_infinite_cutoff"]):
        #if outside of range let pipeline run with fidcucial params but set flag to indicate -infity likelihood
        print("Error: Outside of supported range for w+wa")

    block[cosmo, 'set_likelihood_minus_infinity'] = False 
    block[cosmo, 'wa'] = wpluswa - w
    print("Hello Lukas")
    print(block[cosmo, 'wa'])

    return 0