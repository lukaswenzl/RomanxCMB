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
    wa = block[cosmo, 'wa']

    # if(w + wa > config["w_infinite_cutoff"]):
    #     #if outside of range let pipeline run with fidcucial params but set flag to indicate -infity likelihood
    #     block[cosmo, 'w'] = -1.
    #     block[cosmo, 'w0'] = -1.
    #     block[cosmo, 'wa'] = 0.
    #     block[cosmo, 'set_likelihood_minus_infinity'] = True 
    # else:
    #     block[cosmo, 'w0'] = w
    #     block[cosmo, 'set_likelihood_minus_infinity'] = False 

    # if(w + wa > config["w_infinite_cutoff"]):
    #     #if outside of range let pipeline run with fidcucial params but set flag to indicate -infity likelihood
    #     block[cosmo, 'set_likelihood_minus_infinity'] = True 
    #     print("w0+wa>-1/3, crashing...")
    #     raise ValueError
        
    # else:
    #     block[cosmo, 'w0'] = w
    #     block[cosmo, 'set_likelihood_minus_infinity'] = False 

    if w + wa > -1./3:
        print("Warning : w0+wa>-1/3")

    block[cosmo, 'set_likelihood_minus_infinity'] = False 
    block[cosmo, 'w0'] = w

    return 0