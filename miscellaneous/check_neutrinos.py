"""
Apply a cutoff to the sum of w0 and wa. 

"""
import numpy as np
from cosmosis.datablock import option_section, names
from cosmosis.runtime.utils import Timer

cosmo = names.cosmological_parameters

def setup(options):
    omnuh2_cutoff = options.get_double(option_section, "omnuh2_cutoff", 0.00003)
    #new_name = options.get_double(option_section, "new_name", "matter_power_nl")

    return {"omnuh2_cutoff":omnuh2_cutoff}




def execute(block, config):

    omnuh2 = block[cosmo, 'omnuh2']

    print("Hello Lukas")
    if(omnuh2 > 0. and omnuh2 < config["omnuh2_cutoff"]):
        #low neutrinos do not work with out code -> set to zero below cutoff
        print("HELLO")
        #block[cosmo, 'omnuh2'] = 0.00003
        block[cosmo, 'omnuh2'] = 0.


    return 0