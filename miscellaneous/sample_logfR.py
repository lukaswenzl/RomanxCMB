"""
Sample f_of_R_fR in logspace instead of linear

"""
import numpy as np
from cosmosis.datablock import option_section, names
from cosmosis.runtime.utils import Timer

cosmo = names.cosmological_parameters
MG = "modified_gravity_parameters"

def setup(options):
    sample_logspace_flag = options.get_bool(option_section, "sample_logspace",True)
    #new_name = options.get_double(option_section, "new_name", "matter_power_nl")

    return {"sample_logspace_flag":sample_logspace_flag}




def execute(block, config):

    #only apply if we actually use modified gravity
    if(block.has_value(MG, 'model') and block[MG, "model"] >= 0):
        #check if the flag is set to true
        if(config["sample_logspace_flag"]):
            if(block.has_value(MG, 'f_of_R_fR') and block[MG, 'f_of_R_fR'] != -1):
                print("Warning: Set to sample f_of_R_fR in logspace, will overwrite linear value!")
            #convert the low space value to linear
            block[MG, 'f_of_R_fR'] = np.power(10.,block[MG, 'logf_of_R_fR'] )
    

    return 0