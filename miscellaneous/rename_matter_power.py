"""
Use the growth factor to extrapolate the P(k) we have to higher redshifts.

"""
import numpy as np
from cosmosis.datablock import option_section, names
from cosmosis.runtime.utils import Timer

def setup(options):
    current_name = options.get_string(option_section, "current_name", "matter_power_lin")
    new_name = options.get_string(option_section, "new_name", "matter_power_nl")

    return {"current_name":current_name, "new_name":new_name}




def execute(block, config):

    z, k, P = block.get_grid(config['current_name'], 'z', 'k_h', 'p_k')
    block.replace_grid(config['new_name'], 'z', z, 'k_h', k, 'p_k', P)

    return 0