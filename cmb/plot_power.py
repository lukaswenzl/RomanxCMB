import numpy as np 
import matplotlib as plt

data_file = self.options.get_string("data_file", default_data_file)
ell, c_ell = np.loadtxt(data_file)[0:2] #.T
te_start, ee_start = np.where(np.diff(ell) < 0)[0] + 1
ell = np.split(ell, [te_start, ee_start])