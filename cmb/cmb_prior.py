from cosmosis.gaussian_likelihood import GaussianLikelihood
from cosmosis.datablock import names
from scipy.interpolate import interp1d
import os
import numpy as np

dirname = os.path.split(__file__)[0]
default_data_file = os.path.join(dirname, "CMB_prior.txt")
default_data_file_names = os.path.join(dirname, "CMB_prior_names.txt")

#default_covmat_file = os.path.join(dirname, "covmat.npy")


class CMBLikelihood(GaussianLikelihood):
    like_name = "cmb_prior"
    #x_section = names.cosmological_parameters
    #y_section = names.cosmological_parameters
    #x_name = "ell"
    cosmo = names.cosmological_parameters

    def build_data(self):
        data_file = self.options.get_string("data_file", default_data_file)
        data_file_names = self.options.get_string("data_file_names", default_data_file_names)

        params_names = np.loadtxt(data_file_names, dtype=str) #.T
        params_best = np.loadtxt(data_file)[0] #.T

        #te_start, ee_start = np.where(np.diff(ell) < 0)[0] + 1
        #ell = np.split(ell, [te_start, ee_start])
        #return ell, c_ell
        return params_names, params_best

    def build_covariance(self):
        # covmat should be stored in npy format for size
        covmat_file = self.options.get_string(
            "data_file", default_data_file)
        covmat = np.loadtxt(covmat_file)[1:]
        return covmat

    def extract_theory_points(self, block):
        # Extract the parameters we are testing in this run
        params_names = self.data_x

        params_test = np.array([block[self.cosmo, i] for i in params_names])
        #to get Omegamh2 maybe I need to somehow access consistency section?
        #print("!!still need to check that is is ok to multiply prior to like!") TODO
        return params_test


setup, execute, cleanup = CMBLikelihood.build_module()
