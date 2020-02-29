# script to translate cosmolike covaraince into cosmosis data vector, will turn into module
# part of 6x 2pt forecast for WFIRST_SO
# Author Lukas Wenzl

import numpy as np
from astropy.io import fits
from cosmosis.datablock import names, option_section


def setup(options):
    #This function is called once per processor per chain.
    #It is a chance to read any fixed options from the configuration file,
    #load any data, or do any calculations that are fixed once.

    #Use this syntax to get a single parameter from the ini file section
    #for this module.  There is no type checking here - you get whatever the user
    #put in.
    config = {}
    config["FILE_2pt"] = options.get_string(option_section, "FILE_2pt", "6x2pt_WFIRST_SO.fits")
    # DATA_SETS_2pt = "shear_cl galaxy_shear_cl galaxy_cl shear_cmblensing_cl galaxy_cmblensing_cl cmblensing_cl"
    # DATA_SETS_nspectra = "55 32 10 10 10 1"

    config["COVARIANCE"] = options.get_string(option_section, "COVARIANCE","covariances_WFIRST_SO_gold.txt")

    def read_list(key):
        s = options.get_string(option_section, key)
        return s.split()
    config["DATA_SETS_2pt"] = read_list("DATA_SETS_2pt")
    config["DATA_SETS_nspectra"] = options.get(option_section, "DATA_SETS_nspectra")
    print(config["DATA_SETS_2pt"])
    print(config["DATA_SETS_nspectra"])
    config["n_ell"] = options.get_int(option_section, "n_ell", -1)

    config["OVERWRITE"] = options.get_bool(option_section, "OVERWRITE", True)
    config["FILE_2pt_OUTPUT"] = options.get_string(option_section,"FILE_2pt_OUTPUT", config["FILE_2pt"])


    return config


def execute(block, config):
    #This function is called every time you have a new sample of cosmological and other parameters.
    #It is the main workhorse of the code. The block contains the parameters and results of any
    #earlier modules, and the config is what we loaded earlier.
    covariance =  np.loadtxt(config["COVARIANCE"])
    hdul = fits.open(config["FILE_2pt"])
    print(hdul.info())


    #checking consistency
    try:
        if (len(config["DATA_SETS_2pt"]) != len(config["DATA_SETS_nspectra"])):
            print("WARNING: DATA_SETS_2pt and DATA_SETS_nspectra do not have same length")
        total_spectra_points = sum(config["n_ell"]*config["DATA_SETS_nspectra"])
        if( total_spectra_points != covariance.shape[0] or  total_spectra_points != covariance.shape[1]):
            if(config["n_ell"]== -1):
                print("WARNING: Number of samples per spectrum not specified, can not check it the covariance matrix is consistent")
            else:
                print("WARNING: number of spectra not the same as in covaraince")
            #todo: check number of spectra in datavector
        print("Checked consistency")
    except:
        print("consistency check failed")
        print("Something is wrong")

    hdu = fits.ImageHDU(covariance, name="COVMAT")

    hdr = hdu.header
    hdr['COVDATA'] = True
    count_spectra = 0 #the fits file header should contain where each type of spectra starts, the name for this is STRT_i
    for i in range(len(config["DATA_SETS_2pt"])):
        hdr["NAME_"+str(i)] = config["DATA_SETS_2pt"][i]
        hdr["STRT_"+str(i)] = count_spectra
        count_spectra += config["DATA_SETS_nspectra"][i] * config["n_ell"]

    hdu.header = hdr
    hdul.insert(1,hdu)

    print(hdul.info())
    hdul.writeto(config["FILE_2pt_OUTPUT"], overwrite=config["OVERWRITE"])
    print("wrote file to "+config["FILE_2pt_OUTPUT"])
    hdul.close()
    #We tell CosmoSIS that everything went fine by returning zero
    return 0

def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass



#need this in header:
# >>> hdul[1].header
# XTENSION= 'IMAGE   '           / Image extension
# BITPIX  =                  -64 / array data type
# NAXIS   =                    2 / number of array dimensions
# NAXIS1  =                  900
# NAXIS2  =                  900
# PCOUNT  =                    0 / number of parameters
# GCOUNT  =                    1 / number of groups
# COVDATA =                    T
# EXTNAME = 'COVMAT  '
# STRT_0  =                    0
# NAME_0  = 'xip     '
# STRT_1  =                  200
# NAME_1  = 'xim     '
# STRT_2  =                  400
# NAME_2  = 'gammat  '
# STRT_3  =                  800
