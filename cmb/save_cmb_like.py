# save cmb datavector to file. Calculate gaussian covariance matrix

from __future__ import print_function

from builtins import range
from builtins import object
from cosmosis.datablock import option_section, names
import numpy as np
from scipy.interpolate import interp1d
import os

from SO import SO_Noise_Calculator_Public_v3_1_1 as so_models

dirname = os.path.split(__file__)[0]
default_data_file = os.path.join(dirname, "CMB_datavector.txt")
#default_covmat_file = os.path.join(dirname, "covmat.npy")

def get_scales( x_min, x_max, nbins ):
    """
    Get scales
    """
    log_lims = np.linspace( np.log(x_min), np.log(x_max), nbins+1 )
    lims = np.exp(log_lims)
    
    xmin = log_lims[:-1]
    xmax = log_lims[1:]
    log_mids = 0.5*(xmin + xmax)
    mids = np.exp(log_mids)

    return lims, mids

def get_ell_scales( ell_min, ell_max, nbins, ell_as_int=True):
    #ells should be integers. So get scales using get_scales and then 
    #convert
    lims, mids = get_scales( ell_min, ell_max, nbins)
    if (ell_as_int):
        ell_lims = (np.floor(lims)).astype(int)
        ell_mids = (np.floor(mids)).astype(int)
        #np.exp( 0.5 * ( np.log(ell_lims[1:]) + np.log(ell_lims[:-1]) ) )
    else:
        ell_lims = lims 
        ell_mids = mids
    return ell_lims, ell_mids


def setup(options):
    """
    This needs to specify:
    - The section names of the projected power spectra or correlation functions 
    required:
    spectrum_sections = 
    - If you want a covariance, we also need noise parameters. 
    """
    

    def read_list(key):
        s = options.get_string(option_section, key)
        return s.split()

    config = {}
    ell_min = options.get_int(option_section, "ell_min")
    ell_max = options.get_int(option_section, "ell_max")
    
    config["ell_min"] = ell_min
    config["ell_max"] = ell_max

    use_log_spaced = options.get_bool(option_section, "use_log_spaced", True)
    if(use_log_spaced):
        n_ell = options.get_int(option_section, "n_ell")
        ell_lims, ell = get_ell_scales( ell_min, ell_max, n_ell)
        config["n_ell"] = n_ell
        config["ell"] = ell
        config["ell_lims"] = ell_lims
    else:
        config["ell"] = np.arange(ell_min, ell_max+1)
        config["n_ell"] = ell_max-ell_min
        config["ell_lims"] = None
    



    
    #NOISE PARAMETERS: load from SO module or input
    config['fsky'] = options[option_section, "fsky"] 

    config['load_survey_noise'] = options.get_string(option_section, "load_survey_noise", "")
    config['savetnoise'] = options.get_string(option_section, "savetnoise", "") #TODO
    if (config['load_survey_noise'] == ""):
        print("WARNING: this part is not tested. Better use external noise curve")
        config['frequencies'] = options[option_section, "frequencies"] #in GHz
        config['FWHM'] = options[option_section, "FWHM"] #in arcmin
        config['noise'] = options[option_section, "noise"] #in \mu K
        print("assuming given noise is for T, will use sqrt(2) times that for E")
        config['tcmb'] = 2.7255 #K (https://ui.adsabs.harvard.edu/abs/2009ApJ...707..916F/abstract)
        

        config["Pol_lknee"] = options.get_double(option_section, "Pol_lknee")
        config["Pol_alphaknee"] = options.get_double(option_section, "Pol_alphaknee")

        config["Temp_lknee"] = options.get_double(option_section, "Temp_lknee")
        config["Temp_alphaknee"] = options.get_double(option_section, "Temp_alphaknee")
        config["Temp_Nred"] = options[option_section, "Temp_Nred"] 



    
    # name of the output file and whether to overwrite it.
    config['datafile'] = options.get_string(option_section,"data_file", default_data_file)
    #config['covmat_file'] = options.get_string("covmat_file", default_covmat_file)
    #config['overwrite'] = options.get_bool(option_section, "overwrite", False)

    return config

def execute(block, config):
    spectrum_section = "cmb_cl"


    ell_output = config["ell"]


    print("Saving cmb datavector to {}".format(config['datafile']))

    ell = block[spectrum_section, "ell"]
    if (ell[0] > config["ell_min"] or ell[-1] < config["ell_max"]):
        print("The calculated CMB power spectra do not cover the requested range!")
        print(ell[0])
        print(ell[-1])
        return 1
    tmp = (ell >= config["ell_min"])
    tmp2 = (ell <= config["ell_max"])
    mask = np.logical_and(tmp,tmp2)

    ell = ell[mask]
    spectra = {}
    #we want the C_l so remove prefactors
    f = ell * (ell + 1) / (2 * np.pi)
    spectra["tt"] = (block[spectrum_section, "tt"])[mask]  /f
    spectra["te"] = (block[spectrum_section, "te"])[mask]  /f
    spectra["ee"] = (block[spectrum_section, "ee"])[mask]  /f


    noise = {}
    ones = ell/ell 
    if (config['load_survey_noise'] == "SO"):
        
        mode=1 #baseline, 2 for goal
        #suffix='pdf'
        fsky=config["fsky"]
        ellmax=config["ell_max"]+1 #for some reason the SO module cuts off at ell = N-1, so add 1
        el=50. #elevation
        lat = so_models.SOLatV3point1(mode, el=el)
        bands = lat.get_bands()
        print("band centers: ", lat.get_bands(), "[GHz]")
        print("beam sizes: "  , lat.get_beams(), "[arcmin]")
        N_bands = len(bands)
        ell_noise, N_ell_LA_T_full,N_ell_LA_P_full = lat.get_noise_curves(fsky, ellmax, 1, full_covar=True, deconv_beam=True)
        print(ell_noise)
        WN_levels = lat.get_white_noise(fsky)**.5

        #total noise on power spectra given by inverse sum
        N_ell_LA_T  = N_ell_LA_T_full[range(N_bands),range(N_bands)]
        N_ell_LA_P  = N_ell_LA_P_full[range(N_bands),range(N_bands)]

        N_ell_LA_T_total = 1./np.sum(1./N_ell_LA_T, axis = 0)
        N_ell_LA_P_total = 1./np.sum(1./N_ell_LA_P, axis = 0)

        print("white noise levels: "  , WN_levels, "[uK-arcmin]")

        mask = np.where(np.in1d(ell, ell_noise))[0]
        noise["tt"] = N_ell_LA_T_total[mask]
        noise["te"] = 0.*ones #off diagonal!
        noise["ee"] = N_ell_LA_P_total[mask]
    else:
        print("No survey to load found. Will try to estimate noise from parameters directly. NOT TESTED")
        noise["tt"] = get_noise(ell,config['frequencies'], config['FWHM'], config['noise'], config['tcmb'] ,config["Temp_lknee"] , config["Temp_alphaknee"] , config["Temp_Nred"] , config['savetnoise'])
        noise["te"] = ones #off diagonal!
        noise["ee"] = get_noise(ell,config['frequencies'], config['FWHM'], np.sqrt(2.)*config['noise'], config['tcmb'] ,config["Pol_lknee"] , config["Pol_alphaknee"] , None)


    #units still need work!!!

    spectra_hat = {}
    print("tt spectra!")
    print(spectra["tt"]*ell * (ell + 1) / (2 * np.pi))
    print("tt noise")
    print(noise["tt"]*ell * (ell + 1) / (2 * np.pi))

    spectra_hat["tt"] = spectra["tt"]+noise["tt"]
    spectra_hat["te"] = spectra["te"]+noise["te"]
    spectra_hat["et"] = spectra_hat["te"] ###??????????NEED THIS? prob not
    spectra_hat["ee"] = spectra["ee"]+noise["ee"]

    
    types = ["tt", "te", "ee"]
    n_ell = len(ell_output)
    n_cov = n_ell*len(types)
    covmat = np.zeros((n_cov, n_cov))
    for n1,ij in enumerate(types):
        for n2,kl in enumerate(types):
            dig_cov_ijkl = get_cov_diag_ijkl(spectra_hat, ij[0],ij[1],kl[0],kl[1], config["fsky"], ell)
            dig_cov_ijkl_summed = sum_ell(dig_cov_ijkl, ell, ell_output, config["ell_lims"])
            #find where to put the diagonal submatrix
            print(dig_cov_ijkl_summed)
            inds_x = np.arange( n_ell*n1, n_ell*(n1+1) )
            inds_y = np.arange( n_ell*n2, n_ell*(n2+1))
            cov_inds = np.ix_( inds_x, inds_y )
            covmat[ cov_inds ] = np.diag(dig_cov_ijkl_summed)
    
    #Combine davavector (ell,Cl) and Cov to one matrix
    datavector = np.zeros( (n_cov+2, n_cov))
    inds_x = np.arange( 2, n_cov+2)
    inds_y = np.arange( 0, n_cov)
    cov_inds = np.ix_( inds_x, inds_y )
    datavector[ cov_inds ] = covmat
    datavector[0,0:n_ell] = ell_output
    datavector[0,n_ell:n_ell*2] = ell_output
    datavector[0,n_ell*2:n_ell*3] = ell_output

    output_mask = np.where(np.in1d(ell, ell_output))[0]
    if(len(output_mask) < len(ell_output)):
        print("Requested too many log spaced samples, some of them overlap!")
        print(ell)
        print(ell_output)
        print(output_mask)
    datavector[1,0:n_ell] = spectra["tt"][output_mask]
    datavector[1,n_ell:n_ell*2] = spectra["te"][output_mask]
    datavector[1,n_ell*2:n_ell*3] = spectra["ee"][output_mask]

    print(datavector.size)
    np.savetxt(config['datafile'],datavector)

    return 0

def get_noise(ell,f, FWHM, noise_band, Tcmb, lknee, alphaknee, Nred = None, savenoise="" ):
    #not tested, currently not used

    #N_{\ell, c}^{Y Y}=\left(\frac{\sigma_{c}^{Y} \theta_{\mathrm{FWHM}, c}}{T_{\mathrm{CMB}}}\right) e^{\ell(\ell+1) \theta_{\mathrm{FWHM}, \mathrm{c}}^{2} / 8 \ln 2}
    arcmin_to_rad = 2.*np.pi/(360.*60.) # approx 0.000290888
    print(arcmin_to_rad)
    #noise_white_per_band = np.array([ (noise_band[i]* 10**(-6) * FWHM[i]  / Tcmb) * np.exp(ell*(ell+1) * FWHM[i]*FWHM[i] * arcmin_to_rad*arcmin_to_rad /(8.*np.log(2.))) for i in range(len(f))])
    #noise_per_band = np.array([ (noise[i]* FWHM[i] *arcmin_to_rad / Tcmb) * np.exp(ell*(ell+1) * FWHM[i]*FWHM[i] * arcmin_to_rad*arcmin_to_rad /(8.*np.log(2.))) for i in range(len(f))])
    #noise_white_per_band = np.array([ (noise_band[i] * arcmin_to_rad )**2. * np.exp(ell*(ell+1) * FWHM[i]*FWHM[i] * arcmin_to_rad*arcmin_to_rad /(8.*np.log(2.))) for i in range(len(f))])
    noise_per_band = np.array([ (noise_band[i] * arcmin_to_rad )**2.  for i in range(len(f))])
    
    print("You have to be a bit carefull here. Check the units of noise you are using. ")
    print(np.exp(ell*(ell+1) * FWHM[0]*FWHM[0] * arcmin_to_rad*arcmin_to_rad /(8.*np.log(2.))))
    #N_{\ell}=N_{\mathrm{red}}\left(\frac{\ell}{\ell_{\mathrm{knee}}}\right)^{\alpha_{\mathrm{knee}}}+N_{\mathrm{white}}
    # if(Nred is None):
    #     noise_red_per_band = noise_white_per_band
    # else:
    #     noise_red_per_band = Nred 
    # noise_per_band = np.array([ noise_red_per_band[i] * (ell/lknee)**alphaknee+ noise_white_per_band[i] for i in range(len(f))   ])


    #print((noise_band[0] * FWHM[0] / Tcmb))
    print("focus only on the l=1000 numbers")
    print(ell)
    test =noise_per_band[0]*ell * (ell + 1) / (2 * np.pi)
    print(test)
    print(test[970])
    print(f[-1])
    test = noise_per_band[-1]*ell * (ell + 1) / (2 * np.pi)
    print(test)
    print(test[970])

    #N_{\ell}^{Y Y}=\left[\sum_{c}\left(N_{\ell, c}^{Y Y}\right)^{-1}\right]^{-1}
    noise = 1./np.sum(1./noise_per_band , axis=0)



    
    if(savenoise != ""):
        np.savetxt(savenoise+"_ell.txt", ell)
        np.savetxt(savenoise+"_noise.txt", noise)
    
    print(noise*ell * (ell + 1) / (2 * np.pi))
    return noise 

def get_cov_diag_ijkl(spectra_hat, i,j,k,l, fsky, ell):
    #assuming no bin averaging in l
    prefactor = 1./(fsky*(2*ell+1))
    print(i+j)
    print(k+l)
    diga_cov_ijkl = prefactor * (spectra_hat[i+k]*spectra_hat[j+l] + spectra_hat[i+l]*spectra_hat[j+k])
    print(np.sum(diga_cov_ijkl < 0.))
    return diga_cov_ijkl

def sum_ell(dig_cov, ell, ell_output, ell_lims):
    #don't output full datavector. need to sum. assume that ell already has the right total range.
    #\text{Var}(\bar C_l) = \frac{\sum_{i=n}^m(2i+1)^2 \text{Var}(C_i)}{(\sum_{n=i}^m (2i+1))^2}
    if(len(ell) == len(ell_output)):
        #do nothing in case the log spaced binning is turned off
        return dig_cov

    dig_cov_l2plus1square = (2*ell+1)**2 *dig_cov
    l2plus1 = (2*ell+1)
    
    dig_cov_l2plus1square_summed = bin_array(dig_cov_l2plus1square, bins = ell_lims-ell_lims[0])
    l2plus1_summmed = bin_array(l2plus1, bins = ell_lims-ell_lims[0])
    diag_cov_summed = dig_cov_l2plus1square_summed / l2plus1_summmed**2
    return diag_cov_summed

def bin_array(arr, bins):
    idx = range(len(arr))
    return np.array([np.sum(arr[idx[bins[i]]:idx[bins[i+1]]]) for i in range(len(bins)-1)])








