#
"""
quick script to calculate the scale cuts once for our given fiducial cosmology.
Only needs to run once and then we have a file with the scale cuts

"""
import numpy as np
from cosmosis.datablock import option_section, names
from cosmosis.runtime.utils import Timer

cosmo = names.cosmological_parameters

def setup(options):
    galaxy_scale_cutoff_mpcoverh = options.get_double(option_section, "galaxy_scale_cutoff_mpcoverh", -1)
    ell_max = options.get_double(option_section, "ell_max", -1)
    ell_min = options.get_double(option_section, "ell_min", -1)
    if(galaxy_scale_cutoff_mpcoverh == -1 or ell_max == -1 or ell_min == -1):
        print("Not all configuration parameters given, results will be nonsense or crash.")

    return {"galaxy_scale_cutoff_mpcoverh":galaxy_scale_cutoff_mpcoverh, "ell_min":ell_min, "ell_max":ell_max}




def execute(block, config):
    Rmax = config["galaxy_scale_cutoff_mpcoverh"]
    ell_min = config["ell_min"]
    ell_max = config["ell_max"]


    z = block["distances", 'z']
    chi = block["distances", 'd_m']

    #naming confusing: clustering <-> lenses
    nbin = block["nz_lens", "nbin"] #not sure how to access the values
    #print(nbin)
    clustering_bins = [block["nz_lens", "bin_"+str(i)] for i in range(1,nbin+1)]
    z_bins = block["nz_lens", "z"]

    mean_z = np.zeros(10)
    for i, bin_i in enumerate(clustering_bins):
        norm = np.sum(bin_i)
        bin_i_normed = bin_i/norm
        mean_z[i] = np.sum(bin_i_normed * z_bins)


    print(mean_z)
    chi_zmean = np.interp(mean_z, z, chi)
    print(chi_zmean)

    

    #finding lmax using limber approx: k = l/chi(z mean)
    #be careful with 1/h units!!
    h =  block[cosmo, 'h0']
    Rmax = Rmax * h
    kmax = 2*np.pi/Rmax
    ell_cutoff = kmax * chi_zmean
    print(ell_cutoff)

    #loop through all relevant power spectra
    print("limiting clustering to Rmax = {} Mpc".format(Rmax))
    print("You will need to copy the following scale cuts:")
    cuts = []

    section = "galaxy_cl"
    for i in range(0,nbin):
        for j in range(0,nbin):
            if block.has_value(section, "bin_{}_{}".format(i+1, j+1)):
                if(ell_cutoff[i] < ell_max or ell_cutoff[j] < ell_max ):
                    limiting_ell = np.min([ell_cutoff[i], ell_cutoff[j]])
                    cuts.append("angle_range_{}_{}_{}= {:.1f} {:.1f}".format(section, i+1, j+1, ell_min, limiting_ell))
                    print(cuts[-1])

    #only doing the cutoff for clustering
    section = "galaxy_shear_cl"
    for i in range(0,nbin):
        for j in range(0,nbin):
            if block.has_value(section, "bin_{}_{}".format(i+1, j+1)):
                if(ell_cutoff[i] < ell_max):
                    limiting_ell = ell_cutoff[i]
                    cuts.append("angle_range_{}_{}_{}= {:.1f} {:.1f}".format(section, i+1, j+1, ell_min, limiting_ell))
                    print(cuts[-1])

    section = "galaxy_cmbkappa_cl"
    for i in range(0,nbin):
        j = 0 # there is only one cmb kappa bin
        if block.has_value(section, "bin_{}_{}".format(i+1, j+1)):
            if(ell_cutoff[i] < ell_max):
                limiting_ell = ell_cutoff[i]
                cuts.append("angle_range_{}_{}_{}= {:.1f} {:.1f}".format(section, i+1, j+1, ell_min, limiting_ell))
                print(cuts[-1])

    print("")
    print("")
    print("reached end of scale cuts script. Copy paste them into the params file")
    
    

    
    # galaxy_shear_cl 
    # galaxy_cl
    # galaxy_cmbkappa_cl

    # if(w + wa > config["w_infinite_cutoff"]):
    #     #if outside of range let pipeline run with fidcucial params but set flag to indicate -infity likelihood
    #     block[cosmo, 'w'] = -1.
    #     block[cosmo, 'w0'] = -1.
    #     block[cosmo, 'wa'] = 0.
    #     block[cosmo, 'set_likelihood_minus_infinity'] = True 
    # else:
    #     block[cosmo, 'w0'] = w
    #     block[cosmo, 'set_likelihood_minus_infinity'] = False 

    return 0