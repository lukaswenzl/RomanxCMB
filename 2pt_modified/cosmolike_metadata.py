# script to match cosmosis and cosmolike
# part of 6x 2pt forecast for Roman_SO
# determine bin cuts
# Author Lukas Wenzl


import numpy as np

def load_covariance_metadata(filename="modules/RomanxCMB/cosmolike_data/cov_indices_testruns_apr9.txt"):
    #index starts at 0 for cosmlike
    info = np.genfromtxt(filename, delimiter=" ", dtype=None, names=("cov_idx","ell_idx","ell", "bin1","bin2","tracers"))

    #need to reverse ks -> sk
    idx = np.argwhere(info["tracers"]=="ks")
    tmp = info[idx] 
    tmp["tracers"] = "sk"
    tmp2 = np.copy(tmp["bin1"])
    tmp["bin1"] = np.copy(tmp["bin2"])
    tmp["bin2"] = tmp2 
    np.put(info, idx, tmp)

    #calculate the identifier string
    tracers = {"ss":"shear_cl", "ls":"galaxy_shear_cl", "ll":"galaxy_cl", "lk":"galaxy_cmbkappa_cl", "sk":"shear_cmbkappa_cl", "kk":"cmbkappa_cl"}
    metadata_identifiers = np.array([tracers[i["tracers"]]+str(i["bin1"]+1)+str(i["bin2"]+1)+str(i["ell_idx"])  for i in info])
    return info, metadata_identifiers

def build_metadata_for_cosmosis(spectra, n_ell,ignore_ells=False):
    combinations = []
    identifiers = []
    #identifiers_reverse = []
    for spectrum in spectra:
        for b in spectrum.bin_pairs:
            if(ignore_ells):
                    identifiers.append(spectrum.name+str(b[0])+str(b[1])+str(0))
                    combinations.append((spectrum.name, b[0], b[1]))
            else:
                for i in range(n_ell): #ell index starting from 0
                    #galaxy_cl110 
                    identifiers.append(spectrum.name+str(b[0])+str(b[1])+str(i))
                    #identifiers_reverse.append(spectrum.name+str(b[1])+str(b[0])+str(i))
                    combinations.append((spectrum.name, b[0], b[1]))

    return combinations, identifiers

def give_cuts(metadata_filename, spectra, n_ell):
    #index starts at 1 for cosmosis?
    #format: (spectrum.name, b1, b2)
    info, metadata_identifiers = load_covariance_metadata(metadata_filename)
    
    combinations, identifiers = build_metadata_for_cosmosis(spectra, n_ell, ignore_ells=True)        
    to_cut = np.isin(identifiers, metadata_identifiers, invert=True) 
    #b = np.isin(identifiers_reverse, metadata_identifiers, invert=True) #do i need the inverse test?
    #to_cut = np.logical_and(a,b) 
    # cuts = np.array(cuts)[to_cut] 
    # #needs to be a list              
    cuts = []
    for i, cut in enumerate(combinations):
        if(to_cut[i]):
            cuts.append(cut)
    return cuts


def load_covariance(cov_filename):
    covariance =  np.loadtxt(cov_filename)
    print(covariance)
    print(covariance.shape)
    return covariance

def rearrange_cov(cov_filename, spectra, metadata_filename, n_ell):
    info, metadata_identifiers = load_covariance_metadata(filename=metadata_filename)
    combinations, cosmosis_identifiers = build_metadata_for_cosmosis(spectra,n_ell, ignore_ells=False)

    if(len(cosmosis_identifiers) != len(metadata_identifiers)):
        print("WARNING: cosmosis data and cosmolike covariance matrix do not have same number of indices. Something went wrong")

    #match by sorting and unsorting
    sort_cosmosis = np.argsort(cosmosis_identifiers)
    unsort_cosmosis = np.argsort(sort_cosmosis)
    sort_metadata = np.argsort(metadata_identifiers)
    #unsort_metadata = np.argsort(sort_metadata)

    reshape = sort_metadata[unsort_cosmosis]
    cov = load_covariance(cov_filename)
    cov = cov[reshape, :]#tested with gaussian covariance matrices and it works!
    cov = cov[:, reshape]
    return cov


# info = load_covariance_metadata()
# print(info)

# print(rearange_cov())