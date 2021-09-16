"""
Replace P(k) by P(k)*(D(k,z)/D(k=kmin, z))^2. To be used in combination with eisenstein_hu_cdm.

"""
import numpy as np
from scipy.interpolate.interpolate import interp2d, interp1d

do_plot = False # for debugging only
growth_section = "GROWTH_PARAMETERS"
MG = 'modified_gravity_parameters'

def setup(options):
    return {}

def execute(block, config):

    if block[MG, "model"] == 2.0: 

        z, k, P = block.get_grid('matter_power_lin', 'z', 'k_h', 'p_k')

        d_z = get_d_z_splined(block, z)
        d_z_k = get_d_z_k_splined(block, z, k)

        P_new = P * (d_z_k / d_z[:,np.newaxis])**2
        block.replace_grid('matter_power_lin', 'z', z, 'k_h', k, 'p_k', P_new)

        if do_plot == True:
            mkdir_p('./tmp')
            pfname = './tmp/plot_matter_power_vs_GR.png'
            plot_check_p(z, k, P_new, pfname=pfname)
            plot_check_p_from_file(pfname=pfname)

        return 0
    else: 
        return 0

def get_d_z_splined(block, z):

    z_for_d_z = block[growth_section, 'z']
    d_z = block[growth_section, 'd_z']

    func = interp1d(z_for_d_z, np.log(d_z), kind='cubic')
    d_z_splined = np.exp(func(z))

    return d_z_splined

def get_d_z_k_splined(block, z, k):

    k2, z2, d_k_z = block.get_grid(growth_section, 'k_for_d_k_z', 'z_for_d_k_z', 'd_k_z')
    d_z_k = np.transpose(d_k_z)

    func = interp2d(np.log(k2), z2, np.log(d_z_k), kind='cubic')
    d_z_k_splined = np.exp(func(np.log(k), z))

    if do_plot == True:
        plot_check_interp_d(z2, k2, d_z_k, z, k, d_z_k)

    return d_z_k_splined


def plot_check_interp_d(z2, k2, d_z_k, z, k, d_z_k_splined):

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()

    for iz in range(0, z.size, 49):
        z_tmp = z[iz]
        ind = np.min(np.where(z2 > z_tmp)[0])

        ax.loglog(k2, d_z_k[ind-1,:], ls='-', color='black', label='z=%.2f'%z2[ind-1])
        ax.loglog(k, d_z_k_splined[iz,:], ls=':', lw=2, label='z=%.2f, splined'%z_tmp)
        ax.loglog(k2, d_z_k[ind,:], ls='--', color='grey', label='z=%.2f'%z2[ind])
        
    ax.legend()
    ax.set_xlabel(r'$k [h/Mpc]$')
    ax.set_ylabel(r'$D(k;z) [(Mpc/h)^3]$')
    
    mkdir_p('./tmp')
    pfname = './tmp/plot_check_interp_d_k_z.png'
    fig.savefig(pfname)
    print('Saved plot: {}'.format(pfname))

def plot_check_interp_p(z, k, P, P_new):

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    fig = plt.figure()
    gspec = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, \
        wspace=0.0, hspace=0.0, height_ratios=[4,1.5])
    ax1 = fig.add_subplot(gspec[0, 0])
    ax2 = fig.add_subplot(gspec[1, 0])

    ax1.loglog([], [], ls='-', color='k', lw=2, label='scale-indep.')
    ax1.loglog([], [], ls=':', color='k', lw=1, label='scale-dep.')

    for iz in range(0, z.size, 20):
        z_tmp = z[iz]
        ax1.loglog(k, P[iz,:], ls='-', lw=2, label='z=%.2f'%z_tmp)
        ax1.loglog(k, P_new[iz,:], ls=':', color='k', lw=1)
    
    ax1.legend(ncol=2)
    ax1.set_ylabel(r'$P(k) [(Mpc/h)^3]$')
    
    for iz in range(0, z.size, 20):
        z_tmp = z[iz]
        frac_diff = (P_new[iz,:] - P[iz,:])/P[iz,:]
        ax2.semilogx(k, frac_diff, ls='-', lw=1)
    
    ax2.set_xlabel(r'$k [h/Mpc]$')
    ax2.set_ylabel(r'$\Delta P/P$')
    ax2.set_ylim([-0.05, 0.3])
    
    mkdir_p('./tmp')
    pfname = './tmp/plot_check_interp_matter_power.png'
    fig.savefig(pfname)
    print('Saved plot: {}'.format(pfname))

def load_data(dir_name, fname):
    import os
    fn = os.path.join(dir_name, fname)
    return np.genfromtxt(fn, skip_header=1)

def plot_check_p(z, k, P_new, pfname):

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    dir_GR = './6x2pt_Roman_SO_GR_log_k/matter_power_lin/'
    P_GR = load_data(dir_GR, 'p_k.txt')
    k_GR = load_data(dir_GR, 'k_h.txt')
    z_GR = load_data(dir_GR, 'z.txt')

    func_P_GR = interp2d(np.log(k_GR), z_GR, np.log(P_GR), kind='cubic')
    P_GR_splined = np.exp(func_P_GR(np.log(k), z))

    fig = plt.figure()
    gspec = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, \
        wspace=0.0, hspace=0.0, height_ratios=[4,1.5])
    ax1 = fig.add_subplot(gspec[0, 0])
    ax2 = fig.add_subplot(gspec[1, 0])

    ax1.loglog([], [], ls='-', color='k', lw=2, label='f(R)')
    ax1.loglog([], [], ls=':', color='k', lw=1, label='GR')

    for iz in range(0, 100, 20): #z.size,
        z_tmp = z[iz]
        ax1.loglog(k, P_new[iz,:], ls='-', lw=2, label='z=%.2f'%z_tmp)
        ax1.loglog(k, P_GR_splined[iz,:], ls=':', color='k', lw=1)
    
    ax1.legend(ncol=2)
    ax1.set_ylabel(r'$P(k) [(Mpc/h)^3]$')
    ax1.set_ylim([1e-3, 1e5])
    
    for iz in range(0, 100, 20): #z.size
        z_tmp = z[iz]
        frac_diff = (P_new[iz,:]/P_GR_splined[iz,:])-1
        ax2.semilogx(k, frac_diff, ls='-', lw=1)
        ax2.grid(True)
    
    ax2.set_xlabel(r'$k [h/Mpc]$')
    ax2.set_ylabel(r'$\Delta P/P$')
    ax2.set_ylim([-0.05, 0.3])

    ax1.set_xlim([1e-5, 10])
    ax2.set_xlim([1e-5, 10])
    
    fig.savefig(pfname)
    print('Saved plot: {}'.format(pfname))

def plot_check_p_from_file(pfname):
    dir_2 = './6x2pt_Roman_SO_fof_R_n_1_fR_1e-4_modified_growth_log_k/matter_power_lin/'
    #dir_2 = './6x2pt_Roman_SO_fof_R_n_1_fR_1e-4_modified_growth_log_k_no_sigma8_rescale/matter_power_lin/'
    P_2 = load_data(dir_2, 'p_k.txt')
    k_2 = load_data(dir_2, 'k_h.txt')
    z_2 = load_data(dir_2, 'z.txt')
    pfname = './tmp/plot_matter_power_from_file_vs_GR.png'
    plot_check_p(z_2, k_2, P_2, pfname=pfname)

def mkdir_p(path):
    import os, errno
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise