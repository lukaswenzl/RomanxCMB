from cosmosis.datablock import names, option_section
from linear_alignments import kirk_rassat_host_bridle_power
from linear_alignments import bridle_king
from linear_alignments import bridle_king_corrected
from linear_alignments import linear
from linear_alignments import krause_eifler_blazek
from luminosity_dependance import setup_luminosity_dependance
from luminosity_dependance import resample
import numpy as np


def setup(options):
    method = options[option_section, "method"].lower()
    grid_mode = options.get_bool(option_section, "grid_mode", default=False)
    gal_intrinsic_power = options.get_bool(
        option_section, "do_galaxy_intrinsic", False)
    name = options.get_string(option_section, "name", default="").lower()
    if name:
        suffix = "_" + name
    else:
        suffix = ""

    if method not in ["krhb", "bk", "bk_corrected", "linear", "keb"]:
        raise ValueError('The method in the linear alignment module must'
                         'be either "KRHB" (for Kirk, Rassat, Host, Bridle) or BK for '
                         'Bridle & King or "BK_corrected" for the corrected version of that or keb for krause_eifler_blazek 2016')

    if method == "keb":
        pathtok_e_corr = options[option_section, "k_e_correction_folder"]
        #load k correction files for Limiting magnitude calculation
        k_corr, e_corr = setup_luminosity_dependance(band = "r", pathtok_e_corr=pathtok_e_corr)
    else:
        k_corr = None
        e_corr = None

    return method, gal_intrinsic_power, grid_mode, suffix,k_corr, e_corr


def execute(block, config):
    method, gal_intrinsic_power, grid_mode, suffix,k_corr, e_corr = config

    # load z_lin, k_lin, P_lin, z_nl, k_nl, P_nl, C1, omega_m, H0
    lin = names.matter_power_lin
    nl = names.matter_power_nl
    ia = names.intrinsic_alignment_parameters + suffix
    ia_ii = names.intrinsic_power + suffix
    ia_gi = names.galaxy_intrinsic_power + suffix
    ia_mi = names.matter_intrinsic_power + suffix
    gm = "matter_galaxy_power" + suffix
    cosmo = names.cosmological_parameters

    z_lin, k_lin, p_lin = block.get_grid(lin, "z", "k_h", "p_k")
    z_nl, k_nl, p_nl = block.get_grid(nl, "z", "k_h", "p_k")

    omega_m = block[cosmo, "omega_m"]
    A = block[ia, "A"]

    # run computation and write to datablock
    if method == 'krhb':
        P_II, P_GI, b_I, r_I, k_I, z_I = kirk_rassat_host_bridle_power(
            z_lin, k_lin, p_lin, z_nl, k_nl, p_nl, A, omega_m)
    elif method == 'bk':
        P_II, P_GI, b_I, r_I, k_I, z_I = bridle_king(
            z_nl, k_nl, p_nl, A, omega_m)
    elif method == 'bk_corrected':
        P_II, P_GI, b_I, r_I, k_I, z_I = bridle_king_corrected(
            z_nl, k_nl, p_nl, A, omega_m)
    elif method == "linear":
        P_II, P_GI, b_I, r_I, k_I, z_I = linear(
            z_lin, k_lin, p_lin, A, omega_m)
    elif method == "keb":
        h = block[cosmo, "h0"]
        mlim = block[ia, "mlim"] #fixed for survey
        z0_IA = block[ia, "z0"] #fixed
        z1_IA = block[ia, "z1"] #fixed
        M0 = block[ia, "M0"] #fixed

        beta = block[ia, "beta"]
        eta = block[ia, "eta"]
        eta_highz = block[ia, "eta_highz"]

        #todo: exception when luminosity distance not found
        D_L_values = block["distances", "d_l"]
        z_luminosity_dist = block["distances", "z"]

        #cutoff on redsift: IA beyond redshift 10 is irrelevant for us, set to 0
        cutoff_redshift = block[ia, "cutoff_redshift"]
        z_nl_cutoff = z_nl[z_nl<= cutoff_redshift]

        #resample we only need the z_nl_cutoff sampling
        D_L_values = resample(D_L_values, z_luminosity_dist, z_nl_cutoff, extrapolate_linearly_at_highz=False)
        #if(z_nl[0]==0):
        #    D_L_values[0] = 1000#this might not be a good idea
        D_L = {"D_L":D_L_values, "z":z_nl_cutoff}

        P_II, P_GI, b_I, r_I, k_I, z_I = krause_eifler_blazek(
            z_nl,z_nl_cutoff, k_nl, p_nl, A, omega_m, h, z0_IA,z1_IA,M0,beta,eta,eta_highz,mlim, D_L,k_corr, e_corr)

    if grid_mode:
        block.put_grid(ia, "z", z_I, "k_h", k_I,  "b_I" + suffix, b_I)
        block.put_grid(ia, "z", z_I, "k_h", k_I, "r_I" + suffix, r_I)
    else:
        block.put_grid(ia_mi, "z", z_I, "k_h", k_I,  "p_k", P_GI)
        block.put_grid(ia_ii, "z", z_I, "k_h", k_I, "p_k", P_II)

    # This is a bit of hack...scale GI power spectrum (which is really matter-intrinsic
    # power spectrum) by P_gal_matter/P_delta_delta
    if gal_intrinsic_power:
        z, k, p_gm = block.get_grid(gm, "z", "k_h", "p_k")
        P_gI = P_GI * p_gm / p_nl
        block.put_grid(ia_gi, "z", z, "k_h", k, "p_k", P_gI)
    return 0


def cleanup(config):
    pass
