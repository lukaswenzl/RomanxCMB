from __future__ import print_function
from builtins import str
from cosmosis.runtime.config import Inifile
from cosmosis.runtime.pipeline import LikelihoodPipeline
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d

output_folder = "modules/RomanxCMB/mead/plots/"

whattodo = "power_input_zdep,neutrino_test,halofit_comparison,wdependence "
#whattodo = ",halofit_comparison "
#whattodo = "wdependence"



if ("power_input_zdep" in whattodo):
    #test wheather p_lin(k,z) is used or z=0 is extrapolated with internal growth function. Need to turn off neutrinos for this test

    # The easiest way to start a pipeline it from a parameter file.
    ini = Inifile("modules/RomanxCMB/mead/params_gammagrowth.ini")

    # You can modify things in the ini file object after loading.
    # In this case we will switch off some verbose output
    #ini.set("pipeline", "values",  "%(ROMANxCMB_SRC_DIR)s/values.ini")


    # Make the pipeline itself
    pipeline = LikelihoodPipeline(ini)

    # You can modify which parameters you vary now
    pipeline.set_varied("cosmological_parameters", "gamma0", 0.2, 0.9)
    pipeline.set_varied("cosmological_parameters", "gammaa", -0.5, 0.5)

    #pipeline.set_fixed("cosmological_parameters", "h0", 0.72)

    # You can also override these properties if useful
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False



    params_names = pipeline.varied_params

    index = params_names.index("cosmological_parameters--omnuh2")
    print("the index is")
    print(index)
    vector = pipeline.start_vector()
    vector[index] = 0.
    index = params_names.index("cosmological_parameters--gamma0")
    vector[index] = 0.6#0.545
    index = params_names.index("cosmological_parameters--gammaa")
    vector[index] = 0.3

    fig_growth = plt.figure()  # create a figure object
    ax_growth = fig_growth.add_subplot(1, 1, 1)  # create an axes object in the figure

    plt.figure()
    for i in ["T", "F"]:
        ini.set("mead", "power_input_zdep",  i)
        # Make the pipeline itself
        pipeline = LikelihoodPipeline(ini)
        pipeline.set_varied("cosmological_parameters", "gamma0", 0.2, 0.9)
        pipeline.set_varied("cosmological_parameters", "gammaa", -0.5, 0.5)

        data = pipeline.run_parameters(vector)

        #test
        #ax_growth.plot(data["growth_parameters", "a"], data["growth_parameters", "f_z"])


        k = data['matter_power_lin', 'k_h']
        z = data['matter_power_lin', 'z']
        p_lin = data['matter_power_lin', 'p_k'] 
        Pk_lin_resamplez = interp1d(z, p_lin, axis=0)
        # # Make a plot for this value
        # z_sample = 0.
        # plt.loglog(k, Pk_resamplez(z_sample), "--", color="grey")
        # z_sample = 1.
        # plt.loglog(k, Pk_resamplez(z_sample), "--", color="grey")
        # z_sample = 4.
        # plt.loglog(k, Pk_resamplez(z_sample), "--", color="grey")
        k = data['matter_power_nl', 'k_h']
        z = data['matter_power_nl', 'z']
        p_nl = data['matter_power_nl', 'p_k'] 
        
        Pk_resamplez = interp1d(z, p_nl, axis=0)
        # Make a plot for this value
        # z_sample = 0.
        # plt.loglog(k, Pk_resamplez(z_sample), label=str(i)+ " with z= "+str(z_sample))
        # z_sample = 1.
        # plt.loglog(k, Pk_resamplez(z_sample), label=str(i)+ " with z= "+str(z_sample))
        # z_sample = 4.
        # plt.loglog(k, Pk_resamplez(z_sample), label=str(i)+ " with z= "+str(z_sample))

        z_sample = 0.
        plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez(z_sample),label=str(i)+ " with z= "+str(z_sample) )
        z_sample = 1.
        plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez(z_sample),label=str(i)+ " with z= "+str(z_sample) )
        z_sample = 4.
        plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez(z_sample),label=str(i)+ " with z= "+str(z_sample) )

        print("Done ", i)
    plt.title("power_input_zdep")
    plt.xscale("log")
    plt.xlim(10.**(-3), 100.)
    plt.xlabel("k/h")
    plt.ylim(0.7,1.3)
    plt.ylabel("P / P_lin")

    # Save our plot.
    plt.legend()
    plt.savefig(output_folder+"power_spectrum_input_full_or_z0_RAW.pdf")
    #note: for modified growth we get the wrong normalization for the nl power spectrum, indep of whether we use all z or
    #z= 0 only, thats not a good sign
    #it does work with wCDM


if ("neutrino_test" in whattodo):
    #test wheather p_lin(k,z) is used or z=0 is extrapolated with internal growth function. Need to turn off neutrinos for this test

    # The easiest way to start a pipeline it from a parameter file.
    ini = Inifile("modules/RomanxCMB/mead/params.ini")

    # You can modify things in the ini file object after loading.
    # In this case we will switch off some verbose output
    #ini.set("pipeline", "values",  "%(ROMANxCMB_SRC_DIR)s/values.ini")


    # Make the pipeline itself
    pipeline = LikelihoodPipeline(ini)

    # You can modify which parameters you vary now
    #pipeline.set_fixed("cosmological_parameters", "h0", 0.72)

    # You can also override these properties if useful
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False



    params_names = pipeline.varied_params

    index = params_names.index("cosmological_parameters--omnuh2")
    print("the index is")
    print(index)
    vector = pipeline.start_vector()

    plt.figure()

    vector[index] = 0.
    data = pipeline.run_parameters(vector)
    k = data['matter_power_lin', 'k_h']
    z = data['matter_power_lin', 'z']
    p_lin = data['matter_power_lin', 'p_k'] 
    Pk_lin_resamplez = interp1d(z, p_lin, axis=0)

    for i in [0.,0.0003, 0.0006155, 0.0012]:
        # Make the pipeline itself
        vector[index] = i

        data = pipeline.run_parameters(vector)


        # k = data['matter_power_lin', 'k_h']
        # z = data['matter_power_lin', 'z']
        # p_lin = data['matter_power_lin', 'p_k'] 
        # Pk_lin_resamplez = interp1d(z, p_lin, axis=0)

        k = data['matter_power_nl', 'k_h']
        z = data['matter_power_nl', 'z']
        p_nl = data['matter_power_nl', 'p_k'] 
        
        Pk_resamplez = interp1d(z, p_nl, axis=0)

        z_sample = 0.
        plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez(z_sample),label=str(i)+ " with z= "+str(z_sample) )
        z_sample = 1.
        plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez(z_sample),label=str(i)+ " with z= "+str(z_sample) )
        z_sample = 4.
        plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez(z_sample),label=str(i)+ " with z= "+str(z_sample) )

        print("Done ", i)
    plt.title("neutrino effects")
    plt.xscale("log")
    plt.xlim(10.**(-3), 100.)
    plt.xlabel("k/h")
    plt.ylim(0.7,1.3)
    plt.ylabel("P / P_lin (lin without neutrinos)")

    # Save our plot.
    plt.legend()
    plt.savefig(output_folder+"power_spectrum_neutrino_effectsRAW.pdf")


if ("wdependence" in whattodo):
    #test wheather p_lin(k,z) is used or z=0 is extrapolated with internal growth function. Need to turn off neutrinos for this test

    # The easiest way to start a pipeline it from a parameter file.
    ini = Inifile("modules/RomanxCMB/mead/params.ini")

    # You can modify things in the ini file object after loading.
    # In this case we will switch off some verbose output
    #ini.set("pipeline", "values",  "%(ROMANxCMB_SRC_DIR)s/values.ini")


    # Make the pipeline itself
    pipeline = LikelihoodPipeline(ini)

    # You can modify which parameters you vary now
    #pipeline.set_fixed("cosmological_parameters", "h0", 0.72)

    # You can also override these properties if useful
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False



    params_names = pipeline.varied_params

    index = params_names.index("cosmological_parameters--w")
    index2 = params_names.index("cosmological_parameters--w")

    vector = pipeline.start_vector()

    plt.figure()

    # vector[index] = 0.
    # data = pipeline.run_parameters(vector)
    # k = data['matter_power_lin', 'k_h']
    # z = data['matter_power_lin', 'z']
    # p_lin = data['matter_power_lin', 'p_k'] 
    # Pk_lin_resamplez = interp1d(z, p_lin, axis=0)

    for i in [(-1,0.), (-1.2,0), (-0.8,0), (-1, -0.3), (-0.8, 0.2)]:
        # Make the pipeline itself
        vector[index] = i[0]
        vector[index2] = i[1]


        data = pipeline.run_parameters(vector)


        k = data['matter_power_lin', 'k_h']
        z = data['matter_power_lin', 'z']
        p_lin = data['matter_power_lin', 'p_k'] 
        Pk_lin_resamplez = interp1d(z, p_lin, axis=0)

        k = data['matter_power_nl', 'k_h']
        z = data['matter_power_nl', 'z']
        p_nl = data['matter_power_nl', 'p_k'] 
        
        Pk_resamplez = interp1d(z, p_nl, axis=0)

        z_sample = 0.
        plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez(z_sample),label="w={}, wa={} with z= {}".format(i[0], i[1], z_sample) )
        z_sample = 1.
        plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez(z_sample),label="w={}, wa={} with z= {}".format(i[0], i[1], z_sample) )
        z_sample = 4.
        plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez(z_sample),label="w={}, wa={} with z= {}".format(i[0], i[1], z_sample) )

        print("Done ", i)
    plt.title("w,wa dependence")
    plt.xscale("log")
    plt.xlim(10.**(-3), 100.)
    plt.xlabel("k/h")
    plt.ylim(0.7,1.3)
    plt.ylabel("P / P_lin")

    # Save our plot.
    plt.legend()
    plt.savefig(output_folder+"power_spectrum_wdependence_effectsRAW.pdf")



if ("halofit_comparison" in whattodo):
    #compare mead and halofit
    # The easiest way to start a pipeline it from a parameter file.
    

    plt.figure()
    filenames = ["modules/RomanxCMB/mead/params.ini", "modules/RomanxCMB/mead/params.ini", "modules/RomanxCMB/mead/params_halofit.ini"]
    turn_off_mead_feedback = [False, True, False]
    labels = ["mead with feedback", "mead without feedback", "halofit"]
    #for i in range(3):
    ini = Inifile("modules/RomanxCMB/mead/params.ini")
    ini.set("mead", "feedback",  "F")
    pipeline = LikelihoodPipeline(ini)
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False
    params_names = pipeline.varied_params
    # Make the pipeline itself
    pipeline = LikelihoodPipeline(ini)
    vector = pipeline.start_vector()
    data = pipeline.run_parameters(vector)

    k_lin = data['matter_power_lin', 'k_h']
    z = data['matter_power_lin', 'z']
    p_lin = data['matter_power_lin', 'p_k'] 
    Pk_lin_resamplez = interp1d(z, p_lin, axis=0)
    k_mead = data['matter_power_nl', 'k_h']
    z = data['matter_power_nl', 'z']
    p_nl = data['matter_power_nl', 'p_k'] 
    Pk_mead_resamplez = interp1d(z, p_nl, axis=0)

    #with feedback
    ini.set("mead", "feedback",  "T")
    pipeline = LikelihoodPipeline(ini)
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False

    params_names = pipeline.varied_params
    # Make the pipeline itself
    pipeline = LikelihoodPipeline(ini)
    vector = pipeline.start_vector()
    data = pipeline.run_parameters(vector)
    k_mead_feedback = data['matter_power_nl', 'k_h']
    z = data['matter_power_nl', 'z']
    p_nl = data['matter_power_nl', 'p_k'] 
    Pk_mead_feedback_resamplez = interp1d(z, p_nl, axis=0)

    #halofit
    ini = Inifile("modules/RomanxCMB/mead/params_halofit.ini")
    ini.set("mead", "feedback",  "F")
    pipeline = LikelihoodPipeline(ini)
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False
    params_names = pipeline.varied_params
    # Make the pipeline itself
    pipeline = LikelihoodPipeline(ini)
    vector = pipeline.start_vector()
    data = pipeline.run_parameters(vector)

    k_halofit = data['matter_power_nl', 'k_h']
    z = data['matter_power_nl', 'z']
    p_nl = data['matter_power_nl', 'p_k'] 
    Pk_halofit_resamplez = interp1d(z, p_nl, axis=0)


    z_sample = 0.
    plt.plot(k_halofit,Pk_halofit_resamplez(z_sample) /np.interp(k_halofit, k_mead, Pk_mead_resamplez(z_sample)),color="black", label="halofit at z= "+str(z_sample) )
    plt.plot(k_mead_feedback,Pk_mead_feedback_resamplez(z_sample) /np.interp(k_mead_feedback, k_mead, Pk_mead_resamplez(z_sample)),color="orange",label="mead with feedback at z= "+str(z_sample) )
    
    z_sample = 1.
    plt.plot(k_halofit,Pk_halofit_resamplez(z_sample) /np.interp(k_halofit, k_mead, Pk_mead_resamplez(z_sample)),"--", color="black",label="halofit at z= "+str(z_sample) )
    plt.plot(k_mead_feedback,Pk_mead_feedback_resamplez(z_sample) /np.interp(k_mead_feedback, k_mead, Pk_mead_resamplez(z_sample)),"--", color="orange",label="mead with feedback at z= "+str(z_sample) )
    
    z_sample = 4.
    plt.plot(k_halofit,Pk_halofit_resamplez(z_sample) /np.interp(k_halofit, k_mead, Pk_mead_resamplez(z_sample)),":",color="black",label="halofit at z= "+str(z_sample) )
    plt.plot(k_mead_feedback,Pk_mead_feedback_resamplez(z_sample) /np.interp(k_mead_feedback, k_mead, Pk_mead_resamplez(z_sample)),":", color="orange",label="mead with feedback at z= "+str(z_sample) )
    
    plt.title("fiducial cosmology - halofit vs mead")
    plt.xscale("log")
    plt.xlim(10.**(-3), 100.)
    plt.xlabel("k/h")
    plt.ylim(0.7,1.3)
    plt.ylabel("$P(k) / P_{mead} (k)$")
    plt.axvline(10)
    plt.text(11, 0.72, "extrapolated")

    # Save our plot.
    plt.legend()
    plt.savefig(output_folder+"power_spectrum_halofit_comparison_RAW.pdf")


    

