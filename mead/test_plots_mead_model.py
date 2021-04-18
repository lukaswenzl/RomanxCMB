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

whattodo = "power_input_zdep,neutrino_test,halofit_comparison,wdependence,growth_test,growth_gamma,mead_with_modified_growth,sampling_accuracy"
#whattodo = ",halofit_comparison "
#whattodo = "wdependence"
#whattodo = "power_input_zdep,growth_test,growth_gamma,mead_with_modified_growth"
#whattodo = "power_input_zdep"#halofit_comparison"
whattodo = "mead_with_modified_growth"
whattodo = "sampling_accuracy"
whattodo = "sample_only_above_cutoff"
#whattodo = "halofit_comparison"




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
    vector[index] = -0.05

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
    plt.xlim(10.**(-5), 100.)
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
    plt.xlim(10.**(-5), 100.)
    plt.xlabel("k/h")
    plt.ylim(0.7,1.3)
    plt.ylabel("$P(k) / P_{mead} (k)$")
    plt.axvline(10)
    plt.text(11, 0.72, "extrapolated")

    # Save our plot.
    plt.legend()
    plt.savefig(output_folder+"power_spectrum_halofit_comparison_RAW.pdf")




if ("growth_test" in whattodo):

    fig_growth = plt.figure()  # create a figure object
    ax_growth = fig_growth.add_subplot(1, 1, 1)  # create an axes object in the figure

    value_w = -0.9
    plt.figure()

    i="cosmosis growth module"
    ini = Inifile("modules/RomanxCMB/mead/params.ini")
    pipeline = LikelihoodPipeline(ini)
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False

    params_names = pipeline.varied_params
    vector = pipeline.start_vector()
    index = params_names.index("cosmological_parameters--w")
    vector[index] = value_w

    #ini.set("growth_gamma", "gamma_parametrization",  "F")

    # Make the pipeline itself
    pipeline = LikelihoodPipeline(ini)

    data = pipeline.run_parameters(vector)
    k = data['matter_power_lin', 'k_h']
    z = data['matter_power_lin', 'z']
    p_lin = data['matter_power_lin', 'p_k'] 
    Pk_lin_resamplez_baseline = interp1d(z, p_lin, axis=0)
    k = data['matter_power_nl', 'k_h']
    z = data['matter_power_nl', 'z']
    p_nl = data['matter_power_nl', 'p_k'] 
    Pk_resamplez = interp1d(z, p_nl, axis=0)
    z_sample = 0.
    plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez_baseline(z_sample),label=str(i)+ " with z= "+str(z_sample) )
    z_sample = 1.
    plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez_baseline(z_sample),label=str(i)+ " with z= "+str(z_sample) )
    z_sample = 4.
    plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez_baseline(z_sample),label=str(i)+ " with z= "+str(z_sample) )

    ax_growth.plot(data["growth_parameters", "a"], data["growth_parameters", "d_z"], label="cosmosis module")

    ## my own growth module
    i = "my growth module"

    ini = Inifile("modules/RomanxCMB/mead/params_gammagrowth.ini")
    pipeline = LikelihoodPipeline(ini)
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False

    params_names = pipeline.varied_params
    vector = pipeline.start_vector()
    index = params_names.index("cosmological_parameters--w")
    vector[index] = value_w

    ini.set("growth_gamma", "gamma_parametrization",  "F")

    # Make the pipeline itself
    pipeline = LikelihoodPipeline(ini)

    data = pipeline.run_parameters(vector)
    # k = data['matter_power_lin', 'k_h']
    # z = data['matter_power_lin', 'z']
    # p_lin = data['matter_power_lin', 'p_k'] 
    # Pk_lin_resamplez_baseline = interp1d(z, p_lin, axis=0)
    k = data['matter_power_nl', 'k_h']
    z = data['matter_power_nl', 'z']
    p_nl = data['matter_power_nl', 'p_k'] 
    Pk_resamplez = interp1d(z, p_nl, axis=0)
    z_sample = 0.
    plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez_baseline(z_sample),label=str(i)+ " with z= "+str(z_sample) )
    z_sample = 1.
    plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez_baseline(z_sample),label=str(i)+ " with z= "+str(z_sample) )
    z_sample = 4.
    plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez_baseline(z_sample),label=str(i)+ " with z= "+str(z_sample) )

    ax_growth.plot(data["growth_parameters", "a"], data["growth_parameters", "d_z"], label="my module")

    ## halofit 
    i ="halofit"
    ini = Inifile("modules/RomanxCMB/mead/params_halofit.ini")
    pipeline = LikelihoodPipeline(ini)
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False

    params_names = pipeline.varied_params
    vector = pipeline.start_vector()
    index = params_names.index("cosmological_parameters--w")
    vector[index] = value_w#-0.8

    #ini.set("growth_gamma", "gamma_parametrization",  "F")

    # Make the pipeline itself
    pipeline = LikelihoodPipeline(ini)

    data = pipeline.run_parameters(vector)
    # k = data['matter_power_lin', 'k_h']
    # z = data['matter_power_lin', 'z']
    # p_lin = data['matter_power_lin', 'p_k'] 
    # Pk_lin_resamplez_baseline = interp1d(z, p_lin, axis=0)
    k = data['matter_power_nl', 'k_h']
    z = data['matter_power_nl', 'z']
    p_nl = data['matter_power_nl', 'p_k'] 
    Pk_resamplez = interp1d(z, p_nl, axis=0)
    z_sample = 0.
    plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez_baseline(z_sample),label=str(i)+ " with z= "+str(z_sample) )
    z_sample = 1.
    plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez_baseline(z_sample),label=str(i)+ " with z= "+str(z_sample) )
    z_sample = 4.
    plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez_baseline(z_sample),label=str(i)+ " with z= "+str(z_sample) )

    



    plt.title("growth test mead vs halofit (see large scales)")
    plt.xscale("log")
    plt.xlim(10.**(-3), 100.)
    plt.xlabel("k/h")
    plt.ylim(0.7,1.3)
    plt.ylabel("P / P_lin")

    # Save our plot.
    plt.legend()
    plt.savefig(output_folder+"power_spectrum_growth_test_RAW.pdf")
    
    ax_growth.legend()
    fig_growth.savefig(output_folder+"growth_cosmosis_vs_mymodule.pdf")


if ("growth_gamma" in whattodo):

    fig_growth = plt.figure()  # create a figure object
    ax_growth = fig_growth.add_subplot(1, 1, 1)  # create an axes object in the figure

    value_gamma = 0.55 + 0.1
    plt.figure()

    ## my own growth module
    i = "my growth module"

    ini = Inifile("modules/RomanxCMB/mead/params_gammagrowth.ini")
    pipeline = LikelihoodPipeline(ini)
    pipeline.set_varied("cosmological_parameters", "gamma0", 0.2, 0.9)
    pipeline.set_varied("cosmological_parameters", "gammaa", -0.5, 0.5)
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False

    params_names = pipeline.varied_params
    vector = pipeline.start_vector()
    index = params_names.index("cosmological_parameters--gamma0")
    vector[index] = value_gamma
    index = params_names.index("cosmological_parameters--omnuh2")
    vector[index] = 0.00005

    ini.set("growth_gamma", "gamma_parametrization",  "T")

    # Make the pipeline itself
    pipeline = LikelihoodPipeline(ini)
    pipeline.set_varied("cosmological_parameters", "gamma0", 0.2, 0.9)
    pipeline.set_varied("cosmological_parameters", "gammaa", -0.5, 0.5)

    data = pipeline.run_parameters(vector)
    k = data['matter_power_lin', 'k_h']
    z = data['matter_power_lin', 'z']
    p_lin = data['matter_power_lin', 'p_k'] 
    Pk_lin_resamplez_baseline = interp1d(z, p_lin, axis=0)
    k = data['matter_power_nl', 'k_h']
    z = data['matter_power_nl', 'z']
    p_nl = data['matter_power_nl', 'p_k'] 
    Pk_resamplez = interp1d(z, p_nl, axis=0)
    z_sample = 0.
    plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez_baseline(z_sample),label=str(i)+ " with z= "+str(z_sample) )
    z_sample = 1.
    plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez_baseline(z_sample),label=str(i)+ " with z= "+str(z_sample) )
    z_sample = 4.
    plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez_baseline(z_sample),label=str(i)+ " with z= "+str(z_sample) )

    ax_growth.plot(data["growth_parameters", "a"], data["growth_parameters", "d_z"], label="my module")

    ## halofit 
    i ="halofit"
    ini = Inifile("modules/RomanxCMB/mead/params_halofit_gammagrowth.ini")
    pipeline = LikelihoodPipeline(ini)
    pipeline.set_varied("cosmological_parameters", "gamma0", 0.2, 0.9)
    pipeline.set_varied("cosmological_parameters", "gammaa", -0.5, 0.5)
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False

    params_names = pipeline.varied_params
    vector = pipeline.start_vector()
    index = params_names.index("cosmological_parameters--gamma0")
    vector[index] = value_gamma#-0.8
    index = params_names.index("cosmological_parameters--omnuh2")
    vector[index] = 0.00005

    ini.set("growth_gamma", "gamma_parametrization",  "T")

    # Make the pipeline itself
    pipeline = LikelihoodPipeline(ini)
    pipeline.set_varied("cosmological_parameters", "gamma0", 0.2, 0.9)
    pipeline.set_varied("cosmological_parameters", "gammaa", -0.5, 0.5)

    data = pipeline.run_parameters(vector)
    k = data['matter_power_lin', 'k_h']
    z = data['matter_power_lin', 'z']
    p_lin = data['matter_power_lin', 'p_k'] 
    Pk_lin_resamplez_baseline = interp1d(z, p_lin, axis=0)
    k = data['matter_power_nl', 'k_h']
    z = data['matter_power_nl', 'z']
    p_nl = data['matter_power_nl', 'p_k'] 
    Pk_resamplez = interp1d(z, p_nl, axis=0)
    z_sample = 0.
    plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez_baseline(z_sample),label=str(i)+ " with z= "+str(z_sample) )
    z_sample = 1.
    plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez_baseline(z_sample),label=str(i)+ " with z= "+str(z_sample) )
    z_sample = 4.
    plt.plot(k,Pk_resamplez(z_sample) /Pk_lin_resamplez_baseline(z_sample),label=str(i)+ " with z= "+str(z_sample) )

    



    plt.title("growth test mead vs halofit (see large scales)")
    plt.xscale("log")
    plt.xlim(10.**(-3), 100.)
    plt.xlabel("k/h")
    plt.ylim(0.7,1.3)
    plt.ylabel("P / P_lin")

    # Save our plot.
    plt.legend()
    plt.savefig(output_folder+"power_spectrum_growth_gamma_RAW.pdf")
    
    #ax_growth.legend()
    #fig_growth.savefig(output_folder+"growth2_cosmosis_vs_mymodule.pdf")

    


if ("mead_with_modified_growth" in whattodo):

    fig_z1 = plt.figure()  # create a figure object
    ax_z1 = fig_z1.add_subplot(1, 1, 1)  # create an axes object in the figure

    value_gamma = 0.55 + 0.1
    plt.figure()

    for value_gamma in [0.45,0.475,0.50,0.525,0.55,0.575, 0.6,0.625, 0.65]:
        ## cosmosis growth piped into mead

        ini = Inifile("modules/RomanxCMB/mead/params_gammagrowth.ini")
        ini.set("growth_gamma", "gamma_parametrization",  "T")
        ini.set("mead", "use_cosmosis_growth",  "T")

        pipeline = LikelihoodPipeline(ini)
        pipeline.set_varied("cosmological_parameters", "gamma0", 0.2, 0.9)
        pipeline.set_varied("cosmological_parameters", "gammaa", -0.5, 0.5)
        pipeline.quiet = True
        pipeline.debug = False
        pipeline.timing = False

        params_names = pipeline.varied_params
        vector = pipeline.start_vector()
        index = params_names.index("cosmological_parameters--gamma0")
        vector[index] = value_gamma
        #index = params_names.index("cosmological_parameters--omnuh2")
        #vector[index] = 0.00005
        data = pipeline.run_parameters(vector)
        k = data['matter_power_nl', 'k_h']
        z = data['matter_power_nl', 'z']
        p_nl = data['matter_power_nl', 'p_k'] 
        Pk_resamplez_piped = interp1d(z, p_nl, axis=0)

        #internal growth from mead
        ini = Inifile("modules/RomanxCMB/mead/params_gammagrowth.ini")
        ini.set("growth_gamma", "gamma_parametrization",  "T")
        ini.set("mead", "use_cosmosis_growth",  "F")

        pipeline = LikelihoodPipeline(ini)
        pipeline.set_varied("cosmological_parameters", "gamma0", 0.2, 0.9)
        pipeline.set_varied("cosmological_parameters", "gammaa", -0.5, 0.5)
        pipeline.quiet = True
        pipeline.debug = False
        pipeline.timing = False

        params_names = pipeline.varied_params
        vector = pipeline.start_vector()
        index = params_names.index("cosmological_parameters--gamma0")
        vector[index] = value_gamma
        #index = params_names.index("cosmological_parameters--omnuh2")
        #vector[index] = 0.00005
        data = pipeline.run_parameters(vector)
        k = data['matter_power_nl', 'k_h']
        z = data['matter_power_nl', 'z']
        p_nl = data['matter_power_nl', 'p_k'] 
        Pk_resamplez_mead = interp1d(z, p_nl, axis=0)


        z_sample = 0.
        plt.plot(k,Pk_resamplez_piped(z_sample) /Pk_resamplez_mead(z_sample),label="$\gamma_0 = {}, z= {}$".format(value_gamma,z_sample) )
        z_sample = 1.
        ax_z1.plot(k,Pk_resamplez_piped(z_sample) /Pk_resamplez_mead(z_sample),"--",label="$\gamma_0 = {}, z= {}$".format(value_gamma,z_sample) )
        # z_sample = 2.
        # plt.plot(k,Pk_resamplez_piped(z_sample) /Pk_resamplez_mead(z_sample),":",label="$\gamma_0 = {}, z= {}$".format(value_gamma,z_sample)  )
        # z_sample = 3.
        # plt.plot(k,Pk_resamplez_piped(z_sample) /Pk_resamplez_mead(z_sample),":",label="$\gamma_0 = {}, z= {}$".format(value_gamma,z_sample)  )
        # z_sample = 4.
        # plt.plot(k,Pk_resamplez_piped(z_sample) /Pk_resamplez_mead(z_sample),":",label="$\gamma_0 = {}, z= {}$".format(value_gamma,z_sample) )


    


    plt.title("Adapting equations (21) and (22) for modified growth")
    plt.xscale("log")
    plt.xlim(10.**(-3), 100.)
    plt.xlabel("k/h")
    plt.ylim(0.8,1.2)
    plt.ylabel("$P_{nl} (eq 21,22 modified) / P_{nl} (original)$")

    # Save our plot.
    plt.legend()
    plt.savefig(output_folder+"power_spectrum_mead_with_modified_growth_RAW.pdf")
    

    ax_z1.set_title("Adapting equations (21) and (22) for modified growth")
    ax_z1.set_xscale("log")
    ax_z1.set_xlim(10.**(-3), 100.)
    ax_z1.set_xlabel("k/h")
    ax_z1.set_ylim(0.8,1.2)
    ax_z1.set_ylabel("$P_{nl} (eq 21,22 modified) / P_{nl} (original)$")
    ax_z1.legend()
    fig_z1.savefig(output_folder+"power_spectrum_mead_with_modified_growth_z1_RAW.pdf")

    

if ("sampling_accuracy" in whattodo):
    #compare different samplings

    plt.figure()

    ini = Inifile("modules/RomanxCMB/mead/params.ini")
    ini.set("mead", "nk",  "450")
    ini.set("mead", "na", "250")

    pipeline = LikelihoodPipeline(ini)
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False

    params_names = pipeline.varied_params
    vector = pipeline.start_vector()
    #index = params_names.index("cosmological_parameters--omnuh2")
    #vector[index] = 0.00005
    data = pipeline.run_parameters(vector)

    k_highsampling = data['matter_power_nl', 'k_h']
    z = data['matter_power_nl', 'z']
    p_nl = data['matter_power_nl', 'p_k'] 
    Pk_highsampling = interp1d(z, p_nl, axis=0)

    for i, fraction in enumerate([1,0.75,1./2, 1./4.]):
        ## cosmosis growth piped into mead

        ini = Inifile("modules/RomanxCMB/mead/params.ini")
        ini.set("mead", "nk",  "{}".format(int(200*fraction)))
        ini.set("mead", "na", "{}".format(int(100*fraction)))

        pipeline = LikelihoodPipeline(ini)
        pipeline.quiet = True
        pipeline.debug = False
        pipeline.timing = False

        params_names = pipeline.varied_params
        vector = pipeline.start_vector()
        #index = params_names.index("cosmological_parameters--omnuh2")
        #vector[index] = 0.00005
        data = pipeline.run_parameters(vector)

        k = data['matter_power_nl', 'k_h']
        z = data['matter_power_nl', 'z']
        p_nl = data['matter_power_nl', 'p_k'] 
        Pk = interp1d(z, p_nl, axis=0)

        z_sample = 0.
        plt.plot(k_highsampling,np.interp(k_highsampling,k,Pk(z_sample)) /Pk_highsampling(z_sample) -i/20.,label="fraction = {}, offset = ".format(fraction,-i/20.) )
        z_sample = 1.
        plt.plot(k_highsampling,np.interp(k_highsampling,k,Pk(z_sample)) /Pk_highsampling(z_sample)-i/20.,"--")
        z_sample = 5.
        plt.plot(k_highsampling,np.interp(k_highsampling,k,Pk(z_sample)) /Pk_highsampling(z_sample)-i/20.,":")
        z_sample = 10.
        plt.plot(k_highsampling,np.interp(k_highsampling,k,Pk(z_sample)) /Pk_highsampling(z_sample)-i/20.,":")
        

    


    plt.title("Testing sampling accuracy")
    plt.xscale("log")
    plt.xlim(10.**(-5), 100.)
    plt.xlabel("k/h")
    plt.ylim(0.8,1.2)
    plt.ylabel("$P_{nl} / P_{nl} (high sampling)$")

    # Save our plot.
    plt.legend()
    
    plt.savefig(output_folder+"power_spectrum_sampling_accuracy_RAW.pdf")

if ("sample_only_above_cutoff" in whattodo):
    #compare different samplings

    plt.figure()

    ini = Inifile("modules/RomanxCMB/mead/params.ini")
    ini.set("mead", "optimize_nl_samples",  "F")

    pipeline = LikelihoodPipeline(ini)
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False

    params_names = pipeline.varied_params
    vector = pipeline.start_vector()
    #index = params_names.index("cosmological_parameters--omnuh2")
    #vector[index] = 0.00005
    data = pipeline.run_parameters(vector)

    k_fullsampling = data['matter_power_nl', 'k_h']
    z = data['matter_power_nl', 'z']
    p_nl = data['matter_power_nl', 'p_k'] 
    Pk_fullsampling = interp1d(z, p_nl, axis=0)


    ini = Inifile("modules/RomanxCMB/mead/params.ini")
    ini.set("mead", "optimize_nl_samples",  "T")

    pipeline = LikelihoodPipeline(ini)
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False

    params_names = pipeline.varied_params
    vector = pipeline.start_vector()
    #index = params_names.index("cosmological_parameters--omnuh2")
    #vector[index] = 0.00005
    data = pipeline.run_parameters(vector)

    k_cutoffsampling = data['matter_power_lin', 'k_h']
    z = data['matter_power_lin', 'z']
    p_lin = data['matter_power_lin', 'p_k'] 
    Pk_lin = interp1d(z, p_lin, axis=0)

    k_cutoffsampling = data['matter_power_nl', 'k_h']
    z = data['matter_power_nl', 'z']
    p_nl = data['matter_power_nl', 'p_k'] 
    Pk_cutoffsampling = interp1d(z, p_nl, axis=0)

    

    z_sample = 0.
    plt.plot(k_fullsampling, Pk_cutoffsampling(z_sample) /Pk_lin(z_sample) ,label="z = {}".format(z_sample) )
    z_sample = 1.
    plt.plot(k_fullsampling,Pk_cutoffsampling(z_sample) /Pk_lin(z_sample) ,"--",label="z = {}".format(z_sample))
    z_sample = 5.
    plt.plot(k_fullsampling,Pk_cutoffsampling(z_sample) /Pk_lin(z_sample) ,":",label="z = {}".format(z_sample))
    z_sample = 10.
    plt.plot(k_fullsampling,Pk_cutoffsampling(z_sample) /Pk_lin(z_sample) ,":",label="z = {}".format(z_sample))
        

    


    plt.title("Testing cutoff for non-linear sampling")
    plt.xscale("log")
    plt.xlim(10.**(-5), 100.)
    plt.xlabel("k/h")
    plt.ylim(0.9,1.1)
    plt.ylabel("$P_{nl} / P_{lin}$")

    # Save our plot.
    plt.legend()
    
    plt.savefig(output_folder+"power_spectrum_sample_only_above_cutoff_RAW.pdf")