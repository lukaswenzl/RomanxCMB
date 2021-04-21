from __future__ import print_function
from builtins import str
from cosmosis.runtime.config import Inifile
from cosmosis.runtime.pipeline import LikelihoodPipeline
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d

output_folder = "modules/RomanxCMB/la_model/plots/"

whattodo = "powerspectra,amplitude"
whattodo = "amplitude"

if ("powerspectra" in whattodo):
    # The easiest way to start a pipeline it from a parameter file.
    ini = Inifile("modules/RomanxCMB/la_model/params_IA.ini")
    #ini = Inifile("modules/RomanxCMB/params.ini")


    # You can modify things in the ini file object after loading.
    # In this case we will switch off some verbose output
    #ini.set("pipeline", "values",  "%(ROMANxCMB_SRC_DIR)s/values.ini")

    # Make the pipeline itself
    pipeline = LikelihoodPipeline(ini)


    # You can also override these properties if useful
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False

    #pipeline.set_varied("cosmological_parameters", "omega_m", 0.2, 0.4)

    # Let's look through different values of omega_m
    # and get a Galaxy Galaxy-Lensing spectrum for each of them


    params_names = pipeline.varied_params
    #print(params_names)
    index = params_names.index("intrinsic_alignment_parameters--a")
    print("the index is")
    print(index)
    for A in [5, 5.95, 7]:
    #for A in [0.1, 1., 5.95]:


        # In this method of running the pipeline we
        # pass it a value for each of the parameters 
        # we have told it to vary.
        # We could check what these are by looking at
        #pipeline.varied_params

        vector = pipeline.start_vector()
        vector[index] = A

        data = pipeline.run_parameters(vector)

        # data is a DataBlock - can get things out of it as in any
        # cosmosis module:

        k = data['matter_power_nl', 'k_h']
        z = data['matter_power_nl', 'z']
        p = data['matter_power_nl', 'p_k'] 
        Pk_nl = interp1d(z, p, axis=0)

        k = data['matter_intrinsic_power', 'k_h']
        z = data['matter_intrinsic_power', 'z']
        p = data['matter_intrinsic_power', 'p_k'] 
        Pk_matter_intrinsic = interp1d(z, p, axis=0)

        k = data['intrinsic_power', 'k_h']
        z = data['intrinsic_power', 'z']
        p = data['intrinsic_power', 'p_k'] 
        Pk_intrinsic = interp1d(z, p, axis=0)

        z_sample = 0.5
        plt.plot(k, Pk_nl(z_sample),"--", label="nl, A = {}, z= {}".format(A, z_sample))
        plt.plot(k, Pk_intrinsic(z_sample), label="intrinsic, A = {}, z= {}".format(A, z_sample))


        print("Done ", A)

    # Save our plot.
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.savefig(output_folder+"test_scriptingRAW.pdf")

if ("amplitude" in whattodo):
    # The easiest way to start a pipeline it from a parameter file.
    ini = Inifile("modules/RomanxCMB/la_model/params_IA.ini")
    ini_original = Inifile("modules/RomanxCMB/la_model/params_oldIA.ini")


    # You can modify things in the ini file object after loading.
    # In this case we will switch off some verbose output
    #ini.set("pipeline", "values",  "%(ROMANxCMB_SRC_DIR)s/values.ini")

    # Make the pipeline itself
    pipeline = LikelihoodPipeline(ini)
    pipeline_original = LikelihoodPipeline(ini_original)
    # You can also override these properties if useful
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False
    pipeline_original.quiet = True
    pipeline_original.debug = False
    pipeline_original.timing = False

    #pipeline.set_varied("cosmological_parameters", "omega_m", 0.2, 0.4)

    # Let's look through different values of omega_m
    # and get a Galaxy Galaxy-Lensing spectrum for each of them


    params_names = pipeline.varied_params
    #print(params_names)
    index = params_names.index("intrinsic_alignment_parameters--a")
    print("the index is")
    print(index)
    for A in [1., 5.95]:
    #for A in [0.1, 1., 5.95]:


        # In this method of running the pipeline we
        # pass it a value for each of the parameters 
        # we have told it to vary.
        # We could check what these are by looking at
        #pipeline.varied_params

        vector = pipeline.start_vector()
        vector[index] = A

        data = pipeline.run_parameters(vector)

        vector = pipeline_original.start_vector()
        vector[index] = A

        data_original = pipeline_original.run_parameters(vector)

        # data is a DataBlock - can get things out of it as in any
        # cosmosis module:

        k = data['matter_power_nl', 'k_h']
        z = data['matter_power_nl', 'z']
        p_nl = data['matter_power_nl', 'p_k'] 
        #Pk_nl = interp1d(z, p, axis=0)
        k = data['intrinsic_power', 'k_h']
        z = data['intrinsic_power', 'z']
        p_intrinsic = data['intrinsic_power', 'p_k'] 
        #Pk_intrinsic = interp1d(z, p, axis=0)

        #large k
        plt.plot(z, p_intrinsic[:,0 ]/p_nl[:,0], label="eNLA, A = {}".format(A))

        k = data_original['matter_power_nl', 'k_h']
        z = data_original['matter_power_nl', 'z']
        p_nl = data_original['matter_power_nl', 'p_k'] 
        #Pk_nl = interp1d(z, p, axis=0)
        k = data_original['intrinsic_power', 'k_h']
        z = data_original['intrinsic_power', 'z']
        p_intrinsic = data_original['intrinsic_power', 'p_k'] 
        #Pk_intrinsic = interp1d(z, p, axis=0)

        #large scale
        plt.plot(z, p_intrinsic[:,0 ]/p_nl[:,0],"--", label="simple NLA, A = {}".format(A))


        print("Done ", A)

    # Make the pipeline itself
    pipeline = LikelihoodPipeline(ini)
    pipeline_original = LikelihoodPipeline(ini_original)
    # You can also override these properties if useful
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False
    pipeline_original.quiet = True
    pipeline_original.debug = False
    pipeline_original.timing = False
    params_names = pipeline.varied_params
    #print(params_names)
    index = params_names.index("intrinsic_alignment_parameters--beta")
    print("the index is")
    print(index)

    for beta in [0., 2.]:
    #for A in [0.1, 1., 5.95]:


        # In this method of running the pipeline we
        # pass it a value for each of the parameters 
        # we have told it to vary.
        # We could check what these are by looking at
        #pipeline.varied_params

        vector = pipeline.start_vector()
        vector[index] = beta

        data = pipeline.run_parameters(vector)

        vector = pipeline_original.start_vector()
        vector[index] = beta

        data_original = pipeline_original.run_parameters(vector)

        # data is a DataBlock - can get things out of it as in any
        # cosmosis module:

        k = data['matter_power_nl', 'k_h']
        z = data['matter_power_nl', 'z']
        p_nl = data['matter_power_nl', 'p_k'] 
        #Pk_nl = interp1d(z, p, axis=0)
        k = data['intrinsic_power', 'k_h']
        z = data['intrinsic_power', 'z']
        p_intrinsic = data['intrinsic_power', 'p_k'] 
        #Pk_intrinsic = interp1d(z, p, axis=0)

        #large k
        plt.plot(z, p_intrinsic[:,0 ]/p_nl[:,0], label="eNLA, beta = {}".format(beta))


        print("Done ", beta)


    # Make the pipeline itself
    pipeline = LikelihoodPipeline(ini)
    pipeline_original = LikelihoodPipeline(ini_original)
    # You can also override these properties if useful
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False
    pipeline_original.quiet = True
    pipeline_original.debug = False
    pipeline_original.timing = False
    params_names = pipeline.varied_params
    #print(params_names)
    index = params_names.index("intrinsic_alignment_parameters--eta")
    print("the index is")
    print(index)

    for eta in [0., 1]:
        #for A in [0.1, 1., 5.95]:


        # In this method of running the pipeline we
        # pass it a value for each of the parameters 
        # we have told it to vary.
        # We could check what these are by looking at
        #pipeline.varied_params

        vector = pipeline.start_vector()
        vector[index] = eta

        data = pipeline.run_parameters(vector)

        vector = pipeline_original.start_vector()
        vector[index] = eta

        data_original = pipeline_original.run_parameters(vector)

        # data is a DataBlock - can get things out of it as in any
        # cosmosis module:

        k = data['matter_power_nl', 'k_h']
        z = data['matter_power_nl', 'z']
        p_nl = data['matter_power_nl', 'p_k'] 
        #Pk_nl = interp1d(z, p, axis=0)
        k = data['intrinsic_power', 'k_h']
        z = data['intrinsic_power', 'z']
        p_intrinsic = data['intrinsic_power', 'p_k'] 
        #Pk_intrinsic = interp1d(z, p, axis=0)

        #large k
        plt.plot(z, p_intrinsic[:,0 ]/p_nl[:,0], label="eNLA, eta = {}".format(eta))


        print("Done ", eta)
    
    # Make the pipeline itself
    pipeline = LikelihoodPipeline(ini)
    pipeline_original = LikelihoodPipeline(ini_original)
    # You can also override these properties if useful
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False
    pipeline_original.quiet = True
    pipeline_original.debug = False
    pipeline_original.timing = False

    params_names = pipeline.varied_params
    #print(params_names)
    index = params_names.index("intrinsic_alignment_parameters--eta_highz")
    print("the index is")
    print(index)

    for eta in [-0.5, 0.5]:
        #for A in [0.1, 1., 5.95]:


        # In this method of running the pipeline we
        # pass it a value for each of the parameters 
        # we have told it to vary.
        # We could check what these are by looking at
        #pipeline.varied_params

        vector = pipeline.start_vector()
        vector[index] = eta

        data = pipeline.run_parameters(vector)

        vector = pipeline_original.start_vector()
        vector[index] = eta

        data_original = pipeline_original.run_parameters(vector)

        # data is a DataBlock - can get things out of it as in any
        # cosmosis module:

        k = data['matter_power_nl', 'k_h']
        z = data['matter_power_nl', 'z']
        p_nl = data['matter_power_nl', 'p_k'] 
        #Pk_nl = interp1d(z, p, axis=0)
        k = data['intrinsic_power', 'k_h']
        z = data['intrinsic_power', 'z']
        p_intrinsic = data['intrinsic_power', 'p_k'] 
        #Pk_intrinsic = interp1d(z, p, axis=0)

        #large k
        plt.plot(z, p_intrinsic[:,0 ]/p_nl[:,0], label="eNLA, eta_highz = {}".format(eta))


        print("Done ", eta)

    # Make the pipeline itself
    pipeline = LikelihoodPipeline(ini)
    pipeline_original = LikelihoodPipeline(ini_original)
    # You can also override these properties if useful
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False
    pipeline_original.quiet = True
    pipeline_original.debug = False
    pipeline_original.timing = False
    
    params_names = pipeline_original.varied_params
    #print(params_names)
    index = params_names.index("intrinsic_alignment_parameters--alpha")
    print("the index is")
    print(index)

    index2 = params_names.index("intrinsic_alignment_parameters--a")


    for alpha in [-1., 1]:
    #for A in [0.1, 1., 5.95]:


        # In this method of running the pipeline we
        # pass it a value for each of the parameters 
        # we have told it to vary.
        # We could check what these are by looking at
        #pipeline.varied_params


        vector = pipeline_original.start_vector()
        vector[index] = alpha
        A = 1.
        vector[index2] = A


        data_original = pipeline_original.run_parameters(vector)

        k = data_original['matter_power_nl', 'k_h']
        z = data_original['matter_power_nl', 'z']
        p_nl = data_original['matter_power_nl', 'p_k'] 
        #Pk_nl = interp1d(z, p, axis=0)
        k = data_original['intrinsic_power', 'k_h']
        z = data_original['intrinsic_power', 'z']
        p_intrinsic = data_original['intrinsic_power', 'p_k'] 
        #Pk_intrinsic = interp1d(z, p, axis=0)

        #large scale
        plt.plot(z, p_intrinsic[:,0 ]/p_nl[:,0],"--", label="simple NLA, alpha = {}, A={}".format(alpha, A))


        print("Done ", alpha)

    # Save our plot.
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("z")
    plt.ylabel("P_II /P (large scale) = A^2")
    plt.legend()
    plt.savefig(output_folder+"plot_intrinsic_amplitudeRAW.pdf")
