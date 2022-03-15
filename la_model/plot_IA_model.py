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
#whattodo = "Cls"
whattodo = "red_fraction"#,Cls,amplitude"
whattodo ="la_model_param_dependence"
whattodo = "la_model_param_dependence,Cl_param_dependence,Cl_old_model_param_dependence"
whattodo = "la_model_param_dependence,Cl_param_dependence"
whattodo = "Cl_components"
whattodo = "IAbetaAdegen"

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
    #plt.xscale("log")
    plt.xlim(0,3.5)
    plt.yscale("log")
    plt.xlabel("z")
    plt.ylabel("$P_II /P$ (k fixed) = $f^2_{red} A^2$")
    plt.legend()
    plt.savefig(output_folder+"plot_intrinsic_amplitudeRAW.pdf")


if ("Cls" in whattodo):
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

    index_beta = params_names.index("intrinsic_alignment_parameters--beta")

    plt.figure()
    #for A in [1., 5.95, 10., 100., 1000]:
    for A in [1., 5.95, 10.]:

    #for A in [0.1, 1., 5.95]:


        # In this method of running the pipeline we
        # pass it a value for each of the parameters 
        # we have told it to vary.
        # We could check what these are by looking at
        #pipeline.varied_params

        vector = pipeline.start_vector()
        vector[index] = A
        vector[index_beta] = 0

        data = pipeline.run_parameters(vector)

        vector = pipeline_original.start_vector()
        vector[index] = A

        data_original = pipeline_original.run_parameters(vector)

        # data is a DataBlock - can get things out of it as in any
        # cosmosis module:

        ell = data['shear_cl', 'ell']
        Cl22 = data['shear_cl', 'bin_2_2'] 
        Cl55 = data['shear_cl', 'bin_5_5'] 

        #large k
        plt.plot(ell, Cl22, label="eNLA, A = {}, bin=2,2".format(A))
        plt.plot(ell, Cl55,"--", label="eNLA, A = {}, bin=5,5".format(A))


        ell = data_original['shear_cl', 'ell']
        Cl22 = data_original['shear_cl', 'bin_2_2'] 
        Cl55 = data_original['shear_cl', 'bin_5_5'] 

        #large k
        # plt.plot(ell, Cl22, label="simple NLA, A = {}, bin=2,2".format(A))
        # plt.plot(ell, Cl55, label="simple NLA, A = {}, bin=5,5".format(A))

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

    for beta in [-4, 0., 2., 6]:
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

        ell = data['shear_cl', 'ell']
        Cl22 = data['shear_cl', 'bin_2_2'] 
        Cl55 = data['shear_cl', 'bin_5_5'] 

        #large k
        plt.plot(ell, Cl22, label="eNLA, beta = {}, bin=2,2".format(beta))
        plt.plot(ell, Cl55,"--", label="eNLA, beta = {}, bin=5,5".format(beta))


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

        ell = data['shear_cl', 'ell']
        Cl22 = data['shear_cl', 'bin_2_2'] 
        Cl55 = data['shear_cl', 'bin_5_5'] 

        #large k
        plt.plot(ell, Cl22, label="eNLA, eta = {}, bin=2,2".format(eta))
        plt.plot(ell, Cl55,"--", label="eNLA, eta = {}, bin=5,5".format(eta))


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

        ell = data['shear_cl', 'ell']
        Cl22 = data['shear_cl', 'bin_2_2'] 
        Cl55 = data['shear_cl', 'bin_5_5'] 

        #large k
        plt.plot(ell, Cl22, label="eNLA, eta_highz = {}, bin=2,2".format(eta))
        plt.plot(ell, Cl55,"--", label="eNLA, eta_highz = {}, bin=5,5".format(eta))


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

    # Save our plot.
    #plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("ell")
    plt.ylabel("$C_l$")
    plt.legend()
    plt.savefig(output_folder+"plot_ClRAW.pdf")



if ("red_fraction" in whattodo):
    # The easiest way to start a pipeline it from a parameter file.
    ini = Inifile("modules/RomanxCMB/la_model/params_IA.ini")


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

    vector = pipeline.start_vector()

    data = pipeline.run_parameters(vector)

    z = data["intrinsic_alignment", "z"]
    fred = data["intrinsic_alignment", "fred"]

    plt.figure()
    plt.plot(z, fred)
    plt.ylim(0,0.2)
    plt.xlim(0, 3.5)
    plt.xlabel("z")
    plt.ylabel("fraction of red galaxies")
    

    #plt.legend()
    plt.savefig(output_folder+"plot_f_red_RAW.pdf")

    plt.figure()

    A = data["intrinsic_alignment", "amplitude"]

    plt.figure()
    plt.plot(z, A)
    #plt.ylim(0,0.2)
    plt.xlim(0, 3.5)
    plt.xlabel("z")
    plt.ylabel("A (mlim, z)")
    

    #plt.legend()
    plt.savefig(output_folder+"plot_amplitude_RAW.pdf")




def run_min_max_of_range(pipeline, section, param):
    params_names = pipeline.varied_params
    idx = params_names.index(section+"--"+param.lower())
    vector = pipeline.start_vector()
    vector_max = pipeline.max_vector()
    vector_min = pipeline.min_vector()

    vector[idx] = vector_max[idx]
    data_max = pipeline.run_parameters(vector)

    vector[idx] = vector_min[idx]
    data_min = pipeline.run_parameters(vector)

    return data_max, data_min


    


if ("la_model_param_dependence" in whattodo):
    # The easiest way to start a pipeline it from a parameter file.
    ini = Inifile("modules/RomanxCMB/la_model/params_IA.ini")


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
    plt.figure()    
    section = "intrinsic_alignment_parameters"
    param = "A"
    vector = pipeline.start_vector()
    data = pipeline.run_parameters(vector)
    z = data["intrinsic_alignment", "z"]
    A = np.abs(data["intrinsic_alignment", "amplitude"])

    for param in ["eta", "eta_highz","A", "beta" ]:
        data_max, data_min = run_min_max_of_range(pipeline, section, param)
        z = data_max["intrinsic_alignment", "z"]
        A_max = np.abs(data_max["intrinsic_alignment", "amplitude"])
        A_min = np.abs(data_min["intrinsic_alignment", "amplitude"])

        high = np.maximum(np.maximum(A_max, A_min), A)
        low = np.minimum(np.minimum(A_max, A_min), A)
        plt.fill_between(z, low, high, alpha = 0.3, label=param)

    plt.plot(z, A,color="black", label="fiducial")

    plt.ylim(0,1)
    plt.xlim(0, 3.5)
    plt.xlabel("z")
    plt.ylabel("|A (mlim, z)|")
    

    plt.legend()
    plt.savefig(output_folder+"plot_A_model_dependence_RAW.pdf")


if ("Cl_param_dependence" in whattodo):
    # The easiest way to start a pipeline it from a parameter file.
    ini = Inifile("modules/RomanxCMB/la_model/params_IA.ini")


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
       
    section = "intrinsic_alignment_parameters"
    vector = pipeline.start_vector()
    data = pipeline.run_parameters(vector)


    #tracer = "shear_cl"
    #bin = "bin_2_2"
    for tracer, bin in [("shear_cl", "bin_2_2"),("shear_cl", "bin_3_2"), ("shear_cl", "bin_6_1"),("galaxy_shear_cl", "bin_2_2"), ("galaxy_shear_cl", "bin_2_4"),("shear_cmbkappa_cl", "bin_2_1")]:
    
        plt.figure() 
        ell = data[tracer, 'ell']
        factor = ell*(ell+1)/(2.*np.pi)
        Cl22 = np.abs(data[tracer, bin]) *factor

        plt.title(tracer+" "+bin)
        #Cl55 = data['shear_cl', 'bin_5_5'] 
        

        for param in ["eta", "eta_highz","A", "beta" ]:
            data_max, data_min = run_min_max_of_range(pipeline, section, param)
            ell = data_min[tracer, 'ell']
            observable_max = np.abs(data_max[tracer, bin]) *factor
            observable_min = np.abs(data_min[tracer, bin]) *factor

            high = np.maximum(np.maximum(observable_max, observable_min), Cl22)
            low = np.minimum(np.minimum(observable_max, observable_min), Cl22)
            plt.fill_between(ell, low, high, alpha = 0.3, label=param)

        plt.plot(ell, Cl22,color="black", label="fiducial")

        plt.yscale("log")
        plt.xscale("log")
        #plt.xlim(0, 3.5)
        plt.xlabel("l")
        plt.ylabel("$l(l+1) C_l / (2\pi )$")
        

        plt.legend()
        plt.savefig(output_folder+"plot_Cl_model_dependence_"+tracer+"_"+bin+"_RAW.pdf")

if ("Cl_components" in whattodo):
    # The easiest way to start a pipeline it from a parameter file.
    ini = Inifile("modules/RomanxCMB/la_model/params_IA.ini")


    # You can modify things in the ini file object after loading.
    # In this case we will switch off some verbose output
    #ini.set("pipeline", "values",  "%(ROMANxCMB_SRC_DIR)s/values.ini")

    # Make the pipeline itself
    pipeline = LikelihoodPipeline(ini)
    # You can also override these properties if useful
    pipeline.quiet = True
    pipeline.debug = False
    pipeline.timing = False

       
    section = "intrinsic_alignment_parameters"
    vector = pipeline.start_vector()
    data = pipeline.run_parameters(vector)


    tracer = "shear_cl_gi"
    tracer = "shear_cl_gi"

    #for comparison simple IA model
    ini_old = Inifile("modules/RomanxCMB/la_model/params_oldIA.ini")
    pipeline_old = LikelihoodPipeline(ini_old)
    vector_old = pipeline_old.start_vector()
    data_old = pipeline_old.run_parameters(vector_old)


    #bin = "bin_2_2"
    for tracer, bin in [("shear_cl", "bin_1_1"),("shear_cl", "bin_2_1"),("shear_cl", "bin_2_2"),("shear_cl", "bin_3_2"), ("shear_cl", "bin_5_5"), ("shear_cl", "bin_5_1"),("shear_cl", "bin_8_1")]:
    
        plt.figure() 
        ell = data[tracer, 'ell']
        factor = ell*(ell+1)/(2.*np.pi)
        

        plt.title(tracer+" "+bin)
        #Cl55 = data['shear_cl', 'bin_5_5'] 

        if (tracer == "shear_cl"):
            suff = ["", "_gi","_ii", "_gg"]
        for i in suff:
            if (i == "_gg"):
                Cl = np.abs(data[tracer, bin]) *factor - np.abs(data[tracer+"_gi", bin]) *factor - np.abs(data[tracer+"_ii", bin]) *factor
            else:
                Cl = np.abs(data[tracer+i, bin]) *factor
            plt.plot(ell, Cl, label="full model "+i)

            #ols IA model
            if (i == "_gg"):
                Cl_old = np.abs(data_old[tracer, bin]) *factor - np.abs(data_old[tracer+"_gi", bin]) *factor - np.abs(data_old[tracer+"_ii", bin]) *factor
            else:
                Cl_old = np.abs(data_old[tracer+i, bin]) *factor
            plt.plot(ell, Cl_old,"--", label="simple IA "+i)

        plt.yscale("log")
        plt.xscale("log")
        #plt.xlim(0, 3.5)
        plt.legend()
        plt.xlabel("l")
        plt.ylabel("$l(l+1) C_l / (2\pi )$")
        

        #plt.legend()
        plt.savefig(output_folder+"plot_Cl_components_"+tracer+"_"+bin+"_RAW.pdf")


if ("Cl_old_model_param_dependence" in whattodo):
    # The easiest way to start a pipeline it from a parameter file.
    ini = Inifile("modules/RomanxCMB/la_model/params_oldIA.ini")


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
       
    section = "intrinsic_alignment_parameters"
    vector = pipeline.start_vector()
    data = pipeline.run_parameters(vector)


    #tracer = "shear_cl"
    #bin = "bin_2_2"
    for tracer, bin in [("shear_cl", "bin_2_2"),("shear_cl", "bin_3_2"), ("shear_cl", "bin_6_1"),("galaxy_shear_cl", "bin_2_2"), ("galaxy_shear_cl", "bin_2_4"), ("shear_cmbkappa_cl", "bin_2_1")]:
    
        plt.figure() 
        ell = data[tracer, 'ell']
        factor = ell*(ell+1)/(2.*np.pi)
        Cl22 = np.abs(data[tracer, bin]) *factor

        plt.title(tracer+" "+bin)
        #Cl55 = data['shear_cl', 'bin_5_5'] 
        

        for param in ["A", "alpha" ]:
            data_max, data_min = run_min_max_of_range(pipeline, section, param)
            ell = data_min[tracer, 'ell']
            observable_max = np.abs(data_max[tracer, bin]) *factor
            observable_min = np.abs(data_min[tracer, bin]) *factor

            high = np.maximum(np.maximum(observable_max, observable_min), Cl22)
            low = np.minimum(np.minimum(observable_max, observable_min), Cl22)
            plt.fill_between(ell, low, high, alpha = 0.3, label=param)

        plt.plot(ell, Cl22,color="black", label="fiducial")

        plt.yscale("log")
        plt.xscale("log")
        #plt.xlim(0, 3.5)
        plt.xlabel("l")
        plt.ylabel("$l(l+1) C_l / (2\pi )$")
        

        plt.legend()
        plt.savefig(output_folder+"plot_Cl_old_model_dependence_"+tracer+"_"+bin+"_RAW.pdf")

    
if ("IAbetaAdegen" in whattodo):
    # The easiest way to start a pipeline it from a parameter file.
    ini = Inifile("modules/RomanxCMB/params.ini")
    #ini_original = Inifile("modules/RomanxCMB/la_model/params_oldIA.ini")


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
    index_A = params_names.index("intrinsic_alignment_parameters--a")
    index_beta = params_names.index("intrinsic_alignment_parameters--beta")

    print("the index is")
    print(index_A)
    for A in [ 1., 2., 4., 5.95, 8., 10.]: #
    #for A in [0.1, 1., 5.95]:


        vector = pipeline.start_vector()
        vector[index_A] = A

        data = pipeline.run_parameters(vector)

        A_data = data["intrinsic_alignment", "amplitude"]
        fred = data["intrinsic_alignment", "fred"]
        z = data["intrinsic_alignment", "z"]

        plt.plot(z, -1.*A_data, label="A = {}".format(A))

    for beta in [-0.1, 0.,0.5, 1., 1.5,2., 3.]: #2., 3., 4., 5., 
    #for A in [0.1, 1., 5.95]:


        vector = pipeline.start_vector()
        vector[index_beta] = beta

        data = pipeline.run_parameters(vector)

        A_data = data["intrinsic_alignment", "amplitude"]
        fred = data["intrinsic_alignment", "fred"]
        z = data["intrinsic_alignment", "z"]

        plt.plot(z, -1.*A_data, "--",label="beta = {}".format(beta))

    plt.xlabel("z")
    plt.ylabel("A")
    plt.legend()
    plt.savefig(output_folder+"IAbetaAdegeneracy_RAW.pdf")

