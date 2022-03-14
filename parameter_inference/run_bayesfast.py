#code to run bayesfast
#bayesfast needs to be installed seperately

import numpy as np
from scipy.stats import norm
from scipy.linalg import sqrtm
from cosmosis.runtime.config import Inifile
from cosmosis.runtime.pipeline import LikelihoodPipeline
import dill
import bayesfast as bf
import time
import sys
import os
# from distributed import Client, LocalCluster
# os.chdir('/global/cfs/cdirs/des/shivamp/nl_cosmosis/cosmosis/y3-3x2pt/code')
os.environ['OMP_NUM_THREADS'] = '2'
from multiprocess.pool import Pool
# n_core=2
# from ray.util.multiprocessing import Pool
pool = Pool(processes=14)

def parse_param_files(ini,init_v):

    values_file = ini.get('pipeline', 'values')
    prior_file = ini.get('pipeline', 'priors')

    values_ini = Inifile(values_file)
    prior_ini = Inifile(prior_file)

    pkeys, vkeys = [], []
    init_values, para_range = [], []
    _prior_mu, _prior_sig = [], []
    print(values_ini.keys())
    for k in values_ini.keys():

        param = '{}--{}'.format(k[0][0], k[0][1])
        line = values_ini.get(k[0][0], k[0][1])
        values = line.split()
        if len(values) == 1:
            continue
        elif len(values) == 3:
            vkeys.append(param)
            if init_v is None:
                init_values.append(float(values[1]))
            else:
                init_values.append(init_v[param])
            para_range.append([float(values[0]), float(values[2])])
        else:
            raise("Cannot parse values file, too many quanties specified for {}".format(
                k[0][0], k[0][1]))

    for k in prior_ini.keys():

        param = '{}--{}'.format(k[0][0], k[0][1])
        line = prior_ini.get(k[0][0], k[0][1])
        values = line.split()

        if (len(values) == 1):
            continue
        elif len(values) == 3:
            if values[0] != 'gaussian':
                raise(
                    'Non-gaussian/tophat prior specified. This is not currently supported by BayesFast')

            pkeys.append(param)
            _prior_mu.append(float(values[1]))
            _prior_sig.append(float(values[2]))
        else:
            raise("Cannot parse priors file, too many quanties specified for {}".format(
                k[0][0], k[0][1]))

    init_mu = np.array(init_values)
    print('init_mu is ' , init_mu)
    para_range = np.array(para_range)

    _prior_mu = np.asarray(_prior_mu)
    _prior_sig = np.asarray(_prior_sig)

    pkeys = np.array(pkeys)
    vkeys = np.array(vkeys)

    # only keep parameters in priors that are actually specified in the values file
    pidx = np.in1d(pkeys, vkeys)
    pkeys = pkeys[pidx]
    _prior_mu = _prior_mu[pidx]
    _prior_sig = _prior_sig[pidx]

    init_sig = (para_range[:, 1] - para_range[:, 0]) #/ 1000
    print('parameters are', vkeys)
    idx_dict = dict(zip(vkeys, np.arange(len(vkeys))))

    # nl_params = [p.split(',')[0] for p in ini.get('bayesfast', 'nonlinear-params').split()]
    # _nonlinear_indices = np.array([idx_dict[k] for k in nl_params])
    quadratic_params = [p.split(',')[0] for p in ini.get('bayesfast', 'quadratic-params').split()]
    cubic_params = [p.split(',')[0] for p in ini.get('bayesfast', 'cubic-params').split()]
    _quadratic_indices = [idx_dict[k] for k in quadratic_params]
    _cubic_indices = [idx_dict[k] for k in cubic_params]
    print('quadratic parameters are', quadratic_params)
    print('cubic parameters are', cubic_params)


    useIS = eval(ini.get('bayesfast', 'useIS'))
    n_IS = np.int(ini.get('bayesfast', 'n_IS'))

    n_chain = np.int(ini.get('bayesfast', 'n_chain'))
    n_iter = np.int(ini.get('bayesfast', 'n_iter'))
    n_warmup = np.int(ini.get('bayesfast', 'n_warmup'))
    n_x_0 = np.int(ini.get('bayesfast', 'n_x_0'))

    extend_neutrino = float(ini.get('bayesfast', 'extend_neutrino'))
    extend_IA = float(ini.get('bayesfast', 'extend_IA'))



    # return init_mu, para_range, _prior_mu, _prior_sig, pkeys, vkeys, init_sig,\
    #     _nonlinear_indices, idx_dict, useIS, n_IS, nl_params
    return init_mu, para_range, _prior_mu, _prior_sig, pkeys, vkeys, init_sig, _quadratic_indices, _cubic_indices, idx_dict, useIS, n_IS, n_chain, n_iter, n_warmup, n_x_0, extend_neutrino, extend_IA



def main(ini_string, fname, fnamer0, init_v=None, x_0_cov=None, x_0_fisher=None, init_sig_coeff=100.):
    print("##BF starting main")

    ini = Inifile(ini_string)

    # old_stdout = sys.stdout
    # sys.stdout = open(os.devnull, 'w')

    # init_mu, para_range, _prior_mu, _prior_sig, pkeys, vkeys, init_sig,\
    #     _nonlinear_indices, idx_dict, useIS, n_IS, nl_params = parse_param_files(ini,init_v)
    init_mu, para_range, _prior_mu, _prior_sig, pkeys, vkeys, init_sig, _quadratic_indices, _cubic_indices, idx_dict, useIS, n_IS, n_chain, n_iter, n_warmup, n_x_0, extend_neutrino, extend_IA = parse_param_files(ini,init_v)
    
    print('NL params')
    pipeline = LikelihoodPipeline(ini)
    # sys.stdout = old_stdout

    print("##BF run cosmosis pipeline on start vector")
    start = pipeline.start_vector()
    print(dict(zip(vkeys,start)))
    results = pipeline.run_results(start)

    pkeys = pkeys[np.in1d(np.array(pkeys), vkeys)]
    _prior_indices = np.array([idx_dict[p] for p in pkeys])
    _flat_indices = np.setdiff1d(np.fromiter(
        idx_dict.values(), dtype=int), _prior_indices)

    if len(_prior_indices)>0:
        _prior_norm = (
            -0.5 * np.sum(np.log(2 * np.pi * _prior_sig**2)) - np.sum(np.log(
                norm.cdf(para_range[_prior_indices, 1], _prior_mu, _prior_sig) -
                norm.cdf(para_range[_prior_indices, 0], _prior_mu, _prior_sig))) -
            np.sum(np.log(para_range[_flat_indices, 1] -
                          para_range[_flat_indices, 0]))
        )
    else:
        print('Only flat priors')
        _prior_norm = (
            -0.5 * np.sum(np.log(2 * np.pi * _prior_sig**2))
            - np.sum(np.log(para_range[_flat_indices, 1] -
                          para_range[_flat_indices, 0])))

    _d = results.block['data_vector', '2pt_data']
    nData = _d.shape[0]

    _invC = results.block['data_vector', '2pt_inverse_covariance']
    _invC_r = np.real(sqrtm(_invC))

    _d_diag = _d @ _invC_r
    _norm = results.block['data_vector', '2pt_norm']

    try:
        extra_params = ini.get('pipeline', 'extra_output')
        if (not hasattr(extra_params, '__iter__')) | (type(extra_params) == str):
            extra_params = extra_params.split()
            ep = np.copy(extra_params)
            extra_params = [e.split('/') for e in extra_params]
        else:
            ep = np.copy(extra_params)
            extra_params = [e.split('/') for e in extra_params]            
    except:
            extra_params = []

    nExtraParams = len(extra_params)
    nData += nExtraParams
    print('Also building surrogate for extra parameters: {}'.format(extra_params))
    
    nParams = len(vkeys)

    if len(_prior_indices)>0:

        def des_prior_f(x): # prior chisq + log prior 
            chi2 = -0.5 * np.sum(((x[_prior_indices] - _prior_mu) / _prior_sig)**2)
            return chi2 + _prior_norm

        def des_prior_j(x):  # prior gradient
            foo = np.zeros((1, nParams))
            foo[0, _prior_indices] = - \
                    (x[_prior_indices] - _prior_mu) / _prior_sig**2
            return foo
        
    else:
        def des_prior_f(x): # prior chisq + log prior chi2 = -0.5 *
            return _prior_norm

        def des_prior_j(x):  # prior gradient
            foo = np.zeros((1, nParams))
            return foo

    def des_2pt_theory(x, _invC_r=_invC_r,ini_string=ini_string):

        # run DES pipeline to get data*invCov , theory *invCov
        try:
            import datetime
            import os
            import sys
            os.environ['OMP_NUM_THREADS'] = '2'
            from cosmosis.runtime.config import Inifile
            from cosmosis.runtime.pipeline import LikelihoodPipeline
            old_stdout = sys.stdout
            sys.stdout = open(os.devnull, 'w')
            ini = Inifile(ini_string)
            pipeline = LikelihoodPipeline(ini)
            sys.stdout = old_stdout
            t0 = datetime.datetime.now()
            res = pipeline.run_results(x)
            t0 = datetime.datetime.now() - t0
            print("# Cosmosis pipeline took (seconds) : ", t0.total_seconds())
            try:
                if(res is None):
                    print("Warning: Pipeline returned None")
                    print(datetime.datetime.now())
                    raise
                rres = res.block['data_vector', '2pt_theory'] @ _invC_r
                try:
                    extra_params = ini.get('pipeline', 'extra_output')
                    if (not hasattr(extra_params, '__iter__')) | (type(extra_params) == str):
                        extra_params = extra_params.split()
    #                    print(extra_params)
                        extra_params = [e.split('/') for e in extra_params]
                except:
                    extra_params = None
                if extra_params is not None:
                    try:
                        extra = [res.block[e[0], e[1]] for e in extra_params]
                    except:
                        extra = [1]
                        print(x)
                        pass
                else:
                    extra = [1]
                return np.hstack([rres, extra])
            except:
                print(x)
                return np.nan * np.ones(nData)
                pass

        except Exception as e:
            raise(e)
            print('Failed to run pipeline! Returning nans')
            return np.nan * np.ones(nData)

    def select_2pt(allout, nData=(nData-nExtraParams)):

        return allout[:nData]

    def select_extra(allout, nExtraParams=nExtraParams):

        return allout[-nExtraParams:]

    def chi2_f(m):  # lhood chisq, uses covariance now
        return np.atleast_1d(-0.5 * np.sum((m - _d_diag)**2) + _norm)

    def chi2_fj(m):  # lood chisq gradient
        return (np.atleast_1d(-0.5 * np.sum((m - _d_diag)**2) + _norm),
                -(m - _d_diag)[np.newaxis])

    def select_2pt_jac(allout, nData=(nData - nExtraParams), nExtraParams=nExtraParams):

        jac = np.diag(np.ones(nData + nExtraParams))
        jac[:,-nExtraParams:] = 0
        
        return allout[:nData], jac[:-nExtraParams,:]

    def des_post_f(like, x):  # like+prior
        return like + des_prior_f(x)

    def des_post_fj(like, x):  # like + prior and prior gradients
        return like + des_prior_f(x), np.concatenate(
            (np.ones((1, 1)), des_prior_j(x)), axis=-1)

    print("##BF model definitions")
    # parameters-> theory model
    m_0 = bf.Module(fun=des_2pt_theory, input_vars='x',
                    output_vars=['allout'])

    m_1 = bf.Module(fun=select_2pt, input_vars='allout',
                    fun_and_jac=select_2pt_jac,
                    output_vars=['m'])

    # theory model -> likelihood
    m_2 = bf.Module(fun=chi2_f, fun_and_jac=chi2_fj,
                    input_vars=['m'], output_vars='like')

    # likelihood and parameters -> log posterior
    m_3 = bf.Module(fun=des_post_f, fun_and_jac=des_post_fj,
                    input_vars=['like', 'x'], output_vars='logp')

    # stack modules to go from params -> log posterior
    d_0 = bf.Density(density_name='logp', module_list=[m_0, m_1, m_2, m_3],
                     input_vars='x', input_dims=nParams, input_scales=para_range,
                     hard_bounds=True)

    print('BF posterior density, cosmosis posterior: {}, {}'.format(
        d_0(start), results.post))
    print('These should match closely!!')

    print('Checking size of polynomial inputs. There are {0} model paramters and {1} data points'.format(
        nParams, nData))

    print("##BF surrogate model definitions")
    # approximate theory model with linear model
    s_0 = bf.modules.PolyModel('linear', input_size=nParams, output_size=nData, input_vars='x',
                               output_vars=['allout'], input_scales=para_range)

    # another linear model
    pc_0 = bf.modules.PolyConfig('linear')

    # nonlinear model - quadratic
    # pc_1 = bf.modules.PolyConfig('quadratic', input_mask=_nonlinear_indices)
    # pc_2 = bf.modules.PolyConfig('cubic-2', input_mask=_nonlinear_indices)
    # pc_1 = bf.modules.PolyConfig('quadratic', input_mask=_quadratic_indices)
    # pc_2 = bf.modules.PolyConfig('cubic-2', input_mask=_cubic_indices+_quadratic_indices)
    pc_1 = bf.modules.PolyConfig('quadratic', input_mask=_cubic_indices+_quadratic_indices)
    if len(_cubic_indices)>0:
        pc_2 = bf.modules.PolyConfig('cubic-2', input_mask=_cubic_indices)
        pc = [pc_0, pc_1, pc_2]
    else:
        pc = [pc_0, pc_1]

    # approximate theory model as linear + quadratic
    s_1 = bf.modules.PolyModel(pc,  input_size=nParams, output_size=nData, input_vars='x',
                               output_vars=['allout'], input_scales=para_range)

    # check if in bounds provided by ranges in values.ini
    def _in_bound(xx, bound):
        xxt = np.atleast_2d(xx).T
        return np.product([np.where(xi > bound[i, 0], True, False) *
                           np.where(xi < bound[i, 1], True, False) for i, xi in
                           enumerate(xxt)], axis=0).astype(bool)

    print("##BF surrogate model training points")
    if x_0_cov is not None:
        #raise ValueError
        x_0_cov = np.loadtxt(x_0_cov)
        print("Using covariance file to sample surrogate training points", x_0_cov)
        x_0 = bf.utils.sobol.multivariate_normal(init_mu, x_0_cov, 10000) # can be outside bounds
    elif x_0_fisher is not None:
        #raise ValueError
        x_0_fisher = np.loadtxt(x_0_fisher)

        cov_matrix = np.linalg.inv(x_0_fisher)

        if(extend_neutrino != 1.):
            print("extending the sigma for neutrinos by factor {}".format(extend_neutrino))
            params_names = pipeline.varied_params
            indexomnuh2 = params_names.index("cosmological_parameters--omnuh2")

            cov_matrix[indexomnuh2, :] *=  extend_neutrino
            cov_matrix[:, indexomnuh2] *=  extend_neutrino

        if(extend_IA != 1.):
            print("extending the sigma for IA beta by factor {}".format(extend_neutrino))
            params_names = pipeline.varied_params
            indexA = params_names.index("intrinsic_alignment_parameters--a")
            indexbeta = params_names.index("intrinsic_alignment_parameters--beta")

            print(np.sqrt(cov_matrix[indexbeta, indexbeta]))

            #cov_matrix[indexA, :] *=  extend_IA
            #cov_matrix[:, indexA] *=  extend_IA

            cov_matrix[indexbeta, :] *=  extend_IA
            cov_matrix[:, indexbeta] *=  extend_IA

            print(np.sqrt(cov_matrix[indexbeta, indexbeta]))

        print("Using fisher file to sample surrogate training points", x_0_fisher)
        x_0 = bf.utils.sobol.multivariate_normal(init_mu, cov_matrix, 50000) # can be outside bounds      
        print(x_0[0])  
        print(x_0[1])  
    else:
        x_0 = bf.utils.sobol.multivariate_normal(init_mu, np.diag((init_sig/init_sig_coeff)**2), n_x_0)
        
        
    # bf.utils.parallel.set_backend(n_core)
    bf.utils.parallel.set_backend(pool)
    
    print("##BF checking that points are inside bounds")
    x_0 = x_0[_in_bound(x_0, para_range)]
    x_0 = x_0[:n_x_0,:]
    print(x_0.shape)

    #enforcing w0+wa < 0 to avoid code crashes
    try:
        params_names = pipeline.varied_params
        indexofw0 = params_names.index("cosmological_parameters--w")
        indexofwa = params_names.index("cosmological_parameters--wa")

        w0pluswa = np.array([i[indexofw0] + i[indexofwa] for i in x_0])
        x_0 = x_0[w0pluswa < 0.]
        print("removed any points with w0+wa >= 0. New shape:")
        print(x_0.shape)
    except ValueError as e:
        print ('Warning: Tried to check w0+wa but parameters not found. Ok if these are not varied')

    #import pdb; pdb.set_trace()
    
    sample_trace_0={"n_chain": n_chain,
            "n_iter": n_iter,
            "n_warmup": n_warmup}
    
    print("##BF OptimizeStep")
    sample = [bf.recipe.SampleStep(s_1, alpha_n=4, alpha_min=0.75, sample_trace=sample_trace_0, x_0=x_0),
              bf.recipe.SampleStep(s_1, alpha_n=4, alpha_min=0.75, sample_trace=sample_trace_0, logp_cutoff=False, reuse_samples=1)
             ]
              
    r_0 = bf.recipe.Recipe(density=d_0, sample=sample)

    print("##BF start run")
    print(time.strftime('%H:%M%p %Z on %b %d, %Y'))
    r_0.run()
    
    print("##BF finish run")
    print(time.strftime('%H:%M%p %Z on %b %d, %Y'))
    
    results = r_0.get()
    samples = results.samples
    logq = np.atleast_2d(results.logq).T
    surrogate = r_0._density._surrogate_list[0]
    sigma8 = np.array([surrogate(samples[i,:])[0][-1:] for i in range(len(samples))])

    vkeys = list(vkeys)
    vkeys.extend(['COSMOLOGICAL_PARAMETERS--SIGMA_8','post'])
    chain = np.hstack([samples, sigma8, logq])
    np.savetxt(fname, chain, header=' '.join(vkeys))    
    dill.dump(r_0.get(),open(fnamer0,'wb'))


if __name__ == '__main__':

    print(bf.__file__)

    ini_file = os.path.join(os.environ['RUN_FOLDER'],'params.ini')  

    fname = os.path.join(os.environ['RUN_FOLDER'], 'chain_'+os.environ['RUN_NAME']+'_'+os.environ['DATAFILE']+'.txt')
    fnamer0 = os.path.join(os.environ['RUN_FOLDER'], 'bf_'+os.environ['RUN_NAME']+'_'+os.environ['DATAFILE']+'.dill')

    x_0_fisher = os.path.join(os.environ['RUN_FOLDER'], 'fisher_'+os.environ['RUN_NAME']+'_'+os.environ['DATAFILE']+'.txt')

    main(ini_file, fname, fnamer0, init_sig_coeff=100., x_0_fisher=x_0_fisher)



