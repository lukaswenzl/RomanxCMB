#written for v0.1
#run with 
#python modules/WFIRSTxCMB/script_gaussian_vary_tracer_combinations.py

import os

n_processors = "15"
output_folder = "script_gaussian_vary_tracer_combinations/"
import os
if not os.path.exists(output_folder):
    os.makedirs(output_folder)


#6x2pt
os.system("mpirun -n "+n_processors+" cosmosis modules/WFIRSTxCMB/params_gaussian_covariance.ini -p output.filename='"+output_folder+"6x2pt.txt'")

#3x2pt
os.system("mpirun -n "+n_processors+" cosmosis modules/WFIRSTxCMB/params_gaussian_covariance_only3x2pt.ini -p output.filename='"+output_folder+"3x2pt.txt'")

#shear only
os.system("mpirun -n "+n_processors+" cosmosis modules/WFIRSTxCMB/params_gaussian_covariance_onlyShear.ini -p output.filename='"+output_folder+"shear'.txt")

#galaxy only
os.system("mpirun -n "+n_processors+" cosmosis modules/WFIRSTxCMB/params_gaussian_covariance_onlyClustering.ini -p output.filename='"+output_folder+"galaxy'.txt")

#cmb only
os.system("mpirun -n "+n_processors+" cosmosis modules/WFIRSTxCMB/params_gaussian_covariance_onlyCMBlensing.ini -p output.filename='"+output_folder+"cmblensing'.txt")
