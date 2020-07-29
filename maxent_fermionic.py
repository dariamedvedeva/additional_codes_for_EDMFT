# # # # # # # # # # # # # # # # # # # # # # # #
#         Medvedeva D.S. 2018 - 2020          #
# For questions: medvedeva.ds@gmail.com       #
# Skype: daryacooper                          #
# # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
import subprocess
import os.path
import sys
import os
import time
from scipy import integrate
import iteration_cycle

def create_input_file_maxent(beta, data_file_name_for_maxent, num_of_data_points, PARTICLE_HOLE_SYMMETRY, min_w, max_w, max_iters_for_fitting, NORM):
    file_with_parameters = open("in.param", "w")
    # inverse temperature
    file_with_parameters.write("BETA=")
    file_with_parameters.write(str(beta))
    file_with_parameters.write("\n")

    # 0 || 1
    file_with_parameters.write("PARTICLE_HOLE_SYMMETRY=")
    file_with_parameters.write(str(bool(PARTICLE_HOLE_SYMMETRY)))
    file_with_parameters.write("\n")
    
    # num of data points
    file_with_parameters.write("NDAT=")
    if (PARTICLE_HOLE_SYMMETRY == 0):
        file_with_parameters.write(str(num_of_data_points*2))
    if (PARTICLE_HOLE_SYMMETRY == 1):
        file_with_parameters.write(str(num_of_data_points))
    file_with_parameters.write("\n")

    # num of output frequencies
    file_with_parameters.write("NFREQ=")
    file_with_parameters.write(str(1000))
    file_with_parameters.write("\n")

    # G(iw)
    file_with_parameters.write("DATASPACE=frequency")
    file_with_parameters.write("\n")

    #fermionic|bosonic values
    file_with_parameters.write("KERNEL=fermionic")
    file_with_parameters.write("\n")

    # location of data file
    file_with_parameters.write("DATA=")
    file_with_parameters.write(data_file_name_for_maxent)
    file_with_parameters.write("\n")

    # Minimum frequency
    file_with_parameters.write("OMEGA_MIN=")
    file_with_parameters.write(str(min_w))
    file_with_parameters.write("\n")

    # Maximum frequency
    file_with_parameters.write("OMEGA_MAX=")
    file_with_parameters.write(str(max_w))
    file_with_parameters.write("\n")
    
    # Type of frequency grid (default value: Lorentzian)
    file_with_parameters.write("FREQUENCY_GRID=Quadratic")
    file_with_parameters.write("\n")

#    # log_min for log grid (default value: 0.0001)
#    file_with_parameters.write("LOG_MIN=")
#    file_with_parameters.write(str(0.0001))
#    file_with_parameters.write("\n")

    # Default model for entropy (default value: flat) "Gaussian"
#    file_with_parameters.write("DEFAULT_MODEL=\"double Gaussian\"")     # <== For Susceptibility
    file_with_parameters.write("DEFAULT_MODEL=\"Gaussian\"")           # <== For DOS
    file_with_parameters.write("\n")

    # stddev - For Gaussian models
    file_with_parameters.write("SIGMA=0.5")
    file_with_parameters.write("\n")
    

    # shift of a model (default value: 0)
#    file_with_parameters.write("SHIFT=5.5")
#    file_with_parameters.write("SHIFT=1.0")
    file_with_parameters.write("\n")

    # Maximum Iterations for the fitting routine (default value: 1000)
    file_with_parameters.write("MAX_IT=")
    file_with_parameters.write(str(max_iters_for_fitting))
    file_with_parameters.write("\n")

    # Number of alpha samples (default value: 60)
#    file_with_parameters.write("N_ALPHA=100")
#    file_with_parameters.write("\n")
#
#    file_with_parameters.write("ALPHA_MIN=0.005")
#    file_with_parameters.write("\n")
    
#    file_with_parameters.write("ALPHA_MAX=5")
#    file_with_parameters.write("\n")

    # true to print verbose output (default value: false)
    file_with_parameters.write("VERBOSE=0")
    file_with_parameters.write("\n")
    
#    file_with_parameters.write("NORM=0.374355")
    if(NORM > 0.0):
        file_with_parameters.write("NORM=")
        file_with_parameters.write(str(NORM))
        file_with_parameters.write("\n")
    
    file_with_parameters.close()

def construct_data_file_for_maxent(input_file, output_file):
    w, Gloc = iteration_cycle.read_freq_function(input_file)
    Re_Error = np.zeros(Gloc.shape, np.float)
    Im_Error = np.zeros(Gloc.shape, np.float)
    static_error = 0.0001
    for i in range(len(Re_Error)):
        Re_Error[i] = static_error
        Im_Error[i] = static_error*10
    np.savetxt(output_file, np.column_stack((w.imag, Gloc.real, Re_Error, Gloc.imag, Im_Error)))

def run(path_to_maxent, beta, filename_for_maxent, local, number_of_fermionic_freqs, particle_hole_symm, min_w, max_w, max_iterations_for_fitting, NORM):
    
    create_input_file_maxent(beta, filename_for_maxent, number_of_fermionic_freqs, particle_hole_symm, min_w, max_w, max_iterations_for_fitting, NORM)
    construct_data_file_for_maxent('Gw.dat', filename_for_maxent)
    subprocess.call([path_to_maxent, "in.param"])

def read_real_function(filename):
    # read file
    D = np.loadtxt(filename)
    argument = D[:,0]
    function = D[:,1]
    return argument, function
    
#
##############################################
##                                           #
##            TUNE OF EXECUTION              #
##                                           #
##############################################
#
#server          = False
#path_to_maxent  = '/Users/witcher/workspace/CT_HYB_SEGMENT/Maxent/build2/maxent'
#
## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##                                                             #
##    C T - H Y B    S E G M E N T   C A L C U L A T I O N     #
##                                                             #
## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
## run CT-HYB SEGMENT solver
#number_of_fermionic_freqs               = 1024
#number_of_fermionic_freqs_for_fourier   = 512   # Because of noise we cut the tail of Delta (fermionic hybr. function)
## off and make a Fouriet transform into the \tau - space by the first frequencies with smooth data.
#number_of_bosonic_frequencies           = 1024
#number_of_discrete_tau_points           = 4096  # Friedrich - 4096
#
#start_time = time.time()
#
#print (" ")
#print ("\t++++++++++++++++++++")
#print ("\t    MAXENT  ")
#print ("\t++++++++++++++++++++")
#print (" ")
#
#                # - - - - - - - - - - - - - - - - - #
#                #  (1) S E T  P A R A M E T E R S   #
#                # - - - - - - - - - - - - - - - - - #
#
## Maxent for Gw. The script will create a file 'Gw_for_maxent.dat' and after will use it
## to run maxent.
#filename_for_maxent = 'Gw_for_maxent.dat'
## local = True  --> use Gloc.dat (local    Green's function)
## local = False --> use Gw.dat   (impurity Green's function)
#local = False
#min_w = -5.0
#max_w = 5.0
#max_iterations_for_fitting = 1000000
#
#
#                            # - - - - - -  #
#                            #  (2) R U N   #
#                            # - - - - - -  #
#
#NORM = 1.0 # Look at output and put it here. I am still not sure how it works
#maxent.run(path_to_maxent, beta, filename_for_maxent, local, number_of_fermionic_freqs, particle_hole_symm, min_w, max_w, max_iterations_for_fitting, NORM)
#
#
#                # - - - - - - - - - - - - - - - - - -  #
#                #  (3) D O S  &&  O C C U P A T I O N  #
#                # - - - - - - - - - - - - - - - - - -  #
#
## 'in.out.maxspec.dat' - the name of the file is maxent's output and reserved.
#maxent_filename = 'in.out.maxspec.dat'
#if (os.path.exists(maxent_filename)):
#    w, dos = iteration_cycle.read_real_function(maxent_filename)
#else:
#    print("File >> {} << doesn't exist.".format(maxent_filename))
#    sys.exit()
#
## I check if DOS is normalized on 1. (For one orbital)
#integral = np.round(integrate.trapz(dos, w), 4)
#print("Integral dos = {}".format(integral))
#if (integral == np.round(1.0, 4)):
#    print("DOS is normalized")
#else:
#    dos /= integral
#    # after normalization
#    integral2 = np.round(integrate.trapz(dos, w), 4)
#    print("Integral dos = {} after normalization".format(integral2))
#os.system("mv in.out.maxspec.dat in.out.maxspec_Gw.dat ")
#
#
#print ("************ Calculation is finished. ************")
#print("Time = {} min".format(np.round((time.time() - start_time)/60),2))
#print ("************         *(~.~)*          ************")
