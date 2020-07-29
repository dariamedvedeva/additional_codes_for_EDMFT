import numpy as np
import parameters
import iteration_cycle
import scipy.fftpack as ft
import os
import matplotlib.pyplot as plt

global sqrt32, sqrt3
sqrt32 = np.sqrt(3.)/2.
sqrt3 = np.sqrt(3.)

def func_k(filename):
    # read file
    D = np.loadtxt(filename)
    kx = D[:,0]
    ky = D[:,1]
    V  = D[:,2]
    return kx, ky, V
    
def func_only_real_part(filename):
    # read file
    D = np.loadtxt(filename)
    freq = D[:,0]
    real_part = D[:,1]
    return freq, real_part
    
def func_complex_function(filename):
    # read file
    D = np.loadtxt(filename)
    freq = D[:,0]
    func = D[:,1] + 1j * D[:,2]
    return freq, func

def interaction_GMKG(Nk, lattice_type, param):
    kpoints = np.array([[[x, y] for x in range(Nk)] for y in range(Nk)]).reshape((Nk * Nk, 2))
    kpoints_coordinates = np.zeros((Nk*Nk, 2), dtype = np.float)
    GM_int = np.zeros(Nk * Nk, dtype = np.float)
    MK_int = np.zeros(Nk * Nk, dtype = np.float)
    KG_int = np.zeros(Nk * Nk, dtype = np.float)
    
    b1, b2, b3 = iteration_cycle.get_b_vectors(lattice_type)
    kstep_x, kstep_y = iteration_cycle.get_kstep(lattice_type, Nk)
   
    point_x = []
    point_y = []
    interaction = []
    point_x_MK = []
    point_y_MK = []
    int_MK = []
    
    for k in range(Nk * Nk):
        vec_ = [0.0, 0.0]
        for i in range(2):
            if (lattice_type == "square"):
                kpoints_coordinates[k][i] = kpoints[k][0] * b1[i]/ (Nk - 1) + kpoints[k][1] * b2[i] / (Nk - 1)
            elif (lattice_type == "triangular"):
                kpoints_coordinates[k][i] = kpoints[k][0] * b1[i]/ (Nk - 1) + kpoints[k][1] * b2[i] / (Nk - 1)
    
    x_boundary = 0.0
    y_boundary = 0.0
    i = 0
    for k in kpoints_coordinates:
        # TRIANG
        if lattice_type == 'triangular':
            k_x = k[0]
            k_y = k[1]
            
            x_coord = k_x
            y_coord = k_y
            
            if (x_coord == 0.0) and (y_coord >= 0.0) and (y_coord < 2.*np.pi/sqrt3):
                GM_int[i]  = 2. * param[0] * ( np.cos(x_coord) + np.cos(0.5 * x_coord + sqrt32 * y_coord) + np.cos(0.5 * x_coord - sqrt32 * y_coord) )
                GM_int[i] += 2. * param[1] * ( np.cos(sqrt3 * y_coord) + np.cos(3./2. * x_coord + sqrt32 * y_coord) + np.cos(3./2. * x_coord - sqrt32 * y_coord) )
                GM_int[i] += 2. * param[2] * ( np.cos(2. * x_coord) + np.cos(x_coord + sqrt3 * y_coord) + np.cos(x_coord - sqrt3 * y_coord) )
                interaction.append(GM_int[i])
                point_x.append(x_coord)
                point_y.append(y_coord)
                i += 1

            elif (abs(y_coord - 2.* np.pi/sqrt3) <= 0.001) and (x_coord >= 0.0) and (x_coord <= 2.*np.pi/3.):
                MK_int[i]  = 2. * param[0] * ( np.cos(x_coord) + np.cos(0.5 * x_coord + sqrt32 * y_coord) + np.cos(0.5 * x_coord - sqrt32 * y_coord) )
                MK_int[i] += 2. * param[1] * ( np.cos(sqrt3 * y_coord) + np.cos(3./2. * x_coord + sqrt32 * y_coord) + np.cos(3./2. * x_coord - sqrt32 * y_coord) )
                MK_int[i] += 2. * param[2] * ( np.cos(2. * x_coord) + np.cos(x_coord + sqrt3 * y_coord) + np.cos(x_coord - sqrt3 * y_coord) )
                int_MK.append(MK_int[i])
                point_x_MK.append(x_coord)
                point_y_MK.append(y_coord)
                i += 1

    for j in range(len(point_x_MK)):
        if j !=0 : # it shouldn't be here, but don't touch it (Nb2S)
            point_x.append(point_x_MK[-1-j])
            point_y.append(point_y_MK[-1-j])
            interaction.append(int_MK[-1-j])


    x_coord =  2.*np.pi / 3.
    while x_coord >= 0.0:
        x_coord -= kstep_x
        if (x_coord < 0.0):
            break
        y_coord = x_coord * sqrt3

        KG_int[i]  = 2. * param[0] * ( np.cos(x_coord) + np.cos(0.5 * x_coord + sqrt32 * y_coord) + np.cos(0.5 * x_coord - sqrt32 * y_coord) )
        KG_int[i] += 2. * param[1] * ( np.cos(sqrt3 * y_coord) + np.cos(3./2. * x_coord + sqrt32 * y_coord) + np.cos(3./2. * x_coord - sqrt32 * y_coord) )
        KG_int[i] += 2. * param[2] * ( np.cos(2. * x_coord) + np.cos(x_coord + sqrt3 * y_coord) + np.cos(x_coord - sqrt3 * y_coord) )
        point_x.append(x_coord)
        point_y.append(y_coord)
        interaction.append(KG_int[i])
        i += 1
  
    return point_x, point_y, interaction

     
def read_real_function(filename):
    # read file
    D = np.loadtxt(filename)
    argument = D[:,0]
    function = D[:,1]
    return argument, function
    
def read_imag_part_function(filename):
    # read file
    D = np.loadtxt(filename)
    argument = D[:,0]
    function = 1j * D[:,2]
    return argument, function

"""
    main part
"""

parameters.set_model_parameters()
lattice_type, beta, U, hartree_shift, Nk, num_of_neighbours, t, V, mu, particle_hole_symm, sweeps, time_limit,  delta_mix, lambda_mix, number_of_iterations, start_from_it = \
parameters.get_model_parameters()
#Nk = 727
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#frequencies, impurity_fuction = iteration_cycle.read_real_function(filename)   # for X
filename = 'DOS.dat'
frequencies, spectral_function_A = func_only_real_part(filename)               # for DOS
Im_GF = -np.pi * spectral_function_A
# - Kramers–Kronig relations from Im GF -> Re GF
Re_GF = ft.hilbert(Im_GF) / np.pi
GF = Re_GF + 1j * Im_GF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
kx, ky, int = interaction_GMKG(Nk, lattice_type, t)
np.savetxt("interection_dispersions.dat", np.column_stack((kx, ky, int)))
disp = open("interection_dispersions_value.dat", "w")
for i in range(len(int)):
    disp.write(str(int[i]))
    disp.write("\n")
disp.close()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
filename_real_hybridization = "Delta_r.dat"
freq_r, hybridization = func_complex_function(filename_real_hybridization)

result_file = open("spectral.dat", "w")
for j in range(len(frequencies)):
    point_number = 1
    for i in range(len(int)):
        local_function_val =  1.0 / (1.0 / GF[j] + hybridization[j] - int[i])
#        print(local_function)
        result_file.write(str(point_number))
        result_file.write("\t")
        result_file.write(str(frequencies[j]))
        result_file.write("\t")
        result_file.write(str(local_function_val.imag))
        result_file.write("\t")
        result_file.write(str(kx[i]))
        result_file.write("\t")
        result_file.write(str(ky[i]))
        result_file.write("\n")
        point_number += 1
    result_file.write("\n")
result_file.close()
    
print("int len = {}".format(len(int)))
print("lambda_charge len = {}".format(len(hybridization)))
print("Impurity_fuction len = {}".format(len(GF)))
os.system("gnuplot plot_bands.gnuplot")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
interaction_full_lattice, points = iteration_cycle.interaction_dispersion('_tk_remove_', Nk, lattice_type, t)
local_function_full       = np.zeros(GF.shape, dtype=np.complex128)

for i in range(len(interaction_full_lattice)):
    local_function_full += 1.0 / (1.0 / GF + hybridization - interaction_full_lattice[i])
local_function_full /= Nk**2

np.savetxt("G_real_local.dat", np.column_stack((frequencies, local_function_full.real, local_function_full.imag)))

