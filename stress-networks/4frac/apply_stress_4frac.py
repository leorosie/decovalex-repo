#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 12:04:03 2020

@author: rosie
"""

import os, sys
from pydfnworks import * 
import numpy as np


def compute_normal_stress_projection(sigma_mat, i, normals):

        # Magnitude of normal stress
        sigma_mag = sigma_mat[0][0]*(normals[i][0])**2 + \
                sigma_mat[1][1]*(normals[i][1])**2 + \
                sigma_mat[2][2]*(normals[i][2])**2 + \
                2*(sigma_mat[0][1]*normals[i][0]*normals[i][1] + \
                sigma_mat[1][2]*normals[i][1]*normals[i][2] + \
                sigma_mat[0][2]*normals[i][0]*normals[i][2])

        T_1 = sigma_mat[0][0]*normals[i][0] + \
                sigma_mat[0][1]*normals[i][1] + \
                sigma_mat[0][2]*normals[i][2]

        T_2 = sigma_mat[1][0]*normals[i][0] + \
                sigma_mat[1][1]*normals[i][1] + \
                sigma_mat[1][2]*normals[i][2]

        T_3 = sigma_mat[2][0]*normals[i][0] + \
                sigma_mat[2][1]*normals[i][1] + \
                sigma_mat[2][2]*normals[i][2]

        stress_sqr = (T_1)**2 + (T_2)**2 + (T_3)**2
        tau =  np.sqrt(max(0, stress_sqr - (sigma_mag)**2))

        return sigma_mag, stress_sqr, tau

src_path = os.getcwd()


import sys

case_name = sys.argv[1]
sigma_x = float(sys.argv[2])
sigma_y = float(sys.argv[3])
sigma_z = float(sys.argv[4])
print(f"case name : {case_name}")
print(f"sigma x : {sigma_x}")
print(f"sigma x : {sigma_y}")
print(f"sigma x : {sigma_z}")

jobname = f"{src_path}/{case_name}"

DFN = DFNWORKS(jobname)


DFN.params['domainSize']['value'] = [1000.0, 1000.0, 1000.0]
DFN.params['domainSizeIncrease']['value'] = [100, 100, 100]
DFN.params['h']['value'] = 3
DFN.params['stopCondition']['value'] = 1 #0 for nPoly
DFN.params['nPoly']['value'] = 40 #for reference ignored if p32
DFN.params['boundaryFaces']['value'] = [0,0,1,1,0,0]
#DFN.params['visualizationMode']['value'] = True
DFN.params['seed']['value'] = 200 #seed for random generator 0 seeds off clock
DFN.params['ignoreBoundaryFaces']['value'] = False
DFN.params['boundaryFaces']['value'] = [1,1,0,0,0,0]
DFN.params['disableFram']['value'] = True 

# DFN.add_fracture_family(shape="ell",
#                         distribution="tpl",
#                         probability = 0.219, #unnecessary unless stop condition = 1
#                         p32 = 0.0099,
#                         beta_distribution = 1,
#                         beta = 0,
#                         number_of_points = 12,
#                         kappa=9.41,
#                         theta=1.65*180./3.14,
#                         phi=4.78*180./3.14,
#                         alpha=2.4,
#                         min_radius=30.0,
#                         max_radius=564.0,
#                         hy_variable='transmissivity',
#                         hy_function='correlated',
#                         hy_params={
#                             "alpha": 2.2e-9, 
#                             "beta": 0.8,
#                         })


# DFN.add_fracture_family(shape="ell",
#                         distribution="tpl",
#                         probability = 0.250,
#                         p32 = 0.0113,
#                         beta_distribution = 1,
#                         beta = 0,
#                         number_of_points = 12,
#                         kappa=8.3,
#                         theta=1.57*180./3.14,
#                         phi=3.14*180./3.14,
#                         alpha=2.4,
#                         min_radius=30.0,
#                         max_radius=564.0,
#                         hy_variable='transmissivity',
#                         hy_function='correlated',
#                         hy_params={
#                             "alpha": 2.2e-9, 
#                             "beta": 0.8,
#                         })

# DFN.add_fracture_family(shape="ell",
#                         distribution="tpl",
#                         probability = 0.531,
#                         p32 = 0.0240,
#                         beta_distribution = 1,
#                         beta = 0,
#                         number_of_points = 12,
#                         kappa=5.7,
#                         theta=2.95*180./3.14,
#                         phi=2.62*180./3.14,
#                         alpha=2.4,
#                         min_radius=30.0,
#                         max_radius=564.0,
#                         hy_variable='transmissivity',
#                         hy_function='correlated',
#                         hy_params={
#                             "alpha": 2.2e-9, 
#                             "beta": 0.8,
#                         })

DFN.add_user_fract(shape='ell',
                    radii=600,
                    translation=[-400, 0, 400],
                    normal_vector=[30, 15, 60],
                    number_of_vertices=5,
                    aperture=1.0e-3)

DFN.add_user_fract(shape='ell',
                    radii=1000,
                    translation=[0, 0, 0],
                    normal_vector=[95, 5, 5],
                    # normal_vector=[1, 0, 0],
                    number_of_vertices=5,
                    aperture=1.0e-3)

DFN.add_user_fract(shape='ell',
                    radii=600,
                    aspect_ratio=1,
                    translation=[400, 0, 200],
                    normal_vector=[30, 15, 60],
                    number_of_vertices=5,
                    aperture=1.0e-3)

DFN.add_user_fract(shape='ell',
                    radii=600,
                    aspect_ratio=1,
                    translation=[400, 0, -400],
                    normal_vector=[30, 15, 60],
                    number_of_vertices=5,
                    aperture=5.0e-5)

DFN.make_working_directory(delete=True)
DFN.print_domain_parameters()
DFN.check_input()
DFN.create_network()
DFN.dump_hydraulic_values()
DFN.visual_mode = True
#DFN.output_report() 

#DFN.mesh_network()
#exit()
# create values of the stress tensor
s1 = sigma_x
s2 = sigma_y
s3 = sigma_z
sigma_mat = x = [[s1, 0, 0], [0, s2, 0], [0, 0, s3]]
DFN.dump_hydraulic_values()
# DFN.add_variable_to_mesh("init_aper", "aperture.dat", "reduced_mesh.inp",
#                          "stress.inp")
# modify apertures basedon the stress field
# add final apertures to mesh
# DFN.add_variable_to_mesh("init_aper", "aperture.dat", "full_mesh.inp",
#                          "stress.inp")

k_n0 = 50 * 1e9 # initial normal stiffness (GPa/m)
k_s = 10 * 1e9 # Shear Stiffness GPa/m)
phi = 10 # dialation angle (degrees)
phi_rad = np.radians(phi)
u_r = 0.005 # Residual shear displacement (m)

for ifrac in range(DFN.num_frac):
    a0 = DFN.aperture[ifrac] # initial aperture
    a_r = DFN.aperture[ifrac]* 0.1 # residual aperture 

    v_n_max = a0 - a_r
    normal_stress, normal_stress_sqr, tau = compute_normal_stress_projection(sigma_mat, 
                                                                        ifrac, DFN.normal_vectors)

#    print(f"\n* normal stress: {normal_stress}")
    v_n = normal_stress * v_n_max / (k_n0 * v_n_max + normal_stress)

#    print(f"\n* shear stress: {tau}")
    u_s = tau / k_s
    part1 = min( abs(u_s), u_r)
    part2 =  np.tan(phi_rad)
    v_s = - part1 * part2

    DFN.aperture[ifrac] = a0 - v_n  - v_s
    #print(ifrac + 1, a0, DFN.aperture[ifrac], v_n, v_s)
    # delta = a0 - DFN.aperture[ifrac]
    # print(delta,"\n")
    # if delta < 0:
    #     print("error!")

DFN.perm = DFN.aperture**2/12
DFN.dump_hydraulic_values("stress")
DFN.dump_hydraulic_values(f"../hydraulic_values/{case_name}")

# DFN.add_variable_to_mesh("stress_aper", "aperture.dat", "stress.inp")
# DFN.lagrit2pflotran() 

# mat_perm = 1e-16
# mat_por = 0.1
# cell_size = 20
# DFN.mapdfn_ecpm(mat_perm, mat_por, cell_size)

# DFN.dfnFlow_file = f"{src_path}/cpm_pflotran.in"
# DFN.local_dfnFlow_file = f"cpm_pflotran.in"
# DFN.pflotran()

# # DFN.dfnFlow_file = f"{src_path}/cpm_transport.in"
# # DFN.local_dfnFlow_file = f"cpm_transport.in"
# # DFN.pflotran()
