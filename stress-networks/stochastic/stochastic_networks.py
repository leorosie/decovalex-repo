#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 12:04:03 2020

@author: rosie
"""

import os, sys
from pydfnworks import * 
import numpy as np

src_path = os.getcwd()

import sys

sample_index = int(sys.argv[1])
jobname = f"{src_path}/sample_x{sample_index:02d}"

DFN = DFNWORKS(jobname)


DFN.params['domainSize']['value'] = [1000.0, 1000.0, 1000.0]
DFN.params['domainSizeIncrease']['value'] = [100, 100, 100]
DFN.params['h']['value'] = 3
DFN.params['stopCondition']['value'] = 1 #0 for nPoly
DFN.params['nPoly']['value'] = 40 #for reference ignored if p32
DFN.params['boundaryFaces']['value'] = [0,0,1,1,0,0]
#DFN.params['visualizationMode']['value'] = True
DFN.params['seed']['value'] = 200 * sample_index #seed for random generator 0 seeds off clock
DFN.params['ignoreBoundaryFaces']['value'] = False
DFN.params['boundaryFaces']['value'] = [1,1,0,0,0,0]
DFN.params['disableFram']['value'] = True 

DFN.add_fracture_family(shape="ell",
                        distribution="tpl",
                        probability = 0.219, #unnecessary unless stop condition = 1
                        p32 = 0.0099,
                        beta_distribution = 1,
                        beta = 0,
                        number_of_points = 12,
                        kappa=9.41,
                        theta=1.65*180./3.14,
                        phi=4.78*180./3.14,
                        alpha=2.4,
                        min_radius=30.0,
                        max_radius=564.0,
                        hy_variable='transmissivity',
                        hy_function='correlated',
                        hy_params={
                            "alpha": 2.2e-9, 
                            "beta": 0.8,
                        })


DFN.add_fracture_family(shape="ell",
                        distribution="tpl",
                        probability = 0.250,
                        p32 = 0.0113,
                        beta_distribution = 1,
                        beta = 0,
                        number_of_points = 12,
                        kappa=8.3,
                        theta=1.57*180./3.14,
                        phi=3.14*180./3.14,
                        alpha=2.4,
                        min_radius=30.0,
                        max_radius=564.0,
                        hy_variable='transmissivity',
                        hy_function='correlated',
                        hy_params={
                            "alpha": 2.2e-9, 
                            "beta": 0.8,
                        })

DFN.add_fracture_family(shape="ell",
                        distribution="tpl",
                        probability = 0.531,
                        p32 = 0.0240,
                        beta_distribution = 1,
                        beta = 0,
                        number_of_points = 12,
                        kappa=5.7,
                        theta=2.95*180./3.14,
                        phi=2.62*180./3.14,
                        alpha=2.4,
                        min_radius=30.0,
                        max_radius=564.0,
                        hy_variable='transmissivity',
                        hy_function='correlated',
                        hy_params={
                            "alpha": 2.2e-9, 
                            "beta": 0.8,
                        })

# DFN.add_user_fract(shape='ell',
#                     radii=600,
#                     translation=[-400, 0, 400],
#                     normal_vector=[30, 15, 60],
#                     number_of_vertices=5,
#                     aperture=1.0e-3)

# DFN.add_user_fract(shape='ell',
#                     radii=1000,
#                     translation=[0, 0, 0],
#                     normal_vector=[95, 5, 5],
#                     # normal_vector=[1, 0, 0],
#                     number_of_vertices=5,
#                     aperture=1.0e-3)

# DFN.add_user_fract(shape='ell',
#                     radii=600,
#                     aspect_ratio=1,
#                     translation=[400, 0, 200],
#                     normal_vector=[30, 15, 60],
#                     number_of_vertices=5,
#                     aperture=1.0e-3)

# DFN.add_user_fract(shape='ell',
#                     radii=600,
#                     aspect_ratio=1,
#                     translation=[400, 0, -400],
#                     normal_vector=[30, 15, 60],
#                     number_of_vertices=5,
#                     aperture=5.0e-5)

DFN.make_working_directory(delete=True)
DFN.print_domain_parameters()
DFN.check_input()
DFN.create_network()
DFN.dump_hydraulic_values()
DFN.output_report()

