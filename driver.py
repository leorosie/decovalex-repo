#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 21:28:30 2024

@author: rosie
"""

 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 13:44:47 2024

@author: rleone
"""

from pydfnworks import *
import os

src_path = os.getcwd()

jobname = f"{src_path}/output"
dfnFlow_file = f"{src_path}/dfn_explicit.in"
dfnTrans_file = f"{src_path}/PTDFN_control.dat"

DFN = DFNWORKS(jobname,
               dfnFlow_file=dfnFlow_file,
               dfnTrans_file=dfnTrans_file,
               ncpu=4)

DFN.params['domainSize']['value'] = [1500.0, 1000.0, 200.0]
DFN.params['h']['value'] = 3.0
DFN.params['tripleIntersections']['value'] = True
DFN.params['ignoreBoundaryFaces']['value'] = True
DFN.params['keepIsolatedFractures']['value'] = True
DFN.params['disableFram']['value'] = True

DFN.add_user_fract(shape='rect',
                   radii=510,
                   aspect_ratio=0.694,
                   translation=[0, 0, 5],
                   normal_vector=[0, 0, 1],
                   number_of_vertices=8,
                   aperture=1.0)

DFN.add_user_fract(shape='rect',
                   radii=510,
                   aspect_ratio=0.694,
                   translation=[0, 0, 0],
                   normal_vector=[0, 0, 1],
                   number_of_vertices=8,
                   aperture=1.0)

DFN.add_user_fract(shape='rect',
                   radii=510,
                   aspect_ratio=0.694,
                   translation=[0, 0, -5],
                   normal_vector=[0, 0, 1],
                   number_of_vertices=8,
                   aperture=1.0)


DFN.make_working_directory(delete=True)
DFN.print_domain_parameters()
DFN.check_input()
DFN.create_network()
# DFN.output_report()
DFN.mesh_network(uniform_mesh=True)

DFN.map_to_continuum(l=100, orl=3)
DFN.upscale(mat_perm=1e-15, mat_por=0.01)

DFN.zone2ex(zone_file='all')