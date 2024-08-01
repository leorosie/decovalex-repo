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
import h5py
import os
import numpy as np 
import subprocess 
import sys

def dump_h5_files(DFN, porosity):
    """ Write permeability values to cell ids and permeability values to dfn_properties.h5 file for pflotran. 

    Parameters
    ----------
        self : object
            DFN Class

    Returns
    ---------
        None

    Notes
    ----------
        Hydraulic properties need to attached to the class prior to running this function. Use DFN.assign_hydraulic_properties() to do so. 
    """
    print('*' * 80)
    print("--> Dumping h5 file")
    fp_perm = 'permeability.h5'
    print(f'\n--> Opening HDF5 File {fp_perm}')
    DFN.perm_cell = np.zeros(DFN.num_nodes)
    with h5py.File(fp_perm, mode='w') as h5file:
        print('--> Allocating cell index array')
        print('--> Writing cell indices')
        iarray = np.arange(1,DFN.num_nodes + 1)
        dataset_name = 'Cell Ids'
        h5dset = h5file.create_dataset(dataset_name, data=iarray)
        print('--> Creating permeability array')
        print('--> Note: This script assumes isotropic permeability')
        for i in range(DFN.num_nodes):
            DFN.perm_cell[i] = DFN.perm[DFN.material_ids[i] - 1]
        print('--> Writting Permeability')
        dataset_name = 'Permeability'
        h5dset = h5file.create_dataset(dataset_name, data=DFN.perm_cell)

    fp_porosity = 'porosity.h5'
    print(f'\n--> Opening HDF5 File {fp_porosity}')
    with h5py.File(fp_porosity, mode='w') as h5file:
        print('--> Allocating cell index array')
        print('--> Writing cell indices')
        iarray = np.arange(1,DFN.num_nodes + 1)
        dataset_name = 'Cell Ids'
        h5dset = h5file.create_dataset(dataset_name, data=iarray)
        print('--> Creating porosity array')
        porosity_cell = np.zeros_like(DFN.perm_cell)
        for i in range(DFN.num_nodes):
            porosity_cell[i] = porosity[DFN.material_ids[i] - 1]
        print('--> Writting Porosity')
        dataset_name = 'Porosity'
        h5dset = h5file.create_dataset(dataset_name, data=porosity_cell)

    fp_matid = 'materials.h5'
    print(f'\n--> Opening HDF5 File {fp_matid}')
    with h5py.File(fp_matid, mode='w') as h5file:
        print('--> Allocating cell index array')
        print('--> Writing cell indices')
        dataset_name = 'Materials/Cell Ids'
        h5dset = h5file.create_dataset(dataset_name, data=iarray)
        dataset_name = 'Materials/Material Ids'
        h5dset = h5file.create_dataset(dataset_name, data=DFN.material_ids)


    print("--> Done writting h5 file")
    print('*' * 80)
    print()

src_path = os.getcwd()

jobname = f"{src_path}/output"
dfnFlow_file = f"{src_path}/udfm_repo.in"

DFN = DFNWORKS(jobname,
               dfnFlow_file=dfnFlow_file,
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
DFN.mesh_network()

DFN.map_to_continuum(l=100, orl=3)
DFN.upscale(mat_perm=1e-15, mat_por=0.01)

DFN.zone2ex(zone_file='all')

cmd = f'lagrit < {src_path}/tag_cells.lgi'
subprocess.call(cmd, shell=True)
DFN.uge_file = 'full_mesh.uge'
# DFN.zone2ex(zone_file='boundary_front.zone', face = 'north')
# DFN.zone2ex(zone_file='boundary_bottom.zone', face = 'south')
# DFN.zone2ex(zone_file='boundary_left_w.zone', face = 'west')
DFN.zone2ex(zone_file='repo.zone', face = 'east')

DFN.material_ids = np.genfromtxt('materialid.dat',
                                 skip_header=3).astype(int)

DFN.num_nodes = len(DFN.material_ids)
## Settting perm / porosity parameters 
matrix_perm = 1e-12
repo_perm = 1e-10
matrix_porosity = 0.01
repo_porosity = 0.5 

DFN.perm[0] = matrix_perm
DFN.perm[1] = repo_perm

porosity = [matrix_porosity, repo_porosity]

dump_h5_files(DFN, porosity) 

DFN.ncpu = 16
DFN.pflotran()
DFN.parse_pflotran_vtk_python()
DFN.pflotran_cleanup() 

