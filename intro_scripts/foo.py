from loadmodules import *
import pandas as pd
import time
# import satellite_utilities as sat

target_gas_mass = 3.74806e-06
hubbleparam = 0.67769999999999997
select_halo = 6
nhalos = 30

import shutil
# src = '/cds/kravtsov/Simulations/Auriga/satellite_utilities.py'
# dst =  '/home/lgvanover/intro_scripts/satellite_utilities.py'
# shutil.copyfile(src, dst)
import satellite_utilities as sat


halos_level_3 = [6, 16, 21, 23, 24, 27]


def load_halo(select_halo, sim_level, all_halos=False):
    '''
    creating a df with stats on a select halo's STELLAR mass, level 3 or level 4
    '''
    print(sim_level, select_halo)
    t0 = time.time()
    if sim_level == 3:
        trees_sf1 = '063'
        snap = 63
        # snap = 30
        if all_halos:
            min_snap = 40
        else:
            min_snap = 63
            # min_snap = 30
        path = '/eagle/projects/Auriga/Auriga/level'+str(sim_level)+'/Original/halo_'+str(select_halo)+'/output/'
        treebase = '/eagle/projects/Auriga/Auriga/level'+str(sim_level)+'/Original/mergertrees/Au-'+str(select_halo)+'/trees_sf1_'+str(trees_sf1)

    elif sim_level == 4:
        trees_sf1 = '127'
        snap = 127
        if all_halos:
            min_snap = 30
        else:
            min_snap = 127
        path = '/eagle/projects/Auriga/Auriga/level'+str(sim_level)+'/Original/halo_'+str(select_halo)+'/output/'
        treebase = '/eagle/projects/Auriga/Auriga/level'+str(sim_level)+'/Original/mergertrees/Au-'+str(select_halo)+'/trees_sf1_'+str(trees_sf1)
    
    # snap = 30

    t = load_tree(0,0,base=treebase,verbose=True)
    ifof = 0
    snap_numbers0,redshifts0,subfind_indices0,tree_indices0,ff_tree_indices0,fof_indices0 = t.return_subhalo_first_progenitors(0,snap=snap)

    eps = []
    z = []
    energy = []
    jzmax = []
    j_z = []
    

    #Load data
    attrs = ['pos', 'vel', 'mass', 'age', 'id', 'pot', 'sfr']
    s = gadget_readsnap(snap,snappath=path,hdf5=True,loadonlytype=[4],loadonly=attrs)
    sf = load_subfind(snap,dir=path,hdf5=True)
    # print('time elapsed for readsnap and subfind:', elapsed0)

    s.calc_sf_indizes(sf)
    s.select_halo( sf, .1, use_principal_axis=True, use_cold_gas_spin=False, do_rotation=True ,haloid = fof_indices0[ifof])
    # print('select halo time:', elapsed1)

    g = parse_particledata( s, sf, attrs, radialcut = 0.1*sf.data['frc2'][fof_indices0[ifof]] ,docircularities=True)
    # g.prep_data(ageingyr=True)
    return g