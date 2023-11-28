from loadmodules import *
import pandas as pd

target_gas_mass = 3.74806e-06
hubbleparam = 0.67769999999999997
snap = 127
nhalos = 30

import shutil
# src =  '/home/lgvanover/intro_scripts/satellite_utilities.py'
# dst = '/home/lgvanover/arepo-snap-util/satellite_utilities.py'
# shutil.copyfile(src, dst)

import satellite_utilities as sat

def load_halo(select_halo, sim_level):
    '''
    creating a df with stats on a select halo's stellar particles, level 3 or level 4

    input:
        select_halo (int): halo number 
        sim_level (int): level for auriga simulation. 2, 3, or 4.
    output:
        datafram with the median epsilon for each snapshot and its corresponding redshift
    '''
    tf0 = time.time()
    if sim_level == 3:
        trees_sf1 = '063'
        snap = 63
        min_snap = 20
    elif sim_level == 4:
        trees_sf1 = '127'
        snap = 127
        min_snap = 30
    
    path = '/eagle/projects/Auriga/Auriga/level'+str(sim_level)+'/Original/halo_'+str(select_halo)+'/output/'
    treebase = '/eagle/projects/Auriga/Auriga/level'+str(sim_level)+'/Original/mergertrees/Au-'+str(select_halo)+'/trees_sf1_'+str(trees_sf1)
    print(treebase)
    t = load_tree(0,0,base=treebase,verbose=True)
    ifof = 0
    snap_numbers0,redshifts0,subfind_indices0,tree_indices0,ff_tree_indices0,fof_indices0 = t.return_subhalo_first_progenitors(0,snap=snap)

    Z = []
    Median_Eps = []

    for i in list(range(min_snap, snap)):
    #Load data
        attrs = ['pos', 'vel', 'mass', 'age', 'id', 'pot', 'sfr']
        s = gadget_readsnap(i,snappath=path,hdf5=True,loadonlytype=[4],loadonly=attrs)
        sf = load_subfind(i,dir=path,hdf5=True)

        s.calc_sf_indizes(sf)
        s.select_halo( sf, .1, use_principal_axis=True, use_cold_gas_spin=False, do_rotation=True ,haloid = fof_indices0[ifof])

        g = parse_particledata( s, sf, attrs, radialcut = 0.1*sf.data['frc2'][fof_indices0[ifof]] ,docircularities=True)
        # change flag in prep_data to normbycurve=True for new method
        g.prep_data(ageingyr=True, normbycurve=True)

        sdata = g.sgdata['sdata']
        eps2 = sdata['eps2']
        Z.append(sf.redshift)
        Median_Eps.append(median(eps2))
    
    d = {'Redshift': Z, 'Median Epsilon': Median_Eps}
    df = pd.DataFrame(data=d)
    filename = '/home/lgvanover/introscripts/halos/new/level'+str(sim_level)+'_halo'+str(select_halo)+'.hdf5'
    print('saving to', filename)
    df.to_hdf(filename, 'table')
    return df

if __name__ == '__main__':

    halo = sys.argv[1]
    level = sys.argv[2]
    print('loading halo '+str(halo)+' level '+str(level))
    df = load_halo(int(halo), int(level))
    print(df.head())