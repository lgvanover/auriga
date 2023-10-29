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


def load_halo(select_halo, sim_level):
    '''
    creating a df with stats on a select halo's STELLAR mass, level 3 or level 4
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

    x_data = []
    y_data = []
    ratio = []
    med_smass = []
    med_velocity = []
    med_position = []
    med_sfr = []

    for i in list(range(min_snap, snap)):
    #Load data
        t0 = time.time()
        attrs = ['pos', 'vel', 'mass', 'age', 'id', 'pot', 'sfr']
        s = gadget_readsnap(i,snappath=path,hdf5=True,loadonlytype=[4],loadonly=attrs)
        sf = load_subfind(i,dir=path,hdf5=True)
        t1 = time.time()
        elapsed0 = t1 - t0
        print('time elapsed for readsnap and subfind:', elapsed0)

        s.calc_sf_indizes(sf)
        s.select_halo( sf, .1, use_principal_axis=True, use_cold_gas_spin=False, do_rotation=True ,haloid = fof_indices0[ifof])
        t2 = time.time()
        elapsed1 = t2 - t1
        print('select halo time:', elapsed1)

        g = parse_particledata( s, sf, attrs, radialcut = 0.1*sf.data['frc2'][fof_indices0[ifof]] ,docircularities=True)
        g.prep_data(ageingyr=True)
        t3 = time.time()
        elapsed2 = t3 - t2
        print('parsing and prepping time:', elapsed2)

        sdata = g.sgdata['sdata']
        eps2 = sdata['eps2']
        smass = sdata['mass']
        vel = sdata['vel']
        pos = sdata['pos']

        x_data.append(sf.redshift)
        y_data.append(median(eps2))
        med_smass.append(median(smass))
        med_velocity.append(median(vel))
        med_position.append(median(pos))

        index_neg = where(eps2 <0)
        h = histogram(eps2,bins=50,density=False,weights=smass)
        b = histogram(eps2[index_neg],bins=h[1],density=False,weights=smass[index_neg])
        
        delta_x = h[1][1] - h[1][0]
        bulge = 2* sum(delta_x * b[0])
        total = sum(delta_x*h[0])
        
        ratio.append((total-bulge)/total)
        t4 = time.time()
        elapsed3 = t4 - t3
        print('appending n hist time:', elapsed3)
        print(max([elapsed0,elapsed1,elapsed2,elapsed3]))

    d = {'Redshift': x_data, 'Med_Eps': y_data, 'Disk_Ratio': ratio, 'Med_SMass': med_smass, 'Med_Velocity': med_velocity, 'Med_Position': med_position}
    df = pd.DataFrame(data=d)
    tf1 = time.time()
    total_time = tf1 - tf0
    print('total halo time:', total_time)
    return df


def halo_gas(select_halo, sim_level):
    '''
    creating a df with stats on a select halo's GAS mass, level 3 or level 4
    '''
    if sim_level == 3:
        trees_sf1 = '063'
        snap = 63
        min_snap = 62
    elif sim_level == 4:
        trees_sf1 = '127'
        snap = 127
        min_snap = 30
    
    path = '/eagle/projects/Auriga/Auriga/level'+str(sim_level)+'/Original/halo_'+str(select_halo)+'/output/'
    treebase = '/eagle/projects/Auriga/Auriga/level'+str(sim_level)+'/Original/mergertrees/Au-'+str(select_halo)+'/trees_sf1_'+str(trees_sf1)
    t = load_tree(0,0,base=treebase,verbose=True)
    snap_numbers0,redshifts0,subfind_indices0,tree_indices0,ff_tree_indices0,fof_indices0 = t.return_subhalo_first_progenitors(0,snap=snap)
    ifof = 0

    z = []
    ratio = []
    gmass = []
    med_pos = []
    med_vel = []
    med_sfr = []

    for i in list(range(min_snap, snap)):
    #Load data

        attrs = ['pos', 'vel', 'mass', 'age', 'id', 'pot', 'sfr']
        s = gadget_readsnap(i,snappath=path,hdf5=True,loadonlytype=[0, 4],loadonly=attrs)
        sf = load_subfind(i,dir=path,hdf5=True)

        s.calc_sf_indizes(sf)
        s.select_halo( sf, 3., use_principal_axis=True, use_cold_gas_spin=False, do_rotation=True ,haloid = fof_indices0[ifof])

        g = parse_particledata( s, sf, attrs, radialcut = 0.1*sf.data['frc2'][fof_indices0[ifof]] ,docircularities=True)
        g.prep_data()

        gdata = g.sgdata['gdata']
        sfr = gdata['sfr']
        sfr_nonzero = where(sfr > 0)
        sfr = sfr[sfr_nonzero]

        z.append(sf.redshift)
        gmass.append(median(gdata['mass'][sfr_nonzero]))
        med_pos.append(median(gdata['pos'][sfr_nonzero]))
        med_vel.append(median(gdata['vel'][sfr_nonzero]))
        med_sfr.append(median(gdata['sfr'][sfr_nonzero]))

    d = {'Redshift': z, 'Med_GMass': gmass, 'Med_Pos': med_pos, 'Med_Vel': med_vel, 'Med_SFR': med_sfr}
    df = pd.DataFrame(data=d)
    return df


def halo_df_to_hdf(select_halo, sim_level, type, filter_age=None):
    '''
    saves halo df to hdf5
    '''
    tf0 = time.time()
    if filter_age is None:
        if type == 0:
            df = halo_gas(select_halo, sim_level)
            if sim_level == 3:
                df.to_hdf("halos/gas/level3_gas/halo_"+str(select_halo)+"_level3_gas.hdf5", "table")
            elif sim_level == 4:
                df.to_hdf("halos/gas/level4_gas/halo_"+str(select_halo)+"_gas.hdf5", "table")
        elif type == 4:
            df = load_halo(select_halo, sim_level)
            if sim_level == 3:
                df.to_hdf("halos/level3/halo_"+str(select_halo)+"_level3.hdf5", "table")
            elif sim_level == 4:
                # for testing npm
                df.to_hdf("halos/level4/test/halo_"+str(select_halo)+"test1.hdf5", "table")
    else:
        assert filter_age == 30 or 10 or 1
        df = load_halo_age_filt(select_halo, sim_level, filter_age)
        if sim_level == 3:
            df.to_hdf('halos/level3/age_filt/age'+str(filter_age)+'/halo_'+str(select_halo)+'level3.hdf5', 'table')
        # elif sim_level == 4:
        #     df.to_hdf('halos/level4/age_filt/age'+str(filter_age)+'/halo_'+str(select_halo)+'.hdf5')
    tf1 = time.time()
    total_time = tf1 - tf0
    print('total halo time:', total_time)
    return


halo_df_to_hdf(24, 4, 4)

def load_halo_age_filt(select_halo, sim_level, filter_age):
    '''
    creating a df with stats on a select halo's STELLAR mass, level 3 or level 4
    '''
    
    if sim_level == 3:
        trees_sf1 = '063'
        snap = 63
        min_snap = 25
    elif sim_level == 4:
        trees_sf1 = '127'
        snap = 127
        min_snap = 30

    path = '/eagle/projects/Auriga/Auriga/level'+str(sim_level)+'/Original/halo_'+str(select_halo)+'/output/'
    treebase = '/eagle/projects/Auriga/Auriga/level'+str(sim_level)+'/Original/mergertrees/Au-'+str(select_halo)+'/trees_sf1_'+str(trees_sf1)
    t = load_tree(0,0,base=treebase,verbose=True)
    ifof = 0
    snap_numbers0,redshifts0,subfind_indices0,tree_indices0,ff_tree_indices0,fof_indices0 = t.return_subhalo_first_progenitors(0,snap=snap)

    x_data = []
    y_data = []
    ratio = []
    med_smass = []
    med_velocity = []
    med_position = []
    med_sfr = []

    for i in list(range(min_snap, snap)):
    #Load data
        attrs = ['pos', 'vel', 'mass', 'age', 'id', 'pot', 'sfr']
        s = gadget_readsnap(i,snappath=path,hdf5=True,loadonlytype=[4],loadonly=attrs)
        sf = load_subfind(i,dir=path,hdf5=True)

        s.calc_sf_indizes(sf)
        s.select_halo( sf, 3., use_principal_axis=True, use_cold_gas_spin=False, do_rotation=True ,haloid = fof_indices0[ifof])

        g = parse_particledata( s, sf, attrs, radialcut = 0.1*sf.data['frc2'][fof_indices0[ifof]] ,docircularities=True)
        g.prep_data(ageingyr=True)

        sdata = g.sgdata['sdata']
        snap_age = sat.return_lookbacktime_from_a((np.array([s.time])+1.0)**(-1.0))[0]
        index = where((sdata['age'] - snap_age) < filter_age)

        eps2 = sdata['eps2'][index]
        smass = sdata['mass'][index]
        vel = sdata['vel'][index]
        pos = sdata['pos'][index]

        x_data.append(sf.redshift)
        y_data.append(median(eps2))
        med_smass.append(median(smass))
        med_velocity.append(median(vel))
        med_position.append(median(pos))

        index_neg = where(eps2 <0)
        h = histogram(eps2,bins=50,density=False,weights=smass)
        b = histogram(eps2[index_neg],bins=h[1],density=False,weights=smass[index_neg])
        
        delta_x = h[1][1] - h[1][0]
        bulge = 2* sum(delta_x * b[0])
        total = sum(delta_x*h[0])
        
        ratio.append((total-bulge)/total)

    d = {'Redshift': x_data, 'Med_Eps': y_data, 'Disk_Ratio': ratio, 'Med_SMass': med_smass, 'Med_Velocity': med_velocity, 'Med_Position': med_position}
    df = pd.DataFrame(data=d)
    return df
