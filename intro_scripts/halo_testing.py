from loadmodules import *
import pandas as pd
import time
import sys
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
            min_snap = 22
        else:
            min_snap = 63
            # min_snap = 30
        path = '/eagle/projects/Auriga/Auriga/level'+str(sim_level)+'/Original/halo_'+str(select_halo)+'/output/'
        treebase = '/eagle/projects/Auriga/Auriga/level'+str(sim_level)+'/Original/mergertrees/Au-'+str(select_halo)+'/trees_sf1_'+str(trees_sf1)


    elif sim_level == 4:
        trees_sf1 = '127'
        snap = 127
        if all_halos:
            min_snap = 0
        else:
            min_snap = 127
            
        path = '/eagle/projects/Auriga/Auriga/level'+str(sim_level)+'/Original/halo_'+str(select_halo)+'/output/'
        print(path)
        treebase = '/eagle/projects/Auriga/Auriga/level'+str(sim_level)+'/Original/mergertrees/Au-'+str(select_halo)+'/trees_sf1_'+str(trees_sf1)
        print(treebase)

    elif sim_level == 2:
        trees_sf1 = '127'
        snap = 127
        if all_halos:
            min_snap = 25
        else:
            min_snap = 127
            
        path = '/eagle/projects/Auriga/Auriga/level'+str(sim_level)+'/Original/halo_'+str(select_halo)+'/output/'
        print(path)
        treebase = '/eagle/projects/Auriga/Auriga/level'+str(sim_level)+'/Original/mergertrees/halo_'+str(select_halo)+'/trees_sf1_'+str(trees_sf1)
        print(treebase)
    

    t = load_tree(0,0,base=treebase,verbose=True)
    ifof = 0
    snap_numbers0,redshifts0,subfind_indices0,tree_indices0,ff_tree_indices0,fof_indices0 = t.return_subhalo_first_progenitors(0,snap=snap)
    print("KEYS")
    print(t.data.keys())

    smas = []

    eps = []
    z = []
    energy = []
    jzmax = []
    j_z = []
    # for i in list(range(min_snap, snap)):


    for i in list(range(min_snap, snap + 1)):

        #Load data
        attrs = ['pos', 'vel', 'mass', 'age', 'id', 'pot', 'sfr']
        s = gadget_readsnap(i,snappath=path,hdf5=True,loadonlytype=[4],loadonly=attrs)
        sf = load_subfind(i,dir=path,hdf5=True)
        # print('time elapsed for readsnap and subfind:', elapsed0)

        s.calc_sf_indizes(sf)
        s.select_halo( sf, .1, use_principal_axis=True, use_cold_gas_spin=False, do_rotation=True ,haloid = fof_indices0[ifof])
        # print('select halo time:', elapsed1)

        g = parse_particledata( s, sf, attrs, radialcut = 0.1*sf.data['frc2'][fof_indices0[ifof]] ,docircularities=True)
        # g.prep_data(ageingyr=True)

        istars, = np.where( (g.s.data['type'] == 4) & (g.s.data['age'] > 0) & (g.s.r() < g.radialcut) )

        nstars = size(istars)
        eps2  = pylab.zeros( nstars )
        jcmax = pylab.zeros( nstars )
        spec_energy = pylab.zeros( nstars )

        smas.append(sf.data['smas'])

        print("KEYS2")
        print(g.s.data.keys())
        pos = g.s.data['pos']
        vel = g.s.data['vel']
        pot = g.s.data['pot']
        # smas.append(t.data['smit'])
        stellarm = sf.data['smit'][:,4]

        j  = pylab.cross( pos[istars,:], vel[istars,:] )
        jz = j[:,0]

        E = 0.5 * (vel[istars,:]**2).sum(axis=1) + pot[istars]
        E = E/1e5

        isort = np.argsort(E)
        jz_sorted = jz[isort]
        E_sorted = E[isort]

        nslices = 10
        print("len esorted", len(E_sorted))
        deltaE = (max(E) - min(E))/nslices
        Emin = min(E)

        E_est=[]
        jz_est=[]
        split_index = []

        E_est.append(E_sorted[0])
        jz_est.append(jz_sorted[0])
        split_index.append(0)
        nselect = nstars//nslices

        for i in range(nslices):

            max_jz = max(abs(jz_sorted[i*nselect:(i+1)*nselect]))
            imax = np.where(abs(jz_sorted)==max_jz)[0][0]
            split_index.append(imax)
            E_est.append(E_sorted[imax])
            jz_est.append(abs(max_jz))

        Emax = max(E)
        imax_final = np.where(E_sorted == Emax)
        E_est.append(E_sorted[imax_final][0])
        jz_est.append(jz_sorted[imax_final])
        split_index.append(imax_final[0][0])

        npart =len(E)
        slope_E = []
        y_int = []

        for i in range(1, len(E_est)):

            slope = (E_est[i]-E_est[i-1])/(jz_est[i]-jz_est[i-1])
            # if slope <= 0:
            #     slope = slope_E[-1]
            #     yint = y_int[-1]
            # else:
            #     yint = E_est[i] - slope*jz_est[i]
            yint = E_est[i] - slope*jz_est[i]
            
            slope_E.append(slope)
            y_int.append(yint)

        # jz_max = []
        epsilon = []
        # jz_max = np.zeros_like(E_sorted)

        from parallel_decorators import vectorize_parallel

        @vectorize_parallel(method='processes', num_procs=11)
        def slice_normalization(i, esorted, slope, intercept):

            slice_slope = slope[i]
            slice_intercept = intercept[i]

            jz_max = (esorted[split_index[i]:split_index[i+1]] - slice_intercept) / slice_slope
            return jz_max 
        
        print("split index", split_index)
        print("esorted", E_sorted)
        print("slope", slope_E)
        print("yint", y_int)

        results = slice_normalization(arange(len(split_index[0:-1])), E_sorted, slope_E, y_int)
        # results = array(slice_normalization(arange(len(split_index[0:-1])), E_sorted, slope_E, y_int))
        # print(results)
        # print(size(results))
        # results[nslices*nselect:] = (E_sorted[split_index[-1]] - y_int[-1]) / slope_E[-1]
        
        # print(len(jz_sorted))
        # print(results)
        results = np.concatenate(results).ravel()
        results = np.append(results, (E_sorted[split_index[-1]] - y_int[-1]) / slope_E[-1])
        # print(len(results))
        # print(results)
        epsilon = jz_sorted / results
        # print(epsilon)

        # for i in range(len(split_index[0:-1])):

        #     slice_slope = slope_E[i]
        #     slice_intercept = y_int[i]

        #     jz_max[split_index[i]:split_index[i+1]] = (E_sorted[split_index[i]:split_index[i+1]] - slice_intercept) / slice_slope

        # jz_max[-1] = (E_sorted[split_index[11]] - slice_intercept) / slice_slope
        # epsilon = jz_sorted / jz_max

        isort = np.argsort(E)

        # if all_halos:
        #     eps.append(median(epsilon))
        #     z.append(sf.redshift)
        #     energy.append(median(energy))
        #     jzmax.append(median(results))

        # else:
        #     eps.append(epsilon[isort])
        #     z.append(sf.redshift)
        #     energy.append(E[isort] * 1e5)
        #     jzmax.append(results)
        # #     j_z.append(jz)
        eps.append(median(epsilon))
        z.append(sf.redshift)
        energy.append(median(E_sorted * 1e5))
        jzmax.append(median(results))
        j_z.append(median(jz_sorted))

        # print('Eps', eps[0])
        # print('Eps m', median(eps))
        # print('Jz', j_z[0])
        # print('Jz m', median(j_z[0]))
        # print('Jz max', jzmax[0])
        # print('Jz max m', median(jzmax[0]))
        # print('Energy', energy[0])
        # print('Energy m', median(energy[0]))

    print('SMAS:', smas)
    d = {'Redshift': z, 'Epsilon': eps, 'Spec Energy': energy, 'Jz': j_z, 'Jz Max': jzmax, 'smas': smas}
    df = pd.DataFrame(data=d)
    t1 = time.time()
    elapsed = t1 - t0 
    print('TOTAL TIME:', elapsed, '\n')
    return df

# load_halo(24, 4)

def actual_halo(select_halo, sim_level, all_halos=False):
    '''
    creating a df with stats on a select halo's STELLAR mass, level 3 or level 4
    '''
    print('test\n')
    t0 = time.time()
    if sim_level == 3:
        trees_sf1 = '063'
        snap = 63
        # snap = 30
        if all_halos:
            min_snap = 22
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

    t = load_tree(0,0,base=treebase,verbose=True)
    ifof = 0
    snap_numbers0,redshifts0,subfind_indices0,tree_indices0,ff_tree_indices0,fof_indices0 = t.return_subhalo_first_progenitors(0,snap=snap)

    eps = []
    z = []
    spec_energy = []
    j_z = []
    jzmax = []

    # for i in list(range(min_snap, snap)):
    for i in list(range(min_snap, snap + 1)):
        
        #Load data
        attrs = ['pos', 'vel', 'mass', 'age', 'id', 'pot', 'sfr']
        s = gadget_readsnap(i,snappath=path,hdf5=True,loadonlytype=[4],loadonly=attrs)
        sf = load_subfind(i,dir=path,hdf5=True)
        # print('time elapsed for readsnap and subfind:', elapsed0)

        s.calc_sf_indizes(sf)
        s.select_halo( sf, .1, use_principal_axis=True, use_cold_gas_spin=False, do_rotation=True ,haloid = fof_indices0[ifof])
        # print('select halo time:', elapsed1)

        g = parse_particledata( s, sf, attrs, radialcut = 0.1*sf.data['frc2'][fof_indices0[ifof]] ,docircularities=True)
        g.prep_data(ageingyr=True)

        sdata = g.sgdata['sdata']
        eps2 = sdata['eps2']
        se = sdata['spec_en']
        isort = np.argsort(se)

        
        eps.append(median(eps2))
        z.append(sf.redshift)
        # print(np.sort(se))
        spec_energy.append(median(se[isort]))
        # j_z.append(median(sdata['jz']))
        # jzmax.append(median(sdata['jz_max']))
        # spec_energy.append(np.sort(se))
        print(spec_energy)
        print(eps)
    
    d = {'Redshift': z, 'Epsilon': eps, 'Spec Energy': spec_energy}
    df = pd.DataFrame(data=d)
    t1 = time.time()
    elapsed = t1 - t0
    print('TOTAL TIME:', elapsed, '\n')
    
    return df

# actual_halo(24, 4)

def halo_df_to_hdf(select_halo, sim_level, method, all_h=False):
    if method == 'old':
        print('old\n')
        df = actual_halo(select_halo, sim_level, all_halos=all_h)
        if all_h:
            df.to_hdf('halos/new/level'+str(sim_level)+'_halo'+str(select_halo)+'_all_old.hdf5', 'table')
        else:
            df.to_hdf('halos/new/level'+str(sim_level)+'_halo'+str(select_halo)+'_old_snap.hdf5', 'table')
    elif method == 'new':
        print('new\n')
        df = load_halo(select_halo, sim_level, all_halos=all_h)
        if all_h:
            df.to_hdf('halos/new/level'+str(sim_level)+'_halo'+str(select_halo)+'_all_new.hdf5', 'table')
        else:
            df.to_hdf('halos/new/level'+str(sim_level)+'_halo'+str(select_halo)+'_new_snap.hdf5', 'table')
    return

if __name__ == '__main__':

    halo = sys.argv[1]
    level = sys.argv[2]
    m = sys.argv[3]
    # if len(sys.argv) > 3:
    #     allh = True
    # else:
    #     allh = False
    halo_df_to_hdf(int(halo), int(level), m, all_h=True)
    # load_halo(int(halo), int(level), all_halos=False)

# halo_df_to_hdf(6, 2, 'new', all_h=True)