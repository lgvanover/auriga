from loadmodules import *
import pandas as pd
import time
import sys

target_gas_mass = 3.74806e-06
hubbleparam = 0.67769999999999997
select_halo = 6
nhalos = 30

import shutil
import satellite_utilities as sat

def load_halo(select_halo, sim_level, snap):
    '''
    creating a df with stats on a select halo's STELLAR mass, level 3 or level 4, specific snapshot
    '''
    print(select_halo, sim_level, snap)
    t0 = time.time()
    if sim_level == 3:
        trees_sf1 = '063'
        d = 'Au-'+str(select_halo)
    elif sim_level == 4:
        trees_sf1 = '127'
        d = 'Au-'+str(select_halo)
    elif sim_level == 2:
        trees_sf1 = '127'
        d = 'halo_'+str(select_halo)

            
    path = '/eagle/projects/Auriga/Auriga/level'+str(sim_level)+'/Original/halo_'+str(select_halo)+'/output/'
    print(path)
    treebase = '/eagle/projects/Auriga/Auriga/level'+str(sim_level)+'/Original/mergertrees/'+d+'/trees_sf1_'+str(trees_sf1)
    print(treebase)
    

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

    s.calc_sf_indizes(sf)
    s.select_halo( sf, .1, use_principal_axis=True, use_cold_gas_spin=False, do_rotation=True ,haloid = fof_indices0[ifof])

    g = parse_particledata( s, sf, attrs, radialcut = 0.1*sf.data['frc2'][fof_indices0[ifof]] ,docircularities=True)

    istars, = np.where( (g.s.data['type'] == 4) & (g.s.data['age'] > 0) & (g.s.r() < g.radialcut) )

    nstars = size(istars)
    eps2  = pylab.zeros( nstars )
    jcmax = pylab.zeros( nstars )
    spec_energy = pylab.zeros( nstars )

    pos = g.s.data['pos']
    vel = g.s.data['vel']
    pot = g.s.data['pot']
    stellarm = sf.data['smit'][:,4]

    j  = pylab.cross( pos[istars,:], vel[istars,:] )
    jz = j[:,0]

    E = 0.5 * (vel[istars,:]**2).sum(axis=1) + pot[istars]
    E = E/1e5

    isort = np.argsort(E)
    jz_sorted = jz[isort]
    E_sorted = E[isort]

    nslices = 10
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

    epsilon = []

    from parallel_decorators import vectorize_parallel

    @vectorize_parallel(method='processes', num_procs=11)
    def slice_normalization(i, esorted, slope, intercept):

        slice_slope = slope[i]
        slice_intercept = intercept[i]

        jz_max = (esorted[split_index[i]:split_index[i+1]] - slice_intercept) / slice_slope
        return jz_max 

    results = slice_normalization(arange(len(split_index[0:-1])), E_sorted, slope_E, y_int)
    results = np.concatenate(results).ravel()
    results = np.append(results, (E_sorted[split_index[-1]] - y_int[-1]) / slope_E[-1])
    epsilon = jz_sorted / results
    isort = np.argsort(E)

    eps.append(median(epsilon))
    z.append(sf.redshift)
    energy.append(median(E_sorted * 1e5))
    jzmax.append(median(results))
    j_z.append(median(jz_sorted))

    d = {'Redshift': z, 'Epsilon': eps, 'Spec Energy': energy, 'Jz': j_z, 'Jz Max': jzmax}
    df = pd.DataFrame(data=d)
    t1 = time.time()
    elapsed = t1 - t0 
    print('TOTAL TIME:', elapsed, '\n')
    return df

if __name__ == '__main__':

    halo = sys.argv[1]
    level = sys.argv[2]
    s = sys.argv[3]
    load_halo(int(halo), int(level), int(s))
