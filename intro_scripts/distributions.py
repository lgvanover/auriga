from loadmodules import *
import pandas as pd
import time

target_gas_mass = 3.74806e-06
hubbleparam = 0.67769999999999997
select_halo = 6
nhalos = 30
import sys

import shutil
import satellite_utilities as sat


halos_level_3 = [6, 16, 21, 23, 24, 27]


def load_halo(select_halo, sim_level, all_halos=False, age_filt=None):
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
    
    t = load_tree(0,0,base=treebase,verbose=True)
    ifof = 0
    snap_numbers0,redshifts0,subfind_indices0,tree_indices0,ff_tree_indices0,fof_indices0 = t.return_subhalo_first_progenitors(0,snap=snap)

    
    # smass = []

    for isnap in list(range(min_snap, snap + 1)):
        eps = []
        z = []

        #Load data
        attrs = ['pos', 'vel', 'mass', 'age', 'id', 'pot', 'sfr']
        s = gadget_readsnap(isnap,snappath=path,hdf5=True,loadonlytype=[4],loadonly=attrs)
        sf = load_subfind(isnap,dir=path,hdf5=True)

        s.calc_sf_indizes(sf)
        s.select_halo( sf, .1, use_principal_axis=True, use_cold_gas_spin=False, do_rotation=True ,haloid = fof_indices0[ifof])

        # g = parse_particledata( s, sf, attrs, radialcut = 0.1*sf.data['frc2'][fof_indices0[ifof]] ,docircularities=True)
        g = parse_particledata( s, sf, attrs, radialcut = 0.1*sf.data['frc2'][fof_indices0[ifof]] ,docircularities=False)
        g.prep_data(ageingyr=True)


        istars, = np.where( (g.s.data['type'] == 4) & (g.s.data['age'] > 0) & (g.s.r() < g.radialcut) )

        nstars = size(istars)
        eps2  = pylab.zeros( nstars )
        jcmax = pylab.zeros( nstars )
        spec_energy = pylab.zeros( nstars )
        pos = g.s.data['pos']
        vel = g.s.data['vel']
        pot = g.s.data['pot']
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
            if slope <= 0:
                slope = slope_E[-1]
                yint = y_int[-1]
            else:
                yint = E_est[i] - slope*jz_est[i]
            slope_E.append(slope)
            y_int.append(yint)

        # jz_max = []
        epsilon = []
        # jz_max = np.zeros_like(E_sorted)

        from parallel_decorators import vectorize_parallel

        print("\nDEBUGGING\n")

        @vectorize_parallel(method='processes', num_procs=11)
        def slice_normalization(i, esorted, slope, intercept):

            slice_slope = slope[i]
            slice_intercept = intercept[i]

            jz_max = (esorted[split_index[i]:split_index[i+1]] - slice_intercept) / slice_slope
            print(i, jz_max)
            return jz_max 
        print("\nDEBUGGING\n")

        results = array(slice_normalization(arange(len(split_index[0:-1])), E_sorted, slope_E, y_int))
        results = np.concatenate(results).ravel()
        results = np.append(results, (E_sorted[split_index[-1]] - y_int[-1]) / slope_E[-1])
        epsilon = jz_sorted / results

        isort = np.argsort(E)
        sdata = g.sgdata['sdata']

        if age_filt:
            # snap_age = sat.return_lookbacktime_from_a((np.array([s.time])+1.0)**(-1.0))[0]
            snap_age = sat.return_lookbacktime_from_a(np.array([s.time]))[0]
            index = where((sdata['age'] - snap_age) < age_filt)
            epsilon = epsilon[index]
            print(max(sdata['age']), min(sdata['age']))
            print(snap_age)
            print(s.time)

        # print("KEYS")
        # print(sdata.keys())
        # smass = sdata['mass']
        print('AGE FILT E:', epsilon)

        z.append(sf.redshift)
        # smass.append(smass)
        # weights = np.ones_like(epsilon) / len(epsilon)
        h = hist(epsilon,bins=50,histtype='step',density=True)
        y = h[0]
        edges = h[1]
        
        a = np.zeros_like(edges)
        a[:len(y)] = y
        # ed = np.zeros_like(epsilon)
        # ed[:len(edges)] = edges 
        red = np.zeros_like(edges)
        red[:len(z)] = z


        d = {'Redshift': red, 'Bin_Y': a, 'Bin_Edges': edges}
        df = pd.DataFrame(data=d)
        # df = pd.DataFrame.from_dict(d, orient='index')
        # df = df.transpose()
        if age_filt:
            filename = 'halos/new/distr/age_filt/age'+str(age_filt)+'/halo'+str(select_halo)+'_snap'+str(isnap)+'.hdf5'
        else:
            filename = 'halos/new/distr/halo'+str(select_halo)+'_snap'+str(isnap)+'.hdf5'
        df.to_hdf(filename,'table')
    return

if __name__ == '__main__':

    halo = sys.argv[1]
    level = sys.argv[2]
    # load_halo(int(halo), int(level), all_halos=True, age_filt=10)
    df = load_halo(int(halo), int(level))
