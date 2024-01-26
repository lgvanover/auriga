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

# import satellite_utilities as sat

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
        # g.prep_data(ageingyr=True, normbycurve=True)
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

        spec_energy[:] = 0.5 * (vel[istars,:]**2).sum(axis=1) + pot[istars]

        spec_energy = spec_energy/1e5

        isort = np.argsort(spec_energy)
        jz_sorted = jz[isort]
        E_sorted = spec_energy[isort]
        nslices = 10
        delta_E = (max(spec_energy) - min(spec_energy)) / nslices
        Emin = min(spec_energy)

        E_est = []
        jz_est = []
        isplit = []

        E_est.append(E_sorted[0])
        jz_est.append(jz_sorted[0])
        isplit.append(0)
        nselect = nstars // nslices
        for idx in range(nslices):

                max_jz = max(abs(jz_sorted[idx*nselect:(idx+1)*nselect]))
                imax = np.where(abs(jz_sorted)==max_jz)[0][0]
                isplit.append(imax)
                E_est.append(E_sorted[imax])
                jz_est.append(abs(max_jz))

        Emax = max(spec_energy)
        imax_final = np.where(E_sorted == Emax)
        E_est.append(E_sorted[imax_final][0])
        jz_est.append(jz_sorted[imax_final][0])
        isplit.append(imax_final[0][0])

        slope_E = []
        y_int = []

        for idx in range(1, len(E_est)):
            div = jz_est[idx]-jz_est[idx-1]
            if div == 0:
              div = jz_est[idx]
            slope = (E_est[idx]-E_est[idx-1])/div
            if slope <= 0:
                if len(slope_E)==0:
                  slope = (E_est[idx+1]-E_est[idx]/div)
                  yint = E_est[idx]
                else:
                  slope = slope_E[-1]
                  yint = y_int[-1]
            else:
                yint = E_est[idx] - slope*jz_est[idx]
            slope_E.append(slope)
            y_int.append(yint)


        from parallel_decorators import vectorize_parallel

        @vectorize_parallel(method='processes', num_procs=11)
        def slice_normalization(i, esorted, slope, intercept):

            slice_slope = slope[i]
            slice_intercept = intercept[i]
            jz_max = (esorted[isplit[i]:isplit[i+1]]) - slice_intercept / slice_slope
            return jz_max

        results = slice_normalization(arange(len(isplit[0:-1])), E_sorted, slope_E, y_int)
        results = np.concatenate(results).ravel()
        results = np.append(results, (E_sorted[isplit[-1]] - y_int[-1]) / slope_E[-1])
        eps2 = jz_sorted / results
        '''end'''
    
        '''plotting'''
        idx=0
        while idx<len(isplit)-1:
            # plt.plot(jz_sorted[i*nselect:(i+1)*nselect], E_sorted[i*nselect:(i+1)*nselect], '.', ms=1, alpha=0.6)
            plt.plot(jz_sorted[isplit[idx]:isplit[idx+1]], E_sorted[isplit[idx]:isplit[idx+1]], '.', ms=1, alpha=.6)
            idx+=1

        plt.plot(results[isplit],E_sorted[isplit],'or',lw=2,label='Esitmation')
        plt.xlabel('JZ',fontsize=13)
        plt.ylabel('Spec Energy',fontsize=13)
        plt.title('Z Component Angular Velocity vs Spec Energy Snap'+str(i), fontsize=15, pad=15)
        # plt.legend()

        plt.savefig('halo'+str(select_halo)+'_snap'+str(i)+'.png')
    #     Z.append(sf.redshift)
    #     Median_Eps.append(median(eps2))
    
    # d = {'Redshift': Z, 'Median Epsilon': Median_Eps}
    # df = pd.DataFrame(data=d)
    # filename = '/home/lgvanover/intro_scripts/halos/new/level'+str(sim_level)+'_halo'+str(select_halo)+'.hdf5'
    # print('saving to', filename)
    # df.to_hdf(filename, 'table')
    # return df

if __name__ == '__main__':

    halo = sys.argv[1]
    level = sys.argv[2]
    print('loading halo '+str(halo)+' level '+str(level))
    df = load_halo(int(halo), int(level))
