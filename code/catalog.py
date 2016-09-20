''' 
Construct catalogs for conformity measurements 

'''
import numpy as np 
import h5py
import pickle
import time 
import scipy.spatial as scispace
import astropy.cosmology as astrocosmo
from pydl.pydlutils.spheregroup import spherematch

# --- Local ---  
import util as UT
from ChangTools.fitstables import mrdfits


class ConformCatalog(object):
    def __init__(self, Mrcut=18, 
            primary_delv=500., primary_rperp=0.5, 
            neighbor_delv=1000., neighbor_rperp=5.):
        ''' Galaxy Catalog with primaries and secondaries identified 
        for conformity measurements
        '''
        self.Mrcut = Mrcut      # absolute magnitude cutoff in Tinker et al. (2011) sample
        self.M_cut = Tinker_Masscut(self.Mrcut) # mass cut
    
        # Delta v and r_perp for primary identification
        self.primary_delv = primary_delv
        self.primary_rperp = primary_rperp

        # Delta v and r_perp for neighbor identification
        self.neighbor_delv = neighbor_delv
        self.neighbor_rperp = neighbor_rperp

        #def ReadJackknife(self, n_jack, RADec_bins=[5,5]): 
        #    ''' Read jackknife catalog 
        #    '''
        #    jack_file = ''.join([(self.File()).rsplit('.p',1)[0], '.jackknife', 
        #        str(n_jack), 'of', str(RADec_bins[0]), 'x', str(RADec_bins[1]), '.p'])
        #    catalog = pickle.load(open(jack_file, 'rb'))
        #    return catalog

    def Jackknife(self, catalog, n_jack, RADec_bins=[5,5]): 
        ''' Remove n_jack jackknife field from the catalog. 
        Numbering of the jackknife field goes across Dec first. 
        e.g. RADec_bins=[3,4] 

                Dec 
            1 | 2 | 3 | 4
        RA  5 | 6 | 7 | 8 
            9 | 10| 11| 12 

        Parameters
        ----------
        * n_jack : int
            jackknife field number 
        * RADec_bins : list
            2 element list that specifies the total number of RA and 
            Dec bins.
        '''
        if n_jack > (RADec_bins[0] * RADec_bins[1]): 
            raise ValueError
        #catalog = self.Read()

        RA_percents = [100./np.float(RADec_bins[0])*np.float(i) for i in range(1,RADec_bins[0])]
        Dec_percents = [100./np.float(RADec_bins[0])*np.float(i) for i in range(1,RADec_bins[0])]

        RA_limits = np.percentile(catalog['ra'], RA_percents)
        Dec_limits = np.percentile(catalog['dec'], Dec_percents)

        # find RA and Dec limits of jackknife bin 
        i_dec = np.mod(n_jack, RADec_bins[1])
        if i_dec == 0: 
            i_dec = RADec_bins[1]
            i_ra = (n_jack - np.mod(n_jack, RADec_bins[1]))/RADec_bins[1]
        else: 
            i_ra = (n_jack - np.mod(n_jack, RADec_bins[1]))/RADec_bins[1]+1
        
        # jackknife conditions for all galaxies 
        if i_ra == 1: 
            cut_ra = catalog['ra'] > RA_limits[0]
        elif i_ra == RADec_bins[0]: 
            cut_ra = catalog['ra'] < RA_limits[-1]
        else: 
            cut_ra = (catalog['ra'] < RA_limits[i_ra-2]) | (catalog['ra'] > RA_limits[i_ra-1])

        if i_dec == 1: 
            cut_dec = catalog['dec'] > Dec_limits[0] 
        elif i_dec == RADec_bins[1]: 
            cut_dec = catalog['dec'] < Dec_limits[-1] 
        else: 
            cut_dec = (catalog['dec'] < Dec_limits[i_dec-2]) | \
                    (catalog['dec'] > Dec_limits[i_dec-1])
        cut_jack = np.where(cut_ra | cut_dec)[0]
        
        n_gal = len(catalog['ra'])
        jack_catalog = {} 
        for key in catalog.keys(): 
            if isinstance(catalog[key], list): 
                jack_catalog[key] = [catalog[key][i_cut_jack] for i_cut_jack in cut_jack]
            else: 
                jack_catalog[key] = catalog[key][cut_jack]
        #jack_file = ''.join([(self.File()).rsplit('.p',1)[0], '.jackknife', 
        #    str(n_jack), 'of', str(RADec_bins[0]), 'x', str(RADec_bins[1]), '.p'])
        #print jack_file
        #pickle.dump(catalog, open(jack_file, 'wb')) 
        return jack_catalog 

    def Read(self): 
        ''' Read in the conformity catalog
        '''
        catalog = pickle.load(open(self.File(), 'rb'))
        return catalog

    def Build(self, clobber=False): 
        ''' Build conformity catalog. 

        1. identify primaries using VAGC values
        2. identify primaries using MPA-JHU values 
        3. identify neighbors of VAGC primaries 
        4. identify neighbors of MPA-JHU primaries 
        5. save to pickle file
        '''
        if clobber: 
            Build_TinkerCatalog(Mrcut=self.Mrcut)
            Build_MPAJHU_TinkerCatalog(Mrcut=self.Mrcut)
        catalog = MPAJHU_TinkerCatalog(Mrcut=self.Mrcut)

        # identify primaries based on VAGC values
        catalog = IdentifyPrimaries(catalog, 
                del_v_cut=self.primary_delv, r_perp_cut=self.primary_rperp, 
                mpajhu=False)
        # identify primaries based on MPA-JHU values
        catalog = IdentifyPrimaries(catalog, 
                del_v_cut=self.primary_delv, r_perp_cut=self.primary_rperp, 
                mpajhu=True)

        # identify neighbors of VAGC primaries 
        catalog = IdentifyNeighbors(catalog, 
                del_v_cut=self.neighbor_delv, r_perp_cut=self.neighbor_rperp, 
                mpajhu=False)
        # identify neighbors of MPAJHU primaries 
        catalog = IdentifyNeighbors(catalog, 
                del_v_cut=self.neighbor_delv, r_perp_cut=self.neighbor_rperp, 
                mpajhu=True)
        pickle.dump(catalog, open(self.File(), 'wb')) 
        return None 

    def File(self):  
        ''' Conformity catalog file name 
        '''
        conform_file = ''.join([
            UT.dir_dat(), 'conform_catalog/',
            'MPAJHU_TinkerGroupCat.Mr', str(self.Mrcut), '.Mass', str(self.M_cut), 
            self._FileSpec(), '.p']) 
        return conform_file 

    def _FileSpec(self):
        # string that specifies the choices of the primaries and neighbors 
        spec_str = ''.join([
            '.primary', 
            '_delv', str(self.primary_delv),
            '_rperp', str(self.primary_rperp), 
            '.neighbor', 
            '_delv', str(self.neighbor_delv), 
            '_rperp', str(self.neighbor_rperp)])
        return spec_str


def Tinker_Masscut(Mrcut): 
    # given Mr cut return the matching M_cut 
    if Mrcut == 18: 
        return 9.4
    elif Mrcut == 19:
        return 9.8
    elif Mrcut == 20: 
        return 10.2 
    else: 
        raise ValueError


def TinkerCatalog(Mrcut=18): 
    ''' Tinker et al. (2011) group catalog combined into a 
    volume-limited galaxy catalog and return a dictionary with
    all the value. 
    '''
    M_cut = Tinker_Masscut(Mrcut) 
    # read in h5py file 
    tinker_file = ''.join([UT.dir_dat(), 'tinker2011catalogs/',
        'GroupCat.Mr', str(Mrcut), '.Mass', str(M_cut), '.D360.hdf5']) 
    
    catalog = {} 
    f = h5py.File(tinker_file, 'r')
    grp = f['data']
    for col in grp.keys(): 
        catalog[col] = grp[col].value  

    f.close() 
    return catalog


def Build_TinkerCatalog(Mrcut=18): 
    ''' Preprocess the group catalog data into a more python friendly format
    with appropriate *little h* corrections!
    '''
    h = 0.7
    M_cut = Tinker_Masscut(Mrcut) 
    # Read Group Catalog GalData 
    galdata_file = ''.join([UT.dir_dat(), 'tinker2011catalogs/', 
        'clf_groups_M', str(Mrcut), '_', str(M_cut), '_D360.', 'galdata_corr.fits']) 
    gal_data = mrdfits(galdata_file) 

    catalog = {} 
    for column in gal_data.__dict__.keys(): 
        column_data = getattr(gal_data, column)
        if column == 'stellmass':       
            # stellmass is in units of Msol/h^2
            # why jeremy why?/?
            column_data = column_data / h**2
            catalog['mass'] = np.log10(column_data)     # convert to log Mass
        elif column == 'ssfr': 
            column_data += np.log10(h**2)    # little h #*(@*#$
            catalog['ssfr'] = column_data
        elif column == 'cz': # convert to z  
            catalog['z'] = column_data/299792.458
        elif column in ['ra', 'dec']: 
            catalog[column] = column_data * 57.2957795
        else: 
            catalog[column] = column_data 
    catalog['sfr'] = catalog['mass'] + catalog['ssfr']   # calculate SFR form mass and ssfr
    
    # Read Group Catalog probability 
    prob_file = ''.join([UT.dir_dat(), '/tinker2011catalogs/', 
        'clf_groups_M', str(Mrcut), '_', str(M_cut), '_D360.', 'prob.fits']) 
    prob_data = mrdfits(prob_file)            # import probability file 
    for column in prob_data.__dict__.keys(): 
        catalog[column] = getattr(prob_data, column) 
    
    tinker_file = ''.join([UT.dir_dat(), 'tinker2011catalogs/',
        'GroupCat.Mr', str(Mrcut), '.Mass', str(M_cut), '.D360.hdf5']) 
    
    f = h5py.File(tinker_file, 'w')
    grp = f.create_group('data')
    for key in catalog.keys(): 
        grp.create_dataset(key, data=catalog[key])

    f.close() 
    return None 


def MPAJHU_TinkerCatalog(Mrcut=18): 
    ''' Read in the Tinker et al. (2011) matched up to the MPA-JHU catalog
    and return catalog dictionary 
    '''
    M_cut = Tinker_Masscut(Mrcut) 
    mpajhu_tinker_file = ''.join([UT.dir_dat(), 
        'tinker2011catalogs/',
        'GroupCat.Mr', str(Mrcut), '.Mass', str(M_cut), 
        '.D360.MPAJHU.hdf5']) 

    catalog = {} 
    f = h5py.File(mpajhu_tinker_file, 'r')
    grp = f['data']
    for col in grp.keys(): 
        catalog[col] = grp[col].value  

    f.close() 
    return catalog


def Build_MPAJHU_TinkerCatalog(Mrcut=18): 
    ''' Append MPA-JHU SSFR values to the Tinker et al. (2011) catalog.
    The main purpose is to try to reproduce the Kauffmann et al. (2013) results. 
    Galaxies are matched to each other through spherematch. 
    '''
    # import Tinker et al. (2011) catalog with specified Mr cut  
    catalog = TinkerCatalog(Mrcut=Mrcut) 
    
    # import MPA-JHU catalog
    mpajhu_gals = mrdfits(''.join([UT.dir_dat(), 'mpa_jhu/', 'gal_info_dr7_v5_2.fit'])) 
    # SFR total
    mpajhu_sfrtot = mrdfits(''.join([UT.dir_dat(), 'mpa_jhu/', 'gal_totsfr_dr7_v5_2.fits']))
    # SFR fiber
    mpajhu_sfrfib = mrdfits(''.join([UT.dir_dat(), 'mpa_jhu/', 'gal_fibsfr_dr7_v5_2.fits']))
    # SSFR total 
    mpajhu_ssfrtot = mrdfits(''.join([UT.dir_dat(), 'mpa_jhu/', 'gal_totspecsfr_dr7_v5_2.fits']))
    # SSFR fiber
    mpajhu_ssfrfib = mrdfits(''.join([UT.dir_dat(), 'mpa_jhu/', 'gal_fibspecsfr_dr7_v5_2.fits']))
    # stellar mass total 
    mpajhu_masstot = mrdfits(''.join([UT.dir_dat(), 'mpa_jhu/', 'totlgm_dr7_v5_2.fit']))
    # stellar mass fiber 
    mpajhu_massfib = mrdfits(''.join([UT.dir_dat(), 'mpa_jhu/', 'fiblgm_dr7_v5_2.fit']))
    
    t_spherematch = time.time() 
    match = spherematch(
            catalog['ra'], catalog['dec'], mpajhu_gals.ra, mpajhu_gals.dec, 0.000833333) 
    print 'Spherematch with matchlenght = ', 0.000833333
    print 'takes ', time.time() - t_spherematch, 'seconds' 
    print 1.- np.float(len(match[0]))/np.float(len(catalog['ra'])), 'of the VAGC galaxies'
    print 'do not have matches, likely due to fiber collisions'
    if len(match[0]) != len(np.unique(match[0])): 
        raise ValueError
    
    # save the MPAJHU indices, jsut in case
    catalog['mpajhu_index'] = np.repeat(-999, len(catalog['ra'])) 
    catalog['mpajhu_index'][match[0]] = match[1]
    
    # append SFR, SSFR, and mass values to catalog 
    for col in [
            'sfr_tot_mpajhu', 'sfr_fib_mpajhu', 
            'ssfr_tot_mpajhu', 'ssfr_fib_mpajhu', 
            'mass_tot_mpajhu', 'mass_fib_mpajhu']:   # initiate arrays
        catalog[col] = np.repeat(-999., len(catalog['ra']))

    catalog['sfr_tot_mpajhu'][match[0]] = mpajhu_sfrtot.median[match[1]]
    catalog['sfr_fib_mpajhu'][match[0]] = mpajhu_sfrfib.median[match[1]]
    catalog['ssfr_tot_mpajhu'][match[0]] = mpajhu_ssfrtot.median[match[1]]
    catalog['ssfr_fib_mpajhu'][match[0]] = mpajhu_ssfrfib.median[match[1]] 
    catalog['mass_tot_mpajhu'][match[0]] = mpajhu_masstot.median[match[1]] 
    catalog['mass_fib_mpajhu'][match[0]] = mpajhu_massfib.median[match[1]] 

    # trim galaxies without matches
    hasmatch = np.where(catalog['mpajhu_index'] != -999)  
    for key in catalog.keys(): 
        key_val = catalog[key]
        catalog[key] = key_val[hasmatch]
    catalog['mpajhu_tinker_index'] = hasmatch[0] 
    
    M_cut = Tinker_Masscut(Mrcut) 
    mpajhu_tinker_file = ''.join([UT.dir_dat(), 
        'tinker2011catalogs/',
        'GroupCat.Mr', str(Mrcut), '.Mass', str(M_cut), 
        '.D360.MPAJHU.hdf5']) 

    f = h5py.File(mpajhu_tinker_file, 'w')
    grp = f.create_group('data')
    for key in catalog.keys(): 
        grp.create_dataset(key, data=catalog[key])

    f.close() 
    return None


def IdentifyPrimaries(catalog, del_v_cut=500., r_perp_cut=0.5, mpajhu=False): 
    ''' Identify the primary galaxies in the input catalog dictionary using
    the following steps: 

    1. Run a KDTree and identify the neighbors within r_eff of the 'target' galaxy. 
        r_eff = sqrt( r_perp^2 + d_delv,max^2 )
        
    2. Go through each 'target' galaxy and cut by del v and r_perp and impose 
        stellar mass criteria to identify the primaries. 

    Notes
    -----
    * If mpajhu=True then the MPAJHU stellar mass M_* total is used instead of the
        VAGC M* from Tinker et al. (2011) catalog. 
    '''
    c_kms = 299792.458  # speed of light
    # same cosmology as Kauffmann et al. (2013) 
    cosmo_kauff = astrocosmo.FlatLambdaCDM(H0=70., Om0=0.3)     
    
    # x, y, z coordinates of 
    if 'xyz' not in catalog.keys():
        catalog['xyz'] = UT.radecz_to_xyz(catalog['ra'], catalog['dec'], catalog['z'], 
                H0=cosmo_kauff.H0.value, Om0=cosmo_kauff.Om0)

    d_delv_max = cosmo_kauff.comoving_distance(del_v_cut/c_kms).value
    r_eff = np.sqrt(d_delv_max**2 + r_perp_cut**2)
    print 'd_delv_max', d_delv_max
    print 'KDTree radius', r_eff, 'Mpc'
    
    # construct KDTree
    kd_time = time.time() 
    cat_tree = scispace.KDTree(catalog['xyz'])
    kdt_satellite_indices = cat_tree.query_ball_tree(cat_tree, r_eff)
    print 'KDTree takes ', time.time() - kd_time, ' seconds'

    sat_indices = [[] for i in range(len(catalog['z']))] # satellites of primary isolation criteria
    sat_rperp = [[] for i in range(len(catalog['z']))]   # r_perp of satellites 
    n_sat = np.zeros(len(catalog['z'])) # number of satellites  
    isprimary = np.zeros(len(catalog['z']))  
    for i_targ in range(len(catalog['z'])): 
        # calculate Delta v between target and satellites
        del_v_pair = c_kms * (
                catalog['z'][i_targ] - catalog['z'][kdt_satellite_indices[i_targ]]
                )

        # calculate r_perp
        targ_xyz = catalog['xyz'][i_targ]
        sat_xyz = catalog['xyz'][kdt_satellite_indices[i_targ]]
        targ_mag = np.sum(targ_xyz**2) 
        targ_dot_env = np.sum(targ_xyz * sat_xyz, axis=1)
        proj = targ_dot_env / targ_mag
        proj_xyz = np.array([proj[i] * targ_xyz for i in xrange(len(proj))])

        rll = np.sqrt(np.sum((proj_xyz - targ_xyz)**2, axis=1))
        rperp = np.sqrt(np.sum((proj_xyz - sat_xyz)**2, axis=1))

        # So that we don't count the target
        dr = np.sqrt(np.sum((sat_xyz - targ_xyz)**2, axis=1))

        # mass criteria
        if not mpajhu:  
            M_targ = catalog['mass'][i_targ]
            M_sat = catalog['mass'][kdt_satellite_indices[i_targ]]
        else: 
            M_targ = catalog['mass_tot_mpajhu'][i_targ]
            M_sat = catalog['mass_tot_mpajhu'][kdt_satellite_indices[i_targ]]

        keep_sat = np.where(
            (del_v_pair < del_v_cut) & 
            (rperp < r_perp_cut) & 
            (dr > 0.)
            ) 
    
        if len(keep_sat[0]) > 0: 
            n_sat[i_targ] = len(keep_sat[0])
            sat_rperp[i_targ] = list(rperp[keep_sat])
            sat_indices[i_targ] = list(np.array(kdt_satellite_indices[i_targ])[keep_sat])

            if M_sat[keep_sat].max() < (M_targ + np.log10(0.5)):    
                # More massive than twice the most massive neighbor 
                isprimary[i_targ] = 1
        else: 
            # is by default a primary 
            isprimary[i_targ] = 1
    
    if not mpajhu: 
        catalog['n_secondary_vagc'] = n_sat
        catalog['secondary_rperp_vagc'] = sat_rperp
        catalog['secondary_indices_vagc'] = sat_indices
        catalog['primary_vagc'] = np.array(isprimary)
    else:
        catalog['n_secondary_mpajhu'] = n_sat
        catalog['secondary_rperp_mpajhu'] = sat_rperp
        catalog['secondary_indices_mpahju'] = sat_indices
        catalog['primary_mpajhu'] = np.array(isprimary)

    print 'Identified ', np.sum(isprimary), ' primaries' 
    return catalog 


def IdentifyNeighbors(catalog, del_v_cut=500., r_perp_cut=5., mpajhu=False): 
    ''' Identify the neighboring galaxies of primary galaxies in input 
    catalog dictionary using the same steps as for primaries but for different 
    sets of delta v and r_perp cut offs

    1. Run a KDTree and identify the neighbors within r_eff of the 'target' galaxy. 
        r_eff = sqrt( r_perp^2 + d_delv,max^2 )
    '''
    c_kms = 299792.458  # speed of light
    # same cosmology as Kauffmann et al. (2013) 
    cosmo_kauff = astrocosmo.FlatLambdaCDM(H0=70., Om0=0.3)     
    
    if not mpajhu:
        if 'primary_vagc' not in catalog.keys():
            raise ValueError
    else:
        if 'primary_mpajhu' not in catalog.keys():
            raise ValueError

    d_delv_max = cosmo_kauff.comoving_distance(del_v_cut/c_kms).value
    r_eff = np.sqrt(d_delv_max**2 + r_perp_cut**2)
    print 'd_delv_max', d_delv_max
    print 'KDTree radius', r_eff, 'Mpc'
    
    # loop through primaries and identify their neighbors 
    if not mpajhu: 
        is_primary = np.where(catalog['primary_vagc'] == 1)[0]
    else:
        is_primary = np.where(catalog['primary_mpajhu'] == 1)[0]
    primary_xyz = catalog['xyz'][is_primary]

    # construct KDTree
    kd_time = time.time() 
    cat_tree = scispace.KDTree(catalog['xyz'])
    primary_tree = scispace.KDTree(primary_xyz)
    
    kdt_neighbor_indices = primary_tree.query_ball_tree(cat_tree, r_eff)
    print 'KDTree takes ', time.time() - kd_time, ' seconds'

    neigh_indices = [[] for i in range(len(catalog['z']))]  # satellites of primary isolation criteria
    neigh_rperp = [[] for i in range(len(catalog['z']))]    # r_perp of satellites 
    n_neigh = np.zeros(len(catalog['z']))        # number of satellites  
    for i_targ in range(len(is_primary)): 
        # calculate Delta v between target and satellites
        del_v_pair = c_kms * (
                catalog['z'][is_primary[i_targ]] - catalog['z'][kdt_neighbor_indices[i_targ]]
                )

        # calculate r_perp
        targ_xyz = primary_xyz[i_targ]
        sat_xyz = catalog['xyz'][kdt_neighbor_indices[i_targ]]
        targ_mag = np.sum(targ_xyz**2) 
        targ_dot_env = np.sum(targ_xyz * sat_xyz, axis=1)
        proj = targ_dot_env / targ_mag
        proj_xyz = np.array([proj[i] * targ_xyz for i in xrange(len(proj))])

        rll = np.sqrt(np.sum((proj_xyz - targ_xyz)**2, axis=1))
        rperp = np.sqrt(np.sum((proj_xyz - sat_xyz)**2, axis=1))

        # So that we don't count the target
        dr = np.sqrt(np.sum((sat_xyz - targ_xyz)**2, axis=1))

        keep_neigh = np.where(
            (del_v_pair < del_v_cut) & 
            (rperp < r_perp_cut) & 
            (dr > 0.) 
            )
        if not mpajhu: 
            if len(keep_neigh[0]) < catalog['n_secondary_vagc'][is_primary[i_targ]]: 
                raise ValueError
        else:
            if len(keep_neigh[0]) < catalog['n_secondary_mpajhu'][is_primary[i_targ]]: 
                raise ValueError
        n_neigh[is_primary[i_targ]] = len(keep_neigh[0])
        neigh_rperp[is_primary[i_targ]] = list(rperp[keep_neigh])
        neigh_indices[is_primary[i_targ]] = list(np.array(kdt_neighbor_indices[i_targ])[keep_neigh])

    neigh_primary = np.repeat(-999, len(catalog['z'])) 
    neigh_primary[is_primary] = is_primary 
    
    if not mpajhu: 
        catalog['neighbor_primary_vagc'] = neigh_primary  
        catalog['n_neighbor_vagc'] = n_neigh
        catalog['neighbor_rperp_vagc'] = neigh_rperp
        catalog['neighbor_indices_vagc'] = neigh_indices
    else: 
        catalog['neighbor_primary_mpajhu'] = neigh_primary
        catalog['n_neighbor_mpajhu'] = n_neigh
        catalog['neighbor_rperp_mpajhu'] = neigh_rperp
        catalog['neighbor_indices_mpajhu'] = neigh_indices

    return catalog 



"""
    def Jackknife_TinkerCatalog(n_jack, Mrcut=18): 
        ''' Jackknife catalogs of Tinker et al. (2011) group catalog combined into a 
        volume-limited galaxy catalog and return a dictionary with
        all the value. 
        '''
        M_cut = Tinker_Masscut(Mrcut)
        # read in h5py file 
        tinker_file = lambda censat: ''.join([UT.dir_dat(), 'tinker2011catalogs/',
                    'GroupCat.Mr', str(Mrcut), '.Mass', str(M_cut), '.D360.', censat, '.hdf5']) 
        tinker_central = h5py.File(tinker_file('central'), 'r')
        tinker_satellite = h5py.File(tinker_file('satellite'), 'r')

        # now combine the two into a dictionary
        if tinker_central['data'].keys() != tinker_satellite['data'].keys(): 
            # make sure that they have the same columns
            raise ValueError
        tinker_cat = {}
        for col in tinker_central['data'].keys(): 
            tinker_cat[col] = np.concatenate(
                    [tinker_central['data/'+col].value, tinker_satellite['data/'+col].value]
                    )
        #tinker_cat['ra'].min(), tinker_cat['ra'].max()
        #tinker_cat['ra'].min(), tinker_cat['ra'].max()
        #for col in tinker_cat.keys(): 
        return tinker_cat 
"""
