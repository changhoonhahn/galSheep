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
    def __init__(self, catalog_name, catalog_prop={},  
            primary_delv=500., primary_rperp=0.5, 
            neighbor_delv=1000., neighbor_rperp=5.):
        ''' Galaxy Catalog with primaries and secondaries identified 
        for conformity measurements
        '''
        self.catalog_name = catalog_name 
        if self.catalog_name == 'tinker': 
            self.Mrcut = catalog_prop['Mrcut']      # absolute magnitude cutoff in Tinker et al. (2011) sample
            self.M_cut = Tinker_Masscut(self.Mrcut) # mass cut
        elif self.catalog_name == 'tinkauff': 
            self.M_cut = catalog_prop['Mass_cut']
        elif self.catalog_name == 'tinkauff_iso': 
            self.M_cut = catalog_prop['Mass_cut']
        elif self.catalog_name == 'kauff': 
            pass
    
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

    def zSubsample(self, catalog, lowhigh): 
        ''' Divide the catalog into redshift subsample: 
        below the median, above the median
        '''
        if 'z' not in catalog.keys(): 
            raise ValueError
        med_z = np.median(catalog['z'])

        if lowhigh == 'low': 
            z_cut = np.where(catalog['z'] < med_z)[0]
        elif lowhigh == 'high': 
            z_cut = np.where(catalog['z'] >= med_z)[0] 
        else: 
            raise ValueError
        
        z_catalog = {} 
        for key in catalog.keys(): 
            if isinstance(catalog[key], list): 
                z_catalog[key] = [catalog[key][i_zcut] for i_zcut in z_cut]
            else: 
                z_catalog[key] = catalog[key][z_cut]
        return z_catalog, med_z

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
        if self.catalog_name == 'tinker': 
            if clobber: 
                Build_TinkerCatalog(Mrcut=self.Mrcut)
                Build_MPAJHU_TinkerCatalog(Mrcut=self.Mrcut)
            catalog = MPAJHU_TinkerCatalog(Mrcut=self.Mrcut)
            
            # identify primaries based on VAGC values
            vagc_primary_info = IdentifyPrimaries(catalog, masses=catalog['mass'], 
                    del_v_cut=self.primary_delv, r_perp_cut=self.primary_rperp)
            for key in vagc_primary_info.key(): 
                catalog[key+'_vagc'] = vagc_primary_info[key]

            # identify primaries based on MPAJHU values
            mpajhu_primary_info = IdentifyPrimaries(catalog, masses=catalog['mass_tot_mpajhu'],
                    del_v_cut=self.primary_delv, r_perp_cut=self.primary_rperp)
            for key in mpajhu_primary_info.key(): 
                catalog[key+'_mpajhu'] = mpajhu_primary_info[key]

            # identify neighbors of VAGC primaries 
            vagc_neighbor_info = IdentifyNeighbors(catalog, primary_info=vagc_primary_info,
                    del_v_cut=self.neighbor_delv, r_perp_cut=self.neighbor_rperp)
            for key in vagc_neighbor_info.key(): 
                catalog[key+'_vagc'] = vagc_neighbor_info[key]

            # identify neighbors of MPAJHU primaries 
            mpajhu_neighbor_info = IdentifyNeighbors(catalog, primary_info=mpajhu_primary_info,
                    del_v_cut=self.neighbor_delv, r_perp_cut=self.neighbor_rperp)
            for key in mpajhu_neighbor_info.key(): 
                catalog[key+'_vagc'] = mpajhu_neighbor_info[key]

        elif self.catalog_name in ['tinkauff', 'tinkauff_iso', 'kauff']: 
            if self.catalog_name == 'tinkauff':
                if clobber: 
                    Build_TinKauffGroupCat(Mass_cut=self.M_cut)
                catalog = TinKauffGroupCat(Mass_cut=self.M_cut)
            elif self.catalog_name == 'tinkauff_iso':
                if clobber: 
                    Build_TinKauff_IsolationGroupCat(Mass_cut=self.M_cut)
                catalog = TinKauff_IsolationGroupCat(Mass_cut=self.M_cut)
            elif self.catalog_name == 'kauff': 
                if clobber: 
                    Build_KauffmannParent()
                catalog = KauffmannParent()
            
            # identify primaries based on MPAJHU values
            mpajhu_primary_info = IdentifyPrimaries(catalog, masses=catalog['mass_tot_mpajhu'], 
                    del_v_cut=self.primary_delv, r_perp_cut=self.primary_rperp)
            for key in mpajhu_primary_info.keys(): 
                catalog[key+'_mpajhu'] = mpajhu_primary_info[key]

            # identify neighbors of MPAJHU primaries 
            mpajhu_neighbor_info = IdentifyNeighbors(catalog, primary_info=mpajhu_primary_info,
                    del_v_cut=self.neighbor_delv, r_perp_cut=self.neighbor_rperp)
            for key in mpajhu_neighbor_info.keys(): 
                catalog[key+'_mpajhu'] = mpajhu_neighbor_info[key]
        
        pickle.dump(catalog, open(self.File(), 'wb')) 
        return None 

    def File(self):  
        ''' Conformity catalog file name 
        '''
        if self.catalog_name == 'tinker': 
            conform_file = ''.join([
                UT.dir_dat(), 'conform_catalog/',
                'MPAJHU_TinkerGroupCat.Mr', str(self.Mrcut), '.Mass', str(self.M_cut), 
                self._FileSpec(), '.p']) 
        elif self.catalog_name == 'tinkauff':
            conform_file = ''.join([
                UT.dir_dat(), 'conform_catalog/',
                'VAGCdr72brigh34_MPAJHU.GroupCat.Mass', str(self.M_cut), 
                self._FileSpec(), '.p']) 
        elif self.catalog_name == 'tinkauff_iso':
            conform_file = ''.join([
                UT.dir_dat(), 'conform_catalog/',
                'VAGCdr72brigh34_MPAJHU.IsoGroupCat.Mass', str(self.M_cut), 
                self._FileSpec(), '.p']) 
        elif self.catalog_name == 'kauff':
            conform_file = ''.join([
                UT.dir_dat(), 'conform_catalog/',
                'VAGCdr72brigh34_MPAJHU.Kauffmann', self._FileSpec(), '.p']) 
        else: 
            raise ValueError

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


def KauffmannParent(): 
    ''' Read in the the Kauffmannn et al. (2013) parent sample constructed from 
    VAGC dr72bright34 catalog 
    '''
    dr72_file = ''.join([UT.dir_dat(), 'vagc/', 'VAGCdr72.Kauff2013cut.hdf5']) 
    catalog = {} 
    f = h5py.File(dr72_file, 'r')
    grp = f['data']
    for col in grp.keys(): 
        catalog[col] = grp[col].value  
    f.close() 
    return catalog


def Build_KauffmannParent(): 
    ''' Try to create the parent sample of Kauffmann et al.(2013) 
    '''
    # import VAGC dr72bright34
    vagc_dr72 = VAGCdr72bright34_Catalog() 

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
    
    catalog = {} 
    catalog['ra'] = vagc_dr72['ra']
    catalog['dec'] = vagc_dr72['dec']
    catalog['z'] = vagc_dr72['z']
    for i_band, band in enumerate(['u', 'g', 'r', 'i', 'z']):
        catalog['M_'+band] = vagc_dr72['M_'+band]

    # pre cut 
    cut_z = (catalog['z'] > 0.017) & (catalog['z'] < 0.03) 
    pre_cuts = np.where(cut_z)#& cut_stellarmass & cut_absmag)
    for key in catalog.keys(): 
        catalog[key] = catalog[key][pre_cuts]

    t_spherematch = time.time() 
    match = spherematch(catalog['ra'], catalog['dec'], mpajhu_gals.ra, mpajhu_gals.dec, 0.000833333) 
    print 'Spherematch with matchlenght = ', 0.000833333
    print 'takes ', time.time() - t_spherematch, 'seconds' 
    print 1.- np.float(len(match[0]))/np.float(len(catalog['ra'])), 'of the VAGC galaxies'
    print 'do not have matches'
    if len(match[0]) != len(np.unique(match[0])): 
        raise ValueError
    
    # save the MPAJHU indices, jsut in case
    catalog['mpajhu_index'] = np.repeat(-999, len(catalog['ra'])) 
    catalog['mpajhu_index'][match[0]] = match[1]

    # append SFR, SSFR, and mass values to catalog 
    for col in ['sfr_tot_mpajhu', 'sfr_fib_mpajhu', 
            'ssfr_tot_mpajhu', 'ssfr_fib_mpajhu', 
            'mass_tot_mpajhu', 'mass_fib_mpajhu']:   # initiate arrays
        catalog[col] = np.repeat(-999., len(catalog['ra']))

    catalog['sfr_tot_mpajhu'][match[0]] = mpajhu_sfrtot.median[match[1]]
    catalog['sfr_fib_mpajhu'][match[0]] = mpajhu_sfrfib.median[match[1]]
    catalog['ssfr_tot_mpajhu'][match[0]] = mpajhu_ssfrtot.median[match[1]]
    catalog['ssfr_fib_mpajhu'][match[0]] = mpajhu_ssfrfib.median[match[1]]
    catalog['mass_tot_mpajhu'][match[0]] = mpajhu_masstot.median[match[1]]
    catalog['mass_fib_mpajhu'][match[0]] = mpajhu_massfib.median[match[1]]
    
    # kauffmann et al.(2013) cuts
    cut_stellarmass = (catalog['mass_tot_mpajhu'] > 9.25)
    cut_absmag = (catalog['M_r'] < -16.) & (catalog['M_r'] > -24.)
    cut_match = (catalog['mpajhu_index'] != -999)

    final_cuts = np.where(cut_stellarmass & cut_absmag & cut_match)
    for key in catalog.keys(): 
        catalog[key] = catalog[key][final_cuts]
    
    mpajhu_file = ''.join([UT.dir_dat(), 'vagc/', 'VAGCdr72.Kauff2013cut.hdf5']) 

    f = h5py.File(mpajhu_file, 'w')
    grp = f.create_group('data')
    for key in catalog.keys(): 
        grp.create_dataset(key, data=catalog[key])
    f.close() 
    return None


def Build_VAGCdr72_MPAJHU(Ascii=False):
    ''' Build VAGC dr72 with cross referenced MPAJHU stellar masses 
    and SSFRs.
    '''
    # import VAGC dr72bright34
    vagc_dr72 = VAGCdr72bright34_Catalog() 
    print len(vagc_dr72['ra']), ', VAGC dr72bright34 galaxies'

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
    
    catalog = {} 
    catalog['ra'] = vagc_dr72['ra']
    catalog['dec'] = vagc_dr72['dec']
    catalog['z'] = vagc_dr72['z']
    for i_band, band in enumerate(['u', 'g', 'r', 'i', 'z']):
        catalog['M_'+band] = vagc_dr72['M_'+band]

    t_spherematch = time.time() 
    match = spherematch(catalog['ra'], catalog['dec'], mpajhu_gals.ra, mpajhu_gals.dec, 0.000833333) 
    print 'Spherematch with matchlenght = ', 0.000833333
    print 'takes ', time.time() - t_spherematch, 'seconds' 
    print 1.- np.float(len(match[0]))/np.float(len(catalog['ra'])), 'of the VAGC galaxies'
    print 'do not have matches'
    if len(match[0]) != len(np.unique(match[0])): 
        raise ValueError
    
    # save the MPAJHU indices, jsut in case
    catalog['mpajhu_index'] = np.repeat(-999, len(catalog['ra'])) 
    catalog['mpajhu_index'][match[0]] = match[1]

    # append SFR, SSFR, and mass values to catalog 
    for col in ['sfr_tot', 'sfr_fib', 'ssfr_tot', 'ssfr_fib', 'mass_tot', 'mass_fib']:   # initiate arrays
        catalog[col] = np.repeat(-999., len(catalog['ra']))

    catalog['sfr_tot'][match[0]] = mpajhu_sfrtot.median[match[1]]
    catalog['sfr_fib'][match[0]] = mpajhu_sfrfib.median[match[1]]
    catalog['ssfr_tot'][match[0]] = mpajhu_ssfrtot.median[match[1]]
    catalog['ssfr_fib'][match[0]] = mpajhu_ssfrfib.median[match[1]]
    catalog['mass_tot'][match[0]] = mpajhu_masstot.median[match[1]]
    catalog['mass_fib'][match[0]] = mpajhu_massfib.median[match[1]]
    
    mpajhu_file = ''.join([UT.dir_dat(), 'vagc/', 'VAGCdr72.MPAJHU.nocut.hdf5']) 

    f = h5py.File(mpajhu_file, 'w')
    grp = f.create_group('data')
    for key in catalog.keys(): 
        grp.create_dataset(key, data=catalog[key])
    f.close() 

    if Ascii:   # write to Ascii (for jeremy) 
        mpajhu_file = ''.join([UT.dir_dat(), 'vagc/', 'VAGCdr72.MPAJHU.nocut.dat']) 
        column_order = ['ra', 'dec', 'z', 'mass_tot', 'sfr_tot' , 'ssfr_tot', 'mass_fib', 'sfr_fib', 'ssfr_fib']
        data_list = [] 
        data_fmt = ['%10.5f' for i in range(len(column_order))]
        str_header = ''
        for col in column_order: 
            data_list.append(catalog[col])
            if 'mass' in col: 
                str_header += ' '+col+' (Msun),'
            elif 'sfr' in col: 
                if 'ssfr' not in col: 
                    str_header += ' '+col+' (Msun/yr),'
                else: 
                    str_header += ' '+col+','
            else: 
                str_header += ' '+col+','
        np.savetxt(mpajhu_file, (np.vstack(np.array(data_list))).T, 
                fmt=data_fmt, delimiter='\t', header=str_header)
    return None


def VAGCdr72bright34_Catalog(): 
    ''' Read in the VAGC dr72bright34 catalog hdf5 file
    '''
    dr72_file = ''.join([UT.dir_dat(), 'vagc/', 'VAGCdr72bright34.hdf5']) 
    catalog = {} 
    f = h5py.File(dr72_file, 'r')
    grp = f['data']
    for col in grp.keys(): 
        catalog[col] = grp[col].value  
    f.close() 
    return catalog


def Build_VAGCdr72bright34(): 
    ''' Build hdf5 file of VAGC dr72brigh34, which is the
    parent sample of Jeremy's group catalog 
    '''
    # import VAGC dr72bright34
    vagc_photoinfo = np.loadtxt(''.join([UT.dir_dat(), 'vagc/', 'photoinfo_nonan.dr72bright34.dat']), 
            unpack=True, usecols=[0,1,2,3,4,5])
    vagc_lss = np.loadtxt(''.join([UT.dir_dat(), 'vagc/', 'lss.dr72bright34.dat']), 
            unpack=True, usecols=[0,3,4,5])
    if not np.array_equal(vagc_photoinfo[0], vagc_lss[0]):
        raise ValueError

    catalog = {}
    catalog['id'] = vagc_photoinfo[0]
    for i_band, band in enumerate(['u', 'g', 'r', 'i', 'z']): 
        catalog['M_'+band] = vagc_photoinfo[i_band+1]
    catalog['ra'] = vagc_lss[1]
    catalog['dec'] = vagc_lss[2]
    catalog['cz'] = vagc_lss[3]
    catalog['z'] = vagc_lss[3]/299792.458
    print len(catalog['z']), ' total galaxies'
    
    dr72_file = ''.join([UT.dir_dat(), 'vagc/', 'VAGCdr72bright34.hdf5']) 

    f = h5py.File(dr72_file, 'w')
    grp = f.create_group('data')
    for key in catalog.keys(): 
        grp.create_dataset(key, data=catalog[key])
    f.close() 
    return None


def Build_MPAJHU_TinkerCatalog_ASCII(Mrcut=18): 
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
    print mpajhu_massfib.median[match[1]]

    first_cols = ['id_gal', 'ra', 'dec', 'z', 'mass', 'sfr', 'ssfr', 
            'mass_tot_mpajhu', 'mass_fib_mpajhu', 'sfr_tot_mpajhu', 'sfr_fib_mpajhu', 'ssfr_tot_mpajhu', 'ssfr_fib_mpajhu'] 

    data_fmt = []
    data_list = [] 
    for i_key, key in enumerate(first_cols): 
        data_list.append(catalog[key]) 
        if key == 'id_gal': 
            data_fmt.append('%i')
        else: 
            data_fmt.append('%10.5f')

    later_cols = [] 
    for key in catalog.keys(): 
        if key not in first_cols: 
            later_cols.append(key)

    for key in later_cols: 
        data_list.append(catalog[key]) 
        if 'id' in key: 
            data_fmt.append('%i')
        elif 'index' in key: 
            data_fmt.append('%i')
        elif key == 'n_sersic': 
            data_fmt.append('%i')
        elif key == 'stellmass': 
            data_fmt.append('%1.5e')
        else: 
            data_fmt.append('%10.5f')


    str_header = ', '.join(first_cols + later_cols) 

    M_cut = Tinker_Masscut(Mrcut) 
    mpajhu_tinker_file = ''.join([UT.dir_dat(), 
        'tinker2011catalogs/',
        'GroupCat.Mr', str(Mrcut), '.Mass', str(M_cut), 
        '.D360.MPAJHU.dat']) 
    np.savetxt(mpajhu_tinker_file, (np.vstack(np.array(data_list))).T, 
            fmt=data_fmt, delimiter='\t', header=str_header)
    return None


def IdentifyPrimaries(catalog, masses=None, del_v_cut=500., r_perp_cut=0.5): 
    ''' Identify the primary galaxies in the input catalog dictionary using
    the following steps: 

    1. Run a KDTree and identify the neighbors within r_eff of the 'target' galaxy. 
        r_eff = sqrt( r_perp^2 + d_delv,max^2 )
        
    2. Go through each 'target' galaxy and cut by del v and r_perp and impose 
        stellar mass criteria to identify the primaries. 

    Notes
    -----
    '''
    c_kms = 299792.458  # speed of light
    # same cosmology as Kauffmann et al. (2013) 
    cosmo_kauff = astrocosmo.FlatLambdaCDM(H0=70., Om0=0.3)     
    
    # x, y, z coordinates of 
    if 'xyz' not in catalog.keys():
        catalog['xyz'] = UT.radecz_to_xyz(catalog['ra'], catalog['dec'], catalog['z'], 
                H0=cosmo_kauff.H0.value, Om0=cosmo_kauff.Om0)

    if masses is None: 
        raise ValueError 
    else: 
        if len(masses) != len(catalog['z']): 
            raise ValueError

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
        M_targ = masses[i_targ]
        M_sat = masses[kdt_satellite_indices[i_targ]]
        #if not mpajhu:  
        #    M_targ = catalog['mass'][i_targ]
        #    M_sat = catalog['mass'][kdt_satellite_indices[i_targ]]
        #else: 
        #    M_targ = catalog['mass_tot_mpajhu'][i_targ]
        #    M_sat = catalog['mass_tot_mpajhu'][kdt_satellite_indices[i_targ]]

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
   
    out_catalog = {} 
    out_catalog['n_secondary'] = n_sat
    out_catalog['secondary_rperp'] = sat_rperp
    out_catalog['secondary_indices'] = sat_indices
    out_catalog['primary'] = np.array(isprimary)
    print 'Identified ', np.sum(isprimary), ' primaries' 
    return out_catalog 


def IdentifyNeighbors(catalog, primary_info=None, del_v_cut=500., r_perp_cut=5.): 
    ''' Identify the neighboring galaxies of primary galaxies in input 
    catalog dictionary using the same steps as for primaries but for different 
    sets of delta v and r_perp cut offs

    1. Run a KDTree and identify the neighbors within r_eff of the 'target' galaxy. 
        r_eff = sqrt( r_perp^2 + d_delv,max^2 )
    '''
    c_kms = 299792.458  # speed of light
    # same cosmology as Kauffmann et al. (2013) 
    cosmo_kauff = astrocosmo.FlatLambdaCDM(H0=70., Om0=0.3)     
    
    if 'primary' not in primary_info.keys(): 
        raise ValueError
    if len(primary_info['primary']) != len(catalog['z']): 
        raise ValueError

    d_delv_max = cosmo_kauff.comoving_distance(del_v_cut/c_kms).value
    r_eff = np.sqrt(d_delv_max**2 + r_perp_cut**2)
    print 'd_delv_max', d_delv_max
    print 'KDTree radius', r_eff, 'Mpc'
    
    # loop through primaries and identify their neighbors 
    is_primary = np.where(primary_info['primary'] == 1)[0]
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
        if len(keep_neigh[0]) < primary_info['n_secondary'][is_primary[i_targ]]: 
            raise ValueError
        n_neigh[is_primary[i_targ]] = len(keep_neigh[0])
        neigh_rperp[is_primary[i_targ]] = list(rperp[keep_neigh])
        neigh_indices[is_primary[i_targ]] = list(np.array(kdt_neighbor_indices[i_targ])[keep_neigh])

    neigh_primary = np.repeat(-999, len(catalog['z'])) 
    neigh_primary[is_primary] = is_primary 
    
    out_catalog = {} 
    out_catalog['neighbor_primary'] = neigh_primary  
    out_catalog['n_neighbor'] = n_neigh
    out_catalog['neighbor_rperp'] = neigh_rperp
    out_catalog['neighbor_indices'] = neigh_indices
    return out_catalog 


def TinKauffGroupCat(Mass_cut=9.25): 
    ''' Read in the Tinker-Kauffmann Group catalog generated from VAGC 
    dr72bright34 with MPA-JHU galaxy property values. 
    '''
    dr72_file = ''.join([UT.dir_dat(), '/tinkauff/',
        'VAGCdr72_MPAJHU.GroupCat.Mass', str(Mass_cut), '.hdf5']) 
    catalog = {} 
    f = h5py.File(dr72_file, 'r')
    grp = f['data']
    for col in grp.keys(): 
        catalog[col] = grp[col].value  
    f.close() 
    return catalog


def Build_TinKauffGroupCat(Mass_cut=9.25): 
    ''' Compile the outputs of Jeremy's Group Catalog algorithm
    in order to generate catalogs analogous to Kauffmann et al.(2013). 
    Hence TinKauff. 
    ''' 
    # galdata_corr file 
    galdata_file = ''.join([UT.dir_dat(), 'tinkauff/', 
        'clf_groups_JHU_M', str(Mass_cut), '_z0.017_fibcoll.galdata_corr'])
    gal_data = np.loadtxt(galdata_file, unpack=True, usecols=range(1,14))

    catalog = {
            'id': gal_data[0], 
            'ra': gal_data[5] * 57.2957795, 
            'dec': gal_data[6] * 57.2957795,
            'M_r': gal_data[1], 
            'M_g': gal_data[2],
            'z': gal_data[3]/299792.458,
            'mass_tot_mpajhu': np.log10(gal_data[4]), 
            'Dn4000': gal_data[7], 
            'ssfr_tot_mpajhu': gal_data[8], 
            'ssfr_fib_mpajhu': gal_data[9], 
            'sfr_tot_mpajhu': gal_data[10], 
            'sfr_fib_mpajhu': gal_data[11], 
            'mass_fib_mpajhu': gal_data[12], 
            } 

    # prob data 
    prob_file = ''.join([UT.dir_dat(), 'tinkauff/', 
        'clf_groups_JHU_M', str(Mass_cut), '_z0.017_fibcoll.prob'])
    prob_data = np.loadtxt(prob_file, unpack=True, usecols=[1,2,5,12]) 
    if not np.array_equal(catalog['id'], prob_data[0]): 
        raise ValueError
    catalog['p_sat'] = prob_data[2]
    catalog['group_id'] = prob_data[1]
    catalog['angradius_halo'] = prob_data[3]
    # cuts 
    if catalog['ssfr_tot_mpajhu'].min() == -999.: 
        N_cat = len(catalog['p_sat']) 
        nan_cuts = np.where(catalog['ssfr_tot_mpajhu'] != -999.) 
        for key in catalog.keys(): 
            catalog[key] = catalog[key][nan_cuts]
        
        print 'removed = ', N_cat - len(catalog['p_sat']) 

    tinkauff_file = ''.join([UT.dir_dat(), 'tinkauff/',
        'VAGCdr72_MPAJHU.GroupCat.Mass', str(Mass_cut), '.hdf5']) 
    f = h5py.File(tinkauff_file, 'w')
    grp = f.create_group('data')
    
    for column in catalog.keys(): 
        grp.create_dataset(column, data=catalog[column])
    f.close()
    return None


def TinKauff_IsolationGroupCat(Mass_cut=9.25): 
    ''' Read in the Tinker-Kauffmann Group catalog generated from VAGC 
    dr72bright34 with MPA-JHU galaxy property values. Group catalog is 
    generated from Jeremy's new isolation criteria group catalog
    '''
    dr72_file = ''.join([UT.dir_dat(), '/tinkauff/',
        'VAGCdr72_MPAJHU.IsoGroupCat.Mass', str(Mass_cut), '.hdf5']) 
    catalog = {} 
    f = h5py.File(dr72_file, 'r')
    grp = f['data']
    for col in grp.keys(): 
        catalog[col] = grp[col].value  
    f.close() 
    return catalog


def Build_TinKauff_IsolationGroupCat(Mass_cut=9.25): 
    ''' Compile the outputs of Jeremy's new isolation criteria Group 
    Catalog algorithm in order to generate catalogs analogous to Kauffmann et al.(2013). 
    Hence TinKauff. 
    ''' 
    # galdata_corr file 
    galdata_file = ''.join([UT.dir_dat(), 'tinkauff/', 
        'clf_groups_JHU_M', str(Mass_cut), '_z0.017_fibcoll.galdata_corr'])
    gal_data = np.loadtxt(galdata_file, unpack=True, usecols=range(1,14))

    catalog = {
            'id': gal_data[0], 
            'ra': gal_data[5] * 57.2957795, 
            'dec': gal_data[6] * 57.2957795,
            'M_r': gal_data[1], 
            'M_g': gal_data[2],
            'z': gal_data[3]/299792.458,
            'mass_tot_mpajhu': np.log10(gal_data[4]), 
            'Dn4000': gal_data[7], 
            'ssfr_tot_mpajhu': gal_data[8], 
            'ssfr_fib_mpajhu': gal_data[9], 
            'sfr_tot_mpajhu': gal_data[10], 
            'sfr_fib_mpajhu': gal_data[11], 
            'mass_fib_mpajhu': gal_data[12], 
            } 

    # prob data from Jeremy's new isolation criteria group catalog
    prob_file = ''.join([UT.dir_dat(), 'tinkauff/', 
        'clf_groups_JHU_M', str(Mass_cut), '_z0.017_fibcoll.isolation.prob'])
    prob_data = np.loadtxt(prob_file, unpack=True, usecols=[0]) 
    catalog['p_sat'] = prob_data

    # cuts 
    if catalog['ssfr_tot_mpajhu'].min() == -999.: 
        N_cat = len(catalog['p_sat']) 
        nan_cuts = np.where(catalog['ssfr_tot_mpajhu'] != -999.) 
        for key in catalog.keys(): 
            catalog[key] = catalog[key][nan_cuts]
        
        print 'removed = ', N_cat - len(catalog['p_sat']) 

    tinkauff_file = ''.join([UT.dir_dat(), 'tinkauff/',
        'VAGCdr72_MPAJHU.IsoGroupCat.Mass', str(Mass_cut), '.hdf5']) 
    f = h5py.File(tinkauff_file, 'w')
    grp = f.create_group('data')
    
    for column in catalog.keys(): 
        grp.create_dataset(column, data=catalog[column])
    f.close()
    return None
