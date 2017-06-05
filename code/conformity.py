import numpy as np

# --- local --- 
import util as UT
import catalog as clog

from astropy.cosmology import FlatLambdaCDM

# --- plotting ---
import matplotlib.pyplot as plt
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors


def NeighborSSFR_rperp_PrimaryBins(catalog, rperp_bins=np.arange(0., 4.5, 0.5), 
        percentiles=[25, 50, 75, 90], quantiles=None, 
        primary_pipeline='mpajhu', primary_groupid='all', primary_massbin=[10., 10.5], 
        neighbor_pipeline='mpajhu', neighbor_groupid='all', neighbor_massbin=None): 
    ''' Calculate the mean and median SSFR(r_perp) of neighboring galaxies
    for primaries in SSFR percentile bins.
    '''
    # SSFR binning of primaries
    cut_primary = PrimaryIndices(catalog, 
            pipeline=primary_pipeline, 
            group_id=primary_groupid, 
            massbin=primary_massbin)
    
    # SSFR of primaries after final cut
    if primary_pipeline == 'vagc': 
        ssfr_cut_primary = catalog['ssfr'][cut_primary]
    elif primary_pipeline == 'mpajhu': 
        ssfr_cut_primary = catalog['ssfr_tot_mpajhu'][cut_primary]
    
    primary_SSFRbin_limits, primary_SSFRbin_list = SSFR_percentilebins(ssfr_cut_primary, 
            quantiles=quantiles, percentiles=percentiles)
    
    # for each of the primary SSFR bins 
    neighSSFR_rperp_primarybins = []
    for i_ssfrbin, primary_SSFRbin in enumerate(primary_SSFRbin_list): 
        # neighbor indices of primary galaxies within this SSFR bin  
        neigh_inbin = np.concatenate(    
                [np.array(catalog['neighbor_indices_'+neighbor_pipeline][i]) 
                    for i in cut_primary[primary_SSFRbin]]).astype('int') 
        # neighbor r_perp of primary galaxies within this SSFR bin  
        neigh_rperp = np.concatenate(    
                [np.array(catalog['neighbor_rperp_'+neighbor_pipeline][i]) 
                    for i in cut_primary[primary_SSFRbin]])

        if neighbor_pipeline == 'vagc': 
            neigh_ssfr = catalog['ssfr'][neigh_inbin]
            neigh_mass = catalog['mass'][neigh_inbin]
        elif neighbor_pipeline == 'mpajhu': 
            neigh_ssfr = catalog['ssfr_fib_mpajhu'][neigh_inbin]
            neigh_mass = catalog['mass_tot_mpajhu'][neigh_inbin] #neighbor *total* mass 

        if neighbor_groupid != 'all': # include all neighbors 
            neigh_psat = catalog['p_sat'][neigh_inbin]    # satellite probability of neighbors

        # *neighbor* SSFR(r_perp)
        neighborSSFR_rperp = [] 
        for i_rperp in range(len(rperp_bins)-1): 
            cut_nan = (np.isnan(neigh_ssfr) == False)
            cut_rperp = (neigh_rperp > rperp_bins[i_rperp]) & \
                    (neigh_rperp <= rperp_bins[i_rperp+1])
            if neighbor_groupid == 'all': # include all neighbors 
                cut_groupid = np.repeat(True, len(neigh_rperp))
            elif neighbor_groupid == 'centrals': # *central* neighbors 
                cut_groupid = (neigh_psat <= 0.5)
            elif neighbor_groupid == 'pure_centrals': # *pure central* neighbors 
                cut_groupid = (neigh_psat <= 0.01)
            else: 
                raise ValueError
            if neighbor_massbin is None: 
                cut_mass = np.repeat(True, len(neigh_rperp))
            else: 
                cut_mass = (neigh_mass >= neighbor_massbin[0]) & (neigh_mass < neighbor_massbin[1])
    
            cut_tot_neigh = np.where(cut_nan & cut_rperp & cut_groupid & cut_mass) # total cut
            neighborSSFR_rperp.append(neigh_ssfr[cut_tot_neigh])
        neighSSFR_rperp_primarybins.append(neighborSSFR_rperp)

    output_dict = {
            'primary_SSFRbin_limits': primary_SSFRbin_limits,
            'primary_SSFRbin_label': ['0 - 25\%', '25 - 50\%', '50 - 75\%', '> 75\%', '> 90\%'],
            'neighbor_SSFR_rperp_primary_bins': neighSSFR_rperp_primarybins
            }
    return output_dict 


def zSubsample_NeighborSSFR_rperp_PrimaryBins(catalog, lowhigh, 
        rperp_bins=np.arange(0., 4.5, 0.5), percentiles=[25, 50, 75, 90], quantiles=None, 
        primary_pipeline='mpajhu', primary_groupid='all', primary_massbin=[10., 10.5], 
        neighbor_pipeline='mpajhu', neighbor_groupid='all', neighbor_massbin=None): 
    #concat = clog.ConformCatalog()
    concat = clog.ConformCatalog('nothing', catalog_prop={})
    #jack_catalog = concat.Jackknife(catalog, n_jack, RADec_bins=RADec_bins) 
    jack_catalog, z_cut = concat.zSubsample(catalog, lowhigh) 
    print z_cut
    # SSFR binning of primaries
    cut_primary = PrimaryIndices(jack_catalog, 
            pipeline=primary_pipeline, 
            group_id=primary_groupid, 
            massbin=primary_massbin)
    
    # SSFR of primaries after final cut
    if primary_pipeline == 'vagc': 
        ssfr_cut_primary = jack_catalog['ssfr'][cut_primary]
    elif primary_pipeline == 'mpajhu': 
        ssfr_cut_primary = jack_catalog['ssfr_tot_mpajhu'][cut_primary]
    
    primary_SSFRbin_limits, primary_SSFRbin_list = SSFR_percentilebins(ssfr_cut_primary, 
            quantiles=quantiles, percentiles=percentiles)
    
    # for each of the primary SSFR bins 
    neighSSFR_rperp_primarybins = []
    for i_ssfrbin, primary_SSFRbin in enumerate(primary_SSFRbin_list): 
        # neighbor indices of primary galaxies within this SSFR bin  
        neigh_inbin = np.concatenate(    
                [np.array(jack_catalog['neighbor_indices_'+neighbor_pipeline][i]) 
                    for i in cut_primary[primary_SSFRbin]]).astype('int') 
        # neighbor r_perp of primary galaxies within this SSFR bin  
        neigh_rperp = np.concatenate(    
                [np.array(jack_catalog['neighbor_rperp_'+neighbor_pipeline][i]) 
                    for i in cut_primary[primary_SSFRbin]])

        if neighbor_pipeline == 'vagc': 
            neigh_ssfr = catalog['ssfr'][neigh_inbin]
            neigh_mass = catalog['mass'][neigh_inbin]
        elif neighbor_pipeline == 'mpajhu': 
            try: 
                neigh_ssfr = catalog['ssfr_fib_mpajhu'][neigh_inbin]
                neigh_mass = catalog['mass_tot_mpajhu'][neigh_inbin] #neighbor *total* mass 
            except KeyError: 
                neigh_ssfr = catalog['ssfr_fib'][neigh_inbin]
                neigh_mass = catalog['mass_tot'][neigh_inbin] #neighbor *total* mass 
        neigh_psat = catalog['p_sat'][neigh_inbin]    # satellite probability of neighbors
        neigh_z = catalog['z'][neigh_inbin]    # satellite probability of neighbors

        # *neighbor* SSFR(r_perp)
        neighborSSFR_rperp = [] 
        for i_rperp in range(len(rperp_bins)-1): 
            cut_nan = (np.isnan(neigh_ssfr) == False)
            if lowhigh == 'low': 
                cut_z_neigh = (neigh_z < z_cut) 
            elif lowhigh == 'high': 
                cut_z_neigh = (neigh_z >= z_cut) 
            else: 
                raise ValueError
            cut_rperp = (neigh_rperp > rperp_bins[i_rperp]) & \
                    (neigh_rperp <= rperp_bins[i_rperp+1])
            if neighbor_groupid == 'all': # include all neighbors 
                cut_groupid = np.repeat(True, len(neigh_rperp))
            elif neighbor_groupid == 'centrals': # *central* neighbors 
                cut_groupid = (neigh_psat <= 0.5)
            elif neighbor_groupid == 'pure_centrals': # *pure central* neighbors 
                cut_groupid = (neigh_psat <= 0.01)
            elif neighbor_groupid == 'pure_centrals_iso': # *pure central* neighbors for isolation criteria
                # this is Jeremy's new isolation based group catalog (6/2/17) 
                cut_groupid = (neigh_psat <= 0.1)
            else: 
                raise ValueError
            if neighbor_massbin is None: 
                cut_mass = np.repeat(True, len(neigh_rperp))
            else: 
                cut_mass = (neigh_mass >= neighbor_massbin[0]) & (neigh_mass < neighbor_massbin[1])
    
            cut_tot_neigh = np.where(cut_nan & cut_z_neigh & cut_rperp & cut_groupid & cut_mass) # total cut
            neighborSSFR_rperp.append(neigh_ssfr[cut_tot_neigh])
        neighSSFR_rperp_primarybins.append(neighborSSFR_rperp)

    output_dict = {
            'primary_SSFRbin_limits': primary_SSFRbin_limits,
            'primary_SSFRbin_label': ['0 - 25\%', '25 - 50\%', '50 - 75\%', '> 75\%', '> 90\%'],
            'neighbor_SSFR_rperp_primary_bins': neighSSFR_rperp_primarybins
            }
    return output_dict, z_cut


def Jackknife_NeighborSSFR_rperp_PrimaryBins(catalog, n_jack, RADec_bins=[5,5], 
        rperp_bins=np.arange(0., 4.5, 0.5), percentiles=[25, 50, 75, 90], quantiles=None, 
        primary_pipeline='mpajhu', primary_groupid='all', primary_massbin=[10., 10.5], 
        neighbor_pipeline='mpajhu', neighbor_groupid='all', neighbor_massbin=None): 
    concat = clog.ConformCatalog('nothing', catalog_prop={})
    jack_catalog = concat.Jackknife(catalog, n_jack, RADec_bins=RADec_bins) 
    # SSFR binning of primaries
    cut_primary = PrimaryIndices(jack_catalog, 
            pipeline=primary_pipeline, 
            group_id=primary_groupid, 
            massbin=primary_massbin)
    
    # SSFR of primaries after final cut
    if primary_pipeline == 'vagc': 
        ssfr_cut_primary = jack_catalog['ssfr'][cut_primary]
    elif primary_pipeline == 'mpajhu': 
        ssfr_cut_primary = jack_catalog['ssfr_tot_mpajhu'][cut_primary]
    
    primary_SSFRbin_limits, primary_SSFRbin_list = SSFR_percentilebins(ssfr_cut_primary, 
            quantiles=quantiles, percentiles=percentiles)
    
    # for each of the primary SSFR bins 
    neighSSFR_rperp_primarybins = []
    for i_ssfrbin, primary_SSFRbin in enumerate(primary_SSFRbin_list): 
        # neighbor indices of primary galaxies within this SSFR bin  
        neigh_inbin = np.concatenate(    
                [np.array(jack_catalog['neighbor_indices_'+neighbor_pipeline][i]) 
                    for i in cut_primary[primary_SSFRbin]]).astype('int') 
        # neighbor r_perp of primary galaxies within this SSFR bin  
        neigh_rperp = np.concatenate(    
                [np.array(jack_catalog['neighbor_rperp_'+neighbor_pipeline][i]) 
                    for i in cut_primary[primary_SSFRbin]])

        if neighbor_pipeline == 'vagc': 
            neigh_ssfr = catalog['ssfr'][neigh_inbin]
            neigh_mass = catalog['mass'][neigh_inbin]
        elif neighbor_pipeline == 'mpajhu': 
            neigh_ssfr = catalog['ssfr_fib_mpajhu'][neigh_inbin]
            neigh_mass = catalog['mass_tot_mpajhu'][neigh_inbin] #neighbor *total* mass 
        if neighbor_groupid != 'all': 
            neigh_psat = catalog['p_sat'][neigh_inbin]    # satellite probability of neighbors

        # *neighbor* SSFR(r_perp)
        neighborSSFR_rperp = [] 
        for i_rperp in range(len(rperp_bins)-1): 
            cut_nan = (np.isnan(neigh_ssfr) == False)
            cut_rperp = (neigh_rperp > rperp_bins[i_rperp]) & \
                    (neigh_rperp <= rperp_bins[i_rperp+1])
            if neighbor_groupid == 'all': # include all neighbors 
                cut_groupid = np.repeat(True, len(neigh_rperp))
            elif neighbor_groupid == 'centrals': # *central* neighbors 
                cut_groupid = (neigh_psat <= 0.5)
            elif neighbor_groupid == 'pure_centrals': # *pure central* neighbors 
                cut_groupid = (neigh_psat <= 0.01)
            else: 
                raise ValueError
            if neighbor_massbin is None: 
                cut_mass = np.repeat(True, len(neigh_rperp))
            else: 
                cut_mass = (neigh_mass >= neighbor_massbin[0]) & (neigh_mass < neighbor_massbin[1])
    
            cut_tot_neigh = np.where(cut_nan & cut_rperp & cut_groupid & cut_mass) # total cut
            neighborSSFR_rperp.append(neigh_ssfr[cut_tot_neigh])
        neighSSFR_rperp_primarybins.append(neighborSSFR_rperp)

    output_dict = {
            'primary_SSFRbin_limits': primary_SSFRbin_limits,
            'primary_SSFRbin_label': ['0 - 25\%', '25 - 50\%', '50 - 75\%', '> 75\%', '> 90\%'],
            'neighbor_SSFR_rperp_primary_bins': neighSSFR_rperp_primarybins
            }
    return output_dict 


def Plot_NeighborSSFR_rperp_PrimaryBins(cat_dict, rperp_bins=np.arange(0., 4.5, 0.5),
        primary_pipeline='mpajhu', primary_groupid='all', primary_massbin=[10., 10.5], 
        neighbor_pipeline='mpajhu', neighbor_groupid='all', neighbor_massbin=None):
    ''' Plot the SSFR of neighboring galaxies of primary galaxies binned in SSFRs. 
    This is calculated from the function NeighborSSFR_rperp_PrimaryBins.
    '''
    # read conformity catalog based on input catalog dictionary
    if cat_dict['name'] == 'tinker': 
        catalog_prop = {'Mrcut': cat_dict['Mrcut']} 
    elif cat_dict['name'] == 'tinkauff': 
        catalog_prop = {'Mass_cut': cat_dict['Mass_cut']} 
    elif cat_dict['name'] == 'kauff': 
        catalog_prop = {} 

    concat = clog.ConformCatalog(cat_dict['name'], catalog_prop=catalog_prop, 
            primary_delv=cat_dict['primary_delv'], primary_rperp=cat_dict['primary_rperp'],  
            neighbor_delv=cat_dict['neighbor_delv'], neighbor_rperp=cat_dict['neighbor_rperp'])
    #concat = clog.ConformCatalog(Mrcut=cat_dict['Mrcut'], 
    #        primary_delv=cat_dict['primary_delv'], primary_rperp=cat_dict['primary_rperp'],  
    #        neighbor_delv=cat_dict['neighbor_delv'], neighbor_rperp=cat_dict['neighbor_rperp'])
    catalog = concat.Read() 
    concat_file_spec = concat._FileSpec()    # conformity catalog specification string 

    results = NeighborSSFR_rperp_PrimaryBins(catalog, rperp_bins=rperp_bins,
            percentiles=[25, 50, 75, 90], quantiles=None, 
            primary_pipeline=primary_pipeline, 
            primary_groupid=primary_groupid, 
            primary_massbin=primary_massbin, 
            neighbor_pipeline=neighbor_pipeline, 
            neighbor_groupid=neighbor_groupid, 
            neighbor_massbin=neighbor_massbin)
    primary_SSFRbin_limits = results['primary_SSFRbin_limits']
    primary_SSFRbin_label = results['primary_SSFRbin_label']
    neighSSFR_rperp_primarybins = results['neighbor_SSFR_rperp_primary_bins']

    jack_results = []
    jack_bins = [5,5]
    for i_jack in range(1, jack_bins[0]*jack_bins[1]+1): 
        print i_jack
        jack_results_i = Jackknife_NeighborSSFR_rperp_PrimaryBins(
                catalog, n_jack=i_jack, RADec_bins=jack_bins, 
                rperp_bins=rperp_bins, percentiles=None, quantiles=primary_SSFRbin_limits, 
                primary_pipeline=primary_pipeline, 
                primary_groupid=primary_groupid, 
                primary_massbin=primary_massbin, 
                neighbor_pipeline=neighbor_pipeline,
                neighbor_groupid=neighbor_groupid, 
                neighbor_massbin=neighbor_massbin)
        #jack_results_i = Jackknife_NeighborSSFR_rperp_PrimaryBins(
        #        catalog, n_jack=i_jack, RADec_bins=jack_bins, 
        #        rperp_bins=rperp_bins, percentiles=[25,50,75,90], quantiles=None, 
        #        primary_pipeline=primary_pipeline, 
        #        primary_groupid=primary_groupid, 
        #        primary_massbin=primary_massbin, 
        #        neighbor_pipeline=neighbor_pipeline,
        #        neighbor_groupid=neighbor_groupid, 
        #        neighbor_massbin=neighbor_massbin)
        jack_results.append(jack_results_i['neighbor_SSFR_rperp_primary_bins'])

    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure()
    sub = fig.add_subplot(111)

    # Kauffmann et al.(2013) 
    for i_ssfr, ssfrs_rperp in enumerate(neighSSFR_rperp_primarybins):
        kauff_str = ['0to25', '25to50', '50to75', '75plus', '90plus']
        kauff_dat_file = ''.join([UT.dir_dat(), 'literature/', 'kauff2013_', kauff_str[i_ssfr], '.dat']) 
        kauff_rperp, kauff_ssfr = np.loadtxt(kauff_dat_file, unpack=True, usecols=[0,1]) 
        sub.scatter(kauff_rperp, kauff_ssfr, c=pretty_colors[i_ssfr], lw=0) 
        kauffplot, = sub.plot(kauff_rperp, kauff_ssfr, 
                c=pretty_colors[i_ssfr], lw=1, ls='--', label='Kauffmann+(2013)')

    ssfrplots = []
    for i_ssfr, ssfrs_rperp in enumerate(neighSSFR_rperp_primarybins):
        ssfr_tot = np.array([np.median(ssfrs) for ssfrs in ssfrs_rperp]) 
        ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
                c=pretty_colors[i_ssfr], lw=3, ls='-', label=primary_SSFRbin_label[i_ssfr]) 
    
        err_jack = np.zeros(len(ssfr_tot)) 
        for jack_result in jack_results: 
            ssfr_jack = np.array([np.median(ssfrs) for ssfrs in jack_result[i_ssfr]]) 
            err_jack += (ssfr_jack - ssfr_tot)**2
        err_jack *= np.float(jack_bins[0] * jack_bins[1] - 1)/np.float(jack_bins[0] * jack_bins[1])
        sub.errorbar(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, yerr=np.sqrt(err_jack), 
                c=pretty_colors[i_ssfr], elinewidth=2, capsize=5) 
        ssfrplots.append(ssfrplot) 

    if primary_pipeline == 'vagc': 
        label_pipe = 'vagc'
        label_ssfr_neigh = 'vagc'
    elif primary_pipeline == 'mpajhu': 
        label_pipe = 'tot; mpajhu'
        label_ssfr_neigh = 'fib; mpajhu'
    sub.set_title(r"$10.0 < \mathtt{log}\mathcal{M}^\mathtt{("+label_pipe+")}_* < 10.5$ ", fontsize=25) 
    # axes
    sub.set_xlim([0, 4]) 
    sub.set_xticks([0, 1, 2, 3, 4]) 
    sub.set_ylim([-12.25, -9.75])
    sub.set_yticks([-12., -11., -10.]) 
    sub.text(0.1, -11.9, r'Ranked in $\mathtt{log}\mathtt{SSFR^{('+label_pipe+')}}$', fontsize=20) 
    if neighbor_groupid == 'centrals': 
        sub.text(1.3, -10., 'Central Neighbors', fontsize=20)
    elif neighbor_groupid == 'pure_centrals': 
        sub.text(1.3, -10., 'Pure Central Neighbors', fontsize=20)
    if neighbor_massbin is not None:  
        sub.text(1.3, -10.2, 
                str(neighbor_massbin[0])+'$< $log$\mathcal{M}_*^\mathtt{neigh} <$'+str(neighbor_massbin[1]), 
                fontsize=20)
    sub.minorticks_on() 

    sub.set_xlabel(r'$\mathtt{R_{\perp}}$ [Mpc]', fontsize=25) 
    sub.set_ylabel(r'median log($\mathtt{SSFR_{(neigh)}^{('+label_ssfr_neigh+')}}$ [$\mathtt{yr}^{-1}$])', 
            fontsize=25) 
    first_legend = sub.legend(handles=[kauffplot], loc='lower right', handletextpad=0.1)
    ax = plt.gca().add_artist(first_legend)
    sub.legend(handles=ssfrplots, loc='upper right', ncol=2, handletextpad=0.1) 
    
    str_primary_groupid = 'PrimaryAll'       # primary category 
    if primary_groupid == 'centrals':   
        str_primary_groupid = '.PrimaryCentral'
    elif primary_groupid == 'pure_centrals': 
        str_primary_groupid= '.PrimaryPureCentral'
    elif primary_groupid == 'satellites': 
        str_primary_groupid = '.PrimarySatellite'
    if neighbor_groupid == 'all':      # neighbor catalogy
        str_neigh_groupid = '.NeighborAll'
    elif neighbor_groupid == 'centrals': 
        str_neigh_groupid = '.NeighborCentral'
    elif neighbor_groupid == 'pure_centrals': 
        str_neigh_groupid = '.NeighborPureCentral'
    if neighbor_massbin is None:  
        str_neigh_massbin = '' 
    else: 
        str_neigh_massbin = '.'+'_'.join([str(neighbor_massbin[0]), str(neighbor_massbin[1])]) 
    if cat_dict['name'] == 'tinkauff': 
        str_catalog = cat_dict['name']+'_'+str(cat_dict['Mass_cut'])
    else: 
        str_catalog = cat_dict['name']
    fig_file = ''.join([UT.dir_fig(), 
        'neighborSSFR_rprep_primarybins', '.', str_catalog, 
        concat_file_spec, str_primary_groupid, '.', primary_pipeline.upper(),
        str_neigh_groupid, str_neigh_massbin, '.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 


def Plot_Primary_Groups(cat_dict, primary_pipeline='mpajhu', primary_groupid='all', primary_massbin=[10., 10.5]):
    ''' Plot other members of groups that are within the primary 
    '''
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3) # cosmology 

    # read conformity catalog based on input catalog dictionary
    if cat_dict['name'] == 'tinker': 
        catalog_prop = {'Mrcut': cat_dict['Mrcut']} 
    elif cat_dict['name'] == 'tinkauff': 
        catalog_prop = {'Mass_cut': cat_dict['Mass_cut']} 
    elif cat_dict['name'] == 'kauff': 
        catalog_prop = {} 

    concat = clog.ConformCatalog(cat_dict['name'], catalog_prop=catalog_prop, 
            primary_delv=cat_dict['primary_delv'], primary_rperp=cat_dict['primary_rperp'],  
            neighbor_delv=cat_dict['neighbor_delv'], neighbor_rperp=cat_dict['neighbor_rperp'])
    catalog = concat.Read() 
    concat_file_spec = concat._FileSpec()    # conformity catalog specification string 

    cut_primary = PrimaryIndices(catalog, 
            pipeline=primary_pipeline, 
            group_id=primary_groupid, 
            massbin=primary_massbin)

    group_ids = catalog['group_id'].astype('int')  # galaxy id of primary catlaog
    # read group data 
    group_file = ''.join([UT.dir_dat(), '/tinkauff/',
        'clf_groups_JHU_M9.25_z0.017_fibcoll.groups']) 
    group_data = np.loadtxt(group_file, unpack=True, usecols=[1, 3]) 
    
    #rand_group_ids = np.random.choice(range(len(cut_primary)), len(cut_primary), replace=False) 
    #rand_group_ids = range(10)#range(len(cut_primary))
    rand_group_ids = [8, 12, 33, 14]
    
    rand_group_id_list = []
    widths = [] 
    for rand_group_id in rand_group_ids: 
        others = list(np.where(group_ids == group_ids[cut_primary[rand_group_id]])[0])

        if (np.sum(catalog['p_sat'][others] < 0.5) == 1) and (len(others) > 5) and (len(rand_group_id_list) < 4): 
            rand_group_id_list.append(rand_group_id) 

            widths.append(np.ceil(catalog['ra'][others].max() - catalog['ra'][others].min()))
            widths.append(np.ceil(catalog['dec'][others].max() - catalog['dec'][others].min()))
    
    wid_max = np.max(widths)
    print wid_max

    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(24,6))
    bkgd = fig.add_subplot(111, frameon=False)
    bkgd.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    for ii, i_group in enumerate(rand_group_id_list): 
        sub = fig.add_subplot(1, 4, ii+1, aspect='equal')
        
        halo_radius = catalog['angradius_halo'][cut_primary[i_group]]
        halo_radius *= 57.2958
        others = list(np.where(group_ids == group_ids[cut_primary[i_group]])[0])

        if np.sum(catalog['p_sat'][others] < 0.5) != 1: 
            continue 
        if len(others) < 5: 
            continue
        M_group = group_data[1][np.where(group_data[0] == group_ids[cut_primary[i_group]])[0]] 

        central_index = others[np.where(catalog['p_sat'][others] < 0.5)[0]]

        d_M = cosmo.comoving_distance(catalog['z'][cut_primary[i_group]])
        theta = 0.5/d_M.value * 57.2958
        
        # others in the group that's within the primary cut
        c_kms = 299792.458  # speed of light
        
        delv_cut = np.abs(catalog['z'] - catalog['z'][cut_primary[i_group]]) < 500./c_kms
        mass_cut = catalog['mass_tot_mpajhu'] < np.log10(0.5) + catalog['mass_tot_mpajhu'][cut_primary[i_group]]
        circle_cut = ((catalog['ra'] - catalog['ra'][cut_primary[i_group]])**2 + (catalog['dec'] - catalog['dec'][cut_primary[i_group]])**2) < theta**2
        other_inside = list(np.where((group_ids == group_ids[cut_primary[i_group]]) & delv_cut & circle_cut & mass_cut)[0])

        inside_circle_but_outside = np.where(
                (group_ids == group_ids[cut_primary[i_group]]) & circle_cut & mass_cut & 
                (np.abs(catalog['z'] - catalog['z'][cut_primary[i_group]]) >= 500./c_kms))

        if len(inside_circle_but_outside[0]) > 0: 
            print 'primary', catalog['z'][cut_primary[i_group]] # >= 500./c_kms))
            print catalog['z'][inside_circle_but_outside] 
            print c_kms * (catalog['z'][cut_primary[i_group]] - catalog['z'][inside_circle_but_outside]) 

        others.remove(cut_primary[i_group])
        others.remove(central_index)
        for inside in other_inside: 
            others.remove(inside)

        others_size = 50.*(catalog['mass_tot_mpajhu'][others]-9.0)
        if ii == 0: 
            sub.scatter(catalog['ra'][others], catalog['dec'][others], 
                    lw=0, c='k', alpha=0.7, s=others_size, label='Group Members')
        else: 
            sub.scatter(catalog['ra'][others], catalog['dec'][others], 
                    lw=0, c='k', alpha=0.7, s=others_size)
        # central 
        central_size = 50.*(catalog['mass_tot_mpajhu'][central_index]-9.0)
        if ii == 0: 
            sub.scatter(catalog['ra'][central_index], catalog['dec'][central_index], 
                    lw=0, c='y', s=others_size, label='Central')
        else: 
            sub.scatter(catalog['ra'][central_index], catalog['dec'][central_index], 
                    lw=0, c='y', s=others_size)
        halo_circ = plt.Circle((catalog['ra'][central_index], catalog['dec'][central_index]), 
            halo_radius, color=pretty_colors[0], fill=False, lw=2, linestyle='dashed') 
        sub.add_patch(halo_circ)

        primary_size = 50.*(catalog['mass_tot_mpajhu'][cut_primary[i_group]]-9.0)
        if ii == 1: 
            sub.scatter(catalog['ra'][cut_primary[i_group]], catalog['dec'][cut_primary[i_group]], 
                    lw=0, c=pretty_colors[3], s=primary_size, label='Primary')
        else: 
            sub.scatter(catalog['ra'][cut_primary[i_group]], catalog['dec'][cut_primary[i_group]], 
                    lw=0, c=pretty_colors[3], s=primary_size)
        circ = plt.Circle((catalog['ra'][cut_primary[i_group]], catalog['dec'][cut_primary[i_group]]), theta, 
                color=pretty_colors[3], fill=False, lw=2) 
        sub.add_patch(circ)
        
        if len(other_inside) > 0:
            other_inside_size = 50.*(catalog['mass_tot_mpajhu'][other_inside]-9.0)
            sub.scatter(catalog['ra'][other_inside], catalog['dec'][other_inside], 
                    lw=0, c='r', s=other_inside_size) 
        
        if ii == 1: 
            sub.scatter([0], [0], lw=0, c='r', 
                    label=r'$\mathtt{M_{gal} < \frac{1}{2} \; M_{primary}}$')

        # x-axis 
        #xaxis_mid = 0.5*(np.floor(catalog['ra'][others].min()) + np.ceil(catalog['ra'][others].max()))
        #xaxis_wid = np.ceil(catalog['ra'][others].max()) - np.floor(catalog['ra'][others].min()) 
        xaxis_mid = np.median(catalog['ra'][others]) 
        #xaxis_mid = catalog['ra'][central_index]
        #xaxis_wid = np.ceil(catalog['ra'][others].max() - catalog['ra'][others].min()) 
        #yaxis_mid = 0.5*(np.floor(catalog['dec'][others].min()) + np.ceil(catalog['dec'][others].max()))
        #yaxis_wid = np.ceil(catalog['dec'][others].max()) - np.floor(catalog['dec'][others].min()) 
        #yaxis_mid = catalog['dec'][central_index]
        yaxis_mid = np.median(catalog['dec'][others]) 
        #yaxis_wid = np.ceil(catalog['dec'][others].max() - catalog['dec'][others].min()) 
        #wid_max = np.max([xaxis_wid, yaxis_wid]) 
        
        sub.text(0.4, 0.08, '$\mathtt{log\; M_{group} = }$'+str(round(np.log10(M_group),1)), 
                fontsize=25, ha='center', va='center', transform=sub.transAxes) 
        if i_group in [7, 8, 10]: 
            sub.set_xlim([179.2, 183.2])
            sub.set_xticks([180, 181, 182, 183])
            sub.set_ylim([0, 4])
            sub.set_yticks([0, 1, 2, 3, 4])
            
            #sub.text(179.2, 0.2, '$\mathtt{log\; M_{group} = }$'+str(round(np.log10(M_group),1)), fontsize=25)
            #sub.text(182, 3.5, 'ID='+str(i_group), fontsize=20)
        elif i_group in [12]: 
            sub.set_xlim([40.4, 42.4])
            sub.set_xticks([41, 42])
            sub.set_ylim([-9.1, -7.1])
            sub.set_yticks([-9, -8])
            #sub.text(40.4, -9.1, '$\mathtt{log\; M_{group} = }$'+str(round(np.log10(M_group),1)), fontsize=25)
            #sub.text(41.75, -7.4, 'ID='+str(i_group), fontsize=20)
        elif i_group == 14: 
            sub.set_xlim([218.5, 221.5])
            sub.set_xticks([219, 220, 221])
            sub.set_ylim([2, 5])
            sub.set_yticks([2, 3, 4, 5])
            #sub.text(218.7, 2.2, '$\mathtt{log\; M_{group} = }$'+str(round(np.log10(M_group),1)), fontsize=25)
            #sub.text(220.9, 4.6, 'ID='+str(i_group), fontsize=20)
        elif i_group == 33: 
            sub.set_xlim([223., 225.4])
            sub.set_xticks([223, 224, 225])
            sub.set_ylim([8.2, 10.6])
            sub.set_yticks([9, 10])
            #sub.text(223.2, 8.2, '$\mathtt{log\; M_{group} = }$'+str(round(np.log10(M_group),1)), fontsize=25)
            #sub.text(225., 10., 'ID='+str(i_group), fontsize=20)
        elif i_group == 52: 
            sub.set_xlim([199.8, 202.2])
            sub.set_xticks([200, 201, 202])
            sub.set_ylim([12.8, 15.2])
            sub.set_yticks([13, 14, 15])
            #sub.text(200, 13., '$\mathtt{log\; M_{group} = }$'+str(round(np.log10(M_group),1)), fontsize=25)
            #sub.text(201.8, 14.8, 'ID='+str(i_group), fontsize=20)
        elif i_group == 57: 
            sub.set_xlim([243.4, 245.4])
            sub.set_xticks([244, 245])
            sub.set_ylim([34, 36])
            sub.set_yticks([34, 35, 36])
            #sub.text(243.6, 34.2, '$\mathtt{log\; M_{group} = }$'+str(round(np.log10(M_group),1)), fontsize=25)
            #sub.text(244.8, 35.6, 'ID='+str(i_group), fontsize=20)
        else: 
            sub.set_xlim([np.floor(xaxis_mid - 0.5*wid_max), np.ceil(xaxis_mid + 0.5*wid_max)])
            sub.set_xticks(range(int(np.floor(xaxis_mid - 0.5*wid_max)), int(np.ceil(xaxis_mid + 0.5*wid_max))+1))
            # y-axis 
            sub.set_ylim([np.floor(yaxis_mid - 0.5*wid_max), np.ceil(yaxis_mid + 0.5*wid_max)])
            sub.set_yticks(range(int(np.floor(yaxis_mid - 0.5*wid_max)), int(np.ceil(yaxis_mid + 0.5*wid_max))+1))

        if ii in [0,1]: 
            sub.legend(loc='upper left', markerscale=2, scatterpoints=1, handletextpad=0.1)
        sub.minorticks_on() 

    bkgd.set_xticklabels([]) 
    bkgd.set_xlabel('RA', labelpad=10, fontsize=30) 
    bkgd.set_yticklabels([]) 
    bkgd.set_ylabel('Dec', labelpad=30, fontsize=30) 
        
    str_primary_groupid = 'PrimaryAll'       # primary category 
    if primary_groupid == 'centrals':   
        str_primary_groupid = '.PrimaryCentral'
    elif primary_groupid == 'pure_centrals': 
        str_primary_groupid= '.PrimaryPureCentral'
    elif primary_groupid == 'satellites': 
        str_primary_groupid = '.PrimarySatellite'
    if cat_dict['name'] == 'tinkauff': 
        str_catalog = cat_dict['name']+'_'+str(cat_dict['Mass_cut'])
    else: 
        str_catalog = cat_dict['name']
    
    #fig_file = ''.join([UT.dir_fig(), 
    #    'Group_of_primary', '.', str_catalog, 
    #    concat_file_spec, str_primary_groupid, '.', primary_pipeline.upper(), 
    #    '.png']) 
    #fig.savefig(fig_file, bbox_inches='tight') 
    
    fig_file = ''.join([UT.dir_fig(), 
        'Group_of_primary', '.', str_catalog, 
        concat_file_spec, str_primary_groupid, '.', primary_pipeline.upper(), 
        '.pdf']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 


def Plot_NeighborSSFR_rperp_PrimaryBins_zSubsample(cat_dict, rperp_bins=np.arange(0., 4.5, 0.5),
        primary_pipeline='mpajhu', primary_groupid='all', primary_massbin=[10., 10.5], 
        neighbor_pipeline='mpajhu', neighbor_groupid='all', neighbor_massbin=None):
    ''' Plot the SSFR of neighboring galaxies of primary galaxies binned in SSFRs 
    This is calculated from the function NeighborSSFR_rperp_PrimaryBins.
    '''
    # read conformity catalog based on input catalog dictionary
    if cat_dict['name'] == 'tinker': 
        catalog_prop = {'Mrcut': cat_dict['Mrcut']} 
    elif cat_dict['name'] == 'tinkauff': 
        catalog_prop = {'Mass_cut': cat_dict['Mass_cut']} 

    concat = clog.ConformCatalog(cat_dict['name'], catalog_prop=catalog_prop, 
            primary_delv=cat_dict['primary_delv'], primary_rperp=cat_dict['primary_rperp'],  
            neighbor_delv=cat_dict['neighbor_delv'], neighbor_rperp=cat_dict['neighbor_rperp'])
    catalog = concat.Read() 
    concat_file_spec = concat._FileSpec()    # conformity catalog specification string 

    results = NeighborSSFR_rperp_PrimaryBins(catalog, rperp_bins=rperp_bins,
            percentiles=[25, 50, 75, 90], quantiles=None, 
            primary_pipeline=primary_pipeline, 
            primary_groupid=primary_groupid, 
            primary_massbin=primary_massbin, 
            neighbor_pipeline=neighbor_pipeline, 
            neighbor_groupid=neighbor_groupid, 
            neighbor_massbin=neighbor_massbin)
    primary_SSFRbin_limits = results['primary_SSFRbin_limits']
    primary_SSFRbin_label = results['primary_SSFRbin_label']
    neighSSFR_rperp_primarybins = results['neighbor_SSFR_rperp_primary_bins']
    
    lowz_results, low_z_cut = zSubsample_NeighborSSFR_rperp_PrimaryBins(catalog, 'low', rperp_bins=rperp_bins,
            percentiles=[25, 50, 75, 90], quantiles=primary_SSFRbin_limits, 
            primary_pipeline=primary_pipeline, 
            primary_groupid=primary_groupid, 
            primary_massbin=primary_massbin, 
            neighbor_pipeline=neighbor_pipeline, 
            neighbor_groupid=neighbor_groupid, 
            neighbor_massbin=neighbor_massbin)
    highz_results, high_z_cut = zSubsample_NeighborSSFR_rperp_PrimaryBins(catalog, 'high', rperp_bins=rperp_bins,
            percentiles=[25, 50, 75, 90], quantiles=primary_SSFRbin_limits, 
            primary_pipeline=primary_pipeline, 
            primary_groupid=primary_groupid, 
            primary_massbin=primary_massbin, 
            neighbor_pipeline=neighbor_pipeline, 
            neighbor_groupid=neighbor_groupid, 
            neighbor_massbin=neighbor_massbin)
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(15,7))

    for i_res, res in enumerate([lowz_results, highz_results]): 
        sub = fig.add_subplot(1,2,i_res+1)

        # Kauffmann et al.(2013) 
        for i_ssfr, ssfrs_rperp in enumerate(neighSSFR_rperp_primarybins):
            kauff_str = ['0to25', '25to50', '50to75', '75plus', '90plus']
            kauff_dat_file = ''.join([UT.dir_dat(), 'literature/', 'kauff2013_', kauff_str[i_ssfr], '.dat']) 
            kauff_rperp, kauff_ssfr = np.loadtxt(kauff_dat_file, unpack=True, usecols=[0,1]) 
            sub.scatter(kauff_rperp, kauff_ssfr, c=pretty_colors[i_ssfr], lw=0) 
            kauffplot, = sub.plot(kauff_rperp, kauff_ssfr, 
                    c=pretty_colors[i_ssfr], lw=1, ls='--', label='Kauffmann+(2013)')

        ssfrplots = []
        for i_ssfr, ssfrs_rperp in enumerate(res['neighbor_SSFR_rperp_primary_bins']):
            ssfr_tot = np.array([np.median(ssfrs) for ssfrs in ssfrs_rperp]) 
            ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
                    c=pretty_colors[i_ssfr], lw=3, ls='-', label=primary_SSFRbin_label[i_ssfr]) 
            ssfrplots.append(ssfrplot) 

        if primary_pipeline == 'vagc': 
            label_pipe = 'vagc'
            label_ssfr_neigh = 'vagc'
        elif primary_pipeline == 'mpajhu': 
            label_pipe = 'tot; mpajhu'
            label_ssfr_neigh = 'fib; mpajhu'
        # axes
        sub.set_xlim([0, 4]) 
        sub.set_xticks([0, 1, 2, 3, 4]) 
        sub.set_ylim([-12.25, -9.75])
        sub.set_yticks([-12., -11., -10.]) 
        sub.minorticks_on() 

        sub.set_xlabel(r'$\mathtt{R_{\perp}}$ [Mpc]', fontsize=25) 
        if i_res == 0: 
            sub.set_title(r"$\mathtt{z<"+str(round(low_z_cut,3))+"}$", fontsize=25) 
            if neighbor_groupid == 'centrals': 
                sub.text(1.3, -10., 'Central Neighbors', fontsize=20)
            elif neighbor_groupid == 'pure_centrals': 
                sub.text(1.3, -10., 'Pure Central Neighbors', fontsize=20)
            if neighbor_massbin is not None:  
                sub.text(1.3, -10.2, 
                        str(neighbor_massbin[0])+'$< $log$\mathcal{M}_*^\mathtt{neigh} <$'+str(neighbor_massbin[1]), 
                        fontsize=20)
            sub.set_ylabel(r'median log($\mathtt{SSFR_{(neigh)}^{('+label_ssfr_neigh+')}}$ [$\mathtt{yr}^{-1}$])', 
                    fontsize=25) 
            sub.legend(handles=[kauffplot], loc='lower right', handletextpad=0.1)
        else: 
            sub.set_title(r"$\mathtt{z\ge"+str(round(high_z_cut,3))+"}$", fontsize=25) 
            sub.text(0.1, -11.9, r'Ranked in $\mathtt{log}\mathtt{SSFR^{('+label_pipe+')}}$', fontsize=20) 
            sub.legend(handles=ssfrplots, loc='upper right', ncol=2, handletextpad=0.1) 
            sub.set_yticklabels([]) 
    
    str_primary_groupid = 'PrimaryAll'       # primary category 
    if primary_groupid == 'centrals':   
        str_primary_groupid = '.PrimaryCentral'
    elif primary_groupid == 'pure_centrals': 
        str_primary_groupid= '.PrimaryPureCentral'
    elif primary_groupid == 'satellites': 
        str_primary_groupid = '.PrimarySatellite'
    if neighbor_groupid == 'all':      # neighbor catalogy
        str_neigh_groupid = '.NeighborAll'
    elif neighbor_groupid == 'centrals': 
        str_neigh_groupid = '.NeighborCentral'
    elif neighbor_groupid == 'pure_centrals': 
        str_neigh_groupid = '.NeighborPureCentral'
    if neighbor_massbin is None:  
        str_neigh_massbin = '' 
    else: 
        str_neigh_massbin = '.'+'_'.join([str(neighbor_massbin[0]), str(neighbor_massbin[1])]) 
    fig_file = ''.join([UT.dir_fig(), 
        'neighborSSFR_rprep_primarybins_zSubsamples', 
        concat_file_spec, str_primary_groupid, '.', primary_pipeline.upper(),
        str_neigh_groupid, str_neigh_massbin, '.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 


def Plot_PrimaryBins_Pmass(cat_dict,
        primary_pipeline='mpajhu', primary_groupid='all', primary_massbin=[10., 10.5]):
    ''' Plot the mass distribution of primary galaxies within bins of SSFR. 
    '''
    # read conformity catalog based on input catalog dictionary
    if cat_dict['name'] == 'tinker': 
        catalog_prop = {'Mrcut': cat_dict['Mrcut']} 
    elif cat_dict['name'] == 'tinkauff': 
        catalog_prop = {'Mass_cut': cat_dict['Mass_cut']} 
    elif cat_dict['name'] == 'kauff': 
        catalog_prop = {} 

    concat = clog.ConformCatalog(cat_dict['name'], catalog_prop=catalog_prop, 
            primary_delv=cat_dict['primary_delv'], primary_rperp=cat_dict['primary_rperp'],  
            neighbor_delv=cat_dict['neighbor_delv'], neighbor_rperp=cat_dict['neighbor_rperp'])
    #concat = clog.ConformCatalog(Mrcut=cat_dict['Mrcut'], 
    #        primary_delv=cat_dict['primary_delv'], primary_rperp=cat_dict['primary_rperp'],  
    #        neighbor_delv=cat_dict['neighbor_delv'], neighbor_rperp=cat_dict['neighbor_rperp'])
    catalog = concat.Read() 
    concat_file_spec = concat._FileSpec()    # conformity catalog specification string 

    # SSFR binning of primaries
    cut_primary = PrimaryIndices(catalog, 
            pipeline=primary_pipeline, 
            group_id=primary_groupid, 
            massbin=primary_massbin)
    
    # SSFR of primaries after final cut
    if primary_pipeline == 'vagc': 
        ssfr_cut_primary = catalog['ssfr'][cut_primary]
    elif primary_pipeline == 'mpajhu': 
        ssfr_cut_primary = catalog['ssfr_tot_mpajhu'][cut_primary]
    
    primary_SSFRbin_limits, primary_SSFRbin_list = SSFR_percentilebins(ssfr_cut_primary, 
            quantiles=None, percentiles=[25, 50, 75, 90])
    primary_SSFRbin_labels = ['0 - 25\%', '25 - 50\%', '50 - 75\%', '> 75\%', '> 90\%']

    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure()
    sub = fig.add_subplot(111)
    
    mass_bins = np.arange(10.0, 10.6, 0.1)
    for i_ssfrbin, primary_SSFRbin in enumerate(primary_SSFRbin_list): 
        if primary_pipeline == 'vagc': 
            primary_mass = catalog['mass'][cut_primary[primary_SSFRbin]]
        elif primary_pipeline == 'mpajhu': 
            primary_mass = catalog['mass_tot_mpajhu'][cut_primary[primary_SSFRbin]]
        #print primary_SSFRbin_labels[i_ssfrbin]
        #print np.mean(primary_mass)
        pdf, bins = np.histogram(primary_mass, 
                range=[10., 10.5], bins=10) 
        sub.plot(0.5 * (bins[:-1] + bins[1:]), pdf, lw=3,
                color=pretty_colors[i_ssfrbin], label=primary_SSFRbin_labels[i_ssfrbin])

    # axes
    sub.set_xlim([10., 10.5]) 
    sub.set_xticks([10, 10.25, 10.5]) 
    #sub.set_ylim([-12.25, -9.75])
    #sub.set_yticks([-12., -11., -10.]) 
    sub.minorticks_on() 
    sub.legend(loc='upper right', ncol=2, handletextpad=0.1) 

    sub.set_xlabel(r'log$\mathtt{M_{*}}$', fontsize=25) 
    sub.set_ylabel(r'P(log$\mathtt{M_{*}}$)', fontsize=25) 
    
    str_primary_groupid = 'PrimaryAll'       # primary category 
    if primary_groupid == 'centrals':   
        str_primary_groupid = '.PrimaryCentral'
    elif primary_groupid == 'pure_centrals': 
        str_primary_groupid= '.PrimaryPureCentral'
    elif primary_groupid == 'satellites': 
        str_primary_groupid = '.PrimarySatellite'
    if cat_dict['name'] == 'tinkauff': 
        str_catalog = cat_dict['name']+'_'+str(cat_dict['Mass_cut'])
    else: 
        str_catalog = cat_dict['name']
    fig_file = ''.join([UT.dir_fig(), 
        'Pmass_primarybins', '.', str_catalog, 
        concat_file_spec, str_primary_groupid, '.', primary_pipeline.upper(), '.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 


def Plot_PrimaryBins_massbin_SSFR_rperp(cat_dict, rperp_bins=np.arange(0., 4.5, 0.5), 
        primary_pipeline='mpajhu', primary_groupid='all', primary_massbin=[10., 10.5],
        neighbor_pipeline='mpajhu', neighbor_groupid='all', neighbor_massbin=None): 
    ''' This is complicated. First calculate the mass distribuiton of the primaries in 
    percentile bins of their SSFRs. Then, construct a random sample of primaries that 
    reproduces the mass distribution of the 0-25% bin. Then calculate the 
    SSFR_neighbor(r_perp) of that sample to the rest. 
    '''
    # read conformity catalog based on input catalog dictionary
    if cat_dict['name'] == 'tinker': 
        catalog_prop = {'Mrcut': cat_dict['Mrcut']} 
    elif cat_dict['name'] == 'tinkauff': 
        catalog_prop = {'Mass_cut': cat_dict['Mass_cut']} 
    elif cat_dict['name'] == 'kauff': 
        catalog_prop = {} 

    concat = clog.ConformCatalog(cat_dict['name'], catalog_prop=catalog_prop, 
            primary_delv=cat_dict['primary_delv'], primary_rperp=cat_dict['primary_rperp'],  
            neighbor_delv=cat_dict['neighbor_delv'], neighbor_rperp=cat_dict['neighbor_rperp'])
    catalog = concat.Read() 
    concat_file_spec = concat._FileSpec()    # conformity catalog specification string 

    # SSFR binning of primaries
    cut_primary = PrimaryIndices(catalog, 
            pipeline=primary_pipeline, 
            group_id=primary_groupid, 
            massbin=primary_massbin)
    
    # SSFR of primaries after final cut
    if primary_pipeline == 'vagc': 
        ssfr_cut_primary = catalog['ssfr'][cut_primary]
    elif primary_pipeline == 'mpajhu': 
        ssfr_cut_primary = catalog['ssfr_tot_mpajhu'][cut_primary]
    
    primary_SSFRbin_limits, primary_SSFRbin_list = SSFR_percentilebins(ssfr_cut_primary, 
            quantiles=None, percentiles=[25, 50, 75, 90])
    primary_SSFRbin_labels = ['0 - 25\%', '25 - 50\%', '50 - 75\%', '> 75\%', '> 90\%']

    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure()
    sub = fig.add_subplot(111)
    
    # mass distribution of the lowest SSFR bin  
    if primary_pipeline == 'vagc': 
        primary_mass = catalog['mass'][cut_primary]
    elif primary_pipeline == 'mpajhu': 
        primary_mass = catalog['mass_tot_mpajhu'][cut_primary]

    mass_dist, bins = np.histogram(primary_mass[primary_SSFRbin_list[0]], range=[10., 10.5], bins=10) 
    random_indices = [] 
    for i_bin in range(len(bins)-1): 
        within_Mbin = np.where((primary_mass > bins[i_bin]) & (primary_mass <= bins[i_bin+1]))[0]
        random_indices.append(np.random.choice(within_Mbin, size=mass_dist[i_bin], replace=False))
    # indices of random sampel that reproduces P(M*) of lowest SSFR bin
    lowSSFR_random = np.concatenate(random_indices)
    print np.mean(primary_mass[lowSSFR_random])
    
    lowSSFR_random_converse = []
    for i in range(len(cut_primary)): 
        if i not in lowSSFR_random: 
            lowSSFR_random_converse.append(i)
    lowSSFR_random_converse = np.array(lowSSFR_random_converse)
    print np.mean(primary_mass[lowSSFR_random_converse])

    mass_dist, bins = np.histogram(primary_mass[primary_SSFRbin_list[-2]], range=[10., 10.5], bins=10) 
    random_indices = [] 
    for i_bin in range(len(bins)-1): 
        within_Mbin = np.where((primary_mass > bins[i_bin]) & (primary_mass <= bins[i_bin+1]))[0]
        random_indices.append(np.random.choice(within_Mbin, size=mass_dist[i_bin], replace=False))
    highSSFR_random = np.concatenate(random_indices)
    print np.mean(primary_mass[highSSFR_random])

    neighSSFR_rperp_primarybins = [] 
    for i_ssfrbin, primary_SSFRbin in enumerate([lowSSFR_random, lowSSFR_random_converse, highSSFR_random]): 
        # neighbor indices of primary galaxies within this SSFR bin  
        neigh_inbin = np.concatenate(    
                [np.array(catalog['neighbor_indices_'+neighbor_pipeline][i]) 
                    for i in cut_primary[primary_SSFRbin]]).astype('int') 
        # neighbor r_perp of primary galaxies within this SSFR bin  
        neigh_rperp = np.concatenate(    
                [np.array(catalog['neighbor_rperp_'+neighbor_pipeline][i]) 
                    for i in cut_primary[primary_SSFRbin]])

        if neighbor_pipeline == 'vagc': 
            neigh_ssfr = catalog['ssfr'][neigh_inbin]
            neigh_mass = catalog['mass'][neigh_inbin]
        elif neighbor_pipeline == 'mpajhu': 
            neigh_ssfr = catalog['ssfr_fib_mpajhu'][neigh_inbin]
            neigh_mass = catalog['mass_tot_mpajhu'][neigh_inbin] #neighbor *total* mass 

        if neighbor_groupid != 'all': # include all neighbors 
            neigh_psat = catalog['p_sat'][neigh_inbin]    # satellite probability of neighbors

        # *neighbor* SSFR(r_perp)
        neighborSSFR_rperp = [] 
        for i_rperp in range(len(rperp_bins)-1): 
            cut_nan = (np.isnan(neigh_ssfr) == False)
            cut_rperp = (neigh_rperp > rperp_bins[i_rperp]) & \
                    (neigh_rperp <= rperp_bins[i_rperp+1])
            if neighbor_groupid == 'all': # include all neighbors 
                cut_groupid = np.repeat(True, len(neigh_rperp))
            elif neighbor_groupid == 'centrals': # *central* neighbors 
                cut_groupid = (neigh_psat <= 0.5)
            elif neighbor_groupid == 'pure_centrals': # *pure central* neighbors 
                cut_groupid = (neigh_psat <= 0.01)
            else: 
                raise ValueError
            if neighbor_massbin is None: 
                cut_mass = np.repeat(True, len(neigh_rperp))
            else: 
                cut_mass = (neigh_mass >= neighbor_massbin[0]) & (neigh_mass < neighbor_massbin[1])
    
            cut_tot_neigh = np.where(cut_nan & cut_rperp & cut_groupid & cut_mass) # total cut
            neighborSSFR_rperp.append(neigh_ssfr[cut_tot_neigh])
        neighSSFR_rperp_primarybins.append(neighborSSFR_rperp)

    primary_SSFRbin_label = ['Random sample of \n$0-25\%$ mass distribution', 'others',
            'Random sample of \n$>75\%$ mass distribution' ] 
    ssfrplots = []
    for i_ssfr, ssfrs_rperp in enumerate(neighSSFR_rperp_primarybins):
        ssfr_tot = np.array([np.median(ssfrs) for ssfrs in ssfrs_rperp]) 
        ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
                c=pretty_colors[i_ssfr+1], lw=3, ls='-', label=primary_SSFRbin_label[i_ssfr]) 
        ssfrplots.append(ssfrplot) 

    # axes
    sub.set_xlim([0, 4]) 
    sub.set_xticks([0, 1, 2, 3, 4]) 
    sub.set_ylim([-12.25, -9.75])
    sub.set_yticks([-12., -11., -10.]) 
    if neighbor_groupid == 'centrals': 
        sub.text(1.3, -10., 'Central Neighbors', fontsize=20)
    elif neighbor_groupid == 'pure_centrals': 
        sub.text(1.3, -10., 'Pure Central Neighbors', fontsize=20)
    if neighbor_massbin is not None:  
        sub.text(1.3, -10.2, 
                str(neighbor_massbin[0])+'$< $log$\mathcal{M}_*^\mathtt{neigh} <$'+str(neighbor_massbin[1]), 
                fontsize=20)
    sub.minorticks_on() 
    
    if primary_pipeline == 'vagc': 
        label_pipe = 'vagc'
        label_ssfr_neigh = 'vagc'
    elif primary_pipeline == 'mpajhu': 
        label_pipe = 'tot; mpajhu'
        label_ssfr_neigh = 'fib; mpajhu'
    sub.set_xlabel(r'$\mathtt{R_{\perp}}$ [Mpc]', fontsize=25) 
    sub.set_ylabel(r'median log($\mathtt{SSFR_{(neigh)}^{('+label_ssfr_neigh+')}}$ [$\mathtt{yr}^{-1}$])', 
            fontsize=25) 
    sub.legend(handles=ssfrplots, loc='lower right', handletextpad=0.1) 

    str_primary_groupid = 'PrimaryAll'       # primary category 
    if primary_groupid == 'centrals':   
        str_primary_groupid = '.PrimaryCentral'
    elif primary_groupid == 'pure_centrals': 
        str_primary_groupid= '.PrimaryPureCentral'
    elif primary_groupid == 'satellites': 
        str_primary_groupid = '.PrimarySatellite'
    if neighbor_groupid == 'all':      # neighbor catalogy
        str_neigh_groupid = '.NeighborAll'
    elif neighbor_groupid == 'centrals': 
        str_neigh_groupid = '.NeighborCentral'
    elif neighbor_groupid == 'pure_centrals': 
        str_neigh_groupid = '.NeighborPureCentral'
    if neighbor_massbin is None:  
        str_neigh_massbin = '' 
    else: 
        str_neigh_massbin = '.'+'_'.join([str(neighbor_massbin[0]), str(neighbor_massbin[1])]) 
    if cat_dict['name'] == 'tinkauff': 
        str_catalog = cat_dict['name']+'_'+str(cat_dict['Mass_cut'])
    else: 
        str_catalog = cat_dict['name']
    fig_file = ''.join([UT.dir_fig(), 
        'test_convoluted_primarybins_massbin', '.', str_catalog, 
        concat_file_spec, str_primary_groupid, '.', primary_pipeline.upper(),
        str_neigh_groupid, str_neigh_massbin, '.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 


def PlotConformity_Primary_PDF_Rperp_bin(gal_prop, cat_dict, 
        primary_id='mpajhu', primary_massbin=[10., 10.5], cen_sat=False, 
        neighbors='all', neighbor_massbin=None): 
    ''' Divide the primary galaxies into bins of SSFR and then 
    examine the SSFR of their neighboring galaxies. More specifically
    we look at the SSFR distirubiton within a r_perp bin. 
    '''
    # read in conformity catalog 
    concat = clog.ConformCatalog(Mrcut=cat_dict['Mrcut'], 
            primary_delv=cat_dict['primary_delv'], 
            primary_rperp=cat_dict['primary_rperp'],  
            neighbor_delv=cat_dict['neighbor_delv'], 
            neighbor_rperp=cat_dict['neighbor_rperp'])
    catalog = concat.Read() 

    # first sort out all primary galaxies 
    is_primary = np.where(catalog['primary_'+primary_id] == 1)[0]
    if primary_id == 'vagc': 
        ssfr_primary = catalog['ssfr'][is_primary]
        mass_primary = catalog['mass'][is_primary]
        mass_min, mass_max = 9., 12. 
        ssfr_min, ssfr_max = -13.25, -9.25
    elif primary_id == 'mpajhu': 
        ssfr_primary = catalog['ssfr_tot_mpajhu'][is_primary]
        mass_primary = catalog['mass_tot_mpajhu'][is_primary]
        mass_min, mass_max = 9., 12. 
        ssfr_min, ssfr_max = -13., -8.75

    if cen_sat == False:   # all primaries within mass bin with SSFR
        primary_cut = np.where(
                (mass_primary > primary_massbin[0]) & (mass_primary < primary_massbin[1]) & 
                (np.isnan(ssfr_primary) == False))[0]
    elif cen_sat == 'centrals':  # primaries within mass bin with SSFR that are *centrals*
        primary_cut = np.where(
                (mass_primary > primary_massbin[0]) & (mass_primary < primary_massbin[1]) & 
                (np.isnan(ssfr_primary) == False) & (catalog['p_sat'][is_primary] <= 0.5)
                )[0]
    elif cen_sat == 'pure_centrals':  # primaries within mass bin with SSFR that are pure *centrals*
        primary_cut = np.where(
                (mass_primary > primary_massbin[0]) & (mass_primary < primary_massbin[1]) & 
                (np.isnan(ssfr_primary) == False) & (catalog['p_sat'][is_primary] <= 0.01)
                )[0]
    elif cen_sat == 'satellites': # primaries within mass bin with SSFR that are *satellites*
        primary_cut = np.where(
                (mass_primary > primary_massbin[0]) & (mass_primary < primary_massbin[1]) & 
                (np.isnan(ssfr_primary) == False) & (catalog['p_sat'][is_primary] > 0.5)
                )[0]
    else: 
        raise ValueError
    
    ssfr_primary_cut = ssfr_primary[primary_cut]
    print len(primary_cut), 'primaries within cut'
    
    q25, q50, q75, q90 = np.percentile(ssfr_primary_cut, [25, 50, 75, 90])
    #print q25, q50, q75, q90

    # centrals with 0-25th percentile SSFR
    bin_ssfr_0to25 = np.where(ssfr_primary_cut < q25)[0]
    bin_ssfr_25to50 = np.where((ssfr_primary_cut >= q25) & (ssfr_primary_cut < q50))[0]
    bin_ssfr_50to75 = np.where((ssfr_primary_cut >= q50) & (ssfr_primary_cut < q75))[0]
    bin_ssfr_75plus = np.where(ssfr_primary_cut >= q75)[0]
    bin_ssfr_90plus = np.where(ssfr_primary_cut >= q90)[0]
    bin_ssfr_list = [bin_ssfr_0to25, bin_ssfr_25to50, 
            bin_ssfr_50to75, bin_ssfr_75plus, bin_ssfr_90plus]
    bin_ssfr_label = ['0 - 25\%', '25 - 50\%', '50 - 75\%', '> 75\%', '> 90\%']
    PDF_primary_bin = [] 
    
    rperp_bins = np.arange(0., 5., 1.)
    
    mean_primary_bin, median_primary_bin = [], [] 
    for bin_ssfr_cut in bin_ssfr_list: 
        # get SSFR and r_perp of the neighbors of the primaries
        # in each of the SSFR bins. 
        neigh_mass, neigh_ssfr, neigh_rperp = [], [], [] 
        for ip in primary_cut[bin_ssfr_cut]: 
            cen_index = is_primary[ip]
            if cen_sat in ['central', 'pure_central']: 
                if not catalog['p_sat'][cen_index] <= 0.5:  # has to be central 
                    raise ValueError

            neigh_indices = np.array(catalog['neighbor_indices_'+primary_id][ip]).astype('int') 
            if neighbors == 'all':   # keep all neighbors 
                not_sat = np.repeat(True, len(neigh_indices)) 
            elif neighbors == 'centrals':   # remove all satellites 
                neigh_psat = catalog['p_sat'][neigh_indices]
                not_sat = neigh_psat <= 0.5
            elif neighbors == 'pure_centrals':   # remove all satellites 
                neigh_psat = catalog['p_sat'][neigh_indices]
                not_sat = neigh_psat <= 0.01
            else: 
                raise ValueError

            if neighbor_massbin is None: 
                mass_cut = np.repeat(True, len(neigh_indices)) 
            else: 
                if primary_id == 'vagc': 
                    neigh_M = catalog['mass'][neigh_indices] 
                elif primary_id == 'mpajhu': 
                    neigh_M = catalog['mass_tot_mpajhu'][neigh_indices]
                mass_cut = (neigh_M > neighbor_massbin[0]) & (neigh_M < neighbor_massbin[1])
            
            neigh_cut = np.where(not_sat & mass_cut) 
            neigh_rperp.append(np.array(catalog['neighbor_rperp_'+primary_id][ip])[neigh_cut])

            if primary_id == 'vagc': 
                neigh_ssfr.append(catalog['ssfr'][neigh_indices[neigh_cut]]) 
                neigh_mass.append(catalog['mass'][neigh_indices[neigh_cut]]) 
            elif primary_id == 'mpajhu': 
                try: 
                    neigh_ssfr.append(catalog['ssfr_fib_mpajhu'][neigh_indices[neigh_cut]])
                    neigh_mass.append(catalog['mass_tot_mpajhu'][neigh_indices[neigh_cut]]) #neighbor *total* mass 
                except KeyError: 
                    neigh_ssfr.append(catalog['ssfr_fib'][neigh_indices[neigh_cut]])
                    neigh_mass.append(catalog['mass_tot'][neigh_indices[neigh_cut]]) #neighbor *total* mass 

        neigh_rperp = np.concatenate(neigh_rperp)
        neigh_ssfr = np.concatenate(neigh_ssfr) 
        neigh_mass = np.concatenate(neigh_mass) 

        # calculate P(M*) in bins of r_perp 
        PDF_in_rperpbins = []  
        mean_in_rperpbins, median_in_rperpbins = [], [] 
        for i_rperp in range(len(rperp_bins)-1): 
            if gal_prop == 'mass': 
                bin_prop_rperp = np.where( 
                        (neigh_rperp > rperp_bins[i_rperp]) &  
                        (neigh_rperp <= rperp_bins[i_rperp+1]) & 
                        (np.isnan(neigh_mass) == False))
                pdf, bins = np.histogram(neigh_mass[bin_prop_rperp], 
                        range=[mass_min, mass_max], bins=20, normed=True) 
                mean_in_rperpbins.append(np.mean(neigh_mass[bin_prop_rperp]))
                median_in_rperpbins.append(np.median(neigh_mass[bin_prop_rperp]))
            elif gal_prop == 'ssfr': 
                bin_prop_rperp = np.where( 
                        (neigh_rperp > rperp_bins[i_rperp]) &  
                        (neigh_rperp <= rperp_bins[i_rperp+1]) & 
                        (np.isnan(neigh_ssfr) == False))
                pdf, bins = np.histogram(neigh_ssfr[bin_prop_rperp], 
                        range=[ssfr_min, ssfr_max], bins=20, normed=True) 
                mean_in_rperpbins.append(np.mean(neigh_ssfr[bin_prop_rperp]))
                median_in_rperpbins.append(np.median(neigh_ssfr[bin_prop_rperp]))
            PDF_in_rperpbins.append(pdf) 
        PDF_primary_bin.append(PDF_in_rperpbins)
        mean_primary_bin.append(mean_in_rperpbins)
        median_primary_bin.append(median_in_rperpbins)
    
    # plotting 
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(15,4))
    bkgd = fig.add_subplot(111, frameon=False) 
    for i_rp in range(len(PDF_primary_bin[0])):   # rperp bin 
        sub = fig.add_subplot(1,len(PDF_primary_bin[0]),i_rp+1)
        for i_pb in range(len(PDF_primary_bin)):              # primary bin 
            sub.plot(0.5*(bins[:-1]+bins[1:]), (PDF_primary_bin[i_pb])[i_rp], 
                    c=pretty_colors[i_pb], lw=2)
            sub.vlines(mean_primary_bin[i_pb][i_rp], 1.05, 1.2, 
                    color=pretty_colors[i_pb], linewidth=1.5)
            sub.vlines(median_primary_bin[i_pb][i_rp], 0., 0.15, 
                    color=pretty_colors[i_pb], linewidth=1.5)
            if i_rp == 0: 
                if gal_prop == 'mass': 
                    sub.text(mass_min+0.75, 1.025, 'Mean', fontsize=20)
                    sub.text(mass_min+0.55, 0.05, 'Median', fontsize=20) 
                elif gal_prop == 'ssfr': 
                    sub.text(ssfr_min+0.75, 0.725, 'Mean', fontsize=20)
                    sub.text(ssfr_min+0.55, 0.05, 'Median', fontsize=20) 
        sub.set_title(''.join(['$', str(rperp_bins[i_rp]), '< \mathtt{r}_\perp <', str(rperp_bins[i_rp+1]), '$']), 
                fontsize=20) 
        # axes
        if gal_prop == 'mass': 
            sub.set_xlim([mass_min, mass_max]) 
        elif gal_prop == 'ssfr': 
            sub.set_xlim([ssfr_min, ssfr_max]) 
        if i_rp < len(PDF_primary_bin[0])-1: 
            if gal_prop == 'mass': 
                sub.set_xticks([9., 10., 11.]) 
            elif gal_prop == 'ssfr': 
                sub.set_xticks([-13, -12, -11, -10]) 
        else: 
            if gal_prop == 'mass': 
                sub.set_xticks([9., 10., 11., 12.]) 
            elif gal_prop == 'ssfr': 
                sub.set_xticks([-13, -12, -11, -10, -9]) 
        if gal_prop == 'mass': 
            sub.set_ylim([0., 1.2]) 
            sub.set_yticks([0., 0.4, 0.8, 1.2]) 
        elif gal_prop == 'ssfr': 
            sub.set_ylim([0., 0.8]) 
            sub.set_yticks([0., 0.2, 0.4, 0.6, 0.8]) 
        sub.minorticks_on() 
        if i_rp != 0: 
            sub.set_yticklabels([]) 

    bkgd.set_xticklabels([]) 
    bkgd.set_yticklabels([]) 
    if gal_prop == 'mass': 
        if primary_id == 'vagc': 
            bkgd.set_xlabel(r'log($\mathtt{M_*^{vagc}}$ [$\mathtt{M_\odot}$])', fontsize=25, labelpad=25) 
            bkgd.set_ylabel(r'P(log$\mathtt{M_*^{vagc}}$)', fontsize=25, labelpad=25) 
        elif primary_id == 'mpajhu': 
            bkgd.set_xlabel(r'log($\mathtt{M_*^{tot; mpajhu}}$ [$\mathtt{M_\odot}$])', fontsize=25, labelpad=25) 
            bkgd.set_ylabel(r'P(log$\mathtt{M_*^{tot; mpajhu}}$)', fontsize=25, labelpad=25) 
    elif gal_prop == 'ssfr': 
        if primary_id == 'vagc': 
            bkgd.set_xlabel(r'log($\mathtt{sSFR_{vagc}}$ [$\mathtt{yr}^{-1}$])', fontsize=25, labelpad=25) 
            bkgd.set_ylabel(r'P($\mathtt{sSFR_{vagc}}$)', fontsize=25, labelpad=25) 
        elif primary_id == 'mpajhu': 
            bkgd.set_xlabel(r'log($\mathtt{sSFR_{mpajhu}}$ [$\mathtt{yr}^{-1}$])', fontsize=25, labelpad=25) 
            bkgd.set_ylabel(r'P($\mathtt{SSFR^{(fib)}_{mpajhu}}$)', fontsize=25, labelpad=25) 
    bkgd.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    fig.subplots_adjust(wspace=0.0)

    if neighbors == 'all': 
        rm_str = '' 
    elif neighbors == 'centrals': 
        rm_str = '.Neighbors_centrals'
        if gal_prop == 'mass': 
            sub.text(9.25, 1.1, 'Only Central Neighbors', fontsize=15) 
        elif gal_prop == 'ssfr': 
            sub.text(-12.5, 0.7, 'Only Central Neighbors', fontsize=15) 
    elif neighbors == 'pure_centrals':
        rm_str = '.Neighbors_purecentrals'
        if gal_prop == 'mass': 
            sub.text(9.25, 1.1, 'Only Pure Central Neighbors', fontsize=15) 
        elif gal_prop == 'ssfr': 
            sub.text(-12.5, 0.7, 'Only Central Neighbors', fontsize=15) 
    else: 
        raise ValueError

    cen_sat_str = '.AllPrimary'
    if cen_sat == 'centrals': 
        cen_sat_str = '.CentralPrimary'
    elif cen_sat == 'pure_centrals': 
        cen_sat_str = '.PureCentralPrimary'
    elif cen_sat == 'satellites': 
        cen_sat_str = '.SatellitePrimary'
    if neighbor_massbin is None: 
        str_neigh_M = ''
    else: 
        str_neigh_M = ''.join(['.', str(neighbor_massbin[0]), '_', str(neighbor_massbin[1])])
    fig_file = ''.join([UT.dir_fig(), 
        'Conformity', cen_sat_str ,'.P', gal_prop, '_rprep_bin', concat._FileSpec(), '.', 
        primary_id.upper(), rm_str, str_neigh_M, '.png']) 
    
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 


def PlotConformity_Primary_fsat_Rperp_bin(gal_prop, cat_dict, 
        primary_id='mpajhu', primary_massbin=[10., 10.5], cen_sat=False): 
    ''' Divide the primary galaxies into bins of SSFR and then 
    examine the SSFR of their neighboring galaxies. More specifically
    we look at the SSFR distirubiton within a r_perp bin. 
    '''
    # read in conformity catalog 
    concat = clog.ConformCatalog(Mrcut=cat_dict['Mrcut'], 
            primary_delv=cat_dict['primary_delv'], 
            primary_rperp=cat_dict['primary_rperp'],  
            neighbor_delv=cat_dict['neighbor_delv'], 
            neighbor_rperp=cat_dict['neighbor_rperp'])
    catalog = concat.Read() 

    # first sort out all primary galaxies 
    is_primary = np.where(catalog['primary_'+primary_id] == 1)[0]
    if primary_id == 'vagc': 
        ssfr_primary = catalog['ssfr'][is_primary]
        mass_primary = catalog['mass'][is_primary]
        mass_min, mass_max = 9., 12. 
        ssfr_min, ssfr_max = -13.25, -9.25
    elif primary_id == 'mpajhu': 
        ssfr_primary = catalog['ssfr_tot_mpajhu'][is_primary]
        mass_primary = catalog['mass_tot_mpajhu'][is_primary]
        mass_min, mass_max = 9., 12. 
        ssfr_min, ssfr_max = -13., -8.75

    if cen_sat == False:   # all primaries within mass bin with SSFR
        primary_cut = np.where(
                (mass_primary > primary_massbin[0]) & (mass_primary < primary_massbin[1]) & 
                (np.isnan(ssfr_primary) == False))[0]
    elif cen_sat == 'centrals':  # primaries within mass bin with SSFR that are *centrals*
        primary_cut = np.where(
                (mass_primary > primary_massbin[0]) & (mass_primary < primary_massbin[1]) & 
                (np.isnan(ssfr_primary) == False) & (catalog['p_sat'][is_primary] <= 0.5)
                )[0]
    elif cen_sat == 'pure_centrals':  # primaries within mass bin with SSFR that are pure *centrals*
        primary_cut = np.where(
                (mass_primary > primary_massbin[0]) & (mass_primary < primary_massbin[1]) & 
                (np.isnan(ssfr_primary) == False) & (catalog['p_sat'][is_primary] <= 0.01)
                )[0]
    elif cen_sat == 'satellites': # primaries within mass bin with SSFR that are *satellites*
        primary_cut = np.where(
                (mass_primary > primary_massbin[0]) & (mass_primary < primary_massbin[1]) & 
                (np.isnan(ssfr_primary) == False) & (catalog['p_sat'][is_primary] > 0.5)
                )[0]
    else: 
        raise ValueError
    
    ssfr_primary_cut = ssfr_primary[primary_cut]
    print len(primary_cut), 'primaries within cut'
    
    bin_ssfr_list = PrimarySSFR_percentile_bins(ssfr_primary_cut, 
            percentiles=[25, 50, 75, 90]) 
    bin_ssfr_label = ['0 - 25\%', '25 - 50\%', '50 - 75\%', '> 75\%', '> 90\%']
    
    rperp_bins = np.arange(0., 4.5, 0.5)
    fsat_primary_bin = [] 
    
    for bin_ssfr_cut in bin_ssfr_list: 
        # get SSFR and r_perp of the neighbors of the primaries
        # in each of the SSFR bins. 
        neigh_mass, neigh_ssfr, neigh_rperp, neigh_psat = [], [], [], []
        for ip in primary_cut[bin_ssfr_cut]: 
            cen_index = is_primary[ip]
            if cen_sat in ['central', 'pure_central']: 
                if not catalog['p_sat'][cen_index] <= 0.5:  # has to be central 
                    raise ValueError

            neigh_indices = np.array(catalog['neighbor_indices_'+primary_id][ip]).astype('int') 
            neigh_rperp.append(np.array(catalog['neighbor_rperp_'+primary_id][ip]))

            if primary_id == 'vagc': 
                neigh_ssfr.append(catalog['ssfr'][neigh_indices]) 
                neigh_mass.append(catalog['mass'][neigh_indices]) 
            elif primary_id == 'mpajhu': 
                neigh_ssfr.append(catalog['ssfr_fib_mpajhu'][neigh_indices])
                neigh_mass.append(catalog['mass_tot_mpajhu'][neigh_indices]) #neighbor *total* mass 
            neigh_psat.append(catalog['p_sat'][neigh_indices]) #neighbor p_sat 

        neigh_rperp = np.concatenate(neigh_rperp)
        neigh_ssfr = np.concatenate(neigh_ssfr) 
        neigh_mass = np.concatenate(neigh_mass) 
        neigh_psat = np.concatenate(neigh_psat) 

        # calculate f_sat of the neighbors in bins of r_perp 
        fsat_in_rperpbins = []  
        for i_rperp in range(len(rperp_bins)-1): 
            bin_prop_rperp = np.where( 
                    (neigh_rperp > rperp_bins[i_rperp]) &  
                    (neigh_rperp <= rperp_bins[i_rperp+1]) & 
                    (np.isnan(neigh_mass) == False))
            fsat = np.float(np.sum(neigh_psat[bin_prop_rperp] > 0.5))/np.float(len(neigh_psat[bin_prop_rperp]))
            fsat_in_rperpbins.append(fsat) 
        fsat_primary_bin.append(fsat_in_rperpbins)
    
    # plotting 
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure()
    sub = fig.add_subplot(111)
    for i_pb in range(len(fsat_primary_bin)):              # primary bin 
        sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), fsat_primary_bin[i_pb], 
                c=pretty_colors[i_pb], lw=2)
    # axes
    sub.set_xlim([0., 4.]) 
    sub.set_xticks([0,1,2,3,4])
    sub.set_ylim([0., 0.5]) 
    sub.set_yticks([0., 0.2, 0.4]) 
    sub.minorticks_on() 
    sub.set_xlabel(r'$\mathtt{R_{\perp}}$ [Mpc]', fontsize=25) 
    sub.set_ylabel(r'$\mathtt{f_{sat}}$', fontsize=25) 
    
    cen_sat_str = '.AllPrimary'
    if cen_sat == 'centrals': 
        cen_sat_str = '.CentralPrimary'
    elif cen_sat == 'pure_centrals': 
        cen_sat_str = '.PureCentralPrimary'
    elif cen_sat == 'satellites': 
        cen_sat_str = '.SatellitePrimary'
    fig_file = ''.join([UT.dir_fig(), 
        'Conformity', cen_sat_str ,'.fsat_rprep_bin', concat._FileSpec(), '.', 
        primary_id.upper(), '.png']) 
    
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 


def PlotConformity_Primary_meanM_Rperp_bin(gal_prop, cat_dict, 
        primary_id='mpajhu', primary_massbin=[10., 10.5], cen_sat=False): 
    ''' Divide the primary galaxies into bins of SSFR and then 
    examine the SSFR of their neighboring galaxies. More specifically
    we look at the SSFR distirubiton within a r_perp bin. 
    '''
    # read in conformity catalog 
    concat = clog.ConformCatalog(Mrcut=cat_dict['Mrcut'], 
            primary_delv=cat_dict['primary_delv'], 
            primary_rperp=cat_dict['primary_rperp'],  
            neighbor_delv=cat_dict['neighbor_delv'], 
            neighbor_rperp=cat_dict['neighbor_rperp'])
    catalog = concat.Read() 

    # first sort out all primary galaxies 
    is_primary = np.where(catalog['primary_'+primary_id] == 1)[0]
    if primary_id == 'vagc': 
        ssfr_primary = catalog['ssfr'][is_primary]
        mass_primary = catalog['mass'][is_primary]
        mass_min, mass_max = 9., 12. 
        ssfr_min, ssfr_max = -13.25, -9.25
    elif primary_id == 'mpajhu': 
        ssfr_primary = catalog['ssfr_tot_mpajhu'][is_primary]
        mass_primary = catalog['mass_tot_mpajhu'][is_primary]
        mass_min, mass_max = 9., 12. 
        ssfr_min, ssfr_max = -13., -8.75

    if cen_sat == False:   # all primaries within mass bin with SSFR
        primary_cut = np.where(
                (mass_primary > primary_massbin[0]) & (mass_primary < primary_massbin[1]) & 
                (np.isnan(ssfr_primary) == False))[0]
    elif cen_sat == 'centrals':  # primaries within mass bin with SSFR that are *centrals*
        primary_cut = np.where(
                (mass_primary > primary_massbin[0]) & (mass_primary < primary_massbin[1]) & 
                (np.isnan(ssfr_primary) == False) & (catalog['p_sat'][is_primary] <= 0.5)
                )[0]
    elif cen_sat == 'pure_centrals':  # primaries within mass bin with SSFR that are pure *centrals*
        primary_cut = np.where(
                (mass_primary > primary_massbin[0]) & (mass_primary < primary_massbin[1]) & 
                (np.isnan(ssfr_primary) == False) & (catalog['p_sat'][is_primary] <= 0.01)
                )[0]
    elif cen_sat == 'satellites': # primaries within mass bin with SSFR that are *satellites*
        primary_cut = np.where(
                (mass_primary > primary_massbin[0]) & (mass_primary < primary_massbin[1]) & 
                (np.isnan(ssfr_primary) == False) & (catalog['p_sat'][is_primary] > 0.5)
                )[0]
    else: 
        raise ValueError
    
    ssfr_primary_cut = ssfr_primary[primary_cut]
    print len(primary_cut), 'primaries within cut'
    
    bin_ssfr_list = PrimarySSFR_percentile_bins(ssfr_primary_cut, 
            percentiles=[25, 50, 75, 90]) 
    bin_ssfr_label = ['0 - 25\%', '25 - 50\%', '50 - 75\%', '> 75\%', '> 90\%']
    
    rperp_bins = np.arange(0., 4.5, 0.5)
    meanM_primary_bin = [] 
    
    for bin_ssfr_cut in bin_ssfr_list: 
        # get SSFR and r_perp of the neighbors of the primaries
        # in each of the SSFR bins. 
        neigh_mass, neigh_ssfr, neigh_rperp, neigh_psat = [], [], [], []
        for ip in primary_cut[bin_ssfr_cut]: 
            cen_index = is_primary[ip]
            if cen_sat in ['central', 'pure_central']: 
                if not catalog['p_sat'][cen_index] <= 0.5:  # has to be central 
                    raise ValueError

            neigh_indices = np.array(catalog['neighbor_indices_'+primary_id][ip]).astype('int') 
            neigh_rperp.append(np.array(catalog['neighbor_rperp_'+primary_id][ip]))

            if primary_id == 'vagc': 
                neigh_ssfr.append(catalog['ssfr'][neigh_indices]) 
                neigh_mass.append(catalog['mass'][neigh_indices]) 
            elif primary_id == 'mpajhu': 
                neigh_ssfr.append(catalog['ssfr_fib_mpajhu'][neigh_indices])
                neigh_mass.append(catalog['mass_tot_mpajhu'][neigh_indices]) #neighbor *total* mass 
            neigh_psat.append(catalog['p_sat'][neigh_indices]) #neighbor p_sat 

        neigh_rperp = np.concatenate(neigh_rperp)
        neigh_ssfr = np.concatenate(neigh_ssfr) 
        neigh_mass = np.concatenate(neigh_mass) 
        neigh_psat = np.concatenate(neigh_psat) 

        # calculate f_sat of the neighbors in bins of r_perp 
        meanM_in_rperpbins = []  
        for i_rperp in range(len(rperp_bins)-1): 
            bin_prop_rperp = np.where( 
                    (neigh_rperp > rperp_bins[i_rperp]) &  
                    (neigh_rperp <= rperp_bins[i_rperp+1]) & 
                    (np.isnan(neigh_mass) == False))
            meanM = np.mean(neigh_mass[bin_prop_rperp])
            meanM_in_rperpbins.append(meanM) 
        meanM_primary_bin.append(meanM_in_rperpbins)
    
    # plotting 
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure()
    sub = fig.add_subplot(111)
    for i_pb in range(len(meanM_primary_bin)):              # primary bin 
        sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), meanM_primary_bin[i_pb], 
                c=pretty_colors[i_pb], lw=2)
    # axes
    sub.set_xlim([0., 4.]) 
    sub.set_xticks([0,1,2,3,4])
    sub.set_ylim([10.0, 10.5]) 
    sub.set_yticks([10., 10.2, 10.4]) 
    sub.minorticks_on() 
    sub.set_xlabel(r'$\mathtt{R_{\perp}}$ [Mpc]', fontsize=25) 
    sub.set_ylabel(r'mean log($\mathcal{M_{*}} [\mathtt{M}_\odot] $)', fontsize=25) 
    
    cen_sat_str = '.AllPrimary'
    if cen_sat == 'centrals': 
        cen_sat_str = '.CentralPrimary'
    elif cen_sat == 'pure_centrals': 
        cen_sat_str = '.PureCentralPrimary'
    elif cen_sat == 'satellites': 
        cen_sat_str = '.SatellitePrimary'
    fig_file = ''.join([UT.dir_fig(), 
        'Conformity', cen_sat_str ,'.meanM_rprep_bin', concat._FileSpec(), '.', 
        primary_id.upper(), '.png']) 
    
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 


def PrimaryCentral_match(cat_dict, kauff_shenanigans=False): 
    ''' Check what fraction of the Kauffmann et al. (2013) isolation criteria 
    primaries are classified as 'centrals' in Tinker et al. (2011) group catalog. 
    '''
    # read in conformity catalog 
    if cat_dict['name'] == 'tinker':
        concat = clog.ConformCatalog(Mrcut=cat_dict['Mrcut'], 
                primary_delv=cat_dict['primary_delv'], 
                primary_rperp=cat_dict['primary_rperp'],  
                neighbor_delv=cat_dict['neighbor_delv'], 
                neighbor_rperp=cat_dict['neighbor_rperp'], 
                mpajhu=kauff_shenanigans) 
    catalog = concat.Read() 

    p_sat = catalog['p_sat']    # probability that it's a satetllite or central 
    
    # is primary 
    if not kauff_shenanigans: 
        mass_primary = catalog['mass'] 
    else: 
        mass_primary = catalog['mass_tot_mpajhu'] 
    isprimary = np.where(
            (catalog['primary'] == 1) &
            (mass_primary >= 10.) & 
            (mass_primary < 10.5)) 
    print 'Number of primaries = ', len(isprimary[0]) 
    alsocentral = np.where(p_sat[isprimary] <= 0.5) 
    print 'Number of primaries that are also centrals ', len(alsocentral[0]) 
    print 'corresponds to ', np.float(len(alsocentral[0]))/np.float(len(isprimary[0]))
    return None


def PrimaryIndices(catalog, pipeline='mpajhu', group_id='all', massbin=[10., 10.5]):  
    ''' Return indicies of specified primaries within the catalog dictionary.
    
    Parameters
    ----------
    * catalog : dictionary
        dictionary of all the conformity catalog values
    * pipeline : str
        'vagc' or 'mpajhu'
    * group_id : str
        Group catalog identification based on p_sat. 'all', 'centrals', 'pure_centrals'
    * massbin : list
        [lower mass limit, upper mass limit]
    '''
    # first sort out galaxies that are classified as primaries 
    is_primary = np.where(catalog['primary_'+pipeline] == 1)[0]
    if pipeline == 'vagc': 
        ssfr_primary = catalog['ssfr'][is_primary]
        mass_primary = catalog['mass'][is_primary]
    elif pipeline == 'mpajhu': 
        try: 
            ssfr_primary = catalog['ssfr_tot_mpajhu'][is_primary]
            mass_primary = catalog['mass_tot_mpajhu'][is_primary]   # *TOTAL* stellar mass 
        except KeyError: 
            ssfr_primary = catalog['ssfr_tot'][is_primary]
            mass_primary = catalog['mass_tot'][is_primary]   # *TOTAL* stellar mass 
    else:
        raise ValueError
    if group_id != 'all': 
        psat_primary = catalog['p_sat'][is_primary]
    
    cut_nan = (np.isnan(ssfr_primary) == False)
    cut_mass = (mass_primary > massbin[0]) & (mass_primary < massbin[1])
    if group_id == 'all':   # all primaries 
        cut_group = np.repeat(True, len(ssfr_primary))
    elif group_id == 'centrals':  # *central* primaries
        cut_group = (psat_primary <= 0.5) 
    elif group_id == 'pure_centrals':  # pure *central* primaries
        cut_group = (psat_primary <= 0.01) 
    elif group_id == 'pure_centrals_iso':  # pure *central* primaries
        cut_group = (psat_primary <= 0.1) 
    elif group_id == 'satellites': # primaries within mass bin with SSFR that are *satellites*
        cut_group = (psat_primary > 0.5) 
    elif group_id == 'not_pure_centrals': 
        cut_group = (psat_primary > 0.01) 
    elif group_id == 'not_pure_centrals_iso': 
        cut_group = (psat_primary > 0.1) 
    else: 
        print group_id
        raise ValueError
    cut_tot = np.where(cut_nan & cut_mass & cut_group)[0]
    if group_id != 'all':  
        print np.float(len(cut_tot))/np.float(len(np.where(cut_nan & cut_mass)[0])), ' are '+group_id

    return is_primary[cut_tot]


def SSFR_percentilebins(ssfrs, quantiles=None, percentiles=[25, 50, 75, 90]):
    ''' Given an array of SSFRs classify them into percentile bins. 
    Return the indices that correspond to the bins and also the SSFR 
    cut-offs.

    Parameters
    ----------
    * ssfrs : array
        Array of SSFRs
    * percentiles : list
        List that specifies the percentile of the SSFR bins. 
        [25, 50, 75, 90] --> 0-25, 25-50, 50-75, >75, >90

    Return
    ------
    * quantiles : array 
        Numpy array that specifies the SSFR cut-offs of the primary SSFR bins
    * ssfr_bin_indices : list of arrays
        List of arrays where each array corresponds to the indices 
    '''
    # calculate the percentile cut offs
    if quantiles is None: 
        #print percentiles
        quantiles = np.percentile(ssfrs, percentiles) 
    else: 
        percentiles = [25, 50, 75, 90]
        if quantiles is None: 
            raise ValueError
    #print quantiles
    
    ssfr_bin_indices = [] 
    for i_q in range(len(quantiles)+1): 
        if i_q == 0:  
            prim_inbin = np.where(ssfrs < quantiles[0])[0]
        elif percentiles[i_q-1] >= 75: 
            prim_inbin = np.where(ssfrs > quantiles[i_q-1])[0]
        else: 
            prim_inbin = np.where(
                    (ssfrs >= quantiles[i_q-1]) & 
                    (ssfrs < quantiles[i_q]))[0]
        ssfr_bin_indices.append(prim_inbin)

    return quantiles, ssfr_bin_indices



if __name__=='__main__': 
    #Plot_PrimaryBins_Pmass({'name': 'tinkauff', 'Mass_cut': 10.0, 
    #Plot_PrimaryBins_massbin_SSFR_rperp({'name': 'tinkauff', 'Mass_cut': 10.0, 
    #    'primary_delv': 500., 'primary_rperp': 0.5, 
    #    'neighbor_delv': 500., 'neighbor_rperp': 5.},
    #    primary_pipeline='mpajhu', primary_groupid='pure_centrals', primary_massbin=[10., 10.5],
    #    neighbor_pipeline='mpajhu', neighbor_groupid='all', neighbor_massbin=None) 
    
    #Plot_PrimaryBins_massbin_SSFR_rperp({'name': 'tinkauff', 'Mass_cut': 10.0, 
    #    'primary_delv': 500., 'primary_rperp': 0.5, 
    #    'neighbor_delv': 500., 'neighbor_rperp': 5.},
    #    primary_pipeline='mpajhu', primary_groupid='pure_centrals', primary_massbin=[10., 10.5],
    #    neighbor_pipeline='mpajhu', neighbor_groupid='centrals', neighbor_massbin=[10., 10.5])

    #Plot_NeighborSSFR_rperp_PrimaryBins_zSubsample(
    #Plot_NeighborSSFR_rperp_PrimaryBins(
    #        {'name': 'tinker', 'Mrcut':18, 
    #            'primary_delv': 500., 'primary_rperp': 0.5, 
    #            'neighbor_delv': 500., 'neighbor_rperp': 5.}, 
    #        rperp_bins=np.arange(0., 4.5, 0.5), 
    #        primary_pipeline='mpajhu', primary_groupid='all', primary_massbin=[10., 10.5], 
    #        neighbor_pipeline='mpajhu', neighbor_groupid='all', neighbor_massbin=None)
    #Plot_NeighborSSFR_rperp_PrimaryBins({'name': 'tinkauff', 'Mass_cut': 9.25, 
    Plot_Primary_Groups({'name': 'tinkauff', 'Mass_cut': 9.25, 
                'primary_delv': 500., 'primary_rperp': 0.5, 
                'neighbor_delv': 500., 'neighbor_rperp': 5.}, 
                primary_pipeline='mpajhu', primary_groupid='satellites', primary_massbin=[10., 10.5])
    #Plot_NeighborSSFR_rperp_PrimaryBins({'name': 'tinkauff', 'Mass_cut': 10.0, 
    #            'primary_delv': 500., 'primary_rperp': 0.5, 
    #            'neighbor_delv': 500., 'neighbor_rperp': 5.}, 
    #        rperp_bins=np.arange(0., 4.5, 0.5), 
    #        primary_pipeline='mpajhu', primary_groupid='pure_centrals', primary_massbin=[10., 10.5], 
    #        neighbor_pipeline='mpajhu', neighbor_groupid='centrals', neighbor_massbin=[10., 10.5])

    #Plot_NeighborSSFR_rperp_PrimaryBins({'name': 'kauff', 
    #            'primary_delv': 500., 'primary_rperp': 0.5, 
    #            'neighbor_delv': 500., 'neighbor_rperp': 5.}, 
    #        rperp_bins=np.arange(0., 4.5, 0.5), 
    #        primary_pipeline='mpajhu', primary_groupid='all', primary_massbin=[10., 10.5], 
    #        neighbor_pipeline='mpajhu', neighbor_groupid='all', neighbor_massbin=None)

    #PlotConformity_Primary_meanM_Rperp_bin('ssfr', {'name': 'tinker', 'Mrcut':18, 
    #    'primary_delv': 500., 'primary_rperp': 0.5, 
    #    'neighbor_delv': 500., 'neighbor_rperp': 5.}, 
    #    primary_id='mpajhu', cen_sat='pure_centrals')

    for cen_sat in ['pure_centrals']:#, 'satellites']: 
        for primary_id in ['mpajhu']: 
            for neigh in ['centrals']: #'all'
                pass
                #PlotConformity_Primary_SSFR_rperp('median', {'name': 'tinker', 'Mrcut':18, 
                #    'primary_delv': 500., 'primary_rperp': 0.5, 
                #    'neighbor_delv': 500., 'neighbor_rperp': 5.}, 
                #    primary_id=primary_id, cen_sat=cen_sat, kauff2013=True, 
                #    neighbors=neigh, neighbor_massbin=[10., 10.5])
                #PlotConformity_Primary_PDF_Rperp_bin('mass', {'name': 'tinker', 'Mrcut':18, 
                #    'primary_delv': 500., 'primary_rperp': 0.5, 
                #    'neighbor_delv': 500., 'neighbor_rperp': 5.}, 
                #    primary_id=primary_id, cen_sat=cen_sat, 
                #    neighbors=neigh, neighbor_massbin=None) 
                #PlotConformity_Primary_PDF_Rperp_bin('ssfr', {'name': 'tinker', 'Mrcut':18, 
                #    'primary_delv': 500., 'primary_rperp': 0.5, 
                #    'neighbor_delv': 500., 'neighbor_rperp': 5.}, 
                #    primary_id=primary_id, cen_sat=cen_sat, 
                #    neighbors=neigh, neighbor_massbin=[10.0, 10.5]) 
                
                #PlotConformity_Primary_PDF_Rperp_bin('ssfr', {'name': 'tinker', 'Mrcut':18, 
                #    'primary_delv': 500., 'primary_rperp': 0.5, 
                #    'neighbor_delv': 500., 'neighbor_rperp': 5.}, 
                #    primary_id=primary_id, cen_sat=cen_sat, neighbors=neigh) 
            '''
            PlotConformity_Primary_Pssfr_Rperp_bin({'name': 'tinker', 'Mrcut':18, 
                    'primary_delv': 500., 'primary_rperp': 0.5, 
                    'neighbor_delv': 500., 'neighbor_rperp': 5.}, 
                    primary_id=primary_id, cen_sat=cen_sat, remove_neighbor='none')
            PlotConformity_Primary_Pssfr_Rperp_bin({'name': 'tinker', 'Mrcut':18, 
                    'primary_delv': 500., 'primary_rperp': 0.5, 
                    'neighbor_delv': 500., 'neighbor_rperp': 5.}, 
                    primary_id=primary_id, cen_sat=cen_sat, remove_neighbor='all_satellites')
            '''
