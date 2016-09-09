import numpy as np

# --- local --- 
import util as UT
import catalog as clog

# --- plotting ---
import matplotlib.pyplot as plt
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors


def Conformity_Primary_SSFR_rperp(cat_dict, primary_id='mpajhu', primary_massbin=[10., 10.5], cen_sat=False): 
    ''' Calculate the mean and median SSFR(r_perp) of neighboring galaxies
    of primaries in different SSFR percentiles.
    '''
    # read conformity catalog based on input catalog dictionary
    concat = clog.ConformCatalog(Mrcut=cat_dict['Mrcut'], 
            primary_delv=cat_dict['primary_delv'], primary_rperp=cat_dict['primary_rperp'],  
            neighbor_delv=cat_dict['neighbor_delv'], neighbor_rperp=cat_dict['neighbor_rperp'])
    catalog = concat.Read() 

    # first sort out all primary galaxies 
    is_primary = np.where(catalog['primary_'+primary_id] == 1)[0]
    if primary_id == 'vagc': 
        ssfr_primary = catalog['ssfr'][is_primary]
        mass_primary = catalog['mass'][is_primary]
    elif primary_id == 'mpajhu': 
        #mass_primary = catalog['mass'][is_primary]
        ssfr_primary = catalog['ssfr_tot_mpajhu'][is_primary]
        mass_primary = catalog['mass_tot_mpajhu'][is_primary]
    
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

    # primary SSFR percentile bins 
    q25, q50, q75, q90 = np.percentile(ssfr_primary_cut, [25, 50, 75, 90])
    #print q25, q50, q75, q90
    bin_ssfr_0to25 = np.where(ssfr_primary_cut < q25)[0]
    bin_ssfr_25to50 = np.where((ssfr_primary_cut >= q25) & (ssfr_primary_cut < q50))[0]
    bin_ssfr_50to75 = np.where((ssfr_primary_cut >= q50) & (ssfr_primary_cut < q75))[0]
    bin_ssfr_75plus = np.where(ssfr_primary_cut >= q75)[0]
    bin_ssfr_90plus = np.where(ssfr_primary_cut >= q90)[0]
    bin_ssfr_list = [bin_ssfr_0to25, bin_ssfr_25to50, bin_ssfr_50to75, bin_ssfr_75plus, bin_ssfr_90plus]
    bin_ssfr_label = ['0 - 25\%', '25 - 50\%', '50 - 75\%', '> 75\%', '> 90\%']
    Prim_ssfrbin_Neigh_ssfr = []
    rperp_bins = np.arange(0., 4.5, 0.5)
    
    for bin_ssfr_cut in bin_ssfr_list: 
        neigh_inbin = np.concatenate(
                [np.array(catalog['neighbor_indices_'+primary_id][i]) 
                    for i in primary_cut[bin_ssfr_cut]]).astype('int') 
        neigh_rperp = np.concatenate(
                [np.array(catalog['neighbor_rperp_'+primary_id][i]) 
                    for i in primary_cut[bin_ssfr_cut]])
        if primary_id == 'vagc': 
            neigh_ssfr = catalog['ssfr'][neigh_inbin]
        elif primary_id == 'mpajhu': 
            neigh_ssfr = catalog['ssfr_fib_mpajhu'][neigh_inbin]
        
        # calculate the mean and median of *neighbor* SSFR in bins of r_perp 
        ssfr_in_rperpbins = [] 
        for i_rperp in range(len(rperp_bins)-1): 
            bin_ssfr_rperp = np.where( 
                    (neigh_rperp > rperp_bins[i_rperp]) &  
                    (neigh_rperp <= rperp_bins[i_rperp+1]) & 
                    (np.isnan(neigh_ssfr) == False))
            
            ssfr_in_rperpbins.append(neigh_ssfr[bin_ssfr_rperp])

        Prim_ssfrbin_Neigh_ssfr.append(ssfr_in_rperpbins)

    output_dict = {
            'ConCat_spec': concat._FileSpec(), 
            'ssfr_q25': q25, 'ssfr_q50': q50, 'ssfr_q75': q75, 'ssfr_q90': q90,
            'Primary_ssfr_percentile_labels': bin_ssfr_label,
            'rperp_bins': rperp_bins, 
            'percentile_Neigh_ssfr': Prim_ssfrbin_Neigh_ssfr 
            }
    return output_dict 


def PlotConformity_Primary_SSFR_rperp(ssfr_prop, cat_dict, 
        primary_id='mpajhu', primary_massbin=[10., 10.5], 
        cen_sat=False, kauff2013=False):
    ''' Plot the Conformity measurement from the Conformity_Primary_SSFR_rperp 
    function
    '''
    if ssfr_prop not in ['mean', 'median']: 
        raise ValueError

    results = Conformity_Primary_SSFR_rperp(cat_dict, 
            primary_id=primary_id, 
            primary_massbin=primary_massbin, 
            cen_sat=cen_sat)
    concat_file_spec = results['ConCat_spec']
    rperp_bins = results['rperp_bins']
    bin_ssfr_label = results['Primary_ssfr_percentile_labels']
    Prim_ssfrbin_Neigh_ssfr = results['percentile_Neigh_ssfr']
    
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure()
    sub = fig.add_subplot(111)

    for i_ssfr, ssfrs_rperp in enumerate(Prim_ssfrbin_Neigh_ssfr):
        if ssfr_prop == 'mean':   
            sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), [np.mean(ssfrs) for ssfrs in ssfrs_rperp], 
                    c=pretty_colors[i_ssfr], lw=3, ls='-', label=bin_ssfr_label[i_ssfr]) 
        elif ssfr_prop == 'median': 
            sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), [np.median(ssfrs) for ssfrs in ssfrs_rperp], 
                    c=pretty_colors[i_ssfr], lw=3, ls='-', label=bin_ssfr_label[i_ssfr]) 

        if kauff2013:    # include Kauffmann et al. (2013) data points
            kauff_str = ['0to25', '25to50', '50to75', '75plus', '90plus']
            kauff_dat_file = ''.join([UT.dir_dat(), 'literature/', 
                'kauff2013_', kauff_str[i_ssfr], '.dat']) 
            kauff_rperp, kauff_ssfr = np.loadtxt(kauff_dat_file, unpack=True, usecols=[0,1]) 
            if i_ssfr == len(kauff_str)-1: 
                kauff_label = 'Kauffmann+(2013)'
            else: 
                kauff_label = None

            sub.scatter(kauff_rperp, kauff_ssfr, c=pretty_colors[i_ssfr], lw=0) 
            sub.plot(kauff_rperp, kauff_ssfr, c=pretty_colors[i_ssfr], lw=1, ls='--', label=kauff_label) 
    
    if primary_id == 'vagc': 
        sub.set_title(r"$10.0 < \mathtt{log}\mathcal{M}^\mathtt{(vagc)}_* < 10.5$ ", 
                fontsize=25) 
    elif primary_id == 'mpajhu': 
        sub.set_title(
                r"$10.0 < \mathtt{log}\mathcal{M}^\mathtt{(tot)}_* < 10.5$ ", 
                fontsize=25) 
    # axes
    sub.set_xlim([0, 4]) 
    sub.set_xticks([0, 1, 2, 3, 4]) 
    sub.set_ylim([-12.25, -9.75])
    sub.set_yticks([-12., -11., -10.]) 
    if primary_id == 'vagc': 
        sub.text(0.1, -11.9, r'Ranked in $\mathtt{log}\mathtt{SSFR}$', fontsize=20) 
    elif primary_id == 'mpajhu': 
        sub.text(0.1, -11.9, r'Ranked in $\mathtt{log}\,\mathtt{SSFR^{(tot)}}$', fontsize=20) 
    sub.minorticks_on() 

    sub.set_xlabel(r'$\mathtt{R_{\perp}}$ [Mpc]', fontsize=25) 
    if primary_id == 'vagc': 
        sub.set_ylabel(ssfr_prop+r' log SSFR [$\mathtt{yr}^{-1}$]', fontsize=25) 
    elif primary_id == 'mpajhu': 
        sub.set_ylabel(ssfr_prop+r' log $\mathtt{SSFR^{(fib)}}$ [$\mathtt{yr}^{-1}$]', 
                fontsize=25) 
    sub.legend(loc='lower right', handletextpad=0.1) 

    str_kauff = ''
    if kauff2013:
        str_kauff = '.kauff'
    str_onlycentrals = ''
    if cen_sat == 'centrals': 
        str_onlycentrals = '.Centrals'
    elif cen_sat == 'pure_centrals': 
        str_onlycentrals = '.PureCentrals'
    elif cen_sat == 'satellites': 
        str_onlycentrals = '.Satellites'
    
    fig_file = ''.join([UT.dir_fig(), 'Conformity.SSFR_rprep', str_onlycentrals, 
        concat_file_spec, '.', primary_id.upper(), '.', ssfr_prop, str_kauff, '.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 


def PlotConformity_Primary_Pssfr_Rperp_bin(cat_dict, 
        primary_id='mpajhu', primary_massbin=[10., 10.5], cen_sat=False, remove_neighbor='none'): 
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
        ssfr_min, ssfr_max = -13.25, -9.25
    elif primary_id == 'mpajhu': 
        ssfr_primary = catalog['ssfr_tot_mpajhu'][is_primary]
        mass_primary = catalog['mass_tot_mpajhu'][is_primary]
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
    Pssfr_primary_bin = [] 
    
    rperp_bins = np.arange(0., 5., 1.)
    
    mean_ssfr_primary_bin, median_ssfr_primary_bin = [], [] 
    n_removed = []
    for bin_ssfr_cut in bin_ssfr_list: 
        # get SSFR and r_perp of the neighbors of the primaries
        # in each of the SSFR bins. 
        neigh_ssfr, neigh_rperp = [], [] 
        for ip in primary_cut[bin_ssfr_cut]: 
            cen_index = is_primary[ip]
            if cen_sat in ['central', 'pure_central']: 
                if not catalog['p_sat'][cen_index] <= 0.5:  # has to be central 
                    raise ValueError

            neigh_indices = np.array(catalog['neighbor_indices_'+primary_id][ip]).astype('int') 
            if remove_neighbor == 'none':   # keep all neighbors 
                not_sat = range(len(neigh_indices)) 
            elif remove_neighbor == 'all_satellites':   # remove all satellites 
                neigh_psat = catalog['p_sat'][neigh_indices]
                not_sat = np.where(neigh_psat <= 0.5)
                #print len(neigh_indices) - len(not_sat[0]), 'satellites removed'
                n_removed.append(len(neigh_indices) - len(not_sat[0]))
            elif remove_neighbor == 'cent_satellites':
                not_sat = np.where(catalog['id_group'][neigh_indices] != catalog['id_group'][cen_index])
                n_removed.append(len(neigh_indices) - len(not_sat[0])) #, 'satellites removed'
            else: 
                raise ValueError

            neigh_rperp.append(np.array(catalog['neighbor_rperp_'+primary_id][ip])[not_sat])

            if primary_id == 'vagc': 
                neigh_ssfr.append(catalog['ssfr'][neigh_indices[not_sat]]) 
            elif primary_id == 'mpajhu': 
                neigh_ssfr.append(catalog['ssfr_fib_mpajhu'][neigh_indices[not_sat]])
        neigh_rperp = np.concatenate(neigh_rperp)
        neigh_ssfr = np.concatenate(neigh_ssfr) 

        # calculate the mean/median SSFR in bins of r_perp 
        Pssfr_in_rperpbins = []  
        mean_ssfr_in_rperpbins, median_ssfr_in_rperpbins = [], [] 
        for i_rperp in range(len(rperp_bins)-1): 
            bin_ssfr_rperp = np.where( 
                    (neigh_rperp > rperp_bins[i_rperp]) &  
                    (neigh_rperp <= rperp_bins[i_rperp+1]) & 
                    (np.isnan(neigh_ssfr) == False))
            
            p_ssfr, bin_ssfr = np.histogram(neigh_ssfr[bin_ssfr_rperp], 
                    range=[ssfr_min, ssfr_max], bins=20, normed=True) 

            Pssfr_in_rperpbins.append(p_ssfr) 
            mean_ssfr_in_rperpbins.append(np.mean(neigh_ssfr[bin_ssfr_rperp]))
            median_ssfr_in_rperpbins.append(np.median(neigh_ssfr[bin_ssfr_rperp]))
        Pssfr_primary_bin.append(Pssfr_in_rperpbins)
        mean_ssfr_primary_bin.append(mean_ssfr_in_rperpbins)
        median_ssfr_primary_bin.append(median_ssfr_in_rperpbins)
    
    # plotting 
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(15,4))
    bkgd = fig.add_subplot(111, frameon=False) 
    for i_rp in range(len(Pssfr_primary_bin[0])):   # rperp bin 
        sub = fig.add_subplot(1,len(Pssfr_primary_bin[0]),i_rp+1)
        for i_pb in range(len(Pssfr_primary_bin)):              # primary bin 
            sub.plot(0.5*(bin_ssfr[:-1]+bin_ssfr[1:]), (Pssfr_primary_bin[i_pb])[i_rp], 
                    c=pretty_colors[i_pb], lw=2)
            sub.vlines(mean_ssfr_primary_bin[i_pb][i_rp], 0.65, 1., 
                    color=pretty_colors[i_pb], linewidth=1.5)
            sub.vlines(median_ssfr_primary_bin[i_pb][i_rp], 0., 0.15, 
                    color=pretty_colors[i_pb], linewidth=1.5)
            if i_rp == 0: 
                sub.text(ssfr_min+0.75, 0.725, 'Mean', fontsize=20)
                sub.text(ssfr_min+0.55, 0.05, 'Median', fontsize=20) 
        sub.set_title(''.join(['$', str(rperp_bins[i_rp]), '< \mathtt{r}_\perp <', str(rperp_bins[i_rp+1]), '$']), 
                fontsize=20) 
        # axes
        sub.set_xlim([ssfr_min, ssfr_max]) 
        if i_rp < len(Pssfr_primary_bin[0])-1: 
            sub.set_xticks([-13, -12, -11, -10]) 
        else: 
            sub.set_xticks([-13, -12, -11, -10, -9]) 
        sub.set_ylim([0., 0.8]) 
        sub.set_yticks([0., 0.2, 0.4, 0.6, 0.8]) 
        sub.minorticks_on() 
        if i_rp != 0: 
            sub.set_yticklabels([]) 

    bkgd.set_xticklabels([]) 
    bkgd.set_yticklabels([]) 
    if primary_id == 'vagc': 
        bkgd.set_xlabel(r'log($\mathtt{sSFR_{vagc}}$ [$\mathtt{yr}^{-1}$])', fontsize=25, labelpad=25) 
        bkgd.set_ylabel(r'P($\mathtt{sSFR_{vagc}}$)', fontsize=25, labelpad=25) 
    elif primary_id == 'mpajhu': 
        bkgd.set_xlabel(r'log($\mathtt{sSFR_{mpajhu}}$ [$\mathtt{yr}^{-1}$])', fontsize=25, labelpad=25) 
        bkgd.set_ylabel(r'P($\mathtt{SSFR^{(fib)}_{mpajhu}}$)', fontsize=25, labelpad=25) 
    bkgd.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    fig.subplots_adjust(wspace=0.0)

    if remove_neighbor == 'none': 
        rm_str = '' 
    elif remove_neighbor == 'all_satellites': 
        rm_str = '.allsatellites_removed'
        sub.text(-12.5, 0.7, 'All Satellites Removed', fontsize=15) 
    elif remove_neighbor == 'cent_satellites':
        rm_str = '.centsatellites_removed'
        sub.text(-12.5, 1.5, 'Satellites of Primary Removed', fontsize=15) 
    else: 
        raise ValueError

    cen_sat_str = '.AllPrimary'
    if cen_sat == 'centrals': 
        cen_sat_str = '.CentralPrimary'
    elif cen_sat == 'pure_centrals': 
        cen_sat_str = '.PureCentralPrimary'
    elif cen_sat == 'satellites': 
        cen_sat_str = '.SatellitePrimary'

    fig_file = ''.join([UT.dir_fig(), 
        'Conformity', cen_sat_str ,'.Pssfr_rprep_bin', concat._FileSpec(), '.', 
        primary_id.upper(), rm_str, '.png']) 
    
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



if __name__=='__main__': 
    for cen_sat in [False, 'centrals', 'pure_centrals', 'satellites']: 
        for primary_id in ['vagc', 'mpajhu']: 
            '''
            PlotConformity_Primary_SSFR_rperp('median', {'name': 'tinker', 'Mrcut':18, 
                'primary_delv': 500., 'primary_rperp': 0.5, 
                'neighbor_delv': 500., 'neighbor_rperp': 5.}, 
                primary_id=primary_id, cen_sat=cen_sat, kauff2013=True)

            PlotConformity_Primary_Pssfr_Rperp_bin({'name': 'tinker', 'Mrcut':18, 
                    'primary_delv': 500., 'primary_rperp': 0.5, 
                    'neighbor_delv': 500., 'neighbor_rperp': 5.}, 
                    primary_id=primary_id, cen_sat=cen_sat, remove_neighbor='none')
            '''
            PlotConformity_Primary_Pssfr_Rperp_bin({'name': 'tinker', 'Mrcut':18, 
                    'primary_delv': 500., 'primary_rperp': 0.5, 
                    'neighbor_delv': 500., 'neighbor_rperp': 5.}, 
                    primary_id=primary_id, cen_sat=cen_sat, remove_neighbor='all_satellites')

"""
    def Conformity_Centrals_SSFR_rperp(cat_dict, primary_id='mpajhu', meanmedian='mean', 
            remove_neighbor='none'): 
        ''' Conformity signal in the SSFR(r_perp) for central galaxies binned 
        by their SSFR. 
        '''
        # read in conformity catalog 
        if cat_dict['name'] == 'tinker':
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
        elif primary_id == 'mpajhu': 
            ssfr_primary = catalog['ssfr_tot_mpajhu'][is_primary]
            mass_primary = catalog['mass_tot_mpajhu'][is_primary]
        psat_primary = catalog['p_sat'][is_primary] 
        # get primaries within mass bin 
        primary_cut = np.where((mass_primary > 10.) & (mass_primary < 10.5) & 
                (np.isnan(ssfr_primary) == False) & (psat_primary <= 0.5))[0]
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
        bin_ssfr_dist = [] 
        
        rperp_bins = np.arange(0., 4.5, 0.5)
        
        n_removed = []
        for bin_ssfr_cut in bin_ssfr_list: 
            # loop through the neighbors of the central
            neigh_ssfr, neigh_rperp = [], [] 
            for ip in primary_cut[bin_ssfr_cut]: 
                cen_index = is_primary[ip]
                if not catalog['p_sat'][cen_index] <= 0.5:  # has to be central 
                    raise ValueError

                neigh_indices = np.array(catalog['neighbor_indices_'+primary_id][ip]).astype('int') 
                if remove_neighbor == 'none':   # keep all neighbors 
                    not_sat = range(len(neigh_indices)) 
                elif remove_neighbor == 'all_satellites':   # remove all satellites 
                    neigh_psat = catalog['p_sat'][neigh_indices]
                    not_sat = np.where(neigh_psat <= 0.5)
                    #print len(neigh_indices) - len(not_sat[0]), 'satellites removed'
                    n_removed.append(len(neigh_indices) - len(not_sat[0]))
                elif remove_neighbor == 'cent_satellites':
                    not_sat = np.where(catalog['id_group'][neigh_indices] != catalog['id_group'][cen_index])
                    n_removed.append(len(neigh_indices) - len(not_sat[0])) #, 'satellites removed'
                else: 
                    raise ValueError

                neigh_rperp.append(np.array(catalog['neighbor_rperp_'+primary_id][ip])[not_sat])

                #print catalog['id_group'][cen_index] #catalog['p_sat'][cen_index]
                #print catalog['id_group'][neigh_indices] #catalog['p_sat'][neigh_indices]
                #print np.sum(catalog['id_cent'][neigh_indices] == catalog['id_gal'][cen_index])
                #print np.sum(catalog['id_group'][neigh_indices] == catalog['id_group'][cen_index])
                #print np.sum(catalog['p_sat'][neigh_indices] > 0.5) 
        
                if primary_id == 'vagc': 
                    neigh_ssfr.append(catalog['ssfr'][neigh_indices[not_sat]]) 
                elif primary_id == 'mpajhu': 
                    neigh_ssfr.append(catalog['ssfr_fib_mpajhu'][neigh_indices[not_sat]])
            neigh_rperp = np.concatenate(neigh_rperp)
            neigh_ssfr = np.concatenate(neigh_ssfr) 

            # calculate the mean/median SSFR in bins of r_perp 
            ssfr_in_rperpbins = np.zeros(len(rperp_bins)-1)
            for i_rperp in range(len(rperp_bins)-1): 
                bin_ssfr_rperp = np.where( 
                        (neigh_rperp > rperp_bins[i_rperp]) &  
                        (neigh_rperp <= rperp_bins[i_rperp+1]) & 
                        (np.isnan(neigh_ssfr) == False))

                if meanmedian == 'mean': 
                    ssfr_in_rperpbins[i_rperp] = np.mean(neigh_ssfr[bin_ssfr_rperp])
                elif meanmedian == 'median': 
                    ssfr_in_rperpbins[i_rperp] = np.median(neigh_ssfr[bin_ssfr_rperp])
                else: 
                    raise ValueError
            bin_ssfr_dist.append(ssfr_in_rperpbins)
        
        print np.max(n_removed) 
        print np.float(len(np.where(np.array(n_removed) > 0)[0]))/np.float(len(n_removed)), ' fraction removed'
        
        prettyplot()
        pretty_colors = prettycolors()
        fig = plt.figure()
        sub = fig.add_subplot(111)
        for i_ssfr, ssfr_dist in enumerate(bin_ssfr_dist): 
            if i_ssfr in [1,2]: 
                lstyle = '--' 
            else: 
                lstyle = '-'
            sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_dist, 
                    c=pretty_colors[i_ssfr], lw=3, ls=lstyle, label=bin_ssfr_label[i_ssfr]) 
            
        if primary_id == 'vagc': 
            sub.set_title(r"$10.0 < \mathtt{log}\mathcal{M}^\mathtt{(vagc)}_* < 10.5$ ", 
                    fontsize=25) 
        elif primary_id == 'mpajhu': 
            sub.set_title(
                    r"$10.0 < \mathtt{log}\mathcal{M}^\mathtt{(tot)}_* < 10.5$ ", 
                    fontsize=25) 
        # axes
        sub.set_xlim([0, 4]) 
        sub.set_xticks([0, 1, 2, 3, 4]) 
        sub.set_ylim([-12.25, -9.75])
        sub.set_yticks([-12., -11., -10.]) 
        if primary_id == 'vagc': 
            sub.text(0.1, -11.9, r'Ranked in $\mathtt{log}\mathtt{SSFR}$', fontsize=20) 
        elif primary_id == 'mpajhu': 
            sub.text(0.1, -11.9, r'Ranked in $\mathtt{log}\,\mathtt{SSFR^{(tot)}}$', fontsize=20) 
        sub.minorticks_on() 

        sub.set_xlabel(r'$\mathtt{R_{\perp}}$ [Mpc]', fontsize=25) 
        if primary_id == 'vagc': 
            sub.set_ylabel(meanmedian+r' log SSFR [$\mathtt{yr}^{-1}$]', fontsize=25) 
        elif primary_id == 'mpajhu': 
            sub.set_ylabel(meanmedian+r' log $\mathtt{SSFR^{(fib)}}$ [$\mathtt{yr}^{-1}$]', 
                    fontsize=25) 
        sub.legend(loc='lower right', handletextpad=0.1) 
        
        if remove_neighbor == 'none': 
            rm_str = '' 
        elif remove_neighbor == 'all_satellites': 
            rm_str = '.allsatellites_removed'
            sub.text(2.0, -10., 'All Satellites Removed', fontsize=15) 
        elif remove_neighbor == 'cent_satellites':
            rm_str = '.centsatellites_removed'
            sub.text(1.5, -10., 'Satellites of Primary Removed', fontsize=15) 
        else: 
            raise ValueError

        fig_file = ''.join([UT.dir_fig(), 
            'Conformity.CentralsOnly.SSFR_rprep', concat._FileSpec(), '.', 
            primary_id.upper(), '.', meanmedian, rm_str, '.png']) 
        
        fig.savefig(fig_file, bbox_inches='tight') 
        plt.close()
        return None 
"""
