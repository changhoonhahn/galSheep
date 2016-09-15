import numpy as np

# --- local --- 
import util as UT
import catalog as clog

# --- plotting ---
import matplotlib.pyplot as plt
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors


def Conformity_Primary_SSFR_rperp(cat_dict, primary_id='mpajhu', primary_massbin=[10., 10.5], 
        cen_sat=False, neighbors='all', neighbor_massbin=None): 
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
            neigh_mass = catalog['mass'][neigh_inbin]
        elif primary_id == 'mpajhu': 
            neigh_ssfr = catalog['ssfr_fib_mpajhu'][neigh_inbin]
            neigh_mass = catalog['mass_tot_mpajhu'][neigh_inbin] #neighbor *total* mass 
        neigh_psat = catalog['p_sat'][neigh_inbin]    # satellite probability of neighbors

        # calculate the mean and median of *neighbor* SSFR in bins of r_perp 
        ssfr_in_rperpbins = [] 
        for i_rperp in range(len(rperp_bins)-1): 
            bin_ssfr_rperp = (neigh_rperp > rperp_bins[i_rperp]) & \
                    (neigh_rperp <= rperp_bins[i_rperp+1]) & \
                    (np.isnan(neigh_ssfr) == False)
            if neighbors == 'all': # include all neighbors 
                neigh_cut = np.repeat(True, len(bin_ssfr_rperp))
            elif neighbors == 'centrals': # only include *central* neighbors 
                neigh_cut = (neigh_psat <= 0.5)
            elif neighbors == 'pure_centrals': # only include *pure central* neighbors 
                neigh_cut = (neigh_psat <= 0.01)
            else: 
                raise ValueError
            if neighbor_massbin is None: 
                neigh_mass_cut = np.repeat(True, len(bin_ssfr_rperp))
            else: 
                neigh_mass_cut = (neigh_mass >= neighbor_massbin[0]) & (neigh_mass < neighbor_massbin[1])

            tot_cuts = np.where(bin_ssfr_rperp & neigh_cut & neigh_mass_cut) 
            
            ssfr_in_rperpbins.append(neigh_ssfr[tot_cuts])

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
        cen_sat=False, kauff2013=False, neighbors='all', neighbor_massbin=None):
    ''' Plot the Conformity measurement from the Conformity_Primary_SSFR_rperp 
    function
    '''
    if ssfr_prop not in ['mean', 'median']: 
        raise ValueError

    results = Conformity_Primary_SSFR_rperp(cat_dict, 
            primary_id=primary_id, 
            primary_massbin=primary_massbin, 
            cen_sat=cen_sat, 
            neighbors=neighbors, 
            neighbor_massbin=neighbor_massbin)
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
        sub.set_title(r"$10.0 < \mathtt{log}\mathcal{M}^\mathtt{(tot)}_* < 10.5$ ", 
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

    if neighbors == 'centrals': 
        sub.text(1.3, -10., 'Only Central Neighbors', fontsize=20)
    elif neighbors == 'pure_centrals': 
        sub.text(1.3, -10., 'Only Pure Central Neighbors', fontsize=20)
    if neighbor_massbin is not None:  
        sub.text(1.3, -10.2, str(neighbor_massbin[0])+'$< $log$\mathcal{M}_*^\mathtt{neigh} <$'+str(neighbor_massbin[1]), fontsize=20)
    sub.minorticks_on() 

    sub.set_xlabel(r'$\mathtt{R_{\perp}}$ [Mpc]', fontsize=25) 
    if primary_id == 'vagc': 
        sub.set_ylabel(ssfr_prop+r' log SSFR [$\mathtt{yr}^{-1}$]', fontsize=25) 
    elif primary_id == 'mpajhu': 
        sub.set_ylabel(ssfr_prop+r' log $\mathtt{SSFR^{(fib)}}$ [$\mathtt{yr}^{-1}$]', 
                fontsize=25) 
    sub.legend(loc='lower right', handletextpad=0.1) 
    
    str_kauff = ''   # whether or not kauffmann et al. (2013) is included or not 
    if kauff2013:
        str_kauff = '.kauff'
    str_onlycentrals = ''       # primary category 
    if cen_sat == 'centrals':   
        str_onlycentrals = '.Centrals'
    elif cen_sat == 'pure_centrals': 
        str_onlycentrals = '.PureCentrals'
    elif cen_sat == 'satellites': 
        str_onlycentrals = '.Satellites'
    if neighbors == 'all':      # neighbor catalogy
        neigh_str = ''
    elif neighbors == 'centrals': 
        neigh_str = '.Neighbors_centrals'
    elif neighbors == 'pure_centrals': 
        neigh_str = '.Neighbors_purecentrals'
    if neighbor_massbin is None:  
        neigh_massbin = '' 
    else: 
        neigh_massbin = '_'.join([str(neighbor_massbin[0]), str(neighbor_massbin[1])]) 
    fig_file = ''.join([UT.dir_fig(), 'Conformity.SSFR_rprep', str_onlycentrals, 
        concat_file_spec, '.', primary_id.upper(), '.', ssfr_prop, 
        neigh_str, neigh_massbin, str_kauff, '.png']) 
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
                neigh_ssfr.append(catalog['ssfr_fib_mpajhu'][neigh_indices[neigh_cut]])
                neigh_mass.append(catalog['mass_tot_mpajhu'][neigh_indices[neigh_cut]]) #neighbor *total* mass 

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


def PrimarySSFR_percentile_bins(ssfrs, percentiles=[25, 50, 75, 90]):
    ''' Given SSFRs and percentiles, return a list of indices for the SSFR
    that falls into the quantile range 
    '''
    quantiles = np.percentile(ssfrs, percentiles) 
    
    ssfr_bin_indices = [] 
    for i_q in range(len(quantiles)+1): 
        if i_q == 0:  
            ssfr_bin_indices.append(np.where(ssfrs < quantiles[i_q])[0])
        elif percentiles[i_q-1] >= 75: 
            ssfr_bin_indices.append(np.where(ssfrs >= quantiles[i_q-1])[0])
        else: 
            ssfr_bin_indices.append(
                    np.where((ssfrs >= quantiles[i_q-1]) & 
                        (ssfrs < quantiles[i_q]))[0])
    
    return ssfr_bin_indices 



if __name__=='__main__': 

    PlotConformity_Primary_meanM_Rperp_bin('ssfr', {'name': 'tinker', 'Mrcut':18, 
        'primary_delv': 500., 'primary_rperp': 0.5, 
        'neighbor_delv': 500., 'neighbor_rperp': 5.}, 
        primary_id='mpajhu', cen_sat='pure_centrals')

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
