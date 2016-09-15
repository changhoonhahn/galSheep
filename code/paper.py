'''

Figures for the conformity paper 


'''
import numpy as np

# --- local --- 
import util as UT
import catalog as clog
import conformity as conform

# --- plotting ---
import matplotlib.pyplot as plt
from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors


def Fig_NeighborSSFR_in_Primarybins(cat_dict, primary_massbin=[10., 10.5]): 
    ''' Figure plotting the median SSFR of neighbors as a function of r_perp 
    for primary galaxies binned in SSFR. Direct comparison to Kauffmann et al.(2013) 
    figure. This is the most *raw* measurement. It contains *all* primary galaxies 
    and *all* neighbor galaxies.

    Notes
    -----
    * We want an apples to apples comparison, so we use MPA-JHU derived values 
    rather than the VAGC
    '''
    # read conformity catalog based on input catalog dictionary
    concat = clog.ConformCatalog(Mrcut=cat_dict['Mrcut'], 
            primary_delv=cat_dict['primary_delv'], primary_rperp=cat_dict['primary_rperp'],  
            neighbor_delv=cat_dict['neighbor_delv'], neighbor_rperp=cat_dict['neighbor_rperp'])
    catalog = concat.Read() 

    # first sort out all primary galaxies 
    is_primary = np.where(catalog['primary_mpajhu'] == 1)[0]
    ssfr_primary = catalog['ssfr_tot_mpajhu'][is_primary]
    mass_primary = catalog['mass_tot_mpajhu'][is_primary]
    
    # all primaries within mass bin with SSFR
    primary_cut = np.where(
            (mass_primary > primary_massbin[0]) & (mass_primary < primary_massbin[1]) & 
            (np.isnan(ssfr_primary) == False))[0]
    ssfr_primary_cut = ssfr_primary[primary_cut]
    print len(primary_cut), 'primaries within cut'

    # primary SSFR percentile bins 
    bin_ssfr_list = conform.PrimarySSFR_percentile_bins(ssfr_primary_cut, 
            percentiles=[25, 50, 75, 90]) 
    bin_ssfr_label = ['0 - 25\%', '25 - 50\%', '50 - 75\%', '> 75\%', '> 90\%']

    NeighSSFR_primarybin = [] 
    rperp_bins = np.arange(0., 4.5, 0.5)
    
    for bin_ssfr_cut in bin_ssfr_list: 
        neigh_inbin = np.concatenate(    # neighbor in primary bin 
                [np.array(catalog['neighbor_indices_mpajhu'][i]) 
                    for i in primary_cut[bin_ssfr_cut]]).astype('int') 
        neigh_rperp = np.concatenate(
                [np.array(catalog['neighbor_rperp_mpajhu'][i]) 
                    for i in primary_cut[bin_ssfr_cut]])
        neigh_ssfr = catalog['ssfr_fib_mpajhu'][neigh_inbin]
        neigh_mass = catalog['mass_tot_mpajhu'][neigh_inbin] #neighbor *total* mass 
        neigh_psat = catalog['p_sat'][neigh_inbin]    # satellite probability of neighbors

        # calculate the mean and median of *neighbor* SSFR in bins of r_perp 
        neighSSFR_rperpbins = [] 
        for i_rperp in range(len(rperp_bins)-1): 
            rperp_cut = (neigh_rperp > rperp_bins[i_rperp]) & \
                    (neigh_rperp <= rperp_bins[i_rperp+1]) & \
                    (np.isnan(neigh_ssfr) == False)
            tot_cuts = np.where(rperp_cut) 
            neighSSFR_rperpbins.append(neigh_ssfr[tot_cuts])
        NeighSSFR_primarybin.append(neighSSFR_rperpbins)

    # Kauffmann et al.(2013) local files (data thief-ed)
    str_kauff = ['0to25', '25to50', '50to75', '75plus', '90plus']
    kauff_dat_file = lambda spec: ''.join([UT.dir_dat(), 'literature/kauff2013_', spec, '.dat']) 

    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure()
    sub = fig.add_subplot(111)
    ssfrplots = [] 
    for i_ssfr, ssfrs_rperp in enumerate(NeighSSFR_primarybin):
        # Plot Kauffmann et al.(2013) data points 
        kauff_rperp, kauff_ssfr = np.loadtxt(kauff_dat_file(str_kauff[i_ssfr]), 
                unpack=True, usecols=[0,1]) 
        sub.scatter(kauff_rperp, kauff_ssfr, c=pretty_colors[i_ssfr], lw=0) 
        kauffplot, = sub.plot(kauff_rperp, kauff_ssfr, c=pretty_colors[i_ssfr], 
                lw=1, ls='--', label='Kauffmann+(2013)') 

        ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), 
                [np.median(ssfrs) for ssfrs in ssfrs_rperp], 
                c=pretty_colors[i_ssfr], lw=3, ls='-', label=bin_ssfr_label[i_ssfr]) 
        ssfrplots.append(ssfrplot) 
    
    sub.set_title(r"$10.0 < \mathtt{log}\mathcal{M}^\mathtt{(tot)}_* < 10.5$ ", fontsize=25) 
    # axes
    sub.set_xlim([0, 4]) 
    sub.set_xticks([0, 1, 2, 3, 4]) 
    sub.set_ylim([-12.25, -9.75])
    sub.set_yticks([-12., -11., -10.]) 
    sub.text(0.35, -11.8, r'Primaries ranked', fontsize=20) 
    sub.text(0.35, -12., r'in $\mathtt{log}\,\mathtt{SSFR^{(tot)}}$', fontsize=20) 
    sub.minorticks_on() 

    sub.set_xlabel(r'$\mathtt{R_{\perp}}$ [Mpc]', fontsize=25) 
    sub.set_ylabel(r'median log $\mathtt{SSFR^{(fib)}_{neigh}}$ [$\mathtt{yr}^{-1}$]', 
            fontsize=25) 
    first_legend = sub.legend(handles=[kauffplot], loc='lower right', handletextpad=0.1)
    ax = plt.gca().add_artist(first_legend)
    sub.legend(handles=ssfrplots, loc='upper right', ncol=2, handletextpad=0.1) 

    fig_file = ''.join([UT.dir_fig(), 
        'Conformity.NeighborSSFR_in_PrimarySSFR_bins', concat._FileSpec(), '.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 


def Fig_NeighborSSFR_in_PurePrimarybins(cat_dict, primary_massbin=[10., 10.5]): 
    ''' Figure plotting the median SSFR of neighbors as a function of r_perp 
    for *PURE CENTRAL* primary galaxies binned in SSFR along with the 
    "raw measurement".
    
    It contains *pure central* primary galaxies and *all* neighbor galaxies.

    Notes
    -----
    * We want an apples to apples comparison, so we use MPA-JHU derived values 
    rather than the VAGC
    '''
    # read conformity catalog based on input catalog dictionary
    concat = clog.ConformCatalog(Mrcut=cat_dict['Mrcut'], 
            primary_delv=cat_dict['primary_delv'], primary_rperp=cat_dict['primary_rperp'],  
            neighbor_delv=cat_dict['neighbor_delv'], neighbor_rperp=cat_dict['neighbor_rperp'])
    catalog = concat.Read() 

    # first sort out all primary galaxies 
    is_primary = np.where(catalog['primary_mpajhu'] == 1)[0]
    ssfr_primary = catalog['ssfr_tot_mpajhu'][is_primary]
    mass_primary = catalog['mass_tot_mpajhu'][is_primary]
    
    # all primaries within mass bin with SSFR
    primary_cut = np.where(
            (mass_primary > primary_massbin[0]) & (mass_primary < primary_massbin[1]) & 
            (np.isnan(ssfr_primary) == False))[0]
    ssfr_primary_cut = ssfr_primary[primary_cut]
    print len(primary_cut), 'primaries within cut'
    # PURE central primaries within mass bin with SSFR 
    pure_primary_cut = np.where(
            (mass_primary > primary_massbin[0]) & (mass_primary < primary_massbin[1]) & 
            (np.isnan(ssfr_primary) == False) & (catalog['p_sat'][is_primary] <= 0.01))[0]
    ssfr_pure_primary_cut = ssfr_primary[pure_primary_cut]
    print len(pure_primary_cut), 'pure central primaries within cut'

    # primary SSFR percentile bins 
    bin_ssfr_list = conform.PrimarySSFR_percentile_bins(ssfr_primary_cut, 
            percentiles=[25, 50, 75, 90]) 
    bin_pure_ssfr_list = conform.PrimarySSFR_percentile_bins(ssfr_pure_primary_cut, 
            percentiles=[25, 50, 75, 90]) 
    bin_ssfr_label = ['0 - 25\%', '25 - 50\%', '50 - 75\%', '> 75\%', '> 90\%']

    rperp_bins = np.arange(0., 4.5, 0.5)
    NeighSSFR_primarybin, NeighSSFR_pureprimarybin = [], []

    for i_bin, bin_ssfr_cut in enumerate(bin_ssfr_list): 
        neigh_inbin = np.concatenate(    # neighbor in primary bin 
                [np.array(catalog['neighbor_indices_mpajhu'][i]) 
                    for i in primary_cut[bin_ssfr_cut]]).astype('int') 
        neigh_rperp = np.concatenate(
                [np.array(catalog['neighbor_rperp_mpajhu'][i]) 
                    for i in primary_cut[bin_ssfr_cut]])
        neigh_ssfr = catalog['ssfr_fib_mpajhu'][neigh_inbin]
        neigh_mass = catalog['mass_tot_mpajhu'][neigh_inbin] #neighbor *total* mass 
        neigh_psat = catalog['p_sat'][neigh_inbin]    # satellite probability of neighbors

        pure_neigh_inbin = np.concatenate(    # neighbor in primary bin 
                [np.array(catalog['neighbor_indices_mpajhu'][i]) 
                    for i in pure_primary_cut[bin_pure_ssfr_list[i_bin]]]).astype('int') 
        pure_neigh_rperp = np.concatenate(
                [np.array(catalog['neighbor_rperp_mpajhu'][i]) 
                    for i in pure_primary_cut[bin_pure_ssfr_list[i_bin]]])
        pure_neigh_ssfr = catalog['ssfr_fib_mpajhu'][pure_neigh_inbin]
        pure_neigh_mass = catalog['mass_tot_mpajhu'][pure_neigh_inbin] #neighbor *total* mass 
        pure_neigh_psat = catalog['p_sat'][pure_neigh_inbin]    # satellite probability of neighbors

        # calculate the mean and median of *neighbor* SSFR in bins of r_perp 
        neighSSFR_rperpbins, pure_neighSSFR_rperpbins = [], [] 
        for i_rperp in range(len(rperp_bins)-1): 
            rperp_cut = (neigh_rperp > rperp_bins[i_rperp]) & \
                    (neigh_rperp <= rperp_bins[i_rperp+1]) & \
                    (np.isnan(neigh_ssfr) == False)
            tot_cuts = np.where(rperp_cut) 
            neighSSFR_rperpbins.append(neigh_ssfr[tot_cuts])

            pure_rperp_cut = (pure_neigh_rperp > rperp_bins[i_rperp]) & \
                    (pure_neigh_rperp <= rperp_bins[i_rperp+1]) & \
                    (np.isnan(pure_neigh_ssfr) == False)
            pure_tot_cuts = np.where(pure_rperp_cut) 
            pure_neighSSFR_rperpbins.append(pure_neigh_ssfr[pure_tot_cuts])

        NeighSSFR_primarybin.append(neighSSFR_rperpbins)
        NeighSSFR_pureprimarybin.append(pure_neighSSFR_rperpbins)

    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure()
    sub = fig.add_subplot(111)
    ssfrplots = [] 
    for i_ssfr, ssfrs_rperp in enumerate(NeighSSFR_primarybin):
        sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), 
                [np.median(ssfrs) for ssfrs in ssfrs_rperp], 
                c=pretty_colors[i_ssfr], lw=1, ls='--') 
    
        pure_ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), 
                [np.median(ssfrs) for ssfrs in NeighSSFR_pureprimarybin[i_ssfr]], 
                c=pretty_colors[i_ssfr], lw=3, ls='-', label=bin_ssfr_label[i_ssfr]) 
        ssfrplots.append(pure_ssfrplot) 

    sub.set_title(r"$10.0 < \mathtt{log}\mathcal{M}^\mathtt{(tot)}_* < 10.5$ ", fontsize=25) 
    # axes
    sub.set_xlim([0, 4]) 
    sub.set_xticks([0, 1, 2, 3, 4]) 
    sub.set_ylim([-12.25, -9.75])
    sub.set_yticks([-12., -11., -10.]) 
    sub.text(0.35, -11.8, r'Pure primaries', fontsize=20) 
    sub.text(0.35, -12., r'ranked in $\mathtt{log}\,\mathtt{SSFR^{(tot)}}$', fontsize=20) 
    sub.minorticks_on() 

    sub.set_xlabel(r'$\mathtt{R_{\perp}}$ [Mpc]', fontsize=25) 
    sub.set_ylabel(r'median log $\mathtt{SSFR^{(fib)}_{neigh}}$ [$\mathtt{yr}^{-1}$]', 
            fontsize=25) 
    first_legend = sub.legend(handles=ssfrplots, loc='upper right', ncol=2, handletextpad=0.1) 
    fig_file = ''.join([UT.dir_fig(), 
        'Conformity.NeighborSSFR_in_PurePrimarySSFR_bins', concat._FileSpec(), '.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 


def Fig_CentralNeighborSSFR_in_PurePrimarybins(cat_dict, primary_massbin=[10., 10.5]): 
    ''' Figure plotting the median SSFR of *central* neighbors as a function of r_perp 
    for *PURE CENTRAL* primary galaxies binned in SSFR along with the "raw measurement".
    
    It contains *pure central* primary galaxies and *central* neighbor galaxies.

    Notes
    -----
    * We want an apples to apples comparison, so we use MPA-JHU derived values 
    rather than the VAGC
    '''
    # read conformity catalog based on input catalog dictionary
    concat = clog.ConformCatalog(Mrcut=cat_dict['Mrcut'], 
            primary_delv=cat_dict['primary_delv'], primary_rperp=cat_dict['primary_rperp'],  
            neighbor_delv=cat_dict['neighbor_delv'], neighbor_rperp=cat_dict['neighbor_rperp'])
    catalog = concat.Read() 

    # first sort out all primary galaxies 
    is_primary = np.where(catalog['primary_mpajhu'] == 1)[0]
    ssfr_primary = catalog['ssfr_tot_mpajhu'][is_primary]
    mass_primary = catalog['mass_tot_mpajhu'][is_primary]
    
    # all primaries within mass bin with SSFR
    primary_cut = np.where(
            (mass_primary > primary_massbin[0]) & (mass_primary < primary_massbin[1]) & 
            (np.isnan(ssfr_primary) == False))[0]
    ssfr_primary_cut = ssfr_primary[primary_cut]
    print len(primary_cut), 'primaries within cut'
    # PURE central primaries within mass bin with SSFR 
    pure_primary_cut = np.where(
            (mass_primary > primary_massbin[0]) & (mass_primary < primary_massbin[1]) & 
            (np.isnan(ssfr_primary) == False) & (catalog['p_sat'][is_primary] <= 0.01))[0]
    ssfr_pure_primary_cut = ssfr_primary[pure_primary_cut]
    print len(pure_primary_cut), 'pure central primaries within cut'

    # primary SSFR percentile bins 
    bin_ssfr_list = conform.PrimarySSFR_percentile_bins(ssfr_primary_cut, 
            percentiles=[25, 50, 75, 90]) 
    bin_pure_ssfr_list = conform.PrimarySSFR_percentile_bins(ssfr_pure_primary_cut, 
            percentiles=[25, 50, 75, 90]) 
    bin_ssfr_label = ['0 - 25\%', '25 - 50\%', '50 - 75\%', '> 75\%', '> 90\%']

    rperp_bins = np.arange(0., 4.5, 0.5)
    NeighSSFR_primarybin, NeighSSFR_pureprimarybin = [], []

    for i_bin, bin_ssfr_cut in enumerate(bin_ssfr_list): 
        neigh_inbin = np.concatenate(    # neighbor in primary bin 
                [np.array(catalog['neighbor_indices_mpajhu'][i]) 
                    for i in primary_cut[bin_ssfr_cut]]).astype('int') 
        neigh_rperp = np.concatenate(
                [np.array(catalog['neighbor_rperp_mpajhu'][i]) 
                    for i in primary_cut[bin_ssfr_cut]])
        neigh_ssfr = catalog['ssfr_fib_mpajhu'][neigh_inbin]
        neigh_mass = catalog['mass_tot_mpajhu'][neigh_inbin] #neighbor *total* mass 
        neigh_psat = catalog['p_sat'][neigh_inbin]    # satellite probability of neighbors

        pure_neigh_inbin = np.concatenate(    # neighbor in primary bin 
                [np.array(catalog['neighbor_indices_mpajhu'][i]) 
                    for i in pure_primary_cut[bin_pure_ssfr_list[i_bin]]]).astype('int') 
        pure_neigh_rperp = np.concatenate(
                [np.array(catalog['neighbor_rperp_mpajhu'][i]) 
                    for i in pure_primary_cut[bin_pure_ssfr_list[i_bin]]])
        pure_neigh_ssfr = catalog['ssfr_fib_mpajhu'][pure_neigh_inbin]
        pure_neigh_mass = catalog['mass_tot_mpajhu'][pure_neigh_inbin] #neighbor *total* mass 
        pure_neigh_psat = catalog['p_sat'][pure_neigh_inbin]    # satellite probability of neighbors

        # calculate the mean and median of *neighbor* SSFR in bins of r_perp 
        neighSSFR_rperpbins, pure_neighSSFR_rperpbins = [], [] 
        for i_rperp in range(len(rperp_bins)-1): 
            rperp_cut = (neigh_rperp > rperp_bins[i_rperp]) & \
                    (neigh_rperp <= rperp_bins[i_rperp+1]) & \
                    (np.isnan(neigh_ssfr) == False)
            tot_cuts = np.where(rperp_cut) 
            neighSSFR_rperpbins.append(neigh_ssfr[tot_cuts])

            pure_rperp_cut = (pure_neigh_rperp > rperp_bins[i_rperp]) & \
                    (pure_neigh_rperp <= rperp_bins[i_rperp+1]) & \
                    (np.isnan(pure_neigh_ssfr) == False)
            pure_neigh_cut = (pure_neigh_psat <= 0.5) 
            pure_tot_cuts = np.where(pure_rperp_cut & pure_neigh_cut) 
            pure_neighSSFR_rperpbins.append(pure_neigh_ssfr[pure_tot_cuts])

        NeighSSFR_primarybin.append(neighSSFR_rperpbins)
        NeighSSFR_pureprimarybin.append(pure_neighSSFR_rperpbins)

    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure()
    sub = fig.add_subplot(111)
    ssfrplots = [] 
    for i_ssfr, ssfrs_rperp in enumerate(NeighSSFR_primarybin):
        sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), 
                [np.median(ssfrs) for ssfrs in ssfrs_rperp], 
                c=pretty_colors[i_ssfr], lw=1, ls='--') 
    
        pure_ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), 
                [np.median(ssfrs) for ssfrs in NeighSSFR_pureprimarybin[i_ssfr]], 
                c=pretty_colors[i_ssfr], lw=3, ls='-', label=bin_ssfr_label[i_ssfr]) 
        ssfrplots.append(pure_ssfrplot) 

    sub.set_title(r"$10.0 < \mathtt{log}\mathcal{M}^\mathtt{(tot)}_* < 10.5$ ", fontsize=25) 
    # axes
    sub.set_xlim([0, 4]) 
    sub.set_xticks([0, 1, 2, 3, 4]) 
    sub.set_ylim([-12.25, -9.75])
    sub.set_yticks([-12., -11., -10.]) 
    sub.text(0.35, -11.8, r'Pure central primaries', fontsize=20) 
    sub.text(0.35, -12., r'Central neighbors', fontsize=20) 
    #sub.text(0.35, -12., r'ranked in $\mathtt{log}\,\mathtt{SSFR^{(tot)}}$', fontsize=20) 
    sub.minorticks_on() 

    sub.set_xlabel(r'$\mathtt{R_{\perp}}$ [Mpc]', fontsize=25) 
    sub.set_ylabel(r'median log $\mathtt{SSFR^{(fib)}_{neigh}}$ [$\mathtt{yr}^{-1}$]', 
            fontsize=25) 
    first_legend = sub.legend(handles=ssfrplots, loc='upper right', ncol=2, handletextpad=0.1) 
    fig_file = ''.join([UT.dir_fig(), 
        'Conformity.CentralNeighborSSFR_in_PurePrimarySSFR_bins', concat._FileSpec(), '.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 


def Fig_CentralMassbinNeighborSSFR_in_PurePrimarybins(cat_dict, primary_massbin=[10., 10.5]): 
    ''' Figure plotting the median SSFR of *central* neighbors *within the same mass bin* 
    as a function of r_perp for *PURE CENTRAL* primary galaxies binned in SSFR along with 
    the "raw measurement".
    
    It contains *pure central* primary galaxies and *central* neighbor galaxies *within 
    the primary mass bin*.

    Notes
    -----
    * We want an apples to apples comparison, so we use MPA-JHU derived values 
    rather than the VAGC
    '''
    # read conformity catalog based on input catalog dictionary
    concat = clog.ConformCatalog(Mrcut=cat_dict['Mrcut'], 
            primary_delv=cat_dict['primary_delv'], primary_rperp=cat_dict['primary_rperp'],  
            neighbor_delv=cat_dict['neighbor_delv'], neighbor_rperp=cat_dict['neighbor_rperp'])
    catalog = concat.Read() 

    # first sort out all primary galaxies 
    is_primary = np.where(catalog['primary_mpajhu'] == 1)[0]
    ssfr_primary = catalog['ssfr_tot_mpajhu'][is_primary]
    mass_primary = catalog['mass_tot_mpajhu'][is_primary]
    
    # all primaries within mass bin with SSFR
    primary_cut = np.where(
            (mass_primary > primary_massbin[0]) & (mass_primary < primary_massbin[1]) & 
            (np.isnan(ssfr_primary) == False))[0]
    ssfr_primary_cut = ssfr_primary[primary_cut]
    print len(primary_cut), 'primaries within cut'
    # PURE central primaries within mass bin with SSFR 
    pure_primary_cut = np.where(
            (mass_primary > primary_massbin[0]) & (mass_primary < primary_massbin[1]) & 
            (np.isnan(ssfr_primary) == False) & (catalog['p_sat'][is_primary] <= 0.01))[0]
    ssfr_pure_primary_cut = ssfr_primary[pure_primary_cut]
    print len(pure_primary_cut), 'pure central primaries within cut'

    # primary SSFR percentile bins 
    bin_ssfr_list = conform.PrimarySSFR_percentile_bins(ssfr_primary_cut, 
            percentiles=[25, 50, 75, 90]) 
    bin_pure_ssfr_list = conform.PrimarySSFR_percentile_bins(ssfr_pure_primary_cut, 
            percentiles=[25, 50, 75, 90]) 
    bin_ssfr_label = ['0 - 25\%', '25 - 50\%', '50 - 75\%', '> 75\%', '> 90\%']

    rperp_bins = np.arange(0., 4.5, 0.5)
    NeighSSFR_primarybin, NeighSSFR_pureprimarybin = [], []

    for i_bin, bin_ssfr_cut in enumerate(bin_ssfr_list): 
        neigh_inbin = np.concatenate(    # neighbor in primary bin 
                [np.array(catalog['neighbor_indices_mpajhu'][i]) 
                    for i in primary_cut[bin_ssfr_cut]]).astype('int') 
        neigh_rperp = np.concatenate(
                [np.array(catalog['neighbor_rperp_mpajhu'][i]) 
                    for i in primary_cut[bin_ssfr_cut]])
        neigh_ssfr = catalog['ssfr_fib_mpajhu'][neigh_inbin]
        neigh_mass = catalog['mass_tot_mpajhu'][neigh_inbin] #neighbor *total* mass 
        neigh_psat = catalog['p_sat'][neigh_inbin]    # satellite probability of neighbors

        pure_neigh_inbin = np.concatenate(    # neighbor in primary bin 
                [np.array(catalog['neighbor_indices_mpajhu'][i]) 
                    for i in pure_primary_cut[bin_pure_ssfr_list[i_bin]]]).astype('int') 
        pure_neigh_rperp = np.concatenate(
                [np.array(catalog['neighbor_rperp_mpajhu'][i]) 
                    for i in pure_primary_cut[bin_pure_ssfr_list[i_bin]]])
        pure_neigh_ssfr = catalog['ssfr_fib_mpajhu'][pure_neigh_inbin]
        pure_neigh_mass = catalog['mass_tot_mpajhu'][pure_neigh_inbin] #neighbor *total* mass 
        pure_neigh_psat = catalog['p_sat'][pure_neigh_inbin]    # satellite probability of neighbors

        # calculate the mean and median of *neighbor* SSFR in bins of r_perp 
        neighSSFR_rperpbins, pure_neighSSFR_rperpbins = [], [] 
        for i_rperp in range(len(rperp_bins)-1): 
            rperp_cut = (neigh_rperp > rperp_bins[i_rperp]) & \
                    (neigh_rperp <= rperp_bins[i_rperp+1]) & \
                    (np.isnan(neigh_ssfr) == False)
            tot_cuts = np.where(rperp_cut) 
            neighSSFR_rperpbins.append(neigh_ssfr[tot_cuts])

            pure_rperp_cut = (pure_neigh_rperp > rperp_bins[i_rperp]) & \
                    (pure_neigh_rperp <= rperp_bins[i_rperp+1]) & \
                    (np.isnan(pure_neigh_ssfr) == False)
            pure_neigh_cut = (pure_neigh_psat <= 0.5) 
            pure_mass_cut = (pure_neigh_mass > primary_massbin[0]) & (pure_neigh_mass < primary_massbin[1])
            pure_tot_cuts = np.where(pure_rperp_cut & pure_neigh_cut & pure_mass_cut) 
            pure_neighSSFR_rperpbins.append(pure_neigh_ssfr[pure_tot_cuts])

        NeighSSFR_primarybin.append(neighSSFR_rperpbins)
        NeighSSFR_pureprimarybin.append(pure_neighSSFR_rperpbins)

    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure()
    sub = fig.add_subplot(111)
    ssfrplots = [] 
    for i_ssfr, ssfrs_rperp in enumerate(NeighSSFR_primarybin):
        sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), 
                [np.median(ssfrs) for ssfrs in ssfrs_rperp], 
                c=pretty_colors[i_ssfr], lw=1, ls='--') 
    
        pure_ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), 
                [np.median(ssfrs) for ssfrs in NeighSSFR_pureprimarybin[i_ssfr]], 
                c=pretty_colors[i_ssfr], lw=3, ls='-', label=bin_ssfr_label[i_ssfr]) 
        ssfrplots.append(pure_ssfrplot) 

    sub.set_title(r"$10.0 < \mathtt{log}\mathcal{M}^\mathtt{(tot)}_* < 10.5$ ", fontsize=25) 
    # axes
    sub.set_xlim([0, 4]) 
    sub.set_xticks([0, 1, 2, 3, 4]) 
    sub.set_ylim([-12.25, -9.75])
    sub.set_yticks([-12., -11., -10.]) 
    sub.text(0.35, -11.7, r'Pure central primaries', fontsize=20) 
    sub.text(0.35, -11.9, r'Central neighbors within', fontsize=20) 
    sub.text(0.35, -12.15, r'$10.0 < \mathtt{log}\mathcal{M}^\mathtt{(tot)}_* < 10.5$', fontsize=20) 
    #sub.text(0.35, -12., r'ranked in $\mathtt{log}\,\mathtt{SSFR^{(tot)}}$', fontsize=20) 
    sub.minorticks_on() 

    sub.set_xlabel(r'$\mathtt{R_{\perp}}$ [Mpc]', fontsize=25) 
    sub.set_ylabel(r'median log $\mathtt{SSFR^{(fib)}_{neigh}}$ [$\mathtt{yr}^{-1}$]', 
            fontsize=25) 
    first_legend = sub.legend(handles=ssfrplots, loc='upper right', ncol=2, handletextpad=0.1) 
    fig_file = ''.join([UT.dir_fig(), 
        'Conformity.CentralMassbinNeighborSSFR_in_PurePrimarySSFR_bins', concat._FileSpec(), '.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 




if __name__=='__main__': 
    #Fig_NeighborSSFR_in_Primarybins({'name': 'tinker', 'Mrcut':18, 
    #    'primary_delv': 500., 'primary_rperp': 0.5, 
    #    'neighbor_delv': 500., 'neighbor_rperp': 5.},
    #    primary_massbin=[10., 10.5])
    #Fig_NeighborSSFR_in_PurePrimarybins({'name': 'tinker', 'Mrcut':18, 
    #    'primary_delv': 500., 'primary_rperp': 0.5, 
    #    'neighbor_delv': 500., 'neighbor_rperp': 5.},
    #    primary_massbin=[10., 10.5])
    Fig_CentralMassbinNeighborSSFR_in_PurePrimarybins({'name': 'tinker', 'Mrcut':18, 
        'primary_delv': 500., 'primary_rperp': 0.5, 
        'neighbor_delv': 500., 'neighbor_rperp': 5.},
        primary_massbin=[10., 10.5])

    #Fig_CentralNeighborSSFR_in_PurePrimarybins({'name': 'tinker', 'Mrcut':18, 
    #    'primary_delv': 500., 'primary_rperp': 0.5, 
    #    'neighbor_delv': 500., 'neighbor_rperp': 5.},
    #    primary_massbin=[10., 10.5])
