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


def Fig_NeighborSSFR_rperp_PrimaryBins(cat_dict):
    ''' Ultimate plot of the SSFR of neighboring galaxies of primary galaxies 
    binned in SSFRs. We look at 
    - neighbor SSFR(r_perp) for all primaries 
    - neighbor SSFR(r_perp) for pure central primaries.
    - central neighbor SSFR(r_perp) for pure central primaries. 
    - central neighbor SSFR(r_perp) within 10 < logM* < 10.5 for pure central primaries. 

    '''
    rperp_bins=np.arange(0., 4.5, 0.5)
    # read conformity catalog based on input catalog dictionary
    #concat = clog.ConformCatalog(Mrcut=cat_dict['Mrcut'], 
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
    
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(24,6))
    bkgd = fig.add_subplot(111, frameon=False)

    for i_plot in range(4):  
        sub = fig.add_subplot(1,4,i_plot+1)

        if i_plot == 0: 
            primary_groupid = 'all'
            neighbor_groupid = 'all'
            neighbor_massbin = None 
        elif i_plot == 1: 
            primary_groupid = 'pure_centrals'
            neighbor_groupid = 'all'
            neighbor_massbin = None 
        elif i_plot == 2: 
            primary_groupid = 'pure_centrals'
            neighbor_groupid = 'centrals'
            neighbor_massbin = None 
        elif i_plot == 3: 
            primary_groupid = 'pure_centrals'
            neighbor_groupid = 'centrals'
            neighbor_massbin = [10., 10.5] 

        results = conform.NeighborSSFR_rperp_PrimaryBins(catalog, rperp_bins=rperp_bins,
                percentiles=[25, 50, 75, 90], quantiles=None, 
                primary_pipeline='mpajhu', 
                primary_groupid=primary_groupid, 
                primary_massbin=[10., 10.5], 
                neighbor_pipeline='mpajhu',
                neighbor_groupid=neighbor_groupid, 
                neighbor_massbin=neighbor_massbin)
        primary_SSFRbin_limits = results['primary_SSFRbin_limits']
        primary_SSFRbin_label = results['primary_SSFRbin_label']
        neighSSFR_rperp_primarybins = results['neighbor_SSFR_rperp_primary_bins']
        if i_plot == 0: 
            all_all = neighSSFR_rperp_primarybins
    
        # jackknifes 
        jack_results = []
        jack_bins = [5,5]
        for i_jack in range(1, jack_bins[0]*jack_bins[1]+1): 
            jack_results_i = conform.Jackknife_NeighborSSFR_rperp_PrimaryBins(
                    catalog, n_jack=i_jack, RADec_bins=jack_bins, 
                    rperp_bins=rperp_bins, percentiles=None, quantiles=primary_SSFRbin_limits, 
                    primary_pipeline='mpajhu', 
                    primary_groupid=primary_groupid, 
                    primary_massbin=[10., 10.5], 
                    neighbor_pipeline='mpajhu',
                    neighbor_groupid=neighbor_groupid, 
                    neighbor_massbin=neighbor_massbin)
            jack_results.append(jack_results_i['neighbor_SSFR_rperp_primary_bins'])
    
        if i_plot == 0: # plot Kauffmann et al.(2013) 
            for i_ssfr, ssfrs_rperp in enumerate(neighSSFR_rperp_primarybins):
                kauff_str = ['0to25', '25to50', '50to75', '75plus', '90plus']
                kauff_dat_file = ''.join([UT.dir_dat(), 
                    'literature/', 'kauff2013_', kauff_str[i_ssfr], '.dat']) 
                kauff_rperp, kauff_ssfr = np.loadtxt(kauff_dat_file, unpack=True, usecols=[0,1]) 
                sub.scatter(kauff_rperp, kauff_ssfr, c=pretty_colors[i_ssfr], lw=0) 
                kauffplot, = sub.plot(kauff_rperp, kauff_ssfr, 
                        c=pretty_colors[i_ssfr], lw=2, ls=':', label='Kauffmann+(2013)')
        else: 
            for i_ssfr, ssfrs_rperp in enumerate(all_all):
                ssfr_tot = np.array([np.median(ssfrs) for ssfrs in ssfrs_rperp]) 
                ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
                        c=pretty_colors[i_ssfr], lw=2, ls='--') 

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

        label_pipe = 'tot; mpajhu'
        label_ssfr_neigh = 'fib; mpajhu'
        # axes
        sub.set_xlim([0, 4]) 
        sub.set_ylim([-12.25, -9.75])
        sub.set_xticks([0, 1, 2, 3, 4]) 
        sub.set_yticks([-12., -11., -10.]) 
        if i_plot == 1:
            sub.set_xticklabels(['', 1, 2, 3, 4]) 
            sub.set_yticklabels([]) 
        elif i_plot == 2:
            #sub.set_xticklabels([]) 
            sub.set_xticklabels(['', 1, 2, 3, 4]) 
            sub.set_yticklabels([]) 
        elif i_plot == 3: 
            sub.set_yticklabels([]) 
            sub.set_xticklabels(['', 1, 2, 3, 4]) 
        #sub.text(0.1, -11.9, r'Ranked in $\mathtt{log}\mathtt{SSFR^{('+label_pipe+')}}$', fontsize=20) 

        if i_plot == 0: 
            sub.text(2.0, -10.15, 'All Primaries \nAll Neighbors', fontsize=20)
        elif i_plot == 1: 
            sub.text(1.2, -10.15, 'Pure Central Primaries \nAll Neighbors', fontsize=20)
        elif i_plot == 2: 
            sub.text(1.2, -10.15, 'Pure Central Primaries \nCentral Neighbors', fontsize=20)
        elif i_plot == 3: 
            pass
            #sub.text(0.3, -11.8, 'Pure Central Primaries \nCentral Neighbors', fontsize=20)
            #sub.text(0.3, -12., 
            #        str(neighbor_massbin[0])+'$< $log$\mathcal{M}_*^\mathtt{neigh} <$'+str(neighbor_massbin[1]), 
            #        fontsize=20)
        sub.minorticks_on() 

        if i_plot == 1: 
            sub.legend(handles=ssfrplots, loc='lower right', ncol=2, handletextpad=0.1) 
        elif i_plot == 0: 
            sub.legend(handles=[kauffplot], loc='lower right', handletextpad=0.1, markerscale=10)
    
    #fig.title(r"$10.0 < \mathtt{log}\mathcal{M}^\mathtt{("+label_pipe+")}_* < 10.5$ ", fontsize=25) 
    bkgd.set_xticklabels([]) 
    bkgd.set_xlabel(r'$\mathtt{R_{\perp}}$ [Mpc]', labelpad=20, fontsize=30) 
    bkgd.set_yticklabels([]) 
    bkgd.set_ylabel(''.join([
        r'median log($\mathtt{SSFR_{(neigh)}^{(', 
        label_ssfr_neigh, 
        ')}}$ [$\mathtt{yr}^{-1}$])']), 
        labelpad=30, fontsize=30) 
    fig.subplots_adjust(wspace=0., hspace=0.)
    
    fig_file = ''.join([UT.dir_fig(), 
        'neighborSSFR_rprep_primarybins', concat_file_spec, '.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 


def Fig_forTinker_NeighborSSFR_rperp_PrimaryBins():
    ''' Ultimate plot of the SSFR of neighboring galaxies of primary galaxies 
    binned in SSFRs. We look at 
    - neighbor SSFR(r_perp) for all primaries 
    - neighbor SSFR(r_perp) for pure central primaries.
    - central neighbor SSFR(r_perp) for pure central primaries. 
    - central neighbor SSFR(r_perp) within 10 < logM* < 10.5 for pure central primaries. 

    '''
    rperp_bins=np.arange(0., 4.5, 0.5)
    cat_dict = {'name': 'tinkauff', 'primary_delv': 500., 'primary_rperp': 0.5, 
            'neighbor_delv': 500., 'neighbor_rperp': 5.}

    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(24,6))
    bkgd = fig.add_subplot(111, frameon=False)

    for i_plot in range(4):  
        if i_plot != 3: 
            cat_dict['Mass_cut'] = 9.25
        else: 
            cat_dict['Mass_cut'] = 10.0

        # read conformity catalog based on input catalog dictionary
        #concat = clog.ConformCatalog(Mrcut=cat_dict['Mrcut'], 
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
    
        sub = fig.add_subplot(1,4,i_plot+1)

        if i_plot == 0: 
            primary_groupid = 'all'
            neighbor_groupid = 'all'
            neighbor_massbin = None 
        elif i_plot == 1: 
            primary_groupid = 'pure_centrals'
            neighbor_groupid = 'all'
            neighbor_massbin = None 
        elif i_plot == 2: 
            primary_groupid = 'not_pure_centrals'
            neighbor_groupid = 'all'
            neighbor_massbin = None 
        elif i_plot == 3: 
            primary_groupid = 'pure_centrals'
            neighbor_groupid = 'centrals'
            neighbor_massbin = [10., 10.5] 

        results = conform.NeighborSSFR_rperp_PrimaryBins(catalog, rperp_bins=rperp_bins,
                percentiles=[25, 50, 75, 90], quantiles=None, 
                primary_pipeline='mpajhu', 
                primary_groupid=primary_groupid, 
                primary_massbin=[10., 10.5], 
                neighbor_pipeline='mpajhu',
                neighbor_groupid=neighbor_groupid, 
                neighbor_massbin=neighbor_massbin)
        primary_SSFRbin_limits = results['primary_SSFRbin_limits']
        primary_SSFRbin_label = results['primary_SSFRbin_label']
        neighSSFR_rperp_primarybins = results['neighbor_SSFR_rperp_primary_bins']
        if i_plot == 0: 
            all_all = neighSSFR_rperp_primarybins
    
        # jackknifes 
        jack_results = []
        jack_bins = [5,5]
        for i_jack in range(1, jack_bins[0]*jack_bins[1]+1): 
            jack_results_i = conform.Jackknife_NeighborSSFR_rperp_PrimaryBins(
                    catalog, n_jack=i_jack, RADec_bins=jack_bins, 
                    rperp_bins=rperp_bins, percentiles=None, quantiles=primary_SSFRbin_limits, 
                    primary_pipeline='mpajhu', 
                    primary_groupid=primary_groupid, 
                    primary_massbin=[10., 10.5], 
                    neighbor_pipeline='mpajhu',
                    neighbor_groupid=neighbor_groupid, 
                    neighbor_massbin=neighbor_massbin)
            jack_results.append(jack_results_i['neighbor_SSFR_rperp_primary_bins'])
    
        if i_plot == 0: # plot Kauffmann et al.(2013) 
            for i_ssfr, ssfrs_rperp in enumerate(neighSSFR_rperp_primarybins):
                kauff_str = ['0to25', '25to50', '50to75', '75plus', '90plus']
                kauff_dat_file = ''.join([UT.dir_dat(), 
                    'literature/', 'kauff2013_', kauff_str[i_ssfr], '.dat']) 
                kauff_rperp, kauff_ssfr = np.loadtxt(kauff_dat_file, unpack=True, usecols=[0,1]) 
                #sub.scatter(kauff_rperp, kauff_ssfr, c=pretty_colors[i_ssfr], lw=0) 
                if i_ssfr == 0: 
                    kauffplot, = sub.plot(kauff_rperp, kauff_ssfr, 
                            c=pretty_colors[i_ssfr], lw=2, ls=':', label='Kauffmann+(2013)')
                else: 
                    sub.plot(kauff_rperp, kauff_ssfr, c=pretty_colors[i_ssfr], lw=2, ls=':')

        else: 
            for i_ssfr, ssfrs_rperp in enumerate(all_all):
                ssfr_tot = np.array([np.median(ssfrs) for ssfrs in ssfrs_rperp]) 
                ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
                        c=pretty_colors[i_ssfr], lw=2, ls='--') 

        ssfrplots = []
        for i_ssfr, ssfrs_rperp in enumerate(neighSSFR_rperp_primarybins):
            ssfr_tot = np.array([np.median(ssfrs) for ssfrs in ssfrs_rperp]) 
            ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
                    c=pretty_colors[i_ssfr], lw=3, ls='-', label=primary_SSFRbin_label[i_ssfr]) 
            if i_ssfr == 0: 
                leg_ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
                    c=pretty_colors[i_ssfr], lw=3, ls='-', label='Tinker+(2016)')
        
            err_jack = np.zeros(len(ssfr_tot)) 
            for jack_result in jack_results: 
                ssfr_jack = np.array([np.median(ssfrs) for ssfrs in jack_result[i_ssfr]]) 
                err_jack += (ssfr_jack - ssfr_tot)**2
            err_jack *= np.float(jack_bins[0] * jack_bins[1] - 1)/np.float(jack_bins[0] * jack_bins[1])
            sub.errorbar(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, yerr=np.sqrt(err_jack), 
                    c=pretty_colors[i_ssfr], elinewidth=2, capsize=5) 
            ssfrplots.append(ssfrplot) 

        label_pipe = 'tot; mpajhu'
        label_ssfr_neigh = 'fib; mpajhu'
        # axes
        sub.set_xlim([0, 4]) 
        sub.set_ylim([-12.25, -9.75])
        sub.set_xticks([0, 1, 2, 3, 4]) 
        sub.set_yticks([-12., -11., -10.]) 
        if i_plot == 1:
            sub.set_xticklabels(['', 1, 2, 3, 4]) 
            sub.set_yticklabels([]) 
        elif i_plot == 2:
            #sub.set_xticklabels([]) 
            sub.set_xticklabels(['', 1, 2, 3, 4]) 
            sub.set_yticklabels([]) 
        elif i_plot == 3: 
            sub.set_yticklabels([]) 
            sub.set_xticklabels(['', 1, 2, 3, 4]) 
        #sub.text(0.1, -11.9, r'Ranked in $\mathtt{log}\mathtt{SSFR^{('+label_pipe+')}}$', fontsize=20) 

        if i_plot == 0: 
            sub.text(2.0, -10., 'All Primaries', fontsize=20)
            sub.text(2.0, -10.2, 'All Secondaries', fontsize=20)
        elif i_plot == 1: 
            sub.text(1.0, -10., 'Pure Central Primaries', fontsize=20)
            sub.text(1.0, -10.2, 'All Secondaries', fontsize=20)
        elif i_plot == 2: 
            sub.text(0.6, -10., 'Not Pure Central Primaries', fontsize=20)
            sub.text(1., -10.2, 'All Secondaries', fontsize=20)
        elif i_plot == 3: 
            pass
            #sub.text(0.3, -11.8, 'Pure Central Primaries \nCentral Neighbors', fontsize=20)
            #sub.text(0.3, -12., 
            #        str(neighbor_massbin[0])+'$< $log$\mathcal{M}_*^\mathtt{neigh} <$'+str(neighbor_massbin[1]), 
            #        fontsize=20)
        sub.minorticks_on() 

        if i_plot == 1: 
            sub.legend(handles=ssfrplots, loc='lower right', ncol=2, handletextpad=0.1) 
        elif i_plot == 0: 
            kleg = sub.legend(handles=[kauffplot, leg_ssfrplot], loc='lower right', handletextpad=0.1, 
                    prop={'size': 25}, markerscale=50)
            for legobj in kleg.legendHandles:
                legobj.set_linewidth(5.)
    
    #fig.title(r"$10.0 < \mathtt{log}\mathcal{M}^\mathtt{("+label_pipe+")}_* < 10.5$ ", fontsize=25) 
    bkgd.set_xticklabels([]) 
    bkgd.set_xlabel(r'$\mathtt{R_{\perp}}$ [Mpc]', labelpad=20, fontsize=30) 
    bkgd.set_yticklabels([]) 
    bkgd.set_ylabel(''.join([
        r'log($\mathtt{SSFR^{(fib)}}$ [$\mathtt{yr}^{-1}$])']), 
        labelpad=30, fontsize=30) 
    fig.subplots_adjust(wspace=0., hspace=0.)
    
    fig_file = ''.join([UT.dir_fig(), 
        'forTinker_neighborSSFR_rprep_primarybins', concat_file_spec, '.pdf']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 


def Fig_TinKauff_Iso_NeighborSSFR_rperp_PrimaryBins():
    ''' Ultimate plot of the SSFR of neighboring galaxies of primary galaxies 
    binned in SSFRs. We look at 
    - neighbor SSFR(r_perp) for all primaries 
    - neighbor SSFR(r_perp) for pure central primaries.
    - central neighbor SSFR(r_perp) for pure central primaries. 
    - central neighbor SSFR(r_perp) within 10 < logM* < 10.5 for pure central primaries. 

    '''
    rperp_bins=np.arange(0., 4.5, 0.5)
    cat_dict = {'name': 'tinkauff_iso', 'primary_delv': 500., 'primary_rperp': 0.5, 
            'neighbor_delv': 500., 'neighbor_rperp': 5.}

    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(18,6))
    bkgd = fig.add_subplot(111, frameon=False)

    for i_plot in range(3):  
        cat_dict['Mass_cut'] = 9.25

        # read conformity catalog based on input catalog dictionary
        #concat = clog.ConformCatalog(Mrcut=cat_dict['Mrcut'], 
        if cat_dict['name'] == 'tinker': 
            catalog_prop = {'Mrcut': cat_dict['Mrcut']} 
        elif 'tinkauff' in cat_dict['name']: 
            catalog_prop = {'Mass_cut': cat_dict['Mass_cut']} 
        elif cat_dict['name'] == 'kauff': 
            catalog_prop = {} 
        concat = clog.ConformCatalog(cat_dict['name'], catalog_prop=catalog_prop, 
                primary_delv=cat_dict['primary_delv'], primary_rperp=cat_dict['primary_rperp'],  
                neighbor_delv=cat_dict['neighbor_delv'], neighbor_rperp=cat_dict['neighbor_rperp'])
        catalog = concat.Read() 
        concat_file_spec = concat._FileSpec()    # conformity catalog specification string 
    
        sub = fig.add_subplot(1,3,i_plot+1)

        if i_plot == 0: 
            primary_groupid = 'all'
            neighbor_groupid = 'all'
            neighbor_massbin = None 
        elif i_plot == 1: 
            primary_groupid = 'pure_centrals_iso'
            neighbor_groupid = 'all'
            neighbor_massbin = None 
        elif i_plot == 2: 
            primary_groupid = 'not_pure_centrals_iso'
            neighbor_groupid = 'all'
            neighbor_massbin = None 

        results = conform.NeighborSSFR_rperp_PrimaryBins(catalog, rperp_bins=rperp_bins,
                percentiles=[25, 50, 75, 90], quantiles=None, 
                primary_pipeline='mpajhu', 
                primary_groupid=primary_groupid, 
                primary_massbin=[10., 10.5], 
                neighbor_pipeline='mpajhu',
                neighbor_groupid=neighbor_groupid, 
                neighbor_massbin=neighbor_massbin)
        primary_SSFRbin_limits = results['primary_SSFRbin_limits']
        primary_SSFRbin_label = results['primary_SSFRbin_label']
        neighSSFR_rperp_primarybins = results['neighbor_SSFR_rperp_primary_bins']
        if i_plot == 0: 
            all_all = neighSSFR_rperp_primarybins
    
        # jackknifes 
        jack_results = []
        jack_bins = [5,5]
        for i_jack in range(1, jack_bins[0]*jack_bins[1]+1): 
            jack_results_i = conform.Jackknife_NeighborSSFR_rperp_PrimaryBins(
                    catalog, n_jack=i_jack, RADec_bins=jack_bins, 
                    rperp_bins=rperp_bins, percentiles=None, quantiles=primary_SSFRbin_limits, 
                    primary_pipeline='mpajhu', 
                    primary_groupid=primary_groupid, 
                    primary_massbin=[10., 10.5], 
                    neighbor_pipeline='mpajhu',
                    neighbor_groupid=neighbor_groupid, 
                    neighbor_massbin=neighbor_massbin)
            jack_results.append(jack_results_i['neighbor_SSFR_rperp_primary_bins'])
    
        if i_plot == 0: # plot Kauffmann et al.(2013) 
            for i_ssfr, ssfrs_rperp in enumerate(neighSSFR_rperp_primarybins):
                kauff_str = ['0to25', '25to50', '50to75', '75plus', '90plus']
                kauff_dat_file = ''.join([UT.dir_dat(), 
                    'literature/', 'kauff2013_', kauff_str[i_ssfr], '.dat']) 
                kauff_rperp, kauff_ssfr = np.loadtxt(kauff_dat_file, unpack=True, usecols=[0,1]) 
                #sub.scatter(kauff_rperp, kauff_ssfr, c=pretty_colors[i_ssfr], lw=0) 
                if i_ssfr == 0: 
                    kauffplot, = sub.plot(kauff_rperp, kauff_ssfr, 
                            c=pretty_colors[i_ssfr], lw=2, ls=':', label='Kauffmann+(2013)')
                else: 
                    sub.plot(kauff_rperp, kauff_ssfr, c=pretty_colors[i_ssfr], lw=2, ls=':')

        else: 
            for i_ssfr, ssfrs_rperp in enumerate(all_all):
                ssfr_tot = np.array([np.median(ssfrs) for ssfrs in ssfrs_rperp]) 
                ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
                        c=pretty_colors[i_ssfr], lw=2, ls='--') 

        ssfrplots = []
        for i_ssfr, ssfrs_rperp in enumerate(neighSSFR_rperp_primarybins):
            ssfr_tot = np.array([np.median(ssfrs) for ssfrs in ssfrs_rperp]) 
            ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
                    c=pretty_colors[i_ssfr], lw=3, ls='-', label=primary_SSFRbin_label[i_ssfr]) 
            if i_ssfr == 0: 
                leg_ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
                    c=pretty_colors[i_ssfr], lw=3, ls='-', label='Tinker+(2016)')
        
            err_jack = np.zeros(len(ssfr_tot)) 
            for jack_result in jack_results: 
                ssfr_jack = np.array([np.median(ssfrs) for ssfrs in jack_result[i_ssfr]]) 
                err_jack += (ssfr_jack - ssfr_tot)**2
            err_jack *= np.float(jack_bins[0] * jack_bins[1] - 1)/np.float(jack_bins[0] * jack_bins[1])
            sub.errorbar(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, yerr=np.sqrt(err_jack), 
                    c=pretty_colors[i_ssfr], elinewidth=2, capsize=5) 
            ssfrplots.append(ssfrplot) 

        label_pipe = 'tot; mpajhu'
        label_ssfr_neigh = 'fib; mpajhu'
        # axes
        sub.set_xlim([0, 4]) 
        sub.set_ylim([-12.25, -9.75])
        sub.set_xticks([0, 1, 2, 3, 4]) 
        sub.set_yticks([-12., -11., -10.]) 
        if i_plot == 1:
            sub.set_xticklabels(['', 1, 2, 3, 4]) 
            sub.set_yticklabels([]) 
        elif i_plot == 2:
            #sub.set_xticklabels([]) 
            sub.set_xticklabels(['', 1, 2, 3, 4]) 
            sub.set_yticklabels([]) 
        #sub.text(0.1, -11.9, r'Ranked in $\mathtt{log}\mathtt{SSFR^{('+label_pipe+')}}$', fontsize=20) 

        if i_plot == 0: 
            sub.text(2.0, -10., 'All Primaries', fontsize=20)
            sub.text(2.0, -10.2, 'All Secondaries', fontsize=20)
        elif i_plot == 1: 
            sub.text(1.0, -10., 'Pure Central Primaries', fontsize=20)
            sub.text(1.0, -10.2, 'All Secondaries', fontsize=20)
        elif i_plot == 2: 
            sub.text(0.6, -10., 'Not Pure Central Primaries', fontsize=20)
            sub.text(1., -10.2, 'All Secondaries', fontsize=20)
            #sub.text(0.3, -11.8, 'Pure Central Primaries \nCentral Neighbors', fontsize=20)
            #sub.text(0.3, -12., 
            #        str(neighbor_massbin[0])+'$< $log$\mathcal{M}_*^\mathtt{neigh} <$'+str(neighbor_massbin[1]), 
            #        fontsize=20)
        sub.minorticks_on() 

        if i_plot == 1: 
            sub.legend(handles=ssfrplots, loc='lower right', ncol=2, handletextpad=0.1) 
        elif i_plot == 0: 
            kleg = sub.legend(handles=[kauffplot, leg_ssfrplot], loc='lower right', handletextpad=0.1, 
                    prop={'size': 25}, markerscale=50)
            for legobj in kleg.legendHandles:
                legobj.set_linewidth(5.)
    
    #fig.title(r"$10.0 < \mathtt{log}\mathcal{M}^\mathtt{("+label_pipe+")}_* < 10.5$ ", fontsize=25) 
    bkgd.set_xticklabels([]) 
    bkgd.set_xlabel(r'$\mathtt{R_{\perp}}$ [Mpc]', labelpad=20, fontsize=30) 
    bkgd.set_yticklabels([]) 
    bkgd.set_ylabel(''.join([
        r'log($\mathtt{SSFR^{(fib)}}$ [$\mathtt{yr}^{-1}$])']), 
        labelpad=30, fontsize=30) 
    fig.subplots_adjust(wspace=0., hspace=0.)
    
    fig_file = ''.join([UT.dir_fig(), 
        'TinKauff_Iso_neighborSSFR_rprep_primarybins', concat_file_spec, '.pdf']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 


def Fig_TinKauff_NeighborSSFR_rperp_PrimaryBins_IsoComparison():
    ''' plot of the SSFR of neighboring galaxies of primary galaxies 
    binned in SSFRs. We look at 
    - neighbor SSFR(r_perp) for all primaries 
    - neighbor SSFR(r_perp) for group primaries.
    - neighbor SSFR(r_perp) for isolated primaries.

    '''
    rperp_bins=np.arange(0., 4.5, 0.5)
    cat_dict = {'name': 'tinkauff_iso', 'primary_delv': 500., 'primary_rperp': 0.5, 
            'neighbor_delv': 500., 'neighbor_rperp': 5.}

    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(18,6))
    bkgd = fig.add_subplot(111, frameon=False)

    for i_plot in range(3):  
        cat_dict['Mass_cut'] = 9.25

        # read conformity catalog based on input catalog dictionary
        #concat = clog.ConformCatalog(Mrcut=cat_dict['Mrcut'], 
        if cat_dict['name'] == 'tinker': 
            catalog_prop = {'Mrcut': cat_dict['Mrcut']} 
        elif 'tinkauff' in cat_dict['name']: 
            catalog_prop = {'Mass_cut': cat_dict['Mass_cut']} 
        elif cat_dict['name'] == 'kauff': 
            catalog_prop = {} 
        concat = clog.ConformCatalog(cat_dict['name'], catalog_prop=catalog_prop, 
                primary_delv=cat_dict['primary_delv'], primary_rperp=cat_dict['primary_rperp'],  
                neighbor_delv=cat_dict['neighbor_delv'], neighbor_rperp=cat_dict['neighbor_rperp'])
        catalog = concat.Read() 
        concat_file_spec = concat._FileSpec()    # conformity catalog specification string 
    
        sub = fig.add_subplot(1,3,i_plot+1)

        if i_plot == 0: 
            primary_groupid = 'all'
            neighbor_groupid = 'all'
            neighbor_massbin = None 
        elif i_plot == 1: 
            primary_groupid = 'pure_centrals'
            neighbor_groupid = 'all'
            neighbor_massbin = None 
        elif i_plot == 2: 
            primary_groupid = 'pure_centrals_iso'
            neighbor_groupid = 'all'
            neighbor_massbin = None 

        results = conform.NeighborSSFR_rperp_PrimaryBins(catalog, rperp_bins=rperp_bins,
                percentiles=[25, 50, 75, 90], quantiles=None, 
                primary_pipeline='mpajhu', 
                primary_groupid=primary_groupid, 
                primary_massbin=[10., 10.5], 
                neighbor_pipeline='mpajhu',
                neighbor_groupid=neighbor_groupid, 
                neighbor_massbin=neighbor_massbin)
        primary_SSFRbin_limits = results['primary_SSFRbin_limits']
        primary_SSFRbin_label = results['primary_SSFRbin_label']
        neighSSFR_rperp_primarybins = results['neighbor_SSFR_rperp_primary_bins']
        if i_plot == 0: 
            all_all = neighSSFR_rperp_primarybins
    
        # jackknifes 
        jack_results = []
        jack_bins = [5,5]
        for i_jack in range(1, jack_bins[0]*jack_bins[1]+1): 
            jack_results_i = conform.Jackknife_NeighborSSFR_rperp_PrimaryBins(
                    catalog, n_jack=i_jack, RADec_bins=jack_bins, 
                    rperp_bins=rperp_bins, percentiles=None, quantiles=primary_SSFRbin_limits, 
                    primary_pipeline='mpajhu', 
                    primary_groupid=primary_groupid, 
                    primary_massbin=[10., 10.5], 
                    neighbor_pipeline='mpajhu',
                    neighbor_groupid=neighbor_groupid, 
                    neighbor_massbin=neighbor_massbin)
            jack_results.append(jack_results_i['neighbor_SSFR_rperp_primary_bins'])
    
        if i_plot == 0: # plot Kauffmann et al.(2013) 
            for i_ssfr, ssfrs_rperp in enumerate(neighSSFR_rperp_primarybins):
                kauff_str = ['0to25', '25to50', '50to75', '75plus', '90plus']
                kauff_dat_file = ''.join([UT.dir_dat(), 
                    'literature/', 'kauff2013_', kauff_str[i_ssfr], '.dat']) 
                kauff_rperp, kauff_ssfr = np.loadtxt(kauff_dat_file, unpack=True, usecols=[0,1]) 
                #sub.scatter(kauff_rperp, kauff_ssfr, c=pretty_colors[i_ssfr], lw=0) 
                if i_ssfr == 0: 
                    kauffplot, = sub.plot(kauff_rperp, kauff_ssfr, 
                            c=pretty_colors[i_ssfr], lw=2, ls=':', label='Kauffmann+(2013)')
                else: 
                    sub.plot(kauff_rperp, kauff_ssfr, c=pretty_colors[i_ssfr], lw=2, ls=':')

        else: 
            for i_ssfr, ssfrs_rperp in enumerate(all_all):
                ssfr_tot = np.array([np.median(ssfrs) for ssfrs in ssfrs_rperp]) 
                ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
                        c=pretty_colors[i_ssfr], lw=2, ls='--') 

        ssfrplots = []
        for i_ssfr, ssfrs_rperp in enumerate(neighSSFR_rperp_primarybins):
            ssfr_tot = np.array([np.median(ssfrs) for ssfrs in ssfrs_rperp]) 
            ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
                    c=pretty_colors[i_ssfr], lw=3, ls='-', label=primary_SSFRbin_label[i_ssfr]) 
            if i_ssfr == 0: 
                leg_ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
                    c=pretty_colors[i_ssfr], lw=3, ls='-', label='Tinker+(2016)')
        
            err_jack = np.zeros(len(ssfr_tot)) 
            for jack_result in jack_results: 
                ssfr_jack = np.array([np.median(ssfrs) for ssfrs in jack_result[i_ssfr]]) 
                err_jack += (ssfr_jack - ssfr_tot)**2
            err_jack *= np.float(jack_bins[0] * jack_bins[1] - 1)/np.float(jack_bins[0] * jack_bins[1])
            sub.errorbar(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, yerr=np.sqrt(err_jack), 
                    c=pretty_colors[i_ssfr], elinewidth=2, capsize=5) 
            ssfrplots.append(ssfrplot) 

        label_pipe = 'tot; mpajhu'
        label_ssfr_neigh = 'fib; mpajhu'
        # axes
        sub.set_xlim([0, 4]) 
        sub.set_ylim([-12.25, -9.75])
        sub.set_xticks([0, 1, 2, 3, 4]) 
        sub.set_yticks([-12., -11., -10.]) 
        if i_plot == 1:
            sub.set_xticklabels(['', 1, 2, 3, 4]) 
            sub.set_yticklabels([]) 
        elif i_plot == 2:
            #sub.set_xticklabels([]) 
            sub.set_xticklabels(['', 1, 2, 3, 4]) 
            sub.set_yticklabels([]) 
        #sub.text(0.1, -11.9, r'Ranked in $\mathtt{log}\mathtt{SSFR^{('+label_pipe+')}}$', fontsize=20) 

        if i_plot == 0: 
            sub.text(2.0, -10., 'All Primaries', fontsize=20)
            sub.text(2.0, -10.2, 'All Secondaries', fontsize=20)
        elif i_plot == 1: 
            sub.text(1.0, -10., 'Group Primaries', fontsize=20)
            sub.text(1.0, -10.2, 'All Secondaries', fontsize=20)
        elif i_plot == 2: 
            sub.text(0.6, -10., 'Isolated Primaries', fontsize=20)
            sub.text(1., -10.2, 'All Secondaries', fontsize=20)
        sub.minorticks_on() 

        if i_plot == 1: 
            sub.legend(handles=ssfrplots, loc='lower right', ncol=2, handletextpad=0.1) 
        elif i_plot == 0: 
            kleg = sub.legend(handles=[kauffplot, leg_ssfrplot], loc='lower right', handletextpad=0.1, 
                    prop={'size': 25}, markerscale=50)
            for legobj in kleg.legendHandles:
                legobj.set_linewidth(5.)
    
    #fig.title(r"$10.0 < \mathtt{log}\mathcal{M}^\mathtt{("+label_pipe+")}_* < 10.5$ ", fontsize=25) 
    bkgd.set_xticklabels([]) 
    bkgd.set_xlabel(r'$\mathtt{R_{\perp}}$ [Mpc]', labelpad=20, fontsize=30) 
    bkgd.set_yticklabels([]) 
    bkgd.set_ylabel(''.join([
        r'log($\mathtt{SSFR^{(fib)}}$ [$\mathtt{yr}^{-1}$])']), 
        labelpad=30, fontsize=30) 
    fig.subplots_adjust(wspace=0., hspace=0.)
    
    fig_file = ''.join([UT.dir_fig(), 
        'TinKauff_neighborSSFR_rprep_primarybins', concat_file_spec, '_IsoComparison.pdf']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 


def Fig_forProposal():
    ''' Ultimate plot of the SSFR of neighboring galaxies of primary galaxies 
    binned in SSFRs. We look at 
    - neighbor SSFR(r_perp) for all primaries 
    - neighbor SSFR(r_perp) for pure central primaries.
    - central neighbor SSFR(r_perp) for pure central primaries. 
    - central neighbor SSFR(r_perp) within 10 < logM* < 10.5 for pure central primaries. 

    '''
    rperp_bins=np.arange(0., 4.5, 0.5)
    cat_dict = {'name': 'tinkauff', 'primary_delv': 500., 'primary_rperp': 0.5, 
            'neighbor_delv': 500., 'neighbor_rperp': 5.}

    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(12,4.5))
    bkgd = fig.add_subplot(111, frameon=False)

    for i_plot in range(2):  
        cat_dict['Mass_cut'] = 9.25
        
        # read conformity catalog based on input catalog dictionary
        #concat = clog.ConformCatalog(Mrcut=cat_dict['Mrcut'], 
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
    
        sub = fig.add_subplot(1,2,i_plot+1)

        if i_plot == 0: 
            primary_groupid = 'all'
            neighbor_groupid = 'all'
            neighbor_massbin = None 
        elif i_plot == 1: 
            primary_groupid = 'pure_centrals'
            neighbor_groupid = 'all'
            neighbor_massbin = None 
        elif i_plot == 2: 
            primary_groupid = 'pure_centrals'
            neighbor_groupid = 'centrals'
            neighbor_massbin = None 
        elif i_plot == 3: 
            primary_groupid = 'pure_centrals'
            neighbor_groupid = 'centrals'
            neighbor_massbin = [10., 10.5] 

        results = conform.NeighborSSFR_rperp_PrimaryBins(catalog, rperp_bins=rperp_bins,
                percentiles=[25, 50, 75, 90], quantiles=None, 
                primary_pipeline='mpajhu', 
                primary_groupid=primary_groupid, 
                primary_massbin=[10., 10.5], 
                neighbor_pipeline='mpajhu',
                neighbor_groupid=neighbor_groupid, 
                neighbor_massbin=neighbor_massbin)
        primary_SSFRbin_limits = results['primary_SSFRbin_limits']
        primary_SSFRbin_label = results['primary_SSFRbin_label']
        neighSSFR_rperp_primarybins = results['neighbor_SSFR_rperp_primary_bins']
        if i_plot == 0: 
            all_all = neighSSFR_rperp_primarybins
    
        # jackknifes 
        jack_results = []
        jack_bins = [5,5]
        for i_jack in range(1, jack_bins[0]*jack_bins[1]+1): 
            jack_results_i = conform.Jackknife_NeighborSSFR_rperp_PrimaryBins(
                    catalog, n_jack=i_jack, RADec_bins=jack_bins, 
                    rperp_bins=rperp_bins, percentiles=None, quantiles=primary_SSFRbin_limits, 
                    primary_pipeline='mpajhu', 
                    primary_groupid=primary_groupid, 
                    primary_massbin=[10., 10.5], 
                    neighbor_pipeline='mpajhu',
                    neighbor_groupid=neighbor_groupid, 
                    neighbor_massbin=neighbor_massbin)
            jack_results.append(jack_results_i['neighbor_SSFR_rperp_primary_bins'])
    
        if i_plot == 0: # plot Kauffmann et al.(2013) 
            for i_ssfr, ssfrs_rperp in enumerate(neighSSFR_rperp_primarybins):
                kauff_str = ['0to25', '25to50', '50to75', '75plus', '90plus']
                kauff_dat_file = ''.join([UT.dir_dat(), 
                    'literature/', 'kauff2013_', kauff_str[i_ssfr], '.dat']) 
                kauff_rperp, kauff_ssfr = np.loadtxt(kauff_dat_file, unpack=True, usecols=[0,1]) 
                #sub.scatter(kauff_rperp, kauff_ssfr, c=pretty_colors[i_ssfr], lw=0) 
                if i_ssfr in [0, 3]:
                    if i_ssfr == 0: 
                        kauffplot, = sub.plot(kauff_rperp, kauff_ssfr, 
                                c=pretty_colors[i_ssfr+1], lw=3, ls='-', label='Kauffmann+(2013)')
                    elif i_ssfr == 3: 
                        kauffplot, = sub.plot(kauff_rperp, kauff_ssfr, 
                                c=pretty_colors[i_ssfr], lw=3, ls='-', label='Kauffmann+(2013)')
        #else: 
        #    for i_ssfr, ssfrs_rperp in enumerate(all_all):
        #        ssfr_tot = np.array([np.median(ssfrs) for ssfrs in ssfrs_rperp]) 
        #        ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
        #                c=pretty_colors[i_ssfr], lw=2, ls='--') 
        
        if i_plot != 0: 
            ssfrplots = []
            for i_ssfr, ssfrs_rperp in enumerate(neighSSFR_rperp_primarybins):
                ssfr_tot = np.array([np.median(ssfrs) for ssfrs in ssfrs_rperp]) 
                if i_ssfr in [0, 3]:
                    err_jack = np.zeros(len(ssfr_tot)) 
                    for jack_result in jack_results: 
                        ssfr_jack = np.array([np.median(ssfrs) for ssfrs in jack_result[i_ssfr]]) 
                        err_jack += (ssfr_jack - ssfr_tot)**2
                    err_jack *= np.float(jack_bins[0] * jack_bins[1] - 1)/np.float(jack_bins[0] * jack_bins[1])

                    if i_ssfr == 0: 
                        ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
                                c=pretty_colors[i_ssfr+1], lw=3, ls='-', label='Low sSFR Centrals') 
                        sub.errorbar(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, yerr=np.sqrt(err_jack), 
                                c=pretty_colors[i_ssfr+1], elinewidth=2, capsize=5) 
                    elif i_ssfr == 3: 
                        ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
                                c=pretty_colors[i_ssfr], lw=3, ls='-', label='High sSFR Centrals') 
                        sub.errorbar(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, yerr=np.sqrt(err_jack), 
                                c=pretty_colors[i_ssfr], elinewidth=2, capsize=5) 
                    ssfrplots.append(ssfrplot) 
        else: 
            for i_ssfr, ssfrs_rperp in enumerate(neighSSFR_rperp_primarybins):
                kauff_dat_file = ''.join([UT.dir_dat(), 
                    'literature/', 'kauff2013_', kauff_str[i_ssfr], '.dat']) 
                kauff_rperp, kauff_ssfr = np.loadtxt(kauff_dat_file, unpack=True, usecols=[0,1]) 

                ssfr_tot = np.array([np.median(ssfrs) for ssfrs in ssfrs_rperp]) 
                if i_ssfr in [0, 3]:
                    err_jack = np.zeros(len(ssfr_tot)) 
                    for jack_result in jack_results: 
                        ssfr_jack = np.array([np.median(ssfrs) for ssfrs in jack_result[i_ssfr]]) 
                        err_jack += (ssfr_jack - ssfr_tot)**2
                    err_jack *= np.float(jack_bins[0] * jack_bins[1] - 1)/np.float(jack_bins[0] * jack_bins[1])

                    if i_ssfr == 0: 
                        print kauff_rperp-0.5*(rperp_bins[:-1]+rperp_bins[1:])
                        sub.errorbar(0.5*(rperp_bins[:-1]+rperp_bins[1:]), kauff_ssfr, yerr=np.sqrt(err_jack), 
                                c=pretty_colors[i_ssfr+1], elinewidth=2, capsize=5) 
                    elif i_ssfr == 3: 
                        sub.errorbar(0.5*(rperp_bins[:-1]+rperp_bins[1:]), kauff_ssfr, yerr=np.sqrt(err_jack), 
                                c=pretty_colors[i_ssfr], elinewidth=2, capsize=5) 

        label_pipe = 'tot; mpajhu'
        label_ssfr_neigh = 'fib; mpajhu'
        # axes
        sub.set_xlim([0, 4]) 
        sub.set_ylim([-12., -9.75])
        sub.set_xticks([0, 1, 2, 3, 4]) 
        sub.set_yticks([-12., -11., -10.]) 
        if i_plot == 1:
            #sub.set_xticklabels(['', 1, 2, 3, 4]) 
            sub.set_yticklabels([]) 

        #if i_plot == 0: 
        #    sub.text(1.0, -10.15, 'Kauffmann+(2013)', fontsize=30)
        #elif i_plot == 1: 
        #    sub.text(1.2, -10.15, 'Tinker+(in prep)', fontsize=30)
        sub.minorticks_on() 

        if i_plot == 1: 
            #sub.legend(handles=ssfrplots, loc='lower right', ncol=2, handletextpad=0.1, prop={'size':20})
            sub.legend(handles=ssfrplots, loc='lower right', handletextpad=0.1, prop={'size':20})
        #elif i_plot == 0: 
        #    kleg = sub.legend(handles=[kauffplot], loc='lower right', handletextpad=0.1, 
        #            prop={'size': 25}, markerscale=50)
        #    for legobj in kleg.legendHandles:
        #        legobj.set_linewidth(5.)
    
    #fig.title(r"$10.0 < \mathtt{log}\mathcal{M}^\mathtt{("+label_pipe+")}_* < 10.5$ ", fontsize=25) 
    bkgd.set_xticklabels([]) 
    bkgd.set_xlabel(r'$\mathtt{R_{\perp}}$ [Mpc]', labelpad=20, fontsize=30) 
    bkgd.set_yticklabels([]) 
    bkgd.set_ylabel(''.join([
        r'log($\mathtt{SSFR}$ [$\mathtt{yr}^{-1}$])']), 
        labelpad=30, fontsize=30) 
    fig.subplots_adjust(wspace=0.1, hspace=0.0)
    
    fig_file = ''.join([UT.dir_fig(), 
        'proposal', concat_file_spec, '.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 


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
    blah, bin_ssfr_list = conform.SSFR_percentilebins(ssfr_primary_cut, 
            percentiles=[25, 50, 75, 90]) 
    bin_ssfr_label = ['0 - 25\%', '25 - 50\%', '50 - 75\%', '> 75\%', '> 90\%']

    NeighSSFR_primarybin = [] 
    rperp_bins = np.arange(0., 4.5, 0.5)
    
    for bin_ssfr_cut in bin_ssfr_list: 
        neigh_inbin = np.concatenate(    # neighbor in primary bin 
                [np.array(catalog['neighbor_indices_mpajhu'][i]) 
                    for i in primary_cut[bin_ssfr_cut]]).astype('int') 
        print neigh_inbin, len(neigh_inbin) 
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
            print 'rp=', 0.5 * (rperp_bins[i_rperp] + rperp_bins[i_rperp+1])
            print np.median(neigh_ssfr[tot_cuts]), len(neigh_ssfr[tot_cuts]) 
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
    mass_primary_cut = mass_primary[primary_cut]
    print len(primary_cut), 'primaries within cut'
    print 'mean M* of pure central primaries ', np.log10(np.mean(10.**mass_primary_cut))
    print 'median M* of pure central primaries ', np.median(mass_primary_cut)
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
        print 'SSFR percentile - ', bin_ssfr_label[i_bin]
        print 'mean = ', np.log10(np.mean(10**mass_primary_cut[bin_ssfr_cut]))
        print 'median = ', np.median(mass_primary_cut[bin_ssfr_cut])

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


def Fig_forTinkerReview_NeighborSSFR_rperp_PrimaryBins():
    ''' plot of the SSFR of neighboring galaxies of primary galaxies 
    binned in SSFRs. 
    - central neighbor SSFR(r_perp) within 10 < logM* < 10.5 for pure central primaries. 

    '''
    rperp_bins=np.arange(0., 4.5, 0.5)
    cat_dict = {'name': 'tinkauff', 'primary_delv': 500., 'primary_rperp': 0.5, 
            'neighbor_delv': 500., 'neighbor_rperp': 5.}

    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(6,6))

    cat_dict['Mass_cut'] = 9.25

    # read conformity catalog based on input catalog dictionary
    #concat = clog.ConformCatalog(Mrcut=cat_dict['Mrcut'], 
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
    
    sub = fig.add_subplot(1,1,1)
    primary_groupid = 'all'
    neighbor_groupid = 'all'
    neighbor_massbin = None 

    results = conform.NeighborSSFR_rperp_PrimaryBins(catalog, rperp_bins=rperp_bins,
            percentiles=[25, 50, 75, 90], quantiles=None, 
            primary_pipeline='mpajhu', 
            primary_groupid=primary_groupid, 
            primary_massbin=[10., 10.5], 
            neighbor_pipeline='mpajhu',
            neighbor_groupid=neighbor_groupid, 
            neighbor_massbin=neighbor_massbin)
    primary_SSFRbin_limits = results['primary_SSFRbin_limits']
    primary_SSFRbin_label = results['primary_SSFRbin_label']
    neighSSFR_rperp_primarybins = results['neighbor_SSFR_rperp_primary_bins']
    all_all = neighSSFR_rperp_primarybins
    
    cat_dict['Mass_cut'] = 10.0

    if cat_dict['name'] == 'tinker': 
        catalog_prop = {'Mrcut': cat_dict['Mrcut']} 
    elif cat_dict['name'] == 'tinkauff': 
        catalog_prop = {'Mass_cut': cat_dict['Mass_cut']} 
    elif cat_dict['name'] == 'kauff': 
        catalog_prop = {} 
    # read conformity catalog based on input catalog dictionary
    #concat = clog.ConformCatalog(Mrcut=cat_dict['Mrcut'], 
    concat = clog.ConformCatalog(cat_dict['name'], catalog_prop=catalog_prop, 
            primary_delv=cat_dict['primary_delv'], primary_rperp=cat_dict['primary_rperp'],  
            neighbor_delv=cat_dict['neighbor_delv'], neighbor_rperp=cat_dict['neighbor_rperp'])
    catalog = concat.Read() 
    concat_file_spec = concat._FileSpec()    # conformity catalog specification string 

    primary_groupid = 'pure_centrals'
    neighbor_groupid = 'centrals'
    neighbor_massbin = [10., 10.5] 

    results = conform.NeighborSSFR_rperp_PrimaryBins(catalog, rperp_bins=rperp_bins,
            percentiles=[25, 50, 75, 90], quantiles=None, 
            primary_pipeline='mpajhu', 
            primary_groupid=primary_groupid, 
            primary_massbin=[10., 10.5], 
            neighbor_pipeline='mpajhu',
            neighbor_groupid=neighbor_groupid, 
            neighbor_massbin=neighbor_massbin)
    primary_SSFRbin_limits = results['primary_SSFRbin_limits']
    primary_SSFRbin_label = results['primary_SSFRbin_label']
    neighSSFR_rperp_primarybins = results['neighbor_SSFR_rperp_primary_bins']

    # jackknifes 
    jack_results = []
    jack_bins = [5,5]
    for i_jack in range(1, jack_bins[0]*jack_bins[1]+1): 
        jack_results_i = conform.Jackknife_NeighborSSFR_rperp_PrimaryBins(
                catalog, n_jack=i_jack, RADec_bins=jack_bins, 
                rperp_bins=rperp_bins, percentiles=None, quantiles=primary_SSFRbin_limits, 
                primary_pipeline='mpajhu', 
                primary_groupid=primary_groupid, 
                primary_massbin=[10., 10.5], 
                neighbor_pipeline='mpajhu',
                neighbor_groupid=neighbor_groupid, 
                neighbor_massbin=neighbor_massbin)
        jack_results.append(jack_results_i['neighbor_SSFR_rperp_primary_bins'])
    
    for i_ssfr, ssfrs_rperp in enumerate(all_all):
        ssfr_tot = np.array([np.median(ssfrs) for ssfrs in ssfrs_rperp]) 
        print ssfr_tot
        if i_ssfr == 0: 
            kauffplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
                    c=pretty_colors[i_ssfr], lw=2, ls='--', label='Kauffmann+(2013)') 
        else: 
            ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
                    c=pretty_colors[i_ssfr], lw=2, ls='--') 

    ssfrplots = []
    for i_ssfr, ssfrs_rperp in enumerate(neighSSFR_rperp_primarybins):
        ssfr_tot = np.array([np.median(ssfrs) for ssfrs in ssfrs_rperp]) 
        ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
                c=pretty_colors[i_ssfr], lw=3, ls='-', label=primary_SSFRbin_label[i_ssfr]) 
        if i_ssfr == 0: 
            leg_ssfrplot, = sub.plot(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, 
                c=pretty_colors[i_ssfr], lw=3, ls='-', label='Tinker+(2017)')
    
        err_jack = np.zeros(len(ssfr_tot)) 
        for jack_result in jack_results: 
            ssfr_jack = np.array([np.median(ssfrs) for ssfrs in jack_result[i_ssfr]]) 
            err_jack += (ssfr_jack - ssfr_tot)**2
        err_jack *= np.float(jack_bins[0] * jack_bins[1] - 1)/np.float(jack_bins[0] * jack_bins[1])
        sub.errorbar(0.5*(rperp_bins[:-1]+rperp_bins[1:]), ssfr_tot, yerr=np.sqrt(err_jack), 
                c=pretty_colors[i_ssfr], elinewidth=2, capsize=5) 
        ssfrplots.append(ssfrplot) 

    kleg = sub.legend(handles=[kauffplot, leg_ssfrplot], loc='upper right', handletextpad=0.2, 
            prop={'size': 20}, markerscale=50)
    for legobj in kleg.legendHandles:
        legobj.set_linewidth(5.)
    sub.legend(handles=ssfrplots, loc='lower right', ncol=2, handletextpad=0.2) 
    plt.gca().add_artist(kleg)


    label_pipe = 'tot; mpajhu'
    label_ssfr_neigh = 'fib; mpajhu'
    # axes
    sub.set_xlim([0, 4]) 
    sub.set_ylim([-12.25, -9.75])
    sub.set_xticks([0, 1, 2, 3, 4]) 
    sub.set_yticks([-12., -11., -10.]) 
    #sub.text(0.1, -11.9, r'Ranked in $\mathtt{log}\mathtt{SSFR^{('+label_pipe+')}}$', fontsize=20) 

    sub.minorticks_on() 

    #fig.title(r"$10.0 < \mathtt{log}\mathcal{M}^\mathtt{("+label_pipe+")}_* < 10.5$ ", fontsize=25) 
    sub.set_xlabel(r'$\mathtt{R_{\perp}}$ [Mpc]', labelpad=10, fontsize=30) 
    sub.set_ylabel(''.join([r'log($\mathtt{SSFR}$ [$\mathtt{yr}^{-1}$])']), labelpad=10, fontsize=30) 
    
    fig_file = ''.join([UT.dir_fig(), 
        'forTinker_review.pdf'])
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close()
    return None 


if __name__=='__main__': 
    #Fig_forProposal()
    #Fig_TinKauff_NeighborSSFR_rperp_PrimaryBins_IsoComparison():
    Fig_forTinkerReview_NeighborSSFR_rperp_PrimaryBins()
    #Fig_forTinker_NeighborSSFR_rperp_PrimaryBins()
    #Fig_TinKauff_Iso_NeighborSSFR_rperp_PrimaryBins()
    #Fig_NeighborSSFR_rperp_PrimaryBins(
    #        {'name': 'tinkauff', 'Mrcut':18, 'Mass_cut': 9.25, 
    #            'primary_delv': 500., 'primary_rperp': 0.5, 
    #            'neighbor_delv': 500., 'neighbor_rperp': 5.})
    #
    #Fig_NeighborSSFR_in_Primarybins({'name': 'tinker', 'Mrcut':18, 
    #    'primary_delv': 500., 'primary_rperp': 0.5, 
    #    'neighbor_delv': 500., 'neighbor_rperp': 5.},
    #    primary_massbin=[10., 10.5])
    #Fig_NeighborSSFR_in_PurePrimarybins({'name': 'tinker', 'Mrcut':18, 
    #    'primary_delv': 500., 'primary_rperp': 0.5, 
    #    'neighbor_delv': 500., 'neighbor_rperp': 5.},
    #    primary_massbin=[10., 10.5])
    #Fig_CentralMassbinNeighborSSFR_in_PurePrimarybins({'name': 'tinker', 'Mrcut':18, 
    #    'primary_delv': 500., 'primary_rperp': 0.5, 
    #    'neighbor_delv': 500., 'neighbor_rperp': 5.},
    #    primary_massbin=[10., 10.5])

    #Fig_CentralNeighborSSFR_in_PurePrimarybins({'name': 'tinker', 'Mrcut':18, 
    #    'primary_delv': 500., 'primary_rperp': 0.5, 
    #    'neighbor_delv': 500., 'neighbor_rperp': 5.},
    #    primary_massbin=[10., 10.5])
