import numpy as np
import matplotlib.pyplot as plt

# --- local ---
import util as UT
import catalog as clog

from ChangTools.plotting import prettyplot
from ChangTools.plotting import prettycolors

def Test_KauffmannParent_Pssfr(): 
    ''' look at the P(sSFR) of the Kauffmann et al. sample
    '''
    catalog = clog.KauffmannParent() 
    is_notnan = np.isnan(catalog['ssfr_fib']) == False 
    fib_ssfr = catalog['ssfr_fib'][is_notnan]
    print len(catalog['ssfr_fib']) - np.sum(is_notnan), ' out of ', len(catalog['ssfr_fib']), ' galaxies have NaN SSFR' 

    ssfr_min, ssfr_max = -13., -8.75
    pdf, bins = np.histogram(fib_ssfr, 
            range=[ssfr_min, ssfr_max], bins=20, normed=True) 

    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure()
    sub = fig.add_subplot(111)
    sub.plot(0.5*(bins[:-1]+bins[1:]), pdf, c=pretty_colors[3], lw=3)
    sub.text(-12.8, 0.5, 'median sSFR = '+str(round(np.median(fib_ssfr),3)), fontsize=20)
    sub.vlines(np.median(fib_ssfr), 0., 10., 
            color=pretty_colors[1], linewidth=3, linestyle='--')
    sub.set_title('Kauffmann et al.(2013) cuts on dr72bright34')
    # axes
    sub.set_xlabel(r'log(SSFR)', fontsize=25)
    sub.set_xlim([ssfr_min, ssfr_max]) 
    sub.set_xticks([-13, -11, -9]) 
    sub.set_ylabel(r'P(SSFR)', fontsize=25)
    sub.set_ylim([0.0, 0.6])
    sub.set_yticks([0., 0.2, 0.4, 0.6]) 
    sub.minorticks_on()
    fig_file = ''.join([UT.dir_fig(), 'KauffmannParent.Pssfr.png']) 
    fig.savefig(fig_file) 
    plt.close() 
    return None


def Test_Jackknife(n_jack, RADec_bins=[3,3]):  
    ''' Test the Jackknifing
    '''
    cat_dict = {'name': 'tinker', 'Mrcut':18, 
            'primary_delv': 500., 'primary_rperp': 0.5, 
            'neighbor_delv': 500., 'neighbor_rperp': 5.}
    # read conformity catalog based on input catalog dictionary
    concat = clog.ConformCatalog(Mrcut=cat_dict['Mrcut'], 
            primary_delv=cat_dict['primary_delv'], primary_rperp=cat_dict['primary_rperp'],  
            neighbor_delv=cat_dict['neighbor_delv'], neighbor_rperp=cat_dict['neighbor_rperp'])
    catalog = concat.Read() 
    jack_catalog = concat.Jackknife(catalog, n_jack, RADec_bins=RADec_bins) 
    #jack_catalog = concat.ReadJackknife(n_jack, RADec_bins=RADec_bins)

    prettyplot()
    pretty_colors = prettycolors() 
    fig = plt.figure() 
    sub = fig.add_subplot(111)

    sub.scatter(catalog['ra'], catalog['dec'], s=6, lw=0, c='k') 
    sub.scatter(jack_catalog['ra'], jack_catalog['dec'], s=6, lw=0, c=pretty_colors[3]) 
    
    # axes
    sub.set_xlabel('RA', fontsize=25) 
    sub.set_xlim([-50, 400])
    sub.set_ylabel('Dec', fontsize=25) 
    sub.set_ylim([-20, 80])
    sub.minorticks_on() 

    fig_file = ''.join([UT.dir_fig(), 
        'test_jackknife.', str(n_jack), 'of', 
        str(RADec_bins[0]), 'x', str(RADec_bins[1]), '.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close() 
    return None 


def Test_PrimaryIdentify(del_v_cut=500., r_perp_cut=0.5): 
    ''' Test the primary identification criteria within the RA and Dec plane 
    '''
    # import that confomrity catalog
    concat = clog.ConformCatalog(Mrcut=18, 
            primary_delv=del_v_cut, primary_rperp=r_perp_cut) 
    catalog = concat.Read() 
    
    sat_indices = catalog['satellite_indices']  # satellite indicies
    is_primary = np.where(catalog['primary'] == 1) 

    prettyplot()
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(8,8)) 
    bkgd = fig.add_subplot(111, frameon=False)

    for ii in range(4): 
        sub = fig.add_subplot(2,2,ii+1)
        if ii == 0: 
            i_cen = np.argmax(catalog['n_satellite'][is_primary])
        else: 
            i_cen = np.random.choice(len(is_primary[0]))
            if len(sat_indices[is_primary[0][i_cen]]) == 0:  
                while len(sat_indices[is_primary[0][i_cen]]) == 0: 
                    i_cen = np.random.choice(len(is_primary[0]))
    
        sat_sizes = 100. + 100*(catalog['mass'][sat_indices[is_primary[0][i_cen]]] - catalog['mass'][is_primary[0][i_cen]])
        sub.scatter(
                catalog['ra'][sat_indices[is_primary[0][i_cen]]], 
                catalog['dec'][sat_indices[is_primary[0][i_cen]]], 
                c=pretty_colors[1], lw=0, 
                s=sat_sizes)
        sub.scatter(
                np.repeat(catalog['ra'][is_primary[0][i_cen]],2), 
                np.repeat(catalog['dec'][is_primary[0][i_cen]],2), c='r', lw=0, 
                s=100.) 
        sub.set_xlim([
            catalog['ra'][is_primary[0][i_cen]]-.75, 
            catalog['ra'][is_primary[0][i_cen]]+.75
            ])
        sub.set_ylim([
            catalog['dec'][is_primary[0][i_cen]]-.75, 
            catalog['dec'][is_primary[0][i_cen]]+.75
            ])

    bkgd.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    bkgd.set_ylabel('RA', fontsize=25)
    bkgd.set_xlabel('Dec', fontsize=25)
    
    fig_file = ''.join([UT.dir_fig(), 
        'test_PrimaryIdentify.', 'delv', str(del_v_cut), '.rperp', str(r_perp_cut), '.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close() 
    return None 


def Test_NeighborIdentify(del_v_cut=500., r_perp_cut=5.): 
    ''' Test the primary identification criteria
    '''
    concat = clog.ConformCatalog(Mrcut=18, 
            neighbor_delv = del_v_cut, neighbor_rperp=r_perp_cut
            ) 
    catalog = concat.Read() 
    
    sat_indices = catalog['satellite_indices']  
    neigh_indices = catalog['neighbor_indices']
    is_primary = np.where(catalog['primary'] == 1) 

    prettyplot()
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(8,8)) 
    bkgd = fig.add_subplot(111, frameon=False)

    for ii in range(4): 
        sub = fig.add_subplot(2,2,ii+1)
        if ii == 0: 
            i_cen = np.argmax(catalog['n_neighbor'])
        else: 
            i_cen = np.random.choice(len(is_primary[0]))

        neigh_sizes = 100.*(1. + 
                catalog['mass'][neigh_indices[i_cen]] - 
                catalog['mass'][is_primary[0][i_cen]])
        sub.scatter(
                catalog['ra'][neigh_indices[i_cen]], 
                catalog['dec'][neigh_indices[i_cen]], 
                c=pretty_colors[1], lw=0, 
                s=neigh_sizes)
    
        print len(sat_indices[is_primary[0][i_cen]])
        sat_sizes = 100.*(1.+
                catalog['mass'][sat_indices[is_primary[0][i_cen]]] - 
                catalog['mass'][is_primary[0][i_cen]])
        sub.scatter(
                catalog['ra'][sat_indices[is_primary[0][i_cen]]], 
                catalog['dec'][sat_indices[is_primary[0][i_cen]]], 
                c='g', lw=0, 
                s=sat_sizes)

        sub.scatter(
                np.repeat(catalog['ra'][is_primary[0][i_cen]],2), 
                np.repeat(catalog['dec'][is_primary[0][i_cen]],2), c='r', 
                lw=0, s=100) 

        #sub.text(
        #        catalog['ra'][is_primary[0][i_cen]]-3.5, 
        #        catalog['dec'][is_primary[0][i_cen]]+3.5, 
        #        r"$z_{primary} = $"+str(round(catalog['z'][is_primary[0][i_cen]],3)), 
        #        fontsize=20) 
        sub.set_xlim([
            catalog['ra'][is_primary[0][i_cen]]-2.5, 
            catalog['ra'][is_primary[0][i_cen]]+2.5
            ])
        sub.set_ylim([
            catalog['dec'][is_primary[0][i_cen]]-2.5, 
            catalog['dec'][is_primary[0][i_cen]]+2.5
            ])

    bkgd.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    bkgd.set_ylabel('RA', fontsize=25)
    bkgd.set_xlabel('Dec', fontsize=25)

    fig_file = ''.join([UT.dir_fig(), 
        'test_NeighborIdentify.', 'delv', str(del_v_cut), '.rperp', str(r_perp_cut), '.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close() 
    return None 


def Test_KauffTink_PrimaryIdentification(Mrcut=18, del_v_cut=500., r_perp_cut=0.5):
    '''
    '''
    concat = clog.ConformCatalog(Mrcut=18, 
            primary_delv=del_v_cut, primary_rperp=r_perp_cut, 
            neighbor_delv=500., neighbor_rperp=5., 
            mpajhu=False)
    tink = concat.Read()
    
    concat = clog.ConformCatalog(Mrcut=18, 
            primary_delv=del_v_cut, primary_rperp=r_perp_cut, 
            neighbor_delv=500., neighbor_rperp=5., 
            mpajhu=True)
    kauff = concat.Read()

    tink_isprimary = np.where(tink['primary'] == 1)
    tink_primary_indices = tink_isprimary[0]

    kauff_isprimary = np.where(kauff['primary'] == 1)
    kauff_primary_indices = kauff['mpajhu_tinker_index'][kauff_isprimary]

    print 'tink ', len(tink_primary_indices)
    print 'kauff ', len(kauff_primary_indices)

    common_primaries = np.intersect1d(kauff_primary_indices, tink_primary_indices)
    print 'common ', len(common_primaries) 


def MPAJHU_Tinker(Mrcut=18): 
    ''' Compare the stellar masses inferred from MPA-JHU to the K-correct 
    stellar masses of the VAGC catalog used in Tinker et al. (2011) catalogs. 
    '''
    tink = clog.TinkerCatalog(Mrcut=Mrcut) 

    kauff = clog.MPAJHU_TinkerCatalog(Mrcut=Mrcut) 
    kauff_indices = kauff['mpajhu_tinker_index'] 

    prettyplot()
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(8,8)) 
    for i in range(4):  
        if i == 0: 
            prop = 'mass'
            totorfib = 'tot'
        elif i == 1: 
            prop = 'mass'
            totorfib = 'fib'
        elif i == 2: 
            prop = 'ssfr'
            totorfib = 'tot'
        elif i == 3: 
            prop = 'ssfr'
            totorfib = 'fib'

        sub = fig.add_subplot(2,2,i+1) 
        
        if prop == 'mass': 
            sub.scatter(tink['mass'][kauff_indices], kauff['mass_'+totorfib+'_mpajhu'], 
                    s=6, lw=0, c=pretty_colors[1]) 
            sub.plot([9., 12.], [9., 12.], c='k', lw=3, ls='--') 
            sub.set_xlim([9., 12.]) 
            sub.set_xticks([9., 10., 11., 12.]) 
            sub.set_xlabel(r'VAGC $\mathcal{M}_*$', fontsize=15) 
            sub.set_ylim([9., 12.]) 
            sub.set_yticks([9., 10., 11., 12.]) 
            sub.set_ylabel(r'MPA-JHU $\mathcal{M}_*^\mathtt{'+totorfib+'}$', fontsize=15) 
        elif prop == 'ssfr': 
            sub.scatter(tink['ssfr'][kauff_indices], kauff['ssfr_'+totorfib+'_mpajhu'], 
                    s=6, lw=0, c=pretty_colors[1]) 
            sub.plot([-14., -8], [-14., -8.], c='k', lw=3, ls='--') 
            sub.set_xlim([-14., -8.]) 
            sub.set_xticks([-14., -12, -10, -8.]) 
            sub.set_xlabel(r'VAGC $\mathtt{log}\mathtt{SSFR}$', fontsize=15) 
            sub.set_ylim([-14., -8.]) 
            sub.set_yticks([-14., -12, -10, -8.]) 
            sub.set_ylabel(r'MPA-JHU $\mathtt{log}\mathtt{SSFR^{'+totorfib+'}}$', fontsize=15) 

        if i in [1,3]: 
            sub.yaxis.tick_right() 
            sub.yaxis.set_label_position('right') 
        sub.minorticks_on() 
    
    fig.subplots_adjust(wspace=0.15, hspace=0.3)
    fig_file = ''.join([UT.dir_fig(), 'MPHJHU_Tinker.galprop_comparison.png']) 
    fig.savefig(fig_file, bbox_inches='tight') 
    plt.close() 
    return None



if __name__=='__main__': 
    #clog.Build_VAGCdr72_MPAJHU(Ascii=True)
    clog.Build_TinKauffGroupCat(Mass_cut=9.25)
    #clog.Build_TinKauff_IsolationGroupCat(Mass_cut=9.25)
    clog.Build_TinKauffGroupCat(Mass_cut=10.0)
    #Test_KauffmannParent_Pssfr()
    #for n in [1, 8, 10, 22]:
    #    Test_Jackknife(n, RADec_bins=[5,5])
    #clog.Build_MPAJHU_TinkerCatalog(Mrcut=18)
    #clog.Build_KauffmannParent()
    #clog.Build_VAGCdr72bright34()
    #MPAJHU_Tinker(Mrcut=18)
    #Test_PrimaryIdentify(del_v_cut=500., r_perp_cut=0.5)
    #Test_NeighborIdentify(del_v_cut=500., r_perp_cut=5.)
    #Test_KauffTink_PrimaryIdentification(Mrcut=18, del_v_cut=500., r_perp_cut=0.5)
    #clog.Build_TinkerCatalog(Mrcut=18)
    #clog.Build_MPAJHU_TinkerCatalog(Mrcut=18)
    for delv in [500.]:#, 1000., 1500., 2000.]:
        #tink_concat = clog.ConformCatalog('kauff',  
        tink_concat = clog.ConformCatalog('tinkauff', catalog_prop={'Mass_cut': 10.0}, 
                primary_delv=500., primary_rperp=0.5, 
                neighbor_delv=delv, neighbor_rperp=5.)
        tink_concat.Build() 
        tink_concat = clog.ConformCatalog('tinkauff', catalog_prop={'Mass_cut': 9.25}, 
                primary_delv=500., primary_rperp=0.5, 
                neighbor_delv=delv, neighbor_rperp=5.)
        #tink_concat = clog.ConformCatalog('tinkauff_iso', catalog_prop={'Mass_cut': 9.25}, 
        #        primary_delv=500., primary_rperp=0.5, 
        #        neighbor_delv=delv, neighbor_rperp=5.)
        tink_concat.Build() 
