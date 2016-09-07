import os
import numpy as np

import astropy.cosmology as astrocosmo


def dir_code(): 
    # Directory address to the code directory 
    return os.path.dirname(os.path.realpath(__file__))


def dir_dat():
    # Directory address to the dat directory (dat directory within 
    # repo is a symlink 
    return os.path.dirname(os.path.realpath(__file__)).split('code')[0]+'dat/'


def dir_fig():
    # Directory address to the dat directory (dat directory within 
    # repo is a symlink 
    return os.path.dirname(os.path.realpath(__file__)).split('code')[0]+'fig/'


def radecz_to_xyz(ra, dec, red, H0=70., Om0=0.3): 
    ''' convert RA, Dec, z to x, y, z
     default cosmology consistent with Kauffmann et al. (2013) choices
    '''
    cosmos = astrocosmo.FlatLambdaCDM(H0=H0, Om0=Om0)  
    r = cosmos.comoving_distance(red).value
    phi = ra
    theta = 90. - dec

    DRADEG = 180.0/np.pi
    x = r * np.cos(phi / DRADEG) * np.sin(theta / DRADEG)
    y = r * np.sin(phi / DRADEG) * np.sin(theta / DRADEG)
    z = r * np.cos(theta / DRADEG)

    return np.array([x, y, z]).T
