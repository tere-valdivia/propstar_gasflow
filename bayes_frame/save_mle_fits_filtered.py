import numpy as np
import pandas as pd
from astropy.io import fits
import os
import sys
from skimage import morphology
from astropy.wcs import WCS
sys.path.append('../')
from setup import *
from config import *
import glob
from spectral_cube import SpectralCube

def header_flatten(head):
    """ Flattens 3D header into 2D. Surely it exists somewhere already... """
    flathead = head.copy()
    for key in flathead.keys():
        if key.endswith('3'):
            flathead.pop(key)
    flathead['NAXIS'] = 2
    flathead['WCSAXES'] = 2

    return flathead

molecule = 'N2Hp'
ngaussmapfile = 'nested-sampling/npeaks_cut5.fits'

mlex1, mle1head = fits.getdata('nested-sampling/NGC1333-SE-mle-x1.fits', header=True)
mlex2, mle2head = fits.getdata('nested-sampling/NGC1333-SE-mle-x2.fits', header=True)

Kfile = 'nested-sampling/NGC1333-SE-Ks.fits'
Kcut = 5  # the heuristical ln(Z1/Z2) cut for model selection

Ks = fits.getdata(Kfile)
# make the ln(K)>Kcut based map of LoS component numbers
npeaks_map = np.zeros(np.shape(Ks[0]))
npeaks_map[Ks[0] <= Kcut] = 0
npeaks_map[Ks[0] > Kcut] = 1
Karr_clean = morphology.remove_small_objects(npeaks_map.astype(bool), min_size=101).astype(int)
fits.writeto('nested-sampling/mask_5_s101.fits', Karr_clean, header=header_flatten(fits.getheader(Kfile)), overwrite=True)
Karr2_clean = np.zeros(np.shape(Ks[0]))
Karr2_clean[(Ks[0] > Kcut) & (Ks[1] <= Kcut)] = 0
Karr2_clean[(Ks[0] > Kcut) & (Ks[1] > Kcut)] = 1
Karr2_clean = morphology.remove_small_objects(Karr2_clean.astype(bool), min_size=7).astype(int)
fits.writeto('nested-sampling/mask_5_2ndcomp_s7.fits', Karr2_clean, header=header_flatten(fits.getheader(Kfile)), overwrite=True)

npeaks_map[(Karr_clean == 0)] = 0
npeaks_map[np.where(np.isnan(Ks[0]))] = np.nan
npeaks_map[(Ks[0] > Kcut) & (Karr_clean == 1)] = 1
npeaks_map[(Ks[0] > Kcut) & (Ks[1] > Kcut) & (Karr2_clean == 1)] = 2

# for i in range(1, Ks_lines-1):
#     npeaks_map[(Ks[i-1] > Kcut) & (Ks[i] > Kcut)] = i+1
    

hdu_npeaks = fits.PrimaryHDU(npeaks_map,
                             header=header_flatten(fits.getheader(Kfile)))
hdu_npeaks.writeto('nested-sampling/npeaks_cut5_noislands.fits', overwrite=True)

# here we save the MLE results only where they correspond
mask1G = np.where((npeaks_map>0) & (npeaks_map<2), 1, 0)
mask1G = np.array([mask1G.astype(bool)]*6)
mask2G = np.array([Karr2_clean.astype(bool)]*12)
mlex1[~mask1G] = np.nan
mlex2[~mask2G] = np.nan

# now we filter by error



fits.writeto('nested-sampling/NGC1333-SE-mle-x1_filtered.fits', mlex1, mle1head, overwrite=True)
fits.writeto('nested-sampling/NGC1333-SE-mle-x2_filtered.fits', mlex2, mle2head, overwrite=True)




