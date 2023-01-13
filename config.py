# configuration file for the project
import numpy as np
import matplotlib as mpl
from matplotlib import rc
import astropy.units as u
from pathlib import Path

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
rc('font',**{'family':'serif'})#,'sans-serif':['Helvetica']})
rc('text', usetex=True)

mainfolder = str(Path.home()) + '/propstar_gasflow/'

datafolder = mainfolder + 'data/'
firstlookfolder = mainfolder + 'first_look/'
figfolder = mainfolder + 'figures/'

hc3n_10_9_cube = datafolder + 'NGC1333_HC3N_L24-merged'
hc3n_10_9_cube_s = hc3n_10_9_cube +'_small'

distance = 294 * u.pc # pc, Zucker et al 2018?
hc3n_10_9_rms = 1.1508e-2 # Jy/beam
velrange = np.array([-3, 18])
fwhm_to_sigma = 1. / (8 * np.log(2))**0.5

def beam_size(header):
    beam_fwhm_maj = header['BMAJ'] * u.deg
    beam_fwhm_min = header['BMIN'] * u.deg
    pixsize = np.abs(header['CDELT2']) * u.deg
    beam_sigma_maj = beam_fwhm_maj * fwhm_to_sigma
    beam_sigma_min = beam_fwhm_min * fwhm_to_sigma
    omega_B = 2 * np.pi * beam_sigma_maj * beam_sigma_min
    omega_B_pix2 = (omega_B / pixsize**2).value
    return omega_B, omega_B_pix2