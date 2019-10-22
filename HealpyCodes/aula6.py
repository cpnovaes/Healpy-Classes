#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 16:29:14 2017

@author: camila

  -- Aula 6 --
"""

#%%
"""
STEP 1: Import the packages
"""

import healpy as hp
import matplotlib
from matplotlib import pyplot # ***
import numpy as np

fontsize = 20
matplotlib.rcParams.update({'font.size':fontsize})


#%%

# Preparing map to be used in the class:

mapa = hp.read_map('COM_CMB_IQU-smica_1024_R2.02_full.fits', field={0,3})
    
mapa2 = mapa[0]
mapa2n = hp.reorder(mapa2, r2n=True)
mapa2n[1000000:1030000] = hp.UNSEEN
mapa2r = hp.reorder(mapa2n, n2r=True)
hp.mollview(mapa2r, title='CMB map', unit='K')
pyplot.savefig('mollweide_view_bad.png')

hp.gnomview(mapa2r, rot=[40,80,0], reso=8, title='CMB map', unit='K')
pyplot.savefig('gnomonic_view_bad.png')

mapa_save = [mapa2r,mapa[1]]

hp.write_map('CMB_map_Ns1024.fits',mapa_save, coord='G', column_names=list(['I_CMB','Confidence-Mask']), column_units=list(['K',' ']))

#%%

"""
STEP 2: Missing value
"""
m = np.arange(12.)
m[2] = hp.UNSEEN

m_bad = hp.mask_bad(m)

m_good = hp.mask_good(m)

print('m_bad ->', m_bad)
print('m_good ->', m_good)

hp.mollview(m_bad, title='mask_bad')
hp.mollview(m_good, title='mask_good')

#%%

"""
STEP 3: Treating missing pixels
"""

mapa = hp.read_map('CMB_map_Ns1024.fits',field={0,1}, h=True)
cmb = mapa[0]
mask = mapa[1]

bad = hp.mask_bad(cmb)
good = hp.mask_good(cmb)

hp.mollview(bad, title='Bad')
hp.mollview(good, title='Good')

mask_total = mask * good
hp.mollview(mask_total, title='Final mask')
pyplot.savefig('mask_total.png')


#%%

"""
STEP 4: Masking maps
"""

cmb_masked = hp.ma(cmb)
hp.mollview(cmb)

cmb_masked.mask = np.logical_not(mask)
hp.mollview(cmb_masked, title='Masked CMB map with UNSEEN value')
pyplot.savefig('masked_cmb_unseen.png')

cmb_masked0 = cmb_masked.filled()
cmb_masked0

#import matplotlib.pyplot as plt
#plt.hist(cmb_masked0.compressed(), bins = 1000)


#%%

"""
STEP 5: Removing monopole and dipole
"""

mapa = hp.read_map('LFI_CompMap_Foregrounds-smica_1024_R2.00.fits')

hp.mollview(mapa, norm='hist', title='all multipoles')
#pyplot.savefig('foreg_all_multipoles.png')

print('<mapa> =', np.mean(mapa))

####################
# Removing DIPOLE:::

mapa2 = hp.remove_monopole(mapa)

print('<mapa2> =', np.mean(mapa2))

hp.mollview(mapa2, norm='hist', title='l >= 1')
#pyplot.savefig('foreg_lge1.png')

mapa_mono = mapa - mapa2
print(mapa_mono)

hp.mollview(mapa_mono, norm='hist', title='l = 0 (monopole)')
pyplot.savefig('foreg_l1.png')

# with and without gal_cut ... 
# Compare:::

mapa2 = hp.remove_monopole(mapa, gal_cut=20)
# legal fazer isso em um mapa de foreground -> o monopolo varia bastante

####################
# Removing DIPOLE:::

mapa = hp.read_map('LFI_CompMap_Foregrounds-smica_1024_R2.00.fits')
mapa3 = hp.remove_dipole(mapa)
hp.mollview(mapa3, norm='hist', title='l > 1')
pyplot.savefig('foreg_lgt1.png')

mapa = hp.read_map('LFI_CompMap_Foregrounds-smica_1024_R2.00.fits')
mapa_mono_dip = mapa - mapa3
hp.mollview(mapa_mono_dip, norm='hist', title='0 <= l <= 1 (mono+dip)')
pyplot.savefig('foreg_0gelle1.png')

###
mapa3 = hp.remove_dipole(mapa, gal_cut= 20 )
hp.mollview(mapa3, norm='hist', title='l > 1 (gal_cut=20)')
pyplot.savefig('foreg_cut_lgt1.png')

mapa_mono_dip = mapa - mapa3
hp.mollview(mapa_mono_dip, norm='hist', title='0 <= l <= 1 (mono+dip, gal_cut=20)')
pyplot.savefig('foreg_cut_0gelle1.png')

# if fitval = True mapa2[1] e mapa2[2] returns monopole and dipole vectors.

mapa3 = hp.remove_dipole(mapa, fitval = True)
print(mapa2[1])
print(mapa2[2])

#%%
# Important example::::::::
# ##### WITH THE MASK #####

mapa = hp.read_map('CMB_map_Ns1024.fits',field={0,1}, h=True)
cmb = mapa[0]
mask = mapa[1]

mapa1 = hp.remove_dipole(cmb)
hp.mollview(cmb, title='CMB all multipoles')
hp.mollview(mapa1, title='CMB l > 1')
hp.mollview(cmb-mapa1, title='CMB 0 >= l >=1')

###

cmb_masked = hp.ma(cmb)
cmb_masked.mask = np.logical_not(mask)
hp.mollview(cmb_masked, title='Masked-CMB all multipoles')

mapa2 = hp.remove_dipole(cmb_masked)
hp.mollview(mapa2)
hp.mollview(cmb_masked-mapa2, title='Masked-CMB 0 >= l >=1')


#%%

"""
STEP 6: Finding neighbours
"""

nside = 128
theta = np.pi/2.
phi = np.pi

pix = hp.get_all_neighbours(nside, theta, phi)

# if phi is not given or None, theta is interpreted as pixel number.
pix = hp.get_all_neighbours(nside,100)

hp.get_all_neighbours(1024,5432)
# Out: array([5642, 5431, 5225, 5023, 5226, 5433, 5643, 5857])

