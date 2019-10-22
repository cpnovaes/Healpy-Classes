#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 20 14:03:53 2017

@author: camila

            -- AULA 5 --
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
"""
STEP 2:Ring to nest

"""

Nside = 4
pix_ring = 25
pix_nest = hp.ring2nest(Nside, pix_ring)
print('pix_nest =', pix_nest)

all_pix_ring = np.arange(192)
all_pix_nest = hp.ring2nest(Nside, all_pix_ring)
hp.mollview(all_pix_nest, title='Nest scheme', unit='Pix index')
pyplot.savefig('mollweide_view_nest.png')

#%%
"""
STEP 3:Nest to ring

"""

# Nside = 4
print('pix_nest =', pix_nest)
pix_ring = hp.nest2ring(Nside, pix_nest)
print('pix_ring =', pix_ring)

all_pix_ring = hp.nest2ring(Nside, all_pix_nest)
hp.mollview(all_pix_ring, title='Ring scheme', unit='Pix index')
pyplot.savefig('mollweide_view_ring.png')

#%%
"""
STEP 4: Reordering a map

"""
mapa_ring = hp.read_map('COM_CMB_IQU-smica_1024_R2.02_full.fits')
hp.mollview(mapa_ring, title='CMB map - Ring', unit='K')
pyplot.savefig('mollweide_view_cmb.png')


mapa_nest = hp.reorder(mapa_ring, inp='RING', out='NEST')
hp.mollview(mapa_nest, title='CMB map - Nest', unit='K') # Errado, porque?

# OR:::

mapa_nest = hp.reorder(mapa_ring, r2n=True)
hp.mollview(mapa_nest, title='CMB map - Nest', unit='K') # Errado, porque?



#%%
"""
STEP 5: Nside / Npix / Resolution

"""
# Nside <--> Npix

Nside1 = 64 # = 1, 2, 4, 8, 16, 64, ...
Npix1  = 12*Nside**2 # 12 x Nside²
print('(1)', Nside1, '-->', Npix1)

# OR:::

Npix2  = hp.nside2npix(Nside1) # 12 x Nside²
print('(2)', Nside1, '-->', Npix2)

Nside2  = hp.npix2nside(Npix2) # 12 x Nside²
print('(3)', Nside2, '<--', Npix2)

# Nside <--> Resolution

Nside = 4 # = 1, 2, 4, 8, 16, 64, ...
Resolution  = hp.nside2resol(Nside)
print('(1)', Nside, '-->', Resolution, 'radians')

Resolution  = hp.nside2resol(Nside, arcmin=True)
print('(2)', Nside, '-->', Resolution, 'arcmin =', Resolution/60., 'deg')

Area  = hp.nside2pixarea(Nside, degrees=True)
print('(3)', Nside, '-->', Area, 'square deg')

# get_map_size / get_nside

mapa = hp.read_map('COM_CMB_IQU-smica_1024_R2.02_full.fits')

Npix_mapa  = hp.get_map_size(mapa)
Nside_mapa = hp.get_nside(mapa)

print('Nside =', Nside_mapa)
print('Npix =', Npix_mapa)


#%%
"""
STEP 6: Nside / Npix / Resolution

"""

mapa = hp.read_map('COM_CMB_IQU-smica_1024_R2.02_full.fits')
hp.mollview(mapa, title='CMB map - Nside=1024', unit='K')

mapa_deg = hp.ud_grade(mapa, 64, order_in = 'RING', order_out='NEST')
hp.mollview(mapa_deg, title='CMB map - Nside=64', unit='K', nest=True)
pyplot.savefig('mollweide_view_cmb_nside32.png')
