#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 18:08:01 2017

@author: camila

  -- Aula 8 --

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
STEP 2: Angular distances
"""

# Virgo cluster:::

lon = 283.8
lat = 74.4
vec = hp.dir2vec(lon,phi=lat, lonlat = True)


# Distance from Virgo to Coma:::

lon_c = 235.1
lat_c = 73.0

vec_c = hp.dir2vec(lon_c,phi=lat_c, lonlat = True)

distance = hp.rotator.angdist(vec, vec_c, lonlat=True)
print('Distance: ', distance)

#%%
"""
STEP 3: Rotation
"""

#############
# Examples:::
    
r = hp.Rotator(coord=['G','E'])  # Transforms galactic to ecliptic coordinates
theta_gal, phi_gal = np.pi/2., 0.

theta_ecl, phi_ecl = r(theta_gal, phi_gal)  # Apply the conversion
print(theta_ecl,phi_ecl)

theta_ecl, phi_ecl = hp.Rotator(coord='ge')(theta_gal, phi_gal) # In one line

vec_gal = np.array([1, 0, 0]) #Using vectors
vec_ecl = r(vec_gal)
print(vec_ecl)

######################
# For the whole sky:::
    
nside = 1024
npix = hp.nside2npix(nside)
pixels = np.arange(npix) # [1] Array with the pixel indexes.
hp.mollview(pixels,title='Galactic coods.')
#pyplot.savefig('ipix_Galactic.png')

theta_phi = hp.pix2ang(nside, pixels) # [2] theta, phi of all the pixels
theta_ecl, phi_ecl = r(theta_phi[0], theta_phi[1]) # [3] Rotate to find the 
                                                   # new Eclipt. coords.

ipix = hp.ang2pix(nside, theta_ecl, phi_ecl) # [4] Convert them back to pixel indexes
hp.mollview(ipix,title='Ecliptic coords.') # [5] Visualize!
#pyplot.savefig('ipix_Ecliptic.png')

##############
# For a map:::
mapa = hp.read_map('LFI_CompMap_Foregrounds-smica_1024_R2.00.fits')
hp.mollview(mapa, title='Foreground map in Galactic coods.', norm='hist')

mapa2 = []                  # Applying to a map.
for i in ipix:
    mapa2.append(mapa[i])

mapa3 = mapa[ipix] #************

mapa2 = np.array(mapa2)
hp.mollview(mapa2, title='Foreground map in Ecliptic coods.', norm='hist')

########################
# Performing rotation:::
r = hp.Rotator(rot=[50,310,30])
theta_rot, phi_rot = r(theta_phi[0], theta_phi[1])

ipix_rot = hp.ang2pix(nside, theta_rot, phi_rot)
hp.mollview(ipix_rot)

mapa3 = []                  # Applying to a map.
for i in ipix_rot:
    mapa3.append(mapa[i])

mapa3 = np.array(mapa3)
hp.mollview(mapa3, title='Foreground map rotated', norm='hist')

