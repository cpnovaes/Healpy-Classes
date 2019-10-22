#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 27 11:20:22 2017

@author: camila

  -- Aula 7 --

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
STEP 2: Pix index <--> Angles
"""

# For one pixel:::

ang = hp.pix2ang(16, 1440)
print('Pix = 1440 --> (Theta,Phi) =',ang)

# For the whole map:::

nside = 64
npix = hp.nside2npix(nside)
pixels = np.arange(npix)
theta_phi = hp.pix2ang(nside, pixels)

print(len(theta_phi))

theta = np.rad2deg(theta_phi[0])
phi = np.rad2deg(theta_phi[1])

hp.mollview(theta, title='Projection of Theta variation', unit='Degrees')
pyplot.savefig('theta_map.png')

hp.mollview(phi, title='Projection of Phi variation', unit='Degrees')
pyplot.savefig('phi_map.png')

# Inveting:::

theta0 = np.deg2rad(45) # rad
phi0 = np.deg2rad(30)    # rad
ipix0 = hp.ang2pix(4, theta0, phi0)
print('(theta,phi)=(',theta0,',',phi0,') --> ipix =',ipix0)

pixels_calc = hp.ang2pix(nside, np.deg2rad(theta), np.deg2rad(phi))

hp.mollview(pixels_calc, title='Pixels index')
pyplot.savefig('ipix_map.png')

#%%
"""
STEP 3: Pix index <--> Vec
"""

vec1 = hp.pix2vec(16, 1504)
print('Pix = 1504 --> (x,y,z) =',vec1)

vec2 = hp.pix2vec(16, [1440,  427])
print(vec2)

print('Pix = 1440 --> (x,y,z) = (',vec2[0][0], vec2[1][0], vec2[2][0],')')
print('Pix = 427 --> (x,y,z) = (', vec2[0][1], vec2[1][1], vec2[2][1],')')

vec3 = hp.pix2vec([1, 2], 11)
print(vec3)

print('Nside=1 / Pix = 11  --> (x,y,z) = (',vec3[0][0], vec3[1][0], vec3[2][0],')')
print('Nside=2 / Pix = 11  --> (x,y,z) = (',vec3[0][1], vec3[1][1], vec3[2][1],')')

# Inveting:::

ipix2 = hp.vec2pix(16, vec2[0],vec2[1], vec2[2])

# How do x,y,z variate through the sky?

nside = 64
npix = hp.nside2npix(nside)
pixels = np.arange(npix)

vec = hp.pix2vec(nside, pixels)

hp.mollview(vec[0], title='Projection of x variation', unit='dimensionless')
pyplot.savefig('xVec_map.png')

hp.mollview(vec[1], title='Projection of y variation', unit='dimensionless')
pyplot.savefig('yVec_map.png')

hp.mollview(vec[2], title='Projection of z variation', unit='dimensionless')
pyplot.savefig('zVec_map.png')

#%%
"""
STEP 3: Angles <--> Vec
"""
theta = np.deg2rad(45.) # rad
phi   = np.deg2rad(30.) # rad

vec = hp.ang2vec(theta,phi)
vec
print('(theta,phi) = (',theta,',',phi,') --> (x,y,z)=', vec)

theta_phi = hp.vec2ang(vec)
theta_phi
print('(x,y,z)=', vec,' --> (theta,phi) = (',theta_phi,')')


hp.dir2vec(theta, phi=phi)

#%%
"""
STEP 4: Angles <--> Vec ::OR:: Dir <--> Vec
"""

# vec2dir and dir2vec are very similar to vec2ang and ang2vec
# Which are the difference???

# Create a theta,phi list:
nside = 64
npix = hp.nside2npix(nside)
pixels = np.arange(npix)
theta_phi = hp.pix2ang(nside, pixels)

# ang2vec vs. dir2vec
vec = hp.ang2vec(theta_phi[0],theta_phi[1]) # it works
vec
vec = hp.ang2vec(theta_phi)                 # it DOES NOT work


vec = hp.dir2vec(theta_phi[0],phi=theta_phi[1]) # it works!!!
vec
vec = hp.dir2vec(theta_phi)                 # it works!!!
vec

theta = 30.        # Colatitude
lat = 90. - theta  # Latitude  b
phi = 25.          # Longitude l
lon = phi

vec0 = hp.ang2vec(np.deg2rad(theta),np.deg2rad(phi))
vec0

vec1 = hp.dir2vec(lon,phi=lat, lonlat = True)
vec1

# Virgo cluster:::

lon = 283.8
lat = 74.4
vec = hp.dir2vec(lon,phi=lat, lonlat = True)

#%%
"""
STEP 5: Angular distances
"""

# Distance from Virgo to Coma:::

lon_c = 235.1
lat_c = 73.0

vec_c = hp.dir2vec(lon_c,phi=lat_c, lonlat = True)

distance = hp.rotator.angdist(vec, vec_c, lonlat=True)

#%%
"""
STEP 6: Rotation
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
pixels = np.arange(npix)
hp.mollview(pixels,title='Galactic coods.')
pyplot.savefig('ipix_Galactic.png')

theta_phi = hp.pix2ang(nside, pixels)
theta_ecl, phi_ecl = r(theta_phi[0], theta_phi[1])

ipix = hp.ang2pix(nside, theta_ecl, phi_ecl)
hp.mollview(ipix,title='Ecliptic coords.')
pyplot.savefig('ipix_Ecliptic.png')

##############
# For a map:::
mapa = hp.read_map('LFI_CompMap_Foregrounds-smica_1024_R2.00.fits')
hp.mollview(mapa, title='Foreground map in Galactic coods.', norm='hist')

mapa2 = []
for i in ipix:
    mapa2.append(mapa[i])

mapa2 = np.array(mapa2)
hp.mollview(mapa2, title='Foreground map in Ecliptic coods.', norm='hist')

########################
# Performing rotation:::
r = hp.Rotator(rot=[50,310,30])
theta_rot, phi_rot = r(theta_phi[0], theta_phi[1])

ipix_rot = hp.ang2pix(nside, theta_rot, phi_rot)
hp.mollview(ipix_rot)

mapa3 = []
for i in ipix_rot:
    mapa3.append(mapa[i])

mapa3 = np.array(mapa3)
hp.mollview(mapa3, title='Foreground map rotated', norm='hist')

