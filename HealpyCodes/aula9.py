#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 17:09:55 2017

@author: camila

  -- Aula 9 --

"""
#%%
"""
STEP 1: Import the packages
"""

import healpy as hp
import matplotlib
from matplotlib import pyplot # ***
import numpy as np

# The following comands allow us to define the size of the characteres 
# at the mollview visualization.
fontsize = 20
matplotlib.rcParams.update({'font.size':fontsize})

#%%
"""
STEP 2: Spherical harmonic transforms: map --> Cls
"""
###################################
# Calculating Cls from a CMB map:::

# Just for temperature:::
mapa = hp.read_map('COM_CMB_IQU-smica_1024_R2.02_full.fits')

Cls_TT = hp.anafast(mapa)

# For temperature and polarization:::
mapa = hp.read_map('COM_CMB_IQU-smica_1024_R2.02_full.fits', field=[0,1,2])

Cls = hp.anafast(mapa)

Cls_TT = Cls[0]
Cls_EE = Cls[1]
Cls_BB = Cls[2]
Cls_TE = Cls[3]
Cls_EB = Cls[4]
Cls_TB = Cls[5]


#########
# Plot:::

# TT:::
lmax = len(Cls_TT)
ell = np.arange(lmax)

Dls_TT = ell*(ell+1)*Cls_TT/(2.*np.pi)

pyplot.plot(ell, Dls_TT,linewidth=2.0, color="red") #label="CMB Cls"
#
#plt.ylim(-300.,1000.)
#pyplot.xscale('log')
pyplot.title('Calculated Angular Power Spectrum - TT',fontsize=20)
pyplot.xlabel('Multipole, $\ell$',fontsize=20)
pyplot.ylabel('$[\ell(\ell + 1)/2\pi] C_\ell^{TT}$ $[K^2]$',fontsize=20)
#pyplot.legend(loc='best')
#
fig = pyplot.gcf()
fig.set_size_inches(10, 6)
pyplot.savefig("Cls.png")
pyplot.show()

#-----
res = hp.nside2resol(1024, arcmin = True)
print('Resolution=', res, 'arcmin')

theta = (180./lmax)*60.
print('Theta=', theta, 'arcmin')
#-----

# EE:::

pyplot.plot(ell, Cls_EE,linewidth=2.0, color="red") # label="CMB Cls"
#
pyplot.title('Calculated Angular Power Spectrum - EE',fontsize=20)
pyplot.xlabel('Multipole, $\ell$',fontsize=20)
pyplot.ylabel('$C_\ell^{EE}$ $[K^2]$',fontsize=20) # ($\mu K^2$)')
pyplot.legend(loc='best')

fig = pyplot.gcf()
fig.set_size_inches(10, 6)
pyplot.savefig("Cls_EE.png")
pyplot.show()


# BB:::
# TE:::

Dls_TE = ell*(ell+1)*Cls_TE/(2.*np.pi)

pyplot.plot(ell, Dls_TE,linewidth=2.0, color="red") # label="CMB Cls"
#
pyplot.title('Calculated Angular Power Spectrum - EE',fontsize=20)
pyplot.xlabel('Multipole, $\ell$',fontsize=20)
pyplot.ylabel('$[\ell(\ell + 1)/2\pi] C_\ell^{TE}$ $[K^2]$',fontsize=20) # ($\mu K^2$)')
pyplot.legend(loc='best')

fig = pyplot.gcf()
fig.set_size_inches(10, 6)
pyplot.savefig("Cls_TE.png")
pyplot.show()


# EB:::
# TB:::

#%%
"""
STEP 3: Spherical harmonic transforms: Cls --> fits file (writing the Cls in a fits file)
"""
# For one Cls array:
hp.write_cl('calc_Cls.fits', Cls[0])

# For all:
hp.write_cl('calc_Cls.fits', Cls)


#%%
"""
STEP 4: Spherical harmonic transforms: fits file --> Cls (reading Cls from a fits file)
"""
# To open just the FIRST extension:
Cls2 = hp.read_cl('COM_PowerSpect_CMB_R2.02.fits', h=True) #h=True or False does not work.


# To access the header and chose the extension you want to open:
from astropy.io import fits
hdulist = fits.open('COM_PowerSpect_CMB_R2.02.fits')
hdulist.info()

Cls_TTLOLUNB = hp.read_cl(hdulist[1])
Cls_TTHILUNB = hp.read_cl(hdulist[8])

n1 = len(Cls_TTLOLUNB[0])
n2 = len(Cls_TTHILUNB[0])
n = n1 + n2

ell=np.zeros(n)
Dls=np.zeros(n)

ell[0:n1] = Cls_TTLOLUNB[0]
ell[n1:n] = Cls_TTHILUNB[0]
Dls[0:n1] = Cls_TTLOLUNB[1]
Dls[n1:n] = Cls_TTHILUNB[1]

pyplot.plot(ell, Dls,linewidth=2.0, color="red") 
pyplot.title('Angular Power Spectrum',fontsize=20)
pyplot.xlabel('Multipole, $\ell$',fontsize=20)
pyplot.ylabel('$[\ell(\ell + 1)/2\pi] C_\ell^{TT}$ $[\mu K^2]$',fontsize=20) # ($\mu K^2$)')
fig = pyplot.gcf()
fig.set_size_inches(10, 6)
pyplot.savefig("Cls_TT2.png")
pyplot.show()
