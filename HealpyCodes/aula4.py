#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 18:45:48 2017

@author: camila

            -- AULA 4 --
"""

#%%
"""
STEP 1: Import the packages
"""

import healpy as hp
import matplotlib
from matplotlib import pyplot # ***

#%%
"""
STEP 2: Previously ...

First read the fits file.

"""

mapa = hp.read_map('LFI_CompMap_Foregrounds-smica_1024_R2.00.fits', field={0,1,2})

hp.write_map('new_map.fits',mapa)

fontsize = 20
matplotlib.rcParams.update({'font.size':fontsize})

#%%
"""
STEP 3: Include graticule in a mollweide projection

"""
hp.mollview(mapa[0], norm='hist', title='Including graticule', unit='K', coord=['G','E'])
hp.graticule(dpar=10.,dmer=25.,coord=['G','E'], local=True)
pyplot.savefig('/media/camila/Dados/ON/Disciplina_PG_ON_2017/Aulas/Material/mollweide_view5.png')

#%%
"""
STEP 4: Visualize a map in the other 3 projection types

"""
hp.gnomview(mapa[0],rot=[5,10,0], title='Gnomonic proj.', reso=5., norm='hist', unit='K')
hp.graticule(dpar=5., local=True)
pyplot.savefig('/media/camila/Dados/ON/Disciplina_PG_ON_2017/Aulas/Material/gnomonic_view1.png')

hp.gnomview(mapa[0],rot=[5,10,90], title='Gnom. proj.- rotation in psi', reso=5., norm='hist', unit='K')
hp.graticule(dpar=5., local=True)
pyplot.savefig('/media/camila/Dados/ON/Disciplina_PG_ON_2017/Aulas/Material/gnomonic_view2.png')

#******

hp.cartview(mapa[0], norm='hist', unit='K')
hp.graticule(dpar=10.,dmer=25.)
pyplot.savefig('/media/camila/Dados/ON/Disciplina_PG_ON_2017/Aulas/Material/cartesian_view1.png')

hp.cartview(mapa[0], norm='hist',coord=['G','E'], unit='K')
hp.graticule(dpar=10.,dmer=25.,coord=['G','E'])
pyplot.savefig('/media/camila/Dados/ON/Disciplina_PG_ON_2017/Aulas/Material/cartesian_view2.png')

#******

hp.orthview(mapa[0], norm='hist', unit='K')
hp.graticule(dpar=10.,dmer=25.)
pyplot.savefig('/media/camila/Dados/ON/Disciplina_PG_ON_2017/Aulas/Material/orthographic_view1.png')

hp.orthview(mapa[0], norm='hist',coord=['G','E'], unit='K')
hp.graticule(dpar=10.,dmer=25.,coord=['G','E'])
pyplot.savefig('/media/camila/Dados/ON/Disciplina_PG_ON_2017/Aulas/Material/orthographic_view2.png')

#%%
"""
STEP 5: Write a healpix fits map

"""

mapa[0][:] = mapa[0][:] * 1e6

hp.write_map('new_map.fits',mapa, coord='G', column_names=list(['30GHz','44GHz','70GHz']), column_units=list(['microK','K','K']))

mapa2 = hp.read_map('new_map.fits', field={0,1,2}, h=True)


#%%
"""
STEP 6: Reading/Writing generic data of/into an fits file

"""

#####################################
data = hp.mrdfits('COM_PCCS_030_R2.04.fits')

# OU AINDA:
    
from astropy.io import fits
hdulist = fits.open('COM_PCCS_030_R2.04.fits')
hdulist.info()
data = hp.mrdfits(hdulist, hdu=1)

hdulist.writeto('new_table.fits')

#####################################

len(data)

hp.mwrfits('new_table.fits', data)

#####################################
# Writing FITS files.
#The creation of a FITS file pass through 4 steps.
#
#1) Creation of numpy array with the data.
#import pyfits # ou:::
import astropy.io.fits as pyfits
import numpy as np
    
x = data[10] #np.arange(100)
  
#2) Creation of the HDU from the data.

hdu = pyfits.PrimaryHDU(x)
   
#thus created, the hdu has its basic header and the data.
#
#print hdu.header.ascardlist()
#SIMPLE  =                    T / conforms to FITS standard                     
#BITPIX  =                   64 / array data type                               
#NAXIS   =                    1 / number of array dimensions                    
#NAXIS1  =                  100                                                 
#EXTEND  =                    T                                                 

#3) Once all the keywords are ready, the final HDU list have to be created and written to the file:

hdulist = pyfits.HDUList([hdu])
hdulist.writeto('new.fits')
hdulist.close()

#4) To read it again:

hdulist = pyfits.open('new.fits')
data2 = hdulist[0].data

#####################################

#%%

import matplotlib.pyplot as pl

mapa = hp.read_map('LFI_CompMap_Foregrounds-smica_1024_R2.00.fits')
hp.zoomtool.mollzoom(mapa, hold=True)
pl.show()

