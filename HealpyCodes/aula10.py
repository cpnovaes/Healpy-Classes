#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 11:25:28 2017

@author: camila

  -- Aula 10 --

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
STEP 2: Spherical harmonic transforms: Cls --> map
"""

Cls = hp.read_cl('Cls_bestfitLCDM_PLA2_TT_lmax2508.fits')

# Exemplo:
mapa1 = hp.synfast(Cls,1024, alm=False, pol=False, pixwin=False, fwhm=np.deg2rad(0.))
len(mapa1)
mapa1 = hp.synfast(Cls,1024, alm=True, pol=False, pixwin=False, fwhm=np.deg2rad(0.))
len(mapa1)
len(mapa1[0]), len(mapa1[1])

# Ex1:
mapa1 = hp.synfast(Cls,512, alm=False, pol=False, pixwin=False, fwhm=np.deg2rad(0.))
mapa2 = hp.synfast(Cls,512, alm=False, pol=False, pixwin=False, fwhm=np.deg2rad(0.25))
mapa3 = hp.synfast(Cls,512, alm=False, pol=False, pixwin=True, fwhm=np.deg2rad(0.))

hp.mollview(mapa1, title='fwhm = 0')
pyplot.savefig("cmb_fwhm0.png")
hp.mollview(mapa2, title='fwhm = 0.25deg')
pyplot.savefig("cmb_fwhm0p25.png")

Cls1 = hp.anafast(mapa1)
Cls2 = hp.anafast(mapa2)
Cls3 = hp.anafast(mapa3)

lmax = len(Cls1)
ell = np.arange(lmax)

Dls1 = ell*(ell+1)*Cls1/(2.*np.pi)
Dls2 = ell*(ell+1)*Cls2/(2.*np.pi)
Dls3 = ell*(ell+1)*Cls3/(2.*np.pi)

pyplot.plot(ell, Dls1,linewidth=2.0, color="red",label="fwhm=0/pixwin=False") 
pyplot.plot(ell, Dls2,linewidth=2.0, color="blue",label="fwhm=0.25/pixwin=False") 
pyplot.plot(ell, Dls3,linewidth=2.0, color="green",label="fwhm=0/pixwin=True") 

pyplot.title('Angular Power Spectrum',fontsize=20)
pyplot.xlabel('Multipole, $\ell$',fontsize=20)
pyplot.ylabel('$[\ell(\ell + 1)/2\pi] C_\ell^{TT}$ $[\mu K^2]$',fontsize=20) # ($\mu K^2$)')
pyplot.legend(loc='best')
fig = pyplot.gcf()
fig.set_size_inches(10, 6)
pyplot.savefig("Cls_smo0p25Dif.png")
pyplot.show()

# Ex2:
    
mapa1 = hp.synfast(Cls,512)
mapa2 = hp.synfast(Cls,512)
mapa3 = hp.synfast(Cls,512)

Cls1 = hp.anafast(mapa1)
Cls2 = hp.anafast(mapa2)
Cls3 = hp.anafast(mapa3)

lmax = len(Cls1)
ell = np.arange(lmax)

Dls = ell*(ell+1)*Cls[0:lmax]/(2.*np.pi)
Dls1 = ell*(ell+1)*Cls1/(2.*np.pi)
Dls2 = ell*(ell+1)*Cls2/(2.*np.pi)
Dls3 = ell*(ell+1)*Cls3/(2.*np.pi)

pyplot.plot(ell, Dls1,linewidth=1.0, color="red") 
pyplot.plot(ell, Dls2,linewidth=1.0, color="blue") 
pyplot.plot(ell, Dls3,linewidth=1.0, color="green") 
pyplot.plot(ell, Dls,linewidth=2.0, color="black") 
pyplot.xscale('log')

pyplot.title('Angular Power Spectrum',fontsize=20)
pyplot.xlabel('Multipole, $\ell$',fontsize=20)
pyplot.ylabel('$[\ell(\ell + 1)/2\pi] C_\ell^{TT}$ $[\mu K^2]$',fontsize=20) # ($\mu K^2$)')
fig = pyplot.gcf()
fig.set_size_inches(10, 6)
pyplot.savefig("Cls_dif_to_input.png")
pyplot.show()


#%%
"""
STEP 3: Spherical harmonic transforms: map --> alms
"""

mapa = hp.synfast(Cls,512, alm=True)
mapa1 = mapa[0]
alm1 = mapa[1]

alm2 = hp.map2alm(mapa1)

alm2 = alm1

#%%
"""
STEP 4: Spherical harmonic transforms: alms --> map
"""

mapa1 = hp.synfast(Cls,512, alm=True, pol=False, pixwin=False, fwhm=np.deg2rad(0.))
alms = mapa1[1]

mapa2 = hp.alm2map(alms, 512, pixwin=True)

hp.mollview(mapa1[0], title="WITHOUT pixwin")
pyplot.savefig("mapa1.png")
hp.mollview(mapa2, title="WITH pixwin")
pyplot.savefig("mapa1_pixwin.png")

hp.mollview(mapa2-mapa1[0], title="WITH pixwin")
pyplot.savefig("dif_pixwin.png")

#%%
"""
STEP 5: Spherical harmonic transforms: alms --> map + 1st derivatives
"""

mapa1 = hp.synfast(Cls,512, alm=True, pol=False, pixwin=False, fwhm=np.deg2rad(0.))
alms = mapa1[1]

mapa2 = hp.alm2map_der1(alms, 512)

len(mapa2)

hp.mollview(mapa1[0], title="mapa1[0]")
pyplot.savefig("mapa1_ampl_col0.png")

hp.mollview(mapa2[0], title="mapa2[0]")
pyplot.savefig("mapa2_ampl_col0.png")

hp.mollview(mapa2[1], title="der. theta")
pyplot.savefig("mapa2_derTheta_col1.png")

hp.mollview(mapa2[2], title="der. phi")
pyplot.savefig("mapa2_derPhi_col2.png")

#%%
"""
STEP 6: Spherical harmonic transforms: alms --> fits file
"""
# syntax: healpy.fitsfunc.write_alm(filename, alms, out_dtype=None, lmax=-1, mmax=-1, mmax_in=-1)

mapa1 = hp.synfast(Cls,512, alm=True, pol=False, pixwin=False, fwhm=np.deg2rad(0.))
alms = mapa1[1]

hp.write_alm('alm.fits', alms)

#%%
"""
STEP 7: Spherical harmonic transforms: fits file --> alms
"""
# syntax: healpy.fitsfunc.read_alm(filename, hdu=1, return_mmax=False)

alms1 = hp.read_alm('alm.fits', return_mmax=False)
alms2 = hp.read_alm('alm.fits', return_mmax=True)   # return_mmax=True -> include the 
                                                    # lmax=mmax value in the output.