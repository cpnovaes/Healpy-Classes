#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 10:47:55 2017

@author: camila

  -- Aula 11 --

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
STEP 2: Spherical harmonic transforms tools: Alms()
"""

map_in = hp.read_map('COM_CMB_IQU-smica_1024_R2.02_full.fits')
alm = hp.map2alm(map_in)

l_max = hp.Alm.getlmax(len(alm))

size = hp.Alm.getsize(l_max)
print(len(alm), '/', size)


l, m = hp.Alm.getlm(l_max) # getlm(l_max, index)

index = hp.Alm.getidx(l_max) # getidx(l_max, l, m)

########
# Ex1:::
alm = hp.map2alm(map_in)

lmax = hp.Alm.getlmax(len(alm))

l, m = hp.Alm.getlm(lmax)

l1 = 2

# To take the alm's corresponding to just one multipole, i.e., when l = l_min = l_max
alm[np.where(l != l1)] = (0+0j) # alm is set to zero whenever l != l1
                                # 0+0j --> alms are complex numbers

Nside = 128
map_out1 = hp.sphtfunc.alm2map(alm, Nside, verbose=True) # Reconstruct the map according to 
                                                         # the new alm set.

########
# Ex2:::
alm_in = hp.map2alm(map_in)

lmax = hp.Alm.getlmax(len(alm))

l_m = hp.sphtfunc.Alm.getlm(lmax)
l = l_m[0]

l_min = 10
l_max = 150

# To take the alm's corresponding to an interval of multipolez, i.e., when l_min != l_max
alm = alm_in
alm[l < l_min] = (0+0j)  # alm is set to zero whenever l < l_min and l > l_max
alm[l > l_max] = (0+0j)

Nside = 128
map_out2b = hp.sphtfunc.alm2map(alm, Nside, verbose=True) # Reconstruct the map according to 
               
hp.mollview(map_out2b)
pyplot.savefig("cmb_map_10l150.png")

#%%
"""
STEP 3: Spherical harmonic transforms tools: synalm, alm2cl, almxfl
"""

Cls = hp.read_cl('Cls_bestfitLCDM_PLA2_TT_lmax2508.fits')

####
alm1 = hp.synalm(Cls,lmax=1535) # = synfast, but generating the alm's.
mapa1 = hp.alm2map(alm1,512)
hp.mollview(mapa1, title='Mapa1')
pyplot.savefig("mapa1_synalm_ex.png")

####
Cls_calc = hp.alm2cl(alm1) # = anafast, but upon alm's.

lmax = len(Cls_calc)
ell = np.arange(lmax)

Dls = ell*(ell+1)*Cls[0:lmax]/(2.*np.pi)
Dls_calc = ell*(ell+1)*Cls_calc/(2.*np.pi)

pyplot.plot(ell, Dls_calc,linewidth=2.0, color="red",label="Calculated") 
pyplot.plot(ell, Dls,linewidth=2.0, color="black",label="Best-fit") 
pyplot.title('Angular Power Spectrum',fontsize=20)
pyplot.xlabel('Multipole, $\ell$',fontsize=20)
pyplot.ylabel('$[\ell(\ell + 1)/2\pi] C_\ell^{TT}$ $[\mu K^2]$',fontsize=20) # ($\mu K^2$)')
pyplot.legend(loc='best')
fig = pyplot.gcf()
fig.set_size_inches(10, 6)
pyplot.savefig("Cls_alm2cl_ex.png")
pyplot.show()


####
fl = ell[::-1]**7. #*(ell+1)/(2.*np.pi)

alm2 = hp.almxfl(alm1,fl)

mapa2 = hp.alm2map(alm2,512)
hp.mollview(mapa2, title='Mapa2 = Mapa1 x fl')
pyplot.savefig("mapa2_almxfl_ex.png")

#%%
"""
STEP 4: Spherical harmonic transforms tools: smoothing, smoothalm
"""

mapa = hp.read_map('COM_CMB_IQU-smica_1024_R2_02_full.fits')
hp.mollview(mapa, title='CMB map')
pyplot.savefig("mollview_COM_CMB_IQU-smica_1024_R2.02_full.png")

mapa_smo = hp.smoothing(mapa, fwhm=np.deg2rad(0.75))

hp.mollview(mapa_smo, title='Smoothed map')
pyplot.savefig("mapa_smoothing_ex.png")

alm = hp.map2alm(mapa)

alm_smo = hp.smoothalm(alm, fwhm=np.deg2rad(0.75))

mapa2_alm_smo = hp.alm2map(alm_smo, 1024)

hp.mollview(mapa2_alm_smo, title='Smoothed alm')
pyplot.savefig("mapa_smoothalm_ex.png")

hp.mollview(mapa2_alm_smo-mapa_smo, title='Difference')


#%%
"""
STEP 5: Spherical harmonic transforms tools: pixwin
"""

Cls = hp.read_cl('Cls_bestfitLCDM_PLA2_TT_lmax2508.fits')

# Exercise:
mapa2 = hp.synfast(Cls,512, pixwin=True)
Cls_calc2 = hp.anafast(mapa2)

# Example: The pixwin func has the same effect uppon the Cls calculated from a degraded map!
#mapa2 = hp.synfast(Cls,2048)
#mapa2 = hp.ud_grade(mapa2, 512)
#Cls_calc2 = hp.anafast(mapa2)

lmax = len(Cls_calc2)-1

Nside = 512
p_func = hp.pixwin(Nside)

Cls_norm2 = Cls_calc2
Cls_norm2 = Cls_norm2/p_func[:lmax+1]**2.



# Plots:::
ell = np.arange(lmax+1)

#Dls = ell[0:2509]*(ell[0:2509]+1)*Cls/(2.*np.pi)
Dls = ell*(ell+1)*Cls[0:lmax+1]/(2.*np.pi)

Dls_calc2 = ell*(ell+1)*Cls_calc2/(2.*np.pi)
Dls_norm2 = ell*(ell+1)*Cls_norm2/(2.*np.pi)

Dls_pfunc = ell*(ell+1)*p_func[0:lmax+1]/(2.*np.pi)

pyplot.plot(ell, Dls_calc2,linewidth=2.0, color="red",label="Calculated") 
pyplot.plot(ell, Dls_norm2,linewidth=2.0, color="blue",label="Normalized") 
pyplot.plot(ell, Dls,linewidth=2.0, color="black",label="Best-fit") 
pyplot.plot(ell, p_func[0:lmax+1],linewidth=2.0, color="black",label="Best-fit") 
pyplot.legend(loc='best')
pyplot.title('Angular Power Spectrum',fontsize=20)
pyplot.xlabel('Multipole, $\ell$',fontsize=20)
pyplot.ylabel('$[\ell(\ell + 1)/2\pi] C_\ell^{TT}$ $[\mu K^2]$',fontsize=20) # ($\mu K^2$)')
pyplot.legend(loc='best')
fig = pyplot.gcf()
fig.set_size_inches(10, 6)
#pyplot.savefig("Cls_pixwin_ex2.png")
pyplot.show()

#%%
####
mapa1 = hp.synfast(Cls,1024, pixwin=False, fwhm=np.deg2rad(0.25))
Cls_calc1 = hp.anafast(mapa1)
lmax = len(Cls_calc1)-1

fwhm = np.deg2rad(0.25)
g_func = hp.gauss_beam(fwhm, lmax=lmax)
g_func2 = g_func**2.

Cls_norm = Cls_calc1
Cls_norm[g_func2 != 0.] = Cls_norm[g_func2 != 0.]/g_func2[g_func2 != 0.]

Dls_calc = ell*(ell+1)*Cls_calc/(2.*np.pi)
Dls_norm = ell*(ell+1)*Cls_norm/(2.*np.pi)


