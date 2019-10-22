#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 18:45:48 2017

@author: camila

            -- AULA 3 --
"""

#%%
"""
STEP 1: Import the packages
"""

import healpy as hp
import matplotlib
import numpy as np
from matplotlib import pyplot # ***

#%%
"""
STEP 2: Read the map

First read the larger fits file ... then the smaller.

Questions: 
    1. how many elements there are in the fits file?
    2. what are they?
    3. ring or nest?
    
"""

mapa = hp.read_map('LFI_CompMap_Foregrounds-smica_1024_R2.00.fits', field={0,1,2}, h=True, nest=True)

# WHAT IS A MAP IN PYTHON?????

#mapa_ex = list(range(192))
#hp.mollview(mapa_ex)             # it DOES NOT work
#hp.mollview(np.array(mapa_ex))   # it DOES work
#
## Print the firs 10 elements:
#print(list(mapa[0][0:10]))
#
## To print the header:
#print(list(mapa[3]))



#%%
"""
STEP 3: Visualize the maps and save the figure
"""
fontsize = 20
matplotlib.rcParams.update({'font.size':fontsize})

hp.mollview(mapa[0], nest=True, norm='hist', unit='K', coord=['G','E'])
pyplot.savefig('/media/camila/Dados/ON/Disciplina_PG_ON_2017/Aulas/Material/mollweide_view.png')

hp.mollview(mapa[0], nest=True, norm='hist', rot=[0,210,70], title='Rotated map', unit='K', coord=['G','E'])
pyplot.savefig('/media/camila/Dados/ON/Disciplina_PG_ON_2017/Aulas/Material/mollweide_view2.png')

hp.mollview(mapa[0], nest=True, norm='hist', rot=[0,210,70], title='Rotated map', unit='K', coord=['G','E'], min=0)
pyplot.savefig('/media/camila/Dados/ON/Disciplina_PG_ON_2017/Aulas/Material/mollweide_view3.png')

#%%
"""
STEP 4: Manipulate and Write the map to a fits file

"""
mapa2 = mapa[0]#*1e6

hp.mollview(mapa2, nest=True, norm='hist', title='Foreground map', unit='K')
pyplot.savefig('/media/camila/Dados/ON/Disciplina_PG_ON_2017/Aulas/Material/mollweide_view4.png')

