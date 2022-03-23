######this is code from colossus#####

from __future__ import print_function 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from colossus.cosmology import cosmology
cosmo = cosmology.setCosmology('planck15');
import quser_define

from colossus.halo import profile_nfw

cluster=quser_define.cluster_name

#cluster data from https://arxiv.org/pdf/1510.01961.pdf
d1=pd.DataFrame.from_dict({
'A4038'  : [9.62E14,12.46,0.03],
'A1758'  : [1.71E14,3.91,0.279] #south
},
orient='index',columns=['mvir','cvir','z'])

#cluster mvir-cvir and redshift
Mvir = d1.loc[cluster,'mvir']
cvir = d1.loc[cluster,'cvir']
z = d1.loc[cluster,'z']

#cluster parameter definitions
p_nfw = profile_nfw.NFWProfile(M = Mvir, c = cvir, z = z, mdef = 'vir')

#density and radius
r = 10**np.arange(0,4,0.02)
rho_m = cosmo.rho_m(z)
rho_nfw = p_nfw.density(r)

h=0.7     #Hubble Constant m/s/kpc

#rhos=(p_nfw.par['rhos']*6.76977E-29)/(h**2)
rhos=p_nfw.par['rhos']
#rs=p_nfw.par['rs']*h
rs=p_nfw.par['rs']


