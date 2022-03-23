import numpy as np
from scipy import constants
from qconstants import *
import os.path
import pandas as pd
import halodp
import quser_define

############################## Choose a Cluster from list below ###############################################

cluster=quser_define.cluster_name


############################## List of Cluster ##################################################################


d1=pd.DataFrame.from_dict({
#'nan'  : ['A2142','MH',0.0909,418600,23,129,0.57,0.810,5.82,10.0,0.5,6189.32, 510,2000],

'A4038'  : ['A4038','MH',0.02819, 100000, 86,43, 3.11,0.541,0.0174, 7.95E-6,0.75,98.38]
#'A1758S'  : ['A1758S','MH', 0.279, 1428,3.1},
'A1758'  : ['A1758S','MH', 0.279, 876100,3.9,357.4881512390829, 10.40,0.7,0.37E-2,1.0,438.9050},



orient='index',columns=['Label','type','redshift', 'DL','S_nu','r_c_dum','T_gas_dum','beta','n0','B_dum','eta','r'])
############################# Set and decide the constants start ##################################################

#Hubble constants
h		=0.7		#little h.		unit :	-

#Critical density:
sigm		=0.3		#matter density.	unit :	-
sigl		=0.7		#dark energy density.	unit :	-
sig8		=0.829

#Properties of cluster:
redshift	=d1.loc[cluster,'redshift']
DL		=d1.loc[cluster,'DL']*kpc*100.0	# luminosity distance, in unit of cm.	z=0.0302
#DA		=100000*kpc*100.0
S_nu            =d1.loc[cluster,'S_nu']            #mJy
r_c_dum	=d1.loc[cluster,'r_c_dum']		#Core radius.		unit :	h50^-1 kpc, where H0=50 kms^-1Mpc^-1
T_gas_dum	=d1.loc[cluster,'T_gas_dum']		#T_gas in unit of keV
beta		=d1.loc[cluster,'beta']		#Beta parameter.	unit :	-
n0		=d1.loc[cluster,'n0']		#Thermal electron central density	#unit :	cm^-3. A&A 540, A38 (2012)
B_dum		=d1.loc[cluster,'B_dum']	#central magnetic field.		#unit :	Gauss.  A&A 540, A38 (2012)
eta		=d1.loc[cluster,'eta']		#proportional coefficient related to magnetic field profile

############################## data of rhos and rs from nfw profile ##################################################################

rhos=halodp.rhos
rs=halodp.rs

#Temperature of CMB at any redshift z. (want ot know more? See http://www.cv.nrao.edu/course/astr534/CMB.html)
T		=T0*(1.0+redshift)

#conversion if necessary
B0	=B_dum*(10**-6)
H_0	=100.0*h				#hubble constant.	unit :	km s^-1 Mpc^-1
H_0m	=H_0/(1000.0*si_pc)			#hubble constant.	unit :	s^-1
Ez	=np.sqrt(sigm*(1.0+redshift)**3+sigl)
Hz	=H_0m*Ez				#H(z). 			unit :	s^-1
rho_crit=3.0*Hz**2/(8.0*constants.pi*si_G)	#critical density.	unit :	kg m^-3

T_gas	=T_gas_dum*1000.0*si_e/si_k		#T_gas in unit of Kelvin
r_c_dum2=r_c_dum*0.5/h				#Core radius.		unit :	kpc, where h=0.673
r_c	=r_c_dum2*1000.0*si_pc*100.0		#Core radius.		unit :	cm

############################# Particle upper-limit ################################################################
Relic_density_limit =3.0e-26
Calet=3.0e-24
Boudaud=1.0e-24
Ibarra=0.9e-25

############################# Set and decide the constants end ####################################################
############################# start: parameters for nfw profile ####################################################
#readtable="dummy_a2199"
#d0=pd.read_csv(readtable,delim_whitespace=True,header=0,index_col=0)
#d0=pd.read_csv(readtable,delim_whitespace=True,header=True,names=['rho_s / rho_crit','r_s / kpc','redshift'],comment='#')
#s0=pd.DataFrame({'rho_crit/kgm^-3':rho_crit},index=d0.index)
#d0=pd.concat([d0,s0],axis=1)

#characteristic density in unit of kg m^-3	( multiplied by rho_crit )
#rho_s= rho_crit*d1.loc[cluster,'rho_s']

rho_s= rho_crit*rhos

#characteristic radius in unit of meter
r_s= kpc*rs
#r_s= kpc*d1.loc[cluster,'r_s']
r=kpc*d1.loc[cluster,'r']
############################# end: parameters for nfw profile ####################################################
