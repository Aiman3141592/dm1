from scipy import constants
import numpy as np

#constants in SI unit:
GeV	=1.0e9*constants.e	#Giga electron volt.	unit :	Joule
si_msun	=1.99e30		#solar mass.		unit :	kg
kpc	=1000.0*constants.parsec
mu	=0.61			#mean molecular weight.
T0	=2.73			#Temperature of current Universe.	unit :	K

si_e	=constants.e		#electron charge.	unit :	Coulomb
si_G	=constants.G		#Gravitational constant.unit :	m^3 kg^-1 s^-2
si_c	=constants.c		#speed of light		unit :	m s^-1
si_k	=constants.k		#Boltzmann constant.	unit :	m^2 kg s^-2 K^-1 or J K^-1
si_pc	=constants.parsec	#parsec			unit :	m
si_me	=constants.m_e		#electron mass		unit :	kg
si_mp	=constants.m_p		#proton mass		unit :	kg
si_mc2	=constants.m_e*constants.c**2	#mc2	unit :	Joule

si_sigma=constants.sigma	#Stefan's constant.	unit :	W m^-2 K^-4
si_a	=4.0*si_sigma/si_c	#Radiation density constant.	unit :	J m^-3 K^-4
epsilon_0=constants.epsilon_0	#vacuum permittivity	unit :	m^-3 kg^-1 s^4 A^2
mu_0	=constants.mu_0		#vacuum permeability	unit :	m kg s^-2 A^-2



#constant in cgs unit.
dum_e	=np.sqrt(1.0e-5*4*np.pi*constants.epsilon_0*0.01**2) # = 3.336e-10 Coulomb = 1 statC
cgs_e	=constants.e/dum_e	#electron charge.	unit :	statC
cgs_G	=constants.G*1000	#Gravitational constant.	unit :	cm^3 g^-1 s^-2
cgs_c	=constants.c*100	#speed of light.		unit: cm s^-1
cgs_k	=constants.k*1.0e7	#Boltzmann constant.	unit :	cm^2 g s^-2 K^-1 or erg K^-1
cgs_pc	=constants.parsec*100.0	#parsec			unit :	cm
cgs_me	=constants.m_e*1000.0	#electron mass		unit :	g
cgs_mp	=constants.m_p*1000.0	#proton mass		unit :	g
cgs_mc2	=cgs_me*cgs_c**2	#mc2			unit :	erg

cgs_sigma=constants.sigma*1.0e+3	#Stefan's constant.	unit :	erg cm^-2 K^-4 s^-1
cgs_a	=4.0*cgs_sigma/cgs_c		#Radiation density constant.	unit :	erg cm^-3 K^-4




