import numpy as np
import scipy as sp
from qconstants import *

import qproperties
#reload(qproperties)

############################# Set and decide the constants start ##################################################
#chingyee, decide this:
redshift	=qproperties.redshift	#Redshift
r_c		=qproperties.r_c	#Core radius.		unit :	cm, where H0=67.3
n0		=qproperties.n0		#Thermal electron central density	unit :	cm^-3. Ref: Chen at al (2007)
B0		=qproperties.B0		#central magnetic field.		unit :	Gauss. Estimate only.
beta		=qproperties.beta	#Beta parameter.			unit :	-
eta		=qproperties.eta	#proportional coefficient related to magnetic field profile
T		=qproperties.T		#Temperature of CMB at any redshift z.

############################# Set and decide the constants end ####################################################
############################# b(KE) start #########################################################################
def n_gas(r):
	'''Thermal electron density.
	r = position from cluster center, in unit of cm.

	n_gas = Thermal electron density, in unit of cm^-3.'''

#	return n0				#constant
	return n0*(1.0+(r/r_c)**2)**(-1.5*beta)	#Radial dependent

def bsyn(KE,B):
	'''Energy loss rate due to synchroton radiation.
	KE = Kinetic energy of electron, in unit of GeV.
	B  = Magnetic field, in unit of Gauss.
	bsyn(KE,B) = synchrotron energy loss rate, in unit of GeV s^-1'''

	#Thompson cross section	, unit : cm^2
	sigmat=8.0*np.pi/3.0*(cgs_e**2/(cgs_me*cgs_c**2))**2

	#Energy density		, unit : erg cm^-3
	#! B in microGauss. (for input only. Later I will convert it to Gauss.)
	#! 1 Tesla = 1.0e4 Gauss
	#! 1 J m^-3 = 10 erg cm^-3

	#Lorentz factor		, unit : -
	#KE=Kinetic energy of electron or positron
	KE=KE*GeV*1.0e+7	#convert from GeV to Joule then to erg
	lorentz=(cgs_mc2+KE)/cgs_mc2

	#Power radiated by single electron, in relativistic limit, and average over pitch angle.
	#Two options. Uncomment it if you want to use that option. Remember to comment the unused one.

	#Option 1: unit in erg s^-1 / GeV s^-1
	#! 1 J s^-1 = 1.0e+7 erg s^-1
	qerg=4.0/3.0*sigmat*cgs_c*B**2/(8.0*np.pi)*lorentz**2	#unit : erg s^-1
	qgev=qerg/(GeV*1.0e+7)					#unit : GeV s^-1
	return qgev

	#Option 2: unit in s^-1 (because already divided by m_e*c**2)
#	return 4.0/3.0*sigmat/(cgs_me*cgs_c**2)*cgs_c*B**2/(8.0*np.pi)*lorentz**2


def bIC(KE):
	'''Loss rate due to inverse compton scattering of CMB.
	KE = Kinetic energy of electron, in unit of GeV.
	bIC(KE) = IC energy loss rate, in unit of GeV s^-1.

	It also depends on temperature of CMB and redshift z.'''

	#Thompson cross section	, unit : cm^2
	sigmat=8.0*np.pi/3.0*(cgs_e**2/(cgs_me*cgs_c**2))**2

	#Energy density		, unit : erg cm^-3
	#! 1 J m^-3 = 10 erg cm^-3
	Ucmb=cgs_a*(T**4)

	#Lorentz factor		, unit : -
	#KE=Kinetic energy of electron or positron
	KE=KE*GeV*1.0e+7	#convert from GeV to Joule then to erg
	lorentz=(cgs_mc2+KE)/cgs_mc2

	#Ultrarelativistic e- lose energy because low-energy CMB photon hits him.
	#Two options. Uncomment it if you want to use that option. Remember to comment the unused one.

	#Option 1: unit in erg s^-1 / GeV s^-1
	#! 1 J s^-1 = 1.0e+7 erg s^-1
	qerg=4.0/3.0*sigmat*cgs_c*Ucmb*lorentz**2	#unit : erg s^-1
	qgev=qerg/(GeV*1.0e+7)				#unit : GeV s^-1
	return qgev

	#Option 2: unit in s^-1 (because already divided by m_e*c**2)
#	return 4.0/3.0*sigmat/(cgs_me*cgs_c**2)*cgs_c*Ucmb*lorentz**2

def bcoul(KE,r):
	'''Loss rate due to Coulomb loss.
	KE = Kinetic energy of electron, in unit of GeV.
	r = position from cluster center, in unit of cm.
	bIC(KE,r) = Coulomb energy loss rate, in unit of GeV s^-1.
	it's log or ln?! Answer is ln (maybe only)
	1.2 or 6.13??? Answer is 6.13! Well... actually depends on you do b(E) or b(lorentz).'''

	#Lorentz factor		, unit : -
	#KE=Kinetic energy of electron or positron
	KE=KE*GeV		#convert from GeV to Joule
	lorentz=(si_mc2+KE)/si_mc2

	# unit in GeV s^-1 (instead of b(lorentz) which is in unit of s^-1.)
	qg=6.13e-16*n_gas(r)*(1.0+np.log(lorentz/n_gas(r))/75.0)
	# unit in J s^-1
	qg*GeV
	# unit in erg s^-1
	qg*GeV*1.0e7
	# unit in GeV s^-1
	return qg

def bbrem(KE,r):
	'''Loss rate due to bremsstrahlung.
	KE = Kinetic energy of electron, in unit of GeV.
	r = position from cluster center, in unit of cm.
	bbrem(KE,r) = bremsstrahlung energy loss rate, in unit of GeV s^-1.
	it's log or ln?!
	why different paper seems to give different things?'''

	#Lorentz factor		, unit : -
	#KE=Kinetic energy of electron or positron
	KE=KE*GeV		#convert from GeV to Joule
	lorentz=(si_mc2+KE)/si_mc2

	# unit in s^-1 (because this is the b(lorentz), instead of b(E), which means that is already divided by m_e*c**2)
	qg=1.51e-16*n_gas(r)*lorentz*(np.log(lorentz)+0.36)
	# unit in J s^-1
	qg*si_mc2
	# unit in erg s^-1
	qg*cgs_mc2
	# unit in GeV s^-1
	return qg*si_mc2/GeV

def btot(KE,B,r):
	'''Total loss rate (sum of the four loss rate above) in spherical coord.
	KE = Kinetic energy of electron, in unit of GeV.
	B = Magnetic field, in unit of Gauss.
	r = position from cluster center, in unit of cm.

	btot = total energy loss rate, in unit of GeV s^-1.'''

	A=bsyn(KE,B)+bIC(KE)+bcoul(KE,r)+bbrem(KE,r)	#in unit of GeV s^-1
	return A

def btot_sph(KE,B,r):
	'''Total loss rate in spherical coord.
	KE	= Kinetic energy of electron, in unit of GeV.
	B	= Magnetic field, in unit of Gauss.
	r	= spherical radius	, in unit of cm;

	btot_sph = total energy loss rate, in unit of erg s^-1.'''

	A=bsyn(KE,B)+bIC(KE)+bcoul(KE,r)+bbrem(KE,r)	#in unit of GeV s^-1
	return A*GeV*1.0e+7	#from GeV s^-1 to erg s^-1

############################# b(KE) end ###########################################################################
############################# needed in integration start #########################################################
def rayB(r,typ):
	'''Magnetic field in Gauss.
	r	= spherical radius	, in unit of cm;
	typ	= 'const', 'radial'

	rayB	= magnetic field in unit of Gauss.'''

	if typ in ('const','CONST'):
		return B0
	elif typ in ('radial','rad','RADIAL','RAD'):
		k=B0/(n0**eta)			# Proportional coefficient. Ref: feretti et al 2012 observational properties.pdf
		return k*(n_gas(r))**eta	# Estimate only.


def funcR(nu,KE,B):
	'''Function defined in (Ghisellini, Guilbert, & Svensson 1988).
	It is needed in the integration step to get radio power.

	nu = Radio frequency, in unit of GHz.
	KE = Kinetic energy of electron, in unit of GeV.
	B  = Magnetic field, in unit of Gauss.'''

	#Photon frequency
	nu=nu*1.0e9		#convert from GHz to Hz

	#Lorentz factor
	#KE=Kinetic energy of electron or positron
	KE=KE*GeV		#convert from GeV to Joule
	lorentz=(si_mc2+KE)/si_mc2

	#Magnetic field
	#! 1 Gauss = 1.0e-4 Tesla

	#Cyclotron frequency
	vB=cgs_e*B/(2.0*np.pi*cgs_me*cgs_c)

	#Normalized frequency variable
	x=nu/(3.0*lorentz**2*vB)

	#R(x)
	a1=sp.special.kv(4.0/3.0,x)*sp.special.kv(1.0/3.0,x)
	a2=sp.special.kv(4.0/3.0,x)**2-sp.special.kv(1.0/3.0,x)**2
	return (2.0*x**2)*(a1-3.0/5.0*x*a2)

def AA():
	'''Lump of constants, in radio power formula.'''
	#constant in cgs unit.
	return np.sqrt(3.0)*cgs_e**3/cgs_mc2


############################# needed in integration end ###########################################################
############################# start: nfw profile #######################################################
def nfw_rho(r,rhos,rs):
	'''Dark matter density, in spherical coord.
	r 	= spherical radius, in unit of meter;
	nfw_rho	= dark matter density, in unit of kg m^-3.'''
	A=(r/rs)*(1.0+(r/rs))**2
	return rhos/A


def nfw_numdens_kk(r,m_neutralino,rhos,rs):
	'''Dark matter particle number density.
	r		= radius, in unit of cm;
	m_neutralino	= DM particle mass, in unit of GeV;
	nfw_numdens_neu	= DM particle number density at position r, in unit of cm^-3.'''
	m_neutralino=m_neutralino*GeV/(si_c**2)		#convert GeV to kg
	r=r/100.0						#Convert r from cm to meter
	return nfw_rho(r,rhos,rs)/(m_neutralino*100.0**3)	#neutralino number density in unit of cm^-3

############################# end: nfw profile #######################################################
