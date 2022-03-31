'''
=====================================================================================
 Generate e- source spectrum (Q_source) from PYTHIA dN/dE, which is equivalent to 
 dshayield spectrum from Darksusy.

 KK or SUSY DM, up to you (it should be KK DM here..?). As long as it is from PYTHIA.

 Combination of different KK anni channels (see hooper and profumo 2007 - UED KK.pdf)

 The source spectrum is for the special case:
 sigmav=1.0	#sigmav (cm^3 s^-1)
 BR    =1.0	#BR which shouldn't be modified anymore for KK case, since I have considered BR in Pythia dN/dE
 Nx    =1.0	#DM number density (cm^-3).
=====================================================================================
'''
print(__doc__)
import numpy as np
import warnings
import qproperties
import qfunctions
import qconstants
 
def scndsource(i,kkfile,filename):
   warnings.simplefilter('always',Warning)


###############getting DM mass from maddm.out#########################



		
############################### start: define here ###################################

# Reading from the file to get mass of DM
   fileout = open(filename, 'r')
   content = fileout.readlines()

# Iterating through the content
# Of the file
   numbers = []
   for line in content:
	   for word in line.split():
  	       if word.lstrip("-").replace('.', '', 1).replace('E-', '', 1).replace('E', '', 1).isdigit():
      	          numbers.append(float(word))

   m_kk = numbers[2] # in GeV
   	    
############################### end: define here ###################################

   print "Read dN/dE from file '%s'."%kkfile

   x =np.loadtxt(kkfile)[:,0]
   y =np.loadtxt(kkfile)[:,1]		#unit: # annihilation^-1 Gev^-1
   dNdE=(10**y)*m_kk
   Gev=(10**x)*m_kk

   rhos	=qproperties.rho_s
   rs	=qproperties.r_s
   r       =qproperties.r

#Q_source = sigmav * dN/dE * BR * (Nx**2) / 2.0'
#Please give your sigmav and Nx:'
   sigmav=numbers[31]	#unit: cm^3 s^-1   for number[i], mummup=32 or bbx=30 
   BR    =1.0	                #since its 100% so its equal to 1.0.
   Nx    =qfunctions.nfw_numdens_kk(r,m_kk,rhos,rs)	        #Neutralino number density (cm^-3)

   Qsource=np.array(sigmav*dNdE*BR*(Nx**2/2.0))	#unit: # s^-1 GeV^-1 cm^-3

   dummytable =np.array([[m_kk,sigmav]])
   line1='  E                           Q_source\n  GeV                         # s^-1 GeV^-1 cm^-3\n #.....mchi..... ...sigmav...'


   latesttable=np.array([Gev,Qsource]).transpose()
   latesttable=np.concatenate((dummytable,latesttable),axis=0)

   np.savetxt(outfile,latesttable,'%25.15e', header=line1)


#### iterating through all output#####

n_outputfiles= 44 #number of output file
for i in range (1,n_outputfiles+1,1):
     if i<10:
         outfile_0	=".source_run_01_0%s"%i
         outfile_1	="cluster_%s"%qproperties.cluster 
         outfile	=outfile_1+outfile_0
         kkfile		="run_01_0%s/positrons_spectrum_PPPC4DMID_ew.dat"%i
         filename	="run_01_0%s/maddm.out"%i
         
         scndsource(i,kkfile,filename)
         
         print "\nSource spectrum written to file '%s'."%outfile
     else:
         outfile_0	=".source_run_01_%s"%i
         outfile_1	="cluster_%s"%qproperties.cluster 
         outfile	=outfile_1+outfile_0
         kkfile		="run_01_%s/positrons_spectrum_PPPC4DMID_ew.dat"%i
         filename	="run_01_%s/maddm.out"%i
         
         scndsource(i,kkfile,filename)
         
         print "\nSource spectrum written to file '%s'."%outfile    
         
