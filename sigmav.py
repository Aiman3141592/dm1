import numpy as np 
import matplotlib.pyplot as plt 
import qproperties
import qfunctions
import kkrunthiswithoutloss
import math
from scipy import stats
 
 
#input data	
kkfile =kkfile =kkrunthiswithoutloss.outfilename


x = np.loadtxt(kkfile)[:,1]
y = (np.loadtxt(kkfile)[:,6]/np.loadtxt(kkfile)[:,0])

y_ob = qproperties.S_nu
 
#define coefficient
s=kkrunthiswithoutloss.n_outputfiles
b_2=np.random.normal(scale=1.0, size=s)

p = np.polyfit(np.log(x+b_2), np.log(y), 1)


print (p[0],p[1])

#list of eq.
A=np.log(y_ob)	
B=np.log(x+b_2)
C=A-p[0]*B-p[1]
    	
#sigmav-m_kk constraint sigmav=10*(log(observed radio emission)-slope of line*log(x+np.random)-(Y-intercept)
y_upper = (10**(C))


#plotting the actual points as scatter plot
plt.scatter(x,y_upper)

# putting labels 
plt.xlabel(r'Mass (GeV)') 
plt.ylabel(r'Cross-section (cm$^{3}$ s$^{-1}$)')
	
# turn to log
plt.yscale('log')
#plt.xscale('log')

#save
line1= 'mass   upper(sigmavtot    sigmavBB       sigmavUU)'
np.savetxt('Sig_%s.out'%qproperties.cluster, np.c_[x,y_upper],'%25.15e',header=line1)

# Putting title
plt.title(r'$\sigma_v-M_{LKP}$ constraint for %s'%qproperties.cluster)
	
# function to show plot 
plt.legend()
imgname='%s'%kkfile+'_sigmav.png'
plt.savefig(imgname, dpi='figure', format='png')

	
	
