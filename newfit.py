#this is code for log-log model

import qproperties
import qfunctions
import kkrunthiswithoutloss
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#import file
kkfile =kkrunthiswithoutloss.outfilename


#input data
x = np.loadtxt(kkfile)[:,1]
y = ((np.loadtxt(kkfile)[:,6])/(np.loadtxt(kkfile)[:,0]))


#define coefficient
s=kkrunthiswithoutloss.n_outputfiles
b_2 = np.random.normal(scale=1.0, size=s)
p   = np.polyfit(np.log10(x+b_2), np.log10(y), 1)



#fitting eq. (log-log model)
y_fit=10**(p[0]*np.log10(x+b_2)+p[1])


#plot fitting line
plt.plot(x, y_fit, color = "orange", label = r'$\langle\sigma v\rangle_{total}$', linestyle='solid')

#plot 95%
#sns.regplot(x,y, ci=95)


#plot real data
plt.scatter(x, y, color = "orange", marker = "o", s = 30)




#graph title
#plt.title(r'Synchrotron flux/Cross-section vs Mass for %s fitting model'%qproperties.cluster)

#graph label
plt.xlabel(r'$M_{DM}$ (GeV)')
plt.ylabel(r'$\frac{S_{v}}{\langle\sigma_v\rangle}$ (MJy/$cm^{3}$ $s^{-1}$)')


#graph scale
#plt.xscale('log')
plt.yscale('log')
    	
print(p[0],p[1])
		
plt.legend()
imgname='%s'%kkfile+'_constraintline.png'
plt.savefig(imgname, dpi='figure', format='png')








